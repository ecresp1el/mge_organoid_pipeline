"""Build reproducible, barcode-aligned surface-only Hypergate inputs.

The crude core-panel run precedes this expansion and module-label analysis.
"""
from pathlib import Path
import hashlib
import itertools
import json
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_curve

import render_div90_loupe_subcluster_guidance_expression as atlas
from div90_hypergate_paths import PROJECT_ROOT, REPO_ROOT

ROOT = REPO_ROOT
OUT = PROJECT_ROOT / "results/div90_hypergate_sst_pv"
CORE = ["ERBB4", "CXCR4", "ACKR3", "PLXNA2", "NRP1", "NRP2"]
SST_MODULE = ["SST", "NR2F2", "GRIK1"]
PV_MODULE = ["MEF2C", "MAF", "MAFB", "KCNC1", "KCNC2", "GAD1", "GAD2"]
HOUSEKEEPING = {"B2M", "TFRC", "SLC2A1", "ATP1A1", "ATP1B1", "ATP1B3",
                "SLC3A2", "SLC7A5", "LAMP1", "LAMP2", "HSPA5", "HSP90B1", "P4HB"}
ECM_ISOFORM_AMBIGUITY = {"BCAN", "FRAS1"}


def max_fscore(x, y, beta):
    best = 0.
    for sign in [1, -1]:
        precision, recall, _ = precision_recall_curve(y, sign * x)
        f = (1 + beta**2) * precision * recall / np.maximum(beta**2 * precision + recall, 1e-15)
        best = max(best, float(f.max()))
    return best


def make_jobs(panel, definition, stage):
    label_column = "crude_label" if definition == "crude" else "module_label"
    targets = ["non-SST/PV-candidate", "SST-like"] if definition == "crude" else ["PV-like", "SST-like"]
    rows = []
    for k in [1, 2, 3]:
        for combo in itertools.combinations(panel, k):
            for target in targets:
                for beta in [.5, 1., 2.]:
                    rows.append(dict(job_id=f"{stage}_{len(rows):05d}", definition=definition,
                                     label_column=label_column, target=target, beta=beta,
                                     features=";".join(combo), requested_markers=k, stage=stage))
    return pd.DataFrame(rows)


def main():
    cspa_path = OUT / "provenance/CSPA_human_surfaceome.xlsx"
    cspa = pd.read_excel(cspa_path, sheet_name="Table A")
    cspa.columns = cspa.columns.str.strip()
    high = cspa["CSPA category"].str.startswith("1")
    membrane = pd.to_numeric(cspa["Phobius TM predicted yes/no"], errors="coerce").eq(1) | pd.to_numeric(cspa.GPI, errors="coerce").eq(1)
    cspa["eligible_annotation"] = high & membrane
    cspa["gene"] = cspa["ENTREZ gene symbol"].str.strip().replace({"CXCR7": "ACKR3"})
    cspa.to_csv(OUT / "provenance/cspa_annotation_audit.tsv", sep="\t", index=False)
    curated = cspa.loc[cspa.eligible_annotation].dropna(subset=["gene"]).drop_duplicates("gene").set_index("gene")
    genes = sorted(set(curated.index) | set(CORE) | set(SST_MODULE + PV_MODULE) | {"LHX6", "PVALB"})
    expression_path = OUT / "surface_and_module_expression.tsv.gz"
    if expression_path.exists() and (OUT / "provenance/gene_matches.tsv").exists():
        expression = pd.read_csv(expression_path, sep="\t")
        matches = pd.read_csv(OUT / "provenance/gene_matches.tsv", sep="\t")
        assert set(genes).issubset(expression.columns)
    else:
        expression, matches = atlas.extract_marker_expression_from_h5ad(
            atlas.div90_spec(PROJECT_ROOT), expression_path,
            project_root=PROJECT_ROOT, genes=genes,
        )
    matches.to_csv(OUT / "provenance/gene_matches.tsv", sep="\t", index=False)
    base = pd.read_csv(OUT / "core_input.tsv.gz", sep="\t")
    metadata_cols = ["cell_id", "sample", "cluster", "loupe_label", "loupe_x", "loupe_y", "SST", "LHX6", "ERBB4", "crude_label"]
    data = base[metadata_cols].merge(expression[["cell_id", *[g for g in genes if g not in metadata_cols]]],
                                     on="cell_id", validate="one_to_one")
    assert len(data) == 4768 and data.cell_id.is_unique
    assert (data.SST > 0).sum() == 2919
    module_parameters = {}
    for name, module in [("sst", SST_MODULE), ("pv", PV_MODULE)]:
        use = [g for g in module if data[g].notna().all() and data[g].gt(0).any() and data[g].std(ddof=0) > 0]
        means = data[use].mean(); stds = data[use].std(ddof=0)
        data[name + "_module"] = ((data[use] - means) / stds).mean(axis=1)
        module_parameters[name] = dict(genes=use, means=means.to_dict(), std_population=stds.to_dict(),
                                       low=float(data[name + "_module"].quantile(.35)),
                                       high=float(data[name + "_module"].quantile(.65)))
    sm, pm = module_parameters["sst"], module_parameters["pv"]
    sst = data.sst_module.ge(sm["high"]) & data.pv_module.le(pm["low"])
    pv = data.pv_module.ge(pm["high"]) & data.sst_module.le(sm["low"])
    assert not (sst & pv).any()
    data["module_label"] = np.select([sst, pv], ["SST-like", "PV-like"], default="ambiguous")
    assert sst.sum() > 30 and pv.sum() > 30
    labels = data[[*metadata_cols, "sst_module", "pv_module", "module_label"]]
    labels.to_csv(OUT / "target_labels.tsv.gz", sep="\t", index=False)
    allowed = []
    for gene in sorted(set(curated.index) | set(CORE)):
        values = data[gene]
        reason = ""
        if gene in SST_MODULE + PV_MODULE:
            reason = "target_module_gene_excluded_to_prevent_label_leakage"
        elif gene in ECM_ISOFORM_AMBIGUITY:
            reason = "extracellular_matrix_or_surface_isoform_ambiguity"
        elif gene in HOUSEKEEPING or gene.startswith(("RPL", "RPS", "MT-", "HSP")):
            reason = "housekeeping_or_intracellular_exclusion"
        elif values.isna().any():
            reason = "gene_not_matched"
        elif values.gt(0).mean() < .05 or values.nunique() < 3:
            reason = "less_than_5pct_detection_or_insufficient_variation"
        row = dict(gene=gene, allowed=not bool(reason), exclusion_reason=reason,
                   detection_fraction=float(values.gt(0).mean()),
                   source="user_core_panel" if gene in CORE else "CSPA_2015_high_confidence_TM_or_GPI",
                   uniprot=str(curated.loc[gene, "ID_link"]) if gene in curated.index else "",
                   protein_name=str(curated.loc[gene, "UP_Protein_name"]) if gene in curated.index else "",
                   reference="https://doi.org/10.1371/journal.pone.0121314.s003")
        allowed.append(row)
    annotation = pd.DataFrame(allowed)
    annotation.to_csv(OUT / "provenance/surface_feature_filter_audit.tsv", sep="\t", index=False)
    allowed = annotation.loc[annotation.allowed].copy()
    assert set(CORE).issubset(set(allowed.gene))
    allowed.to_csv(OUT / "allowed_surface_markers.tsv", sep="\t", index=False)
    features = allowed.gene.tolist()
    data[[*labels.columns, *[g for g in features if g not in labels]]].to_csv(OUT / "hypergate_input.tsv.gz", sep="\t", index=False)
    scores = []
    panels = {}
    for definition, col in [("crude", "crude_label"), ("module", "module_label")]:
        keep = data[col].ne("ambiguous")
        targets = ["non-SST/PV-candidate", "SST-like"] if definition == "crude" else ["PV-like", "SST-like"]
        ranked = {}
        for target in targets:
            y = data.loc[keep, col].eq(target).to_numpy()
            for beta in [.5, 1., 2.]:
                vals = {g: max_fscore(data.loc[keep, g].to_numpy(), y, beta) for g in features}
                for g, score in vals.items():
                    scores.append(dict(definition=definition, target=target, beta=beta, gene=g, best_univariate_fbeta=score))
                ranked[(target, beta)] = sorted(vals, key=lambda g: (-vals[g], g))
        # Round-robin across both target directions and beta values; add 6 novel genes.
        expanded = list(CORE)
        for rank in range(len(features)):
            for ranking in ranked.values():
                if ranking[rank] not in expanded:
                    expanded.append(ranking[rank])
                if len(expanded) == 12:
                    break
            if len(expanded) == 12:
                break
        panels[definition] = expanded
        jobs = make_jobs(expanded, definition, "expanded_" + definition)
        if definition == "crude":
            jobs = jobs.loc[~jobs.features.map(lambda s: set(s.split(";")).issubset(CORE))]
        jobs.to_csv(OUT / f"expanded_{definition}_jobs.tsv", sep="\t", index=False)
    pd.DataFrame(scores).to_csv(OUT / "surface_univariate_screen.tsv", sep="\t", index=False)
    params = dict(starting_n=len(data), crude_counts=data.crude_label.value_counts().to_dict(),
                  module_counts=data.module_label.value_counts().to_dict(),
                  module_rule="Mean gene-wise z-scores; high >= pooled 65th percentile and opposite low <= pooled 35th percentile; middle excluded",
                  module_parameters=module_parameters, positivity="log1p(CP10K) > 0",
                  excluded_clusters=["6", "7"], excluded_target_genes=SST_MODULE + PV_MODULE,
                  erbb4_excluded_from_module="Detected in 100% of entry population; retained as allowed gate feature",
                  surface_source="https://doi.org/10.1371/journal.pone.0121314.s003",
                  cspa_sha256=hashlib.sha256(cspa_path.read_bytes()).hexdigest(),
                  allowed_surface_gene_count=len(features), searched_panels=panels,
                  search="Core panel exhaustive 1-3 gene subsets first; all allowed features univariate screened; exhaustive 1-3 gene combinations of 12-gene panel per definition",
                  betas=[.5, 1, 2], thresholds="Hypergate permits inclusive lower and upper bounds",
                  evaluation="Globally optimized, in-sample descriptive metrics; same thresholds applied per sample; no separate sample fitting",
                  h5ad=str(PROJECT_ROOT / "results/python_anndata/varela_div90.h5ad"),
                  project_root=str(PROJECT_ROOT), output_dir=str(OUT), code_root=str(ROOT))
    (OUT / "run_parameters.json").write_text(json.dumps(params, indent=2))
    print(json.dumps({k: params[k] for k in ["crude_counts", "module_counts", "allowed_surface_gene_count", "searched_panels"]}, indent=2), flush=True)


if __name__ == "__main__":
    main()
