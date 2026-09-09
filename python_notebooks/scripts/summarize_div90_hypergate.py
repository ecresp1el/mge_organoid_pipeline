"""Validate R fits, benchmark threshold sweeps, and render the Hypergate report."""
from pathlib import Path
import hashlib
import itertools
import json
import shutil
from zipfile import ZipFile, ZIP_DEFLATED

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from div90_hypergate_paths import PROJECT_ROOT, REPO_ROOT

ROOT = REPO_ROOT
OUT = PROJECT_ROOT / "results/div90_hypergate_sst_pv"
FIG = OUT / "figures"
CORE = ["ERBB4", "CXCR4", "ACKR3", "PLXNA2", "NRP1", "NRP2"]
COLORS = {"SST-like": "#CB4777", "non-SST/PV-candidate": "#277DA8", "PV-like": "#277DA8", "ambiguous": "#C9CED2"}


def normalized_rules(rules):
    return sorted([dict(gene=r["gene"], op=r["op"].strip(), threshold=float(r["threshold"])) for r in rules],
                  key=lambda r: (r["gene"], r["op"], r["threshold"]))


def rule_id(rules):
    return hashlib.sha256(json.dumps(normalized_rules(rules), sort_keys=True).encode()).hexdigest()[:16]


def rule_text(rules):
    return " AND ".join(f"{r['gene']} {r['op']} {r['threshold']:.6g}" for r in normalized_rules(rules)) or "All entry cells"


def capture(data, rules):
    keep = np.ones(len(data), dtype=bool)
    for r in rules:
        v = data[r["gene"]].to_numpy()
        if r["op"] == ">=": keep &= v >= r["threshold"]
        elif r["op"] == "<=": keep &= v <= r["threshold"]
        elif r["op"] == ">": keep &= v > r["threshold"]
        else: raise ValueError(r)
    return keep


def metric(y, c):
    tp = int(np.sum(y & c)); fp = int(np.sum(~y & c))
    fn = int(np.sum(y & ~c)); tn = int(np.sum(~y & ~c))
    n = len(y); n_target = int(y.sum()); n_capture = int(c.sum())
    purity = tp / n_capture if n_capture else np.nan
    recall = tp / n_target if n_target else np.nan
    prior = n_target / n if n else np.nan
    r = dict(tp=tp, fp=fp, fn=fn, tn=tn, n_evaluated=n, n_target=n_target,
             n_captured=n_capture, precision=purity, purity=purity, recall=recall, yield_=recall,
             target_fraction_before=prior, target_fraction_after=purity,
             fold_enrichment=purity / prior if prior > 0 else np.nan)
    for beta, name in [(.5, "F0.5"), (1, "F1"), (2, "F2")]:
        denom = (1 + beta**2) * tp + beta**2 * fn + fp
        r[name] = (1 + beta**2) * tp / denom if denom else np.nan
    return r


def evaluation_masks(data, definition, target):
    col = "crude_label" if definition == "crude" else "module_label"
    eligible = data[col].ne("ambiguous").to_numpy()
    return eligible, data[col].eq(target).to_numpy()


def threshold_grid(values):
    pos = values[values > 0]
    q = np.quantile(pos, [.1, .25, .5, .75, .9]) if len(pos) else []
    return sorted(set([0., .5, 1., 2., 3., 4., *map(float, q)]))


def benchmarks(data):
    yield "entry_baseline", []
    grids = {g: threshold_grid(data[g].to_numpy()) for g in CORE}
    for g in ["ERBB4", "CXCR4"]:
        for threshold in grids[g]:
            yield g + "_high", [dict(gene=g, op=">", threshold=threshold)]
    for g in CORE[1:]:
        for direction in (["high", "low"] if g != "NRP1" else ["high"]):
            for erbb4_threshold, other_threshold in itertools.product(grids["ERBB4"], grids[g]):
                yield f"ERBB4_high_AND_{g}_{direction}", [
                    dict(gene="ERBB4", op=">", threshold=erbb4_threshold),
                    dict(gene=g, op=">" if direction == "high" else "<=", threshold=other_threshold)]


def savefig(fig, name):
    fig.savefig(FIG / (name + ".png"), dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def style_umap(ax):
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values(): spine.set_visible(False)


def plot_context(data, background):
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    for ax, col, title in zip(axes, ["SST", "sst_module", "pv_module"],
                             ["SST expression", "SST-like module", "PV-like module"]):
        ax.scatter(background.loupe_x, background.loupe_y, s=2, c="#E6E8EA", linewidths=0)
        values = data[col]
        vmin, vmax = (0, float(values.quantile(.99))) if col == "SST" else tuple(values.quantile([.01, .99]))
        order = np.argsort(values)
        pts = ax.scatter(data.loupe_x.iloc[order], data.loupe_y.iloc[order], c=values.iloc[order],
                         cmap="viridis", s=5, linewidths=0, vmin=vmin, vmax=vmax)
        fig.colorbar(pts, ax=ax, shrink=.65, label="log1p(CP10K)" if col == "SST" else "Mean gene z-score")
        ax.set_title(title); style_umap(ax)
    fig.suptitle("4,768 LHX6+/ERBB4+ cells on the recovered cortical recluster UMAP", fontsize=13)
    fig.tight_layout()
    savefig(fig, "sst_and_module_umaps")
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    for ax, col, title in zip(axes, ["crude_label", "module_label"], ["First pass: SST detection", "Second pass: opposing module extremes"]):
        ax.scatter(background.loupe_x, background.loupe_y, s=2, c="#ECEEEF", linewidths=0)
        for label, group in data.groupby(col, sort=False):
            ax.scatter(group.loupe_x, group.loupe_y, s=5, c=COLORS[label], label=f"{label} (n={len(group):,})", linewidths=0)
        ax.set_title(title); style_umap(ax); ax.legend(fontsize=7, frameon=False, loc="lower center", bbox_to_anchor=(.5, -.16))
    fig.tight_layout(); savefig(fig, "target_label_umaps")


def plot_best(data, background, best, registry, metrics):
    for definition in ["crude", "module"]:
        row = best[definition]
        rules = registry[row["gate_id"]]; genes = sorted({r["gene"] for r in rules})
        assert len(genes) == 2
        col = "crude_label" if definition == "crude" else "module_label"
        fig, axes = plt.subplots(1, 2, figsize=(10, 4.3))
        ax = axes[0]
        for label, frame in data.groupby(col):
            ax.scatter(frame[genes[0]], frame[genes[1]], c=COLORS[label], s=7, alpha=.45,
                       linewidths=0, label=f"{label} (n={len(frame):,})")
        for r in rules:
            (ax.axvline if r["gene"] == genes[0] else ax.axhline)(r["threshold"], ls="--", c="#171717", lw=1.2)
        bounds = {g: [min(-.1, data[g].min()-.1), float(data[g].max())+.1] for g in genes}
        for r in rules:
            if r["op"] in [">", ">="]: bounds[r["gene"]][0] = max(bounds[r["gene"]][0], r["threshold"])
            else: bounds[r["gene"]][1] = min(bounds[r["gene"]][1], r["threshold"])
        ax.add_patch(Rectangle((bounds[genes[0]][0], bounds[genes[1]][0]),
                              bounds[genes[0]][1]-bounds[genes[0]][0], bounds[genes[1]][1]-bounds[genes[1]][0],
                              facecolor="#66BB66", alpha=.10, edgecolor="#187A28", linewidth=2))
        ax.set_xlabel(genes[0] + " log1p(CP10K)"); ax.set_ylabel(genes[1] + " log1p(CP10K)")
        ax.set_title("Virtual FACS: shaded region captured")
        ax.legend(fontsize=6.5, frameon=False)
        c = capture(data, rules)
        axes[1].scatter(background.loupe_x, background.loupe_y, s=3, c="#E2E5E7", linewidths=0)
        axes[1].scatter(data.loupe_x[c], data.loupe_y[c], s=7, c="#277DA8", linewidths=0)
        style_umap(axes[1]); axes[1].set_title(f"Same gate on recluster UMAP\n{c.sum():,} entry cells captured")
        m = metrics[(definition, row["target"], row["gate_id"])]
        fig.suptitle(f"{definition.capitalize()}: {rule_text(rules)}\nPurity {m['purity']:.1%} | yield {m['recall']:.1%}", fontsize=11)
        fig.tight_layout(); savefig(fig, f"best_pair_{definition}")


def main():
    FIG.mkdir(exist_ok=True)
    data = pd.read_csv(OUT / "hypergate_input.tsv.gz", sep="\t")
    labels = pd.read_csv(OUT / "target_labels.tsv.gz", sep="\t")
    assert data.cell_id.tolist() == labels.cell_id.tolist()
    allowed = set(pd.read_csv(OUT / "allowed_surface_markers.tsv", sep="\t").gene)
    params = json.loads((OUT / "run_parameters.json").read_text())
    paths = [OUT / "core_fits.tsv", *sorted(OUT.glob("expanded_*_fits_*.tsv"))]
    assert len(paths) == 9, paths
    fits = pd.concat([pd.read_csv(p, sep="\t", escapechar="\\") for p in paths], ignore_index=True)
    expected = sum(len(pd.read_csv(OUT / name, sep="\t")) for name in ["core_jobs.tsv", "expanded_crude_jobs.tsv", "expanded_module_jobs.tsv"])
    assert len(fits) == expected, (len(fits), expected)
    assert fits.job_id.is_unique
    errors = fits.loc[fits.error.notna() & fits.error.ne("")]
    errors.to_csv(OUT / "provenance/fit_errors.tsv", sep="\t", index=False)
    if len(errors): raise RuntimeError(f"{len(errors)} failed fits must be resolved")
    registry = {}; captures = {}; all_rows = []; computed = {}; sample_rows = []
    masks = {(definition, target): evaluation_masks(data, definition, target)
             for definition, targets in [("crude", ["non-SST/PV-candidate", "SST-like"]), ("module", ["PV-like", "SST-like"])] for target in targets}
    sample_masks = {sample: data["sample"].eq(sample).to_numpy() for sample in sorted(data["sample"].unique())}

    def register(rules):
        rules = normalized_rules(rules)
        assert all(r["gene"] in allowed for r in rules)
        gid = rule_id(rules)
        if gid not in registry:
            registry[gid] = rules; captures[gid] = capture(data, rules)
        return gid

    def evaluate(definition, target, gid):
        key = (definition, target, gid)
        if key in computed: return computed[key]
        eligible, y = masks[(definition, target)]; c = captures[gid]
        r = metric(y[eligible], c[eligible])
        r.update(n_captured_all_entry=int(c.sum()), n_ambiguous_captured=int((c & ~eligible).sum()),
                 definite_target_fraction_all_captured=int((c & y).sum()) / c.sum() if c.sum() else np.nan)
        computed[key] = r
        for sample, sm in sample_masks.items():
            k = eligible & sm
            sample_rows.append(dict(definition=definition, target=target, gate_id=gid, sample=sample,
                                    gate=rule_text(registry[gid]), **metric(y[k], c[k]),
                                    n_captured_all_entry=int((c & sm).sum()), n_ambiguous_captured=int((c & sm & ~eligible).sum())))
        return r

    for row in fits.itertuples(index=False):
        rules = json.loads(row.rules); gid = register(rules)
        m = evaluate(row.definition, row.target, gid)
        # Independently reproduce the R package's confusion counts at full precision.
        assert all(m[name] == int(getattr(row, name)) for name in ["tp", "fp", "fn", "tn"]), row.job_id
        all_rows.append(dict(job_id=row.job_id, method="R_hypergate", stage=row.stage,
                             definition=row.definition, target=row.target, beta=row.beta,
                             searched_features=row.features, requested_markers=row.requested_markers,
                             actual_markers=len({r["gene"] for r in rules}), n_bounds=len(rules),
                             surface_markers_including_entry=len({r["gene"] for r in rules} | {"ERBB4"}),
                             gate_id=gid, gate=rule_text(rules), rules=json.dumps(normalized_rules(rules)), **m))
    for index, (hypothesis, rules) in enumerate(benchmarks(data)):
        gid = register(rules)
        for definition, target in masks:
            all_rows.append(dict(job_id=f"benchmark_{index}_{definition}_{target}", method="threshold_benchmark", stage=hypothesis,
                                 definition=definition, target=target, beta=np.nan,
                                 searched_features=";".join(sorted({r["gene"] for r in rules})), requested_markers=len(rules),
                                 actual_markers=len({r["gene"] for r in rules}), n_bounds=len(rules),
                                 surface_markers_including_entry=len({r["gene"] for r in rules} | {"ERBB4"}),
                                 gate_id=gid, gate=rule_text(rules), rules=json.dumps(rules), **evaluate(definition, target, gid)))
    results = pd.DataFrame(all_rows).rename(columns={"yield_": "yield"})
    results.to_csv(OUT / "hypergate_results_all.tsv", sep="\t", index=False)
    candidates = results.loc[results.actual_markers.between(1,3)].sort_values(
        ["F1", "actual_markers", "n_bounds", "gate_id"], ascending=[False, True, True, True])
    unique = candidates.drop_duplicates(["definition", "target", "gate_id"])
    top = unique.groupby(["definition", "target"], sort=False).head(20)
    top.to_csv(OUT / "top_gates.tsv", sep="\t", index=False)
    best = {}
    for definition, target in [("crude", "non-SST/PV-candidate"), ("module", "PV-like")]:
        subset = unique.loc[unique.definition.eq(definition) & unique.target.eq(target) & unique.actual_markers.eq(2)]
        best[definition] = subset.iloc[0].to_dict()
    # Choose one common compact rule by the lower F1 across the two target definitions.
    fit_gids = results.loc[results.method.eq("R_hypergate") & results.actual_markers.eq(2), "gate_id"].unique()
    shared_rows = []
    for gid in fit_gids:
        crude = evaluate("crude", "non-SST/PV-candidate", gid)
        module = evaluate("module", "PV-like", gid)
        shared_rows.append(dict(gate_id=gid, gate=rule_text(registry[gid]),
                                min_F1=min(crude["F1"], module["F1"]), mean_F1=(crude["F1"]+module["F1"])/2,
                                crude_purity=crude["purity"], crude_yield=crude["recall"], crude_F1=crude["F1"],
                                crude_enrichment=crude["fold_enrichment"], module_purity=module["purity"], module_yield=module["recall"],
                                module_F1=module["F1"], module_enrichment=module["fold_enrichment"],
                                n_captured_all_entry=int(captures[gid].sum())))
    shared = pd.DataFrame(shared_rows).sort_values(["min_F1", "mean_F1", "gate_id"], ascending=[False,False,True])
    shared.to_csv(OUT / "shared_pair_comparison.tsv", sep="\t", index=False)
    shared_best = shared.iloc[0].to_dict()
    pd.DataFrame(sample_rows).rename(columns={"yield_": "yield"}).to_csv(OUT / "per_sample_validation.tsv", sep="\t", index=False)
    selected_gates = {"best_shared_pair": shared_best["gate_id"], **{f"best_pair_{k}": v["gate_id"] for k,v in best.items()}}
    cells = labels.copy()
    for name, gid in selected_gates.items(): cells[name + "_captured"] = captures[gid]
    cells.to_csv(OUT / "best_gate_cells.tsv.gz", sep="\t", index=False)
    (OUT / "gate_rules.json").write_text(json.dumps(registry, indent=2))
    selected_rules = {name: dict(gate_id=gid, rules=registry[gid], gate=rule_text(registry[gid])) for name,gid in selected_gates.items()}
    (OUT / "best_gate_rules.json").write_text(json.dumps(selected_rules, indent=2))
    bench = results.loc[results.method.eq("threshold_benchmark")]
    benchmark_best = bench.sort_values("F1", ascending=False).drop_duplicates(["definition", "target", "stage"])
    benchmark_best.to_csv(OUT / "benchmark_best_by_hypothesis.tsv", sep="\t", index=False)
    background = pd.read_csv(PROJECT_ROOT / "final_figures/div90_guidance_with_sst/tables/joined_with_sst.tsv.gz", sep="\t")
    background = background.loc[background.set_id.eq("cortical_only")]
    if not all((FIG / name).exists() for name in ["sst_and_module_umaps.png", "target_label_umaps.png"]):
        plot_context(data, background)
    plot_best(data, background, best, registry, computed)
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    for ax, ((definition, target), group) in zip(axes.flat, results.groupby(["definition", "target"])):
        for method, g in group.groupby("method"):
            ax.scatter(g.recall, g.purity, s=8, alpha=.20, label=method)
        prior = float(group.target_fraction_before.iloc[0]); ax.axhline(prior, ls="--", color="#777777", lw=1)
        ax.set(xlim=(0,1.02), ylim=(0,1.02), xlabel="Yield / recall", ylabel="Purity / precision", title=f"{definition}: {target}")
        ax.legend(fontsize=7)
    fig.tight_layout(); savefig(fig, "precision_yield_all_gates")
    fig, ax = plt.subplots(figsize=(12, 7))
    rows = top.loc[top.definition.eq("crude") & top.target.eq("non-SST/PV-candidate")].head(20)
    show = [[r.gate, f"{r.purity:.1%}", f"{r.recall:.1%}", f"{r.F1:.3f}", f"{r.fold_enrichment:.2f}"] for r in rows.itertuples()]
    ax.axis("off"); table = ax.table(cellText=show, colLabels=["Gate (log1p CP10K)", "Purity", "Yield", "F1", "Enrichment"],
                                     colWidths=[.66,.085,.085,.085,.085], cellLoc="left", loc="center")
    table.auto_set_font_size(False); table.set_fontsize(7); table.scale(1,1.5)
    ax.set_title("Top 20 first-pass gates for non-SST/PV-candidate cells", pad=15)
    savefig(fig, "top_20_crude_gates")
    sample = pd.DataFrame(sample_rows)
    sg = sample.loc[sample.gate_id.eq(shared_best["gate_id"]) & ((sample.definition.eq("crude") & sample.target.eq("non-SST/PV-candidate")) |
                      (sample.definition.eq("module") & sample.target.eq("PV-like")))]
    fig, axes = plt.subplots(1,2,figsize=(11,4))
    for ax, definition in zip(axes,["crude","module"]):
        frame = sg.loc[sg.definition.eq(definition)].sort_values("sample")
        x = np.arange(len(frame)); ax.bar(x-.18,frame.purity,.36,label="Purity",color="#277DA8")
        ax.bar(x+.18,frame.recall,.36,label="Yield",color="#C47F36")
        ax.set_xticks(x, [f"{r.sample}\nn={r.n_evaluated}" for r in frame.itertuples()],rotation=35,ha="right",fontsize=7)
        ax.set(ylim=(0,1.05),title=definition+": globally learned shared gate");ax.legend(fontsize=8)
    fig.tight_layout();savefig(fig,"shared_gate_per_sample")
    params.update(hypergate_package_version="0.8.5", n_hypergate_fits=len(fits), n_benchmark_evaluations=len(bench),
                  n_unique_gate_rules=len(registry), confusion_counts_verified_against_R=True,
                  shared_gate_selection="2 distinct markers; maximize minimum F1 across crude non-SST and module PV-like, then mean F1",
                  benchmark_threshold_grid="0,0.5,1,2,3,4 plus 10/25/50/75/90th percentiles of positive expression; high > threshold, low <= threshold",
                  module_metric_denominator="Only 1253 unambiguous labeled cells; ambiguous captures and full-population capture counts reported separately")
    (OUT / "run_parameters.json").write_text(json.dumps(params, indent=2))
    report = ["# DIV90 Hypergate: SST versus non-SST first, module extremes second", "",
              f"Actual R hypergate 0.8.5: {len(fits):,} fits; {len(bench):,} benchmark evaluations. No trajectory reconstruction or UMAP recomputation.", "",
              "## Best common two-marker RNA gate", "", f"**{shared_best['gate']}**", "",
              "All thresholds are log1p(CP10K) RNA expression. LHX6>0 and ERBB4>0 are pre-existing entry criteria.", "",
              "| Target definition | Purity | Yield | F1 | Fold-enrichment |", "|---|---:|---:|---:|---:|",
              f"| SST-negative / non-SST/PV-candidate | {shared_best['crude_purity']:.1%} | {shared_best['crude_yield']:.1%} | {shared_best['crude_F1']:.3f} | {shared_best['crude_enrichment']:.2f} |",
              f"| Module-defined PV-like | {shared_best['module_purity']:.1%} | {shared_best['module_yield']:.1%} | {shared_best['module_F1']:.3f} | {shared_best['module_enrichment']:.2f} |", "",
              "The same exact thresholds are used in both rows. The first row uses all 4,768 entry cells (1,849 targets); the second uses 1,253 unambiguous cells (638 PV-like, 615 SST-like). The remaining 3,515 cells are not labeled as negatives in the module analysis.", "",
              "## ERBB4 + CXCR4 benchmark", ""]
    positive_pair_comparisons = []
    for definition, target in [("crude", "non-SST/PV-candidate"), ("module", "PV-like")]:
        row = benchmark_best.loc[benchmark_best.definition.eq(definition) & benchmark_best.target.eq(target) &
                                 benchmark_best.stage.eq("ERBB4_high_AND_CXCR4_high")].iloc[0]
        b = best[definition]
        positive_pair_comparisons.append(dict(definition=definition, best_pair_F1=b["F1"],
                                              erbb4_cxcr4_positive_F1=row.F1,
                                              difference=b["F1"]-row.F1))
        report += [f"- {definition}: best positive ERBB4/CXCR4 sweep = `{row.gate}`; purity {row.purity:.1%}, yield {row.recall:.1%}, F1 {row.F1:.3f}. Best two-marker gate by F1 = `{b['gate']}`; purity {b['purity']:.1%}, yield {b['recall']:.1%}, F1 {b['F1']:.3f}."]
    verdict = "outperformed by another pair" if all(r["difference"] > .02 for r in positive_pair_comparisons) else "partially supported"
    report += ["", f"**ERBB4 + CXCR4: {verdict}.** Here, 'outperformed' requires an absolute F1 advantage greater than 0.02 in both definitions versus the best positive-positive threshold sweep.", "",
               "## Same-gate sample performance", "",
               "| Definition | Sample | Labeled cells | Purity | Yield |", "|---|---|---:|---:|---:|"]
    for r in sg.itertuples():
        yield_text = f"{r.recall:.1%}" if np.isfinite(r.recall) else "N/A (no target cells)"
        purity_text = f"{r.purity:.1%}" if np.isfinite(r.purity) else "N/A (no captured cells)"
        report.append(f"| {r.definition} | {r.sample} | {r.n_evaluated} | {purity_text} | {yield_text} |")
    report += ["", "Purity varies substantially across samples; the pooled result is modest enrichment rather than a clean SST/non-SST separator. The different starting target fractions are retained in per_sample_validation.tsv. Gates relying on undetected transcripts require particular care when translating to protein measurements."]
    decision = dict(shared_best=shared_best, erbb4_cxcr4_verdict=verdict,
                    positive_pair_comparisons=positive_pair_comparisons,
                    best_pairs={k:{field:v[field] for field in ["gate_id","gate","purity","recall","F1","fold_enrichment"]} for k,v in best.items()})
    (OUT / "decision_summary.json").write_text(json.dumps(decision, indent=2))
    report += ["", "## Labels and surface feature selection", "",
               "First pass: SST>0 versus SST==0, both target directions, with the requested six surface markers. The core run completed before the expanded and module Hypergate searches.", "",
               "SST module: SST, NR2F2, GRIK1. PV-like module: MEF2C, MAF, MAFB, KCNC1, KCNC2, GAD1, GAD2. Scores are means of gene-wise z-scores within the 4,768 entry cells. High is >=65th percentile and the opposing score must be <=35th percentile. Exact centering, scaling and thresholds are saved in run_parameters.json. ERBB4 is excluded from module scoring because every entry cell expresses it.", "",
               f"{params['allowed_surface_gene_count']} expressed surface features passed the CSPA high-confidence plus transmembrane/GPI annotation and expression filters, together with the user-specified core panel. Target-module genes, obvious housekeeping genes, and BCAN/FRAS1 extracellular-matrix or surface-isoform ambiguities are excluded. All allowed genes were screened univariately; combinations were exhaustively fit only within the documented 12-gene panel per definition, not across all surfaceome triplets.", "",
               "[Surfaceome annotation source](https://doi.org/10.1371/journal.pone.0121314.s003). [R Hypergate source](https://github.com/ebecht/hypergate).", "",
               "## Interpretation and files", "",
               "These are globally optimized, in-sample RNA gates. Per-sample results use these same gates and are descriptive robustness checks, not held-out validation. RNA zero can reflect non-detection; non-SST/PV-candidate is not established PV fate. Surface membership does not establish an extracellular antibody's availability or protein sorting performance. No antibody reagent has been validated by this analysis.", "",
               "- hypergate_results_all.tsv: every Hypergate fit and benchmark, exact rules, confusion counts and all requested metrics.",
               "- top_gates.tsv: top 20 distinct rules per definition and target direction, ranked by F1.",
               "- per_sample_validation.tsv: per-sample performance of fixed global rules, including low-count samples.",
               "- shared_pair_comparison.tsv: same two-marker rules compared across both definitions.",
               "- target_labels.tsv.gz, allowed_surface_markers.tsv, hypergate_input.tsv.gz: reproducible aligned inputs.",
               "- best_gate_cells.tsv.gz and best_gate_rules.json: captured cell IDs and exact rules.",
               "- figures/: virtual FACS, original-layout overlays, expression/module UMAPs, precision-yield and ranking plots.",
               "", "## Rerun", "", "From the repository root:", "", "```bash",
               "/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/python python_notebooks/scripts/run_div90_hypergate_workflow.py",
               "```", "", "The workflow reuses completed fit checkpoints and refuses to silently overwrite incomplete ones."]
    (OUT / "REPORT.md").write_text("\n".join(report)+"\n")
    provenance = dict(n_fits_checked=len(fits), all_R_confusion_counts_reproduced=True,
                      no_target_module_genes_in_gate_features=not(set(params["excluded_target_genes"]) & allowed),
                      n_cells=len(data), n_SST_positive=int(data.SST.gt(0).sum()),
                      source_coordinates="Recovered Loupe cortical local UMAP, unchanged")
    (OUT / "provenance/validation.json").write_text(json.dumps(provenance, indent=2))
    handoff = ROOT / "python_notebooks/HANDOFF_div90_hypergate_sst_pv.md"
    if handoff.exists():
        shutil.copy2(handoff, OUT / "PHASE2_HANDOFF.md")
        with (OUT / "REPORT.md").open("a") as handle:
            handle.write("\n## Next-chat handoff\n\nRead [PHASE2_HANDOFF.md](PHASE2_HANDOFF.md). Runtime outputs live in PROJECT_ROOT, outside the Git checkout.\n")
    for name in ["div90_hypergate_paths.py", "run_div90_hypergate.R", "prepare_div90_hypergate.py", "run_div90_hypergate_expanded.py", "summarize_div90_hypergate.py", "run_div90_hypergate_workflow.py"]:
        shutil.copy2(ROOT / "python_notebooks/scripts" / name, OUT / "provenance" / name)
    archive = PROJECT_ROOT / "results/div90_hypergate_sst_pv_results.zip"
    with ZipFile(archive, "w", ZIP_DEFLATED) as z:
        for path in sorted(OUT.rglob("*")):
            if path.is_file(): z.write(path, str(path.relative_to(OUT)))
    with ZipFile(archive) as z: assert z.testzip() is None
    print(json.dumps(shared_best, indent=2), flush=True)
    print("Wrote", archive, flush=True)


if __name__ == "__main__":
    main()
