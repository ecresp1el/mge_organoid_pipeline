#!/usr/bin/env python3
"""Validate marker-supported but confounded DIV90 subcluster splits."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import statsmodels.formula.api as smf
from scipy import sparse
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import adjusted_rand_score


PROJECT_ROOT_DEFAULT = "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder"

CLUSTER1_A = ["EDIL3", "VWC2", "LRFN5", "ITGB8", "CRABP1", "NRCAM", "PCDH10", "CDH4", "RYR3"]
CLUSTER1_B = ["ACKR3", "ZEB2", "BTG1", "GPR173", "CSAD", "WLS"]

CLUSTER3_A = ["L3MBTL1", "PEG3", "SYT4", "LRFN5", "RUNX1T1", "MYO16", "MAGEH1", "TMEFF2", "PLXNA4", "TENT5A"]
CLUSTER3_B = ["BEX3", "USP11", "TMSB4X", "NAP1L3", "CRABP1", "MEST", "DCX", "TUBB2B", "NNAT"]

CLUSTER11_A = ["CCSER1", "ZEB2", "ACKR3", "SLITRK2", "MAF", "RASA1", "PRKG1", "HCN1", "WSB1", "PTCHD4"]
CLUSTER11_B = ["CPE", "TUBB2A", "PEG10", "CNTN1", "PCDH10", "DKK3", "MEST", "ATP1B1"]

CLUSTER3_NEIGHBOR_MARKERS = [
    "BEX3", "CRABP1", "DCX", "TUBB2B", "NNAT", "MEST",
    "LHX6", "SOX6", "ERBB4", "MAF", "GAD1", "GAD2",
]

CLUSTER11_CLUSTER3_MARKERS = [
    "DCX", "STMN2", "DCLK1", "TUBB2A", "TUBB2B", "NNAT", "CPE", "CNTN1",
    "ACKR3", "ZEB2", "MAF", "LHX6", "SOX6", "ERBB4", "GAD1", "GAD2", "KCNC1", "KCNC2",
]

STRESS_GENES = [
    "FOS", "JUN", "JUNB", "JUND", "ATF3", "DDIT3", "HSPA1A", "HSPA1B",
    "HSPH1", "HSP90AA1", "HSP90AB1", "DNAJB1", "HMOX1", "TXNIP",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", default=os.environ.get("PROJECT_ROOT", PROJECT_ROOT_DEFAULT))
    parser.add_argument("--h5ad", default=None)
    parser.add_argument("--sample-map", default="metadata/div30_div90_sample_id_to_biolabel_map.tsv")
    parser.add_argument(
        "--outdir",
        default=None,
        help="Default: PROJECT_ROOT/results/div90_failed_interesting_subcluster_validation/div90_failed_interesting_subcluster_validation_v1",
    )
    return parser.parse_args()


def ensure_dirs(outdir: Path) -> tuple[Path, Path]:
    plot_dir = outdir / "plots"
    table_dir = outdir / "tables"
    plot_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)
    return plot_dir, table_dir


def add_cell_line(adata: ad.AnnData, sample_map_path: Path) -> None:
    sample_map = pd.read_csv(sample_map_path, sep="\t")
    div90 = sample_map.loc[sample_map["DIV"].astype(str).str.upper() == "DIV90"].copy()
    lookup = div90.set_index("run_sample_id")["biological_label"].astype(str).to_dict()
    labels = adata.obs["orig.ident"].astype(str).map(lookup).fillna("unknown")
    adata.obs["cell_line"] = labels.str.replace(r"_rep\d+$", "", regex=True)


def numeric_obs(adata: ad.AnnData) -> None:
    for col in ["nCount_RNA", "nFeature_RNA", "percent.mt", "S.Score", "G2M.Score", "CC.Difference"]:
        if col in adata.obs:
            adata.obs[col] = pd.to_numeric(adata.obs[col].astype(str), errors="coerce")


def matrix_values(adata: ad.AnnData, genes: list[str]) -> tuple[np.ndarray, list[str]]:
    present = [g for g in genes if g in adata.var_names]
    if not present:
        return np.full((adata.n_obs, 0), np.nan), []
    mat = adata[:, present].X
    if sparse.issparse(mat):
        mat = mat.toarray()
    return np.asarray(mat, dtype=float), present


def add_module_score(adata: ad.AnnData, genes: list[str], name: str) -> list[str]:
    values, present = matrix_values(adata, genes)
    adata.obs[name] = np.nan if values.shape[1] == 0 else np.nanmean(values, axis=1)
    return present


def local_split(parent: ad.AnnData, resolution: float, key: str) -> ad.AnnData:
    sub = parent.copy()
    sub.layers["log_expr"] = sub.X.copy()
    sc.pp.highly_variable_genes(sub, n_top_genes=min(2000, sub.n_vars), flavor="seurat")
    sub = sub[:, sub.var["highly_variable"].to_numpy()].copy()
    sc.pp.scale(sub, max_value=10)
    sc.tl.pca(sub, n_comps=min(30, sub.n_obs - 2, sub.n_vars - 1), svd_solver="arpack")
    sc.pp.neighbors(sub, n_neighbors=min(30, sub.n_obs - 1), n_pcs=min(30, sub.obsm["X_pca"].shape[1]))
    sc.tl.umap(sub, random_state=0)
    sc.tl.leiden(sub, resolution=resolution, key_added=key, random_state=0)
    return sub


def restore_scores(parent_full: ad.AnnData, split_hvg: ad.AnnData, score_cols: list[str]) -> ad.AnnData:
    out = split_hvg.copy()
    for col in score_cols:
        out.obs[col] = parent_full.obs.loc[out.obs_names, col].to_numpy()
    for col in ["orig.ident", "cell_line", "nCount_RNA", "nFeature_RNA", "percent.mt", "S.Score", "G2M.Score", "Phase", "cluster_number_name"]:
        if col in parent_full.obs:
            out.obs[col] = parent_full.obs.loc[out.obs_names, col].to_numpy()
    return out


def save_scanpy_umap(adata: ad.AnnData, colors: list[str], path: Path, title: str) -> None:
    fig = sc.pl.umap(adata, color=colors, ncols=4, show=False, frameon=False, size=11, return_fig=True)
    fig.suptitle(title, y=1.01)
    fig.savefig(path, dpi=240, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def composition_table(adata: ad.AnnData, group_col: str, split_col: str, score_cols: list[str], n_label: str) -> pd.DataFrame:
    rows = []
    for group, df in adata.obs.groupby(group_col, observed=False):
        counts = df[split_col].astype(str).value_counts(normalize=False)
        row = {group_col: group, n_label: int(df.shape[0])}
        for sub in sorted(adata.obs[split_col].astype(str).unique()):
            row[f"fraction_subcluster_{sub}"] = float(counts.get(sub, 0) / df.shape[0])
        for score in score_cols:
            row[f"mean_{score}"] = float(pd.to_numeric(df[score], errors="coerce").mean())
        rows.append(row)
    return pd.DataFrame(rows).sort_values(group_col)


def boxplot_scores(adata: ad.AnnData, group_col: str, score_cols: list[str], path: Path, title: str) -> None:
    df = adata.obs[[group_col] + score_cols].copy()
    long = df.melt(id_vars=group_col, value_vars=score_cols, var_name="module", value_name="score")
    fig, ax = plt.subplots(figsize=(max(8, 0.45 * long[group_col].nunique() * len(score_cols)), 5))
    sns.boxplot(data=long, x=group_col, y="score", hue="module", ax=ax, fliersize=1)
    ax.set_title(title)
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    fig.savefig(path, dpi=240)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def sample_correlation_table(adata: ad.AnnData, score_a: str, score_b: str, min_cells: int, split_col: str) -> pd.DataFrame:
    rows = []
    coords = pd.DataFrame(adata.obsm["X_umap"], index=adata.obs_names, columns=["umap_1", "umap_2"])
    obs = adata.obs.join(coords)
    for sample, df in obs.groupby("orig.ident", observed=False):
        if df.shape[0] < min_cells:
            continue
        a = pd.to_numeric(df[score_a], errors="coerce")
        b = pd.to_numeric(df[score_b], errors="coerce")
        finite = a.notna() & b.notna()
        pearson = pearsonr(a[finite], b[finite]).statistic if finite.sum() >= 3 else np.nan
        spearman = spearmanr(a[finite], b[finite]).statistic if finite.sum() >= 3 else np.nan
        high_a = df.loc[a >= a.quantile(0.75), ["umap_1", "umap_2"]]
        high_b = df.loc[b >= b.quantile(0.75), ["umap_1", "umap_2"]]
        dist = np.linalg.norm(high_a.mean().to_numpy() - high_b.mean().to_numpy()) if not high_a.empty and not high_b.empty else np.nan
        rows.append({
            "orig.ident": sample,
            "n_cells": int(df.shape[0]),
            "pearson_module_correlation": pearson,
            "spearman_module_correlation": spearman,
            "umap_centroid_distance_top_quartiles": dist,
            "fraction_subcluster_0": float((df[split_col].astype(str) == "0").mean()),
            "fraction_subcluster_1": float((df[split_col].astype(str) == "1").mean()),
        })
    return pd.DataFrame(rows)


def dotplot_from_obs(adata: ad.AnnData, group_col: str, genes: list[str], path: Path, title: str) -> pd.DataFrame:
    rows = []
    values, present = matrix_values(adata, genes)
    expr = pd.DataFrame(values, index=adata.obs_names, columns=present)
    expr[group_col] = adata.obs[group_col].astype(str).to_numpy()
    for group, df in expr.groupby(group_col):
        for gene in present:
            vals = pd.to_numeric(df[gene], errors="coerce").to_numpy(float)
            rows.append({
                group_col: group,
                "gene": gene,
                "mean_expr": float(np.nanmean(vals)),
                "pct_expressed": float(np.mean(vals > 0) * 100),
            })
    summary = pd.DataFrame(rows)
    groups = sorted(summary[group_col].unique())
    fig, ax = plt.subplots(figsize=(max(8, 0.45 * len(present) + 2), max(4, 0.45 * len(groups) + 1.5)))
    x = {g: i for i, g in enumerate(present)}
    y = {g: i for i, g in enumerate(groups)}
    ax.scatter(
        [x[g] for g in summary["gene"]],
        [y[g] for g in summary[group_col]],
        s=np.maximum(8, summary["pct_expressed"].to_numpy(float) * 3.2),
        c=summary["mean_expr"].to_numpy(float),
        cmap="viridis",
        edgecolors="black",
        linewidths=0.2,
    )
    ax.set_xticks(range(len(present)))
    ax.set_xticklabels(present, rotation=45, ha="right")
    ax.set_yticks(range(len(groups)))
    ax.set_yticklabels(groups)
    ax.invert_yaxis()
    ax.set_title(title)
    fig.colorbar(ax.collections[0], ax=ax, label="Mean expression")
    fig.tight_layout()
    fig.savefig(path, dpi=240)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)
    return summary


def add_qc_modules(adata: ad.AnnData) -> tuple[list[str], list[str]]:
    ribosomal = [g for g in adata.var_names if g.startswith("RPL") or g.startswith("RPS")]
    stress_present = [g for g in STRESS_GENES if g in adata.var_names]
    add_module_score(adata, ribosomal, "ribosomal_module")
    add_module_score(adata, stress_present, "stress_module")
    return ribosomal, stress_present


def cluster3_qc_plot(adata: ad.AnnData, split_col: str, path: Path) -> None:
    cols = ["percent.mt", "ribosomal_module", "stress_module", "nCount_RNA", "nFeature_RNA", "S.Score", "G2M.Score"]
    long = adata.obs[[split_col] + cols].melt(id_vars=split_col, value_vars=cols, var_name="metric", value_name="value")
    long["value"] = pd.to_numeric(long["value"], errors="coerce")
    fig, axes = plt.subplots(2, 4, figsize=(14, 7), squeeze=False)
    for ax, metric in zip(axes.ravel(), cols):
        plot_df = long.loc[long["metric"] == metric].dropna(subset=["value"])
        if plot_df.empty:
            ax.text(0.5, 0.5, "No finite values", ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off()
        else:
            sns.boxplot(data=plot_df, x=split_col, y="value", ax=ax, fliersize=1)
        ax.set_title(metric)
    axes.ravel()[-1].axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=240)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def tertile_table(adata: ad.AnnData, covariate: str, split_col: str, score_cols: list[str]) -> pd.DataFrame:
    obs = adata.obs.copy()
    obs[f"{covariate}_tertile"] = pd.qcut(pd.to_numeric(obs[covariate], errors="coerce"), q=3, labels=["low", "middle", "high"], duplicates="drop")
    rows = []
    for tertile, df in obs.groupby(f"{covariate}_tertile", observed=False):
        row = {f"{covariate}_tertile": str(tertile), "n_cells": int(df.shape[0])}
        for sub in sorted(obs[split_col].astype(str).unique()):
            row[f"fraction_subcluster_{sub}"] = float((df[split_col].astype(str) == sub).mean())
        for score in score_cols:
            row[f"mean_{score}"] = float(pd.to_numeric(df[score], errors="coerce").mean())
        rows.append(row)
    return pd.DataFrame(rows)


def tertile_plot(adata: ad.AnnData, covariate: str, split_col: str, score_cols: list[str], path: Path) -> pd.DataFrame:
    obs = adata.obs.copy()
    tertile_col = f"{covariate}_tertile"
    obs[tertile_col] = pd.qcut(pd.to_numeric(obs[covariate], errors="coerce"), q=3, labels=["low", "middle", "high"], duplicates="drop")
    long = obs[[tertile_col, split_col] + score_cols].melt(id_vars=[tertile_col, split_col], value_vars=score_cols, var_name="module", value_name="score")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    comp = pd.crosstab(obs[tertile_col], obs[split_col].astype(str), normalize="index")
    comp.plot(kind="bar", stacked=True, ax=axes[0])
    axes[0].set_title(f"{split_col} composition by {covariate} tertile")
    axes[0].set_ylabel("fraction")
    sns.boxplot(data=long, x=tertile_col, y="score", hue="module", ax=axes[1], fliersize=1)
    axes[1].set_title(f"Module scores by {covariate} tertile")
    fig.tight_layout()
    fig.savefig(path, dpi=240)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)
    return tertile_table(adata, covariate, split_col, score_cols)


def model_cluster11(adata: ad.AnnData, split_col: str, score_cols: list[str]) -> pd.DataFrame:
    df = adata.obs[[split_col, "nFeature_RNA", "nCount_RNA", "percent.mt", "orig.ident", "cell_line"] + score_cols].copy()
    df = df.rename(columns={split_col: "audit_subcluster", "percent.mt": "percent_mt", "orig.ident": "orig_ident"})
    for col in ["nFeature_RNA", "nCount_RNA", "percent_mt"] + score_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    rows = []
    for score in score_cols:
        formula = f"{score} ~ C(audit_subcluster) + nFeature_RNA + nCount_RNA + percent_mt + C(orig_ident) + C(cell_line)"
        model = smf.ols(formula, data=df).fit()
        for term in model.params.index:
            rows.append({
                "module_score": score,
                "term": term,
                "coef": model.params[term],
                "pvalue": model.pvalues[term],
                "rsquared": model.rsquared,
                "nobs": model.nobs,
            })
    return pd.DataFrame(rows)


def write_interpretation(outdir: Path, lines: list[str]) -> None:
    (outdir / "div90_failed_interesting_subcluster_validation_interpretation.md").write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    project_root = Path(args.project_root).expanduser().resolve()
    h5ad = Path(args.h5ad).expanduser().resolve() if args.h5ad else project_root / "results/python_anndata/varela_div90.h5ad"
    outdir = Path(args.outdir).expanduser().resolve() if args.outdir else project_root / "results/div90_failed_interesting_subcluster_validation/div90_failed_interesting_subcluster_validation_v1"
    plot_dir, table_dir = ensure_dirs(outdir)
    sample_map = Path(args.sample_map)
    if not sample_map.is_absolute():
        sample_map = Path.cwd() / sample_map

    adata = sc.read_h5ad(h5ad)
    numeric_obs(adata)
    add_cell_line(adata, sample_map)
    add_qc_modules(adata)
    adata.obs["cluster_id"] = adata.obs["cluster_id"].astype(str)

    interp = [
        "# DIV90 Failed Interesting Subcluster Validation",
        "",
        "This run validates marker-supported splits in clusters 1, 3, and 11.",
        "It recomputes local Leiden splits, module scores, composition tables, and QC/covariate summaries.",
        "",
    ]

    # Cluster 1
    c1_full = adata[adata.obs["cluster_id"] == "1"].copy()
    c1_a_present = add_module_score(c1_full, CLUSTER1_A, "cluster1_EDIL3_CRABP1_axis")
    c1_b_present = add_module_score(c1_full, CLUSTER1_B, "cluster1_ACKR3_ZEB2_axis")
    c1 = restore_scores(c1_full, local_split(c1_full, 0.2, "cluster1_audit_subcluster"), ["cluster1_EDIL3_CRABP1_axis", "cluster1_ACKR3_ZEB2_axis"])
    save_scanpy_umap(c1, ["cluster1_audit_subcluster", "orig.ident", "cell_line", "nCount_RNA", "nFeature_RNA", "percent.mt", "cluster1_EDIL3_CRABP1_axis", "cluster1_ACKR3_ZEB2_axis"], plot_dir / "cluster1_local_umap_marker_axis.png", "Cluster 1 local UMAP marker-axis validation")
    boxplot_scores(c1, "orig.ident", ["cluster1_EDIL3_CRABP1_axis", "cluster1_ACKR3_ZEB2_axis"], plot_dir / "cluster1_module_scores_by_sample.png", "Cluster 1 module scores by sample")
    boxplot_scores(c1, "cell_line", ["cluster1_EDIL3_CRABP1_axis", "cluster1_ACKR3_ZEB2_axis"], plot_dir / "cluster1_module_scores_by_cell_line.png", "Cluster 1 module scores by cell line")
    c1_sample = composition_table(c1, "orig.ident", "cluster1_audit_subcluster", ["cluster1_EDIL3_CRABP1_axis", "cluster1_ACKR3_ZEB2_axis"], "n_cells_in_cluster1")
    c1_cell = composition_table(c1, "cell_line", "cluster1_audit_subcluster", ["cluster1_EDIL3_CRABP1_axis", "cluster1_ACKR3_ZEB2_axis"], "n_cells_in_cluster1")
    c1_corr = sample_correlation_table(c1, "cluster1_EDIL3_CRABP1_axis", "cluster1_ACKR3_ZEB2_axis", 100, "cluster1_audit_subcluster")
    c1_sample.to_csv(table_dir / "cluster1_sample_composition_table.tsv", sep="\t", index=False)
    c1_cell.to_csv(table_dir / "cluster1_cell_line_composition_table.tsv", sep="\t", index=False)
    c1_corr.to_csv(table_dir / "cluster1_within_sample_module_correlation_table.tsv", sep="\t", index=False)
    interp += ["## Cluster 1", "", f"Present EDIL3/CRABP1-axis genes: {', '.join(c1_a_present)}", f"Present ACKR3/ZEB2-axis genes: {', '.join(c1_b_present)}", ""]

    # Cluster 3
    c3_full = adata[adata.obs["cluster_id"] == "3"].copy()
    c3_a_present = add_module_score(c3_full, CLUSTER3_A, "cluster3_L3MBTL1_SYT4_axis")
    c3_b_present = add_module_score(c3_full, CLUSTER3_B, "cluster3_BEX3_CRABP1_DCX_axis")
    c3 = restore_scores(c3_full, local_split(c3_full, 0.2, "cluster3_audit_subcluster"), ["cluster3_L3MBTL1_SYT4_axis", "cluster3_BEX3_CRABP1_DCX_axis", "ribosomal_module", "stress_module"])
    save_scanpy_umap(c3, ["cluster3_audit_subcluster", "orig.ident", "cell_line", "nCount_RNA", "nFeature_RNA", "percent.mt", "cluster3_L3MBTL1_SYT4_axis", "cluster3_BEX3_CRABP1_DCX_axis"], plot_dir / "cluster3_local_umap_marker_axis.png", "Cluster 3 local UMAP marker-axis validation")
    c3_sample = composition_table(c3, "orig.ident", "cluster3_audit_subcluster", ["cluster3_BEX3_CRABP1_DCX_axis", "cluster3_L3MBTL1_SYT4_axis"], "n_cluster3_cells")
    c3_cell = composition_table(c3, "cell_line", "cluster3_audit_subcluster", ["cluster3_BEX3_CRABP1_DCX_axis", "cluster3_L3MBTL1_SYT4_axis"], "n_cluster3_cells")
    c3_sample.to_csv(table_dir / "cluster3_sample_composition_table.tsv", sep="\t", index=False)
    c3_cell.to_csv(table_dir / "cluster3_cell_line_composition_table.tsv", sep="\t", index=False)
    c3_qc = c3.obs.groupby("cluster3_audit_subcluster", observed=False)[["percent.mt", "ribosomal_module", "stress_module", "nCount_RNA", "nFeature_RNA", "S.Score", "G2M.Score"]].mean().reset_index()
    c3_qc.to_csv(table_dir / "cluster3_qc_marker_review_summary.tsv", sep="\t", index=False)
    cluster3_qc_plot(c3, "cluster3_audit_subcluster", plot_dir / "cluster3_qc_marker_review.png")
    neigh = adata[adata.obs["cluster_id"].isin(["1", "2", "3", "11"])].copy()
    neigh.obs["cluster3_neighbor_group"] = neigh.obs["cluster_id"].astype(str)
    neigh.obs.loc[c3.obs_names, "cluster3_neighbor_group"] = "3_sub" + c3.obs["cluster3_audit_subcluster"].astype(str)
    dot = dotplot_from_obs(neigh, "cluster3_neighbor_group", CLUSTER3_NEIGHBOR_MARKERS, plot_dir / "cluster3_neighbor_cluster_marker_dotplot.png", "Cluster 3 split versus neighboring states")
    dot.to_csv(table_dir / "cluster3_neighbor_cluster_marker_dotplot_table.tsv", sep="\t", index=False)
    interp += ["## Cluster 3", "", f"Present L3MBTL1/SYT4-axis genes: {', '.join(c3_a_present)}", f"Present BEX3/CRABP1/DCX-axis genes: {', '.join(c3_b_present)}", ""]

    # Cluster 11
    c11_full = adata[adata.obs["cluster_id"] == "11"].copy()
    c11_a_present = add_module_score(c11_full, CLUSTER11_A, "cluster11_ACKR3_ZEB2_MAF_axis")
    c11_b_present = add_module_score(c11_full, CLUSTER11_B, "cluster11_CPE_TUBB2A_CNTN1_axis")
    c11 = restore_scores(c11_full, local_split(c11_full, 0.4, "cluster11_audit_subcluster"), ["cluster11_ACKR3_ZEB2_MAF_axis", "cluster11_CPE_TUBB2A_CNTN1_axis"])
    save_scanpy_umap(c11, ["cluster11_audit_subcluster", "nCount_RNA", "nFeature_RNA", "percent.mt", "orig.ident", "cell_line", "cluster11_ACKR3_ZEB2_MAF_axis", "cluster11_CPE_TUBB2A_CNTN1_axis"], plot_dir / "cluster11_local_umap_marker_axis.png", "Cluster 11 local UMAP marker-axis validation")
    tert_nf = tertile_plot(c11, "nFeature_RNA", "cluster11_audit_subcluster", ["cluster11_ACKR3_ZEB2_MAF_axis", "cluster11_CPE_TUBB2A_CNTN1_axis"], plot_dir / "cluster11_nfeature_tertile_marker_scores.png")
    tert_nc = tertile_plot(c11, "nCount_RNA", "cluster11_audit_subcluster", ["cluster11_ACKR3_ZEB2_MAF_axis", "cluster11_CPE_TUBB2A_CNTN1_axis"], plot_dir / "cluster11_ncount_tertile_marker_scores.png")
    tert_nf.to_csv(table_dir / "cluster11_nfeature_tertile_marker_scores.tsv", sep="\t", index=False)
    tert_nc.to_csv(table_dir / "cluster11_ncount_tertile_marker_scores.tsv", sep="\t", index=False)
    model = model_cluster11(c11, "cluster11_audit_subcluster", ["cluster11_ACKR3_ZEB2_MAF_axis", "cluster11_CPE_TUBB2A_CNTN1_axis"])
    model.to_csv(table_dir / "cluster11_covariate_adjusted_model_results.tsv", sep="\t", index=False)
    c11_c3 = adata[adata.obs["cluster_id"].isin(["3", "11"])].copy()
    c11_c3.obs["cluster11_cluster3_group"] = c11_c3.obs["cluster_id"].astype(str)
    c11_c3.obs.loc[c11.obs_names, "cluster11_cluster3_group"] = "11_sub" + c11.obs["cluster11_audit_subcluster"].astype(str)
    c11_c3.obs.loc[c3.obs_names, "cluster11_cluster3_group"] = "3_sub" + c3.obs["cluster3_audit_subcluster"].astype(str)
    dot2 = dotplot_from_obs(c11_c3, "cluster11_cluster3_group", CLUSTER11_CLUSTER3_MARKERS, plot_dir / "cluster11_cluster3_comparison_marker_dotplot.png", "Cluster 11 split compared with cluster 3")
    dot2.to_csv(table_dir / "cluster11_cluster3_comparison_marker_dotplot_table.tsv", sep="\t", index=False)
    # Regress nCount/nFeature and recluster as a sensitivity check.
    c11_reg = c11_full.copy()
    c11_reg.layers["log_expr"] = c11_reg.X.copy()
    sc.pp.highly_variable_genes(c11_reg, n_top_genes=min(2000, c11_reg.n_vars), flavor="seurat")
    c11_reg = c11_reg[:, c11_reg.var["highly_variable"].to_numpy()].copy()
    sc.pp.regress_out(c11_reg, ["nCount_RNA", "nFeature_RNA"])
    sc.pp.scale(c11_reg, max_value=10)
    sc.tl.pca(c11_reg, n_comps=min(30, c11_reg.n_obs - 2, c11_reg.n_vars - 1), svd_solver="arpack")
    sc.pp.neighbors(c11_reg, n_neighbors=min(30, c11_reg.n_obs - 1), n_pcs=min(30, c11_reg.obsm["X_pca"].shape[1]))
    sc.tl.umap(c11_reg, random_state=0)
    sc.tl.leiden(c11_reg, resolution=0.4, key_added="cluster11_regressed_subcluster", random_state=0)
    c11_reg.obs["original_audit_subcluster"] = c11.obs.loc[c11_reg.obs_names, "cluster11_audit_subcluster"].astype(str).to_numpy()
    c11_reg.obs["cluster11_ACKR3_ZEB2_MAF_axis"] = c11.obs.loc[c11_reg.obs_names, "cluster11_ACKR3_ZEB2_MAF_axis"].to_numpy()
    c11_reg.obs["cluster11_CPE_TUBB2A_CNTN1_axis"] = c11.obs.loc[c11_reg.obs_names, "cluster11_CPE_TUBB2A_CNTN1_axis"].to_numpy()
    ari = adjusted_rand_score(c11_reg.obs["original_audit_subcluster"].astype(str), c11_reg.obs["cluster11_regressed_subcluster"].astype(str))
    pd.DataFrame([{"ari_original_vs_regressed_subcluster": ari}]).to_csv(table_dir / "cluster11_regressed_reclustering_comparison.tsv", sep="\t", index=False)
    save_scanpy_umap(c11_reg, ["cluster11_regressed_subcluster", "original_audit_subcluster", "cluster11_ACKR3_ZEB2_MAF_axis", "cluster11_CPE_TUBB2A_CNTN1_axis"], plot_dir / "cluster11_regressed_ncount_nfeature_local_umap.png", "Cluster 11 after regressing nCount/nFeature")
    interp += ["## Cluster 11", "", f"Present ACKR3/ZEB2/MAF-axis genes: {', '.join(c11_a_present)}", f"Present CPE/TUBB2A/CNTN1-axis genes: {', '.join(c11_b_present)}", f"ARI original audit split vs nCount/nFeature-regressed split: {ari:.3f}", ""]

    write_interpretation(outdir, interp)
    print(f"Wrote validation outputs to {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
