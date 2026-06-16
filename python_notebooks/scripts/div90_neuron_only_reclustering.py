#!/usr/bin/env python3
"""Recluster the DIV90 neuronal compartment as one object."""

from __future__ import annotations

import argparse
import gzip
import os
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


PROJECT_ROOT_DEFAULT = "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder"

NEURON_PARENT_CLUSTERS = ["0", "1", "2", "3", "5", "8", "11"]
EXCLUDED_PARENT_CLUSTERS = ["4", "6", "7", "9", "10", "12"]

CANONICAL_MARKERS = [
    "LHX6", "SOX6", "ERBB4", "SATB1", "MAF", "MAFB", "SST", "NPY",
    "GAD1", "GAD2", "KCNC1", "KCNC2",
    "DCX", "STMN2", "DCLK1", "ACKR3", "CXCR4", "ZEB2",
    "CRABP1", "BEX3", "CPE", "TUBB2A", "CNTN1",
    "LHX8", "ISL1", "GBX2", "TAC1", "NKX2-1",
]

NUMERIC_OBS = ["nCount_RNA", "nFeature_RNA", "percent.mt", "S.Score", "G2M.Score", "CC.Difference"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", default=os.environ.get("PROJECT_ROOT", PROJECT_ROOT_DEFAULT))
    parser.add_argument("--h5ad", default=None)
    parser.add_argument("--sample-map", default="metadata/div30_div90_sample_id_to_biolabel_map.tsv")
    parser.add_argument(
        "--outdir",
        default=None,
        help="Default: PROJECT_ROOT/results/div90_neuron_only_reclustering/div90_neuron_only_reclustering_v1",
    )
    parser.add_argument("--resolutions", default="0.1,0.2,0.3,0.4,0.6,0.8,1.0")
    parser.add_argument("--n-top-genes", type=int, default=2500)
    parser.add_argument("--n-pcs", type=int, default=40)
    parser.add_argument("--n-neighbors", type=int, default=30)
    parser.add_argument("--marker-top-n", type=int, default=100)
    return parser.parse_args()


def split_csv(value: str) -> list[str]:
    return [x.strip() for x in value.split(",") if x.strip()]


def resolution_key(resolution: float) -> str:
    return f"neuron_leiden_r{str(resolution).replace('.', '_')}"


def natural_sort(values: list[str]) -> list[str]:
    def key(value: str) -> tuple[int, str]:
        try:
            return (0, f"{int(float(str(value))):05d}")
        except ValueError:
            return (1, str(value))

    return sorted(values, key=key)


def to_numeric_obs(adata: ad.AnnData) -> None:
    for col in NUMERIC_OBS:
        if col in adata.obs:
            adata.obs[col] = pd.to_numeric(adata.obs[col].astype(str), errors="coerce")


def add_cell_line(adata: ad.AnnData, sample_map_path: Path) -> None:
    if "cell_line" in adata.obs:
        return
    if not sample_map_path.exists() or "orig.ident" not in adata.obs:
        adata.obs["cell_line"] = "unknown"
        return
    sample_map = pd.read_csv(sample_map_path, sep="\t")
    div90 = sample_map.loc[sample_map["DIV"].astype(str).str.upper() == "DIV90"].copy()
    lookup = div90.set_index("run_sample_id")["biological_label"].astype(str).to_dict()
    labels = adata.obs["orig.ident"].astype(str).map(lookup).fillna("unknown")
    adata.obs["cell_line"] = labels.str.replace(r"_rep\d+$", "", regex=True)


def make_dirs(outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)


def matrix_values(adata: ad.AnnData, genes: list[str]) -> tuple[np.ndarray, list[str]]:
    present = [gene for gene in genes if gene in adata.var_names]
    if not present:
        return np.zeros((adata.n_obs, 0), dtype=float), []
    mat = adata[:, present].X
    if sparse.issparse(mat):
        mat = mat.toarray()
    return np.asarray(mat, dtype=float), present


def write_tsv_gz(df: pd.DataFrame, path: Path) -> None:
    with gzip.open(path, "wt") as handle:
        df.to_csv(handle, sep="\t", index=False)


def save_scanpy_umap(adata: ad.AnnData, colors: list[str], path: Path, title: str, ncols: int = 3) -> None:
    fig = sc.pl.umap(
        adata,
        color=colors,
        ncols=ncols,
        show=False,
        frameon=False,
        size=8,
        return_fig=True,
    )
    fig.suptitle(title, y=1.01)
    fig.savefig(path, dpi=240, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def composition_table(adata: ad.AnnData, group_col: str, covariate_col: str, resolutions: list[float]) -> pd.DataFrame:
    rows = []
    for res in resolutions:
        key = resolution_key(res)
        obs = adata.obs[[key, covariate_col]].copy()
        obs[key] = obs[key].astype(str)
        obs[covariate_col] = obs[covariate_col].astype(str)
        sub_sizes = obs[key].value_counts().to_dict()
        cov_sizes = obs[covariate_col].value_counts().to_dict()
        for (subcluster, covariate), n_cells in obs.groupby([key, covariate_col], observed=False).size().items():
            rows.append(
                {
                    "resolution": res,
                    "subcluster_key": key,
                    "subcluster": str(subcluster),
                    group_col: str(covariate),
                    "n_cells": int(n_cells),
                    "fraction_of_subcluster": float(n_cells / sub_sizes[str(subcluster)]),
                    "fraction_of_group": float(n_cells / cov_sizes[str(covariate)]),
                    "subcluster_size": int(sub_sizes[str(subcluster)]),
                    "group_size": int(cov_sizes[str(covariate)]),
                }
            )
    return pd.DataFrame(rows).sort_values(["resolution", "subcluster", "n_cells"], ascending=[True, True, False])


def parent_composition_table(adata: ad.AnnData, resolutions: list[float]) -> pd.DataFrame:
    rows = []
    label_lookup = (
        adata.obs[["cluster_id", "cluster_number_name"]]
        .drop_duplicates()
        .assign(cluster_id=lambda df: df["cluster_id"].astype(str))
        .set_index("cluster_id")["cluster_number_name"]
        .astype(str)
        .to_dict()
    )
    for res in resolutions:
        key = resolution_key(res)
        obs = adata.obs[[key, "cluster_id"]].copy()
        obs[key] = obs[key].astype(str)
        obs["cluster_id"] = obs["cluster_id"].astype(str)
        sub_sizes = obs[key].value_counts().to_dict()
        parent_sizes = obs["cluster_id"].value_counts().to_dict()
        for (subcluster, parent), n_cells in obs.groupby([key, "cluster_id"], observed=False).size().items():
            rows.append(
                {
                    "resolution": res,
                    "subcluster_key": key,
                    "subcluster": str(subcluster),
                    "parent_cluster_id": str(parent),
                    "parent_cluster_name": label_lookup.get(str(parent), str(parent)),
                    "n_cells": int(n_cells),
                    "fraction_of_subcluster": float(n_cells / sub_sizes[str(subcluster)]),
                    "fraction_of_parent": float(n_cells / parent_sizes[str(parent)]),
                    "subcluster_size": int(sub_sizes[str(subcluster)]),
                    "parent_cluster_size": int(parent_sizes[str(parent)]),
                }
            )
    return pd.DataFrame(rows).sort_values(["resolution", "subcluster", "n_cells"], ascending=[True, True, False])


def marker_count_summary(markers: pd.DataFrame) -> pd.DataFrame:
    if markers.empty:
        return pd.DataFrame(columns=["resolution", "subcluster", "n_positive_markers", "top_marker"])
    df = markers.copy()
    df["logfoldchanges"] = pd.to_numeric(df["logfoldchanges"], errors="coerce")
    df["pvals_adj"] = pd.to_numeric(df["pvals_adj"], errors="coerce")
    positive = df.loc[(df["logfoldchanges"] >= 0.25) & (df["pvals_adj"] <= 0.05)].copy()
    rows = []
    for (res, group), sub in df.groupby(["resolution", "group"], sort=False):
        pos = positive.loc[(positive["resolution"] == res) & (positive["group"].astype(str) == str(group))]
        top = pos.iloc[0] if not pos.empty else sub.sort_values(["pvals_adj", "logfoldchanges"], ascending=[True, False]).iloc[0]
        rows.append(
            {
                "resolution": float(res),
                "subcluster": str(group),
                "n_positive_markers": int(pos.shape[0]),
                "top_marker": str(top["names"]),
                "top_marker_logfc": float(top["logfoldchanges"]) if pd.notna(top["logfoldchanges"]) else np.nan,
                "top_marker_padj": float(top["pvals_adj"]) if pd.notna(top["pvals_adj"]) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def dominant_metric(composition: pd.DataFrame, value_col: str) -> pd.DataFrame:
    rows = []
    for (res, sub), df in composition.groupby(["resolution", "subcluster"], observed=False):
        top = df.sort_values("fraction_of_subcluster", ascending=False).iloc[0]
        rows.append(
            {
                "resolution": float(res),
                "subcluster": str(sub),
                f"dominant_{value_col}": str(top[value_col]),
                f"dominant_{value_col}_fraction": float(top["fraction_of_subcluster"]),
            }
        )
    return pd.DataFrame(rows)


def resolution_summary(
    adata: ad.AnnData,
    resolutions: list[float],
    parent_comp: pd.DataFrame,
    sample_comp: pd.DataFrame,
    cell_line_comp: pd.DataFrame,
    markers: pd.DataFrame,
) -> pd.DataFrame:
    marker_counts = marker_count_summary(markers)
    parent_dom = dominant_metric(parent_comp, "parent_cluster_id")
    sample_dom = dominant_metric(sample_comp, "orig.ident")
    line_dom = dominant_metric(cell_line_comp, "cell_line")
    rows = []
    for res in resolutions:
        key = resolution_key(res)
        counts = adata.obs[key].astype(str).value_counts()
        m = marker_counts.loc[marker_counts["resolution"] == res]
        n_marker_supported = int((m["n_positive_markers"] >= 2).sum()) if not m.empty else 0
        merged = (
            parent_dom.loc[parent_dom["resolution"] == res]
            .merge(sample_dom.loc[sample_dom["resolution"] == res], on=["resolution", "subcluster"], how="outer")
            .merge(line_dom.loc[line_dom["resolution"] == res], on=["resolution", "subcluster"], how="outer")
        )
        rows.append(
            {
                "resolution": res,
                "subcluster_key": key,
                "n_neuronal_subclusters": int(counts.shape[0]),
                "min_subcluster_size": int(counts.min()),
                "median_subcluster_size": float(counts.median()),
                "max_subcluster_size": int(counts.max()),
                "marker_supported_subclusters": n_marker_supported,
                "marker_supported_fraction": float(n_marker_supported / counts.shape[0]),
                "mean_parent_purity": float(merged["dominant_parent_cluster_id_fraction"].mean()),
                "max_parent_purity": float(merged["dominant_parent_cluster_id_fraction"].max()),
                "mean_sample_dominance": float(merged["dominant_orig.ident_fraction"].mean()),
                "max_sample_dominance": float(merged["dominant_orig.ident_fraction"].max()),
                "mean_cell_line_dominance": float(merged["dominant_cell_line_fraction"].mean()),
                "max_cell_line_dominance": float(merged["dominant_cell_line_fraction"].max()),
            }
        )
    summary = pd.DataFrame(rows)
    summary["recommendation_score"] = summary.apply(score_resolution, axis=1)
    best_idx = summary.sort_values(["recommendation_score", "resolution"], ascending=[False, True]).index[0]
    summary["recommended"] = False
    summary.loc[best_idx, "recommended"] = True
    return summary


def score_resolution(row: pd.Series) -> float:
    score = 0.0
    n_clusters = float(row["n_neuronal_subclusters"])
    min_size = float(row["min_subcluster_size"])
    score += 4.0 * float(row["marker_supported_fraction"])
    score += 1.0 if 7 <= n_clusters <= 12 else -abs(n_clusters - 9) * 0.25
    score += min(min_size / 150.0, 1.0)
    score += 0.8 * float(row["mean_parent_purity"])
    score -= max(0.0, float(row["max_cell_line_dominance"]) - 0.85) * 2.0
    score -= max(0.0, float(row["max_sample_dominance"]) - 0.75) * 2.0
    return score


def marker_dotplot_table(adata: ad.AnnData, group_col: str, genes: list[str]) -> pd.DataFrame:
    values, present = matrix_values(adata, genes)
    expr = pd.DataFrame(values, index=adata.obs_names, columns=present)
    expr[group_col] = adata.obs[group_col].astype(str).to_numpy()
    rows = []
    for group, df in expr.groupby(group_col):
        for gene in present:
            vals = pd.to_numeric(df[gene], errors="coerce").to_numpy(float)
            rows.append(
                {
                    group_col: str(group),
                    "gene": gene,
                    "mean_expr": float(np.nanmean(vals)),
                    "pct_expressed": float(np.mean(vals > 0) * 100),
                }
            )
    return pd.DataFrame(rows)


def plot_marker_dotplot(dot: pd.DataFrame, group_col: str, genes: list[str], path: Path, title: str) -> None:
    present = [gene for gene in genes if gene in set(dot["gene"])]
    groups = natural_sort(list(dot[group_col].astype(str).unique()))
    fig, ax = plt.subplots(figsize=(max(12, 0.42 * len(present) + 3), max(4.5, 0.38 * len(groups) + 2)))
    x = {gene: i for i, gene in enumerate(present)}
    y = {group: i for i, group in enumerate(groups)}
    plot_df = dot.loc[dot["gene"].isin(present)].copy()
    ax.scatter(
        [x[g] for g in plot_df["gene"]],
        [y[str(g)] for g in plot_df[group_col]],
        s=np.maximum(8, plot_df["pct_expressed"].to_numpy(float) * 3.0),
        c=plot_df["mean_expr"].to_numpy(float),
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
    ax.set_xlabel("")
    ax.set_ylabel("Neuron-only Leiden subcluster")
    fig.colorbar(ax.collections[0], ax=ax, label="Mean expression")
    fig.tight_layout()
    fig.savefig(path, dpi=240)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def describe_recommended_clusters(
    adata: ad.AnnData,
    recommended_key: str,
    parent_comp: pd.DataFrame,
    cell_line_comp: pd.DataFrame,
    marker_dot: pd.DataFrame,
    marker_counts: pd.DataFrame,
) -> pd.DataFrame:
    parent = parent_comp.loc[parent_comp["subcluster_key"] == recommended_key].copy()
    line = cell_line_comp.loc[cell_line_comp["subcluster_key"] == recommended_key].copy()
    marker_counts = marker_counts.loc[marker_counts["resolution"] == float(recommended_key.split("_r")[-1].replace("_", "."))].copy()
    rows = []
    for subcluster in natural_sort(list(adata.obs[recommended_key].astype(str).unique())):
        p = parent.loc[parent["subcluster"] == subcluster].sort_values("fraction_of_subcluster", ascending=False)
        l = line.loc[line["subcluster"] == subcluster].sort_values("fraction_of_subcluster", ascending=False)
        d = marker_dot.loc[marker_dot[recommended_key].astype(str) == subcluster].copy()
        top_means = d.sort_values("mean_expr", ascending=False).head(8)["gene"].tolist()
        top_pct = d.sort_values("pct_expressed", ascending=False).head(8)["gene"].tolist()
        mean = d.set_index("gene")["mean_expr"].to_dict()
        neuroblast_score = np.nanmean([mean.get(g, np.nan) for g in ["DCX", "STMN2", "DCLK1", "BEX3", "CRABP1", "TUBB2A", "CPE", "CNTN1"]])
        cortical_score = np.nanmean([mean.get(g, np.nan) for g in ["SST", "NPY", "SATB1", "MAF", "ERBB4", "GAD1", "GAD2"]])
        ventral_score = np.nanmean([mean.get(g, np.nan) for g in ["LHX8", "ISL1", "GBX2", "TAC1", "NKX2-1"]])
        pv_score = np.nanmean([mean.get(g, np.nan) for g in ["SOX6", "KCNC1", "KCNC2", "ERBB4", "GAD1", "GAD2"]])
        if neuroblast_score >= max(cortical_score, ventral_score, pv_score) and neuroblast_score > 1.0:
            state_class = "intermediate/neuroblast-like"
        elif cortical_score >= ventral_score and cortical_score >= pv_score:
            state_class = "candidate cortical interneuron terminal/tip"
        elif ventral_score >= cortical_score:
            state_class = "candidate ventral MGE/striatal-GP terminal/tip"
        else:
            state_class = "candidate PV/maturing interneuron state"
        dominant_line_fraction = float(l.iloc[0]["fraction_of_subcluster"]) if not l.empty else np.nan
        if pd.notna(dominant_line_fraction) and dominant_line_fraction >= 0.70:
            state_class = f"{state_class}; line-biased"
        markers_for_sub = marker_counts.loc[marker_counts["subcluster"].astype(str) == subcluster]
        rows.append(
            {
                "subcluster": subcluster,
                "n_cells": int((adata.obs[recommended_key].astype(str) == subcluster).sum()),
                "dominant_parent_cluster_id": p.iloc[0]["parent_cluster_id"] if not p.empty else "",
                "dominant_parent_cluster_name": p.iloc[0]["parent_cluster_name"] if not p.empty else "",
                "dominant_parent_fraction": float(p.iloc[0]["fraction_of_subcluster"]) if not p.empty else np.nan,
                "dominant_cell_line": l.iloc[0]["cell_line"] if not l.empty else "",
                "dominant_cell_line_fraction": dominant_line_fraction,
                "top_marker": markers_for_sub.iloc[0]["top_marker"] if not markers_for_sub.empty else "",
                "n_positive_markers": int(markers_for_sub.iloc[0]["n_positive_markers"]) if not markers_for_sub.empty else 0,
                "high_mean_canonical_markers": ", ".join(top_means),
                "high_pct_canonical_markers": ", ".join(top_pct),
                "state_class": state_class,
            }
        )
    return pd.DataFrame(rows)


def write_report(
    outdir: Path,
    adata: ad.AnnData,
    summary: pd.DataFrame,
    cluster_desc: pd.DataFrame,
    canonical_present: list[str],
    canonical_missing: list[str],
) -> None:
    rec = summary.loc[summary["recommended"]].iloc[0]
    recommended_key = rec["subcluster_key"]
    summary_md = simple_markdown_table(
        summary[
            [
                "resolution",
                "n_neuronal_subclusters",
                "min_subcluster_size",
                "marker_supported_fraction",
                "mean_parent_purity",
                "max_cell_line_dominance",
                "recommended",
            ]
        ]
    )
    cluster_md = simple_markdown_table(cluster_desc)
    lines = [
        "# DIV90 Neuron-Only Reclustering",
        "",
        "This analysis reclustered the DIV90 neuronal compartment as one object.",
        "It did not recluster each parent cluster separately.",
        "",
        "## Inputs",
        "",
        "```text",
        "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/python_anndata/varela_div90.h5ad",
        "```",
        "",
        f"Neuron parent clusters included: {', '.join(NEURON_PARENT_CLUSTERS)}",
        f"Parent clusters excluded: {', '.join(EXCLUDED_PARENT_CLUSTERS)}",
        f"Neuron-only cells: {adata.n_obs:,}",
        "",
        "Here, `parent cluster` means the original DIV90 cluster labels that already existed in the input object,",
        "stored in `cluster_id` and `cluster_number_name`. It does not mean a lineage-tree parent or a URD parent.",
        "The neuron-only Leiden clusters are new clusters computed after subsetting to neuronal parents only.",
        "",
        "## Resolution Recommendation",
        "",
        f"Recommended resolution: `{rec['resolution']}`",
        f"Recommended subcluster key: `{recommended_key}`",
        f"Number of neuronal states: `{int(rec['n_neuronal_subclusters'])}`",
        "",
        "The recommendation prioritizes interpretable marker-supported states, avoids very small clusters,",
        "and checks parent-cluster, sample, and cell-line representation.",
        "",
        "This recommendation is a pragmatic scoring rule, not a formal hypothesis test. Each resolution received:",
        "",
        "```text",
        "score = 4 * marker_supported_fraction",
        "      + cluster_count_bonus",
        "      + min(min_subcluster_size / 150, 1)",
        "      + 0.8 * mean_parent_purity",
        "      - 2 * max(0, max_cell_line_dominance - 0.85)",
        "      - 2 * max(0, max_sample_dominance - 0.75)",
        "",
        "cluster_count_bonus = +1 if 7-12 clusters, otherwise -0.25 * distance from 9 clusters",
        "```",
        "",
        "A subcluster was counted as marker-supported if it had at least two positive markers among the ranked",
        "genes with adjusted P <= 0.05 and log fold-change >= 0.25.",
        "",
        "So `good enough` means: marker-supported, no tiny cluster problem, enough clusters to separate major",
        "neuronal programs, not so many clusters that the map fragments, and no single sample/cell line completely",
        "dominating the result.",
        "",
        "## Resolution Summary",
        "",
        summary_md,
        "",
        "## Recommended Neuron-Only Clusters",
        "",
        cluster_md,
        "",
        "## Interpretation Guide",
        "",
        "Original parent-cluster correspondence is captured by the dominant parent and parent fraction columns.",
        "Candidate terminal/tip states are the clusters classified as cortical interneuron, ventral MGE/striatal-GP,",
        "or PV/maturing interneuron states. Candidate intermediate states are those classified as",
        "intermediate/neuroblast-like, driven by markers such as `DCX`, `STMN2`, `DCLK1`, `BEX3`, `CRABP1`,",
        "`TUBB2A`, `CPE`, and `CNTN1`.",
        "",
        "Line-biased but biologically plausible clusters are marked with `line-biased` in the state-class column.",
        "These should not be discarded automatically because cell line may track real differentiation differences,",
        "but they should be validated before becoming official terminal states.",
        "",
        "## Canonical Marker Availability",
        "",
        f"Present markers: {', '.join(canonical_present)}",
        f"Missing markers: {', '.join(canonical_missing) if canonical_missing else 'none'}",
        "",
        "## Deliverables",
        "",
        "```text",
        "neuron_only_reclustering_report.md",
        "neuron_only_resolution_summary.tsv",
        "neuron_only_subcluster_markers.tsv.gz",
        "neuron_only_subcluster_parent_composition.tsv",
        "neuron_only_subcluster_cell_line_composition.tsv",
        "neuron_only_subcluster_sample_composition.tsv",
        "neuron_only_umap_parent_cluster.png/pdf",
        "neuron_only_umap_reclustered_resolutions.png/pdf",
        "neuron_only_umap_cell_line.png/pdf",
        "neuron_only_umap_orig_ident.png/pdf",
        "neuron_only_marker_dotplot.png/pdf",
        "```",
    ]
    (outdir / "neuron_only_reclustering_report.md").write_text("\n".join(lines) + "\n")


def simple_markdown_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "_No rows._"
    display = df.copy()
    for col in display.columns:
        if pd.api.types.is_float_dtype(display[col]):
            display[col] = display[col].map(lambda x: "" if pd.isna(x) else f"{x:.3g}")
        else:
            display[col] = display[col].map(lambda x: "" if pd.isna(x) else str(x))
    columns = [str(col) for col in display.columns]
    rows = display.astype(str).values.tolist()
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(row) + " |")
    return "\n".join(lines)


def main() -> int:
    args = parse_args()
    project_root = Path(args.project_root).expanduser().resolve()
    h5ad = Path(args.h5ad).expanduser().resolve() if args.h5ad else project_root / "results/python_anndata/varela_div90.h5ad"
    outdir = (
        Path(args.outdir).expanduser().resolve()
        if args.outdir
        else project_root / "results/div90_neuron_only_reclustering/div90_neuron_only_reclustering_v1"
    )
    sample_map = Path(args.sample_map)
    if not sample_map.is_absolute():
        sample_map = Path.cwd() / sample_map
    make_dirs(outdir)

    print(f"Reading {h5ad}")
    adata = sc.read_h5ad(h5ad)
    adata.obs["cluster_id"] = adata.obs["cluster_id"].astype(str)
    to_numeric_obs(adata)
    add_cell_line(adata, sample_map)

    neuron = adata[adata.obs["cluster_id"].isin(NEURON_PARENT_CLUSTERS)].copy()
    neuron.layers["log_expr"] = neuron.X.copy()
    neuron.obs["parent_cluster"] = neuron.obs["cluster_id"].astype(str) + " - " + neuron.obs["cluster_number_name"].astype(str)

    resolutions = [float(x) for x in split_csv(args.resolutions)]
    print(f"Neuron-only subset: {neuron.n_obs} cells x {neuron.n_vars} genes")

    work = neuron.copy()
    sc.pp.highly_variable_genes(work, n_top_genes=min(args.n_top_genes, work.n_vars), flavor="seurat")
    work = work[:, work.var["highly_variable"].to_numpy()].copy()
    sc.pp.scale(work, max_value=10)
    sc.tl.pca(work, n_comps=min(args.n_pcs, work.n_obs - 2, work.n_vars - 1), svd_solver="arpack")
    sc.pp.neighbors(work, n_neighbors=min(args.n_neighbors, work.n_obs - 1), n_pcs=min(args.n_pcs, work.obsm["X_pca"].shape[1]))
    sc.tl.umap(work, random_state=0)

    neuron.obsm["X_umap"] = work.obsm["X_umap"].copy()
    all_markers = []
    for res in resolutions:
        key = resolution_key(res)
        print(f"Leiden resolution {res}")
        sc.tl.leiden(work, resolution=res, key_added=key, random_state=0)
        neuron.obs[key] = work.obs[key].astype(str).to_numpy()
        if neuron.obs[key].nunique() >= 2:
            sc.tl.rank_genes_groups(neuron, groupby=key, method="wilcoxon", pts=True, layer="log_expr", n_genes=args.marker_top_n)
            markers = sc.get.rank_genes_groups_df(neuron, group=None)
            markers["resolution"] = res
            markers["subcluster_key"] = key
            all_markers.append(markers)
    marker_df = pd.concat(all_markers, ignore_index=True) if all_markers else pd.DataFrame()

    parent_comp = parent_composition_table(neuron, resolutions)
    sample_comp = composition_table(neuron, "orig.ident", "orig.ident", resolutions)
    cell_line_comp = composition_table(neuron, "cell_line", "cell_line", resolutions)
    summary = resolution_summary(neuron, resolutions, parent_comp, sample_comp, cell_line_comp, marker_df)
    recommended = summary.loc[summary["recommended"]].iloc[0]
    recommended_key = str(recommended["subcluster_key"])
    marker_counts = marker_count_summary(marker_df)

    canonical_present = [gene for gene in CANONICAL_MARKERS if gene in neuron.var_names]
    canonical_missing = [gene for gene in CANONICAL_MARKERS if gene not in neuron.var_names]
    dot = marker_dotplot_table(neuron, recommended_key, canonical_present)
    cluster_desc = describe_recommended_clusters(neuron, recommended_key, parent_comp, cell_line_comp, dot, marker_counts)

    summary.to_csv(outdir / "neuron_only_resolution_summary.tsv", sep="\t", index=False)
    parent_comp.to_csv(outdir / "neuron_only_subcluster_parent_composition.tsv", sep="\t", index=False)
    cell_line_comp.to_csv(outdir / "neuron_only_subcluster_cell_line_composition.tsv", sep="\t", index=False)
    sample_comp.to_csv(outdir / "neuron_only_subcluster_sample_composition.tsv", sep="\t", index=False)
    cluster_desc.to_csv(outdir / "neuron_only_recommended_subcluster_interpretation.tsv", sep="\t", index=False)
    dot.to_csv(outdir / "neuron_only_marker_dotplot_table.tsv", sep="\t", index=False)
    write_tsv_gz(marker_df, outdir / "neuron_only_subcluster_markers.tsv.gz")

    save_scanpy_umap(neuron, ["parent_cluster"], outdir / "neuron_only_umap_parent_cluster.png", "DIV90 neuron-only UMAP by original parent cluster", ncols=1)
    save_scanpy_umap(neuron, [resolution_key(res) for res in resolutions], outdir / "neuron_only_umap_reclustered_resolutions.png", "DIV90 neuron-only Leiden resolution sweep", ncols=3)
    save_scanpy_umap(neuron, ["cell_line"], outdir / "neuron_only_umap_cell_line.png", "DIV90 neuron-only UMAP by cell line", ncols=1)
    save_scanpy_umap(neuron, ["orig.ident"], outdir / "neuron_only_umap_orig_ident.png", "DIV90 neuron-only UMAP by sample", ncols=1)
    plot_marker_dotplot(dot, recommended_key, canonical_present, outdir / "neuron_only_marker_dotplot.png", f"Canonical markers at neuron-only resolution {recommended['resolution']}")

    write_report(outdir, neuron, summary, cluster_desc, canonical_present, canonical_missing)
    print(f"Wrote neuron-only reclustering outputs to {outdir}")
    print(f"Recommended resolution: {recommended['resolution']} ({int(recommended['n_neuronal_subclusters'])} states)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
