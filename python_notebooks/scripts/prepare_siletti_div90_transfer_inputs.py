#!/usr/bin/env python3
"""Prepare DIV90 neuron manifests and figure-A heatmaps for Siletti mapping."""

from __future__ import annotations

import argparse
import json
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
from scipy import sparse


PROJECT_ROOT_DEFAULT = "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder"
DEFAULT_RUN_LABEL = "siletti_div90_neuron_prep_v1"
EXCLUDE_NAME_PATTERNS = ("astro", "opc", "stressed", "dividing")
CURATED_MARKERS = [
    "DCX",
    "STMN2",
    "DCLK1",
    "GAD1",
    "GAD2",
    "LHX6",
    "SOX6",
    "SATB1",
    "MAF",
    "MAFB",
    "ERBB4",
    "PVALB",
    "PPARGC1A",
    "KCNC1",
    "KCNC2",
    "GPR149",
    "SST",
    "NPY",
    "CORT",
    "NOS1",
    "CRABP1",
    "ANGPT2",
    "LHX8",
    "ISL1",
    "CHAT",
    "SLC5A7",
    "NKX2-1",
    "FOXP1",
    "BCL11B",
    "TAC1",
    "GBX2",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", default=os.environ.get("PROJECT_ROOT", PROJECT_ROOT_DEFAULT))
    parser.add_argument("--h5ad", default=None)
    parser.add_argument("--run-label", default=os.environ.get("SILETTI_DIV90_PREP_RUN_LABEL", DEFAULT_RUN_LABEL))
    parser.add_argument("--outdir", default=None)
    parser.add_argument("--cluster-col", default="cluster_id")
    parser.add_argument("--cluster-name-col", default="cluster_number_name")
    parser.add_argument("--top-de-genes-per-class", type=int, default=8)
    parser.add_argument("--max-heatmap-genes", type=int, default=120)
    parser.add_argument("--de-method", default="wilcoxon")
    parser.add_argument("--seed", type=int, default=0)
    return parser.parse_args()


def natural_key(value: object) -> tuple[int, str]:
    text = str(value)
    try:
        return (0, f"{int(float(text)):05d}")
    except ValueError:
        return (1, text)


def strip_cluster_prefix(name: object) -> str:
    text = str(name)
    if " - " in text:
        return text.split(" - ", 1)[1]
    return text


def include_for_transfer(cluster_name: str) -> bool:
    low = cluster_name.lower()
    return not any(pattern in low for pattern in EXCLUDE_NAME_PATTERNS)


def cortical_subpallial_hypothesis(cluster_name: str) -> str:
    low = cluster_name.lower()
    if "striatal" in low or "gp" in low or "lhx8" in low:
        return "subpallial_hypothesis"
    if "sst" in low or "cortical" in low or "pv precursor" in low or "crabp1" in low:
        return "cortical_or_pv_precursor_hypothesis"
    return "excluded_or_uncertain"


def matrix_column(adata: ad.AnnData, gene: str) -> np.ndarray:
    idx = adata.var_names.get_loc(gene)
    col = adata.X[:, idx]
    if sparse.issparse(col):
        col = col.toarray()
    return np.asarray(col).reshape(-1)


def zscore_rows(frame: pd.DataFrame) -> pd.DataFrame:
    values = frame.to_numpy(dtype=float)
    means = np.nanmean(values, axis=1, keepdims=True)
    sds = np.nanstd(values, axis=1, keepdims=True)
    sds[~np.isfinite(sds) | (sds == 0)] = 1.0
    z = (values - means) / sds
    return pd.DataFrame(z, index=frame.index, columns=frame.columns)


def save_heatmap(frame: pd.DataFrame, path: Path, title: str, *, center: float | None = 0) -> None:
    width = max(8, 0.36 * frame.shape[1] + 4)
    height = max(6, 0.18 * frame.shape[0] + 2.5)
    fig, ax = plt.subplots(figsize=(width, height))
    sns.heatmap(frame, cmap="RdBu_r" if center == 0 else "viridis", center=center, ax=ax, cbar_kws={"shrink": 0.65})
    ax.set_title(title)
    ax.set_xlabel("DIV90 broad neuron class")
    ax.set_ylabel("Gene")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    args = parse_args()
    np.random.seed(args.seed)
    project_root = Path(args.project_root)
    h5ad = Path(args.h5ad) if args.h5ad else project_root / "results/python_anndata/varela_div90.h5ad"
    outdir = (
        Path(args.outdir)
        if args.outdir
        else project_root / "results/siletti_2023_whb_reference_label_transfer" / args.run_label
    )
    table_dir = outdir / "tables"
    plot_dir = outdir / "plots"
    report_dir = outdir / "reports"
    for directory in [table_dir, plot_dir, report_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(h5ad)
    if args.cluster_col not in adata.obs:
        raise ValueError(f"Missing cluster column: {args.cluster_col}")
    if args.cluster_name_col not in adata.obs:
        raise ValueError(f"Missing cluster-name column: {args.cluster_name_col}")

    obs = adata.obs[[args.cluster_col, args.cluster_name_col, "orig.ident"]].copy()
    obs[args.cluster_col] = obs[args.cluster_col].astype(str)
    obs["div90_broad_class"] = obs[args.cluster_name_col].map(strip_cluster_prefix)
    cluster_manifest = (
        obs[[args.cluster_col, args.cluster_name_col, "div90_broad_class"]]
        .drop_duplicates()
        .sort_values(args.cluster_col, key=lambda s: s.map(natural_key))
    )
    cluster_manifest["include_for_siletti_transfer"] = cluster_manifest["div90_broad_class"].map(include_for_transfer)
    cluster_manifest["cortical_or_subpallial_hypothesis"] = cluster_manifest["div90_broad_class"].map(cortical_subpallial_hypothesis)
    cluster_manifest["exclusion_reason"] = np.where(
        cluster_manifest["include_for_siletti_transfer"],
        "",
        "non_neuronal_or_low_quality_cluster_name_pattern",
    )
    counts = obs.groupby(args.cluster_col).size().rename("n_cells").reset_index()
    cluster_manifest = cluster_manifest.merge(counts, on=args.cluster_col, how="left")
    cluster_manifest.to_csv(table_dir / "div90_neuron_broad_class_manifest.tsv", sep="\t", index=False)

    query_obs = obs.merge(
        cluster_manifest[[args.cluster_col, "include_for_siletti_transfer", "cortical_or_subpallial_hypothesis"]],
        on=args.cluster_col,
        how="left",
    )
    query_obs["cell_id"] = adata.obs["cell_id"].astype(str).to_numpy() if "cell_id" in adata.obs else adata.obs_names.astype(str)
    query_obs["obs_name"] = adata.obs_names.astype(str)
    query_obs.to_csv(table_dir / "div90_siletti_query_cell_manifest.tsv.gz", sep="\t", index=False)
    query_cells = query_obs.loc[query_obs["include_for_siletti_transfer"], ["cell_id", "obs_name", args.cluster_col, "div90_broad_class"]]
    query_cells.to_csv(table_dir / "div90_siletti_query_neuron_cells.tsv", sep="\t", index=False)

    present_markers = [gene for gene in CURATED_MARKERS if gene in adata.var_names]
    missing_markers = [gene for gene in CURATED_MARKERS if gene not in adata.var_names]
    gene_match = pd.DataFrame(
        [{"gene": gene, "present": gene in adata.var_names} for gene in CURATED_MARKERS]
    )
    gene_match.to_csv(table_dir / "div90_siletti_marker_gene_match.tsv", sep="\t", index=False)

    neuron_mask = query_obs["include_for_siletti_transfer"].to_numpy(dtype=bool)
    adata_neuron = adata[neuron_mask].copy()
    neuron_obs = query_obs.loc[neuron_mask].copy()
    adata_neuron.obs["div90_broad_class"] = neuron_obs["div90_broad_class"].astype(str).to_numpy()

    class_order = (
        cluster_manifest.loc[cluster_manifest["include_for_siletti_transfer"]]
        .sort_values(args.cluster_col, key=lambda s: s.map(natural_key))["div90_broad_class"]
        .tolist()
    )
    marker_rows = []
    for broad_class in class_order:
        class_mask = adata_neuron.obs["div90_broad_class"].astype(str).to_numpy() == broad_class
        for gene in present_markers:
            values = matrix_column(adata_neuron, gene)[class_mask]
            marker_rows.append(
                {
                    "div90_broad_class": broad_class,
                    "gene": gene,
                    "mean_expr": float(np.nanmean(values)),
                    "pct_expressed": float(np.mean(values > 0) * 100),
                    "n_cells": int(class_mask.sum()),
                }
            )
    marker_summary = pd.DataFrame(marker_rows)
    marker_summary.to_csv(table_dir / "div90_neuron_curated_marker_summary.tsv", sep="\t", index=False)
    marker_mean = marker_summary.pivot(index="gene", columns="div90_broad_class", values="mean_expr").reindex(index=present_markers, columns=class_order)
    marker_z = zscore_rows(marker_mean)
    marker_mean.to_csv(table_dir / "div90_neuron_curated_marker_mean_expression_matrix.tsv", sep="\t")
    marker_z.to_csv(table_dir / "div90_neuron_curated_marker_zscore_matrix.tsv", sep="\t")
    save_heatmap(marker_z, plot_dir / "figure_A_div90_neuron_curated_marker_heatmap_zscore.png", "DIV90 broad neuron classes: curated marker z-score heatmap")

    sc.tl.rank_genes_groups(adata_neuron, groupby="div90_broad_class", method=args.de_method, n_genes=args.max_heatmap_genes)
    de = sc.get.rank_genes_groups_df(adata_neuron, group=None)
    de.to_csv(table_dir / "div90_neuron_broad_class_de_genes.tsv.gz", sep="\t", index=False)
    top_genes = []
    for broad_class, sub in de.groupby("group", sort=False, observed=False):
        for gene in sub.sort_values(["pvals_adj", "scores"], ascending=[True, False])["names"].head(args.top_de_genes_per_class):
            if gene not in top_genes:
                top_genes.append(str(gene))
    top_genes = top_genes[: args.max_heatmap_genes]
    de_rows = []
    for broad_class in class_order:
        class_mask = adata_neuron.obs["div90_broad_class"].astype(str).to_numpy() == broad_class
        for gene in top_genes:
            values = matrix_column(adata_neuron, gene)[class_mask]
            de_rows.append({"div90_broad_class": broad_class, "gene": gene, "mean_expr": float(np.nanmean(values))})
    de_mean = pd.DataFrame(de_rows).pivot(index="gene", columns="div90_broad_class", values="mean_expr").reindex(columns=class_order)
    de_z = zscore_rows(de_mean)
    de_mean.to_csv(table_dir / "div90_neuron_top_de_gene_mean_expression_matrix.tsv", sep="\t")
    de_z.to_csv(table_dir / "div90_neuron_top_de_gene_zscore_matrix.tsv", sep="\t")
    save_heatmap(de_z, plot_dir / "figure_A_div90_neuron_top_de_gene_heatmap_zscore.png", "DIV90 broad neuron classes: top DE gene z-score heatmap")

    run_config = {
        "project_root": str(project_root),
        "h5ad": str(h5ad),
        "outdir": str(outdir),
        "cluster_col": args.cluster_col,
        "cluster_name_col": args.cluster_name_col,
        "excluded_name_patterns": EXCLUDE_NAME_PATTERNS,
        "curated_markers": CURATED_MARKERS,
        "missing_curated_markers": missing_markers,
        "top_de_genes_per_class": args.top_de_genes_per_class,
        "max_heatmap_genes": args.max_heatmap_genes,
        "de_method": args.de_method,
        "seed": args.seed,
    }
    (table_dir / "siletti_div90_prep_run_config.json").write_text(json.dumps(run_config, indent=2, sort_keys=True) + "\n")

    report = [
        "# Siletti DIV90 Transfer Prep",
        "",
        f"Input H5AD: `{h5ad}`",
        f"Output: `{outdir}`",
        "",
        "## Query Manifest",
        "",
        f"- Total DIV90 cells: {adata.n_obs:,}",
        f"- Included neuron-query cells: {int(query_obs['include_for_siletti_transfer'].sum()):,}",
        f"- Excluded cells: {int((~query_obs['include_for_siletti_transfer']).sum()):,}",
        "",
        "Excluded cluster-name patterns:",
        "",
    ]
    report.extend(f"- `{pattern}`" for pattern in EXCLUDE_NAME_PATTERNS)
    report.extend(
        [
            "",
            "## Outputs",
            "",
            "- `tables/div90_neuron_broad_class_manifest.tsv`",
            "- `tables/div90_siletti_query_neuron_cells.tsv`",
            "- `tables/div90_neuron_curated_marker_summary.tsv`",
            "- `tables/div90_neuron_broad_class_de_genes.tsv.gz`",
            "- `plots/figure_A_div90_neuron_curated_marker_heatmap_zscore.png`",
            "- `plots/figure_A_div90_neuron_top_de_gene_heatmap_zscore.png`",
            "",
            "## Notes",
            "",
            "This prep stage does not use Siletti expression yet. It defines the DIV90 query cells/classes and creates figure-A style marker/DE heatmaps.",
        ]
    )
    (report_dir / "siletti_div90_transfer_prep_report.md").write_text("\n".join(report) + "\n")

    print(f"Included query cells: {int(query_obs['include_for_siletti_transfer'].sum()):,}")
    print(report_dir / "siletti_div90_transfer_prep_report.md")


if __name__ == "__main__":
    main()
