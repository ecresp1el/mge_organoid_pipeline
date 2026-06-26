#!/usr/bin/env python3
"""Render all-Siletti-supercluster DIV90 transfer plots."""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_siletti_div90_jia_figure import (
    adult_label_order,
    draw_river,
    find_one,
    prepare_embedding,
    write_tsv_gz,
)


def read_tsv_gz(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt") as handle:
        return pd.read_csv(handle, sep="\t")


def color_map(labels: list[str]) -> dict[str, str]:
    cmap_names = ["tab20", "tab20b", "tab20c", "Set3"]
    colors: list[str] = []
    for name in cmap_names:
        cmap = plt.get_cmap(name)
        colors.extend(matplotlib.colors.to_hex(cmap(i)) for i in range(cmap.N))
    return {label: colors[i % len(colors)] for i, label in enumerate(labels)}


def attach_predictions(query_umap: pd.DataFrame, transfer_dir: Path) -> pd.DataFrame:
    obs_path = find_one("tables/*_query_obs_with_predictions.tsv.gz", transfer_dir)
    pred = read_tsv_gz(obs_path)
    merged = query_umap.merge(
        pred[["cell_id", "predicted.id", "prediction.score.max"]],
        left_on="seurat_cell_id",
        right_on="cell_id",
        how="left",
    )
    merged["predicted.id"] = merged["predicted.id"].fillna("Unassigned").astype(str)
    merged["assignment_status"] = np.where(merged["predicted.id"].eq("Unassigned"), "unassigned", "assigned")
    return merged


def plot_reference_overlay(ref: pd.DataFrame, query: pd.DataFrame, label_column: str, plots_dir: Path) -> None:
    labels = sorted(ref[label_column].dropna().astype(str).unique())
    colors = color_map(labels)
    fig, ax = plt.subplots(figsize=(10, 8))
    for label in labels:
        sub = ref.loc[ref[label_column].astype(str).eq(label)]
        ax.scatter(
            sub["umap_1"],
            sub["umap_2"],
            s=2.0,
            c=colors[label],
            label=label,
            alpha=0.45,
            linewidths=0,
            rasterized=True,
        )
    ax.scatter(query["umap_1"], query["umap_2"], s=4, c="#111111", alpha=0.35, linewidths=0, rasterized=True, label="DIV90 query")
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title("Siletti all-supercluster reference with DIV90 overlay")
    ax.legend(markerscale=4, fontsize=6, frameon=False, loc="center left", bbox_to_anchor=(1.02, 0.5), ncol=1)
    fig.savefig(plots_dir / "figure_B_all_supercluster_reference_div90_overlay.png", dpi=300, bbox_inches="tight")
    fig.savefig(plots_dir / "figure_B_all_supercluster_reference_div90_overlay.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_query_predictions(ref: pd.DataFrame, query: pd.DataFrame, plots_dir: Path) -> None:
    labels = sorted(query["predicted.id"].dropna().astype(str).unique())
    colors = color_map(labels)
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.scatter(ref["umap_1"], ref["umap_2"], s=1.0, c="#d7d7d7", alpha=0.35, linewidths=0, rasterized=True)
    for label in labels:
        sub = query.loc[query["predicted.id"].astype(str).eq(label)]
        ax.scatter(
            sub["umap_1"],
            sub["umap_2"],
            s=7,
            c=colors[label],
            label=label,
            alpha=0.78,
            linewidths=0,
            rasterized=True,
        )
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title("DIV90 predicted Siletti superclusters in adult reference space")
    ax.legend(markerscale=3, fontsize=7, frameon=False, loc="center left", bbox_to_anchor=(1.02, 0.5))
    fig.savefig(plots_dir / "figure_B_div90_predicted_siletti_superclusters.png", dpi=300, bbox_inches="tight")
    fig.savefig(plots_dir / "figure_B_div90_predicted_siletti_superclusters.pdf", bbox_inches="tight")
    plt.close(fig)


def make_edges(query: pd.DataFrame, source_col: str) -> pd.DataFrame:
    edges = (
        query.groupby([source_col, "predicted.id"], dropna=False)
        .size()
        .reset_index(name="n_cells")
        .rename(columns={source_col: "div90_class", "predicted.id": "adult_subtype"})
    )
    return edges


def sample_proportions(query: pd.DataFrame, sample_col: str) -> pd.DataFrame:
    counts = (
        query.groupby([sample_col, "predicted.id"], dropna=False)
        .size()
        .reset_index(name="n_cells")
        .rename(columns={sample_col: "sample", "predicted.id": "predicted_supercluster"})
    )
    counts["fraction"] = counts["n_cells"] / counts.groupby("sample")["n_cells"].transform("sum")
    return counts


def draw_sample_barplot(props: pd.DataFrame, plots_dir: Path) -> None:
    pivot = props.pivot_table(index="sample", columns="predicted_supercluster", values="fraction", fill_value=0)
    cols = sorted(pivot.columns.astype(str))
    pivot = pivot[cols]
    colors = color_map(cols)
    fig, ax = plt.subplots(figsize=(max(9, 0.75 * pivot.shape[0]), 6))
    bottom = np.zeros(pivot.shape[0])
    x = np.arange(pivot.shape[0])
    for col in cols:
        vals = pivot[col].to_numpy()
        ax.bar(x, vals, bottom=bottom, label=col, color=colors[col], edgecolor="white", linewidth=0.35)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(pivot.index.astype(str), rotation=45, ha="right")
    ax.set_ylim(0, 1)
    ax.set_ylabel("Fraction of DIV90 cells")
    ax.set_title("Predicted Siletti superclusters across DIV90 samples")
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7)
    fig.savefig(plots_dir / "figure_D_sample_predicted_siletti_supercluster_proportions.png", dpi=300, bbox_inches="tight")
    fig.savefig(plots_dir / "figure_D_sample_predicted_siletti_supercluster_proportions.pdf", bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bridge-dir", required=True, type=Path)
    parser.add_argument("--transfer-dir", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--label-column", default="source_supercluster")
    parser.add_argument("--source-class-col", default="div90_broad_class")
    parser.add_argument("--sample-col", default="orig.ident")
    parser.add_argument("--nfeatures", type=int, default=3000)
    parser.add_argument("--n-components", type=int, default=50)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--umap-mode", choices=["reference_project", "combined"], default="reference_project")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    plots_dir = args.outdir / "plots"
    tables_dir = args.outdir / "tables"
    reports_dir = args.outdir / "reports"
    for directory in [plots_dir, tables_dir, reports_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    ref, query_umap, diag = prepare_embedding(
        args.bridge_dir,
        args.transfer_dir,
        args.label_column,
        [],
        args.nfeatures,
        args.n_components,
        args.seed,
        args.umap_mode,
    )
    query = attach_predictions(query_umap, args.transfer_dir)

    plot_reference_overlay(ref, query, args.label_column, plots_dir)
    plot_query_predictions(ref, query, plots_dir)

    edges = make_edges(query, args.source_class_col)
    edges.to_csv(tables_dir / "figure_C_div90_class_to_siletti_supercluster_edges.tsv", sep="\t", index=False)
    draw_river(edges, plots_dir / "figure_C_div90_class_to_siletti_supercluster_river")

    props = sample_proportions(query, args.sample_col)
    props.to_csv(tables_dir / "figure_D_sample_predicted_siletti_supercluster_proportions.tsv", sep="\t", index=False)
    draw_sample_barplot(props, plots_dir)

    query.to_csv(tables_dir / "div90_query_cells_with_all_siletti_supercluster_assignments.tsv.gz", sep="\t", index=False, compression="gzip")
    ref[[ "seurat_cell_id", args.label_column, "cell_type", "umap_1", "umap_2" ]].to_csv(
        tables_dir / "siletti_reference_cells_with_all_supercluster_umap.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    config = {
        "bridge_dir": str(args.bridge_dir),
        "transfer_dir": str(args.transfer_dir),
        "outdir": str(args.outdir),
        "label_column": args.label_column,
        "source_class_col": args.source_class_col,
        "sample_col": args.sample_col,
        "nfeatures_requested": args.nfeatures,
        "n_components_requested": args.n_components,
        "seed": args.seed,
        "umap_mode": args.umap_mode,
        **diag,
    }
    (tables_dir / "plot_config.json").write_text(json.dumps(config, indent=2, sort_keys=True) + "\n")
    (reports_dir / "README_all_supercluster_figure.md").write_text(
        "# Siletti All-Supercluster DIV90 Transfer Figure\n\n"
        "These plots transfer Siletti WHB supercluster labels onto DIV90 query cells.\n"
        "The reference is an all-supercluster downsampled bridge for plotting speed.\n"
    )
    print(json.dumps(config, indent=2, sort_keys=True))
    for path in sorted(args.outdir.rglob("*")):
        if path.is_file():
            print(path)


if __name__ == "__main__":
    main()
