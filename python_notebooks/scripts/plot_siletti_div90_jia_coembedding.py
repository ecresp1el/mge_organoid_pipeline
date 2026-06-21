#!/usr/bin/env python
"""Jia-style joint co-embedding of adult Siletti MGE/LLC and DIV90 neurons.

This is a troubleshooting analysis separate from the final adult-reference
projection figures. It asks whether DIV90 neuron classes co-embed with adult
Siletti MGE/LLC inhibitory neurons in a Jia et al. Fig. 4B-like joint manifold.

The workflow is Python-only and uses the same bridge matrices and selected
genes as the fast KNN Siletti transfer:

shared genes -> log-normalize -> selected HVGs -> joint SVD -> joint UMAP
"""

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
import seaborn as sns
import umap
from scipy import sparse
from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import normalize

from siletti_div90_fast_knn_label_transfer import (
    log_normalize,
    read_bridge_counts,
    safe_token,
    sparse_variance,
)


PROJECT_ROOT_DEFAULT = Path("/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder")
DEFAULT_BRIDGE_DIR = (
    PROJECT_ROOT_DEFAULT
    / "results/siletti_2023_whb_reference_label_transfer"
    / "siletti_div90_seurat_bridge_v1/mge_llc/seurat_bridge"
)
DEFAULT_TRANSFER_DIR = (
    PROJECT_ROOT_DEFAULT
    / "results/siletti_2023_whb_reference_label_transfer"
    / "siletti_div90_fast_knn_transfer_full_mge_llc_v1/mge_llc/svd50_k50_ref0_query0"
)
DEFAULT_OUTDIR = (
    PROJECT_ROOT_DEFAULT
    / "results/siletti_2023_whb_reference_label_transfer"
    / "unintegrated_joint_embedding_adult_div90"
)

ADULT_ORDER = [
    "Cortical SST+ LRP neurons",
    "Cortical SST+ nMt neurons",
    "Cortical SST+ Mt neurons",
    "Cortical PV+ basket neurons",
    "Cortical PV+ Chandelier neurons",
    "Subpallial SST+ LRP neurons",
    "Subpallial SST+ neurons",
    "Subpallial PV+ neurons",
    "Subpallial Cholinergic neurons",
]

ADULT_PALETTE = {
    "Cortical SST+ LRP neurons": "#4c78a8",
    "Cortical SST+ nMt neurons": "#72b7b2",
    "Cortical SST+ Mt neurons": "#54a24b",
    "Cortical PV+ basket neurons": "#e45756",
    "Cortical PV+ Chandelier neurons": "#ff9da6",
    "Subpallial SST+ LRP neurons": "#f58518",
    "Subpallial SST+ neurons": "#b279a2",
    "Subpallial PV+ neurons": "#9d755d",
    "Subpallial Cholinergic neurons": "#bab0ac",
    "Unassigned": "#9a9a9a",
}

SIDE_PALETTE = {
    "Cortical": "#6c8ebf",
    "Subpallial": "#d89555",
    "Unassigned": "#b8b8b8",
}


def adult_side(label: str) -> str:
    label = str(label)
    if label.startswith("Cortical"):
        return "Cortical"
    if label.startswith("Subpallial"):
        return "Subpallial"
    return "Unassigned"


def read_tsv_gz(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt") as handle:
        return pd.read_csv(handle, sep="\t")


def write_tsv_gz(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        df.to_csv(handle, sep="\t", index=False)


def feature_indices(genes: list[str], transfer_dir: Path, nfeatures: int, x_ref: sparse.csr_matrix, x_query: sparse.csr_matrix) -> tuple[np.ndarray, list[str], str]:
    features_path = transfer_dir / "fast_knn/selected_transfer_features.tsv"
    gene_to_idx = {gene: i for i, gene in enumerate(genes)}
    if features_path.exists():
        features = pd.read_csv(features_path, sep="\t")["gene"].astype(str).tolist()
        idx = np.array([gene_to_idx[g] for g in features if g in gene_to_idx], dtype=int)
        return idx, [genes[i] for i in idx], str(features_path)

    combined = sparse.vstack([log_normalize(x_ref), log_normalize(x_query)], format="csr")
    var = sparse_variance(combined)
    idx = np.argsort(var)[::-1][: min(nfeatures, x_ref.shape[1])]
    idx.sort()
    return idx, [genes[i] for i in idx], "computed_from_joint_log_normalized_variance"


def load_joint_inputs(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray, np.ndarray, np.ndarray, dict]:
    x_ref, genes_ref, meta_ref = read_bridge_counts("reference", args.bridge_dir)
    x_query, genes_query, meta_query = read_bridge_counts("query", args.bridge_dir)
    if genes_ref != genes_query:
        raise ValueError("Reference and query bridge gene orders differ.")
    if args.label_column not in meta_ref.columns:
        raise ValueError(f"Reference metadata missing label column: {args.label_column}")
    if args.query_class_col not in meta_query.columns:
        raise ValueError(f"Query metadata missing DIV90 class column: {args.query_class_col}")

    meta_ref[args.label_column] = meta_ref[args.label_column].fillna("Unassigned").astype(str)
    if args.exclude_label and args.exclude_label.upper() != "NONE":
        keep = meta_ref[args.label_column] != args.exclude_label
        x_ref = x_ref[keep.to_numpy(), :]
        meta_ref = meta_ref.loc[keep].reset_index(drop=True)

    idx, selected_genes, feature_source = feature_indices(genes_ref, args.transfer_dir, args.nfeatures, x_ref, x_query)
    x_ref = log_normalize(x_ref)[:, idx]
    x_query = log_normalize(x_query)[:, idx]
    joint = sparse.vstack([x_ref, x_query], format="csr")

    n_components = min(args.n_components, joint.shape[0] - 1, joint.shape[1] - 1)
    svd = TruncatedSVD(n_components=n_components, random_state=args.seed)
    joint_pcs = normalize(svd.fit_transform(joint))
    ref_pcs = joint_pcs[: x_ref.shape[0], :]
    query_pcs = joint_pcs[x_ref.shape[0] :, :]

    reducer = umap.UMAP(
        n_neighbors=args.umap_neighbors,
        min_dist=args.umap_min_dist,
        metric=args.umap_metric,
        random_state=args.seed,
        low_memory=True,
    )
    coords = reducer.fit_transform(joint_pcs)

    meta_ref = meta_ref.copy()
    meta_ref["dataset"] = "Siletti adult reference"
    meta_ref["adult_subtype"] = meta_ref[args.label_column].astype(str)
    meta_ref["adult_side"] = meta_ref["adult_subtype"].map(adult_side)
    meta_ref["umap_1"] = coords[: x_ref.shape[0], 0]
    meta_ref["umap_2"] = coords[: x_ref.shape[0], 1]

    meta_query = meta_query.copy()
    meta_query["dataset"] = "DIV90 query"
    meta_query["div90_class"] = meta_query[args.query_class_col].astype(str)
    meta_query["umap_1"] = coords[x_ref.shape[0] :, 0]
    meta_query["umap_2"] = coords[x_ref.shape[0] :, 1]

    diag = {
        "n_reference_cells": int(meta_ref.shape[0]),
        "n_query_cells": int(meta_query.shape[0]),
        "n_shared_genes_bridge": int(len(genes_ref)),
        "n_selected_features": int(len(selected_genes)),
        "selected_feature_source": feature_source,
        "n_components": int(n_components),
        "svd_explained_variance_ratio_sum": float(svd.explained_variance_ratio_.sum()),
        "umap_neighbors": int(args.umap_neighbors),
        "umap_min_dist": float(args.umap_min_dist),
        "umap_metric": args.umap_metric,
        "embedding_fit_description": "Joint SVD and joint UMAP fit on adult Siletti reference plus DIV90 query cells.",
    }
    return meta_ref, meta_query, ref_pcs, query_pcs, coords, diag


def plot_adult_side_with_div90_classes(ref: pd.DataFrame, query: pd.DataFrame, plots_dir: Path) -> None:
    div90_classes = sorted(query["div90_class"].dropna().astype(str).unique())
    div90_palette = dict(zip(div90_classes, sns.color_palette("tab20", n_colors=max(1, len(div90_classes))).as_hex()))

    fig, ax = plt.subplots(figsize=(8.5, 7.2))
    for side in ["Cortical", "Subpallial", "Unassigned"]:
        sub = ref.loc[ref["adult_side"] == side]
        if sub.empty:
            continue
        ax.scatter(
            sub["umap_1"],
            sub["umap_2"],
            s=2,
            c=SIDE_PALETTE.get(side, "#b8b8b8"),
            alpha=0.23,
            linewidths=0,
            rasterized=True,
            label=f"Adult {side}",
        )
    for label in div90_classes:
        sub = query.loc[query["div90_class"] == label]
        ax.scatter(
            sub["umap_1"],
            sub["umap_2"],
            s=7,
            c=div90_palette[label],
            alpha=0.82,
            linewidths=0,
            rasterized=True,
            label=f"DIV90 {label}",
        )
    ax.set_xlabel("Joint UMAP 1")
    ax.set_ylabel("Joint UMAP 2")
    ax.set_title("Jia-style joint embedding: adult side background, DIV90 classes overlaid")
    ax.legend(markerscale=3, fontsize=6.5, frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.savefig(plots_dir / "joint_umap_adult_side_div90_classes_overlay.png", dpi=300, bbox_inches="tight")
    fig.savefig(plots_dir / "joint_umap_adult_side_div90_classes_overlay.pdf", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.5, 7.2))
    ax.scatter(ref["umap_1"], ref["umap_2"], s=1.8, c="#d0d0d0", alpha=0.45, linewidths=0, rasterized=True, label="Adult reference")
    for label in div90_classes:
        sub = query.loc[query["div90_class"] == label]
        ax.scatter(sub["umap_1"], sub["umap_2"], s=7, c=div90_palette[label], alpha=0.82, linewidths=0, rasterized=True, label=label)
    ax.set_xlabel("Joint UMAP 1")
    ax.set_ylabel("Joint UMAP 2")
    ax.set_title("Jia-style joint embedding: adult gray, DIV90 classes overlaid")
    ax.legend(markerscale=3, fontsize=6.5, frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.savefig(plots_dir / "joint_umap_adult_gray_div90_classes_overlay.png", dpi=300, bbox_inches="tight")
    fig.savefig(plots_dir / "joint_umap_adult_gray_div90_classes_overlay.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_adult_subtypes_with_div90_black(ref: pd.DataFrame, query: pd.DataFrame, plots_dir: Path) -> None:
    labels = [x for x in ADULT_ORDER if x in set(ref["adult_subtype"])]
    labels += sorted(set(ref["adult_subtype"]) - set(labels))
    fig, ax = plt.subplots(figsize=(8.5, 7.2))
    for label in labels:
        sub = ref.loc[ref["adult_subtype"] == label]
        ax.scatter(
            sub["umap_1"],
            sub["umap_2"],
            s=2.2,
            c=ADULT_PALETTE.get(label, "#777777"),
            alpha=0.52,
            linewidths=0,
            rasterized=True,
            label=label,
        )
    ax.scatter(query["umap_1"], query["umap_2"], s=5, c="black", alpha=0.36, linewidths=0, rasterized=True, label="DIV90 query")
    ax.set_xlabel("Joint UMAP 1")
    ax.set_ylabel("Joint UMAP 2")
    ax.set_title("Jia-style joint embedding: adult subtypes with DIV90 overlay")
    ax.legend(markerscale=4, fontsize=7, frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.savefig(plots_dir / "joint_umap_adult_subtypes_div90_black_overlay.png", dpi=300, bbox_inches="tight")
    fig.savefig(plots_dir / "joint_umap_adult_subtypes_div90_black_overlay.pdf", bbox_inches="tight")
    plt.close(fig)


def quantify_coembedding(
    ref: pd.DataFrame,
    query: pd.DataFrame,
    ref_pcs: np.ndarray,
    query_pcs: np.ndarray,
    k: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    k = min(k, ref_pcs.shape[0])
    nn = NearestNeighbors(n_neighbors=k, metric="cosine", algorithm="brute", n_jobs=-1)
    nn.fit(ref_pcs)
    distances, indices = nn.kneighbors(query_pcs)
    neighbor_side = ref["adult_side"].to_numpy()[indices]
    neighbor_subtype = ref["adult_subtype"].to_numpy()[indices]

    rows = []
    for i in range(query.shape[0]):
        row = {
            "seurat_cell_id": query["seurat_cell_id"].iloc[i],
            "div90_class": query["div90_class"].iloc[i],
            "nearest_adult_cosine_distance": float(distances[i, 0]),
            "mean_adult_cosine_distance_k": float(distances[i].mean()),
            "fraction_adult_neighbors_cortical": float(np.mean(neighbor_side[i] == "Cortical")),
            "fraction_adult_neighbors_subpallial": float(np.mean(neighbor_side[i] == "Subpallial")),
            "top_neighbor_adult_side": str(pd.Series(neighbor_side[i]).value_counts().idxmax()),
            "top_neighbor_adult_subtype": str(pd.Series(neighbor_subtype[i]).value_counts().idxmax()),
        }
        rows.append(row)
    per_cell = pd.DataFrame(rows)

    edge_rows = []
    for i, qrow in query.reset_index(drop=True).iterrows():
        counts = pd.Series(neighbor_subtype[i]).value_counts()
        for subtype, n in counts.items():
            edge_rows.append({
                "seurat_cell_id": qrow["seurat_cell_id"],
                "div90_class": qrow["div90_class"],
                "adult_subtype": subtype,
                "n_neighbors": int(n),
                "fraction_neighbors": float(n / k),
            })
    per_cell_edges = pd.DataFrame(edge_rows)

    # Adult-territory thresholds: compare query-to-adult nearest distances with
    # adult-to-adult nearest distances. UMAP thresholds are descriptive only.
    adult_nn = NearestNeighbors(n_neighbors=2, metric="cosine", algorithm="brute", n_jobs=-1)
    adult_nn.fit(ref_pcs)
    adult_svd_dists, _ = adult_nn.kneighbors(ref_pcs)
    svd_threshold = float(np.quantile(adult_svd_dists[:, 1], 0.95))

    ref_umap = ref[["umap_1", "umap_2"]].to_numpy()
    query_umap = query[["umap_1", "umap_2"]].to_numpy()
    umap_nn = NearestNeighbors(n_neighbors=1, metric="euclidean", algorithm="auto")
    umap_nn.fit(ref_umap)
    query_umap_dist, _ = umap_nn.kneighbors(query_umap)
    adult_umap_nn = NearestNeighbors(n_neighbors=2, metric="euclidean", algorithm="auto")
    adult_umap_nn.fit(ref_umap)
    adult_umap_dist, _ = adult_umap_nn.kneighbors(ref_umap)
    umap_threshold = float(np.quantile(adult_umap_dist[:, 1], 0.95))

    per_cell["nearest_adult_umap_distance"] = query_umap_dist[:, 0]
    per_cell["within_adult_svd_territory_95pct"] = per_cell["nearest_adult_cosine_distance"] <= svd_threshold
    per_cell["within_adult_umap_territory_95pct"] = per_cell["nearest_adult_umap_distance"] <= umap_threshold

    summary_rows = []
    for div90_class, sub in per_cell.groupby("div90_class", sort=True):
        subtype_counts = (
            per_cell_edges.loc[per_cell_edges["div90_class"] == div90_class]
            .groupby("adult_subtype")["n_neighbors"]
            .sum()
            .sort_values(ascending=False)
        )
        top_subtype = str(subtype_counts.index[0]) if not subtype_counts.empty else "NA"
        frac_svd = float(sub["within_adult_svd_territory_95pct"].mean())
        frac_umap = float(sub["within_adult_umap_territory_95pct"].mean())
        if frac_svd >= 0.5 and frac_umap >= 0.5:
            status = "embedded_within_adult_territory"
        elif frac_svd < 0.2 and frac_umap < 0.2:
            status = "separate_island"
        else:
            status = "partial_or_boundary_embedding"
        summary_rows.append({
            "div90_class": div90_class,
            "n_cells": int(sub.shape[0]),
            "mean_fraction_adult_neighbors_cortical": float(sub["fraction_adult_neighbors_cortical"].mean()),
            "mean_fraction_adult_neighbors_subpallial": float(sub["fraction_adult_neighbors_subpallial"].mean()),
            "median_nearest_adult_cosine_distance": float(sub["nearest_adult_cosine_distance"].median()),
            "median_nearest_adult_umap_distance": float(sub["nearest_adult_umap_distance"].median()),
            "fraction_within_adult_svd_territory_95pct": frac_svd,
            "fraction_within_adult_umap_territory_95pct": frac_umap,
            "top_adult_subtype_by_neighbor_votes": top_subtype,
            "embedding_status_call": status,
            "adult_svd_95pct_nearest_neighbor_threshold": svd_threshold,
            "adult_umap_95pct_nearest_neighbor_threshold": umap_threshold,
        })
    summary = pd.DataFrame(summary_rows)
    return per_cell, per_cell_edges, summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bridge-dir", type=Path, default=DEFAULT_BRIDGE_DIR)
    parser.add_argument("--transfer-dir", type=Path, default=DEFAULT_TRANSFER_DIR)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--label-column", default="candidate_jia_group")
    parser.add_argument("--exclude-label", default="Excluded / not assigned to Jia-style 9 groups")
    parser.add_argument("--query-class-col", default="div90_broad_class")
    parser.add_argument("--nfeatures", type=int, default=3000)
    parser.add_argument("--n-components", type=int, default=50)
    parser.add_argument("--neighbor-k", type=int, default=50)
    parser.add_argument("--umap-neighbors", type=int, default=30)
    parser.add_argument("--umap-min-dist", type=float, default=0.3)
    parser.add_argument("--umap-metric", default="cosine")
    parser.add_argument("--seed", type=int, default=0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    plots_dir = args.outdir / "plots"
    tables_dir = args.outdir / "tables"
    reports_dir = args.outdir / "reports"
    for path in [plots_dir, tables_dir, reports_dir]:
        path.mkdir(parents=True, exist_ok=True)

    ref, query, ref_pcs, query_pcs, _coords, diag = load_joint_inputs(args)
    plot_adult_side_with_div90_classes(ref, query, plots_dir)
    plot_adult_subtypes_with_div90_black(ref, query, plots_dir)

    per_cell, per_cell_edges, summary = quantify_coembedding(ref, query, ref_pcs, query_pcs, args.neighbor_k)
    summary.to_csv(tables_dir / "div90_class_adult_neighbor_coembedding_summary.tsv", sep="\t", index=False)
    per_cell_edges.to_csv(tables_dir / "div90_cell_adult_neighbor_subtype_composition.tsv", sep="\t", index=False)
    write_tsv_gz(per_cell, tables_dir / "div90_cell_adult_neighbor_coembedding_metrics.tsv.gz")
    write_tsv_gz(query, tables_dir / "div90_query_cells_joint_umap.tsv.gz")
    write_tsv_gz(ref, tables_dir / "siletti_reference_cells_joint_umap.tsv.gz")

    config = {
        "bridge_dir": str(args.bridge_dir),
        "transfer_dir": str(args.transfer_dir),
        "outdir": str(args.outdir),
        "label_column": args.label_column,
        "exclude_label": args.exclude_label,
        "query_class_col": args.query_class_col,
        "nfeatures_requested": args.nfeatures,
        "n_components_requested": args.n_components,
        "neighbor_k": args.neighbor_k,
        "seed": args.seed,
        **diag,
    }
    (tables_dir / "coembedding_config.json").write_text(json.dumps(config, indent=2, sort_keys=True) + "\n")

    readme = [
        "# Siletti DIV90 Jia-Style Co-Embedding",
        "",
        "This analysis is separate from the adult-reference projection workflow.",
        "It fits a joint SVD and joint UMAP on adult Siletti MGE/LLC reference cells plus DIV90 neuron-lineage cells.",
        "",
        "Primary question: do DIV90 classes co-embed with adult Siletti MGE/LLC inhibitory neurons,",
        "or do they form a separate island in the joint manifold?",
        "",
        "Generated diagnostics summarize adult-neighbor composition and nearest-adult distances for each DIV90 class.",
    ]
    (reports_dir / "README_siletti_div90_jia_style_coembedding.md").write_text("\n".join(readme) + "\n")

    print(json.dumps(config, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
