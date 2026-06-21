#!/usr/bin/env python
"""Harmony integration sensitivity for Siletti adult MGE/LLC + DIV90 neurons.

This workflow is separate from:

* the adult-reference projection figures, and
* the unintegrated joint embedding diagnostic.

It is intended to mimic Jia et al.'s fetal+adult integrated analysis more
closely by explicitly correcting the joint SVD/PCA coordinates for source
dataset: ``adult_reference`` versus ``DIV90_query``.
"""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import harmonypy as hm
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import umap
from scipy import sparse
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import normalize

from plot_siletti_div90_jia_coembedding import (
    DEFAULT_BRIDGE_DIR,
    DEFAULT_TRANSFER_DIR,
    PROJECT_ROOT_DEFAULT,
    adult_side,
    feature_indices,
    plot_adult_side_with_div90_classes,
    plot_adult_subtypes_with_div90_black,
    quantify_coembedding,
    write_tsv_gz,
)
from siletti_div90_fast_knn_label_transfer import log_normalize, read_bridge_counts


DEFAULT_OUTDIR = (
    PROJECT_ROOT_DEFAULT
    / "results/siletti_2023_whb_reference_label_transfer"
    / "harmony_integrated_adult_div90"
)
DEFAULT_UNINTEGRATED_DIR = (
    PROJECT_ROOT_DEFAULT
    / "results/siletti_2023_whb_reference_label_transfer"
    / "unintegrated_joint_embedding_adult_div90"
)


def load_svd_inputs(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray, np.ndarray, dict]:
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
    joint_svd = normalize(svd.fit_transform(joint))

    meta_ref = meta_ref.copy()
    meta_ref["dataset"] = "Siletti adult reference"
    meta_ref["source_label"] = "adult_reference"
    meta_ref["adult_subtype"] = meta_ref[args.label_column].astype(str)
    meta_ref["adult_side"] = meta_ref["adult_subtype"].map(adult_side)

    meta_query = meta_query.copy()
    meta_query["dataset"] = "DIV90 query"
    meta_query["source_label"] = "DIV90_query"
    meta_query["div90_class"] = meta_query[args.query_class_col].astype(str)

    metadata = pd.concat(
        [
            pd.DataFrame({"source_label": meta_ref["source_label"].to_numpy()}),
            pd.DataFrame({"source_label": meta_query["source_label"].to_numpy()}),
        ],
        ignore_index=True,
    )

    diag = {
        "n_reference_cells": int(meta_ref.shape[0]),
        "n_query_cells": int(meta_query.shape[0]),
        "n_shared_genes_bridge": int(len(genes_ref)),
        "n_selected_features": int(len(selected_genes)),
        "selected_feature_source": feature_source,
        "n_components": int(n_components),
        "svd_explained_variance_ratio_sum": float(svd.explained_variance_ratio_.sum()),
    }
    return meta_ref, meta_query, joint_svd, metadata, diag


def run_harmony(joint_svd: np.ndarray, metadata: pd.DataFrame, args: argparse.Namespace) -> np.ndarray:
    harmony = hm.run_harmony(
        joint_svd,
        metadata,
        vars_use=["source_label"],
        theta=args.harmony_theta,
        lamb=args.harmony_lambda,
        max_iter_harmony=args.harmony_max_iter,
        random_state=args.seed,
        verbose=True,
    )
    corrected = np.asarray(harmony.Z_corr)
    if corrected.shape[0] != joint_svd.shape[0] and corrected.shape[1] == joint_svd.shape[0]:
        corrected = corrected.T
    if corrected.shape[0] != joint_svd.shape[0]:
        raise ValueError(
            "Harmony output has unexpected shape "
            f"{corrected.shape}; expected {joint_svd.shape[0]} cells."
        )
    return normalize(corrected)


def add_integrated_umap(
    ref: pd.DataFrame,
    query: pd.DataFrame,
    corrected: np.ndarray,
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    reducer = umap.UMAP(
        n_neighbors=args.umap_neighbors,
        min_dist=args.umap_min_dist,
        metric=args.umap_metric,
        random_state=args.seed,
        low_memory=True,
    )
    coords = reducer.fit_transform(corrected)
    ref = ref.copy()
    query = query.copy()
    ref["umap_1"] = coords[: ref.shape[0], 0]
    ref["umap_2"] = coords[: ref.shape[0], 1]
    query["umap_1"] = coords[ref.shape[0] :, 0]
    query["umap_2"] = coords[ref.shape[0] :, 1]
    return ref, query


def plot_by_source(ref: pd.DataFrame, query: pd.DataFrame, plots_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 7.0))
    ax.scatter(
        ref["umap_1"],
        ref["umap_2"],
        s=2.2,
        c="#9b9b9b",
        alpha=0.42,
        linewidths=0,
        rasterized=True,
        label="adult_reference",
    )
    ax.scatter(
        query["umap_1"],
        query["umap_2"],
        s=5.5,
        c="#111111",
        alpha=0.44,
        linewidths=0,
        rasterized=True,
        label="DIV90_query",
    )
    ax.set_xlabel("Harmony UMAP 1")
    ax.set_ylabel("Harmony UMAP 2")
    ax.set_title("Harmony-integrated adult + DIV90 embedding by source")
    ax.legend(frameon=False, markerscale=4, loc="best")
    fig.savefig(plots_dir / "harmony_umap_by_source.png", dpi=300, bbox_inches="tight")
    fig.savefig(plots_dir / "harmony_umap_by_source.pdf", bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bridge-dir", type=Path, default=DEFAULT_BRIDGE_DIR)
    parser.add_argument("--transfer-dir", type=Path, default=DEFAULT_TRANSFER_DIR)
    parser.add_argument("--unintegrated-dir", type=Path, default=DEFAULT_UNINTEGRATED_DIR)
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
    parser.add_argument("--harmony-theta", type=float, default=2.0)
    parser.add_argument("--harmony-lambda", type=float, default=1.0)
    parser.add_argument("--harmony-max-iter", type=int, default=20)
    parser.add_argument("--seed", type=int, default=0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    plots_dir = args.outdir / "plots"
    tables_dir = args.outdir / "tables"
    reports_dir = args.outdir / "reports"
    for path in [plots_dir, tables_dir, reports_dir]:
        path.mkdir(parents=True, exist_ok=True)

    ref, query, joint_svd, harmony_meta, diag = load_svd_inputs(args)
    corrected = run_harmony(joint_svd, harmony_meta, args)
    ref, query = add_integrated_umap(ref, query, corrected, args)

    ref_latent = corrected[: ref.shape[0], :]
    query_latent = corrected[ref.shape[0] :, :]

    plot_by_source(ref, query, plots_dir)
    plot_adult_side_with_div90_classes(ref, query, plots_dir)
    plot_adult_subtypes_with_div90_black(ref, query, plots_dir)

    per_cell, per_cell_edges, summary = quantify_coembedding(ref, query, ref_latent, query_latent, args.neighbor_k)
    summary.to_csv(tables_dir / "div90_class_adult_neighbor_harmony_summary.tsv", sep="\t", index=False)
    per_cell_edges.to_csv(tables_dir / "div90_cell_adult_neighbor_subtype_composition_harmony.tsv", sep="\t", index=False)
    write_tsv_gz(per_cell, tables_dir / "div90_cell_adult_neighbor_harmony_metrics.tsv.gz")
    write_tsv_gz(query, tables_dir / "div90_query_cells_harmony_umap.tsv.gz")
    write_tsv_gz(ref, tables_dir / "siletti_reference_cells_harmony_umap.tsv.gz")

    config = {
        "bridge_dir": str(args.bridge_dir),
        "transfer_dir": str(args.transfer_dir),
        "unintegrated_dir": str(args.unintegrated_dir),
        "outdir": str(args.outdir),
        "integration_method": "Harmony",
        "batch_variable": "source_label",
        "source_labels": ["adult_reference", "DIV90_query"],
        "label_column": args.label_column,
        "exclude_label": args.exclude_label,
        "query_class_col": args.query_class_col,
        "nfeatures_requested": args.nfeatures,
        "n_components_requested": args.n_components,
        "neighbor_k": args.neighbor_k,
        "harmony_theta": args.harmony_theta,
        "harmony_lambda": args.harmony_lambda,
        "harmony_max_iter": args.harmony_max_iter,
        "umap_neighbors": args.umap_neighbors,
        "umap_min_dist": args.umap_min_dist,
        "umap_metric": args.umap_metric,
        "seed": args.seed,
        "embedding_fit_description": "Joint SVD followed by Harmony correction for source_label; UMAP fit on Harmony-corrected latent space.",
        **diag,
    }
    (tables_dir / "harmony_integration_config.json").write_text(json.dumps(config, indent=2, sort_keys=True) + "\n")

    readme = [
        "# Siletti DIV90 Integration Sensitivity",
        "",
        "This workflow keeps the unintegrated joint embedding as a diagnostic and adds Harmony integration.",
        "",
        f"Unintegrated diagnostic directory: `{args.unintegrated_dir}`",
        "",
        "Integrated workflow:",
        "shared genes -> same 3,000 transfer genes -> log-normalize -> joint SVD -> Harmony(source_label) -> UMAP",
        "",
        "The nearest-adult diagnostics are computed in Harmony-corrected latent space.",
    ]
    (reports_dir / "README_harmony_integrated_adult_div90.md").write_text("\n".join(readme) + "\n")

    print(json.dumps(config, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
