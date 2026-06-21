#!/usr/bin/env python
"""Forced Jia-style BBKNN integration for Siletti adult MGE/LLC + DIV90 cells.

This visualization workflow intentionally builds a batch-balanced graph after
Harmony correction:

shared genes -> selected 3,000 transfer genes -> joint SVD -> Harmony(source)
-> BBKNN(source) -> UMAP from the BBKNN graph.
"""

from __future__ import annotations

import argparse
import gzip
import importlib.util
import json
from pathlib import Path

import harmonypy as hm
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from anndata import AnnData
from scipy import sparse
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import normalize

from plot_siletti_div90_jia_coembedding import (
    ADULT_ORDER,
    ADULT_PALETTE,
    DEFAULT_BRIDGE_DIR,
    DEFAULT_TRANSFER_DIR,
    PROJECT_ROOT_DEFAULT,
    SIDE_PALETTE,
    adult_side,
    feature_indices,
)
from siletti_div90_fast_knn_label_transfer import log_normalize, read_bridge_counts


DEFAULT_OUTDIR = (
    PROJECT_ROOT_DEFAULT
    / "results/siletti_2023_whb_reference_label_transfer"
    / "forced_jia_style_integration_adult_div90"
)


def write_tsv_gz(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        df.to_csv(handle, sep="\t", index=False)


def parse_trim(value: str) -> int | None:
    if str(value).lower() in {"none", "null", "na"}:
        return None
    return int(value)


def div90_palette(query: pd.DataFrame) -> dict[str, str]:
    classes = sorted(query["div90_class"].dropna().astype(str).unique())
    return dict(zip(classes, sns.color_palette("tab20", n_colors=max(1, len(classes))).as_hex()))


def load_svd_inputs(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray, pd.DataFrame, dict]:
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
    x_svd = normalize(svd.fit_transform(joint))

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
    return meta_ref, meta_query, x_svd, metadata, diag


def run_harmony(x_svd: np.ndarray, metadata: pd.DataFrame, args: argparse.Namespace) -> np.ndarray:
    ho = hm.run_harmony(
        data_mat=x_svd,
        meta_data=metadata,
        vars_use=["source_label"],
        theta=args.harmony_theta,
        lamb=args.harmony_lambda,
        max_iter_harmony=args.harmony_max_iter,
        max_iter_kmeans=args.harmony_max_iter_kmeans,
        epsilon_cluster=args.harmony_epsilon_cluster,
        epsilon_harmony=args.harmony_epsilon_harmony,
        random_state=args.seed,
        verbose=True,
    )
    z = np.asarray(ho.Z_corr)
    if z.shape[0] == x_svd.shape[1] and z.shape[1] == x_svd.shape[0]:
        z = z.T
    if z.shape != x_svd.shape:
        raise ValueError(f"Unexpected Harmony output shape {z.shape}; expected {x_svd.shape}")
    return normalize(z)


def make_adata(ref: pd.DataFrame, query: pd.DataFrame, x_harmony: np.ndarray) -> AnnData:
    obs = pd.concat(
        [
            ref[["seurat_cell_id", "source_label", "adult_subtype", "adult_side"]].assign(div90_class="adult_reference"),
            query[["seurat_cell_id", "source_label", "div90_class"]].assign(
                adult_subtype="DIV90_query",
                adult_side="DIV90_query",
            ),
        ],
        ignore_index=True,
    )
    obs.index = obs["seurat_cell_id"].astype(str).to_numpy()
    adata = AnnData(X=sparse.csr_matrix((obs.shape[0], 1), dtype=np.float32), obs=obs)
    adata.obsm["X_pca"] = x_harmony.astype(np.float32)
    return adata


def run_forced_graph(adata: AnnData, args: argparse.Namespace) -> str:
    if importlib.util.find_spec("bbknn") is not None:
        import bbknn

        bbknn.bbknn(
            adata,
            batch_key="source_label",
            use_rep="X_pca",
            neighbors_within_batch=args.bbknn_neighbors_within_batch,
            n_pcs=args.bbknn_n_pcs,
            metric=args.bbknn_metric,
            trim=args.bbknn_trim,
            computation=args.bbknn_computation,
            annoy_n_trees=args.bbknn_annoy_n_trees,
            pynndescent_random_state=args.seed,
        )
        return "BBKNN"

    if importlib.util.find_spec("scanorama") is not None:
        raise RuntimeError("Scanorama fallback is available but not implemented for graph-forced UMAP in this script.")

    raise RuntimeError("Neither BBKNN nor Scanorama is available in the active Python environment.")


def add_umap_to_metadata(adata: AnnData, ref: pd.DataFrame, query: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    coords = adata.obsm["X_umap"]
    ref = ref.copy()
    query = query.copy()
    ref["umap_1"] = coords[: ref.shape[0], 0]
    ref["umap_2"] = coords[: ref.shape[0], 1]
    query["umap_1"] = coords[ref.shape[0] :, 0]
    query["umap_2"] = coords[ref.shape[0] :, 1]
    return ref, query


def style_axes(ax: plt.Axes) -> None:
    ax.set_xlabel("Forced integrated UMAP 1")
    ax.set_ylabel("Forced integrated UMAP 2")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def save_plot(fig: plt.Figure, plots_dir: Path, name: str) -> None:
    fig.savefig(plots_dir / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(plots_dir / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_source(ref: pd.DataFrame, query: pd.DataFrame, plots_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.4, 7.1))
    ax.scatter(ref["umap_1"], ref["umap_2"], s=1.6, c="#c8c8c8", alpha=0.38, linewidths=0, rasterized=True, label="adult_reference")
    ax.scatter(query["umap_1"], query["umap_2"], s=5.8, c="#111111", alpha=0.58, linewidths=0, rasterized=True, label="DIV90_query")
    style_axes(ax)
    ax.set_title("Forced integrated adult + DIV90 source check")
    ax.legend(frameon=False, markerscale=4, loc="best")
    save_plot(fig, plots_dir, "FINAL_forced_integrated_source_check")


def plot_adult_gray_div90_classes(ref: pd.DataFrame, query: pd.DataFrame, plots_dir: Path, palette: dict[str, str]) -> None:
    fig, ax = plt.subplots(figsize=(8.8, 7.2))
    ax.scatter(ref["umap_1"], ref["umap_2"], s=1.6, c="#d0d0d0", alpha=0.36, linewidths=0, rasterized=True, label="Adult reference")
    for label in sorted(palette):
        sub = query.loc[query["div90_class"] == label]
        ax.scatter(sub["umap_1"], sub["umap_2"], s=8, c=palette[label], alpha=0.86, linewidths=0, rasterized=True, label=label)
    style_axes(ax)
    ax.set_title("Forced Jia-style integration: adult gray, DIV90 classes")
    ax.legend(markerscale=3, fontsize=6.5, frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    save_plot(fig, plots_dir, "FINAL_forced_integrated_adult_gray_div90_classes")


def plot_adult_broad_div90_classes(ref: pd.DataFrame, query: pd.DataFrame, plots_dir: Path, palette: dict[str, str]) -> None:
    fig, ax = plt.subplots(figsize=(8.8, 7.2))
    for side in ["Cortical", "Subpallial", "Unassigned"]:
        sub = ref.loc[ref["adult_side"] == side]
        if sub.empty:
            continue
        ax.scatter(
            sub["umap_1"],
            sub["umap_2"],
            s=1.8,
            c=SIDE_PALETTE.get(side, "#b8b8b8"),
            alpha=0.24,
            linewidths=0,
            rasterized=True,
            label=f"Adult {side}",
        )
    for label in sorted(palette):
        sub = query.loc[query["div90_class"] == label]
        ax.scatter(sub["umap_1"], sub["umap_2"], s=8, c=palette[label], alpha=0.86, linewidths=0, rasterized=True, label=f"DIV90 {label}")
    style_axes(ax)
    ax.set_title("Forced Jia-style integration: adult broad labels, DIV90 classes")
    ax.legend(markerscale=3, fontsize=6.3, frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    save_plot(fig, plots_dir, "FINAL_forced_integrated_adult_broad_div90_classes")


def plot_adult_subtypes_div90_black(ref: pd.DataFrame, query: pd.DataFrame, plots_dir: Path) -> None:
    labels = [x for x in ADULT_ORDER if x in set(ref["adult_subtype"])]
    labels += sorted(set(ref["adult_subtype"]) - set(labels))
    fig, ax = plt.subplots(figsize=(8.8, 7.2))
    for label in labels:
        sub = ref.loc[ref["adult_subtype"] == label]
        ax.scatter(
            sub["umap_1"],
            sub["umap_2"],
            s=1.9,
            c=ADULT_PALETTE.get(label, "#777777"),
            alpha=0.45,
            linewidths=0,
            rasterized=True,
            label=label,
        )
    ax.scatter(query["umap_1"], query["umap_2"], s=5.8, c="#111111", alpha=0.52, linewidths=0, rasterized=True, label="DIV90 query")
    style_axes(ax)
    ax.set_title("Forced Jia-style integration: adult subtypes, DIV90 black")
    ax.legend(markerscale=4, fontsize=6.8, frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    save_plot(fig, plots_dir, "FINAL_forced_integrated_adult_subtypes_div90_black")


def plot_div90_only(query: pd.DataFrame, plots_dir: Path, palette: dict[str, str]) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 6.8))
    for label in sorted(palette):
        sub = query.loc[query["div90_class"] == label]
        ax.scatter(sub["umap_1"], sub["umap_2"], s=8.5, c=palette[label], alpha=0.88, linewidths=0, rasterized=True, label=label)
    style_axes(ax)
    ax.set_title("DIV90 classes on forced integrated coordinates")
    ax.legend(markerscale=3, fontsize=6.5, frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    save_plot(fig, plots_dir, "FINAL_forced_integrated_div90_classes_only")


def quantify_forced_graph(adata: AnnData, ref: pd.DataFrame, query: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    connectivities = adata.obsp["connectivities"].tocsr()
    obs = adata.obs.reset_index(drop=True)
    adult_idx = np.flatnonzero(obs["source_label"].to_numpy() == "adult_reference")
    adult_set = set(adult_idx.tolist())
    n_ref = ref.shape[0]

    rows = []
    for local_i, qrow in query.reset_index(drop=True).iterrows():
        graph_i = n_ref + local_i
        start, end = connectivities.indptr[graph_i], connectivities.indptr[graph_i + 1]
        neigh = connectivities.indices[start:end]
        weights = connectivities.data[start:end]
        adult_mask = np.array([j in adult_set for j in neigh], dtype=bool)
        adult_neigh = neigh[adult_mask]
        adult_weights = weights[adult_mask]
        adult_obs = obs.iloc[adult_neigh] if adult_neigh.size else obs.iloc[[]]

        side_counts = adult_obs["adult_side"].value_counts()
        subtype_counts = adult_obs["adult_subtype"].value_counts()
        n_adult = int(adult_neigh.size)
        n_total = int(neigh.size)
        n_cortical = int(side_counts.get("Cortical", 0))
        n_subpallial = int(side_counts.get("Subpallial", 0))
        rows.append(
            {
                "seurat_cell_id": qrow["seurat_cell_id"],
                "div90_class": qrow["div90_class"],
                "n_graph_neighbors_total": n_total,
                "n_adult_reference_graph_neighbors": n_adult,
                "fraction_graph_neighbors_adult_reference": float(n_adult / n_total) if n_total else np.nan,
                "n_adult_cortical_graph_neighbors": n_cortical,
                "n_adult_subpallial_graph_neighbors": n_subpallial,
                "fraction_adult_graph_neighbors_cortical": float(n_cortical / n_adult) if n_adult else np.nan,
                "fraction_adult_graph_neighbors_subpallial": float(n_subpallial / n_adult) if n_adult else np.nan,
                "top_adult_side_by_graph_neighbors": str(side_counts.idxmax()) if n_adult else "NA",
                "top_adult_subtype_by_graph_neighbors": str(subtype_counts.idxmax()) if n_adult else "NA",
                "mean_adult_graph_connectivity": float(np.mean(adult_weights)) if n_adult else np.nan,
                "max_adult_graph_connectivity": float(np.max(adult_weights)) if n_adult else np.nan,
            }
        )
    per_cell = pd.DataFrame(rows)

    summary_rows = []
    for div90_class, sub in per_cell.groupby("div90_class", sort=True):
        top_counts = sub["top_adult_subtype_by_graph_neighbors"].value_counts()
        summary_rows.append(
            {
                "div90_class": div90_class,
                "n_cells": int(sub.shape[0]),
                "mean_n_adult_reference_graph_neighbors": float(sub["n_adult_reference_graph_neighbors"].mean()),
                "median_n_adult_reference_graph_neighbors": float(sub["n_adult_reference_graph_neighbors"].median()),
                "mean_fraction_graph_neighbors_adult_reference": float(sub["fraction_graph_neighbors_adult_reference"].mean()),
                "mean_fraction_adult_graph_neighbors_cortical": float(sub["fraction_adult_graph_neighbors_cortical"].mean()),
                "mean_fraction_adult_graph_neighbors_subpallial": float(sub["fraction_adult_graph_neighbors_subpallial"].mean()),
                "top_adult_subtype_by_cellwise_graph_votes": str(top_counts.index[0]) if not top_counts.empty else "NA",
            }
        )
    return per_cell, pd.DataFrame(summary_rows)


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
    parser.add_argument("--harmony-theta", type=float, default=2.0)
    parser.add_argument("--harmony-lambda", type=float, default=1.0)
    parser.add_argument("--harmony-max-iter", type=int, default=20)
    parser.add_argument("--harmony-max-iter-kmeans", type=int, default=50)
    parser.add_argument("--harmony-epsilon-cluster", type=float, default=1e-5)
    parser.add_argument("--harmony-epsilon-harmony", type=float, default=1e-5)
    parser.add_argument("--bbknn-neighbors-within-batch", type=int, default=25)
    parser.add_argument("--bbknn-n-pcs", type=int, default=50)
    parser.add_argument("--bbknn-metric", default="euclidean")
    parser.add_argument("--bbknn-trim", type=parse_trim, default=None)
    parser.add_argument("--bbknn-computation", default="annoy")
    parser.add_argument("--bbknn-annoy-n-trees", type=int, default=50)
    parser.add_argument("--umap-min-dist", type=float, default=0.3)
    parser.add_argument("--umap-spread", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    plots_dir = args.outdir / "plots"
    tables_dir = args.outdir / "tables"
    reports_dir = args.outdir / "reports"
    for path in [plots_dir, tables_dir, reports_dir]:
        path.mkdir(parents=True, exist_ok=True)

    sc.settings.verbosity = 2
    ref, query, x_svd, metadata, diag = load_svd_inputs(args)
    x_harmony = run_harmony(x_svd, metadata, args)
    adata = make_adata(ref, query, x_harmony)
    method = run_forced_graph(adata, args)
    sc.tl.umap(adata, min_dist=args.umap_min_dist, spread=args.umap_spread, random_state=args.seed)
    ref, query = add_umap_to_metadata(adata, ref, query)

    palette = div90_palette(query)
    plot_source(ref, query, plots_dir)
    plot_adult_gray_div90_classes(ref, query, plots_dir, palette)
    plot_adult_broad_div90_classes(ref, query, plots_dir, palette)
    plot_adult_subtypes_div90_black(ref, query, plots_dir)
    plot_div90_only(query, plots_dir, palette)

    per_cell, summary = quantify_forced_graph(adata, ref, query)
    per_cell.to_csv(tables_dir / "forced_integration_neighbor_composition.tsv", sep="\t", index=False)
    summary.to_csv(tables_dir / "forced_integration_div90_to_adult_neighbor_summary.tsv", sep="\t", index=False)
    write_tsv_gz(ref, tables_dir / "forced_integrated_adult_reference_umap.tsv.gz")
    write_tsv_gz(query, tables_dir / "forced_integrated_div90_query_umap.tsv.gz")

    method_info = {
        "integration_method": method,
        "graph_method": "BBKNN batch-balanced neighbor graph" if method == "BBKNN" else method,
        "umap_source": "Scanpy UMAP fit from BBKNN connectivities in adata.obsp['connectivities']",
        "batch_key": "source_label",
        "source_labels": ["adult_reference", "DIV90_query"],
        "intent": "Forced Jia-style integrated visualization; not an adult-reference projection.",
    }
    parameters = {
        "bridge_dir": str(args.bridge_dir),
        "transfer_dir": str(args.transfer_dir),
        "outdir": str(args.outdir),
        "label_column": args.label_column,
        "exclude_label": args.exclude_label,
        "query_class_col": args.query_class_col,
        "nfeatures_requested": args.nfeatures,
        "n_components_requested": args.n_components,
        "harmony_theta": args.harmony_theta,
        "harmony_lambda": args.harmony_lambda,
        "harmony_max_iter": args.harmony_max_iter,
        "harmony_max_iter_kmeans": args.harmony_max_iter_kmeans,
        "harmony_epsilon_cluster": args.harmony_epsilon_cluster,
        "harmony_epsilon_harmony": args.harmony_epsilon_harmony,
        "bbknn_neighbors_within_batch": args.bbknn_neighbors_within_batch,
        "bbknn_n_pcs": args.bbknn_n_pcs,
        "bbknn_metric": args.bbknn_metric,
        "bbknn_trim": args.bbknn_trim,
        "bbknn_computation": args.bbknn_computation,
        "bbknn_annoy_n_trees": args.bbknn_annoy_n_trees,
        "umap_min_dist": args.umap_min_dist,
        "umap_spread": args.umap_spread,
        "seed": args.seed,
        **diag,
    }
    for target in [tables_dir / "forced_integration_method.json", args.outdir / "forced_integration_method.json"]:
        target.write_text(json.dumps(method_info, indent=2, sort_keys=True) + "\n")
    for target in [tables_dir / "forced_integration_parameters.json", args.outdir / "forced_integration_parameters.json"]:
        target.write_text(json.dumps(parameters, indent=2, sort_keys=True) + "\n")
    (tables_dir / "forced_integration_div90_palette.json").write_text(json.dumps(palette, indent=2, sort_keys=True) + "\n")

    readme = [
        "# Forced Jia-Style Integration: Adult Siletti + DIV90",
        "",
        "This output intentionally uses a batch-balanced graph for visualization:",
        "shared genes -> 3,000 selected transfer genes -> joint SVD -> Harmony(source_label) -> BBKNN(source_label) -> UMAP.",
        "",
        "The adult-reference projection, unintegrated joint embedding, and Harmony diagnostic outputs are not modified by this workflow.",
    ]
    (reports_dir / "README_forced_jia_style_integration.md").write_text("\n".join(readme) + "\n")

    print(json.dumps({"method": method_info, "parameters": parameters}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
