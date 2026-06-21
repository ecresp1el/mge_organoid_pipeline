#!/usr/bin/env python
"""Verify Siletti/DIV90 unintegrated versus Harmony integration outputs.

This script does not create another biological mapping. It exports diagnostics
that make the current integration behavior explicit:

* UMAP input provenance for Harmony.
* Source-colored unintegrated and Harmony UMAPs.
* Harmony overlay plots used for inspection.
* Nearest-adult distance percentiles.
* Adult-neighbor subtype and cortical/subpallial composition.
* Source-mixing scores per DIV90 class in latent space.
"""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

from plot_siletti_div90_jia_integration_sensitivity import (
    DEFAULT_BRIDGE_DIR,
    DEFAULT_TRANSFER_DIR,
    DEFAULT_UNINTEGRATED_DIR,
    load_svd_inputs,
    run_harmony,
)


PROJECT_ROOT_DEFAULT = Path("/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder")
DEFAULT_HARMONY_DIR = (
    PROJECT_ROOT_DEFAULT
    / "results/siletti_2023_whb_reference_label_transfer/harmony_integrated_adult_div90"
)
DEFAULT_OUTDIR = (
    PROJECT_ROOT_DEFAULT
    / "results/siletti_2023_whb_reference_label_transfer/integration_sensitivity_diagnostics_adult_div90"
)


def read_tsv_gz(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt") as handle:
        return pd.read_csv(handle, sep="\t")


def copy_if_exists(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.exists():
        shutil.copy2(src, dst)


def source_plot(ref: pd.DataFrame, query: pd.DataFrame, title: str, outbase: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 7.0))
    ax.scatter(
        ref["umap_1"],
        ref["umap_2"],
        s=2.0,
        c="#9b9b9b",
        alpha=0.42,
        linewidths=0,
        rasterized=True,
        label="adult_reference",
    )
    ax.scatter(
        query["umap_1"],
        query["umap_2"],
        s=5.0,
        c="#111111",
        alpha=0.44,
        linewidths=0,
        rasterized=True,
        label="DIV90_query",
    )
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title(title)
    ax.legend(frameon=False, markerscale=4, loc="best")
    fig.savefig(outbase.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(outbase.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def side_by_side_source_plot(
    un_ref: pd.DataFrame,
    un_query: pd.DataFrame,
    ha_ref: pd.DataFrame,
    ha_query: pd.DataFrame,
    outbase: Path,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.8))
    panels = [
        (axes[0], un_ref, un_query, "Unintegrated joint SVD/UMAP"),
        (axes[1], ha_ref, ha_query, "Harmony latent UMAP"),
    ]
    for ax, ref, query, title in panels:
        ax.scatter(ref["umap_1"], ref["umap_2"], s=1.8, c="#9b9b9b", alpha=0.42, linewidths=0, rasterized=True, label="adult_reference")
        ax.scatter(query["umap_1"], query["umap_2"], s=4.2, c="#111111", alpha=0.44, linewidths=0, rasterized=True, label="DIV90_query")
        ax.set_title(title)
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
    axes[1].legend(frameon=False, markerscale=4, loc="best")
    fig.savefig(outbase.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(outbase.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def percentiles(df: pd.DataFrame, method: str, distance_col: str) -> pd.DataFrame:
    probs = [0, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1.0]
    rows = []
    groups = [("ALL_DIV90", df)]
    groups.extend(list(df.groupby("div90_class", sort=True)))
    for label, sub in groups:
        values = pd.to_numeric(sub[distance_col], errors="coerce").dropna().to_numpy()
        row = {"method": method, "div90_class": label, "n_cells": int(values.size), "distance_metric": distance_col}
        if values.size:
            qs = np.quantile(values, probs)
            for prob, value in zip(probs, qs):
                row[f"p{int(prob * 100):02d}"] = float(value)
        rows.append(row)
    return pd.DataFrame(rows)


def subtype_composition(edge_df: pd.DataFrame, method: str) -> pd.DataFrame:
    comp = (
        edge_df.groupby(["div90_class", "adult_subtype"], dropna=False)["n_neighbors"]
        .sum()
        .reset_index()
    )
    comp["total_neighbors"] = comp.groupby("div90_class")["n_neighbors"].transform("sum")
    comp["fraction_neighbors"] = comp["n_neighbors"] / comp["total_neighbors"]
    comp.insert(0, "method", method)
    return comp.sort_values(["div90_class", "fraction_neighbors"], ascending=[True, False])


def broad_side_composition(metrics: pd.DataFrame, method: str) -> pd.DataFrame:
    out = (
        metrics.groupby("div90_class", sort=True)
        .agg(
            n_cells=("seurat_cell_id", "size"),
            mean_fraction_adult_neighbors_cortical=("fraction_adult_neighbors_cortical", "mean"),
            mean_fraction_adult_neighbors_subpallial=("fraction_adult_neighbors_subpallial", "mean"),
            median_nearest_adult_cosine_distance=("nearest_adult_cosine_distance", "median"),
        )
        .reset_index()
    )
    out.insert(0, "method", method)
    return out


def source_mixing_scores(latent: np.ndarray, query_classes: pd.Series, n_ref: int, k: int, method: str) -> pd.DataFrame:
    k_use = min(k + 1, latent.shape[0])
    nn = NearestNeighbors(n_neighbors=k_use, metric="cosine", algorithm="brute", n_jobs=-1)
    nn.fit(latent)
    _, idx = nn.kneighbors(latent[n_ref:, :])
    idx = idx[:, 1:] if idx.shape[1] > k else idx[:, :k]
    is_adult = idx < n_ref
    per_cell = pd.DataFrame(
        {
            "div90_class": query_classes.to_numpy(),
            "fraction_all_neighbors_adult_reference": is_adult.mean(axis=1),
            "fraction_all_neighbors_div90_query": 1.0 - is_adult.mean(axis=1),
        }
    )
    out = (
        per_cell.groupby("div90_class", sort=True)
        .agg(
            n_cells=("fraction_all_neighbors_adult_reference", "size"),
            mean_fraction_all_neighbors_adult_reference=("fraction_all_neighbors_adult_reference", "mean"),
            median_fraction_all_neighbors_adult_reference=("fraction_all_neighbors_adult_reference", "median"),
            mean_fraction_all_neighbors_div90_query=("fraction_all_neighbors_div90_query", "mean"),
        )
        .reset_index()
    )
    out.insert(0, "method", method)
    out["neighbor_k"] = k
    return out


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bridge-dir", type=Path, default=DEFAULT_BRIDGE_DIR)
    parser.add_argument("--transfer-dir", type=Path, default=DEFAULT_TRANSFER_DIR)
    parser.add_argument("--unintegrated-dir", type=Path, default=DEFAULT_UNINTEGRATED_DIR)
    parser.add_argument("--harmony-dir", type=Path, default=DEFAULT_HARMONY_DIR)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--label-column", default="candidate_jia_group")
    parser.add_argument("--exclude-label", default="Excluded / not assigned to Jia-style 9 groups")
    parser.add_argument("--query-class-col", default="div90_broad_class")
    parser.add_argument("--nfeatures", type=int, default=3000)
    parser.add_argument("--n-components", type=int, default=50)
    parser.add_argument("--neighbor-k", type=int, default=50)
    parser.add_argument("--source-mixing-k", type=int, default=50)
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

    un_ref = read_tsv_gz(args.unintegrated_dir / "tables/siletti_reference_cells_joint_umap.tsv.gz")
    un_query = read_tsv_gz(args.unintegrated_dir / "tables/div90_query_cells_joint_umap.tsv.gz")
    ha_ref = read_tsv_gz(args.harmony_dir / "tables/siletti_reference_cells_harmony_umap.tsv.gz")
    ha_query = read_tsv_gz(args.harmony_dir / "tables/div90_query_cells_harmony_umap.tsv.gz")

    source_plot(un_ref, un_query, "Unintegrated joint SVD/UMAP by source", plots_dir / "unintegrated_joint_umap_by_source")
    source_plot(ha_ref, ha_query, "Harmony UMAP by source", plots_dir / "harmony_umap_by_source")
    side_by_side_source_plot(un_ref, un_query, ha_ref, ha_query, plots_dir / "source_umap_unintegrated_vs_harmony_side_by_side")

    copy_if_exists(
        args.harmony_dir / "plots/joint_umap_adult_subtypes_div90_black_overlay.png",
        plots_dir / "harmony_umap_adult_subtypes_div90_black_overlay.png",
    )
    copy_if_exists(
        args.harmony_dir / "plots/joint_umap_adult_subtypes_div90_black_overlay.pdf",
        plots_dir / "harmony_umap_adult_subtypes_div90_black_overlay.pdf",
    )
    copy_if_exists(
        args.harmony_dir / "plots/joint_umap_adult_side_div90_classes_overlay.png",
        plots_dir / "harmony_umap_adult_side_div90_classes_overlay.png",
    )
    copy_if_exists(
        args.harmony_dir / "plots/joint_umap_adult_side_div90_classes_overlay.pdf",
        plots_dir / "harmony_umap_adult_side_div90_classes_overlay.pdf",
    )

    un_metrics = read_tsv_gz(args.unintegrated_dir / "tables/div90_cell_adult_neighbor_coembedding_metrics.tsv.gz")
    ha_metrics = read_tsv_gz(args.harmony_dir / "tables/div90_cell_adult_neighbor_harmony_metrics.tsv.gz")
    pd.concat(
        [
            percentiles(un_metrics, "unintegrated_joint_svd", "nearest_adult_cosine_distance"),
            percentiles(ha_metrics, "harmony_latent", "nearest_adult_cosine_distance"),
        ],
        ignore_index=True,
    ).to_csv(tables_dir / "nearest_adult_distance_percentiles.tsv", sep="\t", index=False)

    un_edges = pd.read_csv(args.unintegrated_dir / "tables/div90_cell_adult_neighbor_subtype_composition.tsv", sep="\t")
    ha_edges = pd.read_csv(args.harmony_dir / "tables/div90_cell_adult_neighbor_subtype_composition_harmony.tsv", sep="\t")
    pd.concat(
        [
            subtype_composition(un_edges, "unintegrated_joint_svd"),
            subtype_composition(ha_edges, "harmony_latent"),
        ],
        ignore_index=True,
    ).to_csv(tables_dir / "nearest_adult_subtype_composition_per_div90_class.tsv", sep="\t", index=False)
    pd.concat(
        [
            broad_side_composition(un_metrics, "unintegrated_joint_svd"),
            broad_side_composition(ha_metrics, "harmony_latent"),
        ],
        ignore_index=True,
    ).to_csv(tables_dir / "nearest_adult_cortical_subpallial_composition_per_div90_class.tsv", sep="\t", index=False)

    # Recompute latent coordinates exactly enough for source-mixing summaries.
    latent_args = argparse.Namespace(
        bridge_dir=args.bridge_dir,
        transfer_dir=args.transfer_dir,
        label_column=args.label_column,
        exclude_label=args.exclude_label,
        query_class_col=args.query_class_col,
        nfeatures=args.nfeatures,
        n_components=args.n_components,
        seed=args.seed,
        harmony_theta=args.harmony_theta,
        harmony_lambda=args.harmony_lambda,
        harmony_max_iter=args.harmony_max_iter,
    )
    ref_meta, query_meta, joint_svd, harmony_meta, latent_diag = load_svd_inputs(latent_args)
    harmony_latent = run_harmony(joint_svd, harmony_meta, latent_args)
    n_ref = ref_meta.shape[0]
    pd.concat(
        [
            source_mixing_scores(joint_svd, query_meta["div90_class"], n_ref, args.source_mixing_k, "unintegrated_joint_svd"),
            source_mixing_scores(harmony_latent, query_meta["div90_class"], n_ref, args.source_mixing_k, "harmony_latent"),
        ],
        ignore_index=True,
    ).to_csv(tables_dir / "source_mixing_score_per_div90_class.tsv", sep="\t", index=False)

    harmony_config = json.loads((args.harmony_dir / "tables/harmony_integration_config.json").read_text())
    umap_provenance = {
        "harmony_plot_output_folder": str(args.harmony_dir),
        "umap_was_fit_on_harmony_corrected_latent_coordinates": True,
        "exact_umap_input_matrix_name_in_script": "corrected",
        "exact_umap_input_description": "normalize(harmony.Z_corr) with shape cells x Harmony components, produced from joint_svd corrected for source_label",
        "script": "python_notebooks/scripts/plot_siletti_div90_jia_integration_sensitivity.py",
        "function": "add_integrated_umap(ref, query, corrected, args)",
        "harmony_config": harmony_config,
    }
    (tables_dir / "harmony_umap_input_provenance.json").write_text(json.dumps(umap_provenance, indent=2, sort_keys=True) + "\n")
    (tables_dir / "verification_config.json").write_text(
        json.dumps(
            {
                "outdir": str(args.outdir),
                "unintegrated_dir": str(args.unintegrated_dir),
                "harmony_dir": str(args.harmony_dir),
                "bridge_dir": str(args.bridge_dir),
                "transfer_dir": str(args.transfer_dir),
                "source_mixing_k": args.source_mixing_k,
                **latent_diag,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )

    readme = [
        "# Siletti DIV90 Integration Verification",
        "",
        "This folder verifies the existing unintegrated and Harmony integration outputs.",
        "",
        "Harmony UMAP provenance:",
        "- UMAP input: Harmony-corrected latent matrix named `corrected` in the script.",
        "- `corrected` is `normalize(harmony.Z_corr)` after source-label correction.",
        "- UMAP was not fit on raw SVD for the Harmony plots.",
        "",
        "Interpretation guardrail:",
        "- Adult-nearest territory assignment means DIV90 cells are assigned to the nearest adult subtype/side.",
        "- True adult-territory embedding requires distances/mixing comparable to adult-adult neighborhoods.",
    ]
    (reports_dir / "README_integration_verification.md").write_text("\n".join(readme) + "\n")

    print(json.dumps(umap_provenance, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
