#!/usr/bin/env python
"""Audit whether Harmony is actually changing the Siletti/DIV90 latent space."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import harmonypy as hm
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import umap
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, silhouette_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import normalize

from plot_siletti_div90_jia_integration_sensitivity import (
    DEFAULT_BRIDGE_DIR,
    DEFAULT_TRANSFER_DIR,
    load_svd_inputs,
)


PROJECT_ROOT_DEFAULT = Path("/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder")
DEFAULT_OUTDIR = (
    PROJECT_ROOT_DEFAULT
    / "results/siletti_2023_whb_reference_label_transfer"
    / "harmony_integrated_adult_div90"
    / "audit"
)


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
    parser.add_argument("--theta", type=float, default=2.0)
    parser.add_argument("--lamb", type=float, default=1.0)
    parser.add_argument("--max-iter-harmony", type=int, default=20)
    parser.add_argument("--max-iter-kmeans", type=int, default=4)
    parser.add_argument("--epsilon-cluster", type=float, default=1e-3)
    parser.add_argument("--epsilon-harmony", type=float, default=1e-2)
    parser.add_argument("--aggressive-theta", type=float, default=4.0)
    parser.add_argument("--aggressive-max-iter-harmony", type=int, default=50)
    parser.add_argument("--aggressive-max-iter-kmeans", type=int, default=50)
    parser.add_argument("--aggressive-epsilon-cluster", type=float, default=1e-5)
    parser.add_argument("--aggressive-epsilon-harmony", type=float, default=1e-5)
    parser.add_argument("--near-identical-mean-abs-threshold", type=float, default=1e-6)
    parser.add_argument("--near-identical-median-corr-threshold", type=float, default=0.999999)
    parser.add_argument("--source-separation-sample-size", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=0)
    return parser.parse_args()


def orient_harmony_matrix(z_corr: np.ndarray, x_svd: np.ndarray) -> tuple[np.ndarray, str]:
    if z_corr.shape[0] == x_svd.shape[1] and z_corr.shape[1] == x_svd.shape[0]:
        return z_corr.T, "transposed_from_components_x_cells"
    if z_corr.shape == x_svd.shape:
        return z_corr, "already_cells_x_components"
    raise ValueError(f"Unexpected Harmony shape: {z_corr.shape}; X_svd shape: {x_svd.shape}")


def component_correlations(x_svd: np.ndarray, x_harmony: np.ndarray) -> list[float]:
    n = min(x_svd.shape[1], x_harmony.shape[1])
    out = []
    for i in range(n):
        corr = np.corrcoef(x_svd[:, i], x_harmony[:, i])[0, 1]
        out.append(float(corr))
    return out


def run_harmony_audit(
    x_svd: np.ndarray,
    meta: pd.DataFrame,
    *,
    theta: float,
    lamb: float,
    max_iter_harmony: int,
    max_iter_kmeans: int,
    epsilon_cluster: float,
    epsilon_harmony: float,
    seed: int,
    label: str,
) -> tuple[np.ndarray, dict]:
    print(f"[audit] Running Harmony label={label}", flush=True)
    print(f"[audit] X_svd shape: {x_svd.shape}", flush=True)
    print(f"[audit] metadata/source vector length: {meta.shape[0]}", flush=True)
    assert x_svd.shape[0] == meta.shape[0]
    ho = hm.run_harmony(
        data_mat=x_svd,
        meta_data=meta,
        vars_use=["source_label"],
        theta=theta,
        lamb=lamb,
        max_iter_harmony=max_iter_harmony,
        max_iter_kmeans=max_iter_kmeans,
        epsilon_cluster=epsilon_cluster,
        epsilon_harmony=epsilon_harmony,
        random_state=seed,
        verbose=True,
    )
    z = np.asarray(ho.Z_corr)
    print(f"[audit] ho.Z_corr shape: {z.shape}", flush=True)
    x_harmony_raw, orientation = orient_harmony_matrix(z, x_svd)
    x_harmony = normalize(x_harmony_raw)
    print(f"[audit] X_harmony_raw_oriented shape after orientation: {x_harmony_raw.shape}", flush=True)
    print(f"[audit] X_harmony shape after workflow normalization: {x_harmony.shape}", flush=True)
    print("[audit] UMAP input variable for Harmony sanity plot: X_harmony", flush=True)

    raw_diff = x_svd - x_harmony_raw
    raw_corrs = component_correlations(x_svd, x_harmony_raw)
    diff = x_svd - x_harmony
    corrs = component_correlations(x_svd, x_harmony)
    info = {
        "harmony_label": label,
        "X_svd_shape": "x".join(map(str, x_svd.shape)),
        "ho_Z_corr_shape": "x".join(map(str, z.shape)),
        "X_harmony_raw_oriented_shape_after_orientation": "x".join(map(str, x_harmony_raw.shape)),
        "X_harmony_shape_after_workflow_normalization": "x".join(map(str, x_harmony.shape)),
        "metadata_rows": int(meta.shape[0]),
        "orientation_action": orientation,
        "umap_input_variable": "X_harmony",
        "theta": theta,
        "lamb": lamb,
        "max_iter_harmony": max_iter_harmony,
        "max_iter_kmeans": max_iter_kmeans,
        "epsilon_cluster": epsilon_cluster,
        "epsilon_harmony": epsilon_harmony,
        "raw_mean_abs_diff": float(np.mean(np.abs(raw_diff))),
        "raw_fro_norm_diff": float(np.linalg.norm(raw_diff)),
        "raw_median_component_corr": float(np.nanmedian(raw_corrs)),
        "raw_min_component_corr": float(np.nanmin(raw_corrs)),
        "raw_max_component_corr": float(np.nanmax(raw_corrs)),
        "mean_abs_diff": float(np.mean(np.abs(diff))),
        "fro_norm_diff": float(np.linalg.norm(diff)),
        "median_component_corr": float(np.nanmedian(corrs)),
        "min_component_corr": float(np.nanmin(corrs)),
        "max_component_corr": float(np.nanmax(corrs)),
    }
    print(f"[audit] raw mean abs diff before workflow normalization: {info['raw_mean_abs_diff']}", flush=True)
    print(f"[audit] raw fro norm diff before workflow normalization: {info['raw_fro_norm_diff']}", flush=True)
    print(f"[audit] raw median component corr before workflow normalization: {info['raw_median_component_corr']}", flush=True)
    print(f"[audit] mean abs diff: {info['mean_abs_diff']}", flush=True)
    print(f"[audit] fro norm diff: {info['fro_norm_diff']}", flush=True)
    print(f"[audit] median component corr: {info['median_component_corr']}", flush=True)
    return x_harmony, info


def source_metrics(x: np.ndarray, source: np.ndarray, seed: int, sample_size: int) -> dict:
    y = (source == "DIV90_query").astype(int)
    idx = np.arange(x.shape[0])
    train_idx, test_idx = train_test_split(idx, test_size=0.3, random_state=seed, stratify=y)
    clf = LogisticRegression(max_iter=1000, solver="lbfgs")
    clf.fit(x[train_idx], y[train_idx])
    pred = clf.predict(x[test_idx])
    acc = accuracy_score(y[test_idx], pred)
    if sample_size and x.shape[0] > sample_size:
        rng = np.random.default_rng(seed)
        sample_idx = rng.choice(x.shape[0], size=sample_size, replace=False)
    else:
        sample_idx = idx
    sil = silhouette_score(x[sample_idx], y[sample_idx], metric="cosine")
    return {
        "source_classifier_accuracy": float(acc),
        "source_silhouette_cosine": float(sil),
        "silhouette_n_cells": int(sample_idx.size),
    }


def scatter_latent(x: np.ndarray, source: np.ndarray, title: str, outbase: Path, xlab: str, ylab: str) -> None:
    colors = np.where(source == "adult_reference", "#8f8f8f", "#111111")
    alpha = np.where(source == "adult_reference", 0.35, 0.45)
    fig, ax = plt.subplots(figsize=(7, 6))
    for src, color, a in [("adult_reference", "#8f8f8f", 0.35), ("DIV90_query", "#111111", 0.45)]:
        mask = source == src
        ax.scatter(x[mask, 0], x[mask, 1], s=3, c=color, alpha=a, linewidths=0, rasterized=True, label=src)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(title)
    ax.legend(frameon=False, markerscale=4)
    fig.savefig(outbase.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(outbase.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def umap_plot(x: np.ndarray, source: np.ndarray, title: str, outbase: Path, seed: int) -> None:
    print(f"[audit] Running UMAP for {title}", flush=True)
    coords = umap.UMAP(n_neighbors=30, min_dist=0.3, metric="cosine", random_state=seed).fit_transform(x)
    scatter_latent(coords, source, title, outbase, "UMAP 1", "UMAP 2")


def main() -> None:
    args = parse_args()
    plots_dir = args.outdir / "plots"
    tables_dir = args.outdir / "tables"
    for path in [plots_dir, tables_dir]:
        path.mkdir(parents=True, exist_ok=True)

    latent_args = argparse.Namespace(
        bridge_dir=args.bridge_dir,
        transfer_dir=args.transfer_dir,
        label_column=args.label_column,
        exclude_label=args.exclude_label,
        query_class_col=args.query_class_col,
        nfeatures=args.nfeatures,
        n_components=args.n_components,
        seed=args.seed,
    )
    ref, query, x_svd, meta, diag = load_svd_inputs(latent_args)
    source = meta["source_label"].astype(str).to_numpy()

    print(f"[audit] X_svd.shape = {x_svd.shape}", flush=True)
    print(f"[audit] metadata/source vector length = {meta.shape[0]}", flush=True)
    assert x_svd.shape[0] == meta.shape[0]

    x_harmony, current_info = run_harmony_audit(
        x_svd,
        meta,
        theta=args.theta,
        lamb=args.lamb,
        max_iter_harmony=args.max_iter_harmony,
        max_iter_kmeans=args.max_iter_kmeans,
        epsilon_cluster=args.epsilon_cluster,
        epsilon_harmony=args.epsilon_harmony,
        seed=args.seed,
        label="current",
    )
    matrix_rows = [current_info]

    near_identical = (
        current_info["mean_abs_diff"] <= args.near_identical_mean_abs_threshold
        and current_info["median_component_corr"] >= args.near_identical_median_corr_threshold
    )
    x_harmony_for_plots = x_harmony
    if near_identical:
        x_harmony_aggressive, aggressive_info = run_harmony_audit(
            x_svd,
            meta,
            theta=args.aggressive_theta,
            lamb=args.lamb,
            max_iter_harmony=args.aggressive_max_iter_harmony,
            max_iter_kmeans=args.aggressive_max_iter_kmeans,
            epsilon_cluster=args.aggressive_epsilon_cluster,
            epsilon_harmony=args.aggressive_epsilon_harmony,
            seed=args.seed,
            label="aggressive",
        )
        matrix_rows.append(aggressive_info)
        x_harmony_for_plots = x_harmony_aggressive

    pd.DataFrame(matrix_rows).to_csv(tables_dir / "harmony_matrix_audit.tsv", sep="\t", index=False)

    sep_rows = []
    for label, mat in [("X_svd_uncorrected", x_svd), ("X_harmony_current", x_harmony)]:
        row = {"matrix": label}
        row.update(source_metrics(mat, source, args.seed, args.source_separation_sample_size))
        sep_rows.append(row)
    if near_identical:
        row = {"matrix": "X_harmony_aggressive"}
        row.update(source_metrics(x_harmony_for_plots, source, args.seed, args.source_separation_sample_size))
        sep_rows.append(row)
    pd.DataFrame(sep_rows).to_csv(tables_dir / "harmony_source_separation_before_after.tsv", sep="\t", index=False)

    params = {
        "bridge_dir": str(args.bridge_dir),
        "transfer_dir": str(args.transfer_dir),
        "outdir": str(args.outdir),
        "X_svd_shape": list(x_svd.shape),
        "metadata_rows": int(meta.shape[0]),
        "current_harmony_parameters": {
            "theta": args.theta,
            "lamb": args.lamb,
            "max_iter_harmony": args.max_iter_harmony,
            "max_iter_kmeans": args.max_iter_kmeans,
            "epsilon_cluster": args.epsilon_cluster,
            "epsilon_harmony": args.epsilon_harmony,
        },
        "aggressive_rerun_triggered": bool(near_identical),
        "aggressive_harmony_parameters": {
            "theta": args.aggressive_theta,
            "lamb": args.lamb,
            "max_iter_harmony": args.aggressive_max_iter_harmony,
            "max_iter_kmeans": args.aggressive_max_iter_kmeans,
            "epsilon_cluster": args.aggressive_epsilon_cluster,
            "epsilon_harmony": args.aggressive_epsilon_harmony,
        },
        "latent_input_diagnostics": diag,
    }
    (tables_dir / "harmony_run_parameters.json").write_text(json.dumps(params, indent=2, sort_keys=True) + "\n")

    scatter_latent(x_svd, source, "Uncorrected X_svd: component 1 vs 2 by source", plots_dir / "X_svd_pc1_pc2_by_source", "SVD 1", "SVD 2")
    scatter_latent(x_harmony_for_plots, source, "Harmony latent: component 1 vs 2 by source", plots_dir / "X_harmony_h1_h2_by_source", "Harmony 1", "Harmony 2")
    print("[audit] UMAP input variable for uncorrected sanity plot: X_svd", flush=True)
    umap_plot(x_svd, source, "UMAP from X_svd by source", plots_dir / "umap_from_X_svd_by_source", args.seed)
    print("[audit] UMAP input variable for Harmony sanity plot: X_harmony", flush=True)
    umap_plot(x_harmony_for_plots, source, "UMAP from X_harmony by source", plots_dir / "umap_from_X_harmony_by_source", args.seed)

    print(json.dumps(matrix_rows, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
