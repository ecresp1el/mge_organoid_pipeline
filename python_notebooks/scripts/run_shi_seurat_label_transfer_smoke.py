#!/usr/bin/env python3
"""Run Seurat-based Shi label transfer smoke plots for DIV30."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mge_organoid_python.shi_label_transfer import (
    attach_reference_labels,
    broad_region_from_label,
    developmental_class_from_label,
    get_umap_coordinates,
    load_shi_table_s2_labels,
    natural_sort_key,
    plot_overlaid_density_by_group,
    plot_umap_categorical,
    plot_umap_continuous,
    reference_week_metadata,
    safe_token,
    shi_week_numeric_from_label,
    summarize_predictions_by_cluster,
)
from mge_organoid_python.shi_prediction_schema import (
    PREDICTION_SCORE_COL,
    UNCERTAINTY_SCORE_COL,
    WEEK_PREDICTION_SCORE_COL,
    WEEK_UNCERTAINTY_SCORE_COL,
    sanitize_shi_label_token,
    validate_canonical_prediction_scores,
)
from scipy.stats import gaussian_kde


MARKER = "[ShiSeuratTransfer]"


def log(message: str) -> None:
    print(f"{MARKER} {message}", flush=True)


def env_path(name: str, default: str) -> Path:
    return Path(os.environ.get(name, default)).expanduser()


def env_int(name: str, default: int) -> int:
    return int(os.environ.get(name, str(default)))


def env_bool(name: str, default: bool) -> bool:
    raw = os.environ.get(name)
    if raw is None:
        return default
    return str(raw).strip().lower() in {"1", "true", "yes", "y"}


def sanitize_score_label(label: str) -> str:
    return sanitize_shi_label_token(label)


def prepare_reference_labels(
    reference_h5ad: Path,
    table_s2_xlsx: Path,
    labels_path: Path,
    table_dir: Path,
    min_match_fraction: float,
) -> list[str]:
    log(f"Reading Shi Table S2 labels: {table_s2_xlsx}")
    table_s2 = load_shi_table_s2_labels(table_s2_xlsx)

    log(f"Reading Shi reference AnnData metadata: {reference_h5ad}")
    reference = ad.read_h5ad(reference_h5ad, backed="r")
    try:
        label_obs, join_summary = attach_reference_labels(
            reference,
            table_s2,
            min_match_fraction=min_match_fraction,
        )
        week_obs = reference_week_metadata(reference.obs)
        week_numeric = pd.to_numeric(week_obs["shi_week_numeric"], errors="coerce")
        week_obs["shi_week_numeric"] = week_numeric
        week_obs["shi_week_label"] = week_numeric.map(
            lambda value: f"GW{int(value):02d}" if pd.notna(value) else pd.NA
        ).astype("string")
        label_obs = label_obs.join(week_obs[["shi_week_label", "shi_week_numeric"]])
        for source_col, output_col in [
            ("cell_id", "reference_cell_id"),
            ("raw_cell_id", "reference_raw_cell_id"),
            ("barcode", "reference_barcode"),
        ]:
            if source_col in reference.obs.columns:
                label_obs[output_col] = reference.obs[source_col].astype(str).to_numpy()
            else:
                label_obs[output_col] = ""
    finally:
        reference.file.close()

    labels_path.parent.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)
    label_obs.reset_index(drop=True).to_csv(labels_path, sep="\t", index=False)
    join_summary.to_csv(table_dir / "shi_reference_table_s2_label_join_summary.tsv", sep="\t", index=False)
    counts = (
        label_obs["shi_label"]
        .dropna()
        .astype(str)
        .value_counts()
        .rename_axis("shi_label")
        .reset_index(name="n_reference_cells")
    )
    counts.to_csv(table_dir / "shi_reference_table_s2_label_counts.tsv", sep="\t", index=False)
    log(f"Wrote Seurat label bridge: {labels_path}")
    return counts["shi_label"].astype(str).tolist()


def run_seurat_transfer(
    repo_root: Path,
    rscript_bin: Path,
    reference_seurat: Path,
    query_seurat: Path,
    labels_path: Path,
    seurat_dir: Path,
    output_prefix: str,
    dims: int,
    nfeatures: int,
    npcs: int,
    seed: int,
) -> None:
    script = repo_root / "python_notebooks" / "scripts" / "seurat_shi_label_transfer_export.R"
    cmd = [
        str(rscript_bin),
        str(script),
        "--reference",
        str(reference_seurat),
        "--query",
        str(query_seurat),
        "--labels",
        str(labels_path),
        "--outdir",
        str(seurat_dir),
        "--reference_assay",
        os.environ.get("SHI_SEURAT_REFERENCE_ASSAY", "RNA"),
        "--query_assay",
        os.environ.get("SHI_SEURAT_QUERY_ASSAY", "RNA"),
        "--transfer_name",
        "shi_label",
        "--output_prefix",
        output_prefix,
        "--normalization_method",
        os.environ.get("SHI_SEURAT_NORMALIZATION_METHOD", "LogNormalize"),
        "--dims",
        str(dims),
        "--nfeatures",
        str(nfeatures),
        "--npcs",
        str(npcs),
        "--seed",
        str(seed),
    ]
    log("Running Seurat transfer/export")
    log("Command: " + " ".join(cmd))
    subprocess.run(cmd, check=True)


def load_predictions(prediction_path: Path) -> pd.DataFrame:
    predictions = pd.read_csv(prediction_path, sep="\t")
    if "cell_id" not in predictions.columns:
        raise ValueError(f"Missing cell_id in Seurat predictions: {prediction_path}")
    if "predicted.id" not in predictions.columns:
        raise ValueError(f"Missing predicted.id in Seurat predictions: {prediction_path}")
    if "prediction.score.max" not in predictions.columns:
        raise ValueError(f"Missing prediction.score.max in Seurat predictions: {prediction_path}")
    predictions["cell_id"] = predictions["cell_id"].astype(str)
    return predictions


def load_week_predictions(prediction_path: Path) -> pd.DataFrame:
    predictions = load_predictions(prediction_path)
    predictions = predictions.rename(
        columns={
            "predicted.id": "predicted_week_label",
            "prediction.score.max": "week_prediction_score",
        }
    )
    return predictions


def align_predictions_to_adata(adata: ad.AnnData, predictions: pd.DataFrame) -> pd.DataFrame:
    pred_by_cell = predictions.set_index("cell_id", drop=False)
    obs_names = pd.Index(adata.obs_names.astype(str))
    direct_matches = obs_names.isin(pred_by_cell.index).sum()
    if direct_matches == adata.n_obs:
        log("Joining Seurat predictions by AnnData obs_names")
        return pred_by_cell.loc[obs_names].reset_index(drop=True)

    if "cell_id" in adata.obs.columns:
        cell_ids = pd.Index(adata.obs["cell_id"].astype(str))
        cell_matches = cell_ids.isin(pred_by_cell.index).sum()
        if cell_matches == adata.n_obs:
            log('Joining Seurat predictions by adata.obs["cell_id"]')
            return pred_by_cell.loc[cell_ids].reset_index(drop=True)

    diagnostics = pd.DataFrame(
        [
            {"join_key": "obs_names", "matched_cells": int(direct_matches), "query_cells": int(adata.n_obs)},
            {
                "join_key": "obs.cell_id",
                "matched_cells": int(adata.obs["cell_id"].astype(str).isin(pred_by_cell.index).sum())
                if "cell_id" in adata.obs.columns
                else 0,
                "query_cells": int(adata.n_obs),
            },
        ]
    )
    raise ValueError(
        "Could not fully join Seurat predictions to AnnData cells:\n"
        + diagnostics.to_string(index=False)
    )


def add_predictions_to_obs(adata: ad.AnnData, aligned: pd.DataFrame, label_order: list[str]) -> list[str]:
    prefix = "shi_seurat_full"
    adata.obs[f"{prefix}_predicted_shi_label"] = aligned["predicted.id"].astype(str).to_numpy()
    adata.obs[f"{prefix}_prediction_score"] = pd.to_numeric(
        aligned["prediction.score.max"],
        errors="coerce",
    ).to_numpy()
    adata.obs[f"{prefix}_uncertainty_score"] = 1.0 - adata.obs[f"{prefix}_prediction_score"].astype(float)
    adata.obs[f"{prefix}_broad_region_class"] = adata.obs[f"{prefix}_predicted_shi_label"].map(broad_region_from_label)
    adata.obs[f"{prefix}_developmental_class"] = adata.obs[f"{prefix}_predicted_shi_label"].map(
        developmental_class_from_label
    )

    score_cols = [col for col in aligned.columns if col.startswith("prediction.score.") and col != "prediction.score.max"]
    label_score_cols = []
    records = []
    for col in score_cols:
        label = col.removeprefix("prediction.score.")
        token = sanitize_score_label(label)
        out_col = f"{prefix}_prediction_score_{token}"
        adata.obs[out_col] = pd.to_numeric(aligned[col], errors="coerce").to_numpy()
        label_score_cols.append(out_col)
        records.append({"shi_label": label, "score_column": out_col})

    score_map = pd.DataFrame.from_records(records)
    adata.uns["shi_seurat_full_score_column_map"] = score_map.to_dict(orient="list")
    if label_order:
        adata.obs[f"{prefix}_predicted_shi_label"] = pd.Categorical(
            adata.obs[f"{prefix}_predicted_shi_label"].astype(str),
            categories=[str(label) for label in label_order],
        )
    return label_score_cols


def add_week_predictions_to_obs(adata: ad.AnnData, aligned: pd.DataFrame) -> list[str]:
    prefix = "shi_seurat_full"
    adata.obs[f"{prefix}_predicted_shi_week_label"] = aligned["predicted_week_label"].astype(str).to_numpy()
    adata.obs[f"{prefix}_week_prediction_score"] = pd.to_numeric(
        aligned["week_prediction_score"],
        errors="coerce",
    ).to_numpy()
    adata.obs[f"{prefix}_week_uncertainty_score"] = 1.0 - adata.obs[f"{prefix}_week_prediction_score"].astype(float)
    adata.obs[f"{prefix}_predicted_shi_week_numeric"] = adata.obs[f"{prefix}_predicted_shi_week_label"].map(
        shi_week_numeric_from_label
    )

    score_cols = [col for col in aligned.columns if col.startswith("prediction.score.") and col != "prediction.score.max"]
    week_score_cols = []
    weighted_sum = np.zeros(adata.n_obs, dtype=float)
    weight_total = np.zeros(adata.n_obs, dtype=float)
    records = []
    for col in score_cols:
        label = col.removeprefix("prediction.score.")
        token = sanitize_score_label(label)
        out_col = f"{prefix}_week_prediction_score_{token}"
        values = pd.to_numeric(aligned[col], errors="coerce").fillna(0.0).to_numpy(dtype=float)
        adata.obs[out_col] = values
        week_score_cols.append(out_col)
        week_numeric = shi_week_numeric_from_label(label)
        if np.isfinite(week_numeric):
            weighted_sum += values * week_numeric
            weight_total += values
        records.append({"shi_week_label": label, "score_column": out_col, "shi_week_numeric": week_numeric})

    expected = np.divide(
        weighted_sum,
        weight_total,
        out=np.full_like(weighted_sum, np.nan, dtype=float),
        where=weight_total > 0,
    )
    adata.obs[f"{prefix}_expected_shi_week_numeric"] = expected
    score_map = pd.DataFrame.from_records(records)
    adata.uns["shi_seurat_full_week_score_column_map"] = score_map.to_dict(orient="list")
    return week_score_cols


def plot_score_panel(
    coords: np.ndarray,
    obs: pd.DataFrame,
    score_cols: list[str],
    labels: list[str],
    path: Path,
) -> None:
    fig, axes = plt.subplots(1, len(score_cols), figsize=(5.2 * len(score_cols), 4.8), squeeze=False)
    for ax, col, label in zip(axes.ravel(), score_cols, labels, strict=True):
        values = pd.to_numeric(obs[col], errors="coerce")
        artist = ax.scatter(
            coords[:, 0],
            coords[:, 1],
            c=values,
            s=1.0,
            cmap="Blues",
            vmin=0.0,
            vmax=1.0,
            linewidths=0,
            rasterized=True,
        )
        ax.set_title(label)
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(artist, ax=ax, label="Seurat prediction score")
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


NON_GREY_COLORS = [
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#17becf",
    "#bcbd22",
    "#e41a1c",
    "#377eb8",
    "#4daf4a",
    "#984ea3",
    "#ff7f00",
    "#a65628",
]


def color_for_index(index: int) -> str:
    return NON_GREY_COLORS[index % len(NON_GREY_COLORS)]


def label_color_map(labels: list[str]) -> dict[str, str]:
    return {str(label): color_for_index(idx) for idx, label in enumerate(labels)}


def style_umap_axis(ax: plt.Axes) -> None:
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_xticks([])
    ax.set_yticks([])


def style_small_umap_axis(ax: plt.Axes, coords: np.ndarray) -> None:
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(float(np.nanmin(coords[:, 0])), float(np.nanmax(coords[:, 0])))
    ax.set_ylim(float(np.nanmin(coords[:, 1])), float(np.nanmax(coords[:, 1])))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")


def plot_predicted_label_with_overlay_panels(
    coords: np.ndarray,
    labels: pd.Series,
    label_order: list[str],
    path: Path,
) -> None:
    values = labels.astype(str)
    categories = [str(label) for label in label_order]
    if not categories:
        categories = sorted(values.dropna().unique())
    left_cols = min(4, max(1, len(categories)))
    rows = int(np.ceil(len(categories) / left_cols))
    fig = plt.figure(figsize=(3.8 * left_cols + 6.0, 3.3 * rows))
    gs = fig.add_gridspec(rows, left_cols + 2, width_ratios=[1] * left_cols + [1.35, 1.35])
    colors = label_color_map(categories)

    for idx, label in enumerate(categories):
        row = idx // left_cols
        col = idx % left_cols
        ax = fig.add_subplot(gs[row, col])
        mask = values.to_numpy() == str(label)
        ax.scatter(coords[:, 0], coords[:, 1], c="#d0d0d0", s=0.6, linewidths=0, rasterized=True)
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            c=colors[label],
            s=1.4,
            linewidths=0,
            rasterized=True,
        )
        ax.set_title(f"{label} (n={int(mask.sum())})", fontsize=9)
        style_small_umap_axis(ax, coords)

    for idx in range(len(categories), rows * left_cols):
        row = idx // left_cols
        col = idx % left_cols
        ax = fig.add_subplot(gs[row, col])
        ax.axis("off")

    ax_full = fig.add_subplot(gs[:, left_cols:])
    for label in categories:
        mask = values.to_numpy() == str(label)
        ax_full.scatter(
            coords[mask, 0],
            coords[mask, 1],
            c=colors[label],
            s=1.0,
            linewidths=0,
            rasterized=True,
            label=f"{label} (n={int(mask.sum())})",
        )
    ax_full.set_title("All predicted Shi labels")
    style_umap_axis(ax_full)
    ax_full.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_stacked_bar_with_cluster_umaps(
    coords: np.ndarray,
    obs: pd.DataFrame,
    label_counts: pd.DataFrame,
    cluster_col: str,
    label_col: str,
    label_order: list[str],
    path: Path,
    title: str,
) -> None:
    clusters = sorted(obs[cluster_col].dropna().astype(str).unique(), key=natural_sort_key)
    cluster_cols = 2 if len(clusters) > 1 else 1
    rows = int(np.ceil(len(clusters) / cluster_cols))
    fig = plt.figure(figsize=(17.0, max(7.2, 3.8 * rows)))
    outer = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.2], wspace=0.24)

    ax_bar = fig.add_subplot(outer[0, 0])
    wide = label_counts.pivot(index=cluster_col, columns=label_col, values="fraction").fillna(0.0)
    wide = wide.reindex(columns=[str(label) for label in label_order], fill_value=0.0)
    wide = wide.loc[sorted(wide.index.astype(str), key=natural_sort_key)]
    colors = [label_color_map([str(label) for label in label_order])[label] for label in wide.columns]
    wide.plot(kind="bar", stacked=True, ax=ax_bar, width=0.85, color=colors)
    ax_bar.set_title(title)
    ax_bar.set_xlabel("seurat_clusters")
    ax_bar.set_ylabel("fraction of cells")
    ax_bar.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, title=label_col, fontsize=8)

    right = outer[0, 1].subgridspec(rows, cluster_cols, wspace=0.12, hspace=0.35)
    cluster_values = obs[cluster_col].astype(str).to_numpy()
    for idx, cluster in enumerate(clusters):
        row = idx // cluster_cols
        col = idx % cluster_cols
        ax = fig.add_subplot(right[row, col])
        mask = cluster_values == str(cluster)
        ax.scatter(coords[:, 0], coords[:, 1], c="#d0d0d0", s=0.6, linewidths=0, rasterized=True)
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            c=color_for_index(idx),
            s=1.4,
            linewidths=0,
            rasterized=True,
        )
        ax.set_title(f"cluster {cluster} (n={int(mask.sum())})", fontsize=9)
        style_small_umap_axis(ax, coords)

    for idx in range(len(clusters), rows * cluster_cols):
        row = idx // cluster_cols
        col = idx % cluster_cols
        ax = fig.add_subplot(right[row, col])
        ax.axis("off")

    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_stacked_bar_fixed_labels(
    label_counts: pd.DataFrame,
    cluster_col: str,
    label_col: str,
    label_order: list[str],
    path: Path,
    title: str,
) -> None:
    wide = label_counts.pivot(index=cluster_col, columns=label_col, values="fraction").fillna(0.0)
    wide = wide.reindex(columns=[str(label) for label in label_order], fill_value=0.0)
    wide = wide.loc[sorted(wide.index.astype(str), key=natural_sort_key)]
    colors = [label_color_map([str(label) for label in label_order])[label] for label in wide.columns]
    fig, ax = plt.subplots(figsize=(max(9.0, 0.6 * wide.shape[0]), 5.8))
    wide.plot(kind="bar", stacked=True, ax=ax, width=0.85, color=colors)
    ax.set_title(title)
    ax.set_xlabel("seurat_clusters")
    ax.set_ylabel("fraction of cells")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, title=label_col, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_ridge_by_group(
    data: pd.DataFrame,
    value_col: str,
    group_col: str,
    path: Path,
    title: str,
    x_label: str,
) -> None:
    work = data[[value_col, group_col]].copy()
    work[value_col] = pd.to_numeric(work[value_col], errors="coerce")
    work[group_col] = work[group_col].astype(str)
    work = work.dropna(subset=[value_col, group_col])
    groups = sorted(work[group_col].unique(), key=natural_sort_key)
    if not groups:
        raise ValueError(f"No groups available for ridge plot: {group_col}")

    x_min = float(np.floor(work[value_col].min())) - 0.5
    x_max = float(np.ceil(work[value_col].max())) + 0.5
    x_grid = np.linspace(x_min, x_max, 300)
    fig, ax = plt.subplots(figsize=(9.0, max(5.0, 0.48 * len(groups))))

    for idx, group in enumerate(groups):
        values = work.loc[work[group_col] == group, value_col].to_numpy(dtype=float)
        values = values[np.isfinite(values)]
        if values.size == 0:
            continue
        if np.unique(values).size > 1:
            kde = gaussian_kde(values)
            density = kde(x_grid)
        else:
            center = values[0]
            density = np.exp(-0.5 * ((x_grid - center) / 0.18) ** 2)
        density = density / density.max() * 0.82 if density.max() > 0 else density
        baseline = idx
        color = color_for_index(idx)
        ax.fill_between(x_grid, baseline, baseline + density, color=color, alpha=0.82, linewidth=0)
        ax.plot(x_grid, baseline + density, color=color, linewidth=1.1)
        ax.plot([x_min, x_max], [baseline, baseline], color="#bdbdbd", linewidth=0.6)

    ax.set_yticks(np.arange(len(groups)))
    ax.set_yticklabels(groups)
    ax.set_ylim(-0.2, len(groups) - 0.1 + 0.9)
    ax.set_xlim(x_min, x_max)
    ax.set_xlabel(x_label)
    ax.set_ylabel(group_col)
    ax.set_title(title)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_long_score_table(obs: pd.DataFrame, score_cols: list[str], path: Path) -> None:
    id_cols = ["cell_id_for_join"]
    if "orig.ident" in obs.columns:
        id_cols.append("orig.ident")
    if "seurat_clusters" in obs.columns:
        id_cols.append("seurat_clusters")
    work = obs[id_cols + score_cols].copy()
    long = work.melt(
        id_vars=id_cols,
        value_vars=score_cols,
        var_name="score_column",
        value_name="score",
    )
    long["shi_label_token"] = long["score_column"].str.removeprefix("shi_seurat_full_prediction_score_")
    long.to_csv(path, sep="\t", index=False)


def main() -> None:
    repo_root = env_path("REPO_ROOT", "/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline")
    project_root = env_path("PROJECT_ROOT", "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder")
    conda_env_bin = env_path("CONDA_ENV_BIN", "/home/elcrespo/miniconda3/envs/mge-organoid-python/bin")

    results_dirname = os.environ.get("SHI_SEURAT_RESULTS_DIRNAME", "shi_reference_div30_seurat_label_transfer")
    run_label = os.environ.get("SHI_SEURAT_RUN_LABEL", "shi_reference_div30_seurat_label_transfer_v1")
    run_dir = project_root / "results" / results_dirname / run_label
    seurat_dir = run_dir / "seurat"
    table_dir = run_dir / "tables"
    plot_dir = run_dir / "plots"
    h5ad_dir = run_dir / "h5ad"
    for directory in [seurat_dir, table_dir, plot_dir, h5ad_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    reference_h5ad = env_path("SHI_SEURAT_REFERENCE_H5AD", str(project_root / "results/python_anndata/shi_2019_paper_qc.h5ad"))
    query_h5ad = env_path("SHI_SEURAT_QUERY_H5AD", str(project_root / "results/python_anndata/varela_div30.h5ad"))
    reference_seurat = env_path(
        "SHI_SEURAT_REFERENCE_RDS",
        str(project_root / "results/shi_2019_paper_qc/shi_2019_seurat.rds"),
    )
    query_seurat = env_path(
        "SHI_SEURAT_QUERY_RDS",
        str(project_root / "results/varela_this_paper/varela_this_paper_seurat.rds"),
    )
    table_s2_xlsx = env_path(
        "SHI_SEURAT_TABLE_S2_XLSX",
        str(project_root / "reference/shi_2021_tables_s2_to_s9/science.abj6641_table_s2.xlsx"),
    )
    rscript_bin = env_path("RSCRIPT", str(conda_env_bin / "Rscript"))
    dims = env_int("SHI_SEURAT_DIMS", 50)
    nfeatures = env_int("SHI_SEURAT_NFEATURES", 3000)
    npcs = env_int("SHI_SEURAT_NPCS", 50)
    seed = env_int("SHI_SEURAT_RANDOM_STATE", 0)
    min_match_fraction = float(os.environ.get("SHI_SEURAT_MIN_LABEL_MATCH_FRACTION", "0.95"))
    overwrite = env_bool("SHI_SEURAT_OVERWRITE", False)
    write_h5ad = env_bool("SHI_SEURAT_WRITE_H5AD", True)
    query_label = os.environ.get("SHI_SEURAT_QUERY_LABEL", "DIV30")
    query_slug = os.environ.get("SHI_SEURAT_QUERY_SLUG", safe_token(query_label).lower())

    labels_path = seurat_dir / "shi_reference_labels_for_seurat.tsv"
    prediction_path = seurat_dir / f"{query_slug}_shi_seurat_full_predictions.tsv.gz"
    score_path = seurat_dir / f"{query_slug}_shi_seurat_full_prediction_scores.tsv.gz"
    week_prediction_path = seurat_dir / f"{query_slug}_shi_seurat_full_week_predictions.tsv.gz"
    week_score_path = seurat_dir / f"{query_slug}_shi_seurat_full_week_prediction_scores.tsv.gz"

    log(f"Run label: {run_label}")
    log(f"Run dir: {run_dir}")
    label_order = prepare_reference_labels(
        reference_h5ad=reference_h5ad,
        table_s2_xlsx=table_s2_xlsx,
        labels_path=labels_path,
        table_dir=table_dir,
        min_match_fraction=min_match_fraction,
    )

    if overwrite or not prediction_path.exists() or not score_path.exists() or not week_prediction_path.exists() or not week_score_path.exists():
        run_seurat_transfer(
            repo_root=repo_root,
            rscript_bin=rscript_bin,
            reference_seurat=reference_seurat,
            query_seurat=query_seurat,
            labels_path=labels_path,
            seurat_dir=seurat_dir,
            output_prefix=query_slug,
            dims=dims,
            nfeatures=nfeatures,
            npcs=npcs,
            seed=seed,
        )
    else:
        log(f"Using existing Seurat predictions: {prediction_path}")

    log(f"Reading {query_label} AnnData for plotting: {query_h5ad}")
    adata = ad.read_h5ad(query_h5ad)
    predictions = load_predictions(prediction_path)
    aligned = align_predictions_to_adata(adata, predictions)
    score_cols = add_predictions_to_obs(adata, aligned, label_order=label_order)
    week_predictions = load_week_predictions(week_prediction_path)
    aligned_week = align_predictions_to_adata(adata, week_predictions)
    week_score_cols = add_week_predictions_to_obs(adata, aligned_week)
    validate_canonical_prediction_scores(
        adata.obs,
        score_cols,
        max_score_col=PREDICTION_SCORE_COL,
        uncertainty_col=UNCERTAINTY_SCORE_COL,
        context="Shi Seurat smoke label scores",
    )
    validate_canonical_prediction_scores(
        adata.obs,
        week_score_cols,
        max_score_col=WEEK_PREDICTION_SCORE_COL,
        uncertainty_col=WEEK_UNCERTAINTY_SCORE_COL,
        context="Shi Seurat smoke week scores",
    )
    if "cell_id" in adata.obs.columns:
        adata.obs["cell_id_for_join"] = adata.obs["cell_id"].astype(str)
    else:
        adata.obs["cell_id_for_join"] = adata.obs_names.astype(str)

    prefix = "shi_seurat_full"
    obs_cols = [
        "cell_id_for_join",
        "orig.ident",
        "seurat_clusters",
        f"{prefix}_predicted_shi_label",
        f"{prefix}_prediction_score",
        f"{prefix}_uncertainty_score",
        f"{prefix}_broad_region_class",
        f"{prefix}_developmental_class",
        f"{prefix}_predicted_shi_week_label",
        f"{prefix}_predicted_shi_week_numeric",
        f"{prefix}_week_prediction_score",
        f"{prefix}_week_uncertainty_score",
        f"{prefix}_expected_shi_week_numeric",
    ]
    obs_cols = [col for col in obs_cols if col in adata.obs.columns]
    adata.obs[obs_cols + score_cols + week_score_cols].to_csv(
        table_dir / f"{query_slug}_shi_seurat_label_transfer_obs.tsv.gz",
        sep="\t",
        index=True,
    )
    write_long_score_table(
        adata.obs,
        score_cols,
        table_dir / f"{query_slug}_shi_seurat_label_transfer_label_scores_long.tsv.gz",
    )
    write_long_score_table(
        adata.obs,
        week_score_cols,
        table_dir / f"{query_slug}_shi_seurat_label_transfer_week_scores_long.tsv.gz",
    )

    cluster_tables = summarize_predictions_by_cluster(
        adata.obs,
        label_col=f"{prefix}_predicted_shi_label",
        score_col=f"{prefix}_prediction_score",
        uncertainty_col=f"{prefix}_uncertainty_score",
        cluster_col="seurat_clusters",
        label_categories=label_order,
    )
    cluster_tables["cluster_summary"].to_csv(
        table_dir / f"{query_slug}_shi_seurat_full_cluster_summaries.tsv",
        sep="\t",
        index=False,
    )
    cluster_tables["label_counts"].to_csv(
        table_dir / f"{query_slug}_shi_seurat_full_label_fractions_by_cluster.tsv",
        sep="\t",
        index=False,
    )

    coords, umap_key = get_umap_coordinates(adata)
    log(f"Using UMAP key: {umap_key}")
    plot_umap_categorical(
        coords,
        adata.obs[f"{prefix}_predicted_shi_label"],
        f"{query_label} Seurat-transfer Shi label",
        plot_dir / f"{query_slug}_umap_shi_seurat_full_predicted_shi_label.png",
        point_size=1.0,
        category_order=label_order,
    )
    plot_predicted_label_with_overlay_panels(
        coords,
        adata.obs[f"{prefix}_predicted_shi_label"],
        label_order=label_order,
        path=plot_dir / f"{query_slug}_umap_shi_seurat_full_predicted_shi_label_with_subtype_overlays.png",
    )
    plot_umap_continuous(
        coords,
        adata.obs[f"{prefix}_prediction_score"],
        f"{query_label} Seurat-transfer max prediction score",
        plot_dir / f"{query_slug}_umap_shi_seurat_full_prediction_score.png",
        point_size=1.0,
        cmap="viridis",
    )

    panel_labels = ["MGE", "LGE", "CGE"]
    panel_cols = [f"{prefix}_prediction_score_{sanitize_score_label(label)}" for label in panel_labels]
    missing_panel = [col for col in panel_cols if col not in adata.obs.columns]
    if missing_panel:
        raise ValueError("Missing expected MGE/LGE/CGE score columns: " + ", ".join(missing_panel))
    plot_score_panel(
        coords,
        adata.obs,
        panel_cols,
        panel_labels,
        plot_dir / f"{query_slug}_umap_shi_seurat_full_mge_lge_cge_score_panel.png",
    )

    if "orig.ident" in adata.obs.columns:
        plot_overlaid_density_by_group(
            adata.obs,
            value_col=f"{prefix}_prediction_score",
            group_col="orig.ident",
            title=f"Seurat-transfer max prediction score by {query_label} sample",
            path=plot_dir / f"{query_slug}_shi_seurat_full_prediction_score_density_by_sample.png",
            bins=np.linspace(0, 1, 41),
            x_label="Seurat max prediction score",
        )

    plot_stacked_bar_fixed_labels(
        cluster_tables["label_counts"],
        cluster_col="seurat_clusters",
        label_col=f"{prefix}_predicted_shi_label",
        label_order=label_order,
        path=plot_dir / f"{query_slug}_shi_seurat_full_shi_label_stacked_bar_by_seurat_clusters.png",
        title=f"Seurat-transfer Shi label fractions by {query_label} seurat_clusters",
    )
    plot_stacked_bar_with_cluster_umaps(
        coords,
        adata.obs,
        cluster_tables["label_counts"],
        cluster_col="seurat_clusters",
        label_col=f"{prefix}_predicted_shi_label",
        label_order=label_order,
        path=plot_dir / f"{query_slug}_shi_seurat_full_shi_label_stacked_bar_by_seurat_clusters_with_cluster_umaps.png",
        title=f"Seurat-transfer Shi label fractions by {query_label} seurat_clusters",
    )

    if "orig.ident" in adata.obs.columns:
        plot_ridge_by_group(
            adata.obs,
            value_col=f"{prefix}_expected_shi_week_numeric",
            group_col="orig.ident",
            path=plot_dir / f"{query_slug}_shi_seurat_full_expected_shi_gw_ridge_by_sample.png",
            title=f"Seurat-transfer expected Shi gestational week by {query_label} sample",
            x_label="Expected Shi gestational week",
        )

    if write_h5ad:
        out_h5ad = h5ad_dir / f"{query_slug}_shi_seurat_label_transfer_predictions.h5ad"
        log(f"Writing annotated AnnData: {out_h5ad}")
        adata.write_h5ad(out_h5ad, compression="gzip")

    manifest = []
    for root in [seurat_dir, table_dir, plot_dir, h5ad_dir]:
        for path in sorted(root.glob("*")):
            if path.is_file():
                manifest.append(
                    {
                        "path": str(path),
                        "relative_path": str(path.relative_to(run_dir)),
                        "size_bytes": path.stat().st_size,
                    }
                )
    pd.DataFrame(manifest).to_csv(table_dir / "shi_seurat_label_transfer_output_manifest.tsv", sep="\t", index=False)

    completion = pd.DataFrame(
        [
            {
                "run_label": run_label,
                "status": "complete",
                "n_query_cells": int(adata.n_obs),
                "n_score_columns": int(len(score_cols)),
                "n_plots": len(list(plot_dir.glob("*.png"))),
            }
        ]
    )
    completion.to_csv(table_dir / "shi_seurat_label_transfer_complete.tsv", sep="\t", index=False)
    log("complete")


if __name__ == "__main__":
    main()
