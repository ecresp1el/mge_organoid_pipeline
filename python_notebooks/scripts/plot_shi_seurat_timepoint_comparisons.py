#!/usr/bin/env python3
"""Plot DIV30/DIV90 comparisons from completed Shi Seurat-transfer tables."""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

from mge_organoid_python.shi_label_transfer import natural_sort_key


SHI_GW_VALUES = np.array([9.0, 12.0, 13.0, 16.0, 18.0], dtype=float)
SHI_GW_LABELS = [f"GW{int(value):02d}" for value in SHI_GW_VALUES]
SHI_GW_EVEN_POSITIONS = np.arange(SHI_GW_VALUES.shape[0], dtype=float)


def env_path(name: str, default: str) -> Path:
    return Path(os.environ.get(name, default)).expanduser()


def env_float(name: str, default: float) -> float:
    return float(os.environ.get(name, str(default)))


def plot_axis_pad(ticks: list[float] | np.ndarray, fraction: float) -> tuple[float, float]:
    tick_values = np.asarray(ticks, dtype=float)
    if tick_values.size < 2:
        return 0.0, 0.0
    left_interval = tick_values[1] - tick_values[0]
    right_interval = tick_values[-1] - tick_values[-2]
    return (
        max(0.0, float(left_interval) * fraction),
        max(0.0, float(right_interval) * fraction),
    )


def load_obs(path: Path, timepoint: str) -> pd.DataFrame:
    cols = [
        "orig.ident",
        "seurat_clusters",
        "shi_seurat_full_predicted_shi_label",
        "shi_seurat_full_prediction_score",
        "shi_seurat_full_prediction_score_MGE",
        "shi_seurat_full_expected_shi_week_numeric",
    ]
    obs = pd.read_csv(path, sep="\t", usecols=cols)
    obs["timepoint"] = timepoint
    obs["sample"] = obs["orig.ident"].astype(str)
    obs["sample_timepoint"] = timepoint + " " + obs["sample"]
    obs["shi_seurat_full_expected_shi_week_numeric"] = pd.to_numeric(
        obs["shi_seurat_full_expected_shi_week_numeric"],
        errors="coerce",
    )
    obs["shi_seurat_full_prediction_score"] = pd.to_numeric(
        obs["shi_seurat_full_prediction_score"],
        errors="coerce",
    )
    obs["shi_seurat_full_prediction_score_MGE"] = pd.to_numeric(
        obs["shi_seurat_full_prediction_score_MGE"],
        errors="coerce",
    )
    return obs


def fixed_gw_axis(data: pd.Series) -> tuple[float, float, list[int]]:
    return 9.0, 18.0, [9, 12, 13, 16, 18]


def expected_gw_to_even_position(values: np.ndarray | pd.Series) -> np.ndarray:
    numeric = np.asarray(pd.to_numeric(pd.Series(values), errors="coerce"), dtype=float)
    return np.interp(
        numeric,
        SHI_GW_VALUES,
        SHI_GW_EVEN_POSITIONS,
        left=SHI_GW_EVEN_POSITIONS[0],
        right=SHI_GW_EVEN_POSITIONS[-1],
    )


def density_for_values(
    values: np.ndarray,
    x_grid: np.ndarray,
    bw_method: float,
) -> np.ndarray:
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.zeros_like(x_grid)
    if np.unique(values).size > 1:
        density = gaussian_kde(values, bw_method=bw_method)(x_grid)
    else:
        density = np.exp(-0.5 * ((x_grid - values[0]) / max(bw_method, 0.05)) ** 2)
    return density / density.max() if density.max() > 0 else density


def density_and_local_score(
    values: np.ndarray,
    scores: np.ndarray,
    x_grid: np.ndarray,
    bw_method: float,
) -> tuple[np.ndarray, np.ndarray]:
    mask = np.isfinite(values) & np.isfinite(scores)
    values = values[mask]
    scores = scores[mask]
    if values.size == 0:
        return np.zeros_like(x_grid), np.full_like(x_grid, np.nan)
    if np.unique(values).size == 1:
        density = np.exp(-0.5 * ((x_grid - values[0]) / max(bw_method, 0.05)) ** 2)
        local_score = np.full_like(x_grid, float(np.nanmean(scores)))
    else:
        kde = gaussian_kde(values, bw_method=bw_method)
        density = kde(x_grid)
        score_weights = np.clip(scores, 1e-8, None)
        score_kde = gaussian_kde(values, bw_method=bw_method, weights=score_weights)
        local_score = float(np.nanmean(score_weights)) * score_kde(x_grid) / np.maximum(density, 1e-12)
    density = density / density.max() if density.max() > 0 else density
    return density, np.clip(local_score, 0.0, 1.0)


def add_score_colored_ridge(
    ax: plt.Axes,
    x_grid: np.ndarray,
    baseline: float,
    density: np.ndarray,
    score: np.ndarray,
    cmap: str,
    norm: Normalize,
) -> PolyCollection:
    vertices = []
    segment_scores = []
    for idx in range(x_grid.size - 1):
        vertices.append(
            [
                (x_grid[idx], baseline),
                (x_grid[idx], baseline + density[idx]),
                (x_grid[idx + 1], baseline + density[idx + 1]),
                (x_grid[idx + 1], baseline),
            ]
        )
        segment_scores.append(float(np.nanmean(score[idx : idx + 2])))
    collection = PolyCollection(vertices, cmap=cmap, norm=norm, edgecolors="none", alpha=0.9)
    collection.set_array(np.asarray(segment_scores, dtype=float))
    ax.add_collection(collection)
    return collection


def plot_expected_gw_ridge_by_sample(
    data: pd.DataFrame,
    path: Path,
    title: str,
    bw_method: float,
    axis_pad_fraction: float,
    even_gw_axis: bool = False,
) -> None:
    value_col = "shi_seurat_full_expected_shi_week_numeric"
    if even_gw_axis:
        x_min, x_max, ticks = 0.0, float(SHI_GW_EVEN_POSITIONS[-1]), SHI_GW_EVEN_POSITIONS.tolist()
    else:
        x_min, x_max, ticks = fixed_gw_axis(data[value_col])
    left_pad, right_pad = plot_axis_pad(ticks, axis_pad_fraction)
    x_plot_min = x_min - left_pad
    x_plot_max = x_max + right_pad
    x_grid = np.linspace(x_plot_min, x_plot_max, 350)
    groups = []
    for timepoint in ["DIV30", "DIV90"]:
        samples = sorted(
            data.loc[data["timepoint"] == timepoint, "sample"].dropna().astype(str).unique(),
            key=natural_sort_key,
        )
        groups.extend([(timepoint, sample, f"{timepoint} {sample}") for sample in samples])

    fig, ax = plt.subplots(figsize=(5.2, max(5.6, 0.5 * len(groups))))
    colors = {"DIV30": "#377eb8", "DIV90": "#e6550d"}
    for idx, (timepoint, sample, label) in enumerate(groups):
        values = data.loc[
            (data["timepoint"] == timepoint) & (data["sample"] == sample),
            value_col,
        ].to_numpy(dtype=float)
        if even_gw_axis:
            values = expected_gw_to_even_position(values)
        density = density_for_values(values, x_grid, bw_method=bw_method) * 0.82
        ax.fill_between(
            x_grid,
            idx,
            idx + density,
            color=colors[timepoint],
            alpha=0.78,
            linewidth=0,
        )
        ax.plot(
            x_grid,
            idx + density,
            color=colors[timepoint],
            linewidth=1.05,
        )
        ax.plot([x_plot_min, x_plot_max], [idx, idx], color="#c7c7c7", linewidth=0.55)

    ax.set_yticks(np.arange(len(groups)))
    ax.set_yticklabels([label for _, _, label in groups])
    ax.set_ylim(-0.2, len(groups) - 0.1 + 0.9)
    ax.set_xlim(x_plot_min, x_plot_max)
    ax.set_xticks(ticks)
    if even_gw_axis:
        ax.set_xticklabels(SHI_GW_LABELS)
        ax.set_xlabel("Shi GW (even spacing)")
    else:
        ax.set_xlabel("Expected Shi GW")
    ax.set_ylabel("sample")
    ax.set_title(title)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_expected_gw_score_colored_ridge_by_sample(
    data: pd.DataFrame,
    path: Path,
    title: str,
    score_col: str,
    score_label: str,
    bw_method: float,
    axis_pad_fraction: float,
    even_gw_axis: bool = False,
) -> None:
    value_col = "shi_seurat_full_expected_shi_week_numeric"
    if even_gw_axis:
        x_min, x_max, ticks = 0.0, float(SHI_GW_EVEN_POSITIONS[-1]), SHI_GW_EVEN_POSITIONS.tolist()
    else:
        x_min, x_max, ticks = fixed_gw_axis(data[value_col])
    left_pad, right_pad = plot_axis_pad(ticks, axis_pad_fraction)
    x_plot_min = x_min - left_pad
    x_plot_max = x_max + right_pad
    x_grid = np.linspace(x_plot_min, x_plot_max, 350)
    groups = []
    for timepoint in ["DIV30", "DIV90"]:
        samples = sorted(
            data.loc[data["timepoint"] == timepoint, "sample"].dropna().astype(str).unique(),
            key=natural_sort_key,
        )
        groups.extend([(timepoint, sample, f"{timepoint} {sample}") for sample in samples])

    fig, ax = plt.subplots(figsize=(5.2, max(5.6, 0.5 * len(groups))))
    norm = Normalize(vmin=0.0, vmax=1.0)
    mappable = None
    for idx, (timepoint, sample, label) in enumerate(groups):
        group = data.loc[(data["timepoint"] == timepoint) & (data["sample"] == sample), [value_col, score_col]]
        values = group[value_col].to_numpy(dtype=float)
        scores = group[score_col].to_numpy(dtype=float)
        if even_gw_axis:
            values = expected_gw_to_even_position(values)
        density, local_score = density_and_local_score(values, scores, x_grid, bw_method=bw_method)
        density = density * 0.82
        ax.plot([x_plot_min, x_plot_max], [idx, idx], color="#c7c7c7", linewidth=0.55)
        mappable = add_score_colored_ridge(
            ax,
            x_grid=x_grid,
            baseline=float(idx),
            density=density,
            score=local_score,
            cmap="viridis",
            norm=norm,
        )
        ax.plot(x_grid, idx + density, color="#3f3f3f", linewidth=0.75)

    ax.set_yticks(np.arange(len(groups)))
    ax.set_yticklabels([label for _, _, label in groups])
    ax.set_ylim(-0.2, len(groups) - 0.1 + 0.9)
    ax.set_xlim(x_plot_min, x_plot_max)
    ax.set_xticks(ticks)
    if even_gw_axis:
        ax.set_xticklabels(SHI_GW_LABELS)
        ax.set_xlabel("Shi GW (even spacing)")
    else:
        ax.set_xlabel("Expected Shi GW")
    ax.set_ylabel("sample")
    ax.set_title(title)
    if mappable is not None:
        colorbar = fig.colorbar(mappable, ax=ax, fraction=0.032, pad=0.02)
        colorbar.set_label(score_label)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_expected_gw_density_by_timepoint(
    data: pd.DataFrame,
    path: Path,
    title: str,
    bw_method: float,
    axis_pad_fraction: float,
    even_gw_axis: bool = False,
) -> None:
    value_col = "shi_seurat_full_expected_shi_week_numeric"
    if even_gw_axis:
        x_min, x_max, ticks = 0.0, float(SHI_GW_EVEN_POSITIONS[-1]), SHI_GW_EVEN_POSITIONS.tolist()
    else:
        x_min, x_max, ticks = fixed_gw_axis(data[value_col])
    left_pad, right_pad = plot_axis_pad(ticks, axis_pad_fraction)
    x_plot_min = x_min - left_pad
    x_plot_max = x_max + right_pad
    x_grid = np.linspace(x_plot_min, x_plot_max, 350)
    colors = {"DIV30": "#377eb8", "DIV90": "#e6550d"}
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    for timepoint in ["DIV30", "DIV90"]:
        values = data.loc[data["timepoint"] == timepoint, value_col].to_numpy(dtype=float)
        if even_gw_axis:
            values = expected_gw_to_even_position(values)
        density = density_for_values(values, x_grid, bw_method=bw_method)
        ax.fill_between(x_grid, density, color=colors[timepoint], alpha=0.25)
        ax.plot(
            x_grid,
            density,
            color=colors[timepoint],
            linewidth=2.0,
            label=f"{timepoint} (n={values.size})",
        )
    ax.set_xlim(x_plot_min, x_plot_max)
    ax.set_xticks(ticks)
    if even_gw_axis:
        ax.set_xticklabels(SHI_GW_LABELS)
        ax.set_xlabel("Shi GW (even spacing)")
    else:
        ax.set_xlabel("Expected Shi GW")
    ax.set_ylabel("density")
    ax.set_title(title)
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_summary_tables(data: pd.DataFrame, table_dir: Path, suffix: str = "") -> None:
    value_col = "shi_seurat_full_expected_shi_week_numeric"
    summary = (
        data.groupby(["timepoint", "sample"], observed=True)
        .agg(
            n_cells=(value_col, "count"),
            expected_gw_mean=(value_col, "mean"),
            expected_gw_median=(value_col, "median"),
            expected_gw_q25=(value_col, lambda x: float(np.nanquantile(x, 0.25))),
            expected_gw_q75=(value_col, lambda x: float(np.nanquantile(x, 0.75))),
            mean_prediction_score=("shi_seurat_full_prediction_score", "mean"),
            mean_mge_score=("shi_seurat_full_prediction_score_MGE", "mean"),
        )
        .reset_index()
    )
    summary.to_csv(table_dir / f"div30_div90_expected_shi_gw_summary_by_sample{suffix}.tsv", sep="\t", index=False)

    label_counts = (
        data.groupby(["timepoint", "shi_seurat_full_predicted_shi_label"], observed=True)
        .size()
        .reset_index(name="n_cells")
    )
    totals = data.groupby("timepoint", observed=True).size().rename("n_cells_timepoint").reset_index()
    label_counts = label_counts.merge(totals, on="timepoint")
    label_counts["fraction"] = label_counts["n_cells"] / label_counts["n_cells_timepoint"]
    label_counts.to_csv(table_dir / f"div30_div90_shi_label_counts_by_timepoint{suffix}.tsv", sep="\t", index=False)


def make_expected_gw_plots(
    data: pd.DataFrame,
    plot_dir: Path,
    table_dir: Path,
    suffix: str,
    title_suffix: str,
    score_col: str,
    score_label: str,
    bw_method: float,
    axis_pad_fraction: float,
) -> None:
    write_summary_tables(data, table_dir, suffix=suffix)
    plot_expected_gw_ridge_by_sample(
        data,
        plot_dir / f"div30_div90_shi_seurat_full_expected_shi_gw_ridge_by_sample_shared_axis{suffix}.png",
        f"Expected Shi GW by sample{title_suffix}",
        bw_method=bw_method,
        axis_pad_fraction=axis_pad_fraction,
    )
    plot_expected_gw_ridge_by_sample(
        data,
        plot_dir / f"div30_div90_shi_seurat_full_expected_shi_gw_ridge_by_sample_shared_axis_even_gw_spacing{suffix}.png",
        f"Expected Shi GW by sample{title_suffix}; even spacing",
        bw_method=bw_method,
        axis_pad_fraction=axis_pad_fraction,
        even_gw_axis=True,
    )
    plot_expected_gw_score_colored_ridge_by_sample(
        data,
        plot_dir / f"div30_div90_shi_seurat_full_expected_shi_gw_ridge_by_sample_shared_axis_prediction_score_color{suffix}.png",
        f"Expected Shi GW by sample{title_suffix}; score color",
        score_col=score_col,
        score_label=score_label,
        bw_method=bw_method,
        axis_pad_fraction=axis_pad_fraction,
    )
    plot_expected_gw_score_colored_ridge_by_sample(
        data,
        plot_dir
        / f"div30_div90_shi_seurat_full_expected_shi_gw_ridge_by_sample_shared_axis_even_gw_spacing_prediction_score_color{suffix}.png",
        f"Expected Shi GW by sample{title_suffix}; score color, even spacing",
        score_col=score_col,
        score_label=score_label,
        bw_method=bw_method,
        axis_pad_fraction=axis_pad_fraction,
        even_gw_axis=True,
    )
    plot_expected_gw_density_by_timepoint(
        data,
        plot_dir / f"div30_div90_shi_seurat_full_expected_shi_gw_density_by_timepoint_shared_axis{suffix}.png",
        f"Expected Shi GW density{title_suffix}",
        bw_method=bw_method,
        axis_pad_fraction=axis_pad_fraction,
    )
    plot_expected_gw_density_by_timepoint(
        data,
        plot_dir / f"div30_div90_shi_seurat_full_expected_shi_gw_density_by_timepoint_shared_axis_even_gw_spacing{suffix}.png",
        f"Expected Shi GW density{title_suffix}; even spacing",
        bw_method=bw_method,
        axis_pad_fraction=axis_pad_fraction,
        even_gw_axis=True,
    )


def main() -> None:
    project_root = env_path("PROJECT_ROOT", "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder")
    div30_obs = env_path(
        "SHI_SEURAT_DIV30_OBS",
        str(
            project_root
            / "results/shi_reference_div30_seurat_label_transfer/shi_reference_div30_seurat_label_transfer_v1/tables/div30_shi_seurat_label_transfer_obs.tsv.gz"
        ),
    )
    div90_obs = env_path(
        "SHI_SEURAT_DIV90_OBS",
        str(
            project_root
            / "results/shi_reference_div90_seurat_label_transfer/shi_reference_div90_seurat_label_transfer_v1/tables/div90_shi_seurat_label_transfer_obs.tsv.gz"
        ),
    )
    results_dirname = os.environ.get(
        "SHI_SEURAT_COMPARISON_RESULTS_DIRNAME",
        "shi_reference_div30_div90_seurat_label_transfer_comparison",
    )
    run_label = os.environ.get("SHI_SEURAT_COMPARISON_RUN_LABEL", "shi_reference_div30_div90_seurat_comparison_v1")
    kde_bw_method = env_float("SHI_SEURAT_COMPARISON_KDE_BW_METHOD", 0.08)
    kde_axis_pad_fraction = env_float("SHI_SEURAT_COMPARISON_KDE_AXIS_PAD_FRACTION", 0.25)
    run_dir = project_root / "results" / results_dirname / run_label
    plot_dir = run_dir / "plots"
    table_dir = run_dir / "tables"
    plot_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)

    data = pd.concat(
        [
            load_obs(div30_obs, "DIV30"),
            load_obs(div90_obs, "DIV90"),
        ],
        ignore_index=True,
    )
    data.to_csv(table_dir / "div30_div90_shi_seurat_label_transfer_obs_subset.tsv.gz", sep="\t", index=False)
    make_expected_gw_plots(
        data,
        plot_dir,
        table_dir,
        suffix="",
        title_suffix="",
        score_col="shi_seurat_full_prediction_score",
        score_label="prediction score",
        bw_method=kde_bw_method,
        axis_pad_fraction=kde_axis_pad_fraction,
    )

    mge = data.loc[data["shi_seurat_full_predicted_shi_label"].astype(str) == "MGE"].copy()
    mge.to_csv(table_dir / "div30_div90_shi_seurat_label_transfer_obs_subset_predicted_mge.tsv.gz", sep="\t", index=False)
    make_expected_gw_plots(
        mge,
        plot_dir,
        table_dir,
        suffix="_predicted_mge",
        title_suffix=" (MGE)",
        score_col="shi_seurat_full_prediction_score_MGE",
        score_label="MGE prediction score",
        bw_method=kde_bw_method,
        axis_pad_fraction=kde_axis_pad_fraction,
    )
    completion = pd.DataFrame(
        [
            {
                "run_label": run_label,
                "status": "complete",
                "n_cells": int(data.shape[0]),
                "n_plots": len(list(plot_dir.glob("*.png"))),
                "kde_bw_method": kde_bw_method,
                "kde_axis_pad_fraction": kde_axis_pad_fraction,
                "numeric_gw_range": "9,18",
                "even_gw_range": "0,4",
            }
        ]
    )
    completion.to_csv(table_dir / "div30_div90_shi_seurat_comparison_complete.tsv", sep="\t", index=False)
    print(f"complete: {run_dir}")


if __name__ == "__main__":
    main()
