#!/usr/bin/env python3
"""Render cross-study receptor capture and in-silico MACS enrichment figure."""

from __future__ import annotations

import argparse
from datetime import datetime
import hashlib
from pathlib import Path
import shutil
import subprocess

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
        "font.size": 8,
        "axes.linewidth": 0.75,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }
)

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
import numpy as np
import pandas as pd

from mge_organoid_python.cross_study_marker_expression import (
    _apply_internal_umap_plot_filters,
    add_marker_plot_coordinates,
    default_cross_study_marker_specs,
)
from render_cross_study_lhx6_expression_final_panel import (
    add_sample_metadata,
    ordered_studies,
)


PROJECT_ROOT_DEFAULT = Path("/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder")
BASE_RUN_DEFAULT = "cross_study_marker_expression_v12_pv_precursors_final_candidate_plus_vipr2"
RECEPTOR_RUN_DEFAULT = "cross_study_marker_expression_v13_receptor_macs"
FINAL_FOLDER_DEFAULT = "fig_cross_study_receptor_macs_enrichment_v1_candidate"
OUTPUT_STEM = "cross_study_receptor_macs_enrichment"

CAPTURE_GENES = ["ERBB4", "CXCR4", "PLXNA2"]
STATE_GENES = ["GAD2", "LHX6", "LHX8", "NKX2-1", "MKI67", "SP8", "PAX6", "NEUROD2"]
BASE_GENES = ["ERBB4", *STATE_GENES]
RECEPTOR_ONLY_GENES = ["CXCR4", "PLXNA2"]
GENE_COLORS = {"ERBB4": "#3B7EA1", "CXCR4": "#8E6BBE", "PLXNA2": "#D05A76"}

STATE_ORDER = [
    "LHX6+ interneuron-like",
    "LHX8+ ventral-like",
    "NKX2.1+ progenitor-like",
    "Dorsal/neurogenic-like",
    "Cycling",
    "Other / undetermined",
]
STATE_COLORS = {
    "LHX6+ interneuron-like": "#2878B5",
    "LHX8+ ventral-like": "#D84A73",
    "NKX2.1+ progenitor-like": "#1B9E77",
    "Dorsal/neurogenic-like": "#E6A23C",
    "Cycling": "#7A5195",
    "Other / undetermined": "#C7CBD1",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT_DEFAULT)
    parser.add_argument("--base-run-label", default=BASE_RUN_DEFAULT)
    parser.add_argument("--receptor-run-label", default=RECEPTOR_RUN_DEFAULT)
    parser.add_argument("--final-dir", type=Path, default=None)
    parser.add_argument("--positive-threshold", type=float, default=0.5)
    parser.add_argument("--capture-sensitivity", type=float, default=0.80)
    parser.add_argument("--background-retention", type=float, default=0.05)
    parser.add_argument("--max-umap-cells-per-study", type=int, default=15000)
    parser.add_argument("--random-state", type=int, default=1729)
    parser.add_argument("--dpi", type=int, default=400)
    parser.add_argument("--fresh", action="store_true")
    return parser.parse_args()


def receptor_table_path(project_root: Path, run_label: str, study_id: str) -> Path:
    return (
        project_root
        / "results"
        / "cross_study_marker_expression"
        / run_label
        / "tables"
        / "per_study"
        / f"{study_id}_marker_expression.tsv.gz"
    )


def load_data(project_root: Path, base_run: str, receptor_run: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    specs = default_cross_study_marker_specs(project_root)
    base_parts = []
    for spec in specs:
        if not spec.include_in_first_plot:
            continue
        path = receptor_table_path(project_root, base_run, spec.study_id)
        if not path.exists():
            raise FileNotFoundError(f"Missing established marker table: {path}")
        usecols = ["cell_id", "study_id", "study_label", "sample", "cluster", "umap_1", "umap_2", *BASE_GENES]
        base_parts.append(
            pd.read_csv(
                path,
                sep="\t",
                usecols=usecols,
                dtype={"cell_id": str, "study_id": str, "study_label": str, "sample": str, "cluster": str},
                low_memory=False,
            )
        )
    base = pd.concat(base_parts, ignore_index=True)
    receptor_parts = []
    availability = []
    for spec in specs:
        if not spec.include_in_first_plot:
            continue
        path = receptor_table_path(project_root, receptor_run, spec.study_id)
        if not path.exists():
            raise FileNotFoundError(f"Missing receptor extraction table: {path}")
        header = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
        missing = [gene for gene in RECEPTOR_ONLY_GENES if gene not in header]
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(missing)}")
        part = pd.read_csv(
            path,
            sep="\t",
            usecols=["cell_id", "study_id", *RECEPTOR_ONLY_GENES],
            dtype={"cell_id": str, "study_id": str},
            low_memory=False,
        )
        for gene in RECEPTOR_ONLY_GENES:
            finite = pd.to_numeric(part[gene], errors="coerce").notna()
            availability.append(
                {
                    "study_id": spec.study_id,
                    "study_label": spec.study_label,
                    "gene": gene,
                    "source_table": str(path),
                    "n_receptor_table_cells": int(part.shape[0]),
                    "n_finite_expression": int(finite.sum()),
                    "expression_available": bool(finite.any()),
                }
            )
        receptor_parts.append(part)
    receptor = pd.concat(receptor_parts, ignore_index=True)
    if receptor.duplicated(["study_id", "cell_id"]).any():
        raise ValueError("Receptor extraction contains duplicated study/cell keys.")
    merged = base.merge(receptor, on=["study_id", "cell_id"], how="left", validate="one_to_one")
    merge_qc = (
        merged.groupby(["study_id", "study_label"], observed=True)
        .agg(
            n_base_table_cells=("cell_id", "size"),
            n_cxcr4_mapped=("CXCR4", lambda x: pd.to_numeric(x, errors="coerce").notna().sum()),
            n_plxna2_mapped=("PLXNA2", lambda x: pd.to_numeric(x, errors="coerce").notna().sum()),
        )
        .reset_index()
    )
    if (merge_qc[["n_cxcr4_mapped", "n_plxna2_mapped"]].min(axis=1) == 0).any():
        failed = merge_qc.loc[merge_qc[["n_cxcr4_mapped", "n_plxna2_mapped"]].min(axis=1) == 0, "study_id"]
        raise ValueError("No receptor values mapped for: " + ", ".join(failed))
    return merged, pd.DataFrame(availability).merge(merge_qc, on=["study_id", "study_label"], how="left")


def positive(data: pd.DataFrame, gene: str, threshold: float) -> pd.Series:
    values = pd.to_numeric(data[gene], errors="coerce")
    return values.ge(threshold) & values.notna()


def assign_states(data: pd.DataFrame, threshold: float) -> pd.DataFrame:
    out = data.copy()
    gad2 = positive(out, "GAD2", threshold)
    lhx6 = positive(out, "LHX6", threshold)
    lhx8 = positive(out, "LHX8", threshold)
    nkx21 = positive(out, "NKX2-1", threshold)
    cycling = positive(out, "MKI67", threshold)
    dorsal = positive(out, "SP8", threshold) | positive(out, "PAX6", threshold) | positive(out, "NEUROD2", threshold)
    state = np.full(out.shape[0], "Other / undetermined", dtype=object)
    state[dorsal.to_numpy()] = "Dorsal/neurogenic-like"
    state[nkx21.to_numpy()] = "NKX2.1+ progenitor-like"
    state[(gad2 & lhx8).to_numpy()] = "LHX8+ ventral-like"
    state[(gad2 & lhx6).to_numpy()] = "LHX6+ interneuron-like"
    state[cycling.to_numpy()] = "Cycling"
    out["cell_state"] = pd.Categorical(state, categories=STATE_ORDER, ordered=True)
    out["independent_target"] = gad2 & (lhx6 | lhx8) & ~cycling
    return out


def model_weights(marker_positive: pd.Series, sensitivity: float, background: float) -> np.ndarray:
    return np.where(marker_positive.to_numpy(dtype=bool), sensitivity, background).astype(float)


def summarize_capture(
    data: pd.DataFrame,
    threshold: float,
    sensitivity: float,
    background: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    summary_rows = []
    composition_rows = []
    sample_rows = []
    for (study_id, study_label), group in data.groupby(["study_id", "study_label"], sort=False, observed=True):
        target = group["independent_target"].to_numpy(dtype=bool)
        baseline = float(target.mean()) if target.size else np.nan
        for gene in CAPTURE_GENES:
            finite = pd.to_numeric(group[gene], errors="coerce").notna().to_numpy()
            if not finite.any():
                continue
            valid = group.loc[finite].copy()
            target_valid = valid["independent_target"].to_numpy(dtype=bool)
            marker_pos = positive(valid, gene, threshold)
            weights = model_weights(marker_pos, sensitivity, background)
            expected_retained = float(weights.sum())
            target_weight = float(weights[target_valid].sum())
            baseline_valid = float(target_valid.mean()) if target_valid.size else np.nan
            purity = target_weight / expected_retained if expected_retained > 0 else np.nan
            recovery = target_weight / float(target_valid.sum()) if target_valid.sum() else np.nan
            summary_rows.append(
                {
                    "study_id": study_id,
                    "study_label": study_label,
                    "capture_gene": gene,
                    "positive_threshold_log1p_cp10k": threshold,
                    "capture_sensitivity": sensitivity,
                    "background_retention": background,
                    "n_cells": int(valid.shape[0]),
                    "n_target_cells": int(target_valid.sum()),
                    "n_marker_positive": int(marker_pos.sum()),
                    "pct_marker_positive": float(marker_pos.mean() * 100.0),
                    "baseline_target_fraction": baseline_valid,
                    "expected_post_sort_yield_fraction": expected_retained / valid.shape[0],
                    "expected_target_recovery": recovery,
                    "expected_post_sort_target_purity": purity,
                    "expected_target_enrichment_fold": purity / baseline_valid if baseline_valid > 0 else np.nan,
                    "target_definition": "GAD2 >= threshold AND (LHX6 >= threshold OR LHX8 >= threshold) AND MKI67 < threshold",
                }
            )
            for state in STATE_ORDER:
                state_mask = valid["cell_state"].astype(str).eq(state).to_numpy()
                state_weight = float(weights[state_mask].sum())
                composition_rows.append(
                    {
                        "study_id": study_id,
                        "study_label": study_label,
                        "capture_gene": gene,
                        "cell_state": state,
                        "expected_retained_weight": state_weight,
                        "expected_post_sort_fraction": state_weight / expected_retained if expected_retained else np.nan,
                    }
                )

            for sample, sample_group in valid.groupby("sample", sort=False, observed=True):
                sample_target = sample_group["independent_target"].to_numpy(dtype=bool)
                sample_pos = positive(sample_group, gene, threshold)
                sample_weights = model_weights(sample_pos, sensitivity, background)
                retained = float(sample_weights.sum())
                target_retained = float(sample_weights[sample_target].sum())
                sample_rows.append(
                    {
                        "study_id": study_id,
                        "study_label": study_label,
                        "sample": str(sample),
                        "capture_gene": gene,
                        "n_cells": int(sample_group.shape[0]),
                        "pct_marker_positive": float(sample_pos.mean() * 100.0),
                        "baseline_target_fraction": float(sample_target.mean()),
                        "expected_post_sort_yield_fraction": retained / sample_group.shape[0],
                        "expected_post_sort_target_purity": target_retained / retained if retained else np.nan,
                    }
                )
    return pd.DataFrame(summary_rows), pd.DataFrame(composition_rows), pd.DataFrame(sample_rows)


def sensitivity_summary(data: pd.DataFrame, threshold: float) -> pd.DataFrame:
    rows = []
    for sensitivity in [0.60, 0.80, 0.95]:
        for background in [0.01, 0.05, 0.10]:
            summary, _, _ = summarize_capture(data, threshold, sensitivity, background)
            summary["scenario"] = f"sensitivity_{sensitivity:g}_background_{background:g}"
            rows.append(summary)
    return pd.concat(rows, ignore_index=True)


def cross_study_rank_summary(summary: pd.DataFrame, minimum_input_target_fraction: float = 0.01) -> pd.DataFrame:
    eligible = summary.loc[summary["baseline_target_fraction"].ge(minimum_input_target_fraction)].copy()
    rank = (
        eligible.groupby("capture_gene", observed=True)
        .agg(
            n_eligible_studies=("study_id", "nunique"),
            median_pct_marker_positive=("pct_marker_positive", "median"),
            median_expected_target_recovery=("expected_target_recovery", "median"),
            median_expected_post_sort_target_purity=("expected_post_sort_target_purity", "median"),
            median_expected_target_enrichment_fold=("expected_target_enrichment_fold", "median"),
            minimum_expected_target_enrichment_fold=("expected_target_enrichment_fold", "min"),
            maximum_expected_target_enrichment_fold=("expected_target_enrichment_fold", "max"),
        )
        .reset_index()
    )
    rank["target_recovery_rank"] = rank["median_expected_target_recovery"].rank(method="min", ascending=False).astype(int)
    rank["target_enrichment_rank"] = rank["median_expected_target_enrichment_fold"].rank(method="min", ascending=False).astype(int)
    rank["eligibility_rule"] = f"baseline independent-target fraction >= {minimum_input_target_fraction:g}"
    return rank.sort_values(["target_recovery_rank", "target_enrichment_rank"], kind="stable").reset_index(drop=True)


def downsample_umap(data: pd.DataFrame, max_cells: int, random_state: int) -> pd.DataFrame:
    if max_cells <= 0:
        return data
    parts = []
    rng = np.random.default_rng(random_state)
    for _, group in data.groupby("study_id", sort=False, observed=True):
        if group.shape[0] <= max_cells:
            parts.append(group)
        else:
            chosen = rng.choice(group.index.to_numpy(), size=max_cells, replace=False)
            parts.append(group.loc[np.sort(chosen)])
    return pd.concat(parts, ignore_index=True)


def sample_palette(samples: list[str]) -> dict[str, tuple[float, float, float, float]]:
    cmap = plt.get_cmap("turbo")
    positions = np.linspace(0.05, 0.95, max(len(samples), 2))
    return {sample: cmap(positions[idx]) for idx, sample in enumerate(samples)}


def plot_sample_umap(ax: plt.Axes, data: pd.DataFrame, sample_key_rows: list[dict[str, object]]) -> None:
    sample_order = (
        data[["sample", "sample_display_label", "sample_age_order", "sample_order_within_study"]]
        .drop_duplicates()
        .sort_values(["sample_age_order", "sample_order_within_study", "sample_display_label"], kind="stable")
    )
    samples = sample_order["sample"].astype(str).tolist()
    colors = sample_palette(samples)
    for idx, row in enumerate(sample_order.itertuples(index=False), start=1):
        sample = str(row.sample)
        subset = data.loc[data["sample"].astype(str).eq(sample)]
        ax.scatter(
            subset["UMAP1_plot"],
            subset["UMAP2_plot"],
            s=0.65,
            color=colors[sample],
            linewidths=0,
            alpha=0.84,
            rasterized=True,
            label=f"S{idx}",
        )
        sample_key_rows.append(
            {
                "study_id": str(data["study_id"].iloc[0]),
                "study_label": str(data["study_label"].iloc[0]),
                "sample_code": f"S{idx}",
                "sample": sample,
                "sample_display_label": str(row.sample_display_label),
                "color_hex": matplotlib.colors.to_hex(colors[sample]),
                "n_cells_after_filter": int((data["sample"].astype(str) == sample).sum()),
            }
        )
    ax.set_title(str(data["study_label"].iloc[0]), fontsize=9.2, fontweight="bold", pad=3)
    ax.text(0.5, 0.99, f"displayed n = {len(data):,}; {len(samples)} samples", transform=ax.transAxes, ha="center", va="top", fontsize=6.4, color="#555B61")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("equal", adjustable="box")
    for spine in ax.spines.values():
        spine.set_visible(False)
    handles, labels = ax.get_legend_handles_labels()
    ncol = 2 if len(samples) > 7 else 1
    ax.legend(
        handles,
        labels,
        loc="upper left",
        bbox_to_anchor=(0.00, -0.01),
        frameon=False,
        fontsize=5.5,
        ncol=ncol,
        handletextpad=0.25,
        columnspacing=0.6,
        markerscale=3.2,
        borderaxespad=0,
    )


def plot_capture_matrix(ax: plt.Axes, summary: pd.DataFrame, studies: pd.DataFrame) -> None:
    y_map = {study: idx for idx, study in enumerate(studies["study_id"].astype(str))}
    x_map = {gene: idx for idx, gene in enumerate(CAPTURE_GENES)}
    cmap = LinearSegmentedColormap.from_list("purity", ["#F2F4F7", "#8CB9D5", "#155A85"])
    finite_purity = summary["expected_post_sort_target_purity"].to_numpy(dtype=float)
    vmax = max(0.30, float(np.nanmax(finite_purity)))
    norm = Normalize(vmin=0, vmax=vmax)
    recovery = summary["expected_target_recovery"].to_numpy(dtype=float)
    sizes = 45 + 430 * np.clip(recovery, 0, 1)
    xs = summary["capture_gene"].map(x_map).to_numpy(dtype=float)
    ys = summary["study_id"].astype(str).map(y_map).to_numpy(dtype=float)
    ax.scatter(xs, ys, s=sizes, c=finite_purity, cmap=cmap, norm=norm, edgecolors="#30353A", linewidths=0.45)
    for x, y, fold, n_target in zip(
        xs,
        ys,
        summary["expected_target_enrichment_fold"],
        summary["n_target_cells"],
        strict=True,
    ):
        label = "low n" if int(n_target) < 25 else "NA" if not np.isfinite(fold) else f"{fold:.1f}×"
        ax.text(x, y, label, ha="center", va="center", fontsize=7.1, fontweight="bold", color="#111416")
    ax.set_xticks(range(len(CAPTURE_GENES)), [rf"$\mathit{{{gene}}}$" for gene in CAPTURE_GENES], fontsize=10.5)
    labels = []
    baseline = summary.groupby("study_id", observed=True)["baseline_target_fraction"].first()
    for row in studies.itertuples(index=False):
        baseline_pct = 100 * baseline.get(row.study_id, np.nan)
        baseline_text = "<0.1%" if 0 < baseline_pct < 0.1 else f"{baseline_pct:.1f}%"
        labels.append(f"{row.study_label}\ninput target {baseline_text}")
    ax.set_yticks(range(len(studies)), labels, fontsize=7.4)
    ax.set_xlim(-0.65, len(CAPTURE_GENES) - 0.35)
    ax.set_ylim(len(studies) - 0.45, -0.55)
    ax.grid(color="#E3E6E8", linewidth=0.55)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(
        0,
        1.015,
        "Dot area = target recovery; color = post-sort target purity; label = enrichment over input",
        transform=ax.transAxes,
        fontsize=7.4,
        color="#555B61",
        ha="left",
        va="bottom",
    )
    cbar = plt.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, fraction=0.035, pad=0.03)
    cbar.set_label("Expected target purity", fontsize=7.5)
    cbar.ax.tick_params(labelsize=7)


def plot_composition_axis(ax: plt.Axes, composition: pd.DataFrame, study_label: str) -> None:
    bottom = np.zeros(len(CAPTURE_GENES), dtype=float)
    for state in STATE_ORDER:
        subset = composition.loc[composition["cell_state"].eq(state)].set_index("capture_gene")
        values = np.array([subset["expected_post_sort_fraction"].get(gene, 0.0) for gene in CAPTURE_GENES], dtype=float)
        ax.bar(
            range(len(CAPTURE_GENES)),
            values,
            bottom=bottom,
            color=STATE_COLORS[state],
            width=0.72,
            edgecolor="white",
            linewidth=0.35,
            label=state,
        )
        bottom += values
    ax.set_title(study_label, fontsize=8.4, fontweight="bold", pad=4)
    ax.set_xticks(range(len(CAPTURE_GENES)), [rf"$\mathit{{{g}}}$" for g in CAPTURE_GENES], rotation=32, ha="right", fontsize=7.2)
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.5, 1.0], ["0", "50", "100"])
    ax.tick_params(axis="y", labelsize=6.5)
    ax.spines[["top", "right"]].set_visible(False)


def render_figure(
    data: pd.DataFrame,
    summary: pd.DataFrame,
    composition: pd.DataFrame,
    studies: pd.DataFrame,
    max_umap_cells: int,
    random_state: int,
    sensitivity: float,
    background: float,
) -> tuple[plt.Figure, pd.DataFrame]:
    figure = plt.figure(figsize=(17.2, 15.8), facecolor="white")
    outer = figure.add_gridspec(
        3,
        1,
        height_ratios=[1.42, 1.02, 1.05],
        left=0.055,
        right=0.985,
        bottom=0.055,
        top=0.875,
        hspace=0.46,
    )
    umap_grid = outer[0, 0].subgridspec(2, 4, wspace=0.30, hspace=0.50)
    matrix_ax = figure.add_subplot(outer[1, 0])
    composition_grid = outer[2, 0].subgridspec(2, 4, wspace=0.30, hspace=0.50)

    umap_data = downsample_umap(data, max_umap_cells, random_state)
    sample_key_rows: list[dict[str, object]] = []
    for idx, study in enumerate(studies.itertuples(index=False)):
        ax = figure.add_subplot(umap_grid[idx // 4, idx % 4])
        subset = umap_data.loc[umap_data["study_id"].astype(str).eq(str(study.study_id))]
        plot_sample_umap(ax, subset, sample_key_rows)
    if len(studies) < 8:
        for idx in range(len(studies), 8):
            ax = figure.add_subplot(umap_grid[idx // 4, idx % 4])
            ax.set_axis_off()

    plot_capture_matrix(matrix_ax, summary, studies)

    legend_handles = None
    legend_labels = None
    for idx, study in enumerate(studies.itertuples(index=False)):
        ax = figure.add_subplot(composition_grid[idx // 4, idx % 4])
        subset = composition.loc[composition["study_id"].astype(str).eq(str(study.study_id))]
        plot_composition_axis(ax, subset, str(study.study_label))
        if idx % 4 == 0:
            ax.set_ylabel("Expected post-sort\ncomposition (%)", fontsize=7.5)
        else:
            ax.tick_params(axis="y", labelleft=False)
        if legend_handles is None:
            legend_handles, legend_labels = ax.get_legend_handles_labels()
    if len(studies) < 8:
        legend_ax = figure.add_subplot(composition_grid[1, 3])
        legend_ax.set_axis_off()
        legend_ax.legend(legend_handles, legend_labels, loc="center left", frameon=False, fontsize=7.4)

    figure.text(0.055, 0.982, "Cross-study receptor-guided enrichment framework", fontsize=16, fontweight="bold", ha="left", va="top")
    figure.text(
        0.055,
        0.956,
        "RNA-expression proxy for comparing ERBB4, CXCR4 and PLXNA2 capture—not surface-protein or antibody validation",
        fontsize=9.2,
        color="#50565B",
        ha="left",
        va="top",
    )
    panel_a_y = outer[0, 0].get_position(figure).y1 + 0.018
    panel_b_y = outer[1, 0].get_position(figure).y1 + 0.018
    panel_c_y = outer[2, 0].get_position(figure).y1 + 0.018
    figure.text(0.035, panel_a_y, "a", fontsize=13, fontweight="bold")
    figure.text(0.055, panel_a_y, "Cross-study UMAPs by sample", fontsize=11.2, fontweight="bold")
    figure.text(0.213, panel_a_y, "S# labels map to full sample names in the sample-key table", fontsize=7.4, color="#555B61")
    figure.text(0.035, panel_b_y, "b", fontsize=13, fontweight="bold")
    figure.text(0.055, panel_b_y, "Theoretical capture performance", fontsize=11.2, fontweight="bold")
    figure.text(0.035, panel_c_y, "c", fontsize=13, fontweight="bold")
    figure.text(
        0.055,
        panel_c_y,
        f"Expected post-MACS cell-state composition ({sensitivity:.0%} marker-positive retention; {background:.0%} background carryover)",
        fontsize=11.2,
        fontweight="bold",
    )
    figure.text(
        0.50,
        0.018,
        "Target: GAD2-positive, LHX6- or LHX8-positive, non-cycling cells. Positivity threshold: 0.5 log1p(CP10K). Studies are empirical protocol-cohort proxies.",
        fontsize=7.3,
        color="#555B61",
        ha="center",
    )
    return figure, pd.DataFrame(sample_key_rows)


def save_figure(fig: plt.Figure, final_dir: Path, dpi: int) -> list[Path]:
    outputs = []
    for ext, ext_dpi in [("png", dpi), ("pdf", 300), ("svg", dpi)]:
        path = final_dir / "figures" / ext / f"{OUTPUT_STEM}.{ext}"
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=ext_dpi, bbox_inches="tight", facecolor="white")
        outputs.append(path)
    plt.close(fig)
    return outputs


def write_readme(final_dir: Path, args: argparse.Namespace) -> Path:
    path = final_dir / "README.md"
    path.write_text(
        f"""# Cross-study receptor-guided enrichment framework

This final-figure candidate compares the theoretical cell populations retained by
ERBB4, CXCR4, or PLXNA2 positive selection across seven organoid study/protocol
cohorts.

## Model contract

- Expression scale: existing log1p(CP10K) cross-study tables.
- Marker-positive threshold: `{args.positive_threshold:g}` log1p(CP10K).
- Independent target: GAD2-positive and LHX6- or LHX8-positive, excluding
  MKI67-positive cycling cells at the same threshold.
- Working in-silico MACS scenario: `{args.capture_sensitivity:.0%}` retention of
  marker-positive cells and `{args.background_retention:.0%}` nonspecific
  retention of marker-negative cells.
- Samarasinghe is controls-only; DIV90 clusters 6 and 7 are excluded as stressed
  cells, matching the existing cross-study visualization contract.

This is an RNA-expression proxy. It does not establish surface localization,
antibody specificity, epitope accessibility, magnetic-bead recovery, or live-cell
compatibility. Protein/flow validation is required before choosing a MACS reagent.

## Figure panels

- a: UMAPs colored by sample, with sample codes mapped in
  `tables/cross_study_sample_key.tsv`.
- b: expected target recovery, post-sort target purity, and fold enrichment.
- c: expected post-sort marker-state composition.

Sensitivity combinations are provided in
`tables/cross_study_receptor_macs_sensitivity_summary.tsv`.

The across-study receptor ranking excludes cohorts with an input independent-
target fraction below 1%, because fold enrichment is unstable when only a
handful of target cells are observed. The full cohort-specific estimates remain
available in the capture-summary table.
"""
    )
    return path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    args = parse_args()
    project_root = args.project_root.resolve()
    final_dir = (args.final_dir or project_root / "final_figures" / FINAL_FOLDER_DEFAULT).resolve()
    if args.fresh and final_dir.exists():
        shutil.rmtree(final_dir)
    for subdir in ["figures/png", "figures/pdf", "figures/svg", "tables", "code", "provenance"]:
        (final_dir / subdir).mkdir(parents=True, exist_ok=True)

    data, availability = load_data(project_root, args.base_run_label, args.receptor_run_label)
    filtered, filter_summary = _apply_internal_umap_plot_filters(data)
    filtered = add_marker_plot_coordinates(filtered)
    specs = default_cross_study_marker_specs(project_root)
    studies = ordered_studies(filtered, specs)
    filtered = add_sample_metadata(filtered, studies)
    filtered = assign_states(filtered, args.positive_threshold)

    summary, composition, sample_summary = summarize_capture(
        filtered,
        args.positive_threshold,
        args.capture_sensitivity,
        args.background_retention,
    )
    sensitivity = sensitivity_summary(filtered, args.positive_threshold)
    rank_summary = cross_study_rank_summary(summary)
    fig, sample_key = render_figure(
        filtered,
        summary,
        composition,
        studies,
        args.max_umap_cells_per_study,
        args.random_state,
        args.capture_sensitivity,
        args.background_retention,
    )
    outputs = save_figure(fig, final_dir, args.dpi)

    tables = {
        "cross_study_receptor_capture_summary.tsv": summary,
        "cross_study_receptor_post_macs_composition.tsv": composition,
        "cross_study_receptor_sample_capture_summary.tsv": sample_summary,
        "cross_study_receptor_macs_sensitivity_summary.tsv": sensitivity,
        "cross_study_receptor_rank_summary.tsv": rank_summary,
        "cross_study_sample_key.tsv": sample_key,
        "cross_study_receptor_expression_availability.tsv": availability,
        "cross_study_receptor_plot_filter_summary.tsv": filter_summary,
    }
    for name, table in tables.items():
        table.to_csv(final_dir / "tables" / name, sep="\t", index=False)

    readme = write_readme(final_dir, args)
    script_path = Path(__file__).resolve()
    shutil.copy2(script_path, final_dir / "code" / script_path.name)
    for source in [
        Path("scripts/12_export_cross_study_marker_expression_tables.R"),
        Path("python_notebooks/scripts/extract_cross_study_receptor_genes_from_h5ad.py"),
        Path("config/cross_study_receptor_macs_genes.tsv"),
        Path("slurm_templates/63_cross_study_receptor_macs_enrichment.sbatch.template"),
        Path("python_notebooks/HANDOFF_cross_study_receptor_macs_enrichment.md"),
    ]:
        if source.exists():
            destination = final_dir / ("code" if source.suffix in {".py", ".R", ".template"} else "provenance") / source.name
            shutil.copy2(source, destination)

    commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=False)
    manifest = pd.DataFrame(
        [
            ("rendered_at", datetime.now().astimezone().isoformat()),
            ("project_root", str(project_root)),
            ("base_run_label", args.base_run_label),
            ("receptor_run_label", args.receptor_run_label),
            ("positive_threshold_log1p_cp10k", f"{args.positive_threshold:g}"),
            ("capture_sensitivity", f"{args.capture_sensitivity:g}"),
            ("background_retention", f"{args.background_retention:g}"),
            ("target_definition", "GAD2+ AND (LHX6+ OR LHX8+) AND not MKI67+"),
            ("git_commit", commit.stdout.strip() if commit.returncode == 0 else "unavailable"),
        ],
        columns=["key", "value"],
    )
    manifest.to_csv(final_dir / "provenance" / "render_manifest.tsv", sep="\t", index=False)
    checksum_files = [readme, *outputs, *sorted((final_dir / "tables").glob("*.tsv")), *sorted((final_dir / "code").glob("*")), *sorted((final_dir / "provenance").glob("*"))]
    checksum_path = final_dir / "provenance" / "sha256_manifest.txt"
    checksum_path.write_text("".join(f"{sha256(path)}  {path.relative_to(final_dir)}\n" for path in checksum_files if path.is_file() and path != checksum_path))

    print(f"final_dir={final_dir}")
    print(f"figure={outputs[0]}")
    print(f"n_cells={len(filtered)}")
    print(f"n_studies={filtered['study_id'].nunique()}")
    print(summary[["study_id", "capture_gene", "pct_marker_positive", "expected_target_recovery", "expected_post_sort_target_purity", "expected_target_enrichment_fold"]].to_string(index=False))


if __name__ == "__main__":
    main()
