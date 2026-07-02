#!/usr/bin/env python3
"""Render integrated publication-style cross-study MGE marker figures."""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
import re
import shutil
import subprocess

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams.update(
    {
        "font.family": "Arial",
        "font.sans-serif": ["Arial", "Nimbus Sans", "Liberation Sans", "DejaVu Sans"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mge_organoid_python.cross_study_marker_expression import (
    PDF_EXPORT_DPI,
    PNG_EXPORT_DPI,
    SVG_EXPORT_DPI,
    _apply_internal_umap_plot_filters,
    _finite_quantile,
    add_marker_plot_coordinates,
    default_cross_study_marker_specs,
    load_marker_expression_tables,
    table_dir,
)
from render_cross_study_lhx6_expression_final_panel import (
    STUDY_COLORS,
    add_sample_metadata,
    normalize_gene,
    ordered_studies,
)


PROJECT_ROOT_DEFAULT = Path("/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder")
RUN_LABEL_DEFAULT = "cross_study_marker_expression_v12_pv_precursors_final_candidate_plus_vipr2"
FINAL_FOLDER_DEFAULT = "fig_cross_study_integrated_mge_marker_expression_v7_overlay_zero_baseline_candidate"
THRESHOLD_DEFAULT = 0.5
FIGURE1_GENES_DEFAULT = ["NKX2.1", "LHX6", "LHX8", "CRABP1", "TAC1", "VIPR2", "SLC6A1", "ZEB2", "SST", "PVALB"]
FIGURE2_GENES_DEFAULT = ["NKX2.1", "LHX6", "LHX8"]
FIGURE1_STEM = "cross_study_canonical_mge_marker_expression_violin_positive_matrix"
FIGURE2_STEM = "this_study_div30_div90_mge_marker_replicate_reproducibility"
FIGURE2_VARIANT_STEM_TEMPLATE = "this_study_div30_div90_mge_marker_replicate_reproducibility_cell_cutoff_{cutoff_token}"

GENE_LABELS = {
    "NKX2-1": "NKX2.1",
    "LHX6": "LHX6",
    "LHX8": "LHX8",
    "CRABP1": "CRABP1",
    "TAC1": "TAC1",
    "VIPR2": "VIPR2",
    "SLC6A1": "SLC6A1",
    "ZEB2": "ZEB2",
    "SST": "SST",
    "PVALB": "PVALB",
}
GENE_COLORS = {
    "NKX2-1": "#0072B2",
    "LHX6": "#009E73",
    "LHX8": "#CC79A7",
    "CRABP1": "#D55E00",
}
GENE_TEXT_COLOR = "#202020"
DIV_COLORS = {
    "varela_div30": "#8FBCE6",
    "varela_div90": "#2F855A",
}
FIGURE_STUDY_ORDER = [
    "varela_div30",
    "varela_div90",
    "siebert_2026",
    "walsh",
    "samarasinghe_2021",
    "bershteyn_2023",
    "bershteyn_2025",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT_DEFAULT)
    parser.add_argument("--run-label", default=RUN_LABEL_DEFAULT)
    parser.add_argument("--final-dir", type=Path, default=None)
    parser.add_argument("--positive-threshold", type=float, default=THRESHOLD_DEFAULT)
    parser.add_argument("--figure1-genes", nargs="+", default=FIGURE1_GENES_DEFAULT)
    parser.add_argument(
        "--figure1-layout",
        choices=["overlay", "paired_rows"],
        default="overlay",
        help="Figure 1 layout. overlay combines percent-positive bars and violins in one row per gene.",
    )
    parser.add_argument("--figure2-genes", nargs="+", default=FIGURE2_GENES_DEFAULT)
    parser.add_argument("--replicate-cell-cutoffs", nargs="+", default=["none", "0.5", "1"])
    parser.add_argument("--fresh", action="store_true", help="Delete and recreate the final figure folder first.")
    parser.add_argument("--random-state", type=int, default=1729)
    parser.add_argument(
        "--max-violin-cells-per-sample",
        type=int,
        default=1800,
        help="Maximum cells used to draw each sample violin. Summary tables and percentages always use all cells.",
    )
    return parser.parse_args()


def gene_label(gene: str) -> str:
    return GENE_LABELS.get(gene, gene)


def finite_values(data: pd.DataFrame, gene: str) -> np.ndarray:
    values = pd.to_numeric(data[gene], errors="coerce").to_numpy(dtype=float)
    return values[np.isfinite(values)]


def save_figure(fig: plt.Figure, final_dir: Path, stem: str) -> list[Path]:
    outputs = []
    for ext, dpi in (("png", PNG_EXPORT_DPI), ("pdf", PDF_EXPORT_DPI), ("svg", SVG_EXPORT_DPI)):
        path = final_dir / "figures" / ext / f"{stem}.{ext}"
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        outputs.append(path)
    plt.close(fig)
    return outputs


def compute_sample_positions(sample_summary: pd.DataFrame) -> pd.Series:
    sample_x: list[float] = []
    current_x = 0.0
    previous_study: str | None = None
    for study_id in sample_summary["study_id"].astype(str).tolist():
        if previous_study is not None and study_id != previous_study:
            current_x += 1.4
        sample_x.append(current_x)
        current_x += 1.0
        previous_study = study_id
    return pd.Series(sample_x, index=sample_summary["sample_plot_id"], name="sample_x")


def summarize_by_sample(data: pd.DataFrame, genes: list[str], threshold: float) -> pd.DataFrame:
    group_cols = [
        "study_id",
        "study_label",
        "study_plot_order",
        "sample",
        "sample_display_label",
        "sample_age_order",
        "sample_order_within_study",
        "sample_plot_id",
    ]
    rows = []
    for keys, group in data.groupby(group_cols, sort=False, observed=True):
        base = dict(zip(group_cols, keys, strict=True))
        for gene in genes:
            values = finite_values(group, gene)
            positive = values >= threshold
            positive_values = values[positive]
            row = dict(base)
            row.update(
                {
                    "gene": gene,
                    "gene_label": gene_label(gene),
                    "positive_threshold_log1p_cp10k": threshold,
                    "n_cells": int(group.shape[0]),
                    "n_finite_expression_values": int(values.size),
                    "n_positive_cells": int(positive.sum()),
                    "fraction_positive_cells": float(positive.mean()) if values.size else np.nan,
                    "pct_positive_cells": float(positive.mean() * 100.0) if values.size else np.nan,
                    "mean_log1p_cp10k": float(np.mean(values)) if values.size else np.nan,
                    "mean_positive_log1p_cp10k": float(np.mean(positive_values)) if positive_values.size else np.nan,
                    "expr_q25": _finite_quantile(values, 0.25),
                    "median_log1p_cp10k": _finite_quantile(values, 0.50),
                    "expr_q75": _finite_quantile(values, 0.75),
                    "expr_q90": _finite_quantile(values, 0.90),
                    "expr_q95": _finite_quantile(values, 0.95),
                    "positive_expr_q25": _finite_quantile(positive_values, 0.25),
                    "positive_median_log1p_cp10k": _finite_quantile(positive_values, 0.50),
                    "positive_expr_q75": _finite_quantile(positive_values, 0.75),
                    "denominator": f"Visualization-filtered cells with finite {gene} expression in each sample",
                }
            )
            rows.append(row)
    out = pd.DataFrame(rows)
    out["gene_order"] = out["gene"].map({gene: idx for idx, gene in enumerate(genes)})
    return out.sort_values(
        [
            "gene_order",
            "study_plot_order",
            "sample_age_order",
            "sample_order_within_study",
            "sample_display_label",
        ],
        kind="stable",
    ).reset_index(drop=True)


def summarize_by_study(data: pd.DataFrame, genes: list[str], threshold: float) -> pd.DataFrame:
    rows = []
    for (study_id, study_label, study_plot_order), group in data.groupby(
        ["study_id", "study_label", "study_plot_order"], sort=False, observed=True
    ):
        for gene in genes:
            values = finite_values(group, gene)
            positive = values >= threshold
            positive_values = values[positive]
            rows.append(
                {
                    "study_id": study_id,
                    "study_label": study_label,
                    "study_plot_order": study_plot_order,
                    "gene": gene,
                    "gene_label": gene_label(gene),
                    "positive_threshold_log1p_cp10k": threshold,
                    "n_cells": int(group.shape[0]),
                    "n_finite_expression_values": int(values.size),
                    "n_positive_cells": int(positive.sum()),
                    "fraction_positive_cells": float(positive.mean()) if values.size else np.nan,
                    "pct_positive_cells": float(positive.mean() * 100.0) if values.size else np.nan,
                    "mean_log1p_cp10k": float(np.mean(values)) if values.size else np.nan,
                    "mean_positive_log1p_cp10k": float(np.mean(positive_values)) if positive_values.size else np.nan,
                    "median_log1p_cp10k": _finite_quantile(values, 0.50),
                    "positive_median_log1p_cp10k": _finite_quantile(positive_values, 0.50),
                    "expr_q90": _finite_quantile(values, 0.90),
                    "expr_q95": _finite_quantile(values, 0.95),
                }
            )
    out = pd.DataFrame(rows)
    out["gene_order"] = out["gene"].map({gene: idx for idx, gene in enumerate(genes)})
    return out.sort_values(["gene_order", "study_plot_order"], kind="stable").reset_index(drop=True)


def apply_figure_study_order(studies: pd.DataFrame) -> pd.DataFrame:
    """Apply the requested visual study order without changing source data logic."""
    out = studies.copy()
    order = {study_id: idx for idx, study_id in enumerate(FIGURE_STUDY_ORDER)}
    fallback = out["study_plot_order"].max() + 1 if "study_plot_order" in out else len(order)
    out["study_plot_order"] = [
        order.get(str(study_id), fallback + idx)
        for idx, study_id in enumerate(out["study_id"].astype(str).tolist())
    ]
    return out.sort_values(["study_plot_order", "study_label"], kind="stable").reset_index(drop=True)


def cutoff_token(cutoff: float | None) -> str:
    if cutoff is None:
        return "none"
    return f"{cutoff:g}".replace(".", "_")


def cutoff_label(cutoff: float | None) -> str:
    if cutoff is None:
        return "No cutoff"
    return f"Cells >= {cutoff:g}"


def parse_cell_cutoffs(values: list[str]) -> list[float | None]:
    cutoffs: list[float | None] = []
    for value in values:
        text = str(value).strip().lower()
        if text in {"none", "no", "all", "all_cells", "uncut", "no_cutoff"}:
            cutoffs.append(None)
        else:
            cutoffs.append(float(text))
    deduped: list[float | None] = []
    seen = set()
    for cutoff in cutoffs:
        key = "none" if cutoff is None else f"{cutoff:g}"
        if key in seen:
            continue
        seen.add(key)
        deduped.append(cutoff)
    return deduped


def summarize_replicate_expression(data: pd.DataFrame, genes: list[str], cell_cutoffs: list[float | None]) -> pd.DataFrame:
    group_cols = [
        "study_id",
        "study_label",
        "study_plot_order",
        "sample",
        "sample_display_label",
        "sample_age_order",
        "sample_order_within_study",
        "sample_plot_id",
    ]
    rows = []
    subset = data.loc[data["study_id"].isin(["varela_div30", "varela_div90"])].copy()
    for keys, group in subset.groupby(group_cols, sort=False, observed=True):
        base = dict(zip(group_cols, keys, strict=True))
        for gene in genes:
            values = finite_values(group, gene)
            for cutoff in cell_cutoffs:
                if cutoff is None:
                    included = np.isfinite(values)
                else:
                    included = values >= cutoff
                included_values = values[included]
                row = dict(base)
                row.update(
                    {
                        "gene": gene,
                        "gene_label": gene_label(gene),
                        "replicate_cell_cutoff_log1p_cp10k": cutoff if cutoff is not None else np.nan,
                        "replicate_cell_cutoff_label": cutoff_label(cutoff),
                        "replicate_cell_cutoff_token": cutoff_token(cutoff),
                        "n_cells": int(group.shape[0]),
                        "n_finite_expression_values": int(values.size),
                        "n_included_cells": int(included_values.size),
                        "fraction_included_cells": float(included_values.size / values.size) if values.size else np.nan,
                        "pct_included_cells": float(included_values.size / values.size * 100.0) if values.size else np.nan,
                        "replicate_plot_value_log1p_cp10k": float(np.mean(included_values)) if included_values.size else np.nan,
                        "included_median_log1p_cp10k": _finite_quantile(included_values, 0.50),
                        "included_expr_q25": _finite_quantile(included_values, 0.25),
                        "included_expr_q75": _finite_quantile(included_values, 0.75),
                        "all_cell_mean_log1p_cp10k": float(np.mean(values)) if values.size else np.nan,
                        "all_cell_median_log1p_cp10k": _finite_quantile(values, 0.50),
                    }
                )
                rows.append(row)
    out = pd.DataFrame(rows)
    out["gene_order"] = out["gene"].map({gene: idx for idx, gene in enumerate(genes)})
    cutoff_order = {cutoff_token(cutoff): idx for idx, cutoff in enumerate(cell_cutoffs)}
    out["cutoff_order"] = out["replicate_cell_cutoff_token"].map(cutoff_order)
    return out.sort_values(
        [
            "cutoff_order",
            "gene_order",
            "study_plot_order",
            "sample_order_within_study",
            "sample_display_label",
        ],
        kind="stable",
    ).reset_index(drop=True)


def ordered_sample_summary(sample_summary: pd.DataFrame, first_gene: str) -> pd.DataFrame:
    samples = sample_summary.loc[sample_summary["gene"].eq(first_gene)].copy()
    return samples.sort_values(
        ["study_plot_order", "sample_age_order", "sample_order_within_study", "sample_display_label"],
        kind="stable",
    ).reset_index(drop=True)


def add_sample_x(data: pd.DataFrame, samples: pd.DataFrame) -> pd.DataFrame:
    positions = compute_sample_positions(samples)
    out = data.copy()
    out["sample_x"] = out["sample_plot_id"].map(positions)
    return out


def apply_axis_style(ax: plt.Axes) -> None:
    ax.grid(axis="y", color="#d8d8d8", linewidth=0.55)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_linewidth(0.7)
    ax.tick_params(axis="both", labelsize=6.3, width=0.6, length=2.5)


def remove_overlay_tick_lines(ax: plt.Axes) -> None:
    ax.grid(False)
    ax.yaxis.grid(False, which="both")
    ax.xaxis.grid(False, which="both")
    for gridline in [*ax.get_ygridlines(), *ax.get_xgridlines()]:
        gridline.set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.tick_params(axis="x", length=0)


def draw_study_separators(ax: plt.Axes, samples: pd.DataFrame, sample_x: np.ndarray, *, ymin: float = 0.0, ymax: float = 1.0) -> None:
    study_values = samples["study_id"].astype(str).to_numpy()
    boundaries = np.flatnonzero(study_values[1:] != study_values[:-1]) + 1
    for boundary in boundaries:
        xpos = float(np.mean([sample_x[boundary - 1], sample_x[boundary]]))
        ax.axvline(xpos, ymin=ymin, ymax=ymax, color="#5b5b5b", linewidth=0.7, clip_on=False)


def draw_study_headers(ax: plt.Axes, samples: pd.DataFrame, sample_x: np.ndarray, *, y: float = 1.10) -> None:
    for _, group in samples.groupby("study_id", sort=False):
        idx = group.index.to_numpy()
        center = float(np.mean(sample_x[idx]))
        label = str(group["study_label"].iloc[0])
        label = label.replace("This Study, ", "This Study\n")
        ax.text(center, y, label, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=7.2)


def plot_figure1_paired_rows(
    annotated: pd.DataFrame,
    sample_summary: pd.DataFrame,
    figure_genes: list[str],
    threshold: float,
    final_dir: Path,
    random_state: int,
    max_violin_cells_per_sample: int,
) -> list[Path]:
    rng = np.random.default_rng(random_state)
    samples = ordered_sample_summary(sample_summary, figure_genes[0])
    samples = add_sample_x(samples, samples)
    plot_data = add_sample_x(annotated, samples)
    sample_x = samples["sample_x"].to_numpy(dtype=float)
    x_min = float(sample_x.min() - 0.8)
    x_max = float(sample_x.max() + 0.8)

    n_genes = len(figure_genes)
    fig_width = max(17.5, 0.31 * samples.shape[0] + 4.4)
    fig_height = 13.4 if n_genes <= 4 else 3.15 * n_genes + 1.2
    fig, axes = plt.subplots(
        n_genes * 2,
        1,
        figsize=(fig_width, fig_height),
        sharex=True,
        gridspec_kw={"height_ratios": [1.08, 0.74] * n_genes, "hspace": 0.16},
        constrained_layout=False,
    )
    violin_axes = axes[0::2]
    bar_axes = axes[1::2]

    expr_max = 0.0
    for gene in figure_genes:
        values = finite_values(plot_data, gene)
        if values.size:
            expr_max = max(expr_max, float(np.quantile(values, 0.995)))
    expr_ylim = max(threshold * 1.65, expr_max * 1.04, 1.0)

    for gene, ax in zip(figure_genes, violin_axes, strict=True):
        values_by_sample = []
        positions = []
        for sample in samples.itertuples(index=False):
            subset = plot_data.loc[plot_data["sample_plot_id"].eq(sample.sample_plot_id)]
            values = finite_values(subset, gene)
            if max_violin_cells_per_sample > 0 and values.size > max_violin_cells_per_sample:
                chosen = rng.choice(values.size, size=max_violin_cells_per_sample, replace=False)
                values = values[np.sort(chosen)]
            values_by_sample.append(values)
            positions.append(float(sample.sample_x))
        violins = ax.violinplot(
            values_by_sample,
            positions=positions,
            widths=0.72,
            showmeans=False,
            showmedians=False,
            showextrema=False,
        )
        for body, study_id in zip(violins["bodies"], samples["study_id"].astype(str).tolist(), strict=True):
            body.set_facecolor(STUDY_COLORS.get(study_id, "#777777"))
            body.set_edgecolor("#242424")
            body.set_linewidth(0.24)
            body.set_alpha(0.76)
        medians_x = []
        medians_y = []
        full_summary = sample_summary.loc[sample_summary["gene"].eq(gene)].set_index("sample_plot_id")
        for xpos, sample in zip(positions, samples.itertuples(index=False), strict=True):
            row = full_summary.loc[sample.sample_plot_id]
            if not np.isfinite(row["median_log1p_cp10k"]):
                continue
            q25 = float(row["expr_q25"])
            median = float(row["median_log1p_cp10k"])
            q75 = float(row["expr_q75"])
            ax.vlines(xpos, q25, q75, color="#202020", linewidth=0.62)
            medians_x.append(xpos)
            medians_y.append(median)
        if medians_x:
            ax.scatter(medians_x, medians_y, s=4, color="#202020", zorder=3)
        ax.axhline(threshold, color="#202020", linewidth=0.75, linestyle=(0, (3, 2)))
        ax.set_ylim(0, expr_ylim)
        ax.set_xlim(x_min, x_max)
        ax.set_ylabel("log1p(CP10K)", fontsize=7.2)
        ax.text(
            -0.045,
            0.50,
            gene_label(gene),
            transform=ax.transAxes,
            ha="right",
            va="center",
            fontsize=10,
            fontweight="bold",
            color=GENE_TEXT_COLOR,
        )
        draw_study_separators(ax, samples, sample_x)
        apply_axis_style(ax)

    draw_study_headers(violin_axes[0], samples, sample_x, y=1.16)
    violin_axes[0].set_title("Expression distributions by sample", fontsize=11, pad=24)

    pct_max = max(108.0, float(np.nanmax(sample_summary["pct_positive_cells"].to_numpy(dtype=float))) + 6.0)
    for gene, ax in zip(figure_genes, bar_axes, strict=True):
        gene_summary = sample_summary.loc[sample_summary["gene"].eq(gene)].copy()
        gene_summary = gene_summary.set_index("sample_plot_id").loc[samples["sample_plot_id"]].reset_index()
        pct_values = gene_summary["pct_positive_cells"].to_numpy(dtype=float)
        colors = [STUDY_COLORS.get(study_id, "#777777") for study_id in gene_summary["study_id"].astype(str)]
        ax.bar(
            sample_x,
            pct_values,
            width=0.80,
            color=colors,
            edgecolor="#202020",
            linewidth=0.22,
        )
        for xpos, pct in zip(sample_x, pct_values, strict=True):
            if not np.isfinite(pct):
                continue
            label_y = min(float(pct) + 1.3, pct_max - 1.6)
            ax.text(
                xpos,
                label_y,
                f"{pct:.0f}",
                ha="center",
                va="bottom",
                fontsize=3.9,
                color="#202020",
                clip_on=False,
            )
        ax.set_ylim(0, pct_max)
        ax.set_ylabel("% >= 0.5", fontsize=7.2)
        ax.text(
            -0.045,
            0.50,
            gene_label(gene),
            transform=ax.transAxes,
            ha="right",
            va="center",
            fontsize=10,
            fontweight="bold",
            color=GENE_TEXT_COLOR,
        )
        draw_study_separators(ax, samples, sample_x)
        apply_axis_style(ax)

    for ax in axes[:-1]:
        ax.set_xticks(sample_x)
        ax.tick_params(axis="x", labelbottom=False)
    axes[-1].set_xticks(sample_x)
    axes[-1].set_xticklabels(samples["sample_display_label"].tolist(), rotation=90, ha="center", va="top", fontsize=5.3)

    fig.suptitle("Canonical MGE marker expression across studies", fontsize=13, y=0.992)
    fig.subplots_adjust(left=0.078, right=0.990, top=0.935, bottom=0.215)
    return save_figure(fig, final_dir, FIGURE1_STEM)


def plot_figure1_overlay(
    annotated: pd.DataFrame,
    sample_summary: pd.DataFrame,
    figure_genes: list[str],
    threshold: float,
    final_dir: Path,
    random_state: int,
    max_violin_cells_per_sample: int,
) -> list[Path]:
    rng = np.random.default_rng(random_state)
    samples = ordered_sample_summary(sample_summary, figure_genes[0])
    samples = add_sample_x(samples, samples)
    plot_data = add_sample_x(annotated, samples)
    sample_x = samples["sample_x"].to_numpy(dtype=float)
    x_min = float(sample_x.min() - 0.8)
    x_max = float(sample_x.max() + 0.8)

    n_genes = len(figure_genes)
    fig_width = max(17.5, 0.31 * samples.shape[0] + 4.4)
    fig_height = max(10.5, 2.15 * n_genes + 2.6)
    fig, axes = plt.subplots(
        n_genes,
        1,
        figsize=(fig_width, fig_height),
        sharex=True,
        gridspec_kw={"hspace": 0.42},
        constrained_layout=False,
    )
    if n_genes == 1:
        axes = [axes]

    expr_max = 0.0
    for gene in figure_genes:
        values = finite_values(plot_data, gene)
        if values.size:
            expr_max = max(expr_max, float(np.quantile(values, 0.995)))
    expr_ylim = max(5.0, threshold * 1.65, expr_max * 1.04)

    right_axes = []
    for gene, ax in zip(figure_genes, axes, strict=True):
        gene_summary = sample_summary.loc[sample_summary["gene"].eq(gene)].copy()
        gene_summary = gene_summary.set_index("sample_plot_id").loc[samples["sample_plot_id"]].reset_index()
        pct_values = gene_summary["pct_positive_cells"].to_numpy(dtype=float)
        colors = [STUDY_COLORS.get(study_id, "#777777") for study_id in gene_summary["study_id"].astype(str)]

        pct_ax = ax.twinx()
        right_axes.append(pct_ax)
        ax.set_zorder(2)
        pct_ax.set_zorder(1)
        ax.patch.set_visible(False)
        pct_ax.bar(
            sample_x,
            pct_values,
            width=0.80,
            color=colors,
            edgecolor="none",
            linewidth=0,
            alpha=0.30,
            zorder=0,
        )
        pct_ax.set_ylim(0, 100)
        pct_ax.set_ylabel("% positive", fontsize=7.2)
        remove_overlay_tick_lines(pct_ax)
        pct_ax.spines["top"].set_visible(False)
        pct_ax.spines["left"].set_visible(False)
        pct_ax.spines["right"].set_linewidth(0.7)
        pct_ax.tick_params(axis="y", labelsize=6.3, width=0.6, length=2.5)
        pct_ax.tick_params(axis="x", bottom=False, labelbottom=False)

        values_by_sample = []
        positions = []
        for sample in samples.itertuples(index=False):
            subset = plot_data.loc[plot_data["sample_plot_id"].eq(sample.sample_plot_id)]
            values = finite_values(subset, gene)
            if max_violin_cells_per_sample > 0 and values.size > max_violin_cells_per_sample:
                chosen = rng.choice(values.size, size=max_violin_cells_per_sample, replace=False)
                values = values[np.sort(chosen)]
            values_by_sample.append(values)
            positions.append(float(sample.sample_x))
        violins = ax.violinplot(
            values_by_sample,
            positions=positions,
            widths=0.72,
            showmeans=False,
            showmedians=False,
            showextrema=False,
        )
        for body, study_id in zip(violins["bodies"], samples["study_id"].astype(str).tolist(), strict=True):
            body.set_facecolor(STUDY_COLORS.get(study_id, "#777777"))
            body.set_edgecolor("none")
            body.set_linewidth(0)
            body.set_alpha(1.00)
            body.set_zorder(3)
        medians_x = []
        medians_y = []
        full_summary = sample_summary.loc[sample_summary["gene"].eq(gene)].set_index("sample_plot_id")
        for xpos, sample in zip(positions, samples.itertuples(index=False), strict=True):
            row = full_summary.loc[sample.sample_plot_id]
            pct = float(row["pct_positive_cells"])
            if np.isfinite(pct):
                ax.text(
                    xpos,
                    1.025,
                    f"{pct:.0f}",
                    transform=ax.get_xaxis_transform(),
                    ha="center",
                    va="bottom",
                    fontsize=4.0,
                    color="#202020",
                    clip_on=False,
                )
            if not np.isfinite(row["median_log1p_cp10k"]):
                continue
            q25 = float(row["expr_q25"])
            median = float(row["median_log1p_cp10k"])
            q75 = float(row["expr_q75"])
            ax.vlines(xpos, q25, q75, color="#202020", linewidth=0.62, zorder=4)
            medians_x.append(xpos)
            medians_y.append(median)
        if medians_x:
            ax.scatter(medians_x, medians_y, s=4, color="#202020", zorder=5)
        ax.axhline(0, color="#202020", linewidth=0.65, zorder=2, clip_on=False)
        ax.axhline(threshold, color="#202020", linewidth=0.75, linestyle=(0, (3, 2)), zorder=2)
        ax.set_ylim(0, expr_ylim)
        ax.set_xlim(x_min, x_max)
        ax.set_ylabel("log1p(CP10K) expression", fontsize=7.2)
        ax.text(
            -0.045,
            0.50,
            gene_label(gene),
            transform=ax.transAxes,
            ha="right",
            va="center",
            fontsize=10,
            fontweight="bold",
            color=GENE_TEXT_COLOR,
        )
        draw_study_separators(ax, samples, sample_x)
        apply_axis_style(ax)
        remove_overlay_tick_lines(ax)

    draw_study_headers(axes[0], samples, sample_x, y=1.16)
    axes[0].set_title("Expression distributions and percent-positive cells by sample", fontsize=11, pad=24)
    for ax in axes[:-1]:
        ax.set_xticks(sample_x)
        ax.tick_params(axis="x", labelbottom=False)
    axes[-1].set_xticks(sample_x)
    axes[-1].set_xticklabels(samples["sample_display_label"].tolist(), rotation=90, ha="center", va="top", fontsize=5.3)
    for pct_ax in right_axes:
        pct_ax.set_xlim(x_min, x_max)

    fig.suptitle("Canonical MGE marker expression across studies", fontsize=13, y=0.992)
    fig.subplots_adjust(left=0.078, right=0.945, top=0.925, bottom=0.235)
    return save_figure(fig, final_dir, FIGURE1_STEM)


def plot_figure1(
    annotated: pd.DataFrame,
    sample_summary: pd.DataFrame,
    figure_genes: list[str],
    threshold: float,
    final_dir: Path,
    random_state: int,
    max_violin_cells_per_sample: int,
    figure1_layout: str,
) -> list[Path]:
    if figure1_layout == "paired_rows":
        return plot_figure1_paired_rows(
            annotated,
            sample_summary,
            figure_genes,
            threshold,
            final_dir,
            random_state,
            max_violin_cells_per_sample,
        )
    return plot_figure1_overlay(
        annotated,
        sample_summary,
        figure_genes,
        threshold,
        final_dir,
        random_state,
        max_violin_cells_per_sample,
    )


def plot_figure2_variant(
    replicate_summary: pd.DataFrame,
    figure_genes: list[str],
    cutoff: float | None,
    final_dir: Path,
    random_state: int,
) -> list[Path]:
    rng = np.random.default_rng(random_state)
    token = cutoff_token(cutoff)
    subset = replicate_summary.loc[replicate_summary["replicate_cell_cutoff_token"].eq(token)].copy()
    n_genes = len(figure_genes)
    fig, axes = plt.subplots(1, n_genes, figsize=(3.15 * n_genes, 4.25), sharey=True, constrained_layout=False)
    if n_genes == 1:
        axes = [axes]

    y_values = subset["replicate_plot_value_log1p_cp10k"].to_numpy(dtype=float)
    finite_y = y_values[np.isfinite(y_values)]
    line_floor = 0.0 if cutoff is None else cutoff * 1.6
    y_max = max(line_floor, float(np.nanmax(finite_y)) * 1.18 if finite_y.size else 1.0, 1.0)
    group_order = ["varela_div30", "varela_div90"]
    labels = ["DIV30", "DIV90"]

    for ax, gene in zip(axes, figure_genes, strict=True):
        values = []
        for study_id in group_order:
            group = subset.loc[subset["gene"].eq(gene) & subset["study_id"].eq(study_id)]
            group = group.sort_values(["sample_order_within_study", "sample_display_label"], kind="stable")
            vals = group["replicate_plot_value_log1p_cp10k"].to_numpy(dtype=float)
            values.append(vals[np.isfinite(vals)])
        bp = ax.boxplot(
            values,
            positions=[1, 2],
            widths=0.52,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "#202020", "linewidth": 1.0},
            boxprops={"facecolor": "none", "edgecolor": "#202020", "linewidth": 0.95},
            whiskerprops={"color": "#202020", "linewidth": 0.8},
            capprops={"color": "#202020", "linewidth": 0.8},
        )
        for patch, study_id in zip(bp["boxes"], group_order, strict=True):
            patch.set_edgecolor(DIV_COLORS[study_id])
            patch.set_facecolor("none")
        for xpos, vals, study_id in zip([1, 2], values, group_order, strict=True):
            jitter = rng.uniform(-0.075, 0.075, size=vals.size)
            ax.scatter(
                np.full(vals.size, xpos, dtype=float) + jitter,
                vals,
                s=29,
                facecolor=DIV_COLORS[study_id],
                edgecolor="#202020",
                linewidth=0.35,
                alpha=0.92,
                zorder=4,
            )
        if cutoff is not None:
            ax.axhline(cutoff, color="#202020", linewidth=0.75, linestyle=(0, (3, 2)))
        ax.set_title(gene_label(gene), fontsize=11, color=GENE_TEXT_COLOR, pad=8, fontweight="bold")
        ax.set_xticks([1, 2])
        ax.set_xticklabels(labels, fontsize=8)
        ax.set_xlim(0.45, 2.55)
        ax.set_ylim(0, y_max)
        apply_axis_style(ax)
    if cutoff is None:
        ylabel = "Sample mean expression, log1p(CP10K)"
        subtitle = "No per-cell cutoff"
    else:
        ylabel = f"Mean expression in cells >= {cutoff:g}, log1p(CP10K)"
        subtitle = f"Per-cell inclusion cutoff: expression >= {cutoff:g}"
    axes[0].set_ylabel(ylabel, fontsize=8.5)
    fig.suptitle(f"This Study DIV30 vs DIV90 replicate reproducibility\n{subtitle}", fontsize=11.5, y=0.985, linespacing=1.15)
    fig.subplots_adjust(left=0.095, right=0.985, top=0.805, bottom=0.145, wspace=0.16)
    stem = FIGURE2_VARIANT_STEM_TEMPLATE.format(cutoff_token=token)
    return save_figure(fig, final_dir, stem)


def plot_figure2_variants(
    replicate_summary: pd.DataFrame,
    figure_genes: list[str],
    cell_cutoffs: list[float | None],
    final_dir: Path,
    random_state: int,
) -> list[Path]:
    outputs = []
    for cutoff in cell_cutoffs:
        outputs.extend(plot_figure2_variant(replicate_summary, figure_genes, cutoff, final_dir, random_state))
    return outputs


def write_readme(
    final_dir: Path,
    project_root: Path,
    run_label: str,
    figure1_genes: list[str],
    figure2_genes: list[str],
    threshold: float,
    replicate_cell_cutoffs: list[float | None],
    figure1_layout: str,
) -> None:
    replicate_lines = []
    for cutoff in replicate_cell_cutoffs:
        token = cutoff_token(cutoff)
        stem = FIGURE2_VARIANT_STEM_TEMPLATE.format(cutoff_token=token)
        replicate_lines.append(f"- `figures/png/{stem}.png`: This Study DIV30 vs DIV90, replicate mean using `{cutoff_label(cutoff)}` cells.")
    text = f"""# Integrated cross-study MGE marker expression figures

Publication-style integrated marker-expression package derived from the existing cross-study LHX6 template.

Source prepared marker-expression run:
`{project_root}/results/cross_study_marker_expression/{run_label}`

Figures:
- `figures/png/{FIGURE1_STEM}.png`: canonical MGE marker comparison across studies and samples.
{chr(10).join(replicate_lines)}

Figure 1 layout:
- `{figure1_layout}`. In overlay mode, each gene has one combined panel with semi-transparent percent-positive bars on the right y-axis behind opaque expression violins on the left y-axis.

Expression scale:
- Cell-level values are the existing `log1p(CP10K)` values from the source marker-expression tables.
- No log10 conversion or rescaling is applied.
- Dashed reference line and positive-cell threshold: expression >= `{threshold:g}` `log1p(CP10K)`.

Figure 1 genes:
{chr(10).join(f"- `{gene_label(gene)}`" for gene in figure1_genes)}

Figure 2 genes:
{chr(10).join(f"- `{gene_label(gene)}`" for gene in figure2_genes)}

Denominator:
Visualization-filtered cells from the marker-expression workflow. Samarasinghe is controls-only, and DIV90 current clusters 6 and 7 are removed as stressed cells for plotting and summaries.

Figure 2 statistic:
Each dot is one This Study replicate/sample. For each replicate panel, the per-cell cutoff is applied first within each sample/gene, then the mean `log1p(CP10K)` expression is calculated from the included cells. Cutoff variants are `{", ".join(cutoff_label(value) for value in replicate_cell_cutoffs)}`.
"""
    final_dir.mkdir(parents=True, exist_ok=True)
    (final_dir / "README.md").write_text(text)


def write_provenance(
    final_dir: Path,
    project_root: Path,
    run_label: str,
    threshold: float,
    outputs: list[Path],
    figure1_genes: list[str],
    figure2_genes: list[str],
    replicate_cell_cutoffs: list[float | None],
    figure1_layout: str,
) -> None:
    provenance = final_dir / "provenance"
    code_dir = final_dir / "code"
    provenance.mkdir(parents=True, exist_ok=True)
    code_dir.mkdir(parents=True, exist_ok=True)
    script_path = Path(__file__).resolve()
    shutil.copy2(script_path, code_dir / script_path.name)
    for path in [
        Path("python_notebooks/scripts/render_cross_study_lhx6_expression_final_panel.py"),
        Path("python_notebooks/HANDOFF_cross_study_marker_synthesis_concept.md"),
        Path("python_notebooks/HANDOFF_transition_to_final_figs.md"),
    ]:
        if path.exists():
            shutil.copy2(path, provenance / path.name)

    git_commit = subprocess.run(["git", "rev-parse", "HEAD"], check=False, capture_output=True, text=True)
    git_status = subprocess.run(["git", "status", "--short"], check=False, capture_output=True, text=True)
    (provenance / "git_commit.txt").write_text(git_commit.stdout if git_commit.returncode == 0 else git_commit.stderr)
    (provenance / "git_status_short.txt").write_text(git_status.stdout if git_status.returncode == 0 else git_status.stderr)
    rows = [
        ("rendered_at", datetime.now().astimezone().isoformat()),
        ("project_root", str(project_root)),
        ("run_label", run_label),
        ("source_table_dir", str(table_dir(project_root, run_label))),
        ("figure1_genes", ",".join(figure1_genes)),
        ("figure1_layout", figure1_layout),
        ("figure1_study_order", ",".join(FIGURE_STUDY_ORDER)),
        ("figure2_genes", ",".join(figure2_genes)),
        ("positive_threshold_log1p_cp10k", f"{threshold:g}"),
        ("replicate_cell_cutoffs_log1p_cp10k", ",".join(cutoff_label(value) for value in replicate_cell_cutoffs)),
        ("expression_scale", "log1p(CP10K); no log10 conversion"),
        ("figure2_dot_statistic", "sample mean after applying per-cell inclusion cutoff within sample/gene"),
        ("png_dpi", str(PNG_EXPORT_DPI)),
        ("pdf_dpi", str(PDF_EXPORT_DPI)),
        ("svg_dpi", str(SVG_EXPORT_DPI)),
        ("div90_visualization_filter", "exclude current clusters 6 and 7 as stressed cells"),
        ("samarasinghe_visualization_filter", "controls only"),
        ("svg_text_policy", "editable text; matplotlib svg.fonttype none; Arial/Nimbus Sans fallback"),
    ]
    rows.extend((f"output_{idx + 1}", str(path)) for idx, path in enumerate(outputs))
    pd.DataFrame(rows, columns=["key", "value"]).to_csv(provenance / "render_manifest.tsv", sep="\t", index=False)

    files = [final_dir / "README.md", code_dir / script_path.name, provenance / "render_manifest.tsv", *outputs]
    files.extend(sorted((final_dir / "tables").glob("*.tsv")))
    sha_lines = []
    for path in files:
        if not path.exists():
            continue
        digest = subprocess.run(["sha256sum", str(path)], check=False, capture_output=True, text=True)
        if digest.returncode == 0:
            sha_lines.append(digest.stdout)
    (provenance / "sha256_manifest.txt").write_text("".join(sha_lines))


def main() -> None:
    args = parse_args()
    project_root = args.project_root
    final_dir = args.final_dir or (project_root / "final_figures" / FINAL_FOLDER_DEFAULT)
    if args.fresh and final_dir.exists():
        shutil.rmtree(final_dir)
    for subdir in ("figures/png", "figures/pdf", "figures/svg", "tables", "logs", "provenance", "code"):
        (final_dir / subdir).mkdir(parents=True, exist_ok=True)

    figure1_genes = [normalize_gene(gene) for gene in args.figure1_genes]
    figure2_genes = [normalize_gene(gene) for gene in args.figure2_genes]
    replicate_cell_cutoffs = parse_cell_cutoffs(args.replicate_cell_cutoffs)
    all_genes = list(dict.fromkeys([*figure1_genes, *figure2_genes]))

    specs = default_cross_study_marker_specs(project_root)
    raw = load_marker_expression_tables(specs, project_root, args.run_label, genes=all_genes)
    filtered, filter_summary = _apply_internal_umap_plot_filters(raw)
    filtered = add_marker_plot_coordinates(filtered)
    finite_plot = np.isfinite(pd.to_numeric(filtered["UMAP1_plot"], errors="coerce")) & np.isfinite(
        pd.to_numeric(filtered["UMAP2_plot"], errors="coerce")
    )
    filtered = filtered.loc[finite_plot].copy()

    studies = apply_figure_study_order(ordered_studies(filtered, specs))
    annotated = add_sample_metadata(filtered, studies)
    sample_summary = summarize_by_sample(annotated, all_genes, args.positive_threshold)
    study_summary = summarize_by_study(annotated, all_genes, args.positive_threshold)
    replicate_summary = summarize_replicate_expression(annotated, figure2_genes, replicate_cell_cutoffs)

    tables_dir = final_dir / "tables"
    filter_summary.to_csv(tables_dir / "integrated_marker_plot_filter_summary.tsv", sep="\t", index=False)
    sample_summary.to_csv(tables_dir / "integrated_marker_sample_summary.tsv", sep="\t", index=False)
    study_summary.to_csv(tables_dir / "integrated_marker_study_summary.tsv", sep="\t", index=False)
    replicate_summary.to_csv(tables_dir / "integrated_marker_replicate_reproducibility_summary.tsv", sep="\t", index=False)

    outputs = []
    outputs.extend(
        plot_figure1(
            annotated,
            sample_summary,
            figure1_genes,
            args.positive_threshold,
            final_dir,
            args.random_state,
            args.max_violin_cells_per_sample,
            args.figure1_layout,
        )
    )
    outputs.extend(plot_figure2_variants(replicate_summary, figure2_genes, replicate_cell_cutoffs, final_dir, args.random_state))

    manifest_rows = []
    for path in outputs:
        figure = path.stem
        replicate_cell_cutoff_token = ""
        replicate_cell_cutoff = np.nan
        if figure.startswith("this_study_div30_div90"):
            match = re.search(r"reproducibility_cell_cutoff_(none|[0-9_]+)$", figure)
            if match:
                replicate_cell_cutoff_token = match.group(1)
                replicate_cell_cutoff = np.nan if replicate_cell_cutoff_token == "none" else float(replicate_cell_cutoff_token.replace("_", "."))
        manifest_rows.append(
            {
                "output": str(path),
                "figure": figure,
                "kind": path.suffix.removeprefix("."),
                "exists": path.exists(),
                "size_bytes": path.stat().st_size if path.exists() else 0,
                "positive_threshold_log1p_cp10k": args.positive_threshold,
                "figure1_layout": args.figure1_layout,
                "replicate_cell_cutoff_token": replicate_cell_cutoff_token,
                "replicate_cell_cutoff_log1p_cp10k": replicate_cell_cutoff,
                "run_label": args.run_label,
                "source_table_dir": str(table_dir(project_root, args.run_label)),
            }
        )
    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(tables_dir / "integrated_marker_output_manifest.tsv", sep="\t", index=False)
    write_readme(
        final_dir,
        project_root,
        args.run_label,
        figure1_genes,
        figure2_genes,
        args.positive_threshold,
        replicate_cell_cutoffs,
        args.figure1_layout,
    )
    write_provenance(
        final_dir,
        project_root,
        args.run_label,
        args.positive_threshold,
        outputs,
        figure1_genes,
        figure2_genes,
        replicate_cell_cutoffs,
        args.figure1_layout,
    )

    print(manifest.to_string(index=False), flush=True)
    print("Final dir:", final_dir, flush=True)


if __name__ == "__main__":
    main()
