#!/usr/bin/env python3
"""Cross-study marker-expression summaries and figures."""
from __future__ import annotations
import argparse
from pathlib import Path
import re
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none", "pdf.fonttype": 42})
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
PNG_EXPORT_DPI = 450
SVG_EXPORT_DPI = 450
PDF_EXPORT_DPI = 300
PREVIEW_DPI = None

CROSS_STUDY_GENES_DEFAULT = ["NKX2.1", "LHX6", "LHX8", "CRABP1", "TAC1", "VIPR2", "SLC6A1", "ZEB2", "SST", "PVALB"]

REPLICATE_GENES_DEFAULT = ["NKX2.1", "LHX6", "LHX8"]

CROSS_STUDY_STEM = "cross_study_canonical_mge_marker_expression_violin_positive_matrix"

REPLICATE_STEM = "this_study_div30_div90_mge_marker_replicate_reproducibility"

REPLICATE_VARIANT_STEM_TEMPLATE = "this_study_div30_div90_mge_marker_replicate_reproducibility_cell_cutoff_{cutoff_token}"

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
STUDY_COLORS = {
    "varela_div30": "#2B6CB0",
    "varela_div90": "#2F855A",
    "siebert_2026": "#805AD5",
    "walsh": "#DD6B20",
    "bershteyn_2025": "#B83280",
    "bershteyn_2023": "#319795",
    "samarasinghe_2021": "#4A5568",
}

BERSHTEYN_2023_SAMPLE_INFO = {
    "D0": ("DIV0 hESC", 0.0, 0.0),
    "D14": ("DIV14 NPC", 14.0, 0.0),
    "MB460": ("DIV42 EOP L1 U", 42.0, 1.1),
    "MB461": ("DIV42 EOP L1 S", 42.0, 1.2),
    "MB279": ("DIV42 EOP L2 U", 42.0, 2.1),
    "MB280": ("DIV42 EOP L2 S", 42.0, 2.2),
    "MB527": ("DIV42 EOP L3 U", 42.0, 3.1),
    "MB528": ("DIV42 EOP L3 S", 42.0, 3.2),
    "MS35mock": ("DIV42 EOP B4 U (MS35mock)", 42.0, 4.0),
    "MS35r41": ("DIV42 EOP B5 U (MS35r41)", 42.0, 5.0),
    "mockv2dw": ("DIV42 EOP B6 U (mockv2dw)", 42.0, 6.0),
    "mockv2ym": ("DIV42 EOP B7 U (mockv2ym)", 42.0, 7.0),
    "r41v2dw": ("DIV42 EOP B8 U (r41v2dw)", 42.0, 8.0),
    "r41v2ym": ("DIV42 EOP B9 U (r41v2ym)", 42.0, 9.0),
}

BERSHTEYN_2025_SAMPLE_INFO = {
    "200319S": ("DIV unknown EOP batch 1 (200319S)", 999.0, 1.0),
    "010519S1": ("DIV unknown EOP batch 2 (010519S1)", 999.0, 2.0),
    "010519S2": ("DIV unknown EOP batch 3 (010519S2)", 999.0, 3.0),
    "111219S": ("DIV unknown EOP batch 4 (111219S)", 999.0, 4.0),
    "251219S": ("DIV unknown EOP batch 5 (251219S)", 999.0, 5.0),
    "150120S": ("DIV unknown EOP batch 6 (150120S)", 999.0, 6.0),
    "200520S1": ("DIV unknown EOP batch 7 (200520S1)", 999.0, 7.0),
    "200520S2": ("DIV unknown EOP batch 8 (200520S2)", 999.0, 8.0),
    "100620S": ("DIV unknown EOP batch 9 (100620S)", 999.0, 9.0),
    "010720S": ("DIV unknown EOP batch 10 (010720S)", 999.0, 10.0),
    "220720S1": ("DIV unknown EOP batch 11 (220720S1)", 999.0, 11.0),
    "220720S2": ("DIV unknown EOP batch 12 (220720S2)", 999.0, 12.0),
    "070120S": ("DIV unknown EOP batch 13 (070120S)", 999.0, 13.0),
    "280120S": ("DIV unknown EOP batch 14 (280120S)", 999.0, 14.0),
}
def gene_label(gene: str) -> str:
    return GENE_LABELS.get(gene, gene)


def finite_values(data: pd.DataFrame, gene: str) -> np.ndarray:
    values = pd.to_numeric(data[gene], errors="coerce").to_numpy(dtype=float)
    return values[np.isfinite(values)]


def save_figure(fig: plt.Figure, final_dir: Path, stem: str) -> list[Path]:
    outputs = []
    formats = (("png", PREVIEW_DPI),) if PREVIEW_DPI else (
        ("png", PNG_EXPORT_DPI), ("pdf", PDF_EXPORT_DPI), ("svg", SVG_EXPORT_DPI)
    )
    for ext, dpi in formats:
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
    """Order studies for display."""
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


def plot_cross_study_expression(
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
    return save_figure(fig, final_dir, CROSS_STUDY_STEM)


def plot_replicate_expression(
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
    stem = REPLICATE_VARIANT_STEM_TEMPLATE.format(cutoff_token=token)
    return save_figure(fig, final_dir, stem)


def plot_replicate_variants(
    replicate_summary: pd.DataFrame,
    figure_genes: list[str],
    cell_cutoffs: list[float | None],
    final_dir: Path,
    random_state: int,
) -> list[Path]:
    outputs = []
    for cutoff in cell_cutoffs:
        outputs.extend(plot_replicate_expression(replicate_summary, figure_genes, cutoff, final_dir, random_state))
    return outputs
def normalize_gene(value: str) -> str:
    text = str(value).strip()
    aliases = {
        "PV": "PVALB",
        "PVALB": "PVALB",
        "NKX2.1": "NKX2-1",
        "NKX2-1": "NKX2-1",
        "NKX2_1": "NKX2-1",
    }
    return aliases.get(text.upper(), text.upper())


def first_number(pattern: str, text: str, default: float = 9999.0) -> float:
    match = re.search(pattern, text, flags=re.IGNORECASE)
    if not match:
        return default
    return float(match.group(1))


def sample_display_and_order(study_id: str, sample: object) -> tuple[str, float, float]:
    sample_text = str(sample)
    if study_id == "bershteyn_2023" and sample_text in BERSHTEYN_2023_SAMPLE_INFO:
        return BERSHTEYN_2023_SAMPLE_INFO[sample_text]
    if study_id == "bershteyn_2025" and sample_text in BERSHTEYN_2025_SAMPLE_INFO:
        return BERSHTEYN_2025_SAMPLE_INFO[sample_text]
    if study_id == "varela_div30":
        return f"DIV30 {sample_text}", 30.0, first_number(r"MW-(\d+)", sample_text)
    if study_id == "varela_div90":
        return f"DIV90 {sample_text}", 90.0, first_number(r"MW-(\d+)", sample_text)
    if study_id == "siebert_2026":
        age_order = 51.0 if sample_text.startswith("Young") else 164.0 if sample_text.startswith("Old") else 9999.0
        div_label = f"likely DIV{int(age_order)}" if age_order < 9999 else "DIV unknown"
        return f"{div_label} {sample_text}", age_order, first_number(r"_(\d+)", sample_text)
    if study_id == "walsh":
        within_order = 1.0 if "9671" in sample_text else 2.0 if "9672" in sample_text else first_number(r"GSM(\d+)", sample_text)
        label = "MEL1 dFB d75" if "9671" in sample_text else "MEL1 vFB d75" if "9672" in sample_text else sample_text
        return f"DIV75 {label}", 75.0, within_order
    if study_id == "samarasinghe_2021":
        div = first_number(r"D(\d+)", sample_text)
        compact = sample_text.replace("_seurat", "").replace("_", " ")
        return f"DIV{int(div)} {compact}" if div < 9999 else f"DIV unknown {compact}", div, first_number(r"_(\d+)_seurat", sample_text)
    return sample_text, first_number(r"D(?:IV)?\s*([0-9]+)", sample_text), first_number(r"(\d+)", sample_text)


def add_sample_metadata(data: pd.DataFrame, studies: pd.DataFrame) -> pd.DataFrame:
    sample_meta = data[["study_id", "study_label", "sample"]].drop_duplicates().copy()
    display = [sample_display_and_order(row.study_id, row.sample) for row in sample_meta.itertuples(index=False)]
    sample_meta["sample_display_label"] = [value[0] for value in display]
    sample_meta["sample_age_order"] = [value[1] for value in display]
    sample_meta["sample_order_within_study"] = [value[2] for value in display]
    sample_meta["sample_plot_id"] = sample_meta["study_id"].astype(str) + "::" + sample_meta["sample"].astype(str)
    sample_meta = sample_meta.merge(studies[["study_id", "study_plot_order"]], on="study_id", how="left")
    return data.merge(sample_meta, on=["study_id", "study_label", "sample"], how="left", validate="many_to_one")
def add_marker_plot_coordinates(data: pd.DataFrame) -> pd.DataFrame:
    """Add plotting-only UMAP coordinates without modifying the source embedding."""
    out = data.copy()
    out["UMAP1_plot"] = pd.to_numeric(out["umap_1"], errors="coerce")
    out["UMAP2_plot"] = pd.to_numeric(out["umap_2"], errors="coerce")
    div90 = out["study_id"].astype(str).eq("varela_div90")
    out.loc[div90, "UMAP2_plot"] = -1.0 * out.loc[div90, "UMAP2_plot"]
    return out


def _finite_quantile(values: np.ndarray, quantile: float) -> float:
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.nan
    return float(np.quantile(values, quantile))


def _apply_internal_umap_plot_filters(data: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Apply hard-coded biological filters that define what each UMAP row shows."""
    summaries = []
    filtered_parts = []
    for study_id, group in data.groupby("study_id", sort=False):
        rules: list[str] = []
        labels: list[str] = []
        kept = group
        if study_id == "samarasinghe_2021":
            rules.append("samarasinghe_2021_controls_only")
            labels.append("Controls only")
            control = group["sample"].astype(str).str.contains("Ctrl", case=False, na=False)
            kept = group.loc[control].copy()
            if kept.empty:
                raise ValueError("Samarasinghe control-only UMAP filter removed all cells.")
            non_control_samples = sorted(kept.loc[~kept["sample"].astype(str).str.contains("Ctrl", case=False, na=False), "sample"].unique())
            if non_control_samples:
                raise ValueError("Samarasinghe control-only UMAP filter retained non-control samples: " + ", ".join(non_control_samples))
        if study_id == "varela_div90":
            rules.append("varela_div90_exclude_current_clusters_6_7_stressed_cells")
            labels.append("Stressed cells removed")
            current_cluster = pd.to_numeric(kept["cluster"], errors="coerce")
            kept = kept.loc[~current_cluster.isin([6, 7])].copy()
            if kept.empty:
                raise ValueError("DIV90 stressed-cell plot filter removed all cells.")
        rule = "+".join(rules) if rules else "none"
        summaries.append(
            {
                "study_id": study_id,
                "study_label": str(group["study_label"].iloc[0]),
                "plot_filter": rule,
                "plot_filter_label": "; ".join(labels),
                "n_cells_before_plot_filter": int(group.shape[0]),
                "n_cells_after_plot_filter": int(kept.shape[0]),
                "n_cells_removed_by_plot_filter": int(group.shape[0] - kept.shape[0]),
                "sample_values_after_plot_filter": ";".join(sorted(kept["sample"].astype(str).unique())),
            }
        )
        filtered_parts.append(kept)
    filtered = pd.concat(filtered_parts, ignore_index=True)
    summary = pd.DataFrame(summaries)
    filtered_total = int(summary["n_cells_after_plot_filter"].sum())
    if filtered_total != filtered.shape[0]:
        raise ValueError(f"Internal UMAP plot filter count mismatch: summary={filtered_total}, data={filtered.shape[0]}")
    return filtered, summary

def main():
    global PREVIEW_DPI
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, nargs="+", required=True, help="Cell tables: TSV or TSV.gz")
    parser.add_argument("--outdir", type=Path, required=True, help="New output folder")
    parser.add_argument("--positive-threshold", type=float, default=0.5)
    parser.add_argument("--genes", nargs="+", default=CROSS_STUDY_GENES_DEFAULT)
    parser.add_argument("--random-state", type=int, default=1729)
    parser.add_argument("--max-violin-cells-per-sample", type=int, default=1800)
    parser.add_argument("--preview", action="store_true", help="Write 80 dpi PNG previews only")
    args = parser.parse_args()
    PREVIEW_DPI = 80 if args.preview else None
    if not np.isfinite(args.positive_threshold) or args.positive_threshold < 0:
        parser.error("--positive-threshold must be finite and nonnegative")
    if args.max_violin_cells_per_sample < 1:
        parser.error("--max-violin-cells-per-sample must be positive")
    if args.outdir.exists():
        parser.error("--outdir must be a new folder")
    genes = list(dict.fromkeys(normalize_gene(g) for g in args.genes))
    data = pd.concat([pd.read_csv(p, sep="\t", dtype={"cell_id": str, "sample": str, "study_id": str})
                      for p in args.input], ignore_index=True)
    required = {"cell_id", "study_id", "study_label", "sample", "cluster", "umap_1", "umap_2", *genes}
    missing = required - set(data)
    if missing:
        parser.error("Missing input columns: " + ", ".join(sorted(missing)))
    if data.empty or data[["cell_id", "study_id", "study_label", "sample"]].isna().any().any():
        parser.error("Input must have cells and nonmissing identifiers")
    if data.duplicated(["study_id", "cell_id"]).any():
        parser.error("Duplicate cell IDs within a study")
    if (data.groupby("study_id")["study_label"].nunique() > 1).any():
        parser.error("Each study_id must have one study_label")
    data, filters = _apply_internal_umap_plot_filters(data)
    data = add_marker_plot_coordinates(data)
    finite = np.isfinite(data["UMAP1_plot"]) & np.isfinite(data["UMAP2_plot"])
    geometry = data.assign(finite_embedding=finite).groupby("study_id").agg(
        n_before_finite_embedding=("finite_embedding", "size"),
        n_after_finite_embedding=("finite_embedding", "sum")).reset_index()
    data = data.loc[finite].copy()
    if data.empty:
        parser.error("No cells remain after filtering and requiring finite embeddings")
    studies = data[["study_id", "study_label"]].drop_duplicates().copy()
    studies["study_plot_order"] = np.arange(len(studies))
    studies = apply_figure_study_order(studies)
    data = add_sample_metadata(data, studies)
    samples = summarize_by_sample(data, genes, args.positive_threshold)
    study_summary = summarize_by_study(data, genes, args.positive_threshold)
    tables = args.outdir / "tables"
    tables.mkdir(parents=True)
    filters.merge(geometry, on="study_id", how="left").to_csv(tables / "cell_filters.tsv", sep="\t", index=False)
    samples.to_csv(tables / "sample_expression.tsv", sep="\t", index=False)
    study_summary.to_csv(tables / "study_expression.tsv", sep="\t", index=False)
    plot_cross_study_expression(data, samples, genes, args.positive_threshold, args.outdir,
                         args.random_state, args.max_violin_cells_per_sample)
    replicate_genes = [normalize_gene(g) for g in REPLICATE_GENES_DEFAULT if normalize_gene(g) in genes]
    if replicate_genes and data["study_id"].isin(["varela_div30", "varela_div90"]).any():
        cutoffs = [None, 0.5, 1.0]
        replicates = summarize_replicate_expression(data, replicate_genes, cutoffs)
        replicates.to_csv(tables / "replicate_expression.tsv", sep="\t", index=False)
        plot_replicate_variants(replicates, replicate_genes, cutoffs, args.outdir, args.random_state)
    pd.DataFrame({"parameter": ["positive_threshold", "random_state", "max_violin_cells_per_sample"],
                  "value": [args.positive_threshold, args.random_state, args.max_violin_cells_per_sample]}).to_csv(
                      tables / "figure_parameters.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
