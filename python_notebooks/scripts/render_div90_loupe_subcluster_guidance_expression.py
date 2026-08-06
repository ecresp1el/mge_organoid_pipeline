#!/usr/bin/env python3
"""Add pooled guidance-gene expression to DIV90 cortical/subpallial Loupe subclusters."""

from __future__ import annotations

import argparse
from datetime import datetime
import hashlib
import itertools
from pathlib import Path
import shutil
import subprocess

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams.update(
    {
        "font.family": "Arial",
        "font.sans-serif": ["Arial", "Nimbus Sans", "Liberation Sans", "DejaVu Sans"],
        "font.size": 7,
        "axes.linewidth": 0.65,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }
)

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from scipy import stats

from mge_organoid_python.cross_study_marker_expression import (
    EXPRESSION_COLOR_FLOOR,
    extract_marker_expression_from_h5ad,
)
from render_div90_loupe_recluster_final_figure import (
    PANEL_LABEL_COLORS,
    RECLUSTER_SPECS,
    find_one,
    read_loupe_csv,
    short_loupe_label,
)
from render_varela_div30_lhx6_lhx8_plxna2_final_figure import (
    expression_floor_colormap,
    format_colorbar_tick,
    round_scale_limit,
)
from mge_organoid_python.cross_study_marker_expression import default_cross_study_marker_specs


PROJECT_ROOT_DEFAULT = Path("/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder")
FIGURE_ID_DEFAULT = "fig_div90_loupe_recluster_annotations_v1"
CSV_DIR_RELATIVE = Path("data/raw/loupe_exports/div90_recluster_annotations")
MEMBERSHIP_RELATIVE = Path(
    "final_figures/fig_div90_loupe_recluster_annotations_v1/"
    "tables/div90_loupe_recluster_membership.tsv.gz"
)
SCRIPT_NAME = "render_div90_loupe_subcluster_guidance_expression.py"
TABLE_PREFIX = "div90_loupe_subcluster_guidance"
GENES = [
    "NKX2-1", "LHX6", "LHX8", "PLXNA2", "NRP2", "NRP1",
    "ROBO1", "ERBB4", "CXCR4", "ACKR3", "SEMA3A", "SEMA3F",
]
GENE_DISPLAY_NAMES = {"ACKR3": "CXCR7 (ACKR3)"}
SETS = ["cortical_only", "subpallial_only"]
SET_TITLES = {
    "cortical_only": "Cortical-only Loupe subclusters",
    "subpallial_only": "Subpallial-only Loupe subclusters",
}
ATLAS_STEMS = {
    "cortical_only": "div90_loupe_cortical_guidance_feature_umaps",
    "subpallial_only": "div90_loupe_subpallial_guidance_feature_umaps",
}
DOTPLOT_STEM = "div90_loupe_guidance_subcluster_dotplots"
PAIRWISE_STEM = "div90_loupe_guidance_pairwise_coexpression"
PLXNA2_STEM = "div90_loupe_plxna2_subcluster_coexpression"


def gene_display_name(gene: str) -> str:
    return GENE_DISPLAY_NAMES.get(gene, gene)


def gene_axis_label(gene: str) -> str:
    return "CXCR7\n(ACKR3)" if gene == "ACKR3" else gene


def gene_title(gene: str) -> str:
    if gene == "NKX2-1":
        return r"$\it{NKX2}$-1"
    if gene == "ACKR3":
        return r"$\it{CXCR7}$ ($\it{ACKR3}$)"
    return rf"$\it{{{gene}}}$"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT_DEFAULT)
    parser.add_argument("--figure-id", default=FIGURE_ID_DEFAULT)
    parser.add_argument("--csv-dir", type=Path, default=None)
    parser.add_argument("--membership", type=Path, default=None)
    parser.add_argument("--dpi", type=int, default=600)
    return parser.parse_args()


def div90_spec(project_root: Path):
    matches = [
        spec
        for spec in default_cross_study_marker_specs(project_root)
        if spec.study_id == "varela_div90"
    ]
    if len(matches) != 1:
        raise ValueError(f"Expected one Varela DIV90 specification; found {len(matches)}")
    return matches[0]


def setup_dirs(outdir: Path) -> None:
    for relative in ["figures/png", "figures/pdf", "figures/svg", "tables", "code", "provenance"]:
        (outdir / relative).mkdir(parents=True, exist_ok=True)


def load_joined_data(
    project_root: Path,
    outdir: Path,
    csv_dir: Path,
    membership_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    expression_path = outdir / "tables" / f"{TABLE_PREFIX}_div90_expression.tsv.gz"
    expression, gene_matches = extract_marker_expression_from_h5ad(
        div90_spec(project_root), expression_path, project_root=project_root, genes=GENES
    )
    expression = expression[["cell_id", *GENES]].copy()
    expression["cell_id"] = expression["cell_id"].astype(str)
    if expression["cell_id"].duplicated().any():
        raise ValueError("DIV90 expression table has duplicate cell identifiers.")

    membership = pd.read_csv(membership_path, sep="\t", dtype={"cell_id": str})
    membership = membership.loc[membership["set_id"].isin(SETS)].copy()
    required = {
        "cell_id", "set_id", "set_title", "loupe_label_index", "loupe_label",
        "loupe_color", "loupe_color_original",
    }
    missing = required.difference(membership.columns)
    if missing:
        raise ValueError(f"Membership table is missing columns: {sorted(missing)}")
    if membership.duplicated(["cell_id", "set_id"]).any():
        raise ValueError("Membership table has duplicate cell/set rows.")

    spec_by_set = {str(spec["set_id"]): spec for spec in RECLUSTER_SPECS}
    frames = []
    for set_id in SETS:
        spec = spec_by_set[set_id]
        local_path = find_one(csv_dir, str(spec["recluster_pattern"]))
        local = read_loupe_csv(local_path)
        member = membership.loc[membership["set_id"] == set_id]
        joined = local.merge(member, on="cell_id", how="left", validate="one_to_one")
        joined = joined.merge(expression, on="cell_id", how="left", validate="one_to_one")
        joined["set_id"] = set_id
        joined["set_title"] = SET_TITLES[set_id]
        joined["expression_available"] = joined[GENES].notna().all(axis=1)
        frames.append(joined)
    combined = pd.concat(frames, ignore_index=True)
    qc = (
        combined.groupby(["set_id", "set_title"], sort=False)
        .agg(
            n_loupe_cells=("cell_id", "size"),
            n_unique_cells=("cell_id", "nunique"),
            n_loupe_subclusters=("loupe_label", "nunique"),
            n_cells_with_expression=("expression_available", "sum"),
        )
        .reset_index()
    )
    qc["n_cells_missing_expression"] = qc["n_loupe_cells"] - qc["n_cells_with_expression"]
    qc["missing_expression_action"] = "reported and excluded from expression summaries/plots"
    return combined, gene_matches, qc


def common_vmax(data: pd.DataFrame) -> float:
    available = data.loc[data["expression_available"]]
    q99 = [float(np.quantile(available[gene], 0.99)) for gene in GENES]
    return round_scale_limit(max(EXPRESSION_COLOR_FLOOR, max(q99)))


def coordinate_bounds(frame: pd.DataFrame) -> tuple[tuple[float, float], tuple[float, float]]:
    def one(values: pd.Series) -> tuple[float, float]:
        low, high = np.quantile(values.to_numpy(dtype=float), [0.005, 0.995])
        pad = max((high - low) * 0.04, 0.1)
        return float(low - pad), float(high + pad)

    return one(frame["loupe_x"]), one(frame["loupe_y"])


def style_umap(ax: plt.Axes, bounds: tuple[tuple[float, float], tuple[float, float]]) -> None:
    ax.set_xlim(*bounds[0])
    ax.set_ylim(*bounds[1])
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.10, 1.05, label, transform=ax.transAxes, fontsize=9, fontweight="bold", va="top")


def make_umap_colorbar_axes(cell, fig: plt.Figure) -> tuple[plt.Axes, plt.Axes]:
    subgrid = cell.subgridspec(2, 1, height_ratios=[4.0, 0.50], hspace=0.20)
    ax = fig.add_subplot(subgrid[0, 0])
    colorbar_slot = fig.add_subplot(subgrid[1, 0])
    colorbar_slot.set_axis_off()
    cax = colorbar_slot.inset_axes([0.25, 0.34, 0.50, 0.27])
    return ax, cax


def plot_context(
    ax: plt.Axes,
    key_ax: plt.Axes,
    frame: pd.DataFrame,
    set_id: str,
    bounds: tuple[tuple[float, float], tuple[float, float]],
) -> None:
    palette = PANEL_LABEL_COLORS[set_id]
    for label, subset in frame.groupby("loupe_label", sort=False):
        ax.scatter(
            subset["loupe_x"], subset["loupe_y"], s=1.4,
            color=palette.get(str(label), "#BDBDBD"), linewidths=0,
            alpha=0.86, rasterized=True,
        )
    style_umap(ax, bounds)
    ax.set_title(SET_TITLES[set_id], fontsize=7.8, fontweight="bold", pad=4)
    centers = (
        frame.groupby(["loupe_label", "loupe_label_index"], dropna=False)[["loupe_x", "loupe_y"]]
        .median()
        .reset_index()
    )
    for row in centers.itertuples(index=False):
        ax.text(
            row.loupe_x, row.loupe_y, short_loupe_label(row.loupe_label, row.loupe_label_index),
            fontsize=6.0, fontweight="bold", ha="center", va="center",
            bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="#888888", lw=0.35, alpha=0.92),
        )
    key_ax.set_axis_off()
    rows = (
        frame.groupby(["loupe_label_index", "loupe_label"], dropna=False)
        .size().reset_index(name="n_cells").sort_values("loupe_label_index")
    )
    nrows = int(np.ceil(len(rows) / 2))
    for index, row in enumerate(rows.itertuples(index=False)):
        col = index // nrows
        line = index % nrows
        color = palette.get(str(row.loupe_label), "#777777")
        rgb = tuple(int(color.lstrip("#")[i : i + 2], 16) for i in (0, 2, 4))
        text_color = "#666666" if sum(rgb) / 3 > 220 else color
        key_ax.text(
            0.01 + col * 0.51,
            0.94 - line * (0.88 / max(nrows - 1, 1)),
            f"{short_loupe_label(row.loupe_label, row.loupe_label_index)}  {row.loupe_label}",
            transform=key_ax.transAxes, color=text_color, fontsize=4.7,
            fontweight="bold", ha="left", va="center",
        )


def plot_feature(
    fig: plt.Figure,
    ax: plt.Axes,
    cax: plt.Axes,
    frame: pd.DataFrame,
    gene: str,
    bounds: tuple[tuple[float, float], tuple[float, float]],
    vmax: float,
) -> None:
    data = frame.loc[frame["expression_available"]].sort_values(gene, kind="stable")
    points = ax.scatter(
        data["loupe_x"], data["loupe_y"], c=data[gene], s=1.35,
        cmap=expression_floor_colormap(vmax), vmin=0, vmax=vmax,
        linewidths=0, rasterized=True,
    )
    style_umap(ax, bounds)
    ax.set_title(gene_title(gene), fontsize=7.5, pad=4)
    colorbar = fig.colorbar(points, cax=cax, orientation="horizontal", ticks=[0, vmax])
    colorbar.set_ticklabels([format_colorbar_tick(0), format_colorbar_tick(vmax)])
    colorbar.ax.tick_params(labelsize=4.5, length=1.2, width=0.35, pad=0.8)
    colorbar.outline.set_linewidth(0.35)
    colorbar.set_label("log1p(CP10K)", fontsize=4.3, labelpad=0.3)


def save_figure(fig: plt.Figure, outdir: Path, stem: str, dpi: int) -> list[Path]:
    outputs = []
    for extension in ["png", "pdf", "svg"]:
        output = outdir / "figures" / extension / f"{stem}.{extension}"
        fig.savefig(output, dpi=dpi, bbox_inches="tight", facecolor="white")
        outputs.append(output)
    plt.close(fig)
    return outputs


def render_feature_atlas(
    frame: pd.DataFrame,
    set_id: str,
    outdir: Path,
    vmax: float,
    dpi: int,
) -> list[Path]:
    bounds = coordinate_bounds(frame)
    fig = plt.figure(figsize=(10.2, 6.15), facecolor="white")
    grid = fig.add_gridspec(
        3, 5, left=0.043, right=0.988, bottom=0.07, top=0.95,
        wspace=0.22, hspace=0.36,
    )
    context_grid = grid[0, 0].subgridspec(2, 1, height_ratios=[4, 1.85], hspace=0.05)
    ax_context = fig.add_subplot(context_grid[0, 0])
    ax_key = fig.add_subplot(context_grid[1, 0])
    plot_context(ax_context, ax_key, frame, set_id, bounds)
    panel_label(ax_context, "a")
    for index, gene in enumerate(GENES, start=1):
        row, col = divmod(index, 5)
        ax, cax = make_umap_colorbar_axes(grid[row, col], fig)
        plot_feature(fig, ax, cax, frame, gene, bounds, vmax)
        panel_label(ax, chr(ord("a") + index))
    for index in range(len(GENES) + 1, 15):
        row, col = divmod(index, 5)
        ax = fig.add_subplot(grid[row, col])
        ax.set_axis_off()
    return save_figure(fig, outdir, ATLAS_STEMS[set_id], dpi)


def subcluster_order(frame: pd.DataFrame) -> pd.DataFrame:
    return (
        frame[["loupe_label_index", "loupe_label"]]
        .drop_duplicates().sort_values("loupe_label_index").reset_index(drop=True)
    )


def expression_summary(data: pd.DataFrame) -> pd.DataFrame:
    rows = []
    available = data.loc[data["expression_available"]]
    for (set_id, label_index, label), frame in available.groupby(
        ["set_id", "loupe_label_index", "loupe_label"], sort=False
    ):
        for gene in GENES:
            values = frame[gene].to_numpy(dtype=float)
            rows.append(
                {
                    "set_id": set_id,
                    "set_title": SET_TITLES[set_id],
                    "loupe_label_index": label_index,
                    "loupe_label": label,
                    "gene": gene,
                    "n_cells": len(values),
                    "mean_log1p_CP10K": float(np.mean(values)),
                    "median_log1p_CP10K": float(np.median(values)),
                    "pct_detected": float(100 * np.mean(values > 0)),
                    "q25_log1p_CP10K": float(np.quantile(values, 0.25)),
                    "q75_log1p_CP10K": float(np.quantile(values, 0.75)),
                    "descriptive_only": True,
                }
            )
    return pd.DataFrame(rows)


def render_dotplots(
    data: pd.DataFrame,
    summary: pd.DataFrame,
    outdir: Path,
    vmax: float,
    dpi: int,
) -> list[Path]:
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 5.0), facecolor="white")
    scatter = None
    for panel, (ax, set_id) in enumerate(zip(axes, SETS, strict=True)):
        frame = data.loc[data["set_id"] == set_id]
        order = subcluster_order(frame)
        y_map = {label: index for index, label in enumerate(order["loupe_label"])}
        subset = summary.loc[summary["set_id"] == set_id].copy()
        subset["x"] = subset["gene"].map({gene: i for i, gene in enumerate(GENES)})
        subset["y"] = subset["loupe_label"].map(y_map)
        sizes = 5 + 65 * subset["pct_detected"].to_numpy(dtype=float) / 100
        scatter = ax.scatter(
            subset["x"], subset["y"], s=sizes, c=subset["mean_log1p_CP10K"],
            cmap="Blues", vmin=0, vmax=vmax, edgecolor="#555555", linewidth=0.22,
        )
        ax.set_xticks(
            np.arange(len(GENES)), [gene_axis_label(gene) for gene in GENES],
            rotation=50, ha="right", fontsize=5.8,
        )
        counts = frame.groupby("loupe_label").size()
        labels = [f"{label} (n={int(counts[label]):,})" for label in order["loupe_label"]]
        ax.set_yticks(np.arange(len(labels)), labels, fontsize=5.6)
        ax.set_ylim(len(labels) - 0.5, -0.5)
        ax.set_xlim(-0.6, len(GENES) - 0.4)
        ax.set_title(SET_TITLES[set_id], fontsize=7.8, fontweight="bold", pad=6)
        ax.grid(False)
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
        palette = PANEL_LABEL_COLORS[set_id]
        for tick, label in zip(ax.get_yticklabels(), order["loupe_label"], strict=True):
            tick.set_color(palette.get(str(label), "#222222"))
            tick.set_fontweight("bold")
        panel_label(ax, chr(ord("a") + panel))
    colorbar = fig.colorbar(scatter, ax=axes, orientation="vertical", fraction=0.020, pad=0.025)
    colorbar.set_label("Mean expression, log1p(CP10K)", fontsize=6.2)
    colorbar.ax.tick_params(labelsize=5.5, length=2)
    size_handles = [
        axes[1].scatter(
            [], [], s=5 + 65 * pct / 100, color="#D9EAF7",
            edgecolor="#555555", linewidth=0.22, label=f"{pct}%",
        )
        for pct in [25, 50, 75, 100]
    ]
    fig.legend(
        handles=size_handles, labels=["25%", "50%", "75%", "100%"],
        title="Cells with detectable expression", frameon=False,
        fontsize=5.4, title_fontsize=5.7, loc="lower center",
        bbox_to_anchor=(0.52, 0.015), ncol=4, handletextpad=0.5, columnspacing=1.1,
    )
    fig.suptitle("DIV90 Loupe subcluster guidance-gene expression", fontsize=9.2, fontweight="bold", y=0.98)
    fig.subplots_adjust(left=0.22, right=0.88, bottom=0.21, top=0.90, wspace=0.54)
    return save_figure(fig, outdir, DOTPLOT_STEM, dpi)


def positive_high_mask(values: np.ndarray) -> tuple[np.ndarray, float]:
    detected = values[values > 0]
    if detected.size == 0:
        return np.zeros(values.size, dtype=bool), np.nan
    threshold = float(np.quantile(detected, 0.75))
    return values >= threshold, threshold


def safe_spearman(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 3 or np.std(a) == 0 or np.std(b) == 0:
        return np.nan
    return float(stats.spearmanr(a, b).statistic)


def coexpression_summary(data: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    family_rows = []
    subcluster_rows = []
    available = data.loc[data["expression_available"]]
    for set_id in SETS:
        family = available.loc[available["set_id"] == set_id]
        groups = [("pooled_family", np.nan, "All Loupe subclusters pooled", family)]
        groups += [
            ("within_subcluster", label_index, label, frame)
            for (label_index, label), frame in family.groupby(
                ["loupe_label_index", "loupe_label"], sort=False
            )
        ]
        for scope, label_index, label, frame in groups:
            for gene_a, gene_b in itertools.combinations_with_replacement(GENES, 2):
                a = frame[gene_a].to_numpy(dtype=float)
                b = frame[gene_b].to_numpy(dtype=float)
                high_a, threshold_a = positive_high_mask(a)
                high_b, threshold_b = positive_high_mask(b)
                marginal_a = float(np.mean(high_a))
                marginal_b = float(np.mean(high_b))
                observed = float(np.mean(high_a & high_b))
                expected = marginal_a * marginal_b
                enrichment = observed / expected if expected > 0 else np.nan
                row = {
                    "set_id": set_id,
                    "set_title": SET_TITLES[set_id],
                    "scope": scope,
                    "loupe_label_index": label_index,
                    "loupe_label": label,
                    "gene_a": gene_a,
                    "gene_b": gene_b,
                    "n_cells": len(frame),
                    "spearman_rho_descriptive": safe_spearman(a, b),
                    "gene_a_high_threshold_log1p_CP10K": threshold_a,
                    "gene_b_high_threshold_log1p_CP10K": threshold_b,
                    "gene_a_high_pct": 100 * marginal_a,
                    "gene_b_high_pct": 100 * marginal_b,
                    "both_high_pct": 100 * observed,
                    "both_high_expected_pct_independence": 100 * expected,
                    "high_high_enrichment_observed_over_expected": enrichment,
                    "high_high_log2_enrichment": (
                        float(np.log2(enrichment)) if enrichment > 0 else np.nan
                    ),
                    "high_definition": "top quartile among detected cells within plotted set or subcluster",
                    "p_values_computed": False,
                }
                if scope == "pooled_family":
                    family_rows.append(row)
                else:
                    subcluster_rows.append(row)
    return pd.DataFrame(family_rows), pd.DataFrame(subcluster_rows)


def symmetric_matrix(frame: pd.DataFrame, metric: str) -> np.ndarray:
    matrix = np.full((len(GENES), len(GENES)), np.nan)
    indices = {gene: i for i, gene in enumerate(GENES)}
    for row in frame.itertuples(index=False):
        i, j = indices[row.gene_a], indices[row.gene_b]
        value = float(getattr(row, metric))
        matrix[i, j] = value
        matrix[j, i] = value
    return matrix


def style_matrix_axis(ax: plt.Axes, show_y: bool) -> None:
    labels = [gene_axis_label(gene) for gene in GENES]
    ax.set_xticks(np.arange(len(GENES)), labels, rotation=50, ha="right", fontsize=5.2)
    ax.set_yticks(np.arange(len(GENES)), labels if show_y else [""] * len(GENES), fontsize=5.2)
    ax.tick_params(length=0)
    index = GENES.index("PLXNA2")
    ax.add_patch(Rectangle((-0.5, index - 0.5), len(GENES), 1, fill=False, ec="#111111", lw=0.75))
    ax.add_patch(Rectangle((index - 0.5, -0.5), 1, len(GENES), fill=False, ec="#111111", lw=0.75))
    for spine in ax.spines.values():
        spine.set_visible(False)


def render_pairwise_coexpression(
    family: pd.DataFrame, outdir: Path, dpi: int
) -> list[Path]:
    fig, axes = plt.subplots(2, 2, figsize=(8.3, 6.7), facecolor="white")
    rho_image = None
    high_image = None
    offdiag = family.loc[family["gene_a"] != family["gene_b"], "high_high_log2_enrichment"]
    high_limit = max(1.0, float(np.nanquantile(np.abs(offdiag), 0.98)))
    for col, set_id in enumerate(SETS):
        frame = family.loc[family["set_id"] == set_id]
        rho = symmetric_matrix(frame, "spearman_rho_descriptive")
        rho_image = axes[0, col].imshow(rho, cmap="coolwarm", vmin=-1, vmax=1, interpolation="none")
        axes[0, col].set_title(SET_TITLES[set_id], fontsize=7.6, fontweight="bold", pad=6)
        style_matrix_axis(axes[0, col], show_y=col == 0)
        high = symmetric_matrix(frame, "high_high_log2_enrichment")
        np.fill_diagonal(high, np.nan)
        high_image = axes[1, col].imshow(
            high, cmap="coolwarm", vmin=-high_limit, vmax=high_limit, interpolation="none"
        )
        style_matrix_axis(axes[1, col], show_y=col == 0)
        panel_label(axes[0, col], chr(ord("a") + col))
        panel_label(axes[1, col], chr(ord("c") + col))
    axes[0, 0].set_ylabel("Pooled-cell Spearman correlation", fontsize=6.5)
    axes[1, 0].set_ylabel("High–high log2(observed/expected)", fontsize=6.5)
    cbar1 = fig.colorbar(rho_image, ax=axes[0, :], orientation="vertical", fraction=0.025, pad=0.025)
    cbar1.set_label("Spearman ρ", fontsize=6)
    cbar1.ax.tick_params(labelsize=5.5, length=2)
    cbar2 = fig.colorbar(high_image, ax=axes[1, :], orientation="vertical", fraction=0.025, pad=0.025)
    cbar2.set_label("log2 enrichment", fontsize=6)
    cbar2.ax.tick_params(labelsize=5.5, length=2)
    fig.suptitle("DIV90 Loupe guidance-gene co-expression", fontsize=9.2, fontweight="bold", y=0.98)
    fig.subplots_adjust(left=0.13, right=0.88, bottom=0.13, top=0.91, wspace=0.28, hspace=0.34)
    return save_figure(fig, outdir, PAIRWISE_STEM, dpi)


def plxna2_matrix(
    frame: pd.DataFrame, set_id: str, metric: str
) -> tuple[np.ndarray, list[str], list[str]]:
    partners = [gene for gene in GENES if gene != "PLXNA2"]
    set_frame = frame.loc[
        (frame["set_id"] == set_id)
        & ((frame["gene_a"] == "PLXNA2") | (frame["gene_b"] == "PLXNA2"))
        & (frame["gene_a"] != frame["gene_b"])
    ].copy()
    set_frame["partner"] = np.where(
        set_frame["gene_a"] == "PLXNA2", set_frame["gene_b"], set_frame["gene_a"]
    )
    order = (
        set_frame[["loupe_label_index", "loupe_label"]]
        .drop_duplicates().sort_values("loupe_label_index")
    )
    matrix = np.full((len(order), len(partners)), np.nan)
    row_map = {label: i for i, label in enumerate(order["loupe_label"])}
    col_map = {gene: i for i, gene in enumerate(partners)}
    for row in set_frame.itertuples(index=False):
        matrix[row_map[row.loupe_label], col_map[row.partner]] = float(getattr(row, metric))
    labels = list(order["loupe_label"])
    return matrix, labels, partners


def style_plxna2_axis(
    ax: plt.Axes,
    labels: list[str],
    partners: list[str],
    set_id: str,
) -> None:
    ax.set_xticks(
        np.arange(len(partners)), [gene_axis_label(gene) for gene in partners],
        rotation=50, ha="right", fontsize=5.2,
    )
    ax.set_yticks(np.arange(len(labels)), labels, fontsize=5.4)
    ax.tick_params(length=0)
    palette = PANEL_LABEL_COLORS[set_id]
    for tick, label in zip(ax.get_yticklabels(), labels, strict=True):
        tick.set_color(palette.get(str(label), "#222222"))
        tick.set_fontweight("bold")
    for spine in ax.spines.values():
        spine.set_visible(False)


def render_plxna2_subcluster_coexpression(
    subcluster: pd.DataFrame, outdir: Path, dpi: int
) -> list[Path]:
    fig, axes = plt.subplots(2, 2, figsize=(8.5, 6.3), facecolor="white")
    rho_values = subcluster.loc[
        ((subcluster["gene_a"] == "PLXNA2") | (subcluster["gene_b"] == "PLXNA2"))
        & (subcluster["gene_a"] != subcluster["gene_b"]),
        "spearman_rho_descriptive",
    ].to_numpy(dtype=float)
    rho_limit = max(0.25, float(np.nanquantile(np.abs(rho_values), 0.98)))
    high_values = subcluster.loc[
        ((subcluster["gene_a"] == "PLXNA2") | (subcluster["gene_b"] == "PLXNA2"))
        & (subcluster["gene_a"] != subcluster["gene_b"]),
        "high_high_log2_enrichment",
    ].to_numpy(dtype=float)
    high_limit = max(1.0, float(np.nanquantile(np.abs(high_values), 0.98)))
    rho_image = None
    high_image = None
    for col, set_id in enumerate(SETS):
        rho, labels, partners = plxna2_matrix(subcluster, set_id, "spearman_rho_descriptive")
        rho_image = axes[0, col].imshow(
            rho, cmap="coolwarm", vmin=-rho_limit, vmax=rho_limit,
            interpolation="none", aspect="auto",
        )
        axes[0, col].set_title(SET_TITLES[set_id], fontsize=7.6, fontweight="bold", pad=6)
        style_plxna2_axis(axes[0, col], labels, partners, set_id)
        high, labels, partners = plxna2_matrix(subcluster, set_id, "high_high_log2_enrichment")
        high_image = axes[1, col].imshow(
            high, cmap="coolwarm", vmin=-high_limit, vmax=high_limit,
            interpolation="none", aspect="auto",
        )
        style_plxna2_axis(axes[1, col], labels, partners, set_id)
        panel_label(axes[0, col], chr(ord("a") + col))
        panel_label(axes[1, col], chr(ord("c") + col))
    axes[0, 0].set_ylabel("PLXNA2–partner Spearman ρ", fontsize=6.5)
    axes[1, 0].set_ylabel("PLXNA2–partner high–high enrichment", fontsize=6.5)
    cbar1 = fig.colorbar(rho_image, ax=axes[0, :], orientation="vertical", fraction=0.024, pad=0.025)
    cbar1.set_label("Spearman ρ", fontsize=6)
    cbar1.ax.tick_params(labelsize=5.5, length=2)
    cbar2 = fig.colorbar(high_image, ax=axes[1, :], orientation="vertical", fraction=0.024, pad=0.025)
    cbar2.set_label("log2(observed/expected)", fontsize=6)
    cbar2.ax.tick_params(labelsize=5.5, length=2)
    fig.suptitle(
        "PLXNA2 co-expression within DIV90 Loupe subclusters",
        fontsize=9.2, fontweight="bold", y=0.98,
    )
    fig.subplots_adjust(left=0.20, right=0.88, bottom=0.13, top=0.90, wspace=0.48, hspace=0.38)
    return save_figure(fig, outdir, PLXNA2_STEM, dpi)


def write_documentation(outdir: Path, vmax: float, qc: pd.DataFrame) -> list[Path]:
    qc_lines = "\n".join(
        f"- {row.set_title}: {int(row.n_cells_with_expression):,}/{int(row.n_loupe_cells):,} cells mapped; {int(row.n_cells_missing_expression)} missing-expression cells excluded."
        for row in qc.itertuples(index=False)
    )
    readme = f"""# DIV90 Loupe subcluster guidance-gene expression

This module extends the existing finalized Loupe recluster package without changing its original figures. It joins the cortical-only and subpallial-only Loupe local UMAP coordinates and recovered `.cloupe` cell-track labels to the existing DIV90 log1p(CP10K) expression matrix.

Genes: {', '.join(gene_display_name(gene) for gene in GENES)}.

## Mapping QC

{qc_lines}

## Figures

- `figures/png/{ATLAS_STEMS['cortical_only']}.png`: cortical Loupe subclusters followed by twelve feature plots.
- `figures/png/{ATLAS_STEMS['subpallial_only']}.png`: subpallial Loupe subclusters followed by twelve feature plots.
- `figures/png/{DOTPLOT_STEM}.png`: pooled subcluster mean expression and detection frequency.
- `figures/png/{PAIRWISE_STEM}.png`: pooled cortical/subpallial pairwise continuous and high–high co-expression.
- `figures/png/{PLXNA2_STEM}.png`: PLXNA2-focused co-expression calculated separately within every Loupe subcluster.

All feature UMAPs use unchanged Loupe local coordinates and the same 0 to {vmax:g} log1p(CP10K) scale. No sample-split panels are generated. Co-expression and subcluster summaries pool cells as requested and are descriptive; no pooled-cell inferential P values are reported. High expression is defined as the top quartile among detected cells within the plotted family or individual subcluster.
"""
    readme_path = outdir / "DIV90_LOUPE_SUBCLUSTER_GUIDANCE_README.md"
    readme_path.write_text(readme)
    caption = f"""DIV90 cortical-only and subpallial-only Loupe recluster guidance-gene analysis.

The feature atlases show the recovered Loupe subcluster annotations followed by {', '.join(gene_display_name(gene) for gene in GENES)} on the corresponding local Loupe recluster coordinates. Every feature plot uses the same 0 to {vmax:g} log1p(CP10K) color range. Dot plots show pooled subcluster mean expression by color and percentage of cells with detectable expression by point size. Pairwise co-expression panels report descriptive pooled-cell Spearman correlations and enrichment of double-high cells relative to independence. The PLXNA2-focused figure repeats these summaries separately within each recovered Loupe subcluster. No sample breakdown or pooled-cell inferential test is shown.
"""
    caption_path = outdir / "div90_loupe_subcluster_guidance_figure_captions.txt"
    caption_path.write_text(caption)
    return [readme_path, caption_path]


def git_value() -> str:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"], check=True, capture_output=True, text=True
        ).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unavailable"


def write_manifest(
    outdir: Path,
    project_root: Path,
    csv_dir: Path,
    membership_path: Path,
    outputs: list[Path],
    vmax: float,
) -> None:
    manifest = pd.DataFrame(
        [
            ("rendered_at_local", datetime.now().isoformat(timespec="seconds")),
            ("source_h5ad", str(div90_spec(project_root).resolved_h5ad_path(project_root))),
            ("source_loupe_csv_dir", str(csv_dir)),
            ("source_membership", str(membership_path)),
            ("genes", ",".join(GENES)),
            ("display_aliases", "ACKR3=CXCR7 (ACKR3)"),
            ("expression_scale", "existing log1p(CP10K)"),
            ("shared_feature_vmax", f"{vmax:g}"),
            ("sample_breakdown", "none"),
            ("coexpression_inference", "none; descriptive pooled-cell summaries"),
            ("git_commit", git_value()),
        ],
        columns=["key", "value"],
    )
    path = outdir / "provenance" / f"{TABLE_PREFIX}_render_manifest.tsv"
    manifest.to_csv(path, sep="\t", index=False)
    outputs.append(path)
    checksum = outdir / "provenance" / f"{TABLE_PREFIX}_sha256_manifest.txt"
    checksum.write_text(
        "".join(
            f"{hashlib.sha256(output.read_bytes()).hexdigest()}  {output.relative_to(outdir)}\n"
            for output in outputs
        )
    )


def main() -> None:
    args = parse_args()
    project_root = args.project_root.resolve()
    outdir = project_root / "final_figures" / args.figure_id
    csv_dir = (args.csv_dir or project_root / CSV_DIR_RELATIVE).resolve()
    membership_path = (args.membership or project_root / MEMBERSHIP_RELATIVE).resolve()
    setup_dirs(outdir)
    data, gene_matches, qc = load_joined_data(
        project_root, outdir, csv_dir, membership_path
    )
    vmax = common_vmax(data)
    expression = expression_summary(data)
    family_coexpression, subcluster_coexpression = coexpression_summary(data)

    outputs: list[Path] = []
    for set_id in SETS:
        outputs += render_feature_atlas(
            data.loc[data["set_id"] == set_id], set_id, outdir, vmax, args.dpi
        )
    outputs += render_dotplots(data, expression, outdir, vmax, args.dpi)
    outputs += render_pairwise_coexpression(family_coexpression, outdir, args.dpi)
    outputs += render_plxna2_subcluster_coexpression(
        subcluster_coexpression, outdir, args.dpi
    )

    tables = {
        f"{TABLE_PREFIX}_gene_matches.tsv": gene_matches,
        f"{TABLE_PREFIX}_mapping_qc.tsv": qc,
        f"{TABLE_PREFIX}_joined_per_cell.tsv.gz": data[
            [
                "cell_id", "set_id", "set_title", "loupe_label_index", "loupe_label",
                "loupe_color", "loupe_x", "loupe_y", "expression_available", *GENES,
            ]
        ],
        f"{TABLE_PREFIX}_subcluster_expression_summary.tsv": expression,
        f"{TABLE_PREFIX}_family_pairwise_coexpression.tsv": family_coexpression,
        f"{TABLE_PREFIX}_subcluster_pairwise_coexpression.tsv": subcluster_coexpression,
    }
    for name, frame in tables.items():
        path = outdir / "tables" / name
        frame.to_csv(path, sep="\t", index=False)
        outputs.append(path)
    outputs += write_documentation(outdir, vmax, qc)
    source = Path(__file__).resolve()
    copied = outdir / "code" / SCRIPT_NAME
    shutil.copy2(source, copied)
    outputs.append(copied)
    write_manifest(outdir, project_root, csv_dir, membership_path, outputs, vmax)

    print(f"Output folder: {outdir}")
    print(f"Genes: {', '.join(gene_display_name(gene) for gene in GENES)}")
    print(qc.to_string(index=False))
    print(f"Shared feature range: 0 to {vmax:g} log1p(CP10K)")
    for stem in [*ATLAS_STEMS.values(), DOTPLOT_STEM, PAIRWISE_STEM, PLXNA2_STEM]:
        print(outdir / "figures" / "png" / f"{stem}.png")


if __name__ == "__main__":
    main()
