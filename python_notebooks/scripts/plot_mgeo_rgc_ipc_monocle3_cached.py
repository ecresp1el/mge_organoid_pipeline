#!/usr/bin/env python3
"""Render Monocle3 MGEO RGC/IPC validation plots from cached flat files.

This script intentionally does not call Monocle, recompute UMAP, learn a graph,
or order cells. It reads plot-ready CSVs created by
scripts/export_mgeo_rgc_ipc_monocle3_plot_cache.R.
"""

from __future__ import annotations

import math
import os
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


PROJECT_ROOT = Path(
    os.environ.get(
        "PROJECT_ROOT",
        "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder",
    )
)
RUN_ROOT = PROJECT_ROOT / "results" / "mgeo_rgc_ipc_monocle3"
INTERACTIVE_DIR = RUN_ROOT / "interactive"
CACHE_DIR = INTERACTIVE_DIR / "plot_cache"
PLOT_DIR = INTERACTIVE_DIR / "plots"
PLOT_DIR.mkdir(parents=True, exist_ok=True)

GENES_TO_PLOT = [
    "NES",
    "VIM",
    "HES1",
    "ASCL1",
    "DLX2",
    "DCX",
    "STMN2",
    "FABP7",
    "HES5",
    "SOX2",
    "DLX1",
    "NKX2-1",
]

# Visualization threshold only. Values below this remain gray, matching the
# plotting style used in 01_seurat_to_anndata.ipynb.
FEATURE_VMIN_EXPRESSION = float(os.environ.get("FEATURE_VMIN_EXPRESSION", "2.0"))
VMAX_PERCENTILE = float(os.environ.get("VMAX_PERCENTILE", "99"))


def rounded_vmax(values: np.ndarray, percentile: float = 99, minimum: float = 1) -> float:
    finite_values = values[np.isfinite(values)]
    positive_values = finite_values[finite_values > 0]
    if positive_values.size == 0:
        return float(minimum)

    vmax = float(np.ceil(np.percentile(positive_values, percentile)))
    return float(max(vmax, minimum))


def make_thresholded_colorbar_cmap(
    expression_cmap: mpl.colors.Colormap,
    background_color: str = "#d9d9d9",
    n_colors: int = 256,
    gray_fraction: float = 0.18,
) -> mpl.colors.ListedColormap:
    background_rgba = mpl.colors.to_rgba(background_color)
    colors = []
    gray_fraction = min(max(gray_fraction, 0.0), 0.95)

    for i in range(n_colors):
        t = i / (n_colors - 1)
        if t < gray_fraction and gray_fraction > 0:
            local_t = t / gray_fraction
            expr_color = expression_cmap(0.0)
            color = tuple(
                (1 - local_t) * background_rgba[j] + local_t * expr_color[j]
                for j in range(4)
            )
        else:
            expr_t = (t - gray_fraction) / (1 - gray_fraction) if gray_fraction < 1 else 1.0
            color = expression_cmap(min(max(expr_t, 0.0), 1.0))
        colors.append(color)

    return mpl.colors.ListedColormap(colors)


def clean_axis(ax: plt.Axes) -> None:
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(False)


def draw_umap_background(
    ax: plt.Axes,
    umap: np.ndarray,
    background_color: str,
    point_size_background: float,
    background_alpha: float,
) -> None:
    ax.scatter(
        umap[:, 0],
        umap[:, 1],
        s=point_size_background,
        c=background_color,
        alpha=background_alpha,
        linewidths=0,
        edgecolors="none",
        rasterized=True,
    )


def add_vertical_colorbar(
    fig: plt.Figure,
    ax: plt.Axes,
    cmap: mpl.colors.Colormap,
    vmax: float,
    colorbar_height: str,
    colorbar_width: str,
    colorbar_labelsize: int,
    threshold: float | None = None,
) -> None:
    norm = mpl.colors.Normalize(vmin=0, vmax=vmax)
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    cax = inset_axes(
        ax,
        width=colorbar_width,
        height=colorbar_height,
        loc="center right",
        bbox_to_anchor=(0.06, 0, 1, 1),
        bbox_transform=ax.transAxes,
        borderpad=0,
    )
    cbar = fig.colorbar(sm, cax=cax, orientation="vertical")
    cbar.set_ticks([0, vmax])
    cbar.set_ticklabels(["0", f"{int(vmax)}"])
    cbar.ax.tick_params(labelsize=colorbar_labelsize, length=2, width=0.5, pad=1)
    cbar.outline.set_linewidth(0.4)

    if threshold is not None:
        cbar.ax.axhline(threshold, color="black", linewidth=0.4, alpha=0.6)


def plot_gene_panel(
    fig: plt.Figure,
    ax: plt.Axes,
    umap: np.ndarray,
    expression: np.ndarray,
    title: str,
    *,
    vmax_percentile: float = 99,
    vmin_expression: float = 2.0,
    cmap: str = "YlOrRd",
    background_color: str = "#d9d9d9",
    point_size_background: float = 0.8,
    point_size_expression: float = 1.0,
    background_alpha: float = 0.35,
    expression_alpha: float = 0.9,
    colorbar_height: str = "33%",
    colorbar_width: str = "3%",
    colorbar_labelsize: int = 14,
) -> None:
    expression_cmap = plt.get_cmap(cmap)
    vmax = rounded_vmax(
        expression,
        percentile=vmax_percentile,
        minimum=max(1, vmin_expression),
    )
    if vmax <= vmin_expression:
        vmax = float(np.ceil(vmin_expression + 1))

    draw_umap_background(
        ax,
        umap,
        background_color,
        point_size_background,
        background_alpha,
    )

    expressing = np.isfinite(expression) & (expression >= vmin_expression)
    if expressing.any():
        ax.scatter(
            umap[expressing, 0],
            umap[expressing, 1],
            s=point_size_expression,
            c=expression[expressing],
            cmap=expression_cmap,
            norm=mpl.colors.Normalize(vmin=vmin_expression, vmax=vmax),
            alpha=expression_alpha,
            linewidths=0,
            edgecolors="none",
            rasterized=True,
        )

    ax.set_title(title.upper(), fontstyle="italic", fontweight="normal", fontsize=11, pad=4)
    clean_axis(ax)

    threshold_fraction = min(max(vmin_expression / vmax, 0.0), 0.95)
    colorbar_cmap = make_thresholded_colorbar_cmap(
        expression_cmap,
        background_color=background_color,
        gray_fraction=threshold_fraction,
    )
    add_vertical_colorbar(
        fig,
        ax,
        colorbar_cmap,
        vmax,
        colorbar_height,
        colorbar_width,
        colorbar_labelsize,
        threshold=vmin_expression,
    )


def plot_pseudotime_panel(
    fig: plt.Figure,
    ax: plt.Axes,
    umap: np.ndarray,
    pseudotime: np.ndarray,
    *,
    vmax_percentile: float = 99,
    cmap: str = "YlOrRd",
    background_color: str = "#d9d9d9",
    point_size_background: float = 0.8,
    point_size_expression: float = 1.0,
    background_alpha: float = 0.35,
    expression_alpha: float = 0.9,
    colorbar_height: str = "33%",
    colorbar_width: str = "3%",
    colorbar_labelsize: int = 14,
) -> None:
    expression_cmap = plt.get_cmap(cmap)
    finite = np.isfinite(pseudotime)
    vmax = rounded_vmax(pseudotime[finite], percentile=vmax_percentile, minimum=1)

    draw_umap_background(
        ax,
        umap,
        background_color,
        point_size_background,
        background_alpha,
    )

    if finite.any():
        ax.scatter(
            umap[finite, 0],
            umap[finite, 1],
            s=point_size_expression,
            c=np.clip(pseudotime[finite], 0, vmax),
            cmap=expression_cmap,
            norm=mpl.colors.Normalize(vmin=0, vmax=vmax),
            alpha=expression_alpha,
            linewidths=0,
            edgecolors="none",
            rasterized=True,
        )

    ax.set_title("pseudotime", fontweight="normal", fontsize=11, pad=4)
    clean_axis(ax)
    add_vertical_colorbar(
        fig,
        ax,
        expression_cmap,
        vmax,
        colorbar_height,
        colorbar_width,
        colorbar_labelsize,
        threshold=None,
    )


def load_partition_frame(partition_label: str) -> pd.DataFrame:
    plot_ready_path = CACHE_DIR / f"mgeo_rgc_ipc_monocle3_plot_ready_partition_{partition_label}.csv"
    expression_path = CACHE_DIR / "mgeo_rgc_ipc_marker_expression_log_normalized.csv"
    if not plot_ready_path.exists():
        raise FileNotFoundError(plot_ready_path)
    if not expression_path.exists():
        raise FileNotFoundError(expression_path)

    plot_df = pd.read_csv(plot_ready_path)
    expr_df = pd.read_csv(expression_path)
    merged = plot_df.merge(expr_df, on="cell_id", how="left", validate="one_to_one")
    return merged


def make_umap_grid(
    plot_df: pd.DataFrame,
    genes: list[str],
    output_prefix: Path,
    *,
    include_pseudotime: bool,
    ncols: int = 3,
    figsize_per_panel: tuple[float, float] = (3.0, 2.8),
) -> None:
    panel_count = len(genes) + int(include_pseudotime)
    nrows = int(math.ceil(panel_count / ncols))
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(figsize_per_panel[0] * ncols, figsize_per_panel[1] * nrows),
        facecolor="white",
        constrained_layout=False,
    )
    axes = np.asarray(axes).reshape(-1)
    for ax in axes[panel_count:]:
        ax.axis("off")

    umap = plot_df[["monocle3_umap_1", "monocle3_umap_2"]].to_numpy()
    axis_index = 0

    if include_pseudotime:
        pseudotime = plot_df["rgc_ipc_pseudotime"].to_numpy(dtype=float)
        plot_pseudotime_panel(
            fig,
            axes[axis_index],
            umap,
            pseudotime,
            vmax_percentile=VMAX_PERCENTILE,
        )
        axis_index += 1

    for gene in genes:
        if gene not in plot_df.columns:
            raise KeyError(f"{gene} not found in marker expression cache")
        expression = plot_df[gene].to_numpy(dtype=float)
        plot_gene_panel(
            fig,
            axes[axis_index],
            umap,
            expression,
            gene,
            vmax_percentile=VMAX_PERCENTILE,
            vmin_expression=FEATURE_VMIN_EXPRESSION,
        )
        axis_index += 1

    fig.subplots_adjust(left=0.02, right=0.96, top=0.94, bottom=0.03, wspace=0.18, hspace=0.22)
    png_path = output_prefix.with_suffix(".png")
    pdf_path = output_prefix.with_suffix(".pdf")
    fig.savefig(png_path, dpi=300)
    fig.savefig(pdf_path)
    plt.close(fig)
    print(f"Saved: {png_path}")
    print(f"Saved: {pdf_path}")


def make_pseudotime_expression_grid(
    plot_df: pd.DataFrame,
    genes: list[str],
    output_prefix: Path,
    ncols: int = 3,
    figsize_per_panel: tuple[float, float] = (3.2, 2.8),
) -> None:
    finite = np.isfinite(plot_df["rgc_ipc_pseudotime"].to_numpy(dtype=float))
    df = plot_df.loc[finite].copy()
    df = df.sort_values("rgc_ipc_pseudotime")

    nrows = int(math.ceil(len(genes) / ncols))
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(figsize_per_panel[0] * ncols, figsize_per_panel[1] * nrows),
        facecolor="white",
        constrained_layout=False,
    )
    axes = np.asarray(axes).reshape(-1)
    for ax in axes[len(genes):]:
        ax.axis("off")

    x = df["rgc_ipc_pseudotime"].to_numpy(dtype=float)
    for ax, gene in zip(axes, genes):
        y = df[gene].to_numpy(dtype=float)
        ax.scatter(x, y, s=0.15, color="#bdbdbd", alpha=0.05, linewidths=0, rasterized=True)

        bins = np.linspace(np.nanmin(x), np.nanmax(x), 180)
        bin_index = np.digitize(x, bins)
        x_smooth = []
        y_smooth = []
        for idx in range(1, len(bins)):
            in_bin = bin_index == idx
            if in_bin.sum() < 20:
                continue
            x_smooth.append(float(np.nanmedian(x[in_bin])))
            y_smooth.append(float(np.nanmean(y[in_bin])))
        if len(x_smooth) > 3:
            ax.plot(x_smooth, y_smooth, color="firebrick", linewidth=1.1)

        ax.set_title(gene.upper(), fontstyle="italic", fontweight="normal", fontsize=11, pad=4)
        ax.grid(False)
        ax.tick_params(axis="both", labelsize=8, length=2, width=0.5)
        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)

    fig.supxlabel("Monocle3 pseudotime", fontsize=11)
    fig.supylabel("log-normalized expression", fontsize=11)
    fig.subplots_adjust(left=0.08, right=0.98, top=0.93, bottom=0.08, wspace=0.28, hspace=0.34)
    png_path = output_prefix.with_suffix(".png")
    pdf_path = output_prefix.with_suffix(".pdf")
    fig.savefig(png_path, dpi=300)
    fig.savefig(pdf_path)
    plt.close(fig)
    print(f"Saved: {png_path}")
    print(f"Saved: {pdf_path}")


def main() -> None:
    print(f"PROJECT_ROOT={PROJECT_ROOT}")
    print(f"CACHE_DIR={CACHE_DIR}")
    print(f"PLOT_DIR={PLOT_DIR}")
    print(f"FEATURE_VMIN_EXPRESSION={FEATURE_VMIN_EXPRESSION}")
    print(f"VMAX_PERCENTILE={VMAX_PERCENTILE}")

    for partition_label in ["true", "false"]:
        print(f"Rendering partition_{partition_label}")
        plot_df = load_partition_frame(partition_label)
        available_genes = [gene for gene in GENES_TO_PLOT if gene in plot_df.columns]
        if not available_genes:
            raise RuntimeError("No requested genes were available in the marker expression cache")

        make_umap_grid(
            plot_df,
            available_genes,
            PLOT_DIR / f"monocle3_marker_umap_grid_python_style_partition_{partition_label}",
            include_pseudotime=True,
        )
        make_umap_grid(
            plot_df,
            available_genes,
            PLOT_DIR / f"monocle3_marker_gene_umap_grid_python_style_partition_{partition_label}",
            include_pseudotime=False,
        )
        make_umap_grid(
            plot_df,
            [],
            PLOT_DIR / f"monocle3_pseudotime_umap_python_style_partition_{partition_label}",
            include_pseudotime=True,
        )
        make_pseudotime_expression_grid(
            plot_df,
            available_genes,
            PLOT_DIR / f"monocle3_marker_pseudotime_grid_python_style_partition_{partition_label}",
        )


if __name__ == "__main__":
    main()
