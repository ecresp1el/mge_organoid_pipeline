"""Plot helpers for Notebook 00 with deterministic save locations."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


@dataclass(frozen=True)
class PlotConfig:
    """Notebook plot output behavior."""

    output_dir: Path
    show: bool = True
    save: bool = True
    dpi: int = 160

    @classmethod
    def from_root(cls, data_root: Path | str, run_label: str, show: bool = True, save: bool = True) -> "PlotConfig":
        output_dir = Path(data_root).expanduser().resolve() / "results" / "notebook00" / run_label / "plots"
        return cls(output_dir=output_dir, show=show, save=save)


def _safe_name(value: object) -> str:
    safe = str(value).strip().replace(" ", "_").replace("/", "_")
    return "".join(ch for ch in safe if ch.isalnum() or ch in {"_", "-", "."})


def save_current_figure(plot_config: PlotConfig, name: str):
    """Save the current Matplotlib figure and optionally display it."""
    out_path = None
    if plot_config.save:
        plot_config.output_dir.mkdir(parents=True, exist_ok=True)
        out_path = plot_config.output_dir / f"{_safe_name(name)}.png"
        plt.gcf().savefig(out_path, dpi=plot_config.dpi, bbox_inches="tight")

    if plot_config.show:
        plt.show()
    else:
        plt.close()
    return out_path


def plot_source_availability(source_table: pd.DataFrame, plot_config: PlotConfig, name: str = "source_availability"):
    """Plot requested sample source availability."""
    plot_df = source_table.copy()
    plot_df["is_available"] = plot_df["load_status"].eq("available")
    colors = plot_df["is_available"].map({True: "#2F855A", False: "#C53030"}).tolist()

    fig_width = max(7, 1.15 * len(plot_df))
    fig, ax = plt.subplots(figsize=(fig_width, 3.8))
    ax.bar(plot_df["run_sample_id"], [1] * len(plot_df), color=colors)
    ax.set_ylim(0, 1.25)
    ax.set_ylabel("source status")
    ax.set_xlabel("run_sample_id")
    ax.set_title(f"{plot_df['data_source'].iloc[0]} source availability")
    ax.set_yticks([])

    for idx, row in enumerate(plot_df.itertuples(index=False)):
        label = "available" if row.load_status == "available" else row.skip_reason
        ax.text(idx, 1.03, label, ha="center", va="bottom", rotation=35, fontsize=8)

    ax.spines[["left", "right", "top"]].set_visible(False)
    fig.tight_layout()
    return save_current_figure(plot_config, name)


def plot_qc_violin(
    adata: ad.AnnData,
    plot_config: PlotConfig,
    keys: Sequence[str] = ("total_counts", "pct_counts_mt"),
    groupby: str = "run_sample_id",
    name: str = "qc_violin_by_sample",
):
    """Plot QC metric violins by sample."""
    sc.pl.violin(
        adata,
        keys=list(keys),
        groupby=groupby,
        rotation=45,
        stripplot=False,
        show=False,
    )
    return save_current_figure(plot_config, name)


def plot_qc_scatter(
    adata: ad.AnnData,
    plot_config: PlotConfig,
    x: str = "total_counts",
    y: str = "pct_counts_mt",
    color: str = "run_sample_id",
    name: str = "qc_scatter_total_counts_pct_counts_mt",
):
    """Plot a QC scatter plot."""
    sc.pl.scatter(
        adata,
        x=x,
        y=y,
        color=color,
        title=f"{x} vs {y}",
        show=False,
    )
    return save_current_figure(plot_config, name)


def plot_embedding(
    adata: ad.AnnData,
    plot_config: PlotConfig,
    basis: str = "umap",
    color: Sequence[str] | str = ("run_sample_id",),
    name: str | None = None,
    **kwargs,
):
    """Plot a Scanpy embedding and save/show it."""
    colors = [color] if isinstance(color, str) else list(color)
    if basis == "umap":
        sc.pl.umap(adata, color=colors, show=False, **kwargs)
    else:
        sc.pl.embedding(adata, basis=basis, color=colors, show=False, **kwargs)
    return save_current_figure(plot_config, name or f"{basis}_{'_'.join(map(_safe_name, colors))}")


def plot_marker_panel(
    adata: ad.AnnData,
    plot_config: PlotConfig,
    markers: Sequence[str],
    basis: str = "umap",
    color_map: str = "magma",
    ncols: int = 3,
    name: str = "marker_panel",
):
    """Plot markers present in the AnnData object."""
    markers_present = [gene for gene in markers if gene in adata.var_names]
    if not markers_present:
        print("No requested markers were present in adata.var_names.")
        return None

    if basis == "umap":
        sc.pl.umap(
            adata,
            color=markers_present,
            color_map=color_map,
            ncols=ncols,
            show=False,
        )
    else:
        sc.pl.embedding(
            adata,
            basis=basis,
            color=markers_present,
            color_map=color_map,
            ncols=ncols,
            show=False,
        )
    return save_current_figure(plot_config, name)


def plot_sample_counts(source_table: pd.DataFrame, plot_config: PlotConfig, name: str = "loaded_vs_skipped_counts"):
    """Plot available versus missing sample counts."""
    summary = (
        source_table.groupby("load_status")
        .size()
        .reindex(["available", "missing_source"], fill_value=0)
        .rename("n_samples")
    )
    colors = ["#2F855A", "#C53030"]
    fig, ax = plt.subplots(figsize=(5.5, 3.8))
    ax.bar(summary.index, summary.values, color=colors)
    ax.set_ylabel("n samples")
    ax.set_title("Requested sample source status")
    ax.spines[["right", "top"]].set_visible(False)
    for idx, value in enumerate(summary.values):
        ax.text(idx, value + 0.05, str(int(value)), ha="center", va="bottom")
    fig.tight_layout()
    return save_current_figure(plot_config, name)
