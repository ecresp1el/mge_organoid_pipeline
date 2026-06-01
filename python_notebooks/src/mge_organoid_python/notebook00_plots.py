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


def plot_highest_expr_genes(
    adata: ad.AnnData,
    plot_config: PlotConfig,
    n_top: int = 20,
    name: str = "highest_expr_genes_top20",
):
    """Plot the highest expressed genes with Scanpy's helper."""
    sc.pl.highest_expr_genes(adata, n_top=n_top, show=False)
    return save_current_figure(plot_config, name)


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


def plot_manual_ec_qc_violin(
    adata: ad.AnnData,
    plot_config: PlotConfig,
    keys: Sequence[str] = ("n_genes_by_counts", "total_counts", "pct_counts_mt"),
    name: str = "manual_ec_qc_violin",
):
    """Plot the requested manual_ec QC metric violins."""
    sc.pl.violin(
        adata,
        keys=list(keys),
        jitter=0.4,
        multi_panel=True,
        show=False,
    )
    return save_current_figure(plot_config, name)


def plot_manual_ec_qc_scatter(
    adata: ad.AnnData,
    plot_config: PlotConfig,
    x: str = "total_counts",
    y: str = "pct_counts_mt",
    name: str | None = None,
):
    """Plot one requested manual_ec QC scatter plot."""
    sc.pl.scatter(
        adata,
        x=x,
        y=y,
        title=f"{x} vs {y}",
        show=False,
    )
    return save_current_figure(plot_config, name or f"manual_ec_scatter_{x}_{y}")


def _sample_grid(
    adata_names: Sequence[str],
    adata_list: Sequence[ad.AnnData],
    width_per_sample: float,
    height: float,
    sharex: bool = False,
    sharey: bool = False,
):
    """Create a 1xN subplot grid for per-sample Scanpy plots."""
    names = [str(name) for name in adata_names]
    sample_adatas = list(adata_list)
    if len(names) != len(sample_adatas):
        raise ValueError(
            f"adata_names has {len(names)} entries but adata_list has {len(sample_adatas)} entries."
        )
    if not sample_adatas:
        return names, sample_adatas, None, []

    fig_width = max(6.5, width_per_sample * len(sample_adatas))
    fig, axes = plt.subplots(
        1,
        len(sample_adatas),
        figsize=(fig_width, height),
        squeeze=False,
        sharex=sharex,
        sharey=sharey,
    )
    return names, sample_adatas, fig, axes.ravel().tolist()


def plot_manual_ec_per_sample_highest_expr_genes(
    adata_names: Sequence[str],
    adata_list: Sequence[ad.AnnData],
    plot_config: PlotConfig,
    n_top: int = 20,
    name: str = "manual_ec_per_sample_highest_expr_genes_top20",
):
    """Plot Scanpy highest-expressed genes as a 1xN per-sample grid."""
    sample_names, sample_adatas, fig, axes = _sample_grid(
        adata_names,
        adata_list,
        width_per_sample=4.2,
        height=4.8,
    )
    if fig is None:
        return None

    for sample_name, sample_adata, ax in zip(sample_names, sample_adatas, axes):
        sc.pl.highest_expr_genes(sample_adata, n_top=n_top, show=False, ax=ax)
        ax.set_title(sample_name, fontsize=10)

    fig.tight_layout()
    return save_current_figure(plot_config, name)


def plot_manual_ec_per_sample_qc_violins(
    adata_names: Sequence[str],
    adata_list: Sequence[ad.AnnData],
    plot_config: PlotConfig,
    keys: Sequence[str] = ("n_genes_by_counts", "total_counts", "pct_counts_mt"),
    name_prefix: str = "manual_ec_per_sample_qc_violin",
):
    """Plot each requested QC violin metric as a 1xN per-sample grid."""
    saved_paths = {}
    for key in keys:
        sample_names, sample_adatas, fig, axes = _sample_grid(
            adata_names,
            adata_list,
            width_per_sample=3.4,
            height=4.2,
            sharey=True,
        )
        if fig is None:
            saved_paths[str(key)] = None
            continue

        for sample_name, sample_adata, ax in zip(sample_names, sample_adatas, axes):
            sc.pl.violin(
                sample_adata,
                keys=[key],
                jitter=0.4,
                show=False,
                ax=ax,
            )
            ax.set_title(sample_name, fontsize=10)
            ax.set_xlabel("")

        fig.tight_layout()
        saved_paths[str(key)] = save_current_figure(plot_config, f"{name_prefix}_{key}")

    return saved_paths


def plot_manual_ec_per_sample_qc_scatter(
    adata_names: Sequence[str],
    adata_list: Sequence[ad.AnnData],
    plot_config: PlotConfig,
    x: str = "total_counts",
    y: str = "pct_counts_mt",
    name: str | None = None,
):
    """Plot one Scanpy QC scatter as a 1xN per-sample grid."""
    sample_names, sample_adatas, fig, axes = _sample_grid(
        adata_names,
        adata_list,
        width_per_sample=3.8,
        height=4.0,
        sharex=True,
        sharey=True,
    )
    if fig is None:
        return None

    for sample_name, sample_adata, ax in zip(sample_names, sample_adatas, axes):
        sc.pl.scatter(
            sample_adata,
            x=x,
            y=y,
            title=sample_name,
            show=False,
            ax=ax,
        )

    fig.tight_layout()
    return save_current_figure(plot_config, name or f"manual_ec_per_sample_scatter_{x}_{y}")


def plot_manual_ec_source_comparison_cell_counts(
    comparison_df: pd.DataFrame,
    plot_config: PlotConfig,
    name: str = "manual_ec_source_comparison_cell_counts",
):
    """Plot retained manual_ec cell counts by sample and source."""
    required = {"run_label", "run_sample_id", "manual_ec_retained_n_cells"}
    if comparison_df.empty or required.difference(comparison_df.columns):
        return None

    pivot_df = comparison_df.pivot_table(
        index="run_sample_id",
        columns="run_label",
        values="manual_ec_retained_n_cells",
        aggfunc="sum",
        fill_value=0,
    )
    ax = pivot_df.plot(kind="bar", figsize=(max(7, 1.2 * len(pivot_df)), 4.5))
    ax.set_xlabel("run_sample_id")
    ax.set_ylabel("manual_ec retained cells")
    ax.set_title("manual_ec retained cells by source")
    ax.legend(title="run_label", fontsize=8)
    plt.tight_layout()
    return save_current_figure(plot_config, name)


def plot_manual_ec_source_comparison_qc_metrics(
    qc_df: pd.DataFrame,
    plot_config: PlotConfig,
    metric: str = "median_n_genes_by_counts",
    name: str = "manual_ec_source_comparison_qc_metrics",
):
    """Plot one manual_ec QC metric by sample and source."""
    required = {"run_label", "run_sample_id", metric}
    if qc_df.empty or required.difference(qc_df.columns):
        return None

    pivot_df = qc_df.pivot_table(
        index="run_sample_id",
        columns="run_label",
        values=metric,
        aggfunc="median",
    )
    ax = pivot_df.plot(kind="bar", figsize=(max(7, 1.2 * len(pivot_df)), 4.5))
    ax.set_xlabel("run_sample_id")
    ax.set_ylabel(metric)
    ax.set_title(f"manual_ec {metric} by source")
    ax.legend(title="run_label", fontsize=8)
    plt.tight_layout()
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
