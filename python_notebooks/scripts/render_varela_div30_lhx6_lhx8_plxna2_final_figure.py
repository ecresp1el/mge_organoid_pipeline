#!/usr/bin/env python3
"""Render Varela DIV30 cluster-context and neuronal LHX6/LHX8/PLXNA2 figures.

The primary neuronal analysis is restricted to original Seurat clusters 1, 2,
and 4. Radial-glia clusters 0, 3, and 7 and inhibitory-progenitor cluster 6
remain visible only in the full-dataset context row. No high/low populations or
pooled-cell significance tests are generated.
"""

from __future__ import annotations

import argparse
from datetime import datetime
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
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
from scipy import stats

from mge_organoid_python.cross_study_marker_expression import (
    BACKGROUND_POINT_COLOR,
    EXPRESSION_COLOR_FLOOR,
    PDF_EXPORT_DPI,
    PNG_EXPORT_DPI,
    SVG_EXPORT_DPI,
    _apply_internal_umap_plot_filters,
    default_cross_study_marker_specs,
    extract_marker_expression_from_h5ad,
)


PROJECT_ROOT_DEFAULT = Path("/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder")
FINAL_FOLDER_DEFAULT = "fig_varela_div30_lhx6_lhx8_plxna2_v1_candidate"
PRIMARY_STEM = "varela_div30_lhx6_lhx8_populations_plxna2"
SPLIT_FEATURE_STEM = "varela_div30_neuronal_features_split_by_organoid"
SPLIT_RELATIONSHIP_STEM = "varela_div30_plxna2_vs_lhx6_minus_lhx8_by_cluster_organoid"
SPLIT_RATIO_RELATIONSHIP_STEM = "varela_div30_plxna2_vs_lhx6_lhx8_log2_ratio_by_cluster_organoid"

GENES = ["LHX6", "LHX8", "PLXNA2"]
NEURONAL_CLUSTERS = ["1", "4", "2"]
EXCLUDED_PRIMARY_CLUSTERS = ["0", "3", "6", "7"]
SCORE_COL = "LHX6_minus_LHX8_zscore"
SCORE_LABEL = "z(LHX6) − z(LHX8)"
RATIO_COL = "LHX6_LHX8_log2_ratio"
RATIO_LABEL = "log2((LHX6 + ε)/(LHX8 + ε))"
CLUSTER_MAPPING_RELATIVE_PATH = Path(
    "results/div30_paper_cluster_annotations/div30_paper_cluster_annotations_v1/"
    "tables/div30_paper_cluster_annotation_mapping.tsv"
)

PUBLICATION_CLUSTER_BY_SEURAT = {
    "0": "1",
    "3": "1",
    "7": "1",
    "6": "2",
    "1": "3",
    "4": "4",
    "2": "5",
}
PUBLICATION_CLUSTER_LABELS = {
    "1": "Radial Glia",
    "2": "Inhibitory Progenitors",
    "3": "SST+ cIN",
    "4": "PV Neuron Precursor",
    "5": "MGE Subpallial Neuron",
}
PUBLICATION_CLUSTER_COLORS = {
    "1": "#F8766D",
    "2": "#A3A500",
    "3": "#00BF7D",
    "4": "#00B0F6",
    "5": "#C77CFF",
}
# Seurat cluster colors inherit the corresponding publication-cluster color.
CLUSTER_COLORS = {
    cluster: PUBLICATION_CLUSTER_COLORS[publication_cluster]
    for cluster, publication_cluster in PUBLICATION_CLUSTER_BY_SEURAT.items()
}
TEXT_COLOR = "#202020"
LINEAGE_GROUPS = {
    "1": "cortical_interneuron_lineage",
    "4": "cortical_interneuron_lineage",
    "2": "MGE_subpallial_neurons",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT_DEFAULT)
    parser.add_argument("--final-dir", type=Path, default=None)
    parser.add_argument("--random-state", type=int, default=1729)
    parser.add_argument(
        "--ratio-epsilon-cp10k",
        type=float,
        default=1.0,
        help="Pseudocount added on the back-transformed linear CP10K scale for the LHX6/LHX8 log2 ratio.",
    )
    parser.add_argument("--fresh", action="store_true")
    return parser.parse_args()


def div30_spec(project_root: Path):
    matches = [
        spec
        for spec in default_cross_study_marker_specs(project_root)
        if spec.study_id == "varela_div30"
    ]
    if len(matches) != 1:
        raise ValueError(f"Expected one varela_div30 spec, found {len(matches)}.")
    return matches[0]


def cluster_mapping_path(project_root: Path) -> Path:
    return project_root / CLUSTER_MAPPING_RELATIVE_PATH


def add_cluster_annotations(
    data: pd.DataFrame, project_root: Path
) -> tuple[pd.DataFrame, pd.DataFrame]:
    path = cluster_mapping_path(project_root)
    if not path.exists():
        raise FileNotFoundError(f"Missing canonical DIV30 cluster mapping: {path}")
    mapping = pd.read_csv(path, sep="\t", dtype={"seurat_clusters": str})
    required = {"seurat_clusters", "paper_cluster_annotation"}
    missing = sorted(required.difference(mapping.columns))
    if missing:
        raise ValueError(f"DIV30 cluster mapping is missing columns: {missing}")
    if mapping["seurat_clusters"].duplicated().any():
        raise ValueError("DIV30 cluster mapping has duplicate seurat_clusters values.")

    out = data.copy()
    out["cluster"] = out["cluster"].astype(str)
    name_map = mapping.set_index("seurat_clusters")["paper_cluster_annotation"]
    out["paper_cluster_annotation"] = out["cluster"].map(name_map)
    if out["paper_cluster_annotation"].isna().any():
        unmapped = sorted(out.loc[out["paper_cluster_annotation"].isna(), "cluster"].unique())
        raise ValueError(f"Unmapped Varela DIV30 clusters: {unmapped}")
    out["cluster_name_label"] = out["cluster"] + " · " + out["paper_cluster_annotation"].astype(str)
    out["publication_cluster"] = out["cluster"].map(PUBLICATION_CLUSTER_BY_SEURAT)
    out["publication_cluster_annotation"] = out["publication_cluster"].map(
        PUBLICATION_CLUSTER_LABELS
    )
    out["publication_cluster_label"] = (
        out["publication_cluster"] + " · " + out["publication_cluster_annotation"]
    )
    if out["publication_cluster"].isna().any():
        unmapped = sorted(out.loc[out["publication_cluster"].isna(), "cluster"].unique())
        raise ValueError(f"Unmapped publication clusters for Seurat clusters: {unmapped}")
    return out, mapping


def add_neuronal_score(
    data: pd.DataFrame,
    ratio_epsilon_cp10k: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if ratio_epsilon_cp10k <= 0:
        raise ValueError("ratio_epsilon_cp10k must be greater than zero.")
    subset = data.loc[data["cluster"].isin(NEURONAL_CLUSTERS)].copy()
    if subset.empty:
        raise ValueError("No cells remain in neuronal clusters 1, 2, and 4.")
    rows: list[dict[str, object]] = []
    for gene in ["LHX6", "LHX8"]:
        values = pd.to_numeric(subset[gene], errors="coerce").to_numpy(dtype=float)
        mean = float(np.nanmean(values))
        sd = float(np.nanstd(values, ddof=0))
        if not np.isfinite(sd) or sd <= 0:
            raise ValueError(f"Cannot z-score {gene}: SD={sd}")
        subset[f"{gene}_zscore"] = (values - mean) / sd
        rows.append(
            {
                "gene": gene,
                "reference_cells": "Varela DIV30 original Seurat clusters 1, 2, and 4",
                "n_cells": int(np.isfinite(values).sum()),
                "mean_log1p_cp10k": mean,
                "population_sd_log1p_cp10k": sd,
            }
        )
    subset[SCORE_COL] = subset["LHX6_zscore"] - subset["LHX8_zscore"]
    lhx6_cp10k = np.expm1(subset["LHX6"].to_numpy(dtype=float))
    lhx8_cp10k = np.expm1(subset["LHX8"].to_numpy(dtype=float))
    subset[RATIO_COL] = np.log2(
        (lhx6_cp10k + ratio_epsilon_cp10k) / (lhx8_cp10k + ratio_epsilon_cp10k)
    )
    subset["lineage_comparison_group"] = subset["cluster"].map(LINEAGE_GROUPS)
    ratio_parameters = pd.DataFrame(
        [
            {
                "method": "stabilized_LHX6_LHX8_log2_ratio",
                "formula": "log2((expm1(LHX6_log1p_cp10k)+epsilon)/(expm1(LHX8_log1p_cp10k)+epsilon))",
                "input_scale": "linear CP10K back-transformed from log1p(CP10K)",
                "epsilon_cp10k": ratio_epsilon_cp10k,
                "reference_cells": "Varela DIV30 original Seurat clusters 1, 2, and 4",
                "n_cells": int(subset.shape[0]),
            }
        ]
    )
    return subset, pd.DataFrame(rows), ratio_parameters


def cluster_annotation_summary(data: pd.DataFrame) -> pd.DataFrame:
    summary = (
        data.groupby(
            [
                "cluster",
                "paper_cluster_annotation",
                "cluster_name_label",
                "publication_cluster",
                "publication_cluster_annotation",
                "publication_cluster_label",
            ],
            observed=True,
        )
        .size()
        .reset_index(name="n_cells")
    )
    summary["included_in_primary_neuronal_analysis"] = summary["cluster"].isin(NEURONAL_CLUSTERS)
    summary["primary_analysis_role"] = np.where(
        summary["cluster"].isin(["1", "4"]),
        "cortical_interneuron_lineage",
        np.where(summary["cluster"].eq("2"), "MGE_subpallial_neurons", "excluded_developmental_state"),
    )
    return summary.sort_values("cluster", key=lambda values: pd.to_numeric(values, errors="coerce")).reset_index(drop=True)


def publication_cluster_annotation_summary(data: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for publication_cluster in ["1", "2", "3", "4", "5"]:
        subset = data.loc[data["publication_cluster"].eq(publication_cluster)]
        rows.append(
            {
                "publication_cluster": publication_cluster,
                "publication_cluster_annotation": PUBLICATION_CLUSTER_LABELS[publication_cluster],
                "source_seurat_clusters": ",".join(
                    sorted(subset["cluster"].unique(), key=int)
                ),
                "n_cells": int(len(subset)),
                "color": PUBLICATION_CLUSTER_COLORS[publication_cluster],
            }
        )
    return pd.DataFrame(rows)


def expression_summaries(
    neuronal: pd.DataFrame,
    positive_threshold: float = 0.5,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    cluster_rows: list[dict[str, object]] = []
    organoid_rows: list[dict[str, object]] = []
    for cluster in NEURONAL_CLUSTERS:
        cluster_data = neuronal.loc[neuronal["cluster"].eq(cluster)]
        annotation = str(cluster_data["paper_cluster_annotation"].iloc[0])
        lineage_group = LINEAGE_GROUPS[cluster]
        for gene in GENES:
            values = cluster_data[gene].to_numpy(dtype=float)
            cluster_rows.append(
                summarize_values(
                    values,
                    cluster=cluster,
                    annotation=annotation,
                    lineage_group=lineage_group,
                    gene=gene,
                    sample="all",
                    positive_threshold=positive_threshold,
                )
            )
        for sample, sample_data in cluster_data.groupby("sample", sort=True, observed=True):
            for gene in GENES:
                organoid_rows.append(
                    summarize_values(
                        sample_data[gene].to_numpy(dtype=float),
                        cluster=cluster,
                        annotation=annotation,
                        lineage_group=lineage_group,
                        gene=gene,
                        sample=str(sample),
                        positive_threshold=positive_threshold,
                    )
                )
    return pd.DataFrame(cluster_rows), pd.DataFrame(organoid_rows)


def lineage_group_expression_summaries(
    neuronal: pd.DataFrame,
    positive_threshold: float = 0.5,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    pooled_rows: list[dict[str, object]] = []
    organoid_rows: list[dict[str, object]] = []
    group_meta = {
        "cortical_interneuron_lineage": ("1+4", "SST+ cIN and PV neuron precursor"),
        "MGE_subpallial_neurons": ("2", "MGE subpallial neurons"),
    }
    for lineage_group, (cluster_label, annotation) in group_meta.items():
        group_data = neuronal.loc[neuronal["lineage_comparison_group"].eq(lineage_group)]
        for gene in GENES:
            pooled_rows.append(
                summarize_values(
                    group_data[gene].to_numpy(dtype=float),
                    cluster=cluster_label,
                    annotation=annotation,
                    lineage_group=lineage_group,
                    gene=gene,
                    sample="all",
                    positive_threshold=positive_threshold,
                )
            )
        for sample, sample_data in group_data.groupby("sample", sort=True, observed=True):
            for gene in GENES:
                organoid_rows.append(
                    summarize_values(
                        sample_data[gene].to_numpy(dtype=float),
                        cluster=cluster_label,
                        annotation=annotation,
                        lineage_group=lineage_group,
                        gene=gene,
                        sample=str(sample),
                        positive_threshold=positive_threshold,
                    )
                )
    return pd.DataFrame(pooled_rows), pd.DataFrame(organoid_rows)


def summarize_values(
    values: np.ndarray,
    *,
    cluster: str,
    annotation: str,
    lineage_group: str,
    gene: str,
    sample: str,
    positive_threshold: float,
) -> dict[str, object]:
    values = values[np.isfinite(values)]
    publication_cluster = (
        "3+4" if cluster == "1+4" else PUBLICATION_CLUSTER_BY_SEURAT[cluster]
    )
    publication_annotation = (
        "Cortical Interneuron Lineage"
        if cluster == "1+4"
        else PUBLICATION_CLUSTER_LABELS[publication_cluster]
    )
    return {
        "sample": sample,
        "sample_display": "all" if sample == "all" else sample.replace("9583-", ""),
        "cluster": cluster,
        "paper_cluster_annotation": annotation,
        "publication_cluster": publication_cluster,
        "publication_cluster_annotation": publication_annotation,
        "lineage_comparison_group": lineage_group,
        "gene": gene,
        "n_cells": int(values.size),
        "mean_log1p_cp10k": float(np.mean(values)),
        "median_log1p_cp10k": float(np.median(values)),
        "q25_log1p_cp10k": float(np.quantile(values, 0.25)),
        "q75_log1p_cp10k": float(np.quantile(values, 0.75)),
        "positive_threshold_log1p_cp10k": positive_threshold,
        "n_positive": int((values >= positive_threshold).sum()),
        "pct_positive": float((values >= positive_threshold).mean() * 100.0),
    }


def relationship_summaries(
    neuronal: pd.DataFrame,
    predictor: str,
    predictor_definition: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    pooled_rows: list[dict[str, object]] = []
    organoid_rows: list[dict[str, object]] = []
    for cluster in NEURONAL_CLUSTERS:
        cluster_data = neuronal.loc[neuronal["cluster"].eq(cluster)]
        pooled_rows.append(summarize_relationship(cluster_data, cluster, "all", predictor, predictor_definition))
        for sample, sample_data in cluster_data.groupby("sample", sort=True, observed=True):
            organoid_rows.append(
                summarize_relationship(sample_data, cluster, str(sample), predictor, predictor_definition)
            )
    return pd.DataFrame(pooled_rows), pd.DataFrame(organoid_rows)


def summarize_relationship(
    data: pd.DataFrame,
    cluster: str,
    sample: str,
    predictor: str,
    predictor_definition: str,
) -> dict[str, object]:
    x = data[predictor].to_numpy(dtype=float)
    y = data["PLXNA2"].to_numpy(dtype=float)
    finite = np.isfinite(x) & np.isfinite(y)
    x = x[finite]
    y = y[finite]
    pearson_r = float(np.corrcoef(x, y)[0, 1]) if x.size >= 3 and np.std(x) > 0 and np.std(y) > 0 else np.nan
    spearman_rho = float(stats.spearmanr(x, y).statistic) if x.size >= 3 else np.nan
    if x.size >= 2 and np.std(x) > 0:
        slope, intercept = np.polyfit(x, y, 1)
    else:
        slope, intercept = np.nan, np.nan
    annotation = str(data["paper_cluster_annotation"].iloc[0])
    return {
        "sample": sample,
        "sample_display": "all" if sample == "all" else sample.replace("9583-", ""),
        "cluster": cluster,
        "paper_cluster_annotation": annotation,
        "publication_cluster": PUBLICATION_CLUSTER_BY_SEURAT[cluster],
        "publication_cluster_annotation": PUBLICATION_CLUSTER_LABELS[
            PUBLICATION_CLUSTER_BY_SEURAT[cluster]
        ],
        "lineage_comparison_group": LINEAGE_GROUPS[cluster],
        "n_cells": int(x.size),
        "predictor": predictor,
        "predictor_definition": predictor_definition,
        "pearson_r_descriptive": pearson_r,
        "spearman_rho_descriptive": spearman_rho,
        "linear_slope_descriptive": float(slope),
        "linear_intercept_descriptive": float(intercept),
        "p_values_computed": False,
    }


def expression_floor_colormap(vmax: float) -> LinearSegmentedColormap:
    floor_fraction = float(np.clip(EXPRESSION_COLOR_FLOOR / max(vmax, EXPRESSION_COLOR_FLOOR + 1e-6), 0, 1))
    return LinearSegmentedColormap.from_list(
        "grayFloorBlue",
        [(0.0, BACKGROUND_POINT_COLOR), (floor_fraction, BACKGROUND_POINT_COLOR), (1.0, "#0000FF")],
    )


def umap_bounds(data: pd.DataFrame) -> tuple[tuple[float, float], tuple[float, float]]:
    x = data["umap_1"].to_numpy(dtype=float)
    y = data["umap_2"].to_numpy(dtype=float)
    xlo, xhi = np.nanquantile(x, [0.001, 0.999])
    ylo, yhi = np.nanquantile(y, [0.001, 0.999])
    xpad = (xhi - xlo) * 0.04
    ypad = (yhi - ylo) * 0.04
    return (float(xlo - xpad), float(xhi + xpad)), (float(ylo - ypad), float(yhi + ypad))


def style_umap_axis(
    ax: plt.Axes,
    bounds: tuple[tuple[float, float], tuple[float, float]],
) -> None:
    ax.set_xlim(*bounds[0])
    ax.set_ylim(*bounds[1])
    ax.set_aspect("equal", adjustable="box")
    ax.set_axis_off()


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.18, 1.09, label, transform=ax.transAxes, fontsize=9.5, fontweight="bold", va="top")


def style_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", width=0.6, length=2.5, colors=TEXT_COLOR, labelsize=6.2)


def make_umap_and_key_axes(cell, figure: plt.Figure, key_height: float = 1.2):
    subgrid = cell.subgridspec(2, 1, height_ratios=[5.0, key_height], hspace=0.04)
    return figure.add_subplot(subgrid[0, 0]), figure.add_subplot(subgrid[1, 0])


def make_umap_and_colorbar_axes(cell, figure: plt.Figure, bar_height: float = 0.46):
    subgrid = cell.subgridspec(
        2,
        3,
        height_ratios=[5.0, bar_height * 0.5],
        width_ratios=[1.0, 2.0, 1.0],
        hspace=0.18,
        wspace=0.0,
    )
    return figure.add_subplot(subgrid[0, :]), figure.add_subplot(subgrid[1, 1])


def plot_cluster_key_umap(
    ax: plt.Axes,
    key_ax: plt.Axes,
    data: pd.DataFrame,
    clusters: list[str],
    bounds: tuple[tuple[float, float], tuple[float, float]],
    title: str,
) -> None:
    for cluster in clusters:
        subset = data.loc[data["cluster"].eq(cluster)]
        ax.scatter(
            subset["umap_1"], subset["umap_2"], s=0.55, color=CLUSTER_COLORS[cluster],
            alpha=0.78, linewidths=0, rasterized=True,
        )
    style_umap_axis(ax, bounds)
    ax.set_title(title, fontsize=7.1, fontweight="bold", pad=3)
    key_ax.set_axis_off()
    ncols = 2 if len(clusters) > 4 else 1
    nrows = int(np.ceil(len(clusters) / ncols))
    for idx, cluster in enumerate(clusters):
        col = idx // nrows
        row = idx % nrows
        subset = data.loc[data["cluster"].eq(cluster)]
        annotation = str(subset["paper_cluster_annotation"].iloc[0])
        x = 0.02 + col * 0.62
        y = 0.91 - row * (0.82 / max(nrows - 1, 1))
        key_ax.text(
            x, y, f"{cluster}  {annotation}", transform=key_ax.transAxes,
            color=CLUSTER_COLORS[cluster], fontsize=3.9, fontweight="bold", ha="left", va="center",
        )


def plot_publication_cluster_key_umap(
    ax: plt.Axes,
    key_ax: plt.Axes,
    data: pd.DataFrame,
    publication_clusters: list[str],
    bounds: tuple[tuple[float, float], tuple[float, float]],
    title: str,
    key_fontsize: float = 6.3,
) -> None:
    """Plot the five publication populations, combining Seurat 0/3/7 as radial glia."""
    for publication_cluster in publication_clusters:
        subset = data.loc[data["publication_cluster"].eq(publication_cluster)]
        ax.scatter(
            subset["umap_1"],
            subset["umap_2"],
            s=0.55,
            color=PUBLICATION_CLUSTER_COLORS[publication_cluster],
            alpha=0.78,
            linewidths=0,
            rasterized=True,
        )
    style_umap_axis(ax, bounds)
    ax.set_title(title, fontsize=8.2, fontweight="bold", pad=3)
    key_ax.set_axis_off()
    # Keep the enlarged key legible by placing every label on its own line.
    ncols = 1
    nrows = len(publication_clusters)
    for idx, publication_cluster in enumerate(publication_clusters):
        col = idx // nrows
        row = idx % nrows
        x = 0.02
        y = 0.93 - row * (0.86 / max(nrows - 1, 1))
        key_ax.text(
            x,
            y,
            f"{publication_cluster}  {PUBLICATION_CLUSTER_LABELS[publication_cluster]}",
            transform=key_ax.transAxes,
            color=PUBLICATION_CLUSTER_COLORS[publication_cluster],
            fontsize=key_fontsize,
            fontweight="bold",
            ha="left",
            va="center",
        )


def plot_feature_umap(
    fig: plt.Figure,
    ax: plt.Axes,
    cax: plt.Axes,
    data: pd.DataFrame,
    feature: str,
    bounds: tuple[tuple[float, float], tuple[float, float]],
    expression_vmax: float,
    score_vmax: float,
    ratio_vmax: float,
    title: str | None = None,
) -> None:
    if feature == SCORE_COL:
        cmap = "coolwarm"
        vmin, vmax = -score_vmax, score_vmax
        label = SCORE_LABEL
    elif feature == RATIO_COL:
        cmap = "coolwarm"
        vmin, vmax = -ratio_vmax, ratio_vmax
        label = RATIO_LABEL
    else:
        cmap = expression_floor_colormap(expression_vmax)
        vmin, vmax = 0.0, expression_vmax
        label = "log1p(CP10K)"
    plot_data = data.sort_values(feature, kind="stable")
    points = ax.scatter(
        plot_data["umap_1"], plot_data["umap_2"], c=plot_data[feature], s=0.55,
        cmap=cmap, vmin=vmin, vmax=vmax, linewidths=0, rasterized=True,
    )
    style_umap_axis(ax, bounds)
    if title is None:
        title = SCORE_LABEL if feature == SCORE_COL else RATIO_LABEL if feature == RATIO_COL else rf"$\it{{{feature}}}$"
    ax.set_title(title, fontsize=7.1, pad=3)
    ticks = [vmin, 0, vmax] if feature in {SCORE_COL, RATIO_COL} else [0, vmax]
    colorbar = fig.colorbar(points, cax=cax, orientation="horizontal", ticks=ticks)
    colorbar.set_ticklabels([format_colorbar_tick(value) for value in ticks])
    colorbar.ax.tick_params(labelsize=3.8, length=1.1, width=0.35, pad=0.7)
    colorbar.outline.set_linewidth(0.35)
    colorbar.set_label(label, fontsize=3.8, labelpad=0.2)


def round_scale_limit(value: float, increment: float = 0.5) -> float:
    """Round a positive display limit upward to a compact publication value."""
    return float(np.ceil(value / increment) * increment)


def format_colorbar_tick(value: float) -> str:
    """Show zero without decimals and other colorbar ticks with at most one decimal."""
    if np.isclose(value, 0.0):
        return "0"
    rounded = round(float(value), 1)
    return f"{rounded:g}"


def common_limits(neuronal: pd.DataFrame) -> tuple[float, float, float, float]:
    per_gene_q99 = [float(np.nanquantile(neuronal[gene], 0.99)) for gene in GENES]
    expression_vmax = round_scale_limit(max(EXPRESSION_COLOR_FLOOR, max(per_gene_q99)))
    score_vmax = round_scale_limit(
        max(1.0, float(np.nanquantile(np.abs(neuronal[SCORE_COL]), 0.99)))
    )
    ratio_vmax = round_scale_limit(
        max(1.0, float(np.nanquantile(np.abs(neuronal[RATIO_COL]), 0.99)))
    )
    expression_ymax = max(2.0, max(float(np.nanquantile(neuronal[gene], 0.995)) for gene in GENES) * 1.04)
    return expression_vmax, score_vmax, ratio_vmax, expression_ymax


def plot_cluster_expression_panel(
    ax: plt.Axes,
    neuronal: pd.DataFrame,
    organoid_summary: pd.DataFrame,
    gene: str,
    y_max: float,
    random_state: int,
) -> None:
    values = [neuronal.loc[neuronal["cluster"].eq(cluster), gene].to_numpy(dtype=float) for cluster in NEURONAL_CLUSTERS]
    rng = np.random.default_rng(random_state)
    draw_values = [rng.choice(value, size=min(value.size, 10_000), replace=False) for value in values]
    violins = ax.violinplot(draw_values, positions=np.arange(3), widths=0.72, showextrema=False)
    for body, cluster in zip(violins["bodies"], NEURONAL_CLUSTERS, strict=True):
        body.set_facecolor(CLUSTER_COLORS[cluster])
        body.set_edgecolor("none")
        body.set_alpha(0.58)
    for xpos, cluster, value in zip(np.arange(3), NEURONAL_CLUSTERS, values, strict=True):
        q25, median, q75 = np.quantile(value, [0.25, 0.5, 0.75])
        ax.plot([xpos, xpos], [q25, q75], color="white", linewidth=2.5, solid_capstyle="round", zorder=3)
        ax.scatter([xpos], [median], s=13, color=TEXT_COLOR, edgecolor="white", linewidth=0.35, zorder=4)
        sample_means = organoid_summary.loc[
            organoid_summary["cluster"].eq(cluster) & organoid_summary["gene"].eq(gene),
            "mean_log1p_cp10k",
        ].to_numpy(dtype=float)
        jitter = np.linspace(-0.13, 0.13, sample_means.size)
        ax.scatter(
            xpos + jitter, sample_means, s=12, facecolor="white", edgecolor=TEXT_COLOR,
            linewidth=0.45, alpha=0.95, zorder=5,
        )
    ax.axvline(1.5, color="#888888", linewidth=0.55, linestyle=(0, (2, 2)))
    ax.set_xticks(np.arange(3))
    ax.set_xticklabels(
        ["3\nSST+ cIN", "4\nPV\nprecursor", "5\nMGE\nsubpallial"],
        linespacing=0.95,
    )
    ax.set_ylim(0, y_max)
    ax.set_title(rf"$\it{{{gene}}}$", fontsize=7.4, fontweight="bold", pad=5)
    ax.set_ylabel("Expression, log1p(CP10K)", fontsize=6.3)
    style_axis(ax)


def deterministic_sample(data: pd.DataFrame, n: int, random_state: int) -> pd.DataFrame:
    if data.shape[0] <= n:
        return data
    return data.sample(n=n, random_state=random_state)


def plot_relationship_panel(
    ax: plt.Axes,
    data: pd.DataFrame,
    cluster: str,
    predictor: str,
    predictor_label: str,
    predictor_limit: float,
    plx_limit: float,
    random_state: int,
) -> None:
    cluster_data = data.loc[data["cluster"].eq(cluster)]
    draw = deterministic_sample(cluster_data, 7_000, random_state)
    ax.scatter(
        draw[predictor], draw["PLXNA2"], s=2.0, color=CLUSTER_COLORS[cluster],
        alpha=0.12, linewidths=0, rasterized=True,
    )
    x = cluster_data[predictor].to_numpy(dtype=float)
    y = cluster_data["PLXNA2"].to_numpy(dtype=float)
    if x.size >= 2 and np.std(x) > 0:
        slope, intercept = np.polyfit(x, y, 1)
        line_x = np.linspace(-predictor_limit, predictor_limit, 80)
        ax.plot(line_x, slope * line_x + intercept, color=CLUSTER_COLORS[cluster], linewidth=1.0)
    publication_cluster = PUBLICATION_CLUSTER_BY_SEURAT[cluster]
    annotation = PUBLICATION_CLUSTER_LABELS[publication_cluster]
    ax.set_title(
        f"{publication_cluster}  {annotation}",
        fontsize=6.8,
        fontweight="bold",
        color=PUBLICATION_CLUSTER_COLORS[publication_cluster],
        pad=5,
    )
    ax.set_xlim(-predictor_limit, predictor_limit)
    ax.set_ylim(0, plx_limit)
    ax.set_xlabel(predictor_label, fontsize=6.2)
    ax.set_ylabel(r"$\it{PLXNA2}$, log1p(CP10K)", fontsize=6.2)
    style_axis(ax)


def save_figure(fig: plt.Figure, final_dir: Path, stem: str) -> list[Path]:
    outputs: list[Path] = []
    for ext, dpi in (("png", PNG_EXPORT_DPI), ("pdf", PDF_EXPORT_DPI), ("svg", SVG_EXPORT_DPI)):
        path = final_dir / "figures" / ext / f"{stem}.{ext}"
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
        outputs.append(path)
    plt.close(fig)
    return outputs


def plot_primary_figure(
    all_data: pd.DataFrame,
    neuronal: pd.DataFrame,
    organoid_summary: pd.DataFrame,
    final_dir: Path,
    random_state: int,
) -> list[Path]:
    expression_vmax, score_vmax, ratio_vmax, expression_ymax = common_limits(neuronal)
    full_bounds = umap_bounds(all_data)
    neuronal_bounds = umap_bounds(neuronal)
    fig = plt.figure(figsize=(7.2, 6.65), facecolor="white")
    outer = fig.add_gridspec(
        3, 1, left=0.07, right=0.985, bottom=0.07, top=0.97,
        height_ratios=[1.0, 1.0, 1.14], hspace=0.55,
    )

    # A-D: full-dataset context.
    top = outer[0].subgridspec(1, 4, width_ratios=[1.55, 1, 1, 1], wspace=0.30)
    ax_a, key_a = make_umap_and_key_axes(top[0, 0], fig, key_height=1.9)
    plot_publication_cluster_key_umap(
        ax_a,
        key_a,
        all_data,
        ["1", "2", "3", "4", "5"],
        full_bounds,
        "DIV30 MGEOs",
        key_fontsize=6.3,
    )
    panel_label(ax_a, "a")
    for idx, (gene, label) in enumerate(zip(GENES, ["b", "c", "d"], strict=True), start=1):
        ax, cax = make_umap_and_colorbar_axes(top[0, idx], fig, bar_height=0.52)
        plot_feature_umap(fig, ax, cax, all_data, gene, full_bounds, expression_vmax, score_vmax, ratio_vmax)
        panel_label(ax, label)

    # E-H: neuronal-only pooled UMAPs.
    middle = outer[1].subgridspec(1, 4, width_ratios=[1.55, 1, 1, 1], wspace=0.30)
    ax_e, key_e = make_umap_and_key_axes(middle[0, 0], fig, key_height=1.55)
    plot_publication_cluster_key_umap(
        ax_e,
        key_e,
        neuronal,
        ["3", "4", "5"],
        neuronal_bounds,
        "DIV30 MGEO neurons",
        key_fontsize=6.3,
    )
    panel_label(ax_e, "e")
    for idx, (feature, label) in enumerate(zip(GENES, ["f", "g", "h"], strict=True), start=1):
        ax, cax = make_umap_and_colorbar_axes(middle[0, idx], fig, bar_height=0.52)
        plot_feature_umap(
            fig, ax, cax, neuronal, feature, neuronal_bounds,
            expression_vmax, score_vmax, ratio_vmax,
        )
        panel_label(ax, label)

    # I-K: descriptive cluster-level distributions with organoid summaries.
    distributions = outer[2].subgridspec(1, 3, wspace=0.33)
    for idx, (gene, label) in enumerate(zip(GENES, ["i", "j", "k"], strict=True)):
        ax = fig.add_subplot(distributions[0, idx])
        plot_cluster_expression_panel(ax, neuronal, organoid_summary, gene, expression_ymax, random_state + idx)
        panel_label(ax, label)

    return save_figure(fig, final_dir, PRIMARY_STEM)


def plot_split_feature_figure(
    neuronal: pd.DataFrame,
    final_dir: Path,
) -> list[Path]:
    expression_vmax, score_vmax, ratio_vmax, _ = common_limits(neuronal)
    bounds = umap_bounds(neuronal)
    samples = sorted(pd.unique(neuronal["sample"]))
    features = ["LHX6", "LHX8", "PLXNA2", SCORE_COL, RATIO_COL]
    fig = plt.figure(figsize=(9.6, 8.0), facecolor="white")
    grid = fig.add_gridspec(
        5, 7, left=0.065, right=0.965, bottom=0.06, top=0.94,
        width_ratios=[1, 1, 1, 1, 1, 1, 0.055], wspace=0.12, hspace=0.2,
    )
    for row, feature in enumerate(features):
        mappable = None
        for col, sample in enumerate(samples):
            ax = fig.add_subplot(grid[row, col])
            sample_data = neuronal.loc[neuronal["sample"].eq(sample)].sort_values(feature, kind="stable")
            if feature == SCORE_COL:
                cmap, vmin, vmax = "coolwarm", -score_vmax, score_vmax
            elif feature == RATIO_COL:
                cmap, vmin, vmax = "coolwarm", -ratio_vmax, ratio_vmax
            else:
                cmap, vmin, vmax = expression_floor_colormap(expression_vmax), 0.0, expression_vmax
            mappable = ax.scatter(
                sample_data["umap_1"], sample_data["umap_2"], c=sample_data[feature], s=0.75,
                cmap=cmap, vmin=vmin, vmax=vmax, linewidths=0, rasterized=True,
            )
            style_umap_axis(ax, bounds)
            if row == 0:
                ax.set_title(str(sample).replace("9583-", ""), fontsize=6.5, fontweight="bold", pad=2)
            if col == 0:
                row_label = SCORE_LABEL if feature == SCORE_COL else RATIO_LABEL if feature == RATIO_COL else feature
                ax.text(-0.11, 0.5, row_label, transform=ax.transAxes, rotation=90, ha="center", va="center", fontsize=6.5, fontweight="bold")
        cax = fig.add_subplot(grid[row, 6])
        ticks = (
            [-score_vmax, 0, score_vmax]
            if feature == SCORE_COL
            else [-ratio_vmax, 0, ratio_vmax]
            if feature == RATIO_COL
            else [0, expression_vmax]
        )
        colorbar = fig.colorbar(mappable, cax=cax, orientation="vertical", ticks=ticks)
        colorbar.ax.tick_params(labelsize=4.5, length=1.5, width=0.4, pad=1)
        colorbar.outline.set_linewidth(0.35)
        colorbar.set_label(
            SCORE_LABEL if feature == SCORE_COL else RATIO_LABEL if feature == RATIO_COL else "log1p(CP10K)",
            fontsize=4.8,
            labelpad=1,
        )
    fig.suptitle(
        "Publication clusters 3, 4, and 5 · feature expression split by organoid",
        fontsize=9,
        fontweight="bold",
        y=0.975,
    )
    return save_figure(fig, final_dir, SPLIT_FEATURE_STEM)


def plot_split_relationship_figure(
    neuronal: pd.DataFrame,
    final_dir: Path,
    predictor: str,
    predictor_label: str,
    predictor_vmax: float,
    stem: str,
    title: str,
) -> list[Path]:
    plx_limit = max(2.0, float(np.nanquantile(neuronal["PLXNA2"], 0.995)) * 1.04)
    samples = sorted(pd.unique(neuronal["sample"]))
    fig, axes = plt.subplots(3, 6, figsize=(10.5, 6.6), sharex=True, sharey=True)
    for row, cluster in enumerate(NEURONAL_CLUSTERS):
        publication_cluster = PUBLICATION_CLUSTER_BY_SEURAT[cluster]
        cluster_name = PUBLICATION_CLUSTER_LABELS[publication_cluster]
        for col, sample in enumerate(samples):
            ax = axes[row, col]
            subset = neuronal.loc[neuronal["cluster"].eq(cluster) & neuronal["sample"].eq(sample)]
            ax.hexbin(
                subset[predictor], subset["PLXNA2"], gridsize=26, mincnt=1,
                extent=(-predictor_vmax, predictor_vmax, 0, plx_limit), cmap="Blues", linewidths=0, rasterized=True,
            )
            ax.set_xlim(-predictor_vmax, predictor_vmax)
            ax.set_ylim(0, plx_limit)
            if row == 0:
                ax.set_title(str(sample).replace("9583-", ""), fontsize=6.4, fontweight="bold")
            if col == 0:
                ax.set_ylabel(
                    f"{publication_cluster}  {cluster_name}\nPLXNA2",
                    fontsize=5.8,
                    color=PUBLICATION_CLUSTER_COLORS[publication_cluster],
                )
            if row == 2:
                ax.set_xlabel(predictor_label, fontsize=5.3)
            style_axis(ax)
    fig.suptitle(title, fontsize=9, fontweight="bold", y=0.985)
    fig.subplots_adjust(left=0.09, right=0.985, bottom=0.09, top=0.935, wspace=0.16, hspace=0.22)
    return save_figure(fig, final_dir, stem)


def write_captions(
    final_dir: Path,
    expression_vmax: float,
    score_vmax: float,
    ratio_vmax: float,
    ratio_epsilon_cp10k: float,
) -> None:
    primary = f"""Figure. Developmentally restricted evaluation of LHX6, LHX8, and PLXNA2 in Varela DIV30 organoids.

(a) All Varela DIV30 cells on the original Seurat UMAP, colored and labeled as the five DIV30 MGEO populations. Original Seurat clusters 0, 3, and 7 are combined visually as population 1, Radial Glia; Seurat cluster 6 maps to population 2, Inhibitory Progenitors; Seurat clusters 1, 4, and 2 map to populations 3 (SST+ cIN), 4 (PV Neuron Precursor), and 5 (MGE Subpallial Neuron), respectively. (b-d) LHX6, LHX8, and PLXNA2 on the same full-dataset UMAP. (e) Primary neuronal subset restricted to populations 3, 4, and 5; populations 1 and 2 are excluded from panels e-k. (f-h) LHX6, LHX8, and PLXNA2 within the neuronal subset, shown with a shared 0 to {expression_vmax:.1f} log1p(CP10K) expression range. (i-k) Cell-level expression distributions by population; black points and white bars show cell medians and interquartile ranges, and the six organoid-level means remain visible as individual points. No pooled-cell significance tests are performed. Score, ratio, and organoid-resolved relationship analyses remain available in the existing companion figures and exported tables but are intentionally omitted from this compact primary figure.
"""
    companions = """Companion figure 1. Neuronal-subset feature UMAPs split by organoid.

Rows show LHX6, LHX8, PLXNA2, z(LHX6)-z(LHX8), and the stabilized LHX6/LHX8 log2 ratio; columns show MW-1 through MW-6. All panels retain the same neuronal-subset UMAP coordinates and row-specific shared color limits used in the primary figure.

Companion figure 2. PLXNA2 versus z(LHX6)-z(LHX8) split by cluster and organoid.

Rows show publication clusters 3, 4, and 5; columns show MW-1 through MW-6. Hexagon intensity represents cell density. Descriptive slopes and correlations are exported in the underlying table, but fitted lines and p-values are intentionally omitted from this organoid-resolved figure.

Companion figure 3. PLXNA2 versus the stabilized LHX6/LHX8 log2 ratio split by cluster and organoid.

Rows show publication clusters 3, 4, and 5; columns show MW-1 through MW-6. Hexagon intensity represents cell density. The ratio uses linear CP10K values with ε={ratio_epsilon_cp10k:g} CP10K. Descriptive slopes and correlations are exported, but fitted lines and p-values are omitted.
"""
    (final_dir / "figure_caption.txt").write_text(primary)
    (final_dir / "companion_figure_captions.txt").write_text(companions)


def write_readme(
    final_dir: Path,
    project_root: Path,
    expression_vmax: float,
    score_vmax: float,
    ratio_vmax: float,
    ratio_epsilon_cp10k: float,
) -> None:
    text = f"""# Varela DIV30 neuronal LHX6/LHX8/PLXNA2 analysis

This package preserves the full DIV30 cluster and feature UMAPs as context and uses the five-population publication mapping for all visual labels and colors.

Primary comparison:
- Publication cluster 1 (`Radial Glia`) combines original Seurat clusters 0, 3, and 7.
- Publication cluster 2 (`Inhibitory Progenitors`) corresponds to original Seurat cluster 6.
- Publication clusters 3, 4, and 5 correspond to SST+ cIN (Seurat 1), PV Neuron Precursor (Seurat 4), and MGE Subpallial Neuron (Seurat 2).
- The neuronal analysis uses publication clusters 3, 4, and 5; publication clusters 1 and 2 are excluded.
- No LHX6/LHX8 high/low states are defined and no pooled-cell significance tests are run.

Workflow contract:
- Source: `{project_root}/results/python_anndata/varela_div30.h5ad`.
- Cluster names: `{project_root / CLUSTER_MAPPING_RELATIVE_PATH}`.
- Expression: existing log1p(CP10K), extracted through the shared cross-study marker workflow.
- Score: `z(LHX6) - z(LHX8)`; both population z-scores are calculated across all cells in clusters 1, 2, and 4.
- Orthogonal ratio: `log2((expm1(LHX6)+epsilon)/(expm1(LHX8)+epsilon))`, using linear CP10K back-transformed from stored log1p(CP10K), with epsilon = `{ratio_epsilon_cp10k:g}` CP10K.
- Shared neuronal feature range: 0 to `{expression_vmax:.6g}` log1p(CP10K).
- Symmetric score range: `{-score_vmax:.6g}` to `{score_vmax:.6g}`.
- Symmetric log2-ratio range: `{-ratio_vmax:.6g}` to `{ratio_vmax:.6g}`.
- Biological replicate: organoid/sample, n=6.

Figures:
- `figures/png/{PRIMARY_STEM}.png`: primary multi-panel figure.
- `figures/png/{SPLIT_FEATURE_STEM}.png`: five features split by organoid.
- `figures/png/{SPLIT_RELATIONSHIP_STEM}.png`: PLXNA2 versus score by cluster and organoid.
- `figures/png/{SPLIT_RATIO_RELATIONSHIP_STEM}.png`: PLXNA2 versus stabilized LHX6/LHX8 log2 ratio by cluster and organoid.

All figures are also exported as PDF and editable-text SVG. Dense cell layers are rasterized inside the vector exports. Tables include the exact per-cell neuronal dataset, cluster and organoid summaries, score parameters, and descriptive relationship summaries.
"""
    (final_dir / "README.md").write_text(text)


def write_provenance(
    final_dir: Path,
    project_root: Path,
    source_h5ad: Path,
    expression_vmax: float,
    score_vmax: float,
    ratio_vmax: float,
    ratio_epsilon_cp10k: float,
    outputs: list[Path],
) -> None:
    provenance = final_dir / "provenance"
    code_dir = final_dir / "code"
    provenance.mkdir(parents=True, exist_ok=True)
    code_dir.mkdir(parents=True, exist_ok=True)
    script_path = Path(__file__).resolve()
    shutil.copy2(script_path, code_dir / script_path.name)
    git_commit = subprocess.run(["git", "rev-parse", "HEAD"], check=False, capture_output=True, text=True)
    git_status = subprocess.run(["git", "status", "--short"], check=False, capture_output=True, text=True)
    (provenance / "git_commit.txt").write_text(git_commit.stdout if git_commit.returncode == 0 else git_commit.stderr)
    (provenance / "git_status_short.txt").write_text(git_status.stdout if git_status.returncode == 0 else git_status.stderr)
    rows = [
        ("rendered_at", datetime.now().astimezone().isoformat()),
        ("project_root", str(project_root)),
        ("source_h5ad", str(source_h5ad)),
        ("cluster_mapping_source", str(cluster_mapping_path(project_root))),
        (
            "publication_visual_mapping",
            "Seurat 0/3/7->publication 1 Radial Glia; 6->2 Inhibitory Progenitors; "
            "1->3 SST+ cIN; 4->4 PV Neuron Precursor; 2->5 MGE Subpallial Neuron",
        ),
        ("study_id", "varela_div30"),
        ("genes", ",".join(GENES)),
        ("primary_included_clusters", ",".join(NEURONAL_CLUSTERS)),
        ("primary_excluded_clusters", ",".join(EXCLUDED_PRIMARY_CLUSTERS)),
        ("score_definition", "z(LHX6)-z(LHX8); population z-scores across clusters 1,2,4"),
        ("ratio_definition", "log2((expm1(LHX6_log1p_cp10k)+epsilon)/(expm1(LHX8_log1p_cp10k)+epsilon))"),
        ("ratio_epsilon_cp10k", f"{ratio_epsilon_cp10k:.9g}"),
        ("shared_expression_vmax_log1p_cp10k", f"{expression_vmax:.9g}"),
        ("symmetric_score_vmax", f"{score_vmax:.9g}"),
        ("symmetric_ratio_vmax", f"{ratio_vmax:.9g}"),
        ("pooled_cell_significance_tests", "none"),
        ("biological_replicate", "organoid/sample; n=6"),
        ("png_dpi", str(PNG_EXPORT_DPI)),
        ("pdf_dpi", str(PDF_EXPORT_DPI)),
        ("svg_dpi", str(SVG_EXPORT_DPI)),
        ("svg_text_policy", "editable text; dense cell layers rasterized"),
    ]
    rows.extend((f"output_{idx + 1}", str(path)) for idx, path in enumerate(outputs))
    pd.DataFrame(rows, columns=["key", "value"]).to_csv(provenance / "render_manifest.tsv", sep="\t", index=False)
    files = [
        final_dir / "README.md",
        final_dir / "figure_caption.txt",
        final_dir / "companion_figure_captions.txt",
        code_dir / script_path.name,
        *outputs,
    ]
    files.extend(sorted((final_dir / "tables").glob("*")))
    sha_lines: list[str] = []
    for path in files:
        if not path.is_file():
            continue
        digest = subprocess.run(["sha256sum", str(path)], check=False, capture_output=True, text=True)
        if digest.returncode == 0:
            sha_lines.append(digest.stdout)
    (provenance / "sha256_manifest.txt").write_text("".join(sha_lines))


def main() -> None:
    args = parse_args()
    project_root = args.project_root.resolve()
    final_dir = args.final_dir or (project_root / "final_figures" / FINAL_FOLDER_DEFAULT)
    if args.fresh and final_dir.exists():
        shutil.rmtree(final_dir)
    for subdir in ("figures/png", "figures/pdf", "figures/svg", "tables", "logs", "provenance", "code"):
        (final_dir / subdir).mkdir(parents=True, exist_ok=True)

    spec = div30_spec(project_root)
    raw_path = final_dir / "tables" / "varela_div30_lhx6_lhx8_plxna2_per_cell.tsv.gz"
    raw, gene_matches = extract_marker_expression_from_h5ad(spec, raw_path, project_root=project_root, genes=GENES)
    gene_matches.to_csv(final_dir / "tables" / "gene_matches.tsv", sep="\t", index=False)
    filtered, filter_summary = _apply_internal_umap_plot_filters(raw)
    filter_summary.to_csv(final_dir / "tables" / "plot_filter_summary.tsv", sep="\t", index=False)
    annotated, mapping = add_cluster_annotations(filtered, project_root)
    mapping.to_csv(final_dir / "tables" / "cluster_annotation_mapping.tsv", sep="\t", index=False)
    cluster_annotation_summary(annotated).to_csv(
        final_dir / "tables" / "cluster_annotation_summary.tsv", sep="\t", index=False
    )
    publication_cluster_annotation_summary(annotated).to_csv(
        final_dir / "tables" / "publication_cluster_annotation_summary.tsv",
        sep="\t",
        index=False,
    )
    neuronal, score_parameters, ratio_parameters = add_neuronal_score(
        annotated, args.ratio_epsilon_cp10k
    )
    score_parameters.to_csv(final_dir / "tables" / "lhx6_lhx8_zscore_parameters.tsv", sep="\t", index=False)
    ratio_parameters.to_csv(final_dir / "tables" / "lhx6_lhx8_log2_ratio_parameters.tsv", sep="\t", index=False)

    per_cell_cols = [
        "cell_id", "sample", "cluster", "paper_cluster_annotation", "cluster_name_label",
        "publication_cluster", "publication_cluster_annotation", "publication_cluster_label",
        "lineage_comparison_group", "umap_1", "umap_2", "LHX6", "LHX8", "PLXNA2",
        "LHX6_zscore", "LHX8_zscore", SCORE_COL, RATIO_COL,
    ]
    neuronal[per_cell_cols].to_csv(
        final_dir / "tables" / "neuronal_clusters_1_2_4_per_cell.tsv.gz", sep="\t", index=False
    )
    cluster_summary, organoid_summary = expression_summaries(neuronal)
    cluster_summary.to_csv(final_dir / "tables" / "neuronal_cluster_expression_summary.tsv", sep="\t", index=False)
    organoid_summary.to_csv(final_dir / "tables" / "neuronal_cluster_organoid_expression_summary.tsv", sep="\t", index=False)
    lineage_summary, lineage_organoid_summary = lineage_group_expression_summaries(neuronal)
    lineage_summary.to_csv(final_dir / "tables" / "neuronal_lineage_group_expression_summary.tsv", sep="\t", index=False)
    lineage_organoid_summary.to_csv(
        final_dir / "tables" / "neuronal_lineage_group_organoid_expression_summary.tsv", sep="\t", index=False
    )
    pooled_relationship, organoid_relationship = relationship_summaries(
        neuronal,
        SCORE_COL,
        "z(LHX6)-z(LHX8); z-scores computed across clusters 1,2,4 pooled",
    )
    pooled_relationship.to_csv(final_dir / "tables" / "plxna2_score_relationship_by_cluster.tsv", sep="\t", index=False)
    organoid_relationship.to_csv(
        final_dir / "tables" / "plxna2_score_relationship_by_cluster_organoid.tsv", sep="\t", index=False
    )
    pooled_ratio_relationship, organoid_ratio_relationship = relationship_summaries(
        neuronal,
        RATIO_COL,
        (
            "log2((expm1(LHX6_log1p_cp10k)+epsilon)/(expm1(LHX8_log1p_cp10k)+epsilon)); "
            f"epsilon={args.ratio_epsilon_cp10k:g} CP10K"
        ),
    )
    pooled_ratio_relationship.to_csv(
        final_dir / "tables" / "plxna2_ratio_relationship_by_cluster.tsv", sep="\t", index=False
    )
    organoid_ratio_relationship.to_csv(
        final_dir / "tables" / "plxna2_ratio_relationship_by_cluster_organoid.tsv", sep="\t", index=False
    )

    expression_vmax, score_vmax, ratio_vmax, _ = common_limits(neuronal)
    scale_manifest = pd.DataFrame(
        [
            {
                "plot_family": "LHX6_LHX8_PLXNA2_feature_umaps",
                "included_clusters": "1,2,4",
                "vmin": 0.0,
                "vmax": expression_vmax,
                "scale": "log1p(CP10K)",
                "color_rule": "gray through 1.0; blue to shared vmax",
            },
            {
                "plot_family": "LHX6_minus_LHX8_score_umaps",
                "included_clusters": "1,2,4",
                "vmin": -score_vmax,
                "vmax": score_vmax,
                "scale": SCORE_LABEL,
                "color_rule": "continuous coolwarm diverging scale centered at zero",
            },
            {
                "plot_family": "LHX6_LHX8_log2_ratio_umaps",
                "included_clusters": "1,2,4",
                "vmin": -ratio_vmax,
                "vmax": ratio_vmax,
                "scale": RATIO_LABEL,
                "color_rule": "continuous coolwarm diverging scale centered at zero",
            },
        ]
    )
    scale_manifest.to_csv(final_dir / "tables" / "feature_plot_scale_manifest.tsv", sep="\t", index=False)

    outputs: list[Path] = []
    outputs.extend(plot_primary_figure(annotated, neuronal, organoid_summary, final_dir, args.random_state))
    outputs.extend(plot_split_feature_figure(neuronal, final_dir))
    outputs.extend(
        plot_split_relationship_figure(
            neuronal,
            final_dir,
            SCORE_COL,
            SCORE_LABEL,
            score_vmax,
            SPLIT_RELATIONSHIP_STEM,
            "PLXNA2 versus LHX6-minus-LHX8 z-score by cluster and organoid",
        )
    )
    outputs.extend(
        plot_split_relationship_figure(
            neuronal,
            final_dir,
            RATIO_COL,
            RATIO_LABEL,
            ratio_vmax,
            SPLIT_RATIO_RELATIONSHIP_STEM,
            "PLXNA2 versus stabilized LHX6/LHX8 log2 ratio by cluster and organoid",
        )
    )
    write_captions(
        final_dir, expression_vmax, score_vmax, ratio_vmax, args.ratio_epsilon_cp10k
    )
    write_readme(
        final_dir,
        project_root,
        expression_vmax,
        score_vmax,
        ratio_vmax,
        args.ratio_epsilon_cp10k,
    )
    write_provenance(
        final_dir,
        project_root,
        spec.resolved_h5ad_path(project_root),
        expression_vmax,
        score_vmax,
        ratio_vmax,
        args.ratio_epsilon_cp10k,
        outputs,
    )

    print(cluster_annotation_summary(annotated).to_string(index=False), flush=True)
    print(score_parameters.to_string(index=False), flush=True)
    print(pooled_relationship.to_string(index=False), flush=True)
    print(pooled_ratio_relationship.to_string(index=False), flush=True)
    for output in outputs:
        print(f"Rendered: {output}", flush=True)
    print(f"Final dir: {final_dir}", flush=True)


if __name__ == "__main__":
    main()
