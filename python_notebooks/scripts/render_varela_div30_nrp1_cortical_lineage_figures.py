#!/usr/bin/env python3
"""Render NRP1/NRP2 lineage-expression and replicate-aware figures for Varela DIV30 MGEOs."""

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
import numpy as np
import pandas as pd
from scipy import stats

from mge_organoid_python.cross_study_marker_expression import (
    EXPRESSION_COLOR_FLOOR,
    _apply_internal_umap_plot_filters,
    extract_marker_expression_from_h5ad,
)
from render_varela_div30_lhx6_lhx8_plxna2_final_figure import (
    PUBLICATION_CLUSTER_BY_SEURAT,
    PUBLICATION_CLUSTER_COLORS,
    PUBLICATION_CLUSTER_LABELS,
    add_cluster_annotations,
    div30_spec,
    make_umap_and_colorbar_axes,
    make_umap_and_key_axes,
    panel_label,
    plot_feature_umap,
    plot_publication_cluster_key_umap,
    round_scale_limit,
    style_axis,
    style_umap_axis,
    umap_bounds,
)
from render_varela_div30_plxna2_replicate_aware_figure import (
    BOOTSTRAP_ITERATIONS,
    bootstrap_paired_difference,
    paired_wilcoxon_exact,
)


PROJECT_ROOT_DEFAULT = Path("/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder")
FINAL_DIR_DEFAULT = PROJECT_ROOT_DEFAULT / "final_figures" / "fig_varela_div30_lhx6_lhx8_plxna2_v1_candidate"
SCRIPT_NAME = "render_varela_div30_nrp1_cortical_lineage_figures.py"
EXPRESSION_STEM = "varela_div30_nrp1_cortical_lineage_expression"
STATS_STEM = "varela_div30_nrp1_cortical_vs_subpallial_replicate_aware"
TABLE_PREFIX = "nrp1_cortical_lineage"
DUAL_EXPRESSION_STEM = "varela_div30_nrp1_nrp2_cortical_lineage_expression"
DUAL_STATS_STEM = "varela_div30_nrp1_nrp2_cortical_vs_subpallial_replicate_aware"
GENES = ["NRP1", "NRP2"]
EXPORT_DPI = 300
NEURONAL_CLUSTERS = ["1", "4", "2"]
CLUSTER_LABELS = {
    "1": "SST+ cIN",
    "4": "PV neuron precursor",
    "2": "MGE subpallial neurons",
}
LINEAGE_ORDER = ["cortical_interneuron_lineage", "subpallial_MGE_lineage"]
LINEAGE_LABELS = {
    "cortical_interneuron_lineage": "Cortical interneuron lineage",
    "subpallial_MGE_lineage": "Subpallial MGE neurons",
}
LINEAGE_SHORT_LABELS = {
    "cortical_interneuron_lineage": "Cortical interneuron\nlineage",
    "subpallial_MGE_lineage": "Subpallial MGE\nneurons",
}
LINEAGE_COLORS = {
    "cortical_interneuron_lineage": "#009E73",
    "subpallial_MGE_lineage": PUBLICATION_CLUSTER_COLORS["5"],
}
CLUSTER_TO_LINEAGE = {
    "1": "cortical_interneuron_lineage",
    "4": "cortical_interneuron_lineage",
    "2": "subpallial_MGE_lineage",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT_DEFAULT)
    parser.add_argument("--final-dir", type=Path, default=FINAL_DIR_DEFAULT)
    parser.add_argument("--random-state", type=int, default=3901)
    parser.add_argument("--bootstrap-iterations", type=int, default=BOOTSTRAP_ITERATIONS)
    return parser.parse_args()


def load_data(project_root: Path, final_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    raw_path = final_dir / "tables" / "varela_div30_nrp1_nrp2_per_cell.tsv.gz"
    raw, matches = extract_marker_expression_from_h5ad(
        div30_spec(project_root), raw_path, project_root=project_root, genes=GENES
    )
    filtered, filter_summary = _apply_internal_umap_plot_filters(raw)
    annotated, _ = add_cluster_annotations(filtered, project_root)
    annotated["cluster"] = annotated["cluster"].astype(str)
    neuronal = annotated.loc[annotated["cluster"].isin(NEURONAL_CLUSTERS)].copy()
    neuronal["organoid"] = neuronal["sample"].astype(str)
    neuronal["lineage"] = neuronal["cluster"].map(CLUSTER_TO_LINEAGE)
    if neuronal["lineage"].isna().any():
        raise ValueError("A selected neuronal cell lacks a lineage mapping.")
    observed = dict(
        neuronal[["cluster", "paper_cluster_annotation"]]
        .drop_duplicates()
        .set_index("cluster")["paper_cluster_annotation"]
    )
    if observed != CLUSTER_LABELS:
        raise ValueError(f"Canonical neuronal annotations changed: {observed}")
    return annotated, neuronal, pd.concat(
        [matches.assign(record_type="gene_match"), filter_summary.assign(record_type="plot_filter")],
        ignore_index=True,
        sort=False,
    )


def summarize_groups(data: pd.DataFrame, group_col: str, labels: dict[str, str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (organoid, group), frame in data.groupby(["organoid", group_col], sort=True):
        row: dict[str, object] = {
            "organoid": organoid,
            group_col: group,
            f"{group_col}_label": labels[group],
            "n_cells": len(frame),
        }
        for gene in GENES:
            values = frame[gene].to_numpy(dtype=float)
            q25, q75 = np.quantile(values, [0.25, 0.75])
            row.update(
                {
                    f"mean_{gene}_log1p_CP10K": float(np.mean(values)),
                    f"median_{gene}_log1p_CP10K": float(np.median(values)),
                    f"sd_{gene}_log1p_CP10K": float(np.std(values, ddof=1)),
                    f"q25_{gene}_log1p_CP10K": float(q25),
                    f"q75_{gene}_log1p_CP10K": float(q75),
                    f"iqr_{gene}_log1p_CP10K": float(q75 - q25),
                    f"pct_{gene}_detected": float(100 * np.mean(values > 0)),
                }
            )
        rows.append(row)
    return pd.DataFrame(rows)


def completeness(lineage_summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for organoid in sorted(lineage_summary["organoid"].unique()):
        observed = set(lineage_summary.loc[lineage_summary["organoid"] == organoid, "lineage"])
        missing = [lineage for lineage in LINEAGE_ORDER if lineage not in observed]
        rows.append(
            {
                "organoid": organoid,
                "complete_pair": not missing,
                "missing_lineages": ";".join(missing) if missing else "none",
                "paired_analysis_consequence": (
                    "included in paired analysis" if not missing else "excluded; missing lineage reported"
                ),
            }
        )
    return pd.DataFrame(rows)


def analyze_endpoint(
    summary: pd.DataFrame,
    gene: str,
    endpoint: str,
    label: str,
    role: str,
    rng: np.random.Generator,
    bootstrap_iterations: int,
) -> dict[str, object]:
    pivot = summary.pivot(index="organoid", columns="lineage", values=endpoint)
    pivot = pivot.dropna(subset=LINEAGE_ORDER).sort_index()
    cortical = pivot[LINEAGE_ORDER[0]].to_numpy(dtype=float)
    subpallial = pivot[LINEAGE_ORDER[1]].to_numpy(dtype=float)
    result = paired_wilcoxon_exact(cortical, subpallial)
    bootstrap = bootstrap_paired_difference(cortical - subpallial, bootstrap_iterations, rng)
    if result["mean_paired_difference"] > 0:
        direction = "cortical interneuron lineage higher"
    elif result["mean_paired_difference"] < 0:
        direction = "subpallial MGE lineage higher"
    else:
        direction = "no mean difference"
    return {
        "endpoint": endpoint,
        "gene": gene,
        "endpoint_label": label,
        "analysis_role": role,
        "contrast": "cortical_interneuron_lineage_vs_subpallial_MGE_lineage",
        "difference_definition": "cortical interneuron lineage minus subpallial MGE lineage",
        "effect_direction": direction,
        **result,
        **bootstrap,
        "multiplicity_adjustment": (
            "none; NRP1 retains its prespecified contrast; NRP2 is an added exploratory parallel analysis"
        ),
    }


def group_mean_ci(lineage_summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for gene in GENES:
        for lineage in LINEAGE_ORDER:
            values = lineage_summary.loc[
                lineage_summary["lineage"] == lineage, f"mean_{gene}_log1p_CP10K"
            ].to_numpy(dtype=float)
            mean = float(np.mean(values))
            sd = float(np.std(values, ddof=1))
            half_width = float(stats.t.ppf(0.975, len(values) - 1) * sd / np.sqrt(len(values)))
            rows.append(
                {
                    "gene": gene,
                    "lineage": lineage,
                    "lineage_label": LINEAGE_LABELS[lineage],
                    "n_organoids": len(values),
                    "mean_of_organoid_lineage_means_log1p_CP10K": mean,
                    "mean_95CI_low": mean - half_width,
                    "mean_95CI_high": mean + half_width,
                    "ci_method": "two-sided t confidence interval across organoid-level means",
                }
            )
    return pd.DataFrame(rows)


def equal_cluster_weight_sensitivity(
    cluster_summary: pd.DataFrame,
    gene: str,
    rng: np.random.Generator,
    bootstrap_iterations: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    pivot = cluster_summary.pivot(
        index="organoid", columns="cluster", values=f"mean_{gene}_log1p_CP10K"
    ).sort_index()
    values = pd.DataFrame(
        {
            "organoid": pivot.index,
            f"equal_weight_cortical_mean_{gene}_log1p_CP10K": (
                pivot["1"].to_numpy(dtype=float) + pivot["4"].to_numpy(dtype=float)
            )
            / 2,
            f"subpallial_MGE_mean_{gene}_log1p_CP10K": pivot["2"].to_numpy(dtype=float),
        }
    )
    cortical = values[f"equal_weight_cortical_mean_{gene}_log1p_CP10K"].to_numpy(dtype=float)
    subpallial = values[f"subpallial_MGE_mean_{gene}_log1p_CP10K"].to_numpy(dtype=float)
    result = paired_wilcoxon_exact(cortical, subpallial)
    bootstrap = bootstrap_paired_difference(cortical - subpallial, bootstrap_iterations, rng)
    statistics = pd.DataFrame(
        [
            {
                "analysis_role": "pooling-definition sensitivity",
                "gene": gene,
                "cortical_definition": "equal-weight mean of organoid cluster-1 and cluster-4 means",
                "subpallial_definition": "organoid cluster-2 mean",
                **result,
                **bootstrap,
            }
        ]
    )
    return values, statistics


def plot_lineage_umap(ax: plt.Axes, key_ax: plt.Axes, neuronal: pd.DataFrame) -> None:
    for lineage in LINEAGE_ORDER:
        frame = neuronal.loc[neuronal["lineage"] == lineage]
        ax.scatter(
            frame["umap_1"], frame["umap_2"], s=0.55, color=LINEAGE_COLORS[lineage],
            alpha=0.78, linewidths=0, rasterized=True,
        )
    style_umap_axis(ax, umap_bounds(neuronal))
    ax.set_title("DIV30 MGEO lineages", fontsize=8.2, fontweight="bold", pad=3)
    key_ax.set_axis_off()
    labels = ["3 + 4  Cortical interneuron lineage", "5  Subpallial MGE neurons"]
    for row, (lineage, label) in enumerate(zip(LINEAGE_ORDER, labels, strict=True)):
        key_ax.text(
            0.02, 0.78 - row * 0.53, label, transform=key_ax.transAxes,
            color=LINEAGE_COLORS[lineage], fontsize=6.3, fontweight="bold",
            ha="left", va="center",
        )


def plot_violin(
    ax: plt.Axes,
    neuronal: pd.DataFrame,
    summary: pd.DataFrame,
    group_col: str,
    order: list[str],
    labels: list[str],
    colors: dict[str, str],
    random_state: int,
    y_max: float,
    gene: str,
) -> None:
    rng = np.random.default_rng(random_state)
    values = [neuronal.loc[neuronal[group_col] == group, gene].to_numpy(dtype=float) for group in order]
    draw = [rng.choice(v, size=min(v.size, 12_000), replace=False) for v in values]
    positions = np.arange(len(order), dtype=float)
    violins = ax.violinplot(draw, positions=positions, widths=0.72, showextrema=False)
    for body, group in zip(violins["bodies"], order, strict=True):
        body.set_facecolor(colors[group])
        body.set_edgecolor("none")
        body.set_alpha(0.58)
    for xpos, group, group_values in zip(positions, order, values, strict=True):
        q25, median, q75 = np.quantile(group_values, [0.25, 0.5, 0.75])
        ax.plot([xpos, xpos], [q25, q75], color="white", linewidth=2.5, solid_capstyle="round")
        ax.scatter([xpos], [median], s=13, color="#202020", edgecolor="white", linewidth=0.35, zorder=4)
        means = summary.loc[
            summary[group_col] == group, f"mean_{gene}_log1p_CP10K"
        ].to_numpy(dtype=float)
        ax.scatter(
            xpos + np.linspace(-0.13, 0.13, len(means)), means, s=12,
            facecolor="white", edgecolor="#202020", linewidth=0.45, zorder=5,
        )
    ax.set_xticks(positions, labels, linespacing=0.95)
    ax.set_ylim(0, y_max)
    ax.set_ylabel("Expression, log1p(CP10K)", fontsize=6.3)
    ax.set_title(rf"$\it{{{gene}}}$", fontsize=7.4, fontweight="bold", pad=5)
    style_axis(ax)


def save_figure(
    fig: plt.Figure, final_dir: Path, stem: str, aliases: tuple[str, ...] = ()
) -> list[Path]:
    outputs = []
    for extension in ["png", "pdf", "svg"]:
        for output_stem in (stem, *aliases):
            output = final_dir / "figures" / extension / f"{output_stem}.{extension}"
            fig.savefig(output, dpi=EXPORT_DPI, bbox_inches="tight", facecolor="white")
            outputs.append(output)
    plt.close(fig)
    return outputs


def render_expression_figure(
    full: pd.DataFrame,
    neuronal: pd.DataFrame,
    cluster_summary: pd.DataFrame,
    lineage_summary: pd.DataFrame,
    final_dir: Path,
    random_state: int,
) -> tuple[list[Path], float]:
    vmax = round_scale_limit(
        max(EXPRESSION_COLOR_FLOOR, *(float(np.quantile(full[gene], 0.99)) for gene in GENES))
    )
    y_max = {
        gene: max(1.0, float(np.quantile(neuronal[gene], 0.995)) * 1.05)
        for gene in GENES
    }
    fig = plt.figure(figsize=(6.2, 10.3), facecolor="white")
    outer = fig.add_gridspec(
        5, 2, left=0.08, right=0.98, bottom=0.055, top=0.98,
        height_ratios=[1, 1, 1, 1.12, 1.12], wspace=0.30, hspace=0.58,
    )
    ax_a, key_a = make_umap_and_key_axes(outer[0, 0], fig, key_height=1.9)
    plot_publication_cluster_key_umap(
        ax_a, key_a, full, ["1", "2", "3", "4", "5"], umap_bounds(full),
        "DIV30 MGEOs", key_fontsize=6.3,
    )
    panel_label(ax_a, "a")
    ax_b, key_b = make_umap_and_key_axes(outer[0, 1], fig, key_height=1.45)
    plot_lineage_umap(ax_b, key_b, neuronal)
    panel_label(ax_b, "b")
    panel_index = 2
    for row, gene in enumerate(GENES, start=1):
        ax_full, cax_full = make_umap_and_colorbar_axes(outer[row, 0], fig, bar_height=0.52)
        plot_feature_umap(fig, ax_full, cax_full, full, gene, umap_bounds(full), vmax, 1.0, 1.0)
        ax_full.set_title(rf"$\it{{{gene}}}$ · all cells", fontsize=7.4, fontweight="bold", pad=4)
        panel_label(ax_full, chr(ord("a") + panel_index))
        panel_index += 1
        ax_neuronal, cax_neuronal = make_umap_and_colorbar_axes(
            outer[row, 1], fig, bar_height=0.52
        )
        plot_feature_umap(
            fig, ax_neuronal, cax_neuronal, neuronal, gene,
            umap_bounds(neuronal), vmax, 1.0, 1.0,
        )
        ax_neuronal.set_title(
            rf"$\it{{{gene}}}$ · neuronal lineages", fontsize=7.4,
            fontweight="bold", pad=4,
        )
        panel_label(ax_neuronal, chr(ord("a") + panel_index))
        panel_index += 1
    for gene_index, gene in enumerate(GENES):
        row = 3 + gene_index
        ax_cluster = fig.add_subplot(outer[row, 0])
        plot_violin(
            ax_cluster, neuronal, cluster_summary, "cluster", NEURONAL_CLUSTERS,
            ["3\nSST+ cIN", "4\nPV precursor", "5\nMGE subpallial"],
            {
                cluster: PUBLICATION_CLUSTER_COLORS[PUBLICATION_CLUSTER_BY_SEURAT[cluster]]
                for cluster in NEURONAL_CLUSTERS
            },
            random_state + 2 * gene_index, y_max[gene], gene,
        )
        panel_label(ax_cluster, chr(ord("a") + panel_index))
        panel_index += 1
        ax_lineage = fig.add_subplot(outer[row, 1])
        plot_violin(
            ax_lineage, neuronal, lineage_summary, "lineage", LINEAGE_ORDER,
            ["3+4\nCortical interneuron\nlineage", "5\nSubpallial MGE\nneurons"],
            LINEAGE_COLORS, random_state + 2 * gene_index + 1, y_max[gene], gene,
        )
        panel_label(ax_lineage, chr(ord("a") + panel_index))
        panel_index += 1
    return save_figure(
        fig, final_dir, EXPRESSION_STEM, aliases=(DUAL_EXPRESSION_STEM,)
    ), vmax


def render_stats_figure(
    summary: pd.DataFrame,
    ci: pd.DataFrame,
    primary: pd.DataFrame,
    final_dir: Path,
) -> list[Path]:
    samples = sorted(summary["organoid"].unique())
    jitter = dict(zip(samples, np.linspace(-0.055, 0.055, len(samples))))
    group_x = {LINEAGE_ORDER[0]: 0.00, LINEAGE_ORDER[1]: 0.75}
    summary_offset = 0.13
    pivots = {
        gene: summary.pivot(
            index="organoid", columns="lineage", values=f"mean_{gene}_log1p_CP10K"
        )
        for gene in GENES
    }
    y_values = np.concatenate([pivot[LINEAGE_ORDER].to_numpy(dtype=float) for pivot in pivots.values()])
    y_min = max(0.0, float(np.min(y_values)) - 0.18 * float(np.ptp(y_values)))
    y_max = float(np.max(y_values)) + 0.28 * max(float(np.ptp(y_values)), 0.25)
    fig, axes = plt.subplots(1, 2, figsize=(6.2, 3.25), facecolor="white", sharey=True)
    for gene_index, (ax, gene) in enumerate(zip(axes, GENES, strict=True)):
        pivot = pivots[gene]
        ci_by_lineage = ci.loc[ci["gene"] == gene].set_index("lineage")
        gene_primary = primary.loc[primary["gene"] == gene].iloc[0]
        for sample in samples:
            x = np.asarray([group_x[lineage] + jitter[sample] for lineage in LINEAGE_ORDER])
            y = pivot.loc[sample, LINEAGE_ORDER].to_numpy(dtype=float)
            ax.plot(x, y, color="#BDBDBD", linewidth=0.8, alpha=0.7, zorder=1)
            ax.scatter(
                x, y, s=27, color=[LINEAGE_COLORS[lineage] for lineage in LINEAGE_ORDER],
                edgecolor="white", linewidth=0.5, zorder=3,
            )
        for lineage in LINEAGE_ORDER:
            row = ci_by_lineage.loc[lineage]
            x = group_x[lineage] + summary_offset
            ax.plot(
                [x, x], [row["mean_95CI_low"], row["mean_95CI_high"]],
                color="#111111", linewidth=1.8,
            )
            ax.scatter(
                [x], [row["mean_of_organoid_lineage_means_log1p_CP10K"]], marker="D",
                s=50, color="#111111", edgecolor="white", linewidth=0.5, zorder=5,
            )
        n_greater = int(gene_primary["n_first_greater_than_second"])
        n_pairs = int(gene_primary["n_paired_organoids"])
        direction_text = f"{n_greater}/{n_pairs} cortical > subpallial"
        ax.set_xticks(
            [group_x[lineage] for lineage in LINEAGE_ORDER],
            [LINEAGE_SHORT_LABELS[lineage] for lineage in LINEAGE_ORDER],
        )
        ax.set_xlim(-0.20, 1.05)
        ax.set_ylim(y_min, y_max)
        ax.set_title(rf"$\it{{{gene}}}$", fontsize=8, fontweight="bold", pad=4)
        ax.text(
            0.00, -0.34,
            f"Exact paired Wilcoxon P = {gene_primary['p_value_raw_exact_two_sided']:.3f}\n"
            f"Rank-biserial r = {gene_primary['rank_biserial_correlation']:.2f}; {direction_text}",
            transform=ax.transAxes, ha="left", va="top", fontsize=6.2, clip_on=False,
        )
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(False)
        ax.tick_params(width=0.6, length=2.5, color="#222222")
        panel_label(ax, chr(ord("a") + gene_index))
    axes[0].set_ylabel("Mean expression, log1p(CP10K)\nper organoid")
    fig.subplots_adjust(left=0.11, right=0.99, bottom=0.40, top=0.92, wspace=0.25)
    return save_figure(fig, final_dir, STATS_STEM, aliases=(DUAL_STATS_STEM,))


def git_value() -> str:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"], check=True, capture_output=True, text=True
        ).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unavailable"


def write_documentation(
    final_dir: Path,
    primary: pd.DataFrame,
    sensitivity: pd.DataFrame,
    vmax: float,
    equal_weight: pd.DataFrame,
) -> list[Path]:
    caption = f"""NRP1 and NRP2 expression in Varela DIV30 MGEO neuronal lineages.

Panels a-b show the five publication populations and the primary neuronal lineage subset. Panels c-f show NRP1 and NRP2 across all cells and within the neuronal subset using a shared 0 to {vmax:g} log1p(CP10K) scale. Panels g-j show descriptive cell-level distributions in the three established neuronal populations and after pooling SST+ cIN and PV precursor cells as the cortical interneuron lineage. Individual organoid means are overlaid; cells are not treated as independent biological replicates.

The replicate-aware figure shows one mean expression value per organoid and lineage for NRP1 and NRP2. Light-gray lines connect matched organoids, colored points denote established lineages, and displaced black diamonds with bars show the across-organoid mean and 95% confidence interval. Exact paired Wilcoxon results are reported beneath each gene. Clusters 0, 3, 6, and 7 are excluded, and neither receptor is used to redefine cluster identities. The NRP1 contrast retains its prior prespecified status; the added NRP2 contrast is exploratory.
"""
    caption_paths = [
        final_dir / "nrp1_cortical_lineage_figure_caption.txt",
        final_dir / "nrp1_nrp2_cortical_lineage_figure_caption.txt",
    ]
    for path in caption_paths:
        path.write_text(caption)
    sensitivity_lines = "\n".join(
        f"| {row.gene} | {row.endpoint_label} | {row.effect_direction} | {row.p_value_raw_exact_two_sided:.5f} | {row.rank_biserial_correlation:.2f} |"
        for row in sensitivity.itertuples()
    )
    primary_lines = []
    for gene in GENES:
        row = primary.loc[primary["gene"] == gene].iloc[0]
        equal = equal_weight.loc[equal_weight["gene"] == gene].iloc[0]
        direction = "higher" if row["mean_paired_difference"] > 0 else "lower"
        primary_lines.extend(
            [
                f"## {gene}",
                "",
                f"- Exact paired Wilcoxon W = {row['wilcoxon_statistic_W_min']:.1f}",
                f"- Exact two-sided P = {row['p_value_raw_exact_two_sided']:.5f}",
                f"- Rank-biserial r = {row['rank_biserial_correlation']:.2f}",
                f"- Mean paired difference (cortical minus subpallial) = {row['mean_paired_difference']:.6f} log1p(CP10K)",
                f"- Bootstrap 95% CI = [{row['mean_difference_bootstrap_95CI_low']:.6f}, {row['mean_difference_bootstrap_95CI_high']:.6f}]",
                f"- Organoids with cortical > subpallial = {int(row['n_first_greater_than_second'])}/{int(row['n_paired_organoids'])}",
                f"- Cortical lineage was {direction} on the observed organoid-level mean scale",
                f"- Equal-cluster-weight sensitivity: P = {equal['p_value_raw_exact_two_sided']:.5f}, r = {equal['rank_biserial_correlation']:.2f}",
                "",
            ]
        )
    results = f"""# Replicate-aware NRP1 and NRP2 lineage results

{chr(10).join(primary_lines)}
| Gene | Endpoint | Direction | Exact P | Rank-biserial r |
|---|---|---|---:|---:|
{sensitivity_lines}

The organoid is the inferential unit. NRP1 retains its prespecified analysis; NRP2 is an added exploratory parallel comparison. These analyses test association with established annotations and do not establish a causal role for either receptor.
"""
    results_paths = [
        final_dir / "nrp1_cortical_lineage_results_summary.md",
        final_dir / "nrp1_nrp2_cortical_lineage_results_summary.md",
    ]
    for path in results_paths:
        path.write_text(results)
    return [*caption_paths, *results_paths]


def write_manifest(final_dir: Path, outputs: list[Path], project_root: Path) -> None:
    manifest = pd.DataFrame(
        [
            ("rendered_at_local", datetime.now().isoformat(timespec="seconds")),
            ("source_h5ad", str(div30_spec(project_root).resolved_h5ad_path(project_root))),
            ("expression_scale", "existing log1p(CP10K)"),
            ("genes", ",".join(GENES)),
            ("export", "300-dpi PNG; PDF/SVG vector with editable text and rasterized cell layers"),
            ("inferential_unit", "organoid"),
            ("publication_mapping", "Seurat 0/3/7->1; 6->2; 1->3; 4->4; 2->5"),
            ("excluded_from_lineage_test", "Seurat clusters 0,3,6,7"),
            ("primary_test", "exact two-sided paired Wilcoxon signed-rank"),
            ("git_commit", git_value()),
        ],
        columns=["key", "value"],
    )
    manifest_path = final_dir / "provenance" / f"{TABLE_PREFIX}_render_manifest.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)
    outputs.append(manifest_path)
    lines = [
        f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(final_dir)}\n"
        for path in outputs
    ]
    (final_dir / "provenance" / f"{TABLE_PREFIX}_sha256_manifest.txt").write_text("".join(lines))


def main() -> None:
    args = parse_args()
    project_root = args.project_root.resolve()
    final_dir = args.final_dir.resolve()
    for relative in ["figures/png", "figures/pdf", "figures/svg", "tables", "code", "provenance"]:
        (final_dir / relative).mkdir(parents=True, exist_ok=True)
    full, neuronal, extraction_qc = load_data(project_root, final_dir)
    cluster_summary = summarize_groups(neuronal, "cluster", CLUSTER_LABELS)
    lineage_summary = summarize_groups(neuronal, "lineage", LINEAGE_LABELS)
    pair_qc = completeness(lineage_summary)
    if not pair_qc["complete_pair"].all():
        raise ValueError(f"Incomplete organoid lineage pairs:\n{pair_qc.to_string(index=False)}")
    rng = np.random.default_rng(args.random_state)
    endpoint_specs = [
        (gene, f"mean_{gene}_log1p_CP10K", "Mean expression", "primary")
        for gene in GENES
    ] + [
        spec
        for gene in GENES
        for spec in [
            (gene, f"median_{gene}_log1p_CP10K", "Median expression", "sensitivity"),
            (gene, f"pct_{gene}_detected", f"{gene}-detected cells (%)", "sensitivity"),
        ]
    ]
    sensitivity = pd.DataFrame(
        [
            analyze_endpoint(
                lineage_summary, gene, endpoint, label, role, rng, args.bootstrap_iterations
            )
            for gene, endpoint, label, role in endpoint_specs
        ]
    )
    primary = sensitivity.loc[sensitivity["analysis_role"] == "primary"].copy()
    equal_value_frames = []
    equal_statistic_frames = []
    for gene in GENES:
        values, statistics = equal_cluster_weight_sensitivity(
            cluster_summary, gene, rng, args.bootstrap_iterations
        )
        equal_value_frames.append(values)
        equal_statistic_frames.append(statistics)
    equal_values = equal_value_frames[0]
    for frame in equal_value_frames[1:]:
        equal_values = equal_values.merge(frame, on="organoid", validate="one_to_one")
    equal_statistics = pd.concat(equal_statistic_frames, ignore_index=True)
    ci = group_mean_ci(lineage_summary)
    tables = final_dir / "tables"
    table_outputs = []
    table_data = {
        f"{TABLE_PREFIX}_extraction_qc.tsv": extraction_qc,
        f"{TABLE_PREFIX}_neuronal_per_cell.tsv.gz": neuronal[
            [
                "cell_id", "sample", "cluster", "paper_cluster_annotation",
                "publication_cluster", "publication_cluster_annotation", "lineage",
                "umap_1", "umap_2", *GENES,
            ]
        ],
        f"{TABLE_PREFIX}_organoid_cluster_aggregates.tsv": cluster_summary,
        f"{TABLE_PREFIX}_organoid_lineage_aggregates.tsv": lineage_summary,
        f"{TABLE_PREFIX}_pair_completeness_qc.tsv": pair_qc,
        f"{TABLE_PREFIX}_group_mean_95CI_summary.tsv": ci,
        f"{TABLE_PREFIX}_sensitivity_endpoints.tsv": sensitivity,
        f"{TABLE_PREFIX}_primary_paired_wilcoxon.tsv": sensitivity.loc[
            sensitivity["analysis_role"] == "primary"
        ],
        f"{TABLE_PREFIX}_equal_cluster_weight_values.tsv": equal_values,
        f"{TABLE_PREFIX}_equal_cluster_weight_statistics.tsv": equal_statistics,
    }
    for name, frame in table_data.items():
        path = tables / name
        frame.to_csv(path, sep="\t", index=False)
        table_outputs.append(path)
    expression_outputs, vmax = render_expression_figure(
        full, neuronal, cluster_summary, lineage_summary, final_dir, args.random_state
    )
    stats_outputs = render_stats_figure(lineage_summary, ci, primary, final_dir)
    documentation = write_documentation(
        final_dir, primary, sensitivity, vmax, equal_statistics
    )
    script_source = Path(__file__).resolve()
    copied_script = final_dir / "code" / SCRIPT_NAME
    shutil.copy2(script_source, copied_script)
    outputs = expression_outputs + stats_outputs + table_outputs + documentation + [copied_script]
    write_manifest(final_dir, outputs, project_root)
    print(f"Output folder: {final_dir}")
    print(f"NRP1/NRP2 neuronal cells: {len(neuronal)}")
    print(f"Complete paired organoids: {int(pair_qc['complete_pair'].sum())}")
    for row in primary.itertuples(index=False):
        print(
            f"Cortical vs subpallial {row.gene}: "
            f"mean difference={row.mean_paired_difference:.6f}, "
            f"P={row.p_value_raw_exact_two_sided:.6f}, "
            f"r_rb={row.rank_biserial_correlation:.3f}, "
            f"cortical greater={int(row.n_first_greater_than_second)}/"
            f"{int(row.n_paired_organoids)}"
        )
    for output in expression_outputs + stats_outputs:
        print(output)


if __name__ == "__main__":
    main()
