#!/usr/bin/env python3
"""Render a lineage-combined counterpart to the DIV30 LHX6/LHX8/PLXNA2 figure.

The canonical all-cluster context is retained. Within the neuronal analysis,
clusters 1 and 4 are displayed and summarized as the cortical interneuron
lineage, while cluster 2 is the subpallial MGE lineage. Existing source tables,
UMAP coordinates, feature scales, z-score definitions, and stabilized ratio
definitions are reused without rewriting the original workflow or outputs.
"""

from __future__ import annotations

import argparse
from datetime import datetime
import hashlib
from pathlib import Path
import shutil
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

REPOSITORY_ROOT = Path("/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline")
LOCAL_PYTHON_SRC = REPOSITORY_ROOT / "python_notebooks" / "src"
if str(LOCAL_PYTHON_SRC) not in sys.path:
    sys.path.insert(0, str(LOCAL_PYTHON_SRC))

from render_varela_div30_lhx6_lhx8_plxna2_final_figure import (
    CLUSTER_COLORS,
    GENES,
    PUBLICATION_CLUSTER_BY_SEURAT,
    PUBLICATION_CLUSTER_COLORS,
    PUBLICATION_CLUSTER_LABELS,
    RATIO_COL,
    RATIO_LABEL,
    SCORE_COL,
    SCORE_LABEL,
    TEXT_COLOR,
    common_limits,
    deterministic_sample,
    make_umap_and_colorbar_axes,
    make_umap_and_key_axes,
    panel_label,
    plot_feature_umap,
    plot_publication_cluster_key_umap,
    save_figure,
    style_axis,
    style_umap_axis,
    umap_bounds,
)


matplotlib.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Nimbus Sans", "Liberation Sans", "DejaVu Sans"],
        "font.size": 7,
        "axes.linewidth": 0.65,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }
)

FINAL_DIR_DEFAULT = Path(
    "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/final_figures/"
    "fig_varela_div30_lhx6_lhx8_plxna2_v1_candidate"
)
FULL_TABLE_RELATIVE = Path("tables/varela_div30_lhx6_lhx8_plxna2_per_cell.tsv.gz")
NEURONAL_TABLE_RELATIVE = Path("tables/neuronal_clusters_1_2_4_per_cell.tsv.gz")
MAPPING_TABLE_RELATIVE = Path("tables/cluster_annotation_mapping.tsv")
SCRIPT_NAME = "render_varela_div30_lhx6_lhx8_plxna2_lineage_orthogonal_figure.py"
FIGURE_STEM = "varela_div30_lhx6_lhx8_lineages_plxna2_orthogonal"
TABLE_PREFIX = "lineage_orthogonal_lhx6_lhx8_plxna2"

LINEAGE_ORDER = ["cortical_interneuron_lineage", "subpallial_MGE_lineage"]
LINEAGE_LABELS = {
    "cortical_interneuron_lineage": "Cortical interneuron lineage",
    "subpallial_MGE_lineage": "Subpallial MGE neurons",
}
LINEAGE_SOURCE_LABELS = {
    "cortical_interneuron_lineage": "publication clusters 3 + 4",
    "subpallial_MGE_lineage": "publication cluster 5",
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
    parser.add_argument("--final-dir", type=Path, default=FINAL_DIR_DEFAULT)
    parser.add_argument("--full-table", type=Path, default=None)
    parser.add_argument("--neuronal-table", type=Path, default=None)
    parser.add_argument("--cluster-mapping", type=Path, default=None)
    parser.add_argument("--random-state", type=int, default=4217)
    return parser.parse_args()


def load_inputs(
    full_path: Path, neuronal_path: Path, mapping_path: Path
) -> tuple[pd.DataFrame, pd.DataFrame]:
    full = pd.read_csv(full_path, sep="\t", dtype={"cluster": str})
    neuronal = pd.read_csv(neuronal_path, sep="\t", dtype={"cluster": str})
    mapping = pd.read_csv(mapping_path, sep="\t", dtype={"seurat_clusters": str})
    full_required = {"cell_id", "sample", "cluster", "umap_1", "umap_2", *GENES}
    neuronal_required = {
        "cell_id",
        "sample",
        "cluster",
        "paper_cluster_annotation",
        "umap_1",
        "umap_2",
        *GENES,
        SCORE_COL,
        RATIO_COL,
    }
    missing_full = sorted(full_required.difference(full.columns))
    missing_neuronal = sorted(neuronal_required.difference(neuronal.columns))
    if missing_full or missing_neuronal:
        raise ValueError(
            f"Missing source columns; full={missing_full}, neuronal={missing_neuronal}"
        )
    name_map = mapping.set_index("seurat_clusters")["paper_cluster_annotation"]
    full["paper_cluster_annotation"] = full["cluster"].map(name_map)
    if full["paper_cluster_annotation"].isna().any():
        missing = sorted(full.loc[full["paper_cluster_annotation"].isna(), "cluster"].unique())
        raise ValueError(f"Full table contains unmapped clusters: {missing}")
    full["publication_cluster"] = full["cluster"].map(PUBLICATION_CLUSTER_BY_SEURAT)
    full["publication_cluster_annotation"] = full["publication_cluster"].map(
        PUBLICATION_CLUSTER_LABELS
    )
    neuronal["publication_cluster"] = neuronal["cluster"].map(
        PUBLICATION_CLUSTER_BY_SEURAT
    )
    neuronal["publication_cluster_annotation"] = neuronal["publication_cluster"].map(
        PUBLICATION_CLUSTER_LABELS
    )
    neuronal["lineage"] = neuronal["cluster"].map(CLUSTER_TO_LINEAGE)
    if neuronal["lineage"].isna().any():
        raise ValueError("Neuronal table contains cells outside clusters 1, 2, and 4.")
    expected = {
        "1": "SST+ cIN",
        "4": "PV neuron precursor",
        "2": "MGE subpallial neurons",
    }
    observed = (
        neuronal[["cluster", "paper_cluster_annotation"]]
        .drop_duplicates()
        .set_index("cluster")["paper_cluster_annotation"]
        .to_dict()
    )
    if observed != expected:
        raise ValueError(f"Canonical neuronal labels changed: observed={observed}")
    return full, neuronal


def plot_lineage_key_umap(
    ax: plt.Axes,
    key_ax: plt.Axes,
    neuronal: pd.DataFrame,
    bounds: tuple[tuple[float, float], tuple[float, float]],
) -> None:
    for lineage in LINEAGE_ORDER:
        subset = neuronal.loc[neuronal["lineage"] == lineage]
        ax.scatter(
            subset["umap_1"],
            subset["umap_2"],
            s=0.55,
            color=LINEAGE_COLORS[lineage],
            alpha=0.78,
            linewidths=0,
            rasterized=True,
        )
    style_umap_axis(ax, bounds)
    ax.set_title("DIV30 MGEO lineages", fontsize=8.2, fontweight="bold", pad=3)
    key_ax.set_axis_off()
    display_numbers = {
        "cortical_interneuron_lineage": "3 + 4",
        "subpallial_MGE_lineage": "5",
    }
    for row, lineage in enumerate(LINEAGE_ORDER):
        key_ax.text(
            0.02,
            0.82 - row * 0.52,
            f"{display_numbers[lineage]}  {LINEAGE_LABELS[lineage]}",
            transform=key_ax.transAxes,
            color=LINEAGE_COLORS[lineage],
            fontsize=6.3,
            fontweight="bold",
            ha="left",
            va="center",
        )


def lineage_organoid_expression_summary(neuronal: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (sample, lineage), group in neuronal.groupby(["sample", "lineage"], sort=True):
        for gene in GENES:
            values = group[gene].to_numpy(dtype=float)
            q25, q75 = np.quantile(values, [0.25, 0.75])
            rows.append(
                {
                    "sample": sample,
                    "sample_display": sample.replace("9583-", ""),
                    "lineage": lineage,
                    "lineage_label": LINEAGE_LABELS[lineage],
                    "source_clusters": LINEAGE_SOURCE_LABELS[lineage],
                    "gene": gene,
                    "n_cells": len(values),
                    "mean_log1p_cp10k": float(np.mean(values)),
                    "median_log1p_cp10k": float(np.median(values)),
                    "q25_log1p_cp10k": float(q25),
                    "q75_log1p_cp10k": float(q75),
                    "pct_detected": float(100 * np.mean(values > 0)),
                }
            )
    return pd.DataFrame(rows)


def plot_lineage_expression_panel(
    ax: plt.Axes,
    neuronal: pd.DataFrame,
    organoid_summary: pd.DataFrame,
    gene: str,
    y_max: float,
    random_state: int,
) -> None:
    values = [
        neuronal.loc[neuronal["lineage"] == lineage, gene].to_numpy(dtype=float)
        for lineage in LINEAGE_ORDER
    ]
    rng = np.random.default_rng(random_state)
    draw_values = [
        rng.choice(value, size=min(len(value), 12_000), replace=False) for value in values
    ]
    positions = np.asarray([0.0, 1.0])
    violins = ax.violinplot(draw_values, positions=positions, widths=0.72, showextrema=False)
    for body, lineage in zip(violins["bodies"], LINEAGE_ORDER, strict=True):
        body.set_facecolor(LINEAGE_COLORS[lineage])
        body.set_edgecolor("none")
        body.set_alpha(0.58)
    for xpos, lineage, value in zip(positions, LINEAGE_ORDER, values, strict=True):
        q25, median, q75 = np.quantile(value, [0.25, 0.5, 0.75])
        ax.plot(
            [xpos, xpos], [q25, q75], color="white", linewidth=2.5,
            solid_capstyle="round", zorder=3,
        )
        ax.scatter(
            [xpos], [median], s=13, color=TEXT_COLOR, edgecolor="white",
            linewidth=0.35, zorder=4,
        )
        means = organoid_summary.loc[
            (organoid_summary["lineage"] == lineage)
            & (organoid_summary["gene"] == gene),
            "mean_log1p_cp10k",
        ].to_numpy(dtype=float)
        jitter = np.linspace(-0.13, 0.13, len(means))
        ax.scatter(
            xpos + jitter,
            means,
            s=12,
            facecolor="white",
            edgecolor=TEXT_COLOR,
            linewidth=0.45,
            alpha=0.95,
            zorder=5,
        )
    ax.set_xticks(
        positions,
        ["3+4\nCortical interneuron\nlineage", "5\nSubpallial MGE\nneurons"],
        linespacing=0.95,
    )
    ax.set_xlim(-0.55, 1.55)
    ax.set_ylim(0, y_max)
    ax.set_title(rf"$\it{{{gene}}}$", fontsize=7.4, fontweight="bold", pad=5)
    ax.set_ylabel("Expression, log1p(CP10K)", fontsize=6.3)
    style_axis(ax)


def plot_lineage_relationship_panel(
    ax: plt.Axes,
    neuronal: pd.DataFrame,
    lineage: str,
    predictor: str,
    predictor_label: str,
    predictor_limit: float,
    plx_limit: float,
    random_state: int,
) -> None:
    data = neuronal.loc[neuronal["lineage"] == lineage]
    draw = deterministic_sample(data, 10_000, random_state)
    ax.scatter(
        draw[predictor],
        draw["PLXNA2"],
        s=2.0,
        color=LINEAGE_COLORS[lineage],
        alpha=0.12,
        linewidths=0,
        rasterized=True,
    )
    x = data[predictor].to_numpy(dtype=float)
    y = data["PLXNA2"].to_numpy(dtype=float)
    if len(x) >= 2 and np.std(x) > 0:
        slope, intercept = np.polyfit(x, y, 1)
        line_x = np.linspace(-predictor_limit, predictor_limit, 80)
        ax.plot(
            line_x,
            slope * line_x + intercept,
            color=LINEAGE_COLORS[lineage],
            linewidth=1.0,
        )
    ax.set_title(
        f"{LINEAGE_SOURCE_LABELS[lineage]}  {LINEAGE_LABELS[lineage]}",
        fontsize=6.8,
        fontweight="bold",
        color=LINEAGE_COLORS[lineage],
        pad=5,
    )
    ax.set_xlim(-predictor_limit, predictor_limit)
    ax.set_ylim(0, plx_limit)
    ax.set_xlabel(predictor_label, fontsize=6.2)
    ax.set_ylabel(r"$\it{PLXNA2}$, log1p(CP10K)", fontsize=6.2)
    style_axis(ax)


def relationship_summary(neuronal: pd.DataFrame, predictor: str) -> pd.DataFrame:
    rows = []
    for sample in ["all", *sorted(neuronal["sample"].unique())]:
        sample_data = neuronal if sample == "all" else neuronal.loc[neuronal["sample"] == sample]
        for lineage in LINEAGE_ORDER:
            data = sample_data.loc[sample_data["lineage"] == lineage]
            x = data[predictor].to_numpy(dtype=float)
            y = data["PLXNA2"].to_numpy(dtype=float)
            finite = np.isfinite(x) & np.isfinite(y)
            x, y = x[finite], y[finite]
            if len(x) >= 2 and np.std(x) > 0:
                slope, intercept = np.polyfit(x, y, 1)
            else:
                slope, intercept = np.nan, np.nan
            rows.append(
                {
                    "sample": sample,
                    "sample_display": "all" if sample == "all" else sample.replace("9583-", ""),
                    "lineage": lineage,
                    "lineage_label": LINEAGE_LABELS[lineage],
                    "source_clusters": LINEAGE_SOURCE_LABELS[lineage],
                    "predictor": predictor,
                    "n_cells": len(x),
                    "pearson_r_descriptive": (
                        float(np.corrcoef(x, y)[0, 1])
                        if len(x) >= 3 and np.std(x) > 0 and np.std(y) > 0
                        else np.nan
                    ),
                    "spearman_rho_descriptive": (
                        float(stats.spearmanr(x, y).statistic) if len(x) >= 3 else np.nan
                    ),
                    "linear_slope_descriptive": float(slope),
                    "linear_intercept_descriptive": float(intercept),
                    "p_values_computed": False,
                }
            )
    return pd.DataFrame(rows)


def render_primary(
    full: pd.DataFrame,
    neuronal: pd.DataFrame,
    organoid_summary: pd.DataFrame,
    final_dir: Path,
    random_state: int,
) -> list[Path]:
    expression_vmax, score_vmax, ratio_vmax, expression_ymax = common_limits(neuronal)
    full_bounds = umap_bounds(full)
    neuronal_bounds = umap_bounds(neuronal)
    fig = plt.figure(figsize=(7.2, 6.65), facecolor="white")
    outer = fig.add_gridspec(
        3,
        1,
        left=0.07,
        right=0.985,
        bottom=0.07,
        top=0.97,
        height_ratios=[1.0, 1.0, 1.14],
        hspace=0.55,
    )

    top = outer[0].subgridspec(1, 4, width_ratios=[1.55, 1, 1, 1], wspace=0.30)
    ax_a, key_a = make_umap_and_key_axes(top[0, 0], fig, key_height=1.9)
    plot_publication_cluster_key_umap(
        ax_a,
        key_a,
        full,
        ["1", "2", "3", "4", "5"],
        full_bounds,
        "DIV30 MGEOs",
        key_fontsize=6.3,
    )
    panel_label(ax_a, "a")
    for idx, (gene, label) in enumerate(zip(GENES, ["b", "c", "d"], strict=True), start=1):
        ax, cax = make_umap_and_colorbar_axes(top[0, idx], fig, bar_height=0.52)
        plot_feature_umap(
            fig, ax, cax, full, gene, full_bounds,
            expression_vmax, score_vmax, ratio_vmax,
        )
        panel_label(ax, label)

    middle = outer[1].subgridspec(1, 4, width_ratios=[1.55, 1, 1, 1], wspace=0.30)
    ax_e, key_e = make_umap_and_key_axes(middle[0, 0], fig, key_height=1.55)
    plot_lineage_key_umap(ax_e, key_e, neuronal, neuronal_bounds)
    panel_label(ax_e, "e")
    for idx, (feature, label) in enumerate(
        zip(GENES, ["f", "g", "h"], strict=True), start=1
    ):
        ax, cax = make_umap_and_colorbar_axes(middle[0, idx], fig, bar_height=0.52)
        plot_feature_umap(
            fig,
            ax,
            cax,
            neuronal,
            feature,
            neuronal_bounds,
            expression_vmax,
            score_vmax,
            ratio_vmax,
        )
        panel_label(ax, label)

    distributions = outer[2].subgridspec(1, 3, wspace=0.33)
    for idx, (gene, label) in enumerate(zip(GENES, ["i", "j", "k"], strict=True)):
        ax = fig.add_subplot(distributions[0, idx])
        plot_lineage_expression_panel(
            ax, neuronal, organoid_summary, gene, expression_ymax, random_state + idx
        )
        panel_label(ax, label)

    return save_figure(fig, final_dir, FIGURE_STEM)


def write_caption(final_dir: Path, expression_vmax: float, score_vmax: float, ratio_vmax: float) -> None:
    caption = f"""Lineage-combined orthogonal view of LHX6, LHX8, and PLXNA2 in Varela DIV30 organoids.

(a) All Varela DIV30 cells on the original Seurat UMAP, colored and labeled as the five DIV30 MGEO populations. Original Seurat clusters 0, 3, and 7 are combined visually as population 1, Radial Glia; Seurat cluster 6 maps to population 2, Inhibitory Progenitors; Seurat clusters 1, 4, and 2 map to populations 3 (SST+ cIN), 4 (PV Neuron Precursor), and 5 (MGE Subpallial Neuron), respectively. (b-d) LHX6, LHX8, and PLXNA2 on the same full-dataset UMAP. (e) Neuronal subset on the unchanged UMAP coordinates, with populations 3 and 4 combined only for display as the cortical interneuron lineage and population 5 shown as the subpallial MGE lineage. (f-h) LHX6, LHX8, and PLXNA2 within this neuronal subset, with the same shared 0 to {expression_vmax:.1f} log1p(CP10K) feature range used by the population figure. (i-k) Descriptive cell-level expression distributions after pooling populations 3 and 4 as the cortical interneuron lineage; the six organoid-level lineage means remain visible as individual points. These cells are not independent biological replicates, and no cell-level inferential tests are shown. Score, ratio, and relationship analyses remain in the exported tables and preserved companion workflows but are intentionally omitted from this compact primary figure.

Original Seurat cluster identities remain recorded in the source tables, while visual labels follow the publication mapping. Cluster identities are not redefined by LHX6, LHX8, or PLXNA2. The lineage pooling is an analysis grouping only. Replicate-aware inference for cortical interneuron lineage versus subpallial MGE neurons is provided separately in `varela_div30_plxna2_cortical_vs_subpallial_replicate_aware` and uses the organoid as the biological replicate.
"""
    (final_dir / "lineage_orthogonal_lhx6_lhx8_plxna2_figure_caption.txt").write_text(caption)


def git_value(args: list[str]) -> str:
    try:
        return subprocess.run(args, check=True, capture_output=True, text=True).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unavailable"


def write_provenance(
    final_dir: Path,
    full_path: Path,
    neuronal_path: Path,
    mapping_path: Path,
    script_source: Path,
    expression_vmax: float,
    score_vmax: float,
    ratio_vmax: float,
) -> None:
    manifest = pd.DataFrame(
        [
            ("rendered_at_local", datetime.now().isoformat(timespec="seconds")),
            ("full_source_table", str(full_path.resolve())),
            ("neuronal_source_table", str(neuronal_path.resolve())),
            ("canonical_mapping_table", str(mapping_path.resolve())),
            ("publication_mapping", "Seurat 0/3/7->pub1; 6->pub2; 1->pub3; 4->pub4; 2->pub5"),
            ("cortical_lineage_definition", "publication clusters 3 and 4 pooled"),
            ("subpallial_lineage_definition", "publication cluster 5"),
            ("expression_vmax_log1p_CP10K", f"{expression_vmax:.9g}"),
            ("zscore_display_vmax", f"{score_vmax:.9g}"),
            ("ratio_display_vmax", f"{ratio_vmax:.9g}"),
            ("inferential_statistics", "none; descriptive figure"),
            ("script", str(script_source.resolve())),
            ("git_commit", git_value(["git", "rev-parse", "HEAD"])),
        ],
        columns=["key", "value"],
    )
    manifest_path = final_dir / "provenance" / f"{TABLE_PREFIX}_render_manifest.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)
    files = [
        final_dir / "figures" / extension / f"{FIGURE_STEM}.{extension}"
        for extension in ["png", "pdf", "svg"]
    ]
    files += [
        final_dir / "tables" / f"{TABLE_PREFIX}_organoid_expression_summary.tsv",
        final_dir / "tables" / f"{TABLE_PREFIX}_score_relationship_summary.tsv",
        final_dir / "tables" / f"{TABLE_PREFIX}_ratio_relationship_summary.tsv",
        final_dir / "tables" / f"{TABLE_PREFIX}_cluster_to_lineage_mapping.tsv",
        final_dir / "lineage_orthogonal_lhx6_lhx8_plxna2_figure_caption.txt",
        final_dir / "code" / SCRIPT_NAME,
        manifest_path,
    ]
    lines = []
    for path in files:
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        lines.append(f"{digest}  {path.relative_to(final_dir)}\n")
    (final_dir / "provenance" / f"{TABLE_PREFIX}_sha256_manifest.txt").write_text(
        "".join(lines)
    )


def main() -> None:
    args = parse_args()
    final_dir = args.final_dir.resolve()
    full_path = (args.full_table or final_dir / FULL_TABLE_RELATIVE).resolve()
    neuronal_path = (args.neuronal_table or final_dir / NEURONAL_TABLE_RELATIVE).resolve()
    mapping_path = (args.cluster_mapping or final_dir / MAPPING_TABLE_RELATIVE).resolve()
    script_source = Path(__file__).resolve()
    for relative in ["figures/png", "figures/pdf", "figures/svg", "tables", "code", "provenance"]:
        (final_dir / relative).mkdir(parents=True, exist_ok=True)

    full, neuronal = load_inputs(full_path, neuronal_path, mapping_path)
    organoid_summary = lineage_organoid_expression_summary(neuronal)
    score_relationship = relationship_summary(neuronal, SCORE_COL)
    ratio_relationship = relationship_summary(neuronal, RATIO_COL)
    expression_vmax, score_vmax, ratio_vmax, _ = common_limits(neuronal)

    organoid_summary.to_csv(
        final_dir / "tables" / f"{TABLE_PREFIX}_organoid_expression_summary.tsv",
        sep="\t",
        index=False,
    )
    score_relationship.to_csv(
        final_dir / "tables" / f"{TABLE_PREFIX}_score_relationship_summary.tsv",
        sep="\t",
        index=False,
    )
    ratio_relationship.to_csv(
        final_dir / "tables" / f"{TABLE_PREFIX}_ratio_relationship_summary.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(
        [
            {
                "cluster": cluster,
                "canonical_cluster_label": label,
                "publication_cluster": PUBLICATION_CLUSTER_BY_SEURAT[cluster],
                "publication_cluster_label": PUBLICATION_CLUSTER_LABELS[
                    PUBLICATION_CLUSTER_BY_SEURAT[cluster]
                ],
                "lineage": CLUSTER_TO_LINEAGE[cluster],
                "lineage_label": LINEAGE_LABELS[CLUSTER_TO_LINEAGE[cluster]],
            }
            for cluster, label in [
                ("1", "SST+ cIN"),
                ("4", "PV neuron precursor"),
                ("2", "MGE subpallial neurons"),
            ]
        ]
    ).to_csv(
        final_dir / "tables" / f"{TABLE_PREFIX}_cluster_to_lineage_mapping.tsv",
        sep="\t",
        index=False,
    )

    outputs = render_primary(full, neuronal, organoid_summary, final_dir, args.random_state)
    write_caption(final_dir, expression_vmax, score_vmax, ratio_vmax)
    shutil.copy2(script_source, final_dir / "code" / SCRIPT_NAME)
    write_provenance(
        final_dir,
        full_path,
        neuronal_path,
        mapping_path,
        script_source,
        expression_vmax,
        score_vmax,
        ratio_vmax,
    )

    print(f"Output folder: {final_dir}")
    print("Original multipanel figure and workflows: preserved")
    print(f"Lineage-combined cells: {len(neuronal)}")
    print(f"Organoid-by-lineage-by-gene summaries: {len(organoid_summary)}")
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
