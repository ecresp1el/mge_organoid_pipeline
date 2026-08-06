#!/usr/bin/env python3
"""Render DIV30 guidance-gene feature UMAPs and replicate-resolved co-expression maps."""

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
    _apply_internal_umap_plot_filters,
    extract_marker_expression_from_h5ad,
)
from render_varela_div30_lhx6_lhx8_plxna2_final_figure import (
    add_cluster_annotations,
    div30_spec,
    make_umap_and_colorbar_axes,
    panel_label,
    plot_feature_umap,
    round_scale_limit,
    umap_bounds,
)


PROJECT_ROOT_DEFAULT = Path("/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder")
FINAL_DIR_DEFAULT = PROJECT_ROOT_DEFAULT / "final_figures" / "fig_varela_div30_lhx6_lhx8_plxna2_v1_candidate"
SCRIPT_NAME = "render_varela_div30_guidance_gene_feature_coexpression.py"
TABLE_PREFIX = "guidance_gene_feature_coexpression"
FULL_STEM = "varela_div30_guidance_genes_feature_umaps_all_cells"
NEURONAL_STEM = "varela_div30_guidance_genes_feature_umaps_neuronal"
COEXPRESSION_STEM = "varela_div30_guidance_genes_coexpression_by_lineage"
GENES = [
    "NRP1", "NRP2", "PLXNA2", "LHX6", "LHX8", "NKX2-1",
    "SOX6", "ROBO1", "ROBO2", "SEMA3A", "SEMA3C",
]
EXPORT_DPI = 300
NEURONAL_CLUSTERS = ["1", "4", "2"]
CLUSTER_TO_LINEAGE = {
    "1": "cortical_interneuron_lineage",
    "4": "cortical_interneuron_lineage",
    "2": "subpallial_MGE_lineage",
}
STRATA = ["all_neuronal", "cortical_interneuron_lineage", "subpallial_MGE_lineage"]
STRATUM_LABELS = {
    "all_neuronal": "All neuronal populations",
    "cortical_interneuron_lineage": "Cortical interneuron lineage",
    "subpallial_MGE_lineage": "Subpallial MGE neurons",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT_DEFAULT)
    parser.add_argument("--final-dir", type=Path, default=FINAL_DIR_DEFAULT)
    return parser.parse_args()


def load_data(
    project_root: Path, final_dir: Path
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    raw_path = final_dir / "tables" / f"{TABLE_PREFIX}_all_cells.tsv.gz"
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
    if neuronal[GENES].isna().any().any():
        raise ValueError("Missing expression values were found after gene extraction.")
    return annotated, neuronal, matches, filter_summary


def common_expression_limit(full: pd.DataFrame) -> float:
    q99 = [float(np.quantile(full[gene], 0.99)) for gene in GENES]
    return round_scale_limit(max(EXPRESSION_COLOR_FLOOR, max(q99)))


def save_figure(fig: plt.Figure, final_dir: Path, stem: str) -> list[Path]:
    outputs = []
    for extension in ["png", "pdf", "svg"]:
        output = final_dir / "figures" / extension / f"{stem}.{extension}"
        fig.savefig(output, dpi=EXPORT_DPI, bbox_inches="tight", facecolor="white")
        outputs.append(output)
    plt.close(fig)
    return outputs


def render_feature_atlas(
    data: pd.DataFrame,
    final_dir: Path,
    stem: str,
    figure_title: str,
    vmax: float,
) -> list[Path]:
    bounds = umap_bounds(data)
    fig = plt.figure(figsize=(8.2, 5.9), facecolor="white")
    grid = fig.add_gridspec(
        3, 4, left=0.05, right=0.985, bottom=0.065, top=0.92,
        wspace=0.27, hspace=0.36,
    )
    for index, gene in enumerate(GENES):
        row, col = divmod(index, 4)
        ax, cax = make_umap_and_colorbar_axes(grid[row, col], fig, bar_height=0.52)
        title = r"$\it{NKX2}$-1" if gene == "NKX2-1" else None
        plot_feature_umap(
            fig, ax, cax, data, gene, bounds, vmax, 1.0, 1.0, title=title
        )
        panel_label(ax, chr(ord("a") + index))
    for index in range(len(GENES), 12):
        row, col = divmod(index, 4)
        ax = fig.add_subplot(grid[row, col])
        ax.set_axis_off()
    fig.suptitle(figure_title, fontsize=9.2, fontweight="bold", y=0.975)
    return save_figure(fig, final_dir, stem)


def stratum_frame(neuronal: pd.DataFrame, stratum: str) -> pd.DataFrame:
    if stratum == "all_neuronal":
        return neuronal
    return neuronal.loc[neuronal["lineage"] == stratum]


def positive_high_mask(values: np.ndarray) -> tuple[np.ndarray, float]:
    detected = values[values > 0]
    if detected.size == 0:
        return np.zeros(values.size, dtype=bool), np.nan
    threshold = float(np.quantile(detected, 0.75))
    return values >= threshold, threshold


def calculate_coexpression(neuronal: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    organoids = sorted(neuronal["organoid"].unique())
    for stratum in STRATA:
        subset = stratum_frame(neuronal, stratum)
        for organoid in organoids:
            frame = subset.loc[subset["organoid"] == organoid]
            if frame.empty:
                continue
            for gene_a, gene_b in itertools.combinations_with_replacement(GENES, 2):
                a = frame[gene_a].to_numpy(dtype=float)
                b = frame[gene_b].to_numpy(dtype=float)
                rho = float(stats.spearmanr(a, b).statistic)
                high_a, threshold_a = positive_high_mask(a)
                high_b, threshold_b = positive_high_mask(b)
                marginal_a = float(np.mean(high_a))
                marginal_b = float(np.mean(high_b))
                observed = float(np.mean(high_a & high_b))
                expected = marginal_a * marginal_b
                enrichment = observed / expected if expected > 0 else np.nan
                rows.append(
                    {
                        "stratum": stratum,
                        "stratum_label": STRATUM_LABELS[stratum],
                        "organoid": organoid,
                        "gene_a": gene_a,
                        "gene_b": gene_b,
                        "n_cells": len(frame),
                        "spearman_rho_descriptive": rho,
                        "gene_a_detected_pct": 100 * float(np.mean(a > 0)),
                        "gene_b_detected_pct": 100 * float(np.mean(b > 0)),
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
                        "high_definition": "top quartile among detected cells, calculated within organoid and stratum",
                        "p_values_computed": False,
                    }
                )
    return pd.DataFrame(rows)


def summarize_coexpression(per_organoid: pd.DataFrame) -> pd.DataFrame:
    metrics = [
        "spearman_rho_descriptive",
        "both_high_pct",
        "both_high_expected_pct_independence",
        "high_high_enrichment_observed_over_expected",
        "high_high_log2_enrichment",
    ]
    rows = []
    for keys, frame in per_organoid.groupby(["stratum", "stratum_label", "gene_a", "gene_b"]):
        row = dict(zip(["stratum", "stratum_label", "gene_a", "gene_b"], keys, strict=True))
        row["n_organoids"] = frame["organoid"].nunique()
        for metric in metrics:
            values = frame[metric].to_numpy(dtype=float)
            row[f"median_{metric}"] = float(np.nanmedian(values))
            row[f"q25_{metric}"] = float(np.nanquantile(values, 0.25))
            row[f"q75_{metric}"] = float(np.nanquantile(values, 0.75))
        row["n_organoids_positive_spearman"] = int(
            (frame["spearman_rho_descriptive"] > 0).sum()
        )
        row["n_organoids_high_high_enriched"] = int(
            (frame["high_high_log2_enrichment"] > 0).sum()
        )
        rows.append(row)
    return pd.DataFrame(rows)


def symmetric_matrix(summary: pd.DataFrame, stratum: str, metric: str) -> np.ndarray:
    matrix = np.full((len(GENES), len(GENES)), np.nan)
    indices = {gene: index for index, gene in enumerate(GENES)}
    frame = summary.loc[summary["stratum"] == stratum]
    for row in frame.itertuples(index=False):
        i, j = indices[row.gene_a], indices[row.gene_b]
        value = float(getattr(row, metric))
        matrix[i, j] = value
        matrix[j, i] = value
    return matrix


def style_heatmap_axis(ax: plt.Axes, show_ylabels: bool) -> None:
    ax.set_xticks(np.arange(len(GENES)), GENES, rotation=55, ha="right", fontsize=5.2)
    if show_ylabels:
        ax.set_yticks(np.arange(len(GENES)), GENES, fontsize=5.2)
    else:
        ax.set_yticks(np.arange(len(GENES)), [""] * len(GENES))
    ax.tick_params(length=0)
    plx_index = GENES.index("PLXNA2")
    ax.add_patch(
        Rectangle(
            (-0.5, plx_index - 0.5), len(GENES), 1,
            fill=False, edgecolor="#111111", linewidth=0.75,
        )
    )
    ax.add_patch(
        Rectangle(
            (plx_index - 0.5, -0.5), 1, len(GENES),
            fill=False, edgecolor="#111111", linewidth=0.75,
        )
    )
    for spine in ax.spines.values():
        spine.set_visible(False)


def render_coexpression_figure(
    summary: pd.DataFrame, final_dir: Path
) -> list[Path]:
    rho_metric = "median_spearman_rho_descriptive"
    high_metric = "median_high_high_log2_enrichment"
    rho_values = summary[rho_metric].to_numpy(dtype=float)
    rho_limit = max(0.25, float(np.nanquantile(np.abs(rho_values), 0.98)))
    high_values = summary.loc[summary["gene_a"] != summary["gene_b"], high_metric].to_numpy(dtype=float)
    high_limit = max(0.5, float(np.nanquantile(np.abs(high_values), 0.98)))
    fig, axes = plt.subplots(2, 3, figsize=(8.6, 5.9), facecolor="white")
    rho_image = None
    high_image = None
    for col, stratum in enumerate(STRATA):
        rho_matrix = symmetric_matrix(summary, stratum, rho_metric)
        rho_image = axes[0, col].imshow(
            rho_matrix, cmap="coolwarm", vmin=-rho_limit, vmax=rho_limit,
            interpolation="none", rasterized=True,
        )
        axes[0, col].set_title(STRATUM_LABELS[stratum], fontsize=7.4, fontweight="bold", pad=5)
        style_heatmap_axis(axes[0, col], show_ylabels=col == 0)
        high_matrix = symmetric_matrix(summary, stratum, high_metric)
        np.fill_diagonal(high_matrix, np.nan)
        high_image = axes[1, col].imshow(
            high_matrix, cmap="coolwarm", vmin=-high_limit, vmax=high_limit,
            interpolation="none", rasterized=True,
        )
        style_heatmap_axis(axes[1, col], show_ylabels=col == 0)
        panel_label(axes[0, col], chr(ord("a") + col))
        panel_label(axes[1, col], chr(ord("d") + col))
    axes[0, 0].set_ylabel("Median within-organoid\nSpearman correlation", fontsize=7)
    axes[1, 0].set_ylabel("Median log2 enrichment of\nhigh–high co-expression", fontsize=7)
    cbar_rho = fig.colorbar(rho_image, ax=axes[0, :], orientation="vertical", fraction=0.018, pad=0.018)
    cbar_rho.set_label("Median Spearman ρ", fontsize=6.2)
    cbar_rho.ax.tick_params(labelsize=5.5, length=2)
    cbar_high = fig.colorbar(high_image, ax=axes[1, :], orientation="vertical", fraction=0.018, pad=0.018)
    cbar_high.set_label("Median log2(observed/expected)", fontsize=6.2)
    cbar_high.ax.tick_params(labelsize=5.5, length=2)
    fig.suptitle(
        "DIV30 MGEO guidance-gene co-expression across six organoids",
        fontsize=9.2, fontweight="bold", y=0.985,
    )
    fig.subplots_adjust(left=0.10, right=0.89, bottom=0.10, top=0.91, wspace=0.24, hspace=0.36)
    return save_figure(fig, final_dir, COEXPRESSION_STEM)


def write_plxna2_summary(summary: pd.DataFrame, final_dir: Path) -> Path:
    focused = summary.loc[
        ((summary["gene_a"] == "PLXNA2") | (summary["gene_b"] == "PLXNA2"))
        & (summary["gene_a"] != summary["gene_b"])
    ].copy()
    focused["partner_gene"] = np.where(
        focused["gene_a"] == "PLXNA2", focused["gene_b"], focused["gene_a"]
    )
    focused = focused.sort_values(["stratum", "median_spearman_rho_descriptive"], ascending=[True, False])
    path = final_dir / "tables" / f"{TABLE_PREFIX}_plxna2_focused_summary.tsv"
    focused.to_csv(path, sep="\t", index=False)
    lines = [
        "# PLXNA2-focused co-expression summary",
        "",
        "Values are medians of six within-organoid descriptive summaries; no cell-level inferential P values are used.",
        "High expression is the top quartile among detected cells within each organoid and stratum.",
        "",
        "| Stratum | Partner | Median Spearman ρ | Median high–high log2 enrichment | Organoids with positive ρ |",
        "|---|---|---:|---:|---:|",
    ]
    for row in focused.itertuples(index=False):
        lines.append(
            f"| {row.stratum_label} | {row.partner_gene} | "
            f"{row.median_spearman_rho_descriptive:.3f} | "
            f"{row.median_high_high_log2_enrichment:.3f} | "
            f"{row.n_organoids_positive_spearman}/6 |"
        )
    md_path = final_dir / "guidance_gene_PLXNA2_coexpression_summary.md"
    md_path.write_text("\n".join(lines) + "\n")
    return md_path


def git_value() -> str:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"], check=True, capture_output=True, text=True
        ).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unavailable"


def write_caption(final_dir: Path, vmax: float) -> Path:
    text = f"""DIV30 MGEO guidance-gene feature and co-expression analysis.

The all-cell and neuronal feature atlases show {', '.join(GENES)} on unchanged UMAP coordinates. Every feature panel uses the same 0 to {vmax:g} log1p(CP10K) color range, allowing direct visual comparison. The neuronal atlas is restricted to the established SST+ cIN, PV neuron precursor, and MGE subpallial neuron populations; radial glia and inhibitory progenitors remain excluded without reannotation.

The co-expression figure summarizes cells within each organoid before combining results across the six biological replicates. Panels a-c show the median within-organoid Spearman correlation for all neuronal populations, the cortical interneuron lineage, and subpallial MGE neurons. Panels d-f show the median log2 enrichment of high–high co-expression relative to independence. High expression is defined separately for each gene, organoid, and stratum as the top quartile among cells with detectable expression. Positive values indicate more double-high cells than expected from the two marginal high-expression frequencies; negative values indicate depletion. Matrix diagonals are omitted from the enrichment row. PLXNA2 rows and columns are outlined. These are descriptive co-expression summaries without pooled-cell inferential P values.
"""
    path = final_dir / "guidance_gene_feature_coexpression_figure_caption.txt"
    path.write_text(text)
    return path


def write_manifest(final_dir: Path, project_root: Path, outputs: list[Path], vmax: float) -> None:
    manifest = pd.DataFrame(
        [
            ("rendered_at_local", datetime.now().isoformat(timespec="seconds")),
            ("source_h5ad", str(div30_spec(project_root).resolved_h5ad_path(project_root))),
            ("genes", ",".join(GENES)),
            ("export", "300-dpi PNG; PDF/SVG vector with editable text and rasterized cell layers"),
            ("expression_scale", "existing log1p(CP10K)"),
            ("shared_feature_vmax", f"{vmax:g}"),
            ("coexpression_unit", "within-organoid cells; median across six organoids"),
            ("high_definition", "top quartile among detected cells within organoid and stratum"),
            ("inferential_p_values", "none"),
            ("git_commit", git_value()),
        ],
        columns=["key", "value"],
    )
    path = final_dir / "provenance" / f"{TABLE_PREFIX}_render_manifest.tsv"
    manifest.to_csv(path, sep="\t", index=False)
    outputs.append(path)
    checksum_path = final_dir / "provenance" / f"{TABLE_PREFIX}_sha256_manifest.txt"
    checksum_path.write_text(
        "".join(
            f"{hashlib.sha256(output.read_bytes()).hexdigest()}  {output.relative_to(final_dir)}\n"
            for output in outputs
        )
    )


def main() -> None:
    args = parse_args()
    project_root = args.project_root.resolve()
    final_dir = args.final_dir.resolve()
    for relative in ["figures/png", "figures/pdf", "figures/svg", "tables", "code", "provenance"]:
        (final_dir / relative).mkdir(parents=True, exist_ok=True)
    full, neuronal, matches, filter_summary = load_data(project_root, final_dir)
    vmax = common_expression_limit(full)
    outputs = []
    outputs += render_feature_atlas(
        full, final_dir, FULL_STEM, "DIV30 MGEO guidance-gene expression · all cells", vmax
    )
    outputs += render_feature_atlas(
        neuronal, final_dir, NEURONAL_STEM,
        "DIV30 MGEO guidance-gene expression · neuronal populations", vmax,
    )
    per_organoid = calculate_coexpression(neuronal)
    summary = summarize_coexpression(per_organoid)
    outputs += render_coexpression_figure(summary, final_dir)
    tables = {
        f"{TABLE_PREFIX}_gene_matches.tsv": matches,
        f"{TABLE_PREFIX}_plot_filter_summary.tsv": filter_summary,
        f"{TABLE_PREFIX}_neuronal_per_cell.tsv.gz": neuronal[
            [
                "cell_id", "sample", "cluster", "paper_cluster_annotation",
                "publication_cluster", "publication_cluster_annotation", "lineage",
                "umap_1", "umap_2", *GENES,
            ]
        ],
        f"{TABLE_PREFIX}_per_organoid_pairwise.tsv": per_organoid,
        f"{TABLE_PREFIX}_pairwise_summary.tsv": summary,
    }
    for name, frame in tables.items():
        path = final_dir / "tables" / name
        frame.to_csv(path, sep="\t", index=False)
        outputs.append(path)
    outputs.append(write_plxna2_summary(summary, final_dir))
    focused_table = final_dir / "tables" / f"{TABLE_PREFIX}_plxna2_focused_summary.tsv"
    outputs.append(focused_table)
    outputs.append(write_caption(final_dir, vmax))
    script_source = Path(__file__).resolve()
    script_copy = final_dir / "code" / SCRIPT_NAME
    shutil.copy2(script_source, script_copy)
    outputs.append(script_copy)
    write_manifest(final_dir, project_root, outputs, vmax)
    print(f"Output folder: {final_dir}")
    print(f"Genes: {', '.join(GENES)}")
    print(f"All cells: {len(full)}; neuronal cells: {len(neuronal)}")
    print(f"Shared feature limit: 0 to {vmax:g} log1p(CP10K)")
    print(f"Within-organoid pairwise records: {len(per_organoid)}")
    for stem in [FULL_STEM, NEURONAL_STEM, COEXPRESSION_STEM]:
        print(final_dir / "figures" / "png" / f"{stem}.png")


if __name__ == "__main__":
    main()
