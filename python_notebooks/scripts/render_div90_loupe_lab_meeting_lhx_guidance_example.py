#!/usr/bin/env python3
"""Render one focused DIV90 Loupe LHX/guidance-receptor lab-meeting figure."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
        "font.size": 8,
        "axes.linewidth": 0.65,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd


PROJECT_ROOT_DEFAULT = Path("/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder")
FIGURE_ID = "fig_div90_loupe_recluster_annotations_v1"
TABLE_RELATIVE = Path(
    f"final_figures/{FIGURE_ID}/tables/"
    "div90_loupe_subcluster_guidance_joined_per_cell.tsv.gz"
)
OUTPUT_NAME = "div90_loupe_lab_meeting_lhx_guidance_best_example.png"

CORTICAL_GENES = ["LHX6", "ERBB4", "PLXNA2", "NRP2"]
DISPLAY_GENES = ["LHX6", "LHX8", "ERBB4", "PLXNA2", "NRP2"]
SUBPALLIAL_REFERENCE_LABELS = ["1+3 GP projection", "7+9 striatal?"]

GRAY = "#D8DADD"
GRAY_DARK = "#8D939A"
BLUE_LIGHT = "#9ECAE1"
BLUE = "#08519C"
MAGENTA_LIGHT = "#F7B6D2"
MAGENTA = "#AD1457"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT_DEFAULT)
    parser.add_argument("--input", type=Path, default=None)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--dpi", type=int, default=600)
    return parser.parse_args()


def load_data(path: Path) -> pd.DataFrame:
    data = pd.read_csv(path, sep="\t", dtype={"cell_id": str})
    required = {
        "cell_id",
        "set_id",
        "loupe_label",
        "loupe_x",
        "loupe_y",
        "expression_available",
        *DISPLAY_GENES,
    }
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    data = data.loc[data["expression_available"].astype(bool)].copy()
    for gene in DISPLAY_GENES:
        data[gene] = pd.to_numeric(data[gene], errors="raise")
    return data


def cluster_summary(data: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (set_id, label), frame in data.groupby(["set_id", "loupe_label"], sort=False):
        row: dict[str, object] = {
            "set_id": set_id,
            "loupe_label": label,
            "n_cells": len(frame),
        }
        for gene in DISPLAY_GENES:
            row[f"{gene}_mean"] = float(frame[gene].mean())
            row[f"{gene}_pct"] = float(100 * frame[gene].gt(0).mean())
        rows.append(row)
    return pd.DataFrame(rows)


def select_cortical_cluster(summary: pd.DataFrame) -> str:
    candidates = summary.loc[
        (summary["set_id"] == "cortical_only") & (summary["n_cells"] >= 100)
    ].copy()
    rank_columns = []
    for gene in CORTICAL_GENES:
        column = f"{gene}_rank"
        candidates[column] = candidates[f"{gene}_mean"].rank(pct=True)
        rank_columns.append(column)
    candidates["joint_rank"] = candidates[rank_columns].mean(axis=1)
    return str(candidates.sort_values(["joint_rank", "n_cells"], ascending=False).iloc[0]["loupe_label"])


def select_subpallial_clusters(data: pd.DataFrame, summary: pd.DataFrame) -> list[str]:
    reference = data.loc[
        (data["set_id"] == "subpallial_only")
        & data["loupe_label"].isin(SUBPALLIAL_REFERENCE_LABELS)
    ]
    reference_plxna2 = float(reference["PLXNA2"].mean())
    reference_nrp2 = float(reference["NRP2"].mean())
    candidates = summary.loc[
        (summary["set_id"] == "subpallial_only")
        & (summary["n_cells"] >= 150)
        & (summary["LHX8_pct"] >= 15)
        & (summary["PLXNA2_mean"] < reference_plxna2)
        & (summary["NRP2_mean"] < reference_nrp2)
    ].copy()
    candidates = candidates.sort_values("LHX8_pct", ascending=False)
    return candidates["loupe_label"].astype(str).tolist()


def style_umap(ax: plt.Axes, frame: pd.DataFrame) -> None:
    for axis in ["loupe_x", "loupe_y"]:
        low, high = np.quantile(frame[axis], [0.005, 0.995])
        pad = max(0.1, 0.04 * (high - low))
        getattr(ax, f"set_{'x' if axis == 'loupe_x' else 'y'}lim")(low - pad, high + pad)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def panel_label(ax: plt.Axes, label: str, y: float = 1.13) -> None:
    ax.text(
        -0.10,
        y,
        label,
        transform=ax.transAxes,
        fontsize=11,
        fontweight="bold",
        ha="left",
        va="top",
    )


def plot_selection_umap(
    ax: plt.Axes,
    family: pd.DataFrame,
    parent_mask: pd.Series,
    selected_mask: pd.Series,
    parent_color: str,
    selected_color: str,
    title: str,
    subtitle: str,
) -> None:
    ax.scatter(
        family["loupe_x"],
        family["loupe_y"],
        s=1.1,
        c=GRAY,
        linewidths=0,
        alpha=0.52,
        rasterized=True,
    )
    parent = family.loc[parent_mask]
    ax.scatter(
        parent["loupe_x"],
        parent["loupe_y"],
        s=2.0,
        c=parent_color,
        linewidths=0,
        alpha=0.74,
        rasterized=True,
    )
    selected = family.loc[selected_mask]
    ax.scatter(
        selected["loupe_x"],
        selected["loupe_y"],
        s=7.0,
        c=selected_color,
        edgecolors="white",
        linewidths=0.18,
        alpha=0.96,
        rasterized=True,
    )
    style_umap(ax, family)
    ax.set_title(title, fontsize=10.2, fontweight="bold", pad=18)
    ax.text(
        0.5,
        1.015,
        subtitle,
        transform=ax.transAxes,
        fontsize=7.1,
        color="#43484E",
        ha="center",
        va="bottom",
    )


def expression_row(frame: pd.DataFrame, label: str, group: str) -> dict[str, object]:
    row: dict[str, object] = {"label": label, "group": group, "n_cells": len(frame)}
    for gene in DISPLAY_GENES:
        row[f"{gene}_mean"] = float(frame[gene].mean())
        row[f"{gene}_pct"] = float(100 * frame[gene].gt(0).mean())
    return row


def plot_dotplot(ax: plt.Axes, rows: pd.DataFrame) -> None:
    vmax = 3.0
    for y, row in rows.iterrows():
        for x, gene in enumerate(DISPLAY_GENES):
            pct = float(row[f"{gene}_pct"])
            mean = float(row[f"{gene}_mean"])
            ax.scatter(
                x,
                y,
                s=14 + 105 * pct / 100,
                c=mean,
                cmap="Blues",
                vmin=0,
                vmax=vmax,
                edgecolor="#61676D",
                linewidth=0.35,
                zorder=3,
            )
    ax.axhspan(0.5, 1.5, color="#EAF2F8", zorder=0)
    ax.axhspan(2.5, 3.5, color="#FCEBF3", zorder=0)
    ax.axhline(1.5, color="#AEB3B8", lw=0.7)
    ax.set_xlim(-0.55, len(DISPLAY_GENES) - 0.45)
    ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.set_xticks(range(len(DISPLAY_GENES)), DISPLAY_GENES, rotation=35, ha="right", fontsize=7.4)
    labels = [f"{row.label}\n(n={int(row.n_cells):,})" for row in rows.itertuples(index=False)]
    ax.set_yticks(range(len(rows)), labels, fontsize=7.0)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title("Focused expression summary", fontsize=10.2, fontweight="bold", pad=18)
    ax.text(
        0.5,
        1.015,
        "color = mean expression; size = cells with detected expression",
        transform=ax.transAxes,
        fontsize=7.1,
        color="#43484E",
        ha="center",
        va="bottom",
    )
    color_norm = matplotlib.colors.Normalize(vmin=0, vmax=vmax)
    mapper = matplotlib.cm.ScalarMappable(norm=color_norm, cmap="Blues")
    colorbar = plt.colorbar(mapper, ax=ax, orientation="horizontal", fraction=0.045, pad=0.17)
    colorbar.set_ticks([0, 1.5, 3])
    colorbar.set_label("Mean log1p(CP10K)", fontsize=6.5, labelpad=2)
    colorbar.ax.tick_params(labelsize=6, length=1.5, pad=1)


def plot_cortical_summary(
    ax: plt.Axes,
    cortical_cluster: pd.DataFrame,
    cortical_selected: pd.DataFrame,
    cluster_label: str,
) -> None:
    ax.set_axis_off()
    fraction = 100 * len(cortical_selected) / len(cortical_cluster)
    ax.text(0.0, 0.94, "Best cortical example", fontsize=8.5, fontweight="bold", color=BLUE, va="top")
    ax.text(
        0.0,
        0.72,
        f"{cluster_label}: {len(cortical_cluster):,} cells",
        fontsize=8,
        fontweight="bold",
        va="top",
    )
    ax.text(
        0.0,
        0.49,
        f"{len(cortical_selected):,} cells ({fraction:.1f}%) co-detect\n"
        "LHX6 + ERBB4 + PLXNA2 + NRP2",
        fontsize=8,
        linespacing=1.35,
        va="top",
    )
    ax.add_line(Line2D([0.0, fraction / 100], [0.10, 0.10], transform=ax.transAxes, lw=7, color=BLUE))
    ax.add_line(Line2D([fraction / 100, 1.0], [0.10, 0.10], transform=ax.transAxes, lw=7, color="#E7EAED"))
    ax.text(fraction / 100, 0.0, f"{fraction:.1f}%", transform=ax.transAxes, fontsize=6.8, ha="center")


def plot_subpallial_summary(
    ax: plt.Axes,
    aggregate: pd.DataFrame,
    selected: pd.DataFrame,
    labels: list[str],
) -> None:
    ax.set_axis_off()
    lhx8_positive = aggregate.loc[aggregate["LHX8"] > 0]
    fraction_all = 100 * len(selected) / len(aggregate)
    fraction_lhx8 = 100 * len(selected) / len(lhx8_positive)
    contributions = selected.groupby("loupe_label").size().sort_values(ascending=False)
    contribution_text = "; ".join(f"{label}: {int(count)}" for label, count in contributions.items())
    ax.text(0.0, 0.94, "Closest subpallial example", fontsize=8.5, fontweight="bold", color=MAGENTA, va="top")
    ax.text(
        0.0,
        0.72,
        f"Aggregate of {len(labels)} LHX8-biased, lower-guidance subclusters: {len(aggregate):,} cells",
        fontsize=7.6,
        fontweight="bold",
        va="top",
    )
    ax.text(
        0.0,
        0.48,
        f"{len(selected):,} LHX8+ / PLXNA2− / NRP2− cells\n"
        f"({fraction_all:.1f}% of aggregate; {fraction_lhx8:.1f}% of LHX8+ cells)",
        fontsize=8,
        linespacing=1.35,
        va="top",
    )
    ax.text(0.0, 0.10, contribution_text, fontsize=6.4, color="#555B61", va="bottom", wrap=True)


def draw_violin(
    ax: plt.Axes,
    values: pd.Series,
    position: float,
    width: float,
    color: str,
) -> None:
    array = values.to_numpy(dtype=float)
    parts = ax.violinplot(
        [array],
        positions=[position],
        widths=width,
        showmeans=False,
        showmedians=False,
        showextrema=False,
        points=150,
        bw_method=0.22,
    )
    for body in parts["bodies"]:
        body.set_facecolor(color)
        body.set_edgecolor("#4F555B")
        body.set_linewidth(0.45)
        body.set_alpha(0.86)
    q25, median, q75 = np.quantile(array, [0.25, 0.5, 0.75])
    ax.vlines(position, q25, q75, color="#2E3338", lw=2.1, zorder=4)
    ax.scatter(
        [position], [median], s=12, color="white", edgecolor="#2E3338", linewidth=0.5, zorder=5
    )


def style_violin_axis(ax: plt.Axes, genes: list[str], title: str, subtitle: str) -> None:
    ax.set_xlim(-0.65, len(genes) - 0.35)
    ax.set_ylim(-0.08, 4.12)
    ax.set_xticks(range(len(genes)), genes, fontsize=8.2)
    ax.set_ylabel("Expression, log1p(CP10K)", fontsize=8)
    ax.set_yticks([0, 1, 2, 3, 4])
    ax.grid(axis="y", color="#E4E7EA", lw=0.55, zorder=0)
    ax.tick_params(axis="x", length=0)
    ax.tick_params(axis="y", labelsize=7, length=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_title(title, fontsize=10.2, fontweight="bold", pad=20)
    ax.text(
        0.5,
        1.015,
        subtitle,
        transform=ax.transAxes,
        fontsize=7.1,
        color="#43484E",
        ha="center",
        va="bottom",
    )


def plot_cortical_violins(
    ax: plt.Axes, parent: pd.DataFrame, selected: pd.DataFrame
) -> None:
    genes = ["LHX6", "ERBB4", "PLXNA2", "NRP2"]
    offsets = [-0.17, 0.17]
    for x, gene in enumerate(genes):
        draw_violin(ax, parent[gene], x + offsets[0], 0.29, BLUE_LIGHT)
        draw_violin(ax, selected[gene], x + offsets[1], 0.29, BLUE)
    style_violin_axis(
        ax,
        genes,
        "Cortical Cluster 3 expression distributions",
        "white points = medians; thick lines = interquartile ranges",
    )
    ax.legend(
        handles=[
            Patch(facecolor=BLUE_LIGHT, edgecolor="#4F555B", label=f"Cluster 3, all cells (n={len(parent):,})"),
            Patch(facecolor=BLUE, edgecolor="#4F555B", label=f"four-gene+ cells (n={len(selected):,})"),
        ],
        frameon=False,
        fontsize=7,
        loc="upper left",
        ncol=2,
        handlelength=1.2,
        columnspacing=1.0,
    )


def plot_subpallial_violins(
    ax: plt.Axes,
    reference: pd.DataFrame,
    aggregate_lhx8: pd.DataFrame,
    low_subset: pd.DataFrame,
    plxna2_threshold: float,
    nrp2_threshold: float,
) -> None:
    genes = ["LHX8", "PLXNA2", "NRP2"]
    offsets = [-0.25, 0.0, 0.25]
    colors = [GRAY_DARK, MAGENTA_LIGHT, MAGENTA]
    frames = [reference, aggregate_lhx8, low_subset]
    for x, gene in enumerate(genes):
        for offset, color, frame in zip(offsets, colors, frames, strict=True):
            draw_violin(ax, frame[gene], x + offset, 0.21, color)
    style_violin_axis(
        ax,
        genes,
        "LHX8+ subpallial expression distributions",
        "the selected subset allows low positive PLXNA2/NRP2 expression",
    )
    for x, threshold, gene in [
        (1, plxna2_threshold, "PLXNA2"),
        (2, nrp2_threshold, "NRP2"),
    ]:
        ax.hlines(threshold, x - 0.39, x + 0.39, color="#7F0000", lw=0.9, ls=(0, (3, 2)))
        ax.text(
            x + 0.40,
            threshold,
            f"low cutoff {threshold:.2f}",
            fontsize=6.2,
            color="#7F0000",
            va="center",
            ha="left",
        )
    ax.legend(
        handles=[
            Patch(facecolor=GRAY_DARK, edgecolor="#4F555B", label=f"GP/striatal LHX8+ reference (n={len(reference):,})"),
            Patch(facecolor=MAGENTA_LIGHT, edgecolor="#4F555B", label=f"four-cluster LHX8+ aggregate (n={len(aggregate_lhx8):,})"),
            Patch(facecolor=MAGENTA, edgecolor="#4F555B", label=f"lower PLXNA2/NRP2 subset (n={len(low_subset):,})"),
        ],
        frameon=False,
        fontsize=6.8,
        loc="upper left",
        ncol=1,
        handlelength=1.2,
    )


def main() -> None:
    args = parse_args()
    project_root = args.project_root.resolve()
    input_path = (args.input or project_root / TABLE_RELATIVE).resolve()
    output_path = (
        args.output
        or project_root / "final_figures" / FIGURE_ID / "figures" / "png" / OUTPUT_NAME
    ).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    data = load_data(input_path)
    summary = cluster_summary(data)
    cortical_label = select_cortical_cluster(summary)
    subpallial_labels = select_subpallial_clusters(data, summary)
    if not subpallial_labels:
        raise ValueError("No subpallial clusters met the prespecified selection criteria.")

    cortical = data.loc[data["set_id"] == "cortical_only"].copy()
    cortical_parent_mask = cortical["loupe_label"].eq(cortical_label)
    cortical_selected_mask = cortical_parent_mask & cortical[CORTICAL_GENES].gt(0).all(axis=1)
    cortical_parent = cortical.loc[cortical_parent_mask]
    cortical_selected = cortical.loc[cortical_selected_mask]

    subpallial = data.loc[data["set_id"] == "subpallial_only"].copy()
    subpallial_parent_mask = subpallial["loupe_label"].isin(subpallial_labels)
    subpallial_reference_mask = (
        subpallial["loupe_label"].isin(SUBPALLIAL_REFERENCE_LABELS)
        & subpallial["LHX8"].gt(0)
    )
    subpallial_reference = subpallial.loc[subpallial_reference_mask]
    plxna2_low_threshold = float(subpallial_reference["PLXNA2"].median())
    nrp2_low_threshold = float(subpallial_reference["NRP2"].median())
    subpallial_parent_lhx8_mask = subpallial_parent_mask & subpallial["LHX8"].gt(0)
    subpallial_selected_mask = (
        subpallial_parent_lhx8_mask
        & subpallial["PLXNA2"].le(plxna2_low_threshold)
        & subpallial["NRP2"].le(nrp2_low_threshold)
    )
    subpallial_parent = subpallial.loc[subpallial_parent_mask]
    subpallial_parent_lhx8 = subpallial.loc[subpallial_parent_lhx8_mask]
    subpallial_selected = subpallial.loc[subpallial_selected_mask]

    if len(cortical_selected) == 0 or len(subpallial_selected) == 0:
        raise ValueError("A selected comparison population is empty.")

    fig = plt.figure(figsize=(13.2, 8.7), facecolor="white")
    grid = fig.add_gridspec(
        2,
        2,
        height_ratios=[1.0, 0.62],
        left=0.055,
        right=0.975,
        bottom=0.105,
        top=0.855,
        wspace=0.23,
        hspace=0.34,
    )
    ax_cortical = fig.add_subplot(grid[0, 0])
    ax_subpallial = fig.add_subplot(grid[0, 1])
    ax_cortical_violin = fig.add_subplot(grid[1, 0])
    ax_subpallial_violin = fig.add_subplot(grid[1, 1])

    plot_selection_umap(
        ax_cortical,
        cortical,
        cortical_parent_mask,
        cortical_selected_mask,
        BLUE_LIGHT,
        BLUE,
        "Cortical-only Loupe recluster",
        f"{cortical_label} in light blue; {len(cortical_selected):,} four-gene+ cells in dark blue ({100 * len(cortical_selected) / len(cortical_parent):.1f}%)",
    )
    plot_selection_umap(
        ax_subpallial,
        subpallial,
        subpallial_parent_mask,
        subpallial_selected_mask,
        MAGENTA_LIGHT,
        MAGENTA,
        "Subpallial-only Loupe recluster",
        f"four parent subclusters in pink; {len(subpallial_selected):,} LHX8+ lower-expression cells in dark magenta",
    )
    plot_cortical_violins(ax_cortical_violin, cortical_parent, cortical_selected)
    plot_subpallial_violins(
        ax_subpallial_violin,
        subpallial_reference,
        subpallial_parent_lhx8,
        subpallial_selected,
        plxna2_low_threshold,
        nrp2_low_threshold,
    )
    panel_label(ax_cortical, "a", y=1.06)
    panel_label(ax_subpallial, "b", y=1.06)
    panel_label(ax_cortical_violin, "c", y=1.10)
    panel_label(ax_subpallial_violin, "d", y=1.10)

    fig.suptitle(
        "DIV90 cortical LHX6+ guidance state versus LHX8+ subpallial cells with lower PLXNA2/NRP2",
        fontsize=14.0,
        fontweight="bold",
        y=0.955,
    )
    fig.text(
        0.5,
        0.895,
        "UMAP locations and violin distributions from the existing pooled cortical-only and subpallial-only Loupe reclusters",
        fontsize=9.0,
        color="#43484E",
        ha="center",
    )
    fig.text(
        0.5,
        0.025,
        f"Subpallial 'lower expression' requires LHX8 detection and PLXNA2 ≤ {plxna2_low_threshold:.2f} plus NRP2 ≤ {nrp2_low_threshold:.2f} log1p(CP10K): "
        "the medians among LHX8+ GP/striatal reference cells. Zero expression is not required.",
        fontsize=7.1,
        color="#555B61",
        ha="center",
    )

    fig.savefig(output_path, dpi=args.dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    print(f"input={input_path}")
    print(f"output={output_path}")
    print(f"cortical_cluster={cortical_label}")
    print(f"cortical_parent_n={len(cortical_parent)}")
    print(f"cortical_four_gene_positive_n={len(cortical_selected)}")
    print(f"cortical_four_gene_positive_pct={100 * len(cortical_selected) / len(cortical_parent):.3f}")
    print(f"subpallial_parent_labels={' | '.join(subpallial_labels)}")
    print(f"subpallial_parent_n={len(subpallial_parent)}")
    print(f"subpallial_parent_lhx8_positive_n={len(subpallial_parent_lhx8)}")
    print(f"subpallial_reference_lhx8_positive_n={len(subpallial_reference)}")
    print(f"plxna2_low_threshold={plxna2_low_threshold:.6f}")
    print(f"nrp2_low_threshold={nrp2_low_threshold:.6f}")
    print(f"subpallial_selected_n={len(subpallial_selected)}")
    print(f"subpallial_selected_pct_parent={100 * len(subpallial_selected) / len(subpallial_parent):.3f}")
    print(
        "subpallial_selected_pct_lhx8_positive="
        f"{100 * len(subpallial_selected) / len(subpallial_parent_lhx8):.3f}"
    )


if __name__ == "__main__":
    main()
