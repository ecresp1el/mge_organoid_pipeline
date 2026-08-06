#!/usr/bin/env python3
"""Render pooled cortical-versus-subpallial DIV90 subpopulation violins."""

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
        "axes.linewidth": 0.75,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "mathtext.fontset": "dejavusans",
    }
)

import matplotlib.pyplot as plt
import anndata as ad
import numpy as np
import pandas as pd


PROJECT_ROOT_DEFAULT = Path("/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder")
FIGURE_ID = "fig_div90_loupe_recluster_annotations_v1"
INPUT_RELATIVE = Path(
    f"final_figures/{FIGURE_ID}/tables/"
    "div90_loupe_subcluster_guidance_joined_per_cell.tsv.gz"
)
H5AD_RELATIVE = Path("results/python_anndata/varela_div90.h5ad")
MEMBERSHIP_RELATIVE = Path(
    f"final_figures/{FIGURE_ID}/tables/div90_loupe_recluster_membership.tsv.gz"
)
OUTPUT_NAME = "div90_guidance_receptor_subpopulations.png"
GENES = ["LHX6", "LHX8", "ERBB4", "PLXNA2", "NRP2", "NRP1"]
CORTICAL_SELECTION_GENES = ["LHX6", "ERBB4", "PLXNA2", "NRP2"]
SUBPALLIAL_REFERENCE_LABELS = ["1+3 GP projection", "7+9 striatal?"]
DETECTION_CUTOFF = 0.0
EXCLUDED_DIV90_CLUSTERS = {"6", "7"}

BLUE = "#397CA8"
MAGENTA = "#B5486B"
FAMILY_COLORS = {
    "cortical_only": BLUE,
    "progenitor_astro_glia": "#168575",
    "subpallial_only": MAGENTA,
}
FAMILY_TITLES = {
    "cortical_only": "Cortical",
    "progenitor_astro_glia": "Progenitors, astro/glia",
    "subpallial_only": "Subpallial",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT_DEFAULT)
    parser.add_argument("--input", type=Path, default=None)
    parser.add_argument("--h5ad", type=Path, default=None)
    parser.add_argument("--membership", type=Path, default=None)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--dpi", type=int, default=400)
    return parser.parse_args()


def draw_violin(
    ax: plt.Axes,
    values: pd.Series,
    position: float,
    color: str,
) -> None:
    array = values.to_numpy(dtype=float)
    parts = ax.violinplot(
        [array],
        positions=[position],
        widths=0.52,
        showmeans=False,
        showmedians=False,
        showextrema=False,
        points=180,
        bw_method=0.22,
    )
    for body in parts["bodies"]:
        body.set_facecolor(color)
        body.set_edgecolor("#3F454A")
        body.set_linewidth(0.55)
        body.set_alpha(0.88)
    q25, median, q75 = np.quantile(array, [0.25, 0.5, 0.75])
    ax.vlines(position, q25, q75, color="#252A2E", lw=2.4, zorder=4)
    ax.scatter(
        [position],
        [median],
        s=17,
        color="white",
        edgecolor="#252A2E",
        linewidth=0.55,
        zorder=5,
    )


def load_data(path: Path) -> pd.DataFrame:
    data = pd.read_csv(path, sep="\t", dtype={"cell_id": str})
    data = data.loc[data["expression_available"].astype(bool)].copy()
    for gene in GENES:
        data[gene] = pd.to_numeric(data[gene], errors="raise")
    return data


def load_original_umap(path: Path) -> tuple[pd.DataFrame, str]:
    adata = ad.read_h5ad(path, backed="r")
    reduction = "X_umap_seurat" if "X_umap_seurat" in adata.obsm else "X_umap"
    coordinates = np.asarray(adata.obsm[reduction])
    frame = pd.DataFrame(
        {
            "cell_id": adata.obs_names.astype(str),
            "cluster_id": adata.obs["cluster_id"].astype(str).to_numpy(),
            "original_umap_1": coordinates[:, 0],
            "original_umap_2": -coordinates[:, 1],
        }
    )
    adata.file.close()
    return frame, reduction


def select_populations(
    data: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, float, float]:
    cortical = data.loc[
        data["set_id"] == "cortical_only"
    ]
    cortical_selected = cortical.loc[
        cortical[CORTICAL_SELECTION_GENES].gt(DETECTION_CUTOFF).all(axis=1)
    ].copy()

    subpallial = data.loc[data["set_id"] == "subpallial_only"].copy()
    reference = subpallial.loc[
        subpallial["loupe_label"].isin(SUBPALLIAL_REFERENCE_LABELS)
        & subpallial["LHX8"].gt(DETECTION_CUTOFF)
    ]
    plxna2_threshold = float(reference["PLXNA2"].median())
    nrp2_threshold = float(reference["NRP2"].median())
    subpallial_selected = subpallial.loc[
        ~subpallial["loupe_label"].isin(SUBPALLIAL_REFERENCE_LABELS)
        & subpallial["LHX8"].gt(DETECTION_CUTOFF)
    ].copy()
    if cortical_selected.empty or subpallial_selected.empty:
        raise ValueError("A pooled comparison population is empty.")
    return cortical_selected, subpallial_selected, plxna2_threshold, nrp2_threshold


def plot_selection_umap(
    ax: plt.Axes,
    family: pd.DataFrame,
    selected: pd.DataFrame,
    color: str,
    title: str,
    subtitle: str,
    x_column: str = "loupe_x",
    y_column: str = "loupe_y",
    background_size: float = 1.0,
    selected_size: float = 4.1,
) -> None:
    ax.scatter(
        family[x_column],
        family[y_column],
        s=background_size,
        color="#D9DCDF",
        linewidths=0,
        alpha=0.55,
        rasterized=True,
    )
    ax.scatter(
        selected[x_column],
        selected[y_column],
        s=selected_size,
        color=color,
        edgecolors="white",
        linewidths=0.12,
        alpha=0.92,
        rasterized=True,
    )
    for coordinate, setter in [(x_column, ax.set_xlim), (y_column, ax.set_ylim)]:
        low, high = np.quantile(family[coordinate].to_numpy(dtype=float), [0.005, 0.995])
        pad = max(0.1, 0.04 * (high - low))
        setter(low - pad, high + pad)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    if title:
        ax.set_title(title, fontsize=12.2, fontweight="bold", color=color, pad=18)
    if subtitle:
        ax.text(
            0.5,
            1.01,
            subtitle,
            transform=ax.transAxes,
            fontsize=8.8,
            color="#555A5E",
            ha="center",
            va="bottom",
        )


def add_panel_label(
    ax: plt.Axes, label: str, y: float = 1.08, x: float = -0.10
) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        ha="left",
        va="top",
    )


def main() -> None:
    args = parse_args()
    project_root = args.project_root.resolve()
    input_path = (args.input or project_root / INPUT_RELATIVE).resolve()
    h5ad_path = (args.h5ad or project_root / H5AD_RELATIVE).resolve()
    membership_path = (args.membership or project_root / MEMBERSHIP_RELATIVE).resolve()
    output_path = (
        args.output
        or project_root / "final_figures" / FIGURE_ID / "figures" / "png" / OUTPUT_NAME
    ).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    data = load_data(input_path)
    original_umap, original_reduction = load_original_umap(h5ad_path)
    original_umap = original_umap.loc[
        ~original_umap["cluster_id"].isin(EXCLUDED_DIV90_CLUSTERS)
    ].copy()
    retained_cell_ids = set(original_umap["cell_id"])
    membership = pd.read_csv(
        membership_path,
        sep="\t",
        usecols=["cell_id", "set_id"],
        dtype={"cell_id": str},
    ).drop_duplicates(["cell_id", "set_id"])
    membership = membership.loc[membership["cell_id"].isin(retained_cell_ids)].copy()
    membership_original = membership.merge(
        original_umap,
        on="cell_id",
        how="inner",
        validate="many_to_one",
    )
    data = data.loc[data["cell_id"].isin(retained_cell_ids)].copy()
    cortical, subpallial, plxna2_threshold, nrp2_threshold = select_populations(data)
    cortical_family = data.loc[data["set_id"] == "cortical_only"]
    subpallial_family = data.loc[data["set_id"] == "subpallial_only"]
    cortical_original = cortical[["cell_id"]].merge(
        original_umap, on="cell_id", how="inner", validate="one_to_one"
    )
    subpallial_original = subpallial[["cell_id"]].merge(
        original_umap, on="cell_id", how="inner", validate="one_to_one"
    )
    if len(cortical_original) != len(cortical) or len(subpallial_original) != len(subpallial):
        raise ValueError("Not every selected cell mapped to the original DIV90 UMAP.")

    fig = plt.figure(figsize=(13.2, 6.8), facecolor="white")
    outer_grid = fig.add_gridspec(
        2,
        1,
        height_ratios=[0.65, 0.80],
        left=0.075,
        right=0.985,
        bottom=0.105,
        top=0.855,
        hspace=0.27,
    )
    family_grid = outer_grid[0, 0].subgridspec(1, 3, wspace=0.10)
    analysis_grid = outer_grid[1, 0].subgridspec(1, 2, wspace=0.14)
    family_axes = [fig.add_subplot(family_grid[0, index]) for index in range(3)]
    cortical_container = fig.add_subplot(analysis_grid[0, 0])
    subpallial_container = fig.add_subplot(analysis_grid[0, 1])
    for container in [cortical_container, subpallial_container]:
        container.set_axis_off()
    ax_cortical_violin = cortical_container.inset_axes([0.00, 0.00, 1.00, 0.57])
    ax_subpallial_violin = subpallial_container.inset_axes(
        [0.00, 0.00, 1.00, 0.57], sharey=ax_cortical_violin
    )
    ax_cortical_original = cortical_container.inset_axes([0.60, 0.63, 0.38, 0.34])
    ax_subpallial_original = subpallial_container.inset_axes([0.60, 0.63, 0.38, 0.34])

    family_order = ["cortical_only", "progenitor_astro_glia", "subpallial_only"]
    for ax, set_id in zip(family_axes, family_order, strict=True):
        family = membership_original.loc[membership_original["set_id"] == set_id]
        plot_selection_umap(
            ax,
            original_umap,
            family,
            FAMILY_COLORS[set_id],
            FAMILY_TITLES[set_id],
            f"major population membership (n = {family['cell_id'].nunique():,})",
            x_column="original_umap_1",
            y_column="original_umap_2",
        )

    plot_selection_umap(
        ax_cortical_original,
        original_umap,
        cortical_original,
        BLUE,
        "",
        "",
        x_column="original_umap_1",
        y_column="original_umap_2",
        background_size=0.45,
        selected_size=2.2,
    )
    plot_selection_umap(
        ax_subpallial_original,
        original_umap,
        subpallial_original,
        MAGENTA,
        "",
        "",
        x_column="original_umap_1",
        y_column="original_umap_2",
        background_size=0.45,
        selected_size=2.2,
    )

    violin_panels = [
        (ax_cortical_violin, cortical, BLUE),
        (ax_subpallial_violin, subpallial, MAGENTA),
    ]
    for ax, frame, color in violin_panels:
        for x, gene in enumerate(GENES):
            draw_violin(ax, frame[gene], x, color)
        ax.set_xlim(-0.58, len(GENES) - 0.42)
        ax.set_ylim(-0.08, 4.12)
        gene_labels = [rf"$\mathit{{{gene}}}$" for gene in GENES]
        ax.set_xticks(range(len(GENES)), gene_labels, fontsize=13.2)
        ax.set_yticks([0, 1, 2, 3, 4])
        ax.tick_params(axis="x", length=0)
        ax.tick_params(axis="y", labelsize=9.5, length=3.4, width=0.8, direction="out")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color("#24272A")
        ax.spines["bottom"].set_color("#24272A")
    ax_cortical_violin.set_ylabel(
        "Normalized expression\nlog1p(CP10K)", fontsize=10.8, labelpad=9
    )
    ax_subpallial_violin.tick_params(axis="y", labelleft=False)

    for container, title, frame, color in [
        (cortical_container, "Cortical subpopulation", cortical, BLUE),
        (subpallial_container, "Subpallial subpopulation", subpallial, MAGENTA),
    ]:
        container.text(
            0.02,
            0.965,
            title,
            transform=container.transAxes,
            fontsize=12.5,
            fontweight="bold",
            color=color,
            ha="left",
            va="top",
        )
        container.text(
            0.02,
            0.885,
            f"n = {len(frame):,} cells",
            transform=container.transAxes,
            fontsize=9.5,
            color="#555A5E",
            ha="left",
            va="top",
        )

    for ax, label in zip(family_axes, ["a", "b", "c"], strict=True):
        add_panel_label(ax, label, y=1.10, x=-0.17)
    add_panel_label(cortical_container, "d", y=1.00, x=-0.06)
    add_panel_label(subpallial_container, "e", y=1.00, x=-0.06)
    fig.text(
        0.065,
        0.975,
        "DIV90 guidance-receptor subpopulations",
        fontsize=12.5,
        fontweight="bold",
        color="#17191B",
        ha="left",
        va="top",
    )
    fig.savefig(output_path, dpi=args.dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    print(f"input={input_path}")
    print(f"h5ad={h5ad_path}")
    print(f"membership={membership_path}")
    print(f"original_reduction={original_reduction}")
    print(f"output={output_path}")
    print(f"cortical_subpopulation_n={len(cortical)}")
    print(f"subpallial_subpopulation_n={len(subpallial)}")
    print(f"plxna2_low_threshold={plxna2_threshold:.6f}")
    print(f"nrp2_low_threshold={nrp2_threshold:.6f}")


if __name__ == "__main__":
    main()
