#!/usr/bin/env python3
"""Render all-cell cluster UMAP inventories exported by Step 00."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
import re

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams.update(
    {
        "font.family": "Arial",
        "font.sans-serif": ["Arial", "Nimbus Sans", "Liberation Sans", "DejaVu Sans"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


STUDY_ORDER = [
    "varela_div30",
    "varela_div90",
    "walsh",
    "bershteyn_2025",
    "bershteyn_2023",
    "siebert_2026",
]
POINT_SIZE = 1.6


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--umap-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--manifest-outdir", type=Path, required=True)
    parser.add_argument("--png-dpi", type=int, default=300)
    parser.add_argument("--pdf-dpi", type=int, default=300)
    parser.add_argument("--svg-dpi", type=int, default=300)
    parser.add_argument("--make-svg", choices=["true", "false"], default="false")
    return parser.parse_args()


def natural_sort_key(value: object) -> list[object]:
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", str(value))]


def load_tables(umap_dir: Path) -> pd.DataFrame:
    frames = []
    for study_id in STUDY_ORDER:
        path = umap_dir / "per_study" / f"{study_id}_umap_cluster_inventory.tsv.gz"
        if not path.is_file():
            raise FileNotFoundError(f"Missing UMAP audit table: {path}")
        frame = pd.read_csv(path, sep="\t", low_memory=False)
        frame["study_id"] = frame["study_id"].astype(str)
        frames.append(frame)
    data = pd.concat(frames, ignore_index=True)
    if data["study_id"].nunique() != 6:
        raise ValueError("The UMAP inventory must contain exactly six studies.")
    data["UMAP1_plot"] = pd.to_numeric(data["UMAP1_original"], errors="coerce")
    data["UMAP2_plot"] = pd.to_numeric(data["UMAP2_original"], errors="coerce")
    div90 = data["study_id"].eq("varela_div90")
    data.loc[div90, "UMAP2_plot"] = -1.0 * data.loc[div90, "UMAP2_plot"]
    finite = np.isfinite(data["UMAP1_plot"]) & np.isfinite(data["UMAP2_plot"])
    if not finite.all():
        raise ValueError(f"Found {(~finite).sum()} cells without finite UMAP coordinates.")
    return data


def cluster_palette(labels: list[str]) -> dict[str, tuple[float, float, float, float]]:
    colors: list[tuple[float, float, float, float]] = []
    for name in ["tab20", "tab20b", "tab20c", "Set3", "Dark2", "Accent"]:
        cmap = plt.get_cmap(name)
        colors.extend(cmap(i) for i in range(getattr(cmap, "N", 8)))
    return {label: colors[index % len(colors)] for index, label in enumerate(labels)}


def plot_panel(ax: plt.Axes, subset: pd.DataFrame) -> dict[str, object]:
    study_id = str(subset["study_id"].iloc[0])
    display_name = str(subset["display_name"].iloc[0])
    labels = sorted(subset["display_cluster_label"].astype(str).unique(), key=natural_sort_key)
    colors = cluster_palette(labels)

    for label in labels:
        cluster_data = subset.loc[subset["display_cluster_label"].astype(str).eq(label)]
        ax.scatter(
            cluster_data["UMAP1_plot"],
            cluster_data["UMAP2_plot"],
            s=POINT_SIZE,
            c=[colors[label]],
            linewidths=0,
            rasterized=True,
        )

    positions = (
        subset.groupby(["display_cluster_id", "display_cluster_label"], sort=False, observed=True)
        .agg(x=("UMAP1_plot", "median"), y=("UMAP2_plot", "median"))
        .reset_index()
    )
    for row in positions.itertuples(index=False):
        ax.text(
            row.x,
            row.y,
            str(row.display_cluster_id),
            ha="center",
            va="center",
            fontsize=6.5,
            fontweight="bold",
            color="black",
            bbox={
                "boxstyle": "round,pad=0.12",
                "facecolor": "white",
                "edgecolor": "black",
                "linewidth": 0.25,
                "alpha": 0.80,
            },
        )

    n_raw = subset["raw_cluster_id"].astype(str).nunique()
    n_display = len(labels)
    count_text = f"{n_display} displayed clusters"
    if n_raw != n_display:
        count_text += f" from {n_raw} raw clusters"
    orientation = "DIV90 final-figure vertical orientation; all cells retained" if study_id == "varela_div90" else "original UMAP orientation; all cells retained"
    ax.set_title(f"{display_name}\n{count_text}, n={subset.shape[0]:,}\n{orientation}", fontsize=9)
    ax.set_aspect("equal", adjustable="box")
    ax.set_axis_off()

    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=colors[label], markersize=4, label=label)
        for label in labels
    ]
    ax.legend(
        handles=handles,
        title="Cluster number and name",
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        fontsize=5.2,
        title_fontsize=5.8,
        frameon=False,
        borderaxespad=0,
    )
    return {
        "study_id": study_id,
        "display_name": display_name,
        "n_cells": int(subset.shape[0]),
        "n_raw_clusters": int(n_raw),
        "n_display_clusters": int(n_display),
        "display_cluster_labels": " | ".join(labels),
        "plot_transform": str(subset["plot_transform"].iloc[0]),
        "all_cells_retained": True,
    }


def save_figure(
    fig: plt.Figure,
    prefix: Path,
    png_dpi: int,
    pdf_dpi: int,
    make_svg: bool,
    svg_dpi: int,
) -> list[str]:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(prefix.with_suffix(".png"), dpi=png_dpi, bbox_inches="tight")
    fig.savefig(prefix.with_suffix(".pdf"), dpi=pdf_dpi, bbox_inches="tight")
    extensions = ["png", "pdf"]
    if make_svg:
        fig.savefig(prefix.with_suffix(".svg"), dpi=svg_dpi, bbox_inches="tight")
        extensions.append("svg")
    plt.close(fig)
    return extensions


def render_figure(
    data: pd.DataFrame,
    study_ids: list[str],
    prefix: Path,
    title: str,
    png_dpi: int,
    pdf_dpi: int,
    make_svg: bool,
    svg_dpi: int,
    ncols: int,
) -> tuple[list[dict[str, object]], dict[str, object], list[str]]:
    subsets = [data.loc[data["study_id"].eq(study_id)].copy() for study_id in study_ids]
    if any(subset.empty for subset in subsets):
        raise ValueError(f"One or more requested studies are absent for {prefix.name}")
    nrows = int(np.ceil(len(subsets) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(7.0 * ncols, 5.2 * nrows), squeeze=False)
    summaries = []
    for ax, subset in zip(axes.ravel(), subsets):
        summaries.append(plot_panel(ax, subset))
    for ax in axes.ravel()[len(subsets) :]:
        ax.axis("off")
    fig.suptitle(title, fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.97), w_pad=7.0, h_pad=2.0)
    extensions = save_figure(fig, prefix, png_dpi, pdf_dpi, make_svg, svg_dpi)
    manifest = {
        "figure_id": prefix.name,
        "study_ids": ";".join(study_ids),
        "n_cells": int(sum(item["n_cells"] for item in summaries)),
        "png": str(prefix.with_suffix(".png")),
        "pdf": str(prefix.with_suffix(".pdf")),
        "svg": str(prefix.with_suffix(".svg")) if make_svg else "",
        "png_dpi": png_dpi,
        "pdf_dpi": pdf_dpi,
        "svg_enabled": make_svg,
        "svg_dpi": svg_dpi if make_svg else "",
        "all_cells_retained": True,
        "notes": "Cluster colors and centroid numbers; names in side legends; DIV90 UMAP2 plotting-only flip; no cells filtered.",
    }
    return summaries, manifest, extensions


def main() -> None:
    args = parse_args()
    make_svg = args.make_svg == "true"
    for format_name, dpi in {
        "PNG": args.png_dpi,
        "PDF": args.pdf_dpi,
        "SVG": args.svg_dpi,
    }.items():
        if dpi != 300:
            raise ValueError(f"{format_name} rasterized layers must use 300 dpi; got {dpi}")
    data = load_tables(args.umap_dir)
    args.manifest_outdir.mkdir(parents=True, exist_ok=True)
    output_subdirs = ["png", "pdf"] + (["svg"] if make_svg else [])
    for subdir in output_subdirs:
        (args.outdir / subdir).mkdir(parents=True, exist_ok=True)

    figure_specs = [
        (
            STUDY_ORDER,
            "six_study_input_umap_cluster_inventory",
            "Paper 2 registered input UMAP and cluster inventory — all cells",
            3,
        ),
        (
            ["varela_div30", "varela_div90"],
            "varela_div30_div90_input_umap_cluster_inventory",
            "Varela DIV30 and DIV90 input clusters — all cells",
            2,
        ),
    ]
    for study_id in STUDY_ORDER:
        figure_specs.append(
            (
                [study_id],
                f"{study_id}_input_umap_cluster_inventory",
                f"{study_id} input UMAP and cluster inventory — all cells",
                1,
            )
        )

    panel_rows: list[dict[str, object]] = []
    manifest_rows: list[dict[str, object]] = []
    for study_ids, basename, title, ncols in figure_specs:
        # Render once per format directory so paths follow the final-figure contract.
        staging_prefix = args.outdir / basename
        summaries, manifest, extensions = render_figure(
            data,
            study_ids,
            staging_prefix,
            title,
            args.png_dpi,
            args.pdf_dpi,
            make_svg,
            args.svg_dpi,
            ncols,
        )
        for extension in extensions:
            source = staging_prefix.with_suffix(f".{extension}")
            destination = args.outdir / extension / source.name
            source.replace(destination)
            manifest[extension] = str(destination)
        for summary in summaries:
            summary["figure_id"] = basename
            panel_rows.append(summary)
        manifest_rows.append(manifest)

    pd.DataFrame(panel_rows).to_csv(args.manifest_outdir / "umap_figure_panel_summary.tsv", sep="\t", index=False)
    pd.DataFrame(manifest_rows).to_csv(args.manifest_outdir / "umap_figure_manifest.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
