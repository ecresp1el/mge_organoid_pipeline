#!/usr/bin/env python3
"""Render final DIV30/DIV90 maturation-score UMAP overlays from Seurat scores."""

from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import subprocess
import textwrap
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import numpy as np
import pandas as pd

from mge_organoid_python.paths import resolve_project_root


BACKGROUND_GREY = "#d0d0d0"
SCORE_BLUE = "#0000ff"
GREY_BLUE_CMAP = mpl.colors.LinearSegmentedColormap.from_list("greyBlue", [BACKGROUND_GREY, SCORE_BLUE])
JIA_RGC_MEAN_COL = "jia_score_RGC1_RGC2_mean"

JIA_SCORE_COLUMNS = [
    (JIA_RGC_MEAN_COL, "RG scores"),
    ("jia_score_IPC", "IP scores"),
]
PREDEFINED_SCORE_COLUMNS = [
    (JIA_RGC_MEAN_COL, "RG scores"),
    ("jia_score_IPC", "IP scores"),
    ("immature_module_score", "Immature neuron score"),
    ("mature_module_score", "Mature neuron score"),
]

DIV30_PAPER_CLUSTER_MAP = {
    "0": ("1", "Radial glia"),
    "3": ("1", "Radial glia"),
    "7": ("1", "Radial glia"),
    "6": ("2", "Inhibitory progenitors"),
    "1": ("3", "SST+ cIN"),
    "4": ("4", "PV neuron precursor"),
    "2": ("5", "MGE subpallial neurons"),
}

IMMATURE_GENES = [
    "DCX",
    "STMN2",
    "STMN4",
    "SOX11",
    "TUBB3",
    "TUBB2B",
    "ELAVL4",
    "GAP43",
    "CXCR4",
    "ACKR3",
]
MATURE_GENES = [
    "RBFOX3",
    "SNAP25",
    "SYT1",
    "SYN1",
    "SYN2",
    "DLG4",
    "VAMP2",
    "SLC12A5",
    "GAD1",
    "GAD2",
    "SLC6A1",
    "ERBB4",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render final maturation-score UMAP overlays from a Seurat AddModuleScore table."
    )
    parser.add_argument("--project-root", default=os.environ.get("PROJECT_ROOT"))
    parser.add_argument("--marker-csv", default=None)
    parser.add_argument("--score-table", required=True, help="TSV/TSV.GZ exported by scripts/38_export_maturation_module_scores_seurat.R.")
    parser.add_argument("--run-label", default="maturation_scores_v1")
    parser.add_argument("--final-folder", default="maturation_scores")
    parser.add_argument("--ctrl-size", type=int, default=50)
    parser.add_argument("--nbin", type=int, default=24)
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument("--point-size", type=float, default=None)
    return parser.parse_args()


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.5,
            "ytick.major.width": 0.5,
        }
    )


def git_status(repo_root: Path) -> dict[str, str]:
    def run_git(args: list[str]) -> str:
        try:
            return subprocess.check_output(
                ["git", *args],
                cwd=repo_root,
                text=True,
                stderr=subprocess.DEVNULL,
            ).strip()
        except Exception:
            return ""

    return {
        "commit": run_git(["rev-parse", "HEAD"]),
        "status_short": run_git(["status", "--short"]),
    }


def copy_if_different(src: Path, dst: Path) -> None:
    src = src.resolve()
    dst = dst.resolve()
    if src == dst:
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def clear_output_files(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)
    for child in path.iterdir():
        if child.is_file() or child.is_symlink():
            child.unlink()


def read_score_table(path: Path) -> pd.DataFrame:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        data = pd.read_csv(handle, sep="\t", low_memory=False)
    required = {"dataset", "cell_id", "umap_1", "umap_2", "plot_include"}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Score table is missing required column(s): {sorted(missing)}")
    for column, _ in PREDEFINED_SCORE_COLUMNS:
        if column not in data.columns:
            raise ValueError(f"Score table is missing score column: {column}")
    return data


def bool_series(values: pd.Series) -> pd.Series:
    if values.dtype == bool:
        return values.fillna(False)
    text = values.astype(str).str.lower().str.strip()
    return text.isin(["true", "t", "1", "yes", "y"])


def prepare_plot_data(data: pd.DataFrame) -> pd.DataFrame:
    out = data.copy()
    out["plot_include"] = bool_series(out["plot_include"])
    out["umap_1_plot"] = pd.to_numeric(out["umap_1"], errors="coerce")
    out["umap_2_plot"] = pd.to_numeric(out["umap_2"], errors="coerce")
    div90 = out["dataset"].astype(str).eq("DIV90")
    out.loc[div90, "umap_2_plot"] = -1.0 * out.loc[div90, "umap_2_plot"]
    return out


def all_score_columns() -> list[tuple[str, str]]:
    seen = set()
    out = []
    for column, display in [*JIA_SCORE_COLUMNS, *PREDEFINED_SCORE_COLUMNS]:
        if column in seen:
            continue
        seen.add(column)
        out.append((column, display))
    return out


def score_color_norms(
    data: pd.DataFrame,
    columns: list[tuple[str, str]],
    table_dir: Path,
) -> tuple[dict[str, mpl.colors.Normalize], dict[str, dict[str, float]]]:
    plot_data = data.loc[data["plot_include"]]
    rows = []
    limits: dict[str, dict[str, float]] = {}
    norms: dict[str, mpl.colors.Normalize] = {}
    for column, display in columns:
        values = pd.to_numeric(plot_data[column], errors="coerce").to_numpy(dtype=float)
        finite = np.isfinite(values)
        if finite.any():
            raw_min = float(np.nanpercentile(values[finite], 1))
            raw_max = float(np.nanpercentile(values[finite], 99))
        else:
            raw_min, raw_max = 0.0, 1.0
        if raw_max <= raw_min:
            raw_max = raw_min + 1e-6
        norms[column] = mpl.colors.Normalize(vmin=raw_min, vmax=raw_max)
        limits[column] = {"raw_min": raw_min, "raw_max": raw_max}
        rows.append(
            {
                "score_column": column,
                "display": display,
                "plot_transform": "raw Seurat AddModuleScore values, clipped to this score column's 1st-99th percentile range across plotted cells",
                "display_vmin_1pct_plotted_cells": raw_min,
                "display_vmax_99pct_plotted_cells": raw_max,
            }
        )
    pd.DataFrame(rows).to_csv(table_dir / "maturation_score_color_scaling.tsv", sep="\t", index=False)
    return norms, limits


def clean_axis(ax: plt.Axes) -> None:
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_aspect("equal", adjustable="box")
    for spine in ax.spines.values():
        spine.set_visible(False)


def plot_score_layer(
    fig: plt.Figure,
    ax: plt.Axes,
    frame: pd.DataFrame,
    score_column: str,
    title: str,
    norm: mpl.colors.Normalize,
    point_size: float,
) -> mpl.collections.PathCollection:
    coords = frame[["umap_1_plot", "umap_2_plot"]].to_numpy(dtype=float)
    values = pd.to_numeric(frame[score_column], errors="coerce").to_numpy()
    finite = np.isfinite(values) & np.isfinite(coords[:, 0]) & np.isfinite(coords[:, 1])

    ax.scatter(
        coords[:, 0],
        coords[:, 1],
        s=point_size,
        c=BACKGROUND_GREY,
        alpha=1.0,
        linewidths=0,
        rasterized=True,
    )
    scatter = ax.scatter(
        coords[finite, 0],
        coords[finite, 1],
        s=point_size,
        c=values[finite],
        cmap=GREY_BLUE_CMAP,
        norm=norm,
        alpha=0.9,
        linewidths=0,
        rasterized=True,
    )
    ax.set_title(title, fontsize=9, fontweight="normal", pad=4)
    clean_axis(ax)
    return scatter


def cluster_labels_for_frame(frame: pd.DataFrame, dataset: str) -> pd.DataFrame:
    out = frame.copy()
    if dataset == "DIV30":
        raw = out["seurat_clusters"].astype(str)
        mapped = raw.map(DIV30_PAPER_CLUSTER_MAP)
        out["cluster_number"] = mapped.map(lambda value: value[0] if isinstance(value, tuple) else "")
        out["cluster_name"] = mapped.map(lambda value: value[1] if isinstance(value, tuple) else "")
        out["cluster_label"] = np.where(
            out["cluster_number"].astype(str).ne(""),
            out["cluster_number"].astype(str) + " - " + out["cluster_name"].astype(str),
            raw,
        )
    else:
        if "cluster_number_name" in out.columns:
            label = out["cluster_number_name"].astype(str).replace({"nan": ""})
        else:
            label = out["seurat_clusters"].astype(str)
        out["cluster_label"] = label
        out["cluster_number"] = out["cluster_label"].str.extract(r"^\s*([0-9]+)", expand=False).fillna("")
        out["cluster_name"] = out["cluster_label"].str.replace(r"^\s*[0-9]+\s*[-.]\s*", "", regex=True)
    out["cluster_sort"] = pd.to_numeric(out["cluster_number"], errors="coerce")
    return out


def draw_cluster_layer(ax: plt.Axes, frame: pd.DataFrame, dataset: str, point_size: float) -> pd.DataFrame:
    labeled = cluster_labels_for_frame(frame, dataset)
    labels = (
        labeled[["cluster_label", "cluster_sort"]]
        .drop_duplicates()
        .sort_values(["cluster_sort", "cluster_label"], na_position="last")
        ["cluster_label"]
        .tolist()
    )
    palette = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("tab20b").colors) + list(plt.get_cmap("tab20c").colors)
    color_map = {label: palette[idx % len(palette)] for idx, label in enumerate(labels)}
    ax.scatter(
        labeled["umap_1_plot"],
        labeled["umap_2_plot"],
        s=point_size,
        c=labeled["cluster_label"].map(color_map),
        alpha=0.86,
        linewidths=0,
        rasterized=True,
    )
    label_table = (
        labeled.groupby("cluster_label", sort=False)
        .agg(
            n_cells=("cell_id", "size"),
            cluster_number=("cluster_number", "first"),
            cluster_name=("cluster_name", "first"),
        )
        .reset_index()
    )
    label_table["cluster_sort"] = pd.to_numeric(label_table["cluster_number"], errors="coerce")
    label_table = label_table.sort_values(["cluster_sort", "cluster_label"], na_position="last")

    x = pd.to_numeric(labeled["umap_1_plot"], errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(labeled["umap_2_plot"], errors="coerce").to_numpy(dtype=float)
    finite = np.isfinite(x) & np.isfinite(y)
    x_min, x_max = float(np.nanmin(x[finite])), float(np.nanmax(x[finite]))
    y_min, y_max = float(np.nanmin(y[finite])), float(np.nanmax(y[finite]))
    x_range = max(x_max - x_min, 1e-6)
    y_range = max(y_max - y_min, 1e-6)
    legend_pad = 0.95 * x_range if dataset == "DIV90" else 0.70 * x_range
    ax.set_xlim(x_min - 0.04 * x_range, x_max + legend_pad)
    ax.set_ylim(y_min - 0.04 * y_range, y_max + 0.04 * y_range)

    y_positions = np.linspace(y_max - 0.04 * y_range, y_min + 0.04 * y_range, label_table.shape[0])
    x_marker = x_max + 0.08 * x_range
    x_text = x_max + 0.13 * x_range
    for (_, row), y_pos in zip(label_table.iterrows(), y_positions):
        label = textwrap.fill(str(row["cluster_label"]), width=28 if dataset == "DIV90" else 22)
        color = color_map.get(str(row["cluster_label"]), "#333333")
        ax.scatter([x_marker], [y_pos], s=16, c=[color], linewidths=0, clip_on=False, zorder=5)
        ax.text(
            x_text,
            y_pos,
            label,
            ha="left",
            va="center",
            fontsize=4.8 if dataset == "DIV90" else 5.6,
            color="black",
            zorder=5,
            linespacing=0.9,
            clip_on=False,
        )
    ax.set_title(f"{dataset} clusters", fontsize=9, fontweight="normal", pad=4)
    clean_axis(ax)
    label_table.insert(0, "dataset", dataset)
    return label_table


def render_grid(
    data: pd.DataFrame,
    columns: list[tuple[str, str]],
    norms: dict[str, mpl.colors.Normalize],
    output_stem: Path,
    dpi: int,
    point_size: float | None,
    table_dir: Path | None = None,
    include_cluster_column: bool = False,
) -> dict[str, dict[str, float]]:
    datasets = ["DIV30", "DIV90"]
    n_cols = len(columns) + (1 if include_cluster_column else 0)
    width_ratios = [1.55, *([1.0] * len(columns))] if include_cluster_column else [1.0] * len(columns)
    fig, axes = plt.subplots(
        len(datasets),
        n_cols,
        figsize=(2.75 * n_cols, 2.75 * len(datasets)),
        constrained_layout=True,
        squeeze=False,
        gridspec_kw={"width_ratios": width_ratios},
    )
    cluster_rows = []
    score_scatters: dict[str, mpl.collections.PathCollection] = {}

    for row_idx, dataset in enumerate(datasets):
        frame = data.loc[data["dataset"].astype(str).eq(dataset) & data["plot_include"]].copy()
        size = point_size if point_size is not None else (0.25 if len(frame) > 80000 else 0.75)
        offset = 0
        if include_cluster_column:
            cluster_rows.append(draw_cluster_layer(axes[row_idx, 0], frame, dataset, size))
            offset = 1
        for col_idx, (score_column, display) in enumerate(columns):
            score_scatters[score_column] = plot_score_layer(
                fig,
                axes[row_idx, col_idx + offset],
                frame,
                score_column,
                f"{dataset} {display}",
                norms[score_column],
                size,
            )

    offset = 1 if include_cluster_column else 0
    for col_idx, (score_column, _) in enumerate(columns):
        scatter = score_scatters.get(score_column)
        if scatter is None:
            continue
        cbar = fig.colorbar(
            scatter,
            ax=axes[:, col_idx + offset].ravel().tolist(),
            fraction=0.030,
            pad=0.012,
            shrink=0.46,
            aspect=10,
        )
        cbar.ax.tick_params(labelsize=6, length=2, width=0.4, pad=1)
        cbar.outline.set_linewidth(0.4)

    for ext in ["png", "pdf", "svg"]:
        path = output_stem.with_suffix(f".{ext}")
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        print(f"[MaturationScores] wrote {path}", flush=True)
    plt.close(fig)
    if include_cluster_column and table_dir is not None and cluster_rows:
        pd.concat(cluster_rows, ignore_index=True).to_csv(
            table_dir / "predefined_maturation_scores_grid_cluster_label_positions.tsv",
            sep="\t",
            index=False,
        )
    return {column: {"plot_vmin": norms[column].vmin, "plot_vmax": norms[column].vmax} for column, _ in columns}


def render_div90_cluster_reference(
    data: pd.DataFrame,
    output_stem: Path,
    table_dir: Path,
    dpi: int,
    point_size: float | None,
) -> None:
    frame = data.loc[data["dataset"].astype(str).eq("DIV90") & data["plot_include"]].copy()
    if frame.empty:
        return
    if "cluster_number_name" in frame.columns:
        label_col = "cluster_number_name"
    elif "seurat_clusters" in frame.columns:
        label_col = "seurat_clusters"
    else:
        label_col = "div90_visualization_cluster"

    frame["cluster_label"] = frame[label_col].astype(str).replace({"nan": ""})
    frame.loc[frame["cluster_label"].eq(""), "cluster_label"] = frame.loc[frame["cluster_label"].eq(""), "div90_visualization_cluster"].astype(str)
    frame["cluster_number"] = frame["cluster_label"].str.extract(r"^\s*([0-9]+)", expand=False)
    frame["cluster_sort"] = pd.to_numeric(frame["cluster_number"], errors="coerce")
    labels = (
        frame[["cluster_label", "cluster_sort"]]
        .drop_duplicates()
        .sort_values(["cluster_sort", "cluster_label"], na_position="last")
        ["cluster_label"]
        .tolist()
    )
    palette = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("tab20b").colors) + list(plt.get_cmap("tab20c").colors)
    color_map = {label: palette[idx % len(palette)] for idx, label in enumerate(labels)}

    size = point_size if point_size is not None else 0.75
    fig, ax = plt.subplots(figsize=(5.2, 4.7), constrained_layout=True)
    ax.scatter(
        frame["umap_1_plot"],
        frame["umap_2_plot"],
        s=size,
        c=frame["cluster_label"].map(color_map),
        alpha=0.85,
        linewidths=0,
        rasterized=True,
    )
    positions = (
        frame.groupby("cluster_label", sort=False)
        .agg(
            umap_1_plot=("umap_1_plot", "median"),
            umap_2_plot=("umap_2_plot", "median"),
            n_cells=("cell_id", "size"),
            cluster_number=("cluster_number", "first"),
        )
        .reset_index()
    )
    for _, row in positions.iterrows():
        text = ax.text(
            row["umap_1_plot"],
            row["umap_2_plot"],
            str(row["cluster_label"]),
            ha="center",
            va="center",
            fontsize=6.5,
            color="black",
            zorder=5,
        )
        text.set_path_effects([path_effects.withStroke(linewidth=2.5, foreground="white")])

    ax.set_title(f"DIV90 clusters used for score overlays\nclusters 6/7 removed, n={len(frame):,}", fontsize=10)
    clean_axis(ax)
    for ext in ["png", "pdf", "svg"]:
        path = output_stem.with_suffix(f".{ext}")
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        print(f"[MaturationScores] wrote {path}", flush=True)
    plt.close(fig)

    counts = (
        frame.groupby(["cluster_number", "cluster_label"], dropna=False, sort=False)
        .size()
        .reset_index(name="n_plotted_cells")
    )
    counts["cluster_sort"] = pd.to_numeric(counts["cluster_number"], errors="coerce")
    counts = counts.sort_values(["cluster_sort", "cluster_label"], na_position="last").drop(columns=["cluster_sort"])
    counts.to_csv(table_dir / "div90_cluster_number_name_plot_counts.tsv", sep="\t", index=False)


def write_plot_filter_summary(data: pd.DataFrame, table_dir: Path) -> None:
    rows = []
    for dataset, frame in data.groupby("dataset", sort=False):
        include = frame["plot_include"].astype(bool)
        rows.append(
            {
                "dataset": dataset,
                "n_total_cells": int(len(frame)),
                "n_plotted_cells": int(include.sum()),
                "n_excluded_cells": int((~include).sum()),
                "div90_visualization_cluster_col": frame.get("div90_visualization_cluster_col", pd.Series([""])).astype(str).replace("nan", "").iloc[0],
                "rule": "DIV90 excludes current clusters 6 and 7; DIV90 UMAP2 is multiplied by -1 for plotting" if str(dataset) == "DIV90" else "all cells plotted",
            }
        )
    pd.DataFrame(rows).to_csv(table_dir / "maturation_scores_plot_filter_summary.tsv", sep="\t", index=False)


def write_readme(
    final_dir: Path,
    run_dir: Path,
    marker_csv: Path,
    score_table: Path,
    ctrl_size: int,
    nbin: int,
    random_state: int,
) -> None:
    readme = f"""# Maturation Scores

Unified DIV30/DIV90 UMAP overlays for Jia and predefined maturation modules.

## Scoring

- Jia marker source: `{marker_csv}`
- Scoring method for every displayed module: `Seurat::AddModuleScore`
- Parameters: `ctrl={ctrl_size}`, `nbin={nbin}`, `seed={random_state}`, assay `RNA`
- Jia displays: `mean(jia_score_RGC1, jia_score_RGC2)` and `jia_score_IPC`
- Predefined displays: the Jia RGC1/RGC2 mean, Jia IPC, immature module, and mature module
- No mature-minus-immature MGE maturation score is computed or plotted in this final package.

Immature module genes requested:
`{", ".join(IMMATURE_GENES)}`

Mature module genes requested:
`{", ".join(MATURE_GENES)}`

Gene-level coverage for the requested and resolved genes is in
`tables/maturation_score_gene_report.tsv`. The requested module gene list is in
`tables/maturation_score_module_gene_sets_requested.tsv`.

## DIV90 Visualization

DIV90 plotting follows the final-figure handoff rule used elsewhere: current
clusters 6 and 7 are excluded from the UMAP view as stressed cells, and the
plotted DIV90 coordinates use `UMAP2_plot = -1 * UMAP2`. The per-cell score
table still retains all cells and records the plot-inclusion flag.

## Color Scales

Score overlays are rendered with the Shi-style grey-to-blue colormap, but each
score column autoscales to its raw Seurat AddModuleScore 1st-99th percentile
range across plotted DIV30/DIV90 cells. The display ranges are recorded in
`tables/maturation_score_color_scaling.tsv`.

## Inputs

- Seurat score table: `{score_table}`
- Reproducible run directory: `{run_dir}`

## Outputs

- Jia module overlays: `figures/png/maturation_scores_umap_grid.png`,
  `figures/pdf/maturation_scores_umap_grid.pdf`, and
  `figures/svg/maturation_scores_umap_grid.svg`
- Predefined module overlays:
  `figures/png/predefined_maturation_scores_umap_grid.png`,
  `figures/pdf/predefined_maturation_scores_umap_grid.pdf`, and
  `figures/svg/predefined_maturation_scores_umap_grid.svg`. The first column
  is a DIV30/DIV90 cluster-number/name UMAP with the labels listed beside the
  UMAP rather than over the point cloud.
- DIV90 cluster reference UMAP:
  `figures/png/div90_cluster_number_name_reference_umap.png`,
  `figures/pdf/div90_cluster_number_name_reference_umap.pdf`, and
  `figures/svg/div90_cluster_number_name_reference_umap.svg`
- Tables and gene audits are in `tables/`
"""
    (final_dir / "README.md").write_text(readme, encoding="utf-8")


def main() -> int:
    args = parse_args()
    configure_matplotlib()
    repo_root = Path(__file__).resolve().parents[2]
    project_root = resolve_project_root(args.project_root)
    marker_csv = Path(args.marker_csv or project_root / "reference" / "Jia_et_al_2026_Science_3_progs.csv").expanduser().resolve()
    score_table = Path(args.score_table).expanduser().resolve()
    if not score_table.exists():
        raise FileNotFoundError(score_table)

    run_dir = project_root / "results" / "maturation_scores" / args.run_label
    final_dir = project_root / "final_figures" / args.final_folder
    result_plot_dir = run_dir / "figures"
    result_table_dir = run_dir / "tables"
    final_plot_dirs = {
        "png": final_dir / "figures" / "png",
        "pdf": final_dir / "figures" / "pdf",
        "svg": final_dir / "figures" / "svg",
    }
    for path in [
        result_plot_dir,
        result_table_dir,
        final_dir / "code",
        final_dir / "tables",
        final_dir / "logs",
        final_dir / "provenance",
        *final_plot_dirs.values(),
    ]:
        path.mkdir(parents=True, exist_ok=True)

    for path in [final_dir / "tables", final_dir / "provenance", *final_plot_dirs.values()]:
        clear_output_files(path)

    data = prepare_plot_data(read_score_table(score_table))
    copy_if_different(score_table, result_table_dir / score_table.name)
    write_plot_filter_summary(data, result_table_dir)
    score_norms, score_limits = score_color_norms(data, all_score_columns(), result_table_dir)

    jia_limits = render_grid(
        data,
        JIA_SCORE_COLUMNS,
        score_norms,
        result_plot_dir / "maturation_scores_umap_grid",
        dpi=args.dpi,
        point_size=args.point_size,
    )
    predefined_limits = render_grid(
        data,
        PREDEFINED_SCORE_COLUMNS,
        score_norms,
        result_plot_dir / "predefined_maturation_scores_umap_grid",
        dpi=args.dpi,
        point_size=args.point_size,
        table_dir=result_table_dir,
        include_cluster_column=True,
    )
    render_div90_cluster_reference(
        data,
        result_plot_dir / "div90_cluster_number_name_reference_umap",
        result_table_dir,
        dpi=args.dpi,
        point_size=args.point_size,
    )

    for ext, dest_dir in final_plot_dirs.items():
        for stem in [
            "maturation_scores_umap_grid",
            "predefined_maturation_scores_umap_grid",
            "div90_cluster_number_name_reference_umap",
        ]:
            path = result_plot_dir / f"{stem}.{ext}"
            if path.exists():
                shutil.copy2(path, dest_dir / f"{stem}.{ext}")
    for table_path in result_table_dir.iterdir():
        if table_path.is_file():
            shutil.copy2(table_path, final_dir / "tables" / table_path.name)
    shutil.copy2(Path(__file__).resolve(), final_dir / "code" / Path(__file__).name)

    provenance = {
        "run_dir": str(run_dir),
        "final_dir": str(final_dir),
        "marker_csv": str(marker_csv),
        "score_table": str(score_table),
        "scoring_method": "Seurat::AddModuleScore",
        "ctrl_size": args.ctrl_size,
        "nbin": args.nbin,
        "random_state": args.random_state,
        "dpi": args.dpi,
        "div90_visualization_filter": "exclude current clusters 6 and 7 as stressed cells",
        "div90_plot_coordinate_transform": "UMAP1_plot = UMAP1; UMAP2_plot = -1 * UMAP2",
        "score_display_transform": "raw Seurat AddModuleScore values shown with per-score grey-to-blue 1st-99th percentile scaling across plotted cells",
        "score_display_raw_limits": score_limits,
        "color_limits": {"jia_grid": jia_limits, "predefined_grid": predefined_limits},
        "git": git_status(repo_root),
    }
    (run_dir / "maturation_scores_provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    shutil.copy2(run_dir / "maturation_scores_provenance.json", final_dir / "provenance" / "maturation_scores_provenance.json")
    write_readme(final_dir, run_dir, marker_csv, score_table, args.ctrl_size, args.nbin, args.random_state)
    print(f"[MaturationScores] final package: {final_dir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
