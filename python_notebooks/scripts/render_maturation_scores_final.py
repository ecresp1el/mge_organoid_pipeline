#!/usr/bin/env python3
"""Render final DIV30/DIV90 maturation-score UMAP overlays from Seurat scores."""

from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import subprocess
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mge_organoid_python.paths import resolve_project_root


JIA_RGC_MEAN_COL = "jia_score_RGC1_RGC2_mean"

JIA_SCORE_COLUMNS = [
    (JIA_RGC_MEAN_COL, "Jia RGC1/RGC2 mean", "viridis"),
    ("jia_score_IPC", "Jia IPC", "viridis"),
]
PREDEFINED_SCORE_COLUMNS = [
    (JIA_RGC_MEAN_COL, "Jia RGC1/RGC2 mean", "viridis"),
    ("jia_score_IPC", "Jia IPC", "viridis"),
    ("immature_module_score", "Immature module", "magma"),
    ("mature_module_score", "Mature module", "magma"),
]

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
    for column, _, _ in PREDEFINED_SCORE_COLUMNS:
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


def clean_axis(ax: plt.Axes) -> None:
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_aspect("equal", adjustable="box")
    for spine in ax.spines.values():
        spine.set_visible(False)


def shared_norms(data: pd.DataFrame, columns: list[tuple[str, str, str]]) -> dict[str, mpl.colors.Normalize]:
    norms: dict[str, mpl.colors.Normalize] = {}
    plot_data = data.loc[data["plot_include"]]
    for column, _, _ in columns:
        values = pd.to_numeric(plot_data[column], errors="coerce").to_numpy()
        finite = np.isfinite(values)
        if finite.any():
            vmin = float(np.nanpercentile(values[finite], 1))
            vmax = float(np.nanpercentile(values[finite], 99))
        else:
            vmin, vmax = 0.0, 1.0
        if vmax <= vmin:
            vmax = vmin + 1e-6
        norms[column] = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    return norms


def plot_score_layer(
    fig: plt.Figure,
    ax: plt.Axes,
    frame: pd.DataFrame,
    score_column: str,
    title: str,
    cmap: str,
    norm: mpl.colors.Normalize,
    point_size: float,
) -> None:
    coords = frame[["umap_1_plot", "umap_2_plot"]].to_numpy(dtype=float)
    values = pd.to_numeric(frame[score_column], errors="coerce").to_numpy()
    finite = np.isfinite(values) & np.isfinite(coords[:, 0]) & np.isfinite(coords[:, 1])

    ax.scatter(
        coords[:, 0],
        coords[:, 1],
        s=point_size,
        c="#d0d0d0",
        alpha=0.35,
        linewidths=0,
        rasterized=True,
    )
    scatter = ax.scatter(
        coords[finite, 0],
        coords[finite, 1],
        s=point_size,
        c=values[finite],
        cmap=cmap,
        norm=norm,
        alpha=0.9,
        linewidths=0,
        rasterized=True,
    )
    ax.set_title(title, fontsize=9, fontweight="normal", pad=4)
    clean_axis(ax)
    cbar = fig.colorbar(scatter, ax=ax, fraction=0.042, pad=0.01)
    cbar.ax.tick_params(labelsize=6, length=2, width=0.4, pad=1)
    cbar.outline.set_linewidth(0.4)


def render_grid(
    data: pd.DataFrame,
    columns: list[tuple[str, str, str]],
    output_stem: Path,
    dpi: int,
    point_size: float | None,
) -> dict[str, dict[str, float]]:
    datasets = ["DIV30", "DIV90"]
    norms = shared_norms(data, columns)
    fig, axes = plt.subplots(
        len(datasets),
        len(columns),
        figsize=(2.75 * len(columns), 2.75 * len(datasets)),
        constrained_layout=True,
        squeeze=False,
    )

    for row_idx, dataset in enumerate(datasets):
        frame = data.loc[data["dataset"].astype(str).eq(dataset) & data["plot_include"]].copy()
        size = point_size if point_size is not None else (0.25 if len(frame) > 80000 else 0.75)
        for col_idx, (score_column, display, cmap) in enumerate(columns):
            plot_score_layer(
                fig,
                axes[row_idx, col_idx],
                frame,
                score_column,
                f"{dataset} {display}",
                cmap,
                norms[score_column],
                size,
            )

    for ext in ["png", "pdf", "svg"]:
        path = output_stem.with_suffix(f".{ext}")
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        print(f"[MaturationScores] wrote {path}", flush=True)
    plt.close(fig)
    return {column: {"vmin_1pct": norm.vmin, "vmax_99pct": norm.vmax} for column, _, _ in columns for norm in [norms[column]]}


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

Each score uses one shared color scale across DIV30 and DIV90, clipped to the
combined 1st and 99th percentiles of plotted cells for that score.

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
  `figures/svg/predefined_maturation_scores_umap_grid.svg`
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

    jia_limits = render_grid(
        data,
        JIA_SCORE_COLUMNS,
        result_plot_dir / "maturation_scores_umap_grid",
        dpi=args.dpi,
        point_size=args.point_size,
    )
    predefined_limits = render_grid(
        data,
        PREDEFINED_SCORE_COLUMNS,
        result_plot_dir / "predefined_maturation_scores_umap_grid",
        dpi=args.dpi,
        point_size=args.point_size,
    )

    for ext, dest_dir in final_plot_dirs.items():
        for stem in ["maturation_scores_umap_grid", "predefined_maturation_scores_umap_grid"]:
            shutil.copy2(result_plot_dir / f"{stem}.{ext}", dest_dir / f"{stem}.{ext}")
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
