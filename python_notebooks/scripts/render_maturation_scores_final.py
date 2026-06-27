#!/usr/bin/env python3
"""Render unified Jia maturation-score UMAP overlays for DIV30 and DIV90."""

from __future__ import annotations

import argparse
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
import scanpy as sc

from mge_organoid_python.gene_program_scoring import (
    choose_first_existing,
    choose_umap_key,
    match_program_genes,
    programs_from_marker_table,
    read_marker_program_csv,
    safe_token,
    scanpy_score_genes_control_audit,
    score_programs_scanpy,
    select_program_markers,
)
from mge_organoid_python.loader import cached_h5ad_path
from mge_organoid_python.paths import resolve_project_root
from mge_organoid_python.studies import default_studies


PROGRAM_ORDER = ["RGC1", "RGC2", "IPC"]
PROGRAM_DISPLAY = {
    "RGC1": "Jia RGC1",
    "RGC2": "Jia RGC2",
    "IPC": "Jia IPC",
}
MATURATION_COL = "jia_maturation_index_IPC_minus_mean_RGC1_RGC2"
MATURATION_DISPLAY = "IPC - mean(RGC1, RGC2)"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Score Jia RGC1/RGC2/IPC programs in DIV30 and DIV90 and render final UMAP overlays."
    )
    parser.add_argument(
        "--project-root",
        default=os.environ.get("PROJECT_ROOT"),
        help="Runtime project root. Defaults to PROJECT_ROOT or the documented Great Lakes root.",
    )
    parser.add_argument(
        "--marker-csv",
        default=None,
        help="Jia marker CSV. Defaults to PROJECT_ROOT/reference/Jia_et_al_2026_Science_3_progs.csv.",
    )
    parser.add_argument(
        "--run-label",
        default="maturation_scores_v1",
        help="Results run label under PROJECT_ROOT/results/maturation_scores/.",
    )
    parser.add_argument(
        "--final-folder",
        default="maturation_scores",
        help="Final figure folder under PROJECT_ROOT/final_figures/.",
    )
    parser.add_argument("--ctrl-size", type=int, default=50)
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument("--top-n", type=int, default=None)
    parser.add_argument("--min-avg-log2fc", type=float, default=None)
    parser.add_argument("--max-p-val-adj", type=float, default=None)
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument("--point-size", type=float, default=None)
    parser.add_argument(
        "--skip-control-audit",
        action="store_true",
        help="Skip Scanpy control-gene audit tables if only the plots are needed.",
    )
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


def load_study(study_id: str, project_root: Path):
    studies = {study.study_id: study for study in default_studies()}
    if study_id not in studies:
        raise KeyError(f"Unknown study_id: {study_id}")
    h5ad_path = cached_h5ad_path(studies[study_id], project_root=project_root)
    if not h5ad_path.exists():
        raise FileNotFoundError(h5ad_path)
    print(f"[MaturationScores] loading {study_id}: {h5ad_path}", flush=True)
    adata = sc.read_h5ad(h5ad_path)
    print(
        f"[MaturationScores] loaded {study_id}: {adata.n_obs} cells x {adata.n_vars} genes",
        flush=True,
    )
    return adata, h5ad_path


def score_one(
    adata,
    dataset_label: str,
    selected_markers: pd.DataFrame,
    ctrl_size: int,
    random_state: int,
    table_dir: Path,
    skip_control_audit: bool,
) -> dict[str, str]:
    programs = programs_from_marker_table(
        selected_markers,
        gene_col="gene",
        program_col="cluster",
        program_order=PROGRAM_ORDER,
    )
    matched, overlap_summary, overlap_detail = match_program_genes(programs, adata.var_names)
    overlap_summary.insert(0, "dataset", dataset_label)
    overlap_detail.insert(0, "dataset", dataset_label)
    overlap_summary.to_csv(
        table_dir / f"{safe_token(dataset_label)}_jia_program_gene_overlap_summary.tsv",
        sep="\t",
        index=False,
    )
    overlap_detail.to_csv(
        table_dir / f"{safe_token(dataset_label)}_jia_program_gene_overlap_detail.tsv",
        sep="\t",
        index=False,
    )
    print(
        f"[MaturationScores] {dataset_label} gene overlap\n"
        + overlap_summary.to_string(index=False),
        flush=True,
    )

    score_columns = score_programs_scanpy(
        adata,
        matched,
        score_prefix="jia_score_",
        ctrl_size=ctrl_size,
        random_state=random_state,
    )
    adata.obs[MATURATION_COL] = (
        pd.to_numeric(adata.obs[score_columns["IPC"]], errors="coerce")
        - (
            pd.to_numeric(adata.obs[score_columns["RGC1"]], errors="coerce")
            + pd.to_numeric(adata.obs[score_columns["RGC2"]], errors="coerce")
        )
        / 2.0
    )

    if not skip_control_audit:
        control_summary, control_detail, program_bins = scanpy_score_genes_control_audit(
            adata,
            matched,
            score_columns=score_columns,
            ctrl_size=ctrl_size,
            random_state=random_state,
        )
        for frame in (control_summary, control_detail, program_bins):
            frame.insert(0, "dataset", dataset_label)
        control_summary.to_csv(
            table_dir / f"{safe_token(dataset_label)}_scanpy_control_gene_summary.tsv",
            sep="\t",
            index=False,
        )
        control_detail.to_csv(
            table_dir / f"{safe_token(dataset_label)}_scanpy_control_gene_detail.tsv",
            sep="\t",
            index=False,
        )
        program_bins.to_csv(
            table_dir / f"{safe_token(dataset_label)}_scanpy_program_gene_bins.tsv",
            sep="\t",
            index=False,
        )

    return score_columns


def obs_score_table(adata, dataset_label: str, score_columns: dict[str, str], cluster_col: str | None) -> pd.DataFrame:
    cols = []
    for candidate in ["cell_id", "orig.ident", "sample", "cluster_number_name", "seurat_clusters"]:
        if candidate in adata.obs.columns and candidate not in cols:
            cols.append(candidate)
    if cluster_col and cluster_col in adata.obs.columns and cluster_col not in cols:
        cols.append(cluster_col)
    cols.extend([score_columns[p] for p in PROGRAM_ORDER if p in score_columns])
    cols.append(MATURATION_COL)
    out = adata.obs[cols].copy()
    out.insert(0, "dataset", dataset_label)
    out.insert(1, "obs_name", adata.obs_names.astype(str))
    return out


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
    coords: np.ndarray,
    values: np.ndarray,
    title: str,
    cmap: str,
    point_size: float,
    center_zero: bool = False,
) -> None:
    finite = np.isfinite(values)
    ax.scatter(
        coords[:, 0],
        coords[:, 1],
        s=point_size,
        c="#d9d9d9",
        alpha=0.35,
        linewidths=0,
        rasterized=True,
    )
    if center_zero:
        lim = float(np.nanpercentile(np.abs(values[finite]), 99)) if finite.any() else 1.0
        lim = max(lim, 1e-6)
        norm = mpl.colors.TwoSlopeNorm(vmin=-lim, vcenter=0.0, vmax=lim)
    else:
        vmin = float(np.nanpercentile(values[finite], 1)) if finite.any() else 0.0
        vmax = float(np.nanpercentile(values[finite], 99)) if finite.any() else 1.0
        if vmax <= vmin:
            vmax = vmin + 1e-6
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
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


def render_umap_grid(
    scored: list[dict[str, object]],
    output_stem: Path,
    dpi: int,
    point_size: float | None,
) -> None:
    n_rows = len(scored)
    columns = [
        ("RGC1", "RGC1", "viridis", False),
        ("RGC2", "RGC2", "viridis", False),
        ("IPC", "IPC", "viridis", False),
        ("maturation", MATURATION_DISPLAY, "coolwarm", True),
    ]
    fig, axes = plt.subplots(
        n_rows,
        len(columns),
        figsize=(2.85 * len(columns), 2.75 * n_rows),
        constrained_layout=True,
        squeeze=False,
    )
    for row_idx, item in enumerate(scored):
        adata = item["adata"]
        dataset_label = str(item["label"])
        umap_key = str(item["umap_key"])
        score_columns = item["score_columns"]
        coords = np.asarray(adata.obsm[umap_key])
        size = point_size if point_size is not None else (0.25 if adata.n_obs > 80000 else 0.75)
        for col_idx, (program, display, cmap, center_zero) in enumerate(columns):
            if program == "maturation":
                values = pd.to_numeric(adata.obs[MATURATION_COL], errors="coerce").to_numpy()
            else:
                values = pd.to_numeric(adata.obs[score_columns[program]], errors="coerce").to_numpy()
            title = f"{dataset_label} {display}"
            plot_score_layer(fig, axes[row_idx, col_idx], coords, values, title, cmap, size, center_zero)

    for ext in ["png", "pdf", "svg"]:
        path = output_stem.with_suffix(f".{ext}")
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        print(f"[MaturationScores] wrote {path}", flush=True)
    plt.close(fig)


def write_readme(
    final_dir: Path,
    run_dir: Path,
    marker_csv: Path,
    h5ad_paths: dict[str, Path],
    ctrl_size: int,
    random_state: int,
) -> None:
    readme = f"""# Maturation Scores

Unified DIV30/DIV90 Jia program score UMAP overlays.

## Scoring

- Marker source: `{marker_csv}`
- Programs: `RGC1`, `RGC2`, `IPC`
- Method: `scanpy.tl.score_genes` on each AnnData `.X`, `use_raw=False`
- Control genes: expression-binned Scanpy controls, `ctrl_size={ctrl_size}`, `random_state={random_state}`
- Derived maturation index: `jia_score_IPC - mean(jia_score_RGC1, jia_score_RGC2)`

The derived maturation index is only a display summary. The three Jia program
scores are exported separately and should remain the primary score columns.

## Inputs

- DIV30 H5AD: `{h5ad_paths["DIV30"]}`
- DIV90 H5AD: `{h5ad_paths["DIV90"]}`

## Outputs

- Main overlays: `figures/png/maturation_scores_umap_grid.png`,
  `figures/pdf/maturation_scores_umap_grid.pdf`, and
  `figures/svg/maturation_scores_umap_grid.svg`
- Tables and audits are in `tables/`
- Reproducible run outputs are mirrored from `{run_dir}`
"""
    (final_dir / "README.md").write_text(readme, encoding="utf-8")


def main() -> int:
    args = parse_args()
    configure_matplotlib()
    repo_root = Path(__file__).resolve().parents[2]
    project_root = resolve_project_root(args.project_root)
    marker_csv = Path(args.marker_csv or project_root / "reference" / "Jia_et_al_2026_Science_3_progs.csv")
    marker_csv = marker_csv.expanduser().resolve()
    if not marker_csv.exists():
        raise FileNotFoundError(marker_csv)

    run_dir = project_root / "results" / "maturation_scores" / args.run_label
    final_dir = project_root / "final_figures" / args.final_folder
    result_plot_dir = run_dir / "figures"
    result_table_dir = run_dir / "tables"
    final_plot_dirs = {
        "png": final_dir / "figures" / "png",
        "pdf": final_dir / "figures" / "pdf",
        "svg": final_dir / "figures" / "svg",
    }
    for path in [result_plot_dir, result_table_dir, final_dir / "code", final_dir / "tables", final_dir / "logs", final_dir / "provenance", *final_plot_dirs.values()]:
        path.mkdir(parents=True, exist_ok=True)

    markers = read_marker_program_csv(marker_csv, gene_col="gene", program_col="cluster")
    selected_markers = select_program_markers(
        markers,
        gene_col="gene",
        program_col="cluster",
        top_n=args.top_n,
        min_avg_log2fc=args.min_avg_log2fc,
        max_p_val_adj=args.max_p_val_adj,
    )
    selected_markers.to_csv(result_table_dir / "jia_program_markers_selected.tsv", sep="\t", index=False)
    markers.to_csv(result_table_dir / "jia_program_markers_full.tsv", sep="\t", index=False)

    h5ad_paths: dict[str, Path] = {}
    scored = []
    all_obs_tables = []
    for dataset_label, study_id in [("DIV30", "varela_div30"), ("DIV90", "varela_div90")]:
        adata, h5ad_path = load_study(study_id, project_root)
        h5ad_paths[dataset_label] = h5ad_path
        cluster_col = choose_first_existing(
            adata.obs.columns,
            ["cluster_number_name", "paper_cluster_annotation", "seurat_clusters", "RNA_snn_res.0.5", "RNA_snn_res.0.2"],
        )
        umap_key = choose_umap_key(adata)
        score_columns = score_one(
            adata,
            dataset_label,
            selected_markers,
            args.ctrl_size,
            args.random_state,
            result_table_dir,
            args.skip_control_audit,
        )
        all_obs_tables.append(obs_score_table(adata, dataset_label, score_columns, cluster_col))
        scored.append(
            {
                "label": dataset_label,
                "adata": adata,
                "umap_key": umap_key,
                "score_columns": score_columns,
                "cluster_col": cluster_col or "",
            }
        )

    obs_scores = pd.concat(all_obs_tables, ignore_index=True)
    obs_scores.to_csv(result_table_dir / "div30_div90_jia_maturation_scores_obs.tsv.gz", sep="\t", index=False)

    render_umap_grid(
        scored,
        result_plot_dir / "maturation_scores_umap_grid",
        dpi=args.dpi,
        point_size=args.point_size,
    )

    for ext, dest_dir in final_plot_dirs.items():
        shutil.copy2(result_plot_dir / f"maturation_scores_umap_grid.{ext}", dest_dir / f"maturation_scores_umap_grid.{ext}")
    for table_path in result_table_dir.iterdir():
        if table_path.is_file():
            shutil.copy2(table_path, final_dir / "tables" / table_path.name)
    shutil.copy2(Path(__file__).resolve(), final_dir / "code" / Path(__file__).name)

    provenance = {
        "run_dir": str(run_dir),
        "final_dir": str(final_dir),
        "marker_csv": str(marker_csv),
        "h5ad_paths": {key: str(value) for key, value in h5ad_paths.items()},
        "ctrl_size": args.ctrl_size,
        "random_state": args.random_state,
        "top_n": args.top_n,
        "min_avg_log2fc": args.min_avg_log2fc,
        "max_p_val_adj": args.max_p_val_adj,
        "dpi": args.dpi,
        "git": git_status(repo_root),
    }
    (run_dir / "maturation_scores_provenance.json").write_text(
        json.dumps(provenance, indent=2),
        encoding="utf-8",
    )
    shutil.copy2(run_dir / "maturation_scores_provenance.json", final_dir / "provenance" / "maturation_scores_provenance.json")
    write_readme(final_dir, run_dir, marker_csv, h5ad_paths, args.ctrl_size, args.random_state)
    print(f"[MaturationScores] final package: {final_dir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
