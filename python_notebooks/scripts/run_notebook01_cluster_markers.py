#!/usr/bin/env python
"""Run quick Notebook 01 Leiden cluster marker ranking."""

from __future__ import annotations

import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import pandas as pd
import scanpy as sc


def _add_repo_src_to_path() -> Path:
    start = Path(__file__).resolve()
    for candidate in [start.parent, *start.parents]:
        src_dir = candidate / "python_notebooks" / "src"
        if src_dir.exists():
            sys.path.insert(0, str(src_dir))
            return candidate
    raise FileNotFoundError("Could not find python_notebooks/src from script path.")


REPO_ROOT = _add_repo_src_to_path()

from mge_organoid_python.data_sources import resolve_data_root  # noqa: E402
from mge_organoid_python.notebook01_workflow import (  # noqa: E402
    Notebook01InputPaths,
    Notebook01InputSettings,
    Notebook01OutputPaths,
    Notebook01RunSettings,
)


def env_int(name: str, default: int) -> int:
    raw = os.environ.get(name)
    return default if raw is None or not raw.strip() else int(raw)


def main() -> None:
    data_root = resolve_data_root()
    notebook00_run_label = os.environ.get(
        "NOTEBOOK01_NOTEBOOK00_RUN_LABEL",
        "cellranger_filtered_manual_ec_div30_core_samples_freeze",
    )
    run_label = os.environ.get(
        "NOTEBOOK01_RUN_LABEL",
        f"{notebook00_run_label}_ccdifference_seurat_order_pcs10_neighbors20_leiden08_v1",
    )
    groupby = os.environ.get("NOTEBOOK01_MARKER_GROUPBY", "leiden")
    method = os.environ.get("NOTEBOOK01_MARKER_METHOD", "wilcoxon")
    top_n = env_int("NOTEBOOK01_MARKER_TOP_N", 25)
    rank_n_genes = env_int("NOTEBOOK01_MARKER_RANK_N_GENES", 200)
    restrict_to_branch_hvgs = os.environ.get("NOTEBOOK01_MARKER_RESTRICT_TO_BRANCH_HVGS", "0").strip() in {
        "1",
        "true",
        "TRUE",
        "yes",
    }

    input_settings = Notebook01InputSettings(notebook00_run_label=notebook00_run_label)
    run_settings = Notebook01RunSettings(run_label=run_label)
    input_paths = Notebook01InputPaths.from_data_root(data_root, settings=input_settings)
    output_paths = Notebook01OutputPaths.from_data_root(data_root, settings=run_settings)
    output_paths.ensure_dirs()

    branch_h5ad = output_paths.h5ad_dir / "combined" / "regressed_ccdifference" / "analysis_hvg_scaled_umap.h5ad"
    marker_table = output_paths.table_dir / f"notebook01_combined_{groupby}_markers_{method}_top{rank_n_genes}.tsv"
    marker_top_table = output_paths.table_dir / f"notebook01_combined_{groupby}_markers_{method}_top{top_n}.tsv"
    summary_table = output_paths.table_dir / f"notebook01_combined_{groupby}_marker_run_summary.tsv"
    plot_path = output_paths.plot_dir / "combined" / "regressed_ccdifference" / f"rank_genes_groups_{groupby}_{method}_top{top_n}.png"
    plot_path.parent.mkdir(parents=True, exist_ok=True)

    if not input_paths.combined_normalized_log1p.exists():
        raise FileNotFoundError(f"Missing normalized/log1p input: {input_paths.combined_normalized_log1p}")
    if not branch_h5ad.exists():
        raise FileNotFoundError(f"Missing final branch h5ad: {branch_h5ad}")

    print(f"Repo root: {REPO_ROOT}", flush=True)
    print(f"Data root: {data_root}", flush=True)
    print(f"Notebook 00 run label: {notebook00_run_label}", flush=True)
    print(f"Notebook 01 run label: {run_label}", flush=True)
    print(f"Groupby: {groupby}", flush=True)
    print(f"Method: {method}", flush=True)
    print(f"Input expression h5ad: {input_paths.combined_normalized_log1p}", flush=True)
    print(f"Final branch h5ad for clusters: {branch_h5ad}", flush=True)
    print(f"Restrict to branch HVGs: {restrict_to_branch_hvgs}", flush=True)

    branch = sc.read_h5ad(branch_h5ad, backed="r")
    try:
        if groupby not in branch.obs.columns:
            raise KeyError(f"Missing groupby column {groupby!r} in final branch h5ad.")
        branch_obs_names = branch.obs_names.copy()
        branch_var_names = branch.var_names.copy()
        groups = branch.obs[groupby].astype("category").copy()
    finally:
        branch.file.close()

    expr = sc.read_h5ad(input_paths.combined_normalized_log1p)
    gene_index = branch_var_names if restrict_to_branch_hvgs else expr.var_names
    expr = expr[branch_obs_names, gene_index].copy()
    expr.obs[groupby] = groups.loc[expr.obs_names].astype("category")
    expr.uns["notebook01_cluster_marker_context"] = {
        "run_label": run_label,
        "expression_matrix": "Notebook 00 normalized/log1p .X",
        "groupby_source": str(branch_h5ad),
        "groupby": groupby,
        "method": method,
        "restrict_to_branch_hvgs": restrict_to_branch_hvgs,
    }

    print(f"Marker AnnData shape: {expr.shape}", flush=True)
    print(f"Number of {groupby} groups: {expr.obs[groupby].nunique()}", flush=True)

    sc.tl.rank_genes_groups(
        expr,
        groupby=groupby,
        method=method,
        n_genes=min(rank_n_genes, expr.n_vars),
        use_raw=False,
    )
    marker_df = sc.get.rank_genes_groups_df(expr, group=None)
    marker_df.insert(0, "method", method)
    marker_df.insert(0, "groupby", groupby)
    marker_df.insert(0, "run_label", run_label)
    marker_df.to_csv(marker_table, sep="\t", index=False)

    top_df = marker_df.groupby("group", sort=False).head(top_n).copy()
    top_df.to_csv(marker_top_table, sep="\t", index=False)

    ax = sc.pl.rank_genes_groups(expr, n_genes=top_n, sharey=False, show=False)
    fig = ax[0].figure if isinstance(ax, list) else ax.figure
    fig.savefig(plot_path, dpi=160, bbox_inches="tight")

    summary = pd.DataFrame(
        [
            {
                "run_label": run_label,
                "input_expression_h5ad": str(input_paths.combined_normalized_log1p),
                "groupby_source_h5ad": str(branch_h5ad),
                "groupby": groupby,
                "method": method,
                "n_cells": int(expr.n_obs),
                "n_genes_ranked": int(expr.n_vars),
                "n_groups": int(expr.obs[groupby].nunique()),
                "rank_n_genes": int(min(rank_n_genes, expr.n_vars)),
                "top_n": int(top_n),
                "restrict_to_branch_hvgs": bool(restrict_to_branch_hvgs),
                "marker_table": str(marker_table),
                "marker_top_table": str(marker_top_table),
                "plot_path": str(plot_path),
            }
        ]
    )
    summary.to_csv(summary_table, sep="\t", index=False)

    print(f"Wrote marker table: {marker_table}", flush=True)
    print(f"Wrote top marker table: {marker_top_table}", flush=True)
    print(f"Wrote marker plot: {plot_path}", flush=True)
    print(f"Wrote summary: {summary_table}", flush=True)


if __name__ == "__main__":
    main()
