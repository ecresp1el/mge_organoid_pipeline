#!/usr/bin/env python
"""Validate Notebook 01 cell-cycle genes against Notebook 00 checkpoints."""

from __future__ import annotations

import os
import sys
from pathlib import Path

import anndata as ad


def _add_repo_src_to_path() -> Path:
    start = Path(__file__).resolve()
    for candidate in [start.parent, *start.parents]:
        src_dir = candidate / "python_notebooks" / "src"
        if src_dir.exists():
            sys.path.insert(0, str(src_dir))
            return candidate
    raise FileNotFoundError("Could not find python_notebooks/src from script path.")


REPO_ROOT = _add_repo_src_to_path()

from mge_organoid_python.cell_cycle import select_cell_cycle_genes  # noqa: E402
from mge_organoid_python.data_sources import resolve_data_root  # noqa: E402
from mge_organoid_python.notebook01_workflow import (  # noqa: E402
    Notebook01InputPaths,
    Notebook01InputSettings,
    Notebook01OutputPaths,
    Notebook01RunSettings,
)


def main() -> None:
    data_root = resolve_data_root()
    notebook00_run_label = os.environ.get(
        "NOTEBOOK01_NOTEBOOK00_RUN_LABEL",
        "cellranger_filtered_manual_ec_div30_core_samples_freeze",
    )
    run_label = os.environ.get(
        "NOTEBOOK01_RUN_LABEL",
        f"{notebook00_run_label}_ccdifference_seurat_order_v1",
    )

    input_settings = Notebook01InputSettings(notebook00_run_label=notebook00_run_label)
    run_settings = Notebook01RunSettings(run_label=run_label)
    input_paths = Notebook01InputPaths.from_data_root(data_root, settings=input_settings)
    output_paths = Notebook01OutputPaths.from_data_root(data_root, settings=run_settings)
    output_paths.ensure_dirs()

    input_path = input_paths.combined_normalized_log1p
    if not input_path.exists():
        raise FileNotFoundError(f"Missing Notebook 00 checkpoint: {input_path}")

    print(f"Repo root: {REPO_ROOT}")
    print(f"Data root: {data_root}")
    print(f"Notebook 00 run label: {notebook00_run_label}")
    print(f"Notebook 01 run label: {run_label}")
    print(f"Input checkpoint: {input_path}")
    print(f"Output tables: {output_paths.table_dir}")

    backed = ad.read_h5ad(input_path, backed="r")
    try:
        var_names = list(map(str, backed.var_names))
    finally:
        backed.file.close()

    selection = select_cell_cycle_genes(var_names)
    gene_table = selection.gene_table
    present_table = gene_table.loc[gene_table["present_in_adata"]].copy()
    missing_table = gene_table.loc[~gene_table["present_in_adata"]].copy()
    summary_table = selection.summary_table.copy()
    summary_table.insert(0, "input_path", str(input_path))
    summary_table.insert(0, "run_label", run_label)
    summary_table.insert(0, "notebook00_run_label", notebook00_run_label)

    gene_table.to_csv(output_paths.table_dir / "cell_cycle_gene_source.tsv", sep="\t", index=False)
    present_table.to_csv(output_paths.table_dir / "cell_cycle_genes_present.tsv", sep="\t", index=False)
    missing_table.to_csv(output_paths.table_dir / "cell_cycle_genes_missing.tsv", sep="\t", index=False)
    summary_table.to_csv(output_paths.table_dir / "cell_cycle_gene_summary.tsv", sep="\t", index=False)

    print(summary_table.to_string(index=False))
    print(f"S genes present: {len(selection.s_genes)} / 43")
    print(f"G2M genes present: {len(selection.g2m_genes)} / 54")
    print(f"Cell-cycle genes present: {len(selection.cell_cycle_genes)} / 97")

    if len(selection.s_genes) < 10:
        raise AssertionError("Too few S-phase genes present for reliable scoring.")
    if len(selection.g2m_genes) < 10:
        raise AssertionError("Too few G2M-phase genes present for reliable scoring.")


if __name__ == "__main__":
    main()
