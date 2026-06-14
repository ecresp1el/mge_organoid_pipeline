#!/usr/bin/env python3
"""Plot requested DIV90 genes on the exported Seurat UMAP from cached H5AD."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pandas as pd


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _add_src_to_path() -> None:
    src = _repo_root() / "python_notebooks" / "src"
    if str(src) not in sys.path:
        sys.path.insert(0, str(src))


def _split_csv(value: str) -> list[str]:
    return [part.strip() for part in value.split(",") if part.strip()]


def _default_project_root() -> Path:
    return Path(os.environ.get("PROJECT_ROOT", "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder"))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create DIV90 gene-expression UMAP plots from the cached AnnData file. "
            "This does not submit jobs or load the source Seurat RDS."
        )
    )
    parser.add_argument(
        "--genes",
        required=True,
        help="Comma-separated genes to plot, for example: LHX8,ISL1,LHX6,CRABP1,ANGPT2",
    )
    parser.add_argument(
        "--outdir",
        default=None,
        help="Output directory. Default: PROJECT_ROOT/results/div90_gene_umaps_from_h5ad/<auto_label>",
    )
    parser.add_argument(
        "--project-root",
        default=str(_default_project_root()),
        help="Project/data root containing results/python_anndata/varela_div90.h5ad.",
    )
    parser.add_argument(
        "--h5ad",
        default=None,
        help="Optional explicit H5AD path. Default: PROJECT_ROOT/results/python_anndata/varela_div90.h5ad.",
    )
    parser.add_argument(
        "--max-cells",
        type=int,
        default=None,
        help="Optional downsample for plotting only. By default all cells are plotted.",
    )
    parser.add_argument(
        "--vmax-quantile",
        type=float,
        default=0.99,
        help="Positive-expression quantile used as each gene's color scale maximum.",
    )
    parser.add_argument(
        "--point-size",
        type=float,
        default=0.36,
        help="Point size for expressing cells.",
    )
    parser.add_argument(
        "--background-point-size",
        type=float,
        default=0.15,
        help="Point size for all-cell grey background.",
    )
    return parser.parse_args()


def main() -> int:
    _add_src_to_path()
    from mge_organoid_python.cross_study_marker_expression import (
        CrossStudyMarkerSpec,
        extract_marker_expression_from_h5ad,
        plot_marker_umap_grid,
    )

    args = parse_args()
    genes = _split_csv(args.genes)
    if not genes:
        raise SystemExit("--genes did not contain any gene symbols")

    project_root = Path(args.project_root).expanduser().resolve()
    h5ad_path = Path(args.h5ad).expanduser().resolve() if args.h5ad else project_root / "results/python_anndata/varela_div90.h5ad"
    if not h5ad_path.exists():
        raise FileNotFoundError(f"Missing DIV90 H5AD: {h5ad_path}")

    safe_label = "_".join(gene.replace("/", "_") for gene in genes[:8])
    if len(genes) > 8:
        safe_label += f"_plus{len(genes) - 8}"
    outdir = (
        Path(args.outdir).expanduser().resolve()
        if args.outdir
        else project_root / "results/div90_gene_umaps_from_h5ad" / safe_label
    )
    table_dir = outdir / "tables"
    plot_dir = outdir / "plots"
    table_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    spec = CrossStudyMarkerSpec(
        study_id="varela_div90",
        study_label="This Study, DIV 90",
        seurat_path="/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/umap_props_output/clustered_day90_with_cluster_names_2.rds",
        h5ad_path=str(h5ad_path),
        sample_col="orig.ident",
        cluster_col="seurat_clusters",
    )

    marker_table = table_dir / "div90_requested_gene_expression_umap_table.tsv.gz"
    expression, gene_matches = extract_marker_expression_from_h5ad(
        spec=spec,
        output_path=marker_table,
        project_root=project_root,
        genes=genes,
        obsm_keys=("X_umap_seurat", "X_umap"),
    )
    gene_match_path = table_dir / "div90_requested_gene_match_table.tsv"
    gene_matches.to_csv(gene_match_path, sep="\t", index=False)

    grid_path = plot_dir / "div90_requested_gene_umap_grid.png"
    manifest = plot_marker_umap_grid(
        data=expression,
        output_path=grid_path,
        genes=genes,
        specs=[spec],
        title="DIV90 requested gene expression on exported Seurat UMAP",
        max_cells_per_study=args.max_cells,
        point_size=args.point_size,
        background_point_size=args.background_point_size,
        vmax_quantile=args.vmax_quantile,
    )
    manifest_path = table_dir / "div90_requested_gene_umap_plot_manifest.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)

    if "matched" in gene_matches:
        missing = gene_matches.loc[~gene_matches["matched"].astype(bool), "gene"].tolist()
    elif "matched_feature" in gene_matches:
        missing = gene_matches.loc[~gene_matches["matched_feature"].astype(str).str.len().astype(bool), "gene"].tolist()
    else:
        missing = []
    summary = pd.DataFrame(
        [
            {
                "h5ad_path": str(h5ad_path),
                "outdir": str(outdir),
                "n_cells": int(expression.shape[0]),
                "genes_requested": ",".join(genes),
                "genes_missing": ",".join(missing),
                "umap_plot_png": str(grid_path),
                "umap_plot_pdf": str(grid_path.with_suffix(".pdf")),
                "expression_table": str(marker_table),
                "gene_match_table": str(gene_match_path),
                "plot_manifest": str(manifest_path),
            }
        ]
    )
    summary_path = table_dir / "div90_requested_gene_umap_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)

    print(f"Wrote plot: {grid_path}")
    print(f"Wrote PDF:  {grid_path.with_suffix('.pdf')}")
    print(f"Wrote gene match table: {gene_match_path}")
    print(f"Wrote expression table: {marker_table}")
    if missing:
        print("Missing requested genes:", ", ".join(missing))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
