# Cross-study marker expression

- `export_marker_expression.R`: expression and metadata export from Seurat objects.
- `plot_marker_expression.py`: marker-expression panels and sample summaries.
- `studies.example.tsv` and `genes.example.tsv`: input templates.
- `figures/`: cross-study marker-expression PNG figure.

Export dependencies: R with Seurat and Matrix. Plotting dependencies: Python
3.10 or newer with the packages in `requirements.txt`.

## Inputs

A study manifest containing object paths, assay and layer names, and metadata
columns, plus a one-column `gene` table. The example files show the required
format. Relative object paths are resolved against `--project-root`.

The plotter accepts cell tables containing `cell_id`, `study_id`, `study_label`,
`sample`, `cluster`, `umap_1`, `umap_2`, and one expression column per gene.

## Usage

Run from this folder:

```sh
Rscript export_marker_expression.R --study-table studies.tsv \
  --gene-table genes.example.tsv --project-root . --outdir marker_tables
python3 plot_marker_expression.py --input marker_tables/*_marker_expression.tsv.gz \
  --outdir marker_figures
```

Add `--preview` to produce 80 dpi PNGs. Standard output includes PNG, PDF, and SVG.

## Expression summaries

For raw counts, `log1p_cp10k` computes log1p(10000 × gene counts / total counts
across all genes). The `none` option uses the supplied normalized layer, which
must contain log1p(CP10K) values for comparison on this scale. Missing genes are
represented as missing values. Gene matches are exported in separate tables.

Samarasinghe samples containing `Ctrl` are included. DIV90 numeric clusters 6 and
7 are excluded as stressed cells. Cells with nonfinite embedding coordinates are
excluded before summaries.

Percent-positive values use finite gene-expression values within each included
sample as their denominator, with expression >= 0.5 by default. Violin rendering
uses at most 1,800 cells per sample; summaries use all included finite values.
Replicate panels show per-sample means using all finite values or per-cell
cutoffs of 0.5 and 1. Replicate panels report sample-level summaries.

Outputs include expression panels, sample and study summaries, cell-filter
counts, and figure parameters.
