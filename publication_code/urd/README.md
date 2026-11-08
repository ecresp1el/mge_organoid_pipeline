# URD analysis and figures

- `urd_analysis.R`: cell and gene filtering, variable-gene selection, graph construction, flood pseudotime, and lineage reconstruction.
- `plot_urd.R`: DIV90 cluster, pseudotime, and marker-expression figures.
- `render_div90_figure.sh`: rendering command for the DIV90 composite figure.
- `figures/`: DIV90 URD composite PNG.

Dependencies: URD, Matrix, ggplot2, cowplot, and scales. The PCA-kNN option uses
RANN. SVG export can use svglite.

## Inputs

A gene-by-cell raw count matrix with unique gene and cell names, and cell
metadata with matching row names. The `urd_root_candidate` column contains
logical root flags. A separate column contains terminal-state labels. The order
of the supplied tip labels determines numeric tip IDs.

The DIV90 figure renderer uses a URD tree with cluster annotations and UMAP
coordinates in its metadata.

## Usage

Run from this folder:

```r
source("urd_analysis.R")
counts <- readRDS("counts.rds")
metadata <- read.delim("cell_metadata.tsv", row.names = 1, check.names = FALSE)

pt <- urd_pseudotime(counts, metadata, dm_method = "urd_expression")
tree <- urd_lineage(pt$object, annotation_col = "tip_label",
                    tip_labels = readLines("tip_labels.txt"),
                    pseudotime_name = pt$parameters$pseudotime_name)
saveRDS(tree$object, "lineage_tree.rds")
saveRDS(list(pseudotime = pt$parameters, lineage = tree$parameters,
             tip_mapping = tree$tip_mapping), "analysis_parameters.rds")
```

```sh
Rscript plot_urd.R --tree-rds lineage_tree.rds --outdir urd_figures \
  --annotation-col cluster_number_name --pseudotime-name root_pseudotime
```

For the DIV90 composite with identity and migration-guidance marker panels:

```sh
bash render_div90_figure.sh lineage_tree.rds urd_figures
```

## Analysis

Variable genes are ranked by variance in URD-normalized log expression.
Pseudotime is calculated by repeated flooding from the supplied root cells.
Lineage reconstruction uses pseudotime-weighted tip-to-root random walks and
visitation divergence.

The `urd_expression` option uses URD graph construction. The `pca_knn` option
uses a symmetric Gaussian-weighted graph in PC space; its stored coordinates
are scaled PCs rather than diffusion eigenvectors.

Root flags, terminal-state labels, cell inclusion, and analysis parameters are
supplied for the dataset being analyzed. For the DIV90 seven-tip design,
endpoints comprise five neuronal groups, astrocytes, and OPCs.

Marker display values use log1p(CP10K), with gene-specific colour limits.
The figure renderer applies DIV90 cluster recoding and writes plot exclusions,
colour limits, and figure parameters with its outputs.
