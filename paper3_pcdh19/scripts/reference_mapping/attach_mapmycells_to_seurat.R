#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat is required")

parse_args <- function(args) {
  values <- list()
  i <- 1L
  while (i <= length(args)) {
    if (!startsWith(args[[i]], "--") || i == length(args)) stop("Use --name value arguments")
    values[[substring(args[[i]], 3L)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  required <- c("seurat", "mapmycells", "output")
  missing <- required[!required %in% names(values)]
  if (length(missing)) stop("Missing: ", paste(missing, collapse = ", "))
  values
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
query <- readRDS(args$seurat)
if (!inherits(query, "Seurat")) stop("Input is not a Seurat object")
mapping <- read.delim(gzfile(args$mapmycells), check.names = FALSE)
if (anyDuplicated(mapping$cell_id) || !setequal(mapping$cell_id, colnames(query))) {
  stop("MapMyCells/Seurat exact cell identity failed")
}
mapping_columns <- grep("^mmc_", names(mapping), value = TRUE)
if (!length(mapping_columns)) stop("No MapMyCells columns found")
rownames(mapping) <- mapping$cell_id
query <- Seurat::AddMetaData(query, mapping[colnames(query), mapping_columns, drop = FALSE])
query@misc$MapMyCells <- list(
  role = "independent adult Allen WMB comparator; not E15 developmental ground truth",
  hierarchy = "all levels and confidence fields preserved with mmc_ prefixes",
  cell_identity_join = "exact stable sample_id_barcode; 450,788/450,788",
  reclustered = FALSE,
  integrated = FALSE,
  existing_umap_recomputed = FALSE
)

temporary <- paste0(args$output, ".partial")
on.exit(unlink(temporary), add = TRUE)
saveRDS(query, temporary, compress = FALSE)
validated <- readRDS(temporary)
if (!inherits(validated, "Seurat") || ncol(validated) != ncol(query) ||
    !all(mapping_columns %in% colnames(validated@meta.data))) {
  stop("Combined Seurat validation failed")
}
rm(validated)
if (!file.rename(temporary, args$output)) stop("Could not atomically publish combined Seurat RDS")
if (normalizePath(args$seurat, mustWork = TRUE) != normalizePath(args$output, mustWork = TRUE)) {
  unlink(args$seurat)
}
writeLines(c(
  "PASS",
  paste0("cells=", ncol(query)),
  paste0("mapmycells_columns=", length(mapping_columns)),
  "mind_fields_preserved=true",
  "reclustered=false",
  "integrated=false",
  "existing_umap_recomputed=false"
), file.path(dirname(args$output), "COMBINED_SEURAT_SUCCESS.txt"))
