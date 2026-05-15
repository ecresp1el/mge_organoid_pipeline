#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(zellkonverter)
})

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_msg <- function(...) message("[R ", timestamp(), "] ", paste0(..., collapse = ""))

parse_args <- function(args) {
  out <- list(
    seurat = NULL,
    h5ad = NULL,
    assay = "RNA",
    reduction = "umap",
    expression_layer = "data",
    overwrite = "false"
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unknown argument: ", key, call. = FALSE)
    name <- substring(key, 3L)
    if (!(name %in% names(out))) stop("Unknown argument: ", key, call. = FALSE)
    if (i == length(args)) stop("Missing value for argument: ", key, call. = FALSE)
    out[[name]] <- args[[i + 1L]]
    i <- i + 2L
  }
  out
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
overwrite <- tolower(opt$overwrite) %in% c("1", "true", "yes", "y")

if (is.null(opt$seurat) || !nzchar(opt$seurat)) stop("--seurat is required", call. = FALSE)
if (is.null(opt$h5ad) || !nzchar(opt$h5ad)) stop("--h5ad is required", call. = FALSE)
if (!file.exists(opt$seurat)) stop("Seurat source does not exist: ", opt$seurat, call. = FALSE)

if (file.exists(opt$h5ad) && !overwrite) {
  log_msg("Existing H5AD found; skipping conversion: ", opt$h5ad)
  quit(save = "no", status = 0)
}

log_msg("R executable: ", R.home("bin"))
log_msg("R version: ", paste(R.version$major, R.version$minor, sep = "."))
log_msg("Seurat version: ", as.character(utils::packageVersion("Seurat")))
log_msg("SeuratObject version: ", as.character(utils::packageVersion("SeuratObject")))
log_msg("zellkonverter version: ", as.character(utils::packageVersion("zellkonverter")))
log_msg("RETICULATE_PYTHON: ", Sys.getenv("RETICULATE_PYTHON", unset = "<unset>"))
log_msg("RETICULATE_AUTOCONFIGURE: ", Sys.getenv("RETICULATE_AUTOCONFIGURE", unset = "<unset>"))

if (!nzchar(Sys.getenv("RETICULATE_PYTHON", unset = ""))) {
  stop(
    "RETICULATE_PYTHON is not set. Refusing to let reticulate auto-install Python.",
    call. = FALSE
  )
}

log_msg("Reading Seurat RDS: ", opt$seurat)
obj <- readRDS(opt$seurat)
if (!inherits(obj, "Seurat")) {
  stop("Object is not a Seurat object: ", opt$seurat, call. = FALSE)
}
log_msg("Loaded Seurat object with ", ncol(obj), " cells and ", nrow(obj), " features")

if (!(opt$assay %in% Seurat::Assays(obj))) {
  stop(
    "Missing assay '", opt$assay, "'. Available assays: ",
    paste(Seurat::Assays(obj), collapse = ", "),
    call. = FALSE
  )
}
if (!(opt$reduction %in% Seurat::Reductions(obj))) {
  stop(
    "Missing reduction '", opt$reduction, "'. Available reductions: ",
    paste(Seurat::Reductions(obj), collapse = ", "),
    call. = FALSE
  )
}

Seurat::DefaultAssay(obj) <- opt$assay
if (
  inherits(obj[[opt$assay]], "Assay5") &&
  exists("JoinLayers", where = asNamespace("SeuratObject"), inherits = FALSE)
) {
  log_msg("Assay5 detected; joining layers for assay: ", opt$assay)
  obj <- SeuratObject::JoinLayers(obj, assay = opt$assay)
}

log_msg("Converting Seurat object to SingleCellExperiment")
sce <- Seurat::as.SingleCellExperiment(obj, assay = opt$assay)

log_msg("Transferring reduction to reducedDim X_umap: ", opt$reduction)
emb <- Seurat::Embeddings(obj, reduction = opt$reduction)
if (is.null(dim(emb)) || nrow(emb) == 0 || ncol(emb) < 2) {
  stop("Reduction '", opt$reduction, "' does not contain a two-dimensional embedding.", call. = FALSE)
}
missing_emb <- setdiff(colnames(sce), rownames(emb))
if (length(missing_emb) > 0) {
  stop(
    "Reduction '", opt$reduction, "' is missing embeddings for ",
    length(missing_emb), " cells.",
    call. = FALSE
  )
}
emb <- emb[colnames(sce), , drop = FALSE]
emb <- as.matrix(emb[, seq_len(2), drop = FALSE])
colnames(emb) <- c("UMAP_1", "UMAP_2")
SingleCellExperiment::reducedDim(sce, "X_umap") <- emb

assay_names <- names(SummarizedExperiment::assays(sce))
x_name <- NULL
if (identical(opt$expression_layer, "data") && "logcounts" %in% assay_names) {
  x_name <- "logcounts"
} else if (opt$expression_layer %in% assay_names) {
  x_name <- opt$expression_layer
} else if ("counts" %in% assay_names) {
  x_name <- "counts"
} else if (length(assay_names) > 0) {
  x_name <- assay_names[[1]]
}

dir.create(dirname(opt$h5ad), showWarnings = FALSE, recursive = TRUE)
log_msg("Writing H5AD: ", opt$h5ad, " (X_name=", x_name, ")")
zellkonverter::writeH5AD(sce, file = opt$h5ad, X_name = x_name)
log_msg("Finished writing H5AD: ", opt$h5ad)
