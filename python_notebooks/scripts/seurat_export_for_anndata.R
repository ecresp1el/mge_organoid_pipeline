#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_msg <- function(...) message("[R ", timestamp(), "] ", paste0(..., collapse = ""))

parse_args <- function(args) {
  out <- list(
    seurat = NULL,
    outdir = NULL,
    assay = "RNA",
    reduction = "umap",
    expression_layer = "data"
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

get_assay_matrix <- function(obj, assay, expression_layer) {
  if (expression_layer == "data") {
    layer_candidates <- c("data", "logcounts", "counts")
  } else {
    layer_candidates <- c(expression_layer, "data", "counts")
  }

  last_error <- NULL
  for (layer in unique(layer_candidates)) {
    mat <- tryCatch(
      SeuratObject::GetAssayData(obj, assay = assay, layer = layer),
      error = function(e) {
        last_error <<- e
        NULL
      }
    )
    if (!is.null(mat)) {
      return(list(matrix = mat, layer = layer))
    }
  }
  stop("Unable to extract assay matrix. Last error: ", conditionMessage(last_error), call. = FALSE)
}

safe_tsv_value <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (inherits(x, "POSIXt") || inherits(x, "Date")) x <- as.character(x)
  if (is.list(x)) x <- vapply(x, function(v) paste(as.character(v), collapse = ";"), character(1))
  x
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(opt$seurat) || !nzchar(opt$seurat)) stop("--seurat is required", call. = FALSE)
if (is.null(opt$outdir) || !nzchar(opt$outdir)) stop("--outdir is required", call. = FALSE)
if (!file.exists(opt$seurat)) stop("Seurat source does not exist: ", opt$seurat, call. = FALSE)

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

log_msg("R executable: ", R.home("bin"))
log_msg("R version: ", paste(R.version$major, R.version$minor, sep = "."))
log_msg("Seurat version: ", as.character(utils::packageVersion("Seurat")))
log_msg("SeuratObject version: ", as.character(utils::packageVersion("SeuratObject")))
log_msg("Reading Seurat RDS: ", opt$seurat)
obj <- readRDS(opt$seurat)
if (!inherits(obj, "Seurat")) stop("Object is not a Seurat object: ", opt$seurat, call. = FALSE)
log_msg("Loaded Seurat object with ", ncol(obj), " cells and ", nrow(obj), " features")

if (!(opt$assay %in% Seurat::Assays(obj))) {
  stop("Missing assay '", opt$assay, "'. Available assays: ", paste(Seurat::Assays(obj), collapse = ", "), call. = FALSE)
}
if (!(opt$reduction %in% Seurat::Reductions(obj))) {
  stop("Missing reduction '", opt$reduction, "'. Available reductions: ", paste(Seurat::Reductions(obj), collapse = ", "), call. = FALSE)
}

Seurat::DefaultAssay(obj) <- opt$assay
if (
  inherits(obj[[opt$assay]], "Assay5") &&
  exists("JoinLayers", where = asNamespace("SeuratObject"), inherits = FALSE)
) {
  log_msg("Assay5 detected; joining layers for assay: ", opt$assay)
  obj <- SeuratObject::JoinLayers(obj, assay = opt$assay)
}

log_msg("Extracting assay matrix for layer preference: ", opt$expression_layer)
mat_res <- get_assay_matrix(obj, opt$assay, opt$expression_layer)
mat <- mat_res$matrix
layer_used <- mat_res$layer
if (!inherits(mat, "sparseMatrix")) mat <- Matrix::Matrix(mat, sparse = TRUE)
log_msg("Using layer: ", layer_used)
log_msg("Matrix dimensions features x cells: ", nrow(mat), " x ", ncol(mat))

matrix_path <- file.path(opt$outdir, "matrix.mtx")
features_path <- file.path(opt$outdir, "features.tsv")
barcodes_path <- file.path(opt$outdir, "barcodes.tsv")
obs_path <- file.path(opt$outdir, "obs.tsv")
umap_path <- file.path(opt$outdir, "umap.tsv")
manifest_path <- file.path(opt$outdir, "manifest.tsv")

log_msg("Writing sparse Matrix Market: ", matrix_path)
Matrix::writeMM(mat, matrix_path)

log_msg("Writing feature names: ", features_path)
write.table(
  data.frame(feature_id = rownames(mat), stringsAsFactors = FALSE),
  features_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

log_msg("Writing cell barcodes: ", barcodes_path)
write.table(
  data.frame(cell_id = colnames(mat), stringsAsFactors = FALSE),
  barcodes_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

log_msg("Writing cell metadata: ", obs_path)
cell_ids <- colnames(mat)
obs <- obj@meta.data[cell_ids, , drop = FALSE]
obs <- as.data.frame(lapply(obs, safe_tsv_value), stringsAsFactors = FALSE)
obs <- data.frame(cell_id = cell_ids, obs, check.names = FALSE, stringsAsFactors = FALSE)
write.table(obs, obs_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

log_msg("Writing UMAP coordinates: ", umap_path)
emb <- Seurat::Embeddings(obj, reduction = opt$reduction)
missing_emb <- setdiff(colnames(mat), rownames(emb))
if (length(missing_emb) > 0) {
  stop("Reduction '", opt$reduction, "' is missing embeddings for ", length(missing_emb), " cells.", call. = FALSE)
}
emb <- emb[colnames(mat), seq_len(2), drop = FALSE]
umap <- data.frame(
  cell_id = rownames(emb),
  UMAP_1 = emb[, 1],
  UMAP_2 = emb[, 2],
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write.table(umap, umap_path, sep = "\t", quote = FALSE, row.names = FALSE)

log_msg("Writing manifest: ", manifest_path)
manifest <- data.frame(
  key = c("assay", "layer_used", "reduction", "n_features", "n_cells"),
  value = c(opt$assay, layer_used, opt$reduction, as.character(nrow(mat)), as.character(ncol(mat))),
  stringsAsFactors = FALSE
)
write.table(manifest, manifest_path, sep = "\t", quote = FALSE, row.names = FALSE)

log_msg("Finished Seurat export for AnnData: ", opt$outdir)
