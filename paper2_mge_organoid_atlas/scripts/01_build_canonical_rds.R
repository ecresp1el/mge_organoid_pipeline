#!/usr/bin/env Rscript

# Build one current minimal Seurat object from an audited Seurat source.
# The saved RDS contains only an RNA assay, selected metadata, and compact
# provenance: no reductions, graphs, neighbors, commands, tools, SCT models,
# integrated assays, images, or scaled matrix. Matrix Market/TSV files are
# temporary bridges for the independently validated H5AD writer and are
# removed by that writer after success.

suppressPackageStartupMessages({
  library(digest)
  library(Matrix)
  library(Seurat)
  library(SeuratObject)
})

options(stringsAsFactors = FALSE, digits = 17)

log_msg <- function(...) {
  message("[canonical-rds ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(..., collapse = ""))
}

parse_args <- function(args) {
  out <- list(manifest = NULL, `study-id` = NULL, outdir = NULL)
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
  for (name in names(out)) {
    if (is.null(out[[name]]) || !nzchar(out[[name]])) stop("--", name, " is required", call. = FALSE)
  }
  out
}

split_semicolon <- function(x) {
  values <- trimws(strsplit(as.character(x), ";", fixed = TRUE)[[1]])
  unique(values[nzchar(values)])
}

safe_metadata_frame <- function(df) {
  out <- as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)
  for (name in colnames(out)) {
    value <- out[[name]]
    if (is.factor(value)) value <- as.character(value)
    if (inherits(value, "POSIXt") || inherits(value, "Date")) value <- as.character(value)
    if (is.list(value)) value <- vapply(value, function(x) paste(as.character(x), collapse = ";"), character(1))
    out[[name]] <- value
  }
  out
}

get_layer <- function(obj, assay, layer) {
  tryCatch(
    SeuratObject::LayerData(obj, assay = assay, layer = layer),
    error = function(e) NULL
  )
}

as_dgc <- function(x) {
  if (inherits(x, "dgCMatrix")) return(x)
  as(Matrix::Matrix(x, sparse = TRUE), "dgCMatrix")
}

write_tsv <- function(df, path) {
  write.table(
    df,
    path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )
}

matrix_summary <- function(mat, matrix_name) {
  x <- mat@x
  data.frame(
    matrix = matrix_name,
    n_features = nrow(mat),
    n_cells = ncol(mat),
    nnz = length(x),
    value_sum = if (length(x)) sum(x) else 0,
    value_sum_squares = if (length(x)) sum(x * x) else 0,
    value_min_nonzero = if (length(x)) min(x) else NA_real_,
    value_max = if (length(x)) max(x) else 0,
    stringsAsFactors = FALSE
  )
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
manifest <- read.delim(opt$manifest, check.names = FALSE, stringsAsFactors = FALSE)
row <- manifest[manifest$study_id == opt$`study-id`, , drop = FALSE]
if (nrow(row) != 1L) stop("Expected exactly one manifest row for ", opt$`study-id`, call. = FALSE)

study_id <- row$study_id[[1]]
display_name <- row$display_name[[1]]
source_rds <- row$source_rds[[1]]
source_sha256 <- row$source_sha256[[1]]
assay <- row$rna_assay[[1]]
normalized_policy <- row$normalized_policy[[1]]
sample_column <- row$sample_column[[1]]
metadata_columns <- split_semicolon(row$metadata_columns[[1]])
outdir <- normalizePath(opt$outdir, mustWork = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
bridge_dir <- file.path(outdir, ".bridge")
dir.create(bridge_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(source_rds)) stop("Missing source RDS: ", source_rds, call. = FALSE)
observed_sha <- digest::digest(source_rds, algo = "sha256", file = TRUE, serialize = FALSE)
if (!identical(observed_sha, source_sha256)) {
  stop("Source SHA-256 mismatch for ", study_id, ": expected ", source_sha256, ", observed ", observed_sha, call. = FALSE)
}

log_msg("Loading ", study_id, " from ", source_rds)
source <- readRDS(source_rds)
if (!inherits(source, "Seurat")) stop("Source is not a Seurat object", call. = FALSE)
if (!(assay %in% Seurat::Assays(source))) stop("Source is missing assay ", assay, call. = FALSE)

if (inherits(source[[assay]], "Assay5")) {
  log_msg("Joining sample-split ", assay, " layers in memory")
  source <- SeuratObject::JoinLayers(source, assay = assay)
}

counts <- get_layer(source, assay, "counts")
if (is.null(counts)) stop("Required raw counts layer is absent for ", study_id, call. = FALSE)
counts <- as_dgc(counts)
if (any(!is.finite(counts@x)) || any(counts@x < 0)) stop("Counts contain invalid values", call. = FALSE)
if (any(abs(counts@x - round(counts@x)) > 1e-8)) stop("Counts are not integer-valued", call. = FALSE)

normalized <- NULL
normalized_source_layer <- "absent"
if (identical(normalized_policy, "preserve_if_present")) {
  candidate <- get_layer(source, assay, "data")
  if (!is.null(candidate) && nrow(candidate) > 0L && ncol(candidate) > 0L) {
    normalized <- as_dgc(candidate)
    normalized_source_layer <- "RNA:data"
  }
} else if (!identical(normalized_policy, "none")) {
  stop("Unsupported normalized_policy: ", normalized_policy, call. = FALSE)
}

if (!is.null(normalized)) {
  if (!identical(dim(normalized), dim(counts))) stop("Normalized/count dimensions differ", call. = FALSE)
  if (!identical(rownames(normalized), rownames(counts))) stop("Normalized/count feature IDs differ", call. = FALSE)
  if (!identical(colnames(normalized), colnames(counts))) stop("Normalized/count cell IDs differ", call. = FALSE)
  if (any(!is.finite(normalized@x))) stop("Normalized matrix contains non-finite values", call. = FALSE)
}

source_cell_id <- colnames(counts)
source_feature_id <- rownames(counts)
if (is.null(source_cell_id) || anyNA(source_cell_id) || any(!nzchar(source_cell_id)) || anyDuplicated(source_cell_id)) {
  stop("Source cell IDs must be nonempty and unique", call. = FALSE)
}
if (is.null(source_feature_id) || anyNA(source_feature_id) || any(!nzchar(source_feature_id)) || anyDuplicated(source_feature_id)) {
  stop("Source feature IDs must be nonempty and unique", call. = FALSE)
}
canonical_cell_id <- paste0(study_id, "::", source_cell_id)
colnames(counts) <- canonical_cell_id
if (!is.null(normalized)) colnames(normalized) <- canonical_cell_id

missing_metadata <- setdiff(c(metadata_columns, sample_column), colnames(source@meta.data))
if (length(missing_metadata)) {
  stop("Required metadata columns are absent: ", paste(missing_metadata, collapse = ", "), call. = FALSE)
}
cell_metadata <- source@meta.data[source_cell_id, metadata_columns, drop = FALSE]
cell_metadata <- safe_metadata_frame(cell_metadata)
cell_metadata <- data.frame(
  canonical_cell_id = canonical_cell_id,
  source_cell_id = source_cell_id,
  study_id = study_id,
  study_label = display_name,
  canonical_sample_id = as.character(source@meta.data[source_cell_id, sample_column]),
  cell_metadata,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
rownames(cell_metadata) <- canonical_cell_id

source_feature_metadata <- tryCatch(source[[assay]][[]], error = function(e) NULL)
if (is.null(source_feature_metadata) || nrow(source_feature_metadata) != length(source_feature_id)) {
  source_feature_metadata <- data.frame(row.names = source_feature_id)
} else {
  source_feature_metadata <- source_feature_metadata[source_feature_id, , drop = FALSE]
}
source_feature_metadata <- safe_metadata_frame(source_feature_metadata)
identifier_patterns <- c(
  "gene_id", "gene_ids", "gene_symbol", "symbol", "feature_id",
  "feature_name", "feature_type", "genome", "ensembl", "ensembl_id"
)
keep_feature_metadata <- tolower(colnames(source_feature_metadata)) %in% identifier_patterns
source_feature_metadata <- source_feature_metadata[, keep_feature_metadata, drop = FALSE]
generated_feature_columns <- c("canonical_feature_id", "source_feature_id", "feature_symbol")
colliding_feature_columns <- colnames(source_feature_metadata) %in% generated_feature_columns
colnames(source_feature_metadata)[colliding_feature_columns] <- paste0(
  "source_metadata_",
  colnames(source_feature_metadata)[colliding_feature_columns]
)
feature_metadata <- data.frame(
  canonical_feature_id = source_feature_id,
  source_feature_id = source_feature_id,
  feature_symbol = source_feature_id,
  source_feature_metadata,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
rownames(feature_metadata) <- source_feature_id

canonical_provenance <- list(
  schema_name = "paper2_mge_organoid_canonical_seurat",
  schema_version = "1.0.0",
  study_id = study_id,
  study_label = display_name,
  source = list(
    source_rds = source_rds,
    source_sha256 = source_sha256,
    source_assay = assay,
    counts_source_layer = "RNA:counts",
    normalized_source_layer = normalized_source_layer,
    normalized_policy = normalized_policy,
    sample_source_column = sample_column,
    cells_subsetted = FALSE,
    features_subsetted = FALSE,
    normalization_performed = FALSE,
    gene_harmonization_performed = FALSE,
    metadata_harmonization_performed = FALSE,
    biological_content_changed = FALSE
  ),
  software = list(
    R = R.version.string,
    Seurat = as.character(utils::packageVersion("Seurat")),
    SeuratObject = as.character(utils::packageVersion("SeuratObject")),
    Matrix = as.character(utils::packageVersion("Matrix"))
  )
)

log_msg("Creating clean current Seurat object")
canonical <- Seurat::CreateSeuratObject(
  counts = counts,
  assay = "RNA",
  project = study_id,
  meta.data = cell_metadata,
  min.cells = 0,
  min.features = 0
)
if (!is.null(normalized)) {
  SeuratObject::LayerData(canonical, assay = "RNA", layer = "data") <- normalized
}
canonical[["RNA"]] <- SeuratObject::AddMetaData(canonical[["RNA"]], metadata = feature_metadata)
Seurat::VariableFeatures(canonical[["RNA"]]) <- character(0)
canonical@meta.data <- cell_metadata[colnames(canonical), , drop = FALSE]
canonical@active.ident <- factor(rep(study_id, ncol(canonical)), levels = study_id)
names(canonical@active.ident) <- colnames(canonical)
canonical@assays <- canonical@assays["RNA"]
canonical@reductions <- list()
canonical@graphs <- list()
canonical@neighbors <- list()
canonical@images <- list()
canonical@commands <- list()
canonical@tools <- list()
canonical@misc <- list(canonical_provenance = canonical_provenance)
Seurat::DefaultAssay(canonical) <- "RNA"

if (!inherits(canonical, "Seurat")) stop("Canonical RDS is not a Seurat object", call. = FALSE)
if (!identical(Seurat::Assays(canonical), "RNA")) stop("Canonical object must contain only RNA", call. = FALSE)
if (length(Seurat::Reductions(canonical)) || length(SeuratObject::Graphs(canonical)) || length(canonical@neighbors)) {
  stop("Canonical object retained reductions, graphs, or neighbors", call. = FALSE)
}
if (length(canonical@commands) || length(canonical@tools) || length(canonical@images)) {
  stop("Canonical object retained commands, tools, or images", call. = FALSE)
}
canonical_layers <- SeuratObject::Layers(canonical[["RNA"]])
if (any(grepl("^scale\\.data", canonical_layers))) stop("Canonical object retained scale.data", call. = FALSE)
if (any(grepl("integrated|SCT", Seurat::Assays(canonical), ignore.case = TRUE))) {
  stop("Canonical object retained integrated/SCT assay", call. = FALSE)
}

rds_path <- file.path(outdir, paste0(study_id, "_minimal.rds"))
log_msg("Saving minimal RDS: ", rds_path)
saveRDS(canonical, rds_path, compress = "gzip", version = 3)

# Reload before producing the bridge so H5AD is constructed from the exact RDS
# representation rather than from a parallel in-memory branch.
rm(canonical, counts, normalized, source)
invisible(gc())
canonical <- readRDS(rds_path)
if (!inherits(canonical, "Seurat")) stop("Reloaded RDS is not a Seurat object", call. = FALSE)
counts <- as_dgc(SeuratObject::LayerData(canonical, assay = "RNA", layer = "counts"))
reloaded_layers <- SeuratObject::Layers(canonical[["RNA"]])
normalized <- if ("data" %in% reloaded_layers) {
  as_dgc(SeuratObject::LayerData(canonical, assay = "RNA", layer = "data"))
} else {
  NULL
}
cell_metadata <- safe_metadata_frame(canonical@meta.data[colnames(canonical), , drop = FALSE])
feature_metadata <- safe_metadata_frame(canonical[["RNA"]][[]][rownames(canonical), , drop = FALSE])

log_msg("Writing temporary counts bridge")
Matrix::writeMM(counts, file.path(bridge_dir, "counts.mtx"))
if (!is.null(normalized)) {
  log_msg("Writing temporary normalized-expression bridge")
  Matrix::writeMM(normalized, file.path(bridge_dir, "lognorm.mtx"))
}
write_tsv(cell_metadata, file.path(bridge_dir, "obs.tsv"))
write_tsv(feature_metadata, file.path(bridge_dir, "var.tsv"))
metadata_schema <- data.frame(
  table = c(rep("obs", ncol(cell_metadata)), rep("var", ncol(feature_metadata))),
  column = c(colnames(cell_metadata), colnames(feature_metadata)),
  r_class = c(
    vapply(cell_metadata, function(x) class(x)[[1]], character(1)),
    vapply(feature_metadata, function(x) class(x)[[1]], character(1))
  ),
  stringsAsFactors = FALSE
)
write_tsv(metadata_schema, file.path(bridge_dir, "metadata_schema.tsv"))

matrix_rows <- matrix_summary(counts, "counts")
if (!is.null(normalized)) matrix_rows <- rbind(matrix_rows, matrix_summary(normalized, "lognorm"))
write_tsv(matrix_rows, file.path(outdir, "matrix_summary_rds.tsv"))

manifest_out <- data.frame(
  key = c(
    "schema_name", "schema_version", "object_class", "study_id", "study_label", "source_rds",
    "source_sha256", "source_assay", "counts_source_layer", "normalized_source_layer",
    "n_cells", "n_features", "counts_nnz", "lognorm_present", "metadata_columns",
    "assays_present", "layers_present", "reductions_present", "graphs_present",
    "neighbors_present", "commands_present", "tools_present", "images_present",
    "scale_data_present", "integrated_assays_present", "R_version", "Seurat_version",
    "SeuratObject_version"
  ),
  value = c(
    "paper2_mge_organoid_canonical_seurat", "1.0.0", "Seurat", study_id, display_name, source_rds,
    source_sha256, assay, "RNA:counts", normalized_source_layer,
    ncol(counts), nrow(counts), length(counts@x), !is.null(normalized),
    paste(colnames(cell_metadata), collapse = ";"), "RNA", paste(reloaded_layers, collapse = ";"),
    "false", "false", "false", "false", "false", "false", "false", "false",
    R.version.string, as.character(utils::packageVersion("Seurat")),
    as.character(utils::packageVersion("SeuratObject"))
  ),
  stringsAsFactors = FALSE
)
write_tsv(manifest_out, file.path(outdir, "canonical_manifest.tsv"))
writeLines(capture.output(utils::sessionInfo()), file.path(outdir, "r_sessionInfo.txt"))
log_msg("RDS and bridge complete for ", study_id)
