#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

addon_library <- Sys.getenv("R_LIBS_USER", unset = "")
if (nzchar(addon_library)) .libPaths(c(addon_library, .libPaths()))

required_packages <- c("Matrix", "Seurat", "hdf5r")
for (package in required_packages) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " is required")
}

parse_args <- function(args) {
  values <- list()
  i <- 1L
  while (i <= length(args)) {
    if (!startsWith(args[[i]], "--") || i == length(args)) stop("Use --name value arguments")
    values[[substring(args[[i]], 3L)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  required <- c("cellranger-root", "sample-key", "query-metadata", "reference-matrix",
                "reference-join", "output-dir", "threshold", "seed", "cpus")
  missing <- required[!required %in% names(values)]
  if (length(missing) > 0L) stop("Missing: ", paste(missing, collapse = ", "))
  values
}

read_nested_rds <- function(path) {
  payload <- readBin(path, what = "raw", n = file.info(path)$size)
  layers <- 0L
  while (length(payload) >= 2L && identical(as.integer(payload[1:2]), c(31L, 139L))) {
    payload <- memDecompress(payload, type = "gzip")
    layers <- layers + 1L
  }
  connection <- rawConnection(payload, open = "rb")
  on.exit(close(connection), add = TRUE)
  list(object = readRDS(connection), gzip_layers = layers)
}

write_tsv_gz <- function(data, path) {
  temporary <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  connection <- gzfile(temporary, open = "wt")
  write.table(data, connection, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  close(connection)
  if (!file.rename(temporary, path)) stop("Could not publish ", path)
}

prefix_transfer <- function(transfer, prefix) {
  names(transfer) <- paste0(prefix, names(transfer))
  transfer
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
set.seed(as.integer(args$seed))
options(future.globals.maxSize = 200 * 1024^3)
if (requireNamespace("future", quietly = TRUE)) future::plan("sequential")
threshold <- as.numeric(args$threshold)
if (!is.finite(threshold) || threshold < 0 || threshold > 1) stop("Invalid threshold")
output_dir <- normalizePath(args[["output-dir"]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- read.delim(gzfile(args[["query-metadata"]]), check.names = FALSE)
if (anyDuplicated(metadata$cell_id)) stop("Query cell IDs are not unique")
rownames(metadata) <- metadata$cell_id
sample_key <- read.csv(args[["sample-key"]], check.names = FALSE)

objects <- vector("list", nrow(sample_key))
for (i in seq_len(nrow(sample_key))) {
  sample_id <- sample_key$technical_sample_id[[i]]
  matrix_path <- file.path(args[["cellranger-root"]], "per_sample_outs", sample_id,
                           "sample_filtered_feature_bc_matrix.h5")
  counts <- Seurat::Read10X_h5(matrix_path, use.names = TRUE, unique.features = TRUE)
  if (is.list(counts)) counts <- counts[["Gene Expression"]]
  stable_ids <- paste0(sample_id, "_", colnames(counts))
  colnames(counts) <- stable_ids
  if (!identical(stable_ids, metadata[stable_ids, "cell_id"])) {
    stop("Exact query metadata/count order failed for ", sample_id)
  }
  objects[[i]] <- Seurat::CreateSeuratObject(
    counts = counts,
    project = "paper3_e15_mge",
    meta.data = metadata[stable_ids, setdiff(names(metadata), "cell_id"), drop = FALSE],
    min.cells = 0,
    min.features = 0
  )
  rm(counts)
  invisible(gc())
}
query <- if (length(objects) == 1L) objects[[1L]] else merge(objects[[1L]], y = objects[-1L])
rm(objects)
invisible(gc())
if (!setequal(colnames(query), metadata$cell_id)) stop("Merged query cell identity failed")

coordinates <- as.matrix(metadata[colnames(query), c("vendor_umap_1", "vendor_umap_2")])
storage.mode(coordinates) <- "double"
colnames(coordinates) <- c("cellrangerUMAP_1", "cellrangerUMAP_2")
query[["cellranger.umap"]] <- Seurat::CreateDimReducObject(
  embeddings = coordinates, key = "cellrangerUMAP_", assay = "RNA"
)
query@misc$coordinate_warning <- paste(
  "These Cell Ranger UMAP coordinates were computed independently per sample.",
  "Use only within-sample or faceted displays."
)
query@misc$clustering_warning <- paste(
  "existing_cluster is sample-qualified because Cell Ranger graph clusters were",
  "computed independently per sample. No clustering was recomputed."
)

loaded <- read_nested_rds(args[["reference-matrix"]])
reference_counts <- loaded$object
if (!(is.matrix(reference_counts) || inherits(reference_counts, "Matrix"))) {
  stop("Bandler CA301 reference did not resolve to a matrix")
}
join <- read.delim(args[["reference-join"]], check.names = FALSE)
if (nrow(join) != 4481L || any(join$assignment_status != "definitive")) {
  stop("Expected exactly 4,481 definitive CA301 MIND joins")
}
if (anyDuplicated(join$cell_id) || !all(join$cell_id %in% colnames(reference_counts))) {
  stop("Bandler reference join/count identity failed")
}
reference_counts <- reference_counts[, join$cell_id, drop = FALSE]
reference_meta <- data.frame(
  MIND_class = join$later_atlas_class,
  MIND_cluster = join$later_atlas_cluster,
  row.names = join$cell_id,
  check.names = FALSE
)
reference <- Seurat::CreateSeuratObject(
  counts = reference_counts,
  project = "Bandler_CA301_E15.5_MIND",
  meta.data = reference_meta,
  min.cells = 0,
  min.features = 0
)
rm(reference_counts, loaded)
invisible(gc())

reference <- Seurat::NormalizeData(reference, verbose = FALSE)
reference <- Seurat::FindVariableFeatures(reference, selection.method = "vst", nfeatures = 2000,
                                           verbose = FALSE)
transfer_features <- intersect(Seurat::VariableFeatures(reference), rownames(query))
if (length(transfer_features) < 500L) stop("Too few shared variable features: ", length(transfer_features))
reference <- Seurat::ScaleData(reference, features = transfer_features, verbose = FALSE)
reference <- Seurat::RunPCA(reference, features = transfer_features, npcs = 30,
                            seed.use = as.integer(args$seed), verbose = FALSE)
query <- Seurat::NormalizeData(query, verbose = FALSE)
anchors <- Seurat::FindTransferAnchors(
  reference = reference,
  query = query,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  reduction = "pcaproject",
  features = transfer_features,
  dims = 1:30,
  verbose = TRUE
)

class_transfer <- Seurat::TransferData(
  anchorset = anchors, refdata = reference$MIND_class, dims = 1:30, verbose = TRUE
)
cluster_transfer <- Seurat::TransferData(
  anchorset = anchors, refdata = reference$MIND_cluster, dims = 1:30, verbose = TRUE
)
class_transfer <- prefix_transfer(class_transfer, "mind_class_")
cluster_transfer <- prefix_transfer(cluster_transfer, "mind_cluster_")
query <- Seurat::AddMetaData(query, class_transfer)
query <- Seurat::AddMetaData(query, cluster_transfer)

query$MIND_class_raw <- query$mind_class_predicted.id
query$MIND_cluster_raw <- query$mind_cluster_predicted.id
query$MIND_class_confidence <- query$mind_class_prediction.score.max
query$MIND_cluster_confidence <- query$mind_cluster_prediction.score.max
query$MIND_class <- ifelse(query$MIND_class_confidence >= threshold,
                           query$MIND_class_raw, "Unassigned")
query$MIND_cluster <- ifelse(query$MIND_cluster_confidence >= threshold,
                             query$MIND_cluster_raw, "Unassigned")
query@misc$MIND_reference <- list(
  paper = "Bandler et al. 2022 cells with later MIND interactive-atlas labels",
  sample = "CA301 / GSM5684876 / WT MGE E15.5",
  mapped_reference_cells = ncol(reference),
  broad_labels = sort(unique(reference$MIND_class)),
  state_labels = sort(unique(reference$MIND_cluster)),
  unassigned_threshold = threshold,
  threshold_rule = "prediction.score.max < threshold",
  note = "Later MIND 12-state labels are distinct from the original paper's 21-state taxonomy."
)

output_metadata <- query@meta.data
output_metadata$cell_id <- rownames(output_metadata)
output_metadata <- output_metadata[, c("cell_id", setdiff(names(output_metadata), "cell_id")), drop = FALSE]
write_tsv_gz(output_metadata, file.path(output_dir, "mind_label_transfer_per_cell.tsv.gz"))

composition <- as.data.frame(table(
  sample_id = output_metadata$sample_id,
  existing_cluster = output_metadata$existing_cluster,
  MIND_class = output_metadata$MIND_class,
  MIND_cluster = output_metadata$MIND_cluster
), stringsAsFactors = FALSE)
composition <- composition[composition$Freq > 0, ]
write.table(composition, file.path(output_dir, "mind_composition.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
saveRDS(query, file.path(output_dir, "paper3_query_with_mind_labels.rds"), compress = FALSE)
capture.output(sessionInfo(), file = file.path(output_dir, "R_sessionInfo.txt"))
writeLines(c(
  "PASS",
  paste0("query_cells=", ncol(query)),
  paste0("reference_cells=", ncol(reference)),
  paste0("unassigned_threshold=", threshold),
  paste0("seed=", args$seed),
  "reclustered=false",
  "integrated=false",
  "existing_umap_recomputed=false"
), file.path(output_dir, "SUCCESS.txt"))
