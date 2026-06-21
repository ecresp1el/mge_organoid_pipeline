#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

parse_args <- function(args) {
  out <- list(
    bridge_dir = NULL,
    transfer_dir = NULL,
    outdir = NULL,
    max_reference_cells = 3000L,
    max_query_cells = 3000L,
    nfeatures = 3000L,
    dims = 20L,
    seed = 0L,
    subtype_label_column = "candidate_jia_group",
    query_class_col = "div90_broad_class",
    exclude_label = "Excluded / not assigned to Jia-style 9 groups",
    k_anchor = 5L,
    k_filter = "100",
    k_score = 30L,
    max_features = 200L,
    nn_method = "annoy",
    n_trees = 50L,
    seurat_verbose = "true"
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unknown argument: ", key, call. = FALSE)
    name <- gsub("-", "_", substring(key, 3L))
    if (!(name %in% names(out))) stop("Unknown argument: ", key, call. = FALSE)
    if (i == length(args)) stop("Missing value for argument: ", key, call. = FALSE)
    out[[name]] <- args[[i + 1L]]
    i <- i + 2L
  }
  int_names <- c("max_reference_cells", "max_query_cells", "nfeatures", "dims", "seed", "k_anchor", "k_score", "max_features", "n_trees")
  for (name in int_names) out[[name]] <- as.integer(out[[name]])
  if (toupper(as.character(out$k_filter)) == "NA") out$k_filter <- NA_integer_ else out$k_filter <- as.integer(out$k_filter)
  out$seurat_verbose <- tolower(as.character(out$seurat_verbose)) %in% c("true", "t", "1", "yes", "y")
  out
}

read_tsv <- function(path) {
  utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

rss_mb <- function() {
  status_path <- "/proc/self/status"
  if (!file.exists(status_path)) return(NA_real_)
  line <- grep("^VmRSS:", readLines(status_path, warn = FALSE), value = TRUE)
  if (length(line) == 0L) return(NA_real_)
  as.numeric(gsub("[^0-9]", "", line[[1]])) / 1024
}

progress_path <- NULL
validation_path <- NULL
failure_path <- NULL

mark_progress <- function(step, status, detail = "") {
  row <- data.frame(
    timestamp = timestamp(),
    step = step,
    status = status,
    detail = as.character(detail),
    rss_mb = round(rss_mb(), 3),
    stringsAsFactors = FALSE
  )
  if (!is.null(progress_path)) {
    con <- file(progress_path, open = if (file.exists(progress_path)) "at" else "wt")
    on.exit(close(con), add = TRUE)
    utils::write.table(row, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = !file.exists(progress_path))
  }
  message("[R ", row$timestamp, "] ", step, " [", status, "]: ", detail, " rss_mb=", row$rss_mb)
  flush.console()
  invisible(NULL)
}

audit_check <- function(check, status, detail = "") {
  row <- data.frame(
    timestamp = timestamp(),
    check = check,
    status = status,
    detail = as.character(detail),
    stringsAsFactors = FALSE
  )
  con <- file(validation_path, open = if (file.exists(validation_path)) "at" else "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(row, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = !file.exists(validation_path))
  if (!identical(status, "PASS")) {
    writeLines(paste(timestamp(), check, status, detail, sep = "\t"), con = failure_path, sep = "\n", useBytes = TRUE)
    stop("Object integrity check failed: ", check, " :: ", detail, call. = FALSE)
  }
  invisible(TRUE)
}

timed_step <- function(step, expr, detail = "") {
  mark_progress(step, "start", detail)
  elapsed_start <- proc.time()[["elapsed"]]
  value <- force(expr)
  elapsed <- proc.time()[["elapsed"]] - elapsed_start
  mark_progress(step, "end", paste0("elapsed_sec=", round(elapsed, 3)))
  value
}

layer_data <- function(object, assay, layer) {
  tryCatch(
    LayerData(object = object[[assay]], layer = layer),
    error = function(e) NULL
  )
}

nonzero_count <- function(mat) {
  if (is.null(mat)) return(NA_real_)
  if (inherits(mat, "sparseMatrix")) return(length(mat@x))
  sum(mat != 0, na.rm = TRUE)
}

finite_matrix <- function(mat) {
  if (inherits(mat, "sparseMatrix")) return(all(is.finite(mat@x)))
  all(is.finite(mat))
}

validate_one_object <- function(object, object_name, assay, features, dims_use, require_pca) {
  audit_check(paste0(object_name, "_assay_present"), ifelse(assay %in% Assays(object), "PASS", "FAIL"), paste(Assays(object), collapse = ";"))
  audit_check(paste0(object_name, "_default_assay"), ifelse(DefaultAssay(object) == assay, "PASS", "FAIL"), DefaultAssay(object))
  audit_check(paste0(object_name, "_rna_assay_class"), "PASS", paste(class(object[[assay]]), collapse = ";"))
  audit_check(paste0(object_name, "_layers"), "PASS", paste(Layers(object[[assay]]), collapse = ";"))

  counts <- layer_data(object, assay, "counts")
  data <- layer_data(object, assay, "data")
  audit_check(paste0(object_name, "_counts_layer_exists"), ifelse(!is.null(counts), "PASS", "FAIL"), "RNA counts layer")
  audit_check(paste0(object_name, "_data_layer_exists"), ifelse(!is.null(data), "PASS", "FAIL"), "RNA data layer")
  audit_check(paste0(object_name, "_counts_dims"), "PASS", paste(dim(counts), collapse = "x"))
  audit_check(paste0(object_name, "_data_dims"), "PASS", paste(dim(data), collapse = "x"))
  audit_check(paste0(object_name, "_data_nonempty_nonzero"), ifelse(prod(dim(data)) > 0 && nonzero_count(data) > 0, "PASS", "FAIL"), paste0("nonzero=", nonzero_count(data)))

  audit_check(paste0(object_name, "_feature_rownames_unique"), ifelse(!anyDuplicated(rownames(object)), "PASS", "FAIL"), paste0("n_features=", nrow(object)))
  audit_check(paste0(object_name, "_cell_colnames_unique"), ifelse(!anyDuplicated(colnames(object)), "PASS", "FAIL"), paste0("n_cells=", ncol(object)))
  audit_check(paste0(object_name, "_metadata_rownames_match_cells"), ifelse(identical(rownames(object[[]]), colnames(object)), "PASS", "FAIL"), "meta.data rownames vs colnames(object)")
  audit_check(paste0(object_name, "_no_duplicated_genes"), ifelse(!anyDuplicated(rownames(data)), "PASS", "FAIL"), "RNA data rownames")
  audit_check(paste0(object_name, "_no_duplicated_cells"), ifelse(!anyDuplicated(colnames(data)), "PASS", "FAIL"), "RNA data colnames")

  missing_features <- setdiff(features, rownames(object))
  audit_check(paste0(object_name, "_selected_features_present"), ifelse(length(missing_features) == 0L, "PASS", "FAIL"), paste0("missing=", length(missing_features)))
  selected_data <- data[features, , drop = FALSE]
  audit_check(paste0(object_name, "_selected_data_finite"), ifelse(finite_matrix(selected_data), "PASS", "FAIL"), "selected RNA data")
  audit_check(paste0(object_name, "_selected_data_nonzero"), ifelse(nonzero_count(selected_data) > 0, "PASS", "FAIL"), paste0("nonzero=", nonzero_count(selected_data)))

  if (require_pca) {
    audit_check(paste0(object_name, "_pca_exists"), ifelse("pca" %in% Reductions(object), "PASS", "FAIL"), paste(Reductions(object), collapse = ";"))
    emb <- Embeddings(object, "pca")
    audit_check(paste0(object_name, "_pca_rownames_match_cells"), ifelse(identical(rownames(emb), colnames(object)), "PASS", "FAIL"), "PCA embeddings rownames vs cells")
    audit_check(paste0(object_name, "_requested_dims_available"), ifelse(max(dims_use) <= ncol(emb), "PASS", "FAIL"), paste0("requested=", max(dims_use), "; available=", ncol(emb)))
    audit_check(paste0(object_name, "_pca_no_na_nan_inf"), ifelse(all(is.finite(emb[, dims_use, drop = FALSE])), "PASS", "FAIL"), "PCA embeddings requested dims")
  }
  invisible(TRUE)
}

validate_transfer_object <- function(reference, query, assay = "RNA", features, dims_use) {
  audit_check("seurat_version", "PASS", as.character(utils::packageVersion("Seurat")))
  audit_check("seuratobject_version", "PASS", as.character(utils::packageVersion("SeuratObject")))
  audit_check("r_version", "PASS", paste(R.version$major, R.version$minor, sep = "."))
  validate_one_object(reference, "reference", assay, features, dims_use, require_pca = TRUE)
  validate_one_object(query, "query", assay, features, dims_use, require_pca = FALSE)
  shared <- intersect(rownames(reference), rownames(query))
  audit_check("shared_features_length", ifelse(length(shared) > 0L, "PASS", "FAIL"), paste0("shared=", length(shared)))
  audit_check("all_anchor_features_in_both_objects", ifelse(all(features %in% rownames(reference)) && all(features %in% rownames(query)), "PASS", "FAIL"), paste0("features=", length(features)))
  invisible(TRUE)
}

load_bridge_object <- function(prefix, bridge_dir) {
  mat <- Matrix::readMM(file.path(bridge_dir, paste0(prefix, "_counts.mtx")))
  genes <- read_tsv(file.path(bridge_dir, paste0(prefix, "_genes.tsv")))$gene
  barcodes <- read_tsv(file.path(bridge_dir, paste0(prefix, "_barcodes.tsv")))$cell_id
  meta <- read_tsv(file.path(bridge_dir, paste0(prefix, "_metadata.tsv.gz")))
  if (nrow(mat) != length(genes)) stop(prefix, ": matrix rows do not match genes", call. = FALSE)
  if (ncol(mat) != length(barcodes)) stop(prefix, ": matrix columns do not match barcodes", call. = FALSE)
  rownames(mat) <- make.unique(as.character(genes))
  colnames(mat) <- as.character(barcodes)
  rownames(meta) <- as.character(meta$seurat_cell_id)
  CreateSeuratObject(counts = mat, meta.data = meta, min.cells = 0, min.features = 0)
}

read_selected_features <- function(transfer_dir, nfeatures) {
  path <- file.path(transfer_dir, "fast_knn", "selected_transfer_features.tsv")
  if (!file.exists(path)) return(NULL)
  features <- read_tsv(path)$gene
  unique(as.character(features))[seq_len(min(length(unique(features)), nfeatures))]
}

stratified_cells <- function(meta, group_col, max_cells, seed) {
  if (max_cells <= 0L || nrow(meta) <= max_cells) return(rownames(meta))
  set.seed(seed)
  groups <- split(rownames(meta), as.character(meta[[group_col]]), drop = TRUE)
  groups <- groups[order(names(groups))]
  per_group <- ceiling(max_cells / length(groups))
  picked <- unlist(lapply(groups, function(x) if (length(x) <= per_group) x else sample(x, per_group)), use.names = FALSE)
  if (length(picked) > max_cells) picked <- sample(picked, max_cells)
  sort(picked)
}

prepare_for_anchor <- function(object, features, npcs, verbose) {
  VariableFeatures(object) <- features
  object <- NormalizeData(object, normalization.method = "LogNormalize", verbose = verbose)
  object <- ScaleData(object, features = features, verbose = verbose)
  object <- RunPCA(object, features = features, npcs = npcs, verbose = verbose)
  object
}

make_stripped_from_data_layer <- function(object, assay, features, npcs, verbose) {
  data <- layer_data(object, assay, "data")
  data <- data[features, colnames(object), drop = FALSE]
  meta <- object[[]]
  stripped <- CreateSeuratObject(counts = data, assay = assay, meta.data = meta, min.cells = 0, min.features = 0)
  stripped <- SetAssayData(stripped, assay = assay, layer = "data", new.data = data)
  DefaultAssay(stripped) <- assay
  VariableFeatures(stripped) <- features
  stripped <- ScaleData(stripped, features = features, verbose = verbose)
  stripped <- RunPCA(stripped, features = features, npcs = npcs, verbose = verbose)
  stripped
}

main <- function() {
  opt <- parse_args(commandArgs(trailingOnly = TRUE))
  if (is.null(opt$bridge_dir) || is.null(opt$transfer_dir) || is.null(opt$outdir)) {
    stop("--bridge-dir, --transfer-dir, and --outdir are required", call. = FALSE)
  }
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
  progress_path <<- file.path(opt$outdir, "s0_progress.tsv")
  validation_path <<- file.path(opt$outdir, "s0_object_integrity_audit.tsv")
  failure_path <<- file.path(opt$outdir, "s0_object_integrity_failures.log")
  for (path in c(progress_path, validation_path, failure_path)) if (file.exists(path)) file.remove(path)
  utils::capture.output(sessionInfo(), file = file.path(opt$outdir, "s0_session_info.txt"))

  mark_progress("run", "start", opt$outdir)
  reference <- timed_step("load_reference", load_bridge_object("reference", opt$bridge_dir), opt$bridge_dir)
  query <- timed_step("load_query", load_bridge_object("query", opt$bridge_dir), opt$bridge_dir)

  if (!(opt$subtype_label_column %in% colnames(reference[[]]))) stop("Missing adult subtype column: ", opt$subtype_label_column, call. = FALSE)
  if (!(opt$query_class_col %in% colnames(query[[]]))) stop("Missing query class column: ", opt$query_class_col, call. = FALSE)
  if (!is.na(opt$exclude_label) && nzchar(opt$exclude_label) && toupper(opt$exclude_label) != "NONE") {
    before <- ncol(reference)
    keep <- rownames(reference[[]])[as.character(reference[[opt$subtype_label_column]][, 1]) != opt$exclude_label]
    reference <- subset(reference, cells = keep)
    mark_progress("filter_reference_label", "end", paste0("before=", before, "; after=", ncol(reference)))
  }

  ref_cells <- stratified_cells(reference[[]], opt$subtype_label_column, opt$max_reference_cells, opt$seed)
  query_cells <- stratified_cells(query[[]], opt$query_class_col, opt$max_query_cells, opt$seed + 1L)
  reference <- subset(reference, cells = ref_cells)
  query <- subset(query, cells = query_cells)
  mark_progress("stratified_downsample", "end", paste0("reference=", ncol(reference), "; query=", ncol(query)))

  shared <- intersect(rownames(reference), rownames(query))
  selected <- read_selected_features(opt$transfer_dir, opt$nfeatures)
  if (is.null(selected)) selected <- shared
  features <- intersect(selected, shared)
  features <- features[seq_len(min(length(features), opt$nfeatures))]
  npcs <- min(50L, length(features) - 1L, ncol(reference) - 1L, ncol(query) - 1L)
  dims_use <- seq_len(min(opt$dims, npcs))
  write_tsv(data.frame(feature = features), file.path(opt$outdir, "s0_selected_features.tsv"))
  mark_progress("selected_features", "info", paste0("n=", length(features), "; npcs=", npcs, "; dims=1:", max(dims_use)))

  reference <- timed_step("prepare_original_reference", prepare_for_anchor(reference, features, npcs, opt$seurat_verbose))
  query <- timed_step("prepare_original_query", prepare_for_anchor(query, features, npcs, opt$seurat_verbose))
  mark_progress("validate_original_objects", "start", "")
  validate_transfer_object(reference, query, assay = "RNA", features = features, dims_use = dims_use)
  mark_progress("validate_original_objects", "end", "passed")

  stripped_reference <- timed_step("build_stripped_reference", make_stripped_from_data_layer(reference, "RNA", features, npcs, opt$seurat_verbose))
  stripped_query <- timed_step("build_stripped_query", make_stripped_from_data_layer(query, "RNA", features, npcs, opt$seurat_verbose))
  mark_progress("validate_stripped_objects", "start", "")
  validate_transfer_object(stripped_reference, stripped_query, assay = "RNA", features = features, dims_use = dims_use)
  mark_progress("validate_stripped_objects", "end", "passed")

  object_sizes <- data.frame(
    object = c("original_reference", "original_query", "stripped_reference", "stripped_query"),
    size_bytes = c(
      as.numeric(object.size(reference)),
      as.numeric(object.size(query)),
      as.numeric(object.size(stripped_reference)),
      as.numeric(object.size(stripped_query))
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(object_sizes, file.path(opt$outdir, "s0_object_sizes.tsv"))

  anchor_args <- list(
    reference = stripped_reference,
    query = stripped_query,
    normalization.method = "LogNormalize",
    reference.assay = "RNA",
    query.assay = "RNA",
    reduction = "pcaproject",
    reference.reduction = "pca",
    features = features,
    dims = dims_use,
    npcs = npcs,
    k.anchor = opt$k_anchor,
    k.filter = opt$k_filter,
    k.score = opt$k_score,
    max.features = opt$max_features,
    nn.method = opt$nn_method,
    n.trees = opt$n_trees,
    verbose = opt$seurat_verbose
  )
  mark_progress("find_transfer_anchors_stripped_s0", "start", paste0("features=", length(features), "; cells=", ncol(stripped_reference), "x", ncol(stripped_query), "; reduction=pcaproject; k.filter=", ifelse(is.na(opt$k_filter), "NA", opt$k_filter)))
  elapsed_start <- proc.time()[["elapsed"]]
  anchors <- do.call(FindTransferAnchors, anchor_args)
  elapsed <- proc.time()[["elapsed"]] - elapsed_start
  n_anchors <- nrow(as.data.frame(anchors@anchors))
  mark_progress("find_transfer_anchors_stripped_s0", "end", paste0("elapsed_sec=", round(elapsed, 3), "; anchors=", n_anchors))

  summary <- data.frame(
    completed = TRUE,
    n_anchors = n_anchors,
    elapsed_sec = round(elapsed, 3),
    reference_cells = ncol(stripped_reference),
    query_cells = ncol(stripped_query),
    n_features = length(features),
    dims = paste0("1:", max(dims_use)),
    reduction = "pcaproject",
    k_filter = ifelse(is.na(opt$k_filter), "NA", as.character(opt$k_filter)),
    nn_method = opt$nn_method,
    stringsAsFactors = FALSE
  )
  write_tsv(summary, file.path(opt$outdir, "s0_stripped_pcaproject_summary.tsv"))
  mark_progress("run", "end", "completed")
}

main()
