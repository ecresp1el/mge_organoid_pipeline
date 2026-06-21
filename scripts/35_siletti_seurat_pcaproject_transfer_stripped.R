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
    subtype_label_column = "candidate_jia_group",
    broad_label_column = "jia_side",
    query_class_col = "div90_broad_class",
    exclude_label = "Excluded / not assigned to Jia-style 9 groups",
    nfeatures = 3000L,
    npcs = 50L,
    dims = 20L,
    seed = 0L,
    max_reference_cells = 0L,
    max_query_cells = 0L,
    k_anchor = 5L,
    k_filter = "100",
    k_score = 30L,
    max_features = 200L,
    k_weight = 50L,
    nn_method = "annoy",
    n_trees = 50L,
    seurat_verbose = "true",
    save_rds = "false"
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
  int_names <- c("nfeatures", "npcs", "dims", "seed", "max_reference_cells", "max_query_cells", "k_anchor", "k_score", "max_features", "k_weight", "n_trees")
  for (name in int_names) out[[name]] <- as.integer(out[[name]])
  if (toupper(as.character(out$k_filter)) == "NA") out$k_filter <- NA_integer_ else out$k_filter <- as.integer(out$k_filter)
  for (name in c("seurat_verbose", "save_rds")) {
    out[[name]] <- tolower(as.character(out[[name]])) %in% c("true", "t", "1", "yes", "y")
  }
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

adult_broad <- function(x) {
  x <- as.character(x)
  out <- ifelse(grepl("^Cortical", x), "Cortical", ifelse(grepl("^Subpallial", x), "Subpallial", "Unassigned"))
  out[is.na(x) | !nzchar(x)] <- "Unassigned"
  out
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

nonzero_count <- function(mat) {
  if (inherits(mat, "sparseMatrix")) return(length(mat@x))
  sum(mat != 0, na.rm = TRUE)
}

validate_transfer_object <- function(reference, query, assay = "RNA", features, dims_use) {
  for (object_name in c("reference", "query")) {
    object <- if (object_name == "reference") reference else query
    if (!(assay %in% Assays(object))) stop(object_name, ": missing assay ", assay, call. = FALSE)
    if (DefaultAssay(object) != assay) stop(object_name, ": DefaultAssay is not ", assay, call. = FALSE)
    data <- layer_data(object, assay, "data")
    counts <- layer_data(object, assay, "counts")
    if (is.null(counts)) stop(object_name, ": missing RNA counts layer", call. = FALSE)
    if (is.null(data)) stop(object_name, ": missing RNA data layer", call. = FALSE)
    if (prod(dim(data)) == 0L || nonzero_count(data) == 0L) stop(object_name, ": RNA data layer empty/nonzero=0", call. = FALSE)
    if (anyDuplicated(rownames(object))) stop(object_name, ": duplicated genes/features", call. = FALSE)
    if (anyDuplicated(colnames(object))) stop(object_name, ": duplicated cells", call. = FALSE)
    if (!identical(rownames(object[[]]), colnames(object))) stop(object_name, ": metadata rownames do not match cell names", call. = FALSE)
    missing_features <- setdiff(features, rownames(object))
    if (length(missing_features) > 0L) stop(object_name, ": missing selected features: ", length(missing_features), call. = FALSE)
    selected <- data[features, , drop = FALSE]
    if (inherits(selected, "sparseMatrix")) {
      if (!all(is.finite(selected@x))) stop(object_name, ": non-finite selected RNA data", call. = FALSE)
    } else if (!all(is.finite(selected))) {
      stop(object_name, ": non-finite selected RNA data", call. = FALSE)
    }
  }
  if (!("pca" %in% Reductions(reference))) stop("reference: missing PCA reduction", call. = FALSE)
  emb <- Embeddings(reference, "pca")
  if (!identical(rownames(emb), colnames(reference))) stop("reference: PCA rownames do not match cells", call. = FALSE)
  if (max(dims_use) > ncol(emb)) stop("requested dims exceed available reference PCs", call. = FALSE)
  if (!all(is.finite(emb[, dims_use, drop = FALSE]))) stop("reference: PCA embeddings have NA/NaN/Inf", call. = FALSE)
  if (length(intersect(rownames(reference), rownames(query))) < length(features)) stop("shared feature set smaller than selected features", call. = FALSE)
  invisible(TRUE)
}

prepare_original <- function(object, features, npcs, verbose) {
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

prediction_score_cols <- function(df) {
  setdiff(grep("^prediction\\.score\\.", colnames(df), value = TRUE), "prediction.score.max")
}

choose_k_weight <- function(anchors, requested) {
  n_anchors <- nrow(as.data.frame(anchors@anchors))
  if (n_anchors <= requested) max(1L, n_anchors - 1L) else requested
}

main <- function() {
  opt <- parse_args(commandArgs(trailingOnly = TRUE))
  if (is.null(opt$bridge_dir) || is.null(opt$transfer_dir) || is.null(opt$outdir)) {
    stop("--bridge-dir, --transfer-dir, and --outdir are required", call. = FALSE)
  }
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(opt$outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(opt$outdir, "seurat"), recursive = TRUE, showWarnings = FALSE)
  progress_path <<- file.path(opt$outdir, "seurat_pcaproject_transfer_progress.tsv")
  if (file.exists(progress_path)) file.remove(progress_path)
  utils::capture.output(sessionInfo(), file = file.path(opt$outdir, "seurat_pcaproject_session_info.txt"))

  mark_progress("run", "start", opt$outdir)
  reference <- timed_step("load_reference", load_bridge_object("reference", opt$bridge_dir), opt$bridge_dir)
  query <- timed_step("load_query", load_bridge_object("query", opt$bridge_dir), opt$bridge_dir)
  mark_progress("object_sizes_initial", "info", paste0("reference=", ncol(reference), "; query=", ncol(query), "; genes=", nrow(reference)))

  if (!(opt$subtype_label_column %in% colnames(reference[[]]))) stop("Missing subtype label column: ", opt$subtype_label_column, call. = FALSE)
  if (!(opt$broad_label_column %in% colnames(reference[[]]))) stop("Missing broad label column: ", opt$broad_label_column, call. = FALSE)
  if (!(opt$query_class_col %in% colnames(query[[]]))) stop("Missing query class column: ", opt$query_class_col, call. = FALSE)

  reference$adult_subtype_label <- as.character(reference[[opt$subtype_label_column]][, 1])
  reference$adult_broad_label <- adult_broad(reference[[opt$broad_label_column]][, 1])
  query$div90_class <- as.character(query[[opt$query_class_col]][, 1])

  if (!is.na(opt$exclude_label) && nzchar(opt$exclude_label) && toupper(opt$exclude_label) != "NONE") {
    before <- ncol(reference)
    keep <- colnames(reference)[reference$adult_subtype_label != opt$exclude_label]
    reference <- subset(reference, cells = keep)
    mark_progress("filter_reference_label", "end", paste0("before=", before, "; after=", ncol(reference), "; excluded=", opt$exclude_label))
  }

  if (opt$max_reference_cells > 0L) {
    reference <- subset(reference, cells = stratified_cells(reference[[]], "adult_subtype_label", opt$max_reference_cells, opt$seed))
  }
  if (opt$max_query_cells > 0L) {
    query <- subset(query, cells = stratified_cells(query[[]], "div90_class", opt$max_query_cells, opt$seed + 1L))
  }
  mark_progress("cell_scope", "info", paste0("reference=", ncol(reference), "; query=", ncol(query)))

  shared <- intersect(rownames(reference), rownames(query))
  selected <- read_selected_features(opt$transfer_dir, opt$nfeatures)
  if (is.null(selected)) selected <- shared
  features <- intersect(selected, shared)
  features <- features[seq_len(min(length(features), opt$nfeatures))]
  if (length(features) < 50L) stop("Too few selected shared features: ", length(features), call. = FALSE)
  write_tsv(data.frame(feature = features), file.path(opt$outdir, "seurat", "selected_transfer_features.tsv"))

  npcs <- min(opt$npcs, length(features) - 1L, ncol(reference) - 1L, ncol(query) - 1L)
  dims_use <- seq_len(min(opt$dims, npcs))
  mark_progress("feature_scope", "info", paste0("features=", length(features), "; npcs=", npcs, "; dims=1:", max(dims_use)))

  reference <- timed_step("prepare_original_reference", prepare_original(reference, features, npcs, opt$seurat_verbose))
  query <- timed_step("prepare_original_query", prepare_original(query, features, npcs, opt$seurat_verbose))
  timed_step("validate_original_objects", validate_transfer_object(reference, query, assay = "RNA", features = features, dims_use = dims_use))

  stripped_reference <- timed_step("build_stripped_reference", make_stripped_from_data_layer(reference, "RNA", features, npcs, opt$seurat_verbose))
  stripped_query <- timed_step("build_stripped_query", make_stripped_from_data_layer(query, "RNA", features, npcs, opt$seurat_verbose))
  timed_step("validate_stripped_objects", validate_transfer_object(stripped_reference, stripped_query, assay = "RNA", features = features, dims_use = dims_use))

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
  write_tsv(object_sizes, file.path(opt$outdir, "tables", "seurat_pcaproject_object_sizes.tsv"))

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
  mark_progress("find_transfer_anchors", "start", paste0("reduction=pcaproject; features=", length(features), "; dims=1:", max(dims_use), "; k.filter=", ifelse(is.na(opt$k_filter), "NA", opt$k_filter)))
  anchor_start <- proc.time()[["elapsed"]]
  anchors <- do.call(FindTransferAnchors, anchor_args)
  anchor_elapsed <- proc.time()[["elapsed"]] - anchor_start
  anchor_count <- nrow(as.data.frame(anchors@anchors))
  mark_progress("find_transfer_anchors", "end", paste0("elapsed_sec=", round(anchor_elapsed, 3), "; anchors=", anchor_count))

  k_weight_used <- choose_k_weight(anchors, opt$k_weight)
  broad_pred <- timed_step("transfer_broad", as.data.frame(TransferData(
    anchorset = anchors,
    refdata = stripped_reference$adult_broad_label,
    dims = dims_use,
    k.weight = k_weight_used,
    verbose = opt$seurat_verbose
  )))
  subtype_pred <- timed_step("transfer_subtype", as.data.frame(TransferData(
    anchorset = anchors,
    refdata = stripped_reference$adult_subtype_label,
    dims = dims_use,
    k.weight = k_weight_used,
    verbose = opt$seurat_verbose
  )))

  predictions <- data.frame(
    cell_id = colnames(stripped_query),
    div90_class = stripped_query$div90_class,
    predicted_broad = broad_pred$predicted.id,
    prediction_score_broad = broad_pred$prediction.score.max,
    predicted_subtype = subtype_pred$predicted.id,
    prediction_score_subtype = subtype_pred$prediction.score.max,
    stringsAsFactors = FALSE
  )
  write_tsv(predictions, file.path(opt$outdir, "seurat_pcaproject_per_cell_predictions.tsv"))

  class_summary <- aggregate(
    cbind(prediction_score_broad, prediction_score_subtype) ~ div90_class + predicted_broad + predicted_subtype,
    predictions,
    function(x) mean(x, na.rm = TRUE)
  )
  class_counts <- as.data.frame(table(predictions$div90_class, predictions$predicted_broad, predictions$predicted_subtype), stringsAsFactors = FALSE)
  colnames(class_counts) <- c("div90_class", "predicted_broad", "predicted_subtype", "n_cells")
  class_counts <- class_counts[class_counts$n_cells > 0L, , drop = FALSE]
  class_summary <- merge(class_counts, class_summary, by = c("div90_class", "predicted_broad", "predicted_subtype"), all.x = TRUE)
  write_tsv(class_summary, file.path(opt$outdir, "seurat_pcaproject_prediction_scores_by_class.tsv"))

  anchor_summary <- data.frame(
    method = "stripped_rna_data_pcaproject",
    bridge_dir = opt$bridge_dir,
    transfer_dir = opt$transfer_dir,
    reference_cells = ncol(stripped_reference),
    query_cells = ncol(stripped_query),
    n_features = length(features),
    npcs = npcs,
    dims = paste0("1:", max(dims_use)),
    k_anchor = opt$k_anchor,
    k_filter = ifelse(is.na(opt$k_filter), "NA", as.character(opt$k_filter)),
    k_score = opt$k_score,
    max_features = opt$max_features,
    k_weight = k_weight_used,
    nn_method = opt$nn_method,
    n_anchors = anchor_count,
    find_transfer_anchors_elapsed_sec = round(anchor_elapsed, 3),
    stringsAsFactors = FALSE
  )
  write_tsv(anchor_summary, file.path(opt$outdir, "seurat_pcaproject_anchor_summary.tsv"))
  if (opt$save_rds) {
    saveRDS(anchors, file.path(opt$outdir, "seurat", "seurat_pcaproject_anchors.rds"))
  }
  mark_progress("run", "end", "completed")
}

main()
