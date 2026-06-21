#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
})

options(future.globals.maxSize = 8 * 1024^3)
if (requireNamespace("future", quietly = TRUE)) {
  future::plan("sequential")
}

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
    max_reference_cells = 3000L,
    max_query_cells = 3000L,
    seed = 0L,
    k_anchor = 5L,
    k_filter = "100",
    k_score = 30L,
    max_features = 200L,
    k_weight = 50L,
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
  int_names <- c("nfeatures", "npcs", "dims", "max_reference_cells", "max_query_cells", "seed", "k_anchor", "k_score", "max_features", "k_weight", "n_trees")
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
  tryCatch(LayerData(object = object[[assay]], layer = layer), error = function(e) NULL)
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

choose_k_weight <- function(anchors, requested) {
  n_anchors <- nrow(as.data.frame(anchors@anchors))
  if (n_anchors <= requested) max(1L, n_anchors - 1L) else requested
}

palette_for <- function(labels, palette = "Dark 3") {
  labels <- sort(unique(as.character(labels)))
  cols <- grDevices::hcl.colors(max(3, length(labels)), palette = palette)
  stats::setNames(cols[seq_along(labels)], labels)
}

save_plot <- function(p, path, width = 8.5, height = 7) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggsave(path, p, width = width, height = height, dpi = 300, bg = "white")
}

plot_pc_scatter <- function(df, color_col, outfile, title, palette = NULL, alpha = 0.8, size = 0.55) {
  p <- ggplot(df, aes(PC_1, PC_2, color = .data[[color_col]])) +
    geom_point(alpha = alpha, size = size, stroke = 0) +
    coord_equal() +
    theme_classic(base_size = 11) +
    labs(title = title, x = "Reference PCA 1", y = "Reference PCA 2", color = NULL)
  if (!is.null(palette)) p <- p + scale_color_manual(values = palette, na.value = "#999999")
  save_plot(p, outfile)
}

plot_overlay <- function(ref_df, query_df, outfile, mode, div90_pal) {
  if (mode == "adult_gray_div90_classes") {
    p <- ggplot() +
      geom_point(data = ref_df, aes(PC_1, PC_2), color = "#c9c9c9", alpha = 0.28, size = 0.38) +
      geom_point(data = query_df, aes(PC_1, PC_2, color = div90_class), alpha = 0.88, size = 0.8) +
      scale_color_manual(values = div90_pal, name = "DIV90 class")
  } else if (mode == "predicted_broad") {
    pal <- c(Cortical = "#2f6db3", Subpallial = "#d4782a", Unassigned = "#999999")
    p <- ggplot() +
      geom_point(data = ref_df, aes(PC_1, PC_2), color = "#d6d6d6", alpha = 0.18, size = 0.35) +
      geom_point(data = query_df, aes(PC_1, PC_2, color = predicted_broad), alpha = 0.88, size = 0.8) +
      scale_color_manual(values = pal, name = "Predicted broad", na.value = "#999999")
  } else if (mode == "predicted_subtype") {
    pal <- palette_for(query_df$predicted_subtype, "Dark 3")
    p <- ggplot() +
      geom_point(data = ref_df, aes(PC_1, PC_2), color = "#d6d6d6", alpha = 0.18, size = 0.35) +
      geom_point(data = query_df, aes(PC_1, PC_2, color = predicted_subtype), alpha = 0.88, size = 0.8) +
      scale_color_manual(values = pal, name = "Predicted subtype", na.value = "#999999")
  } else {
    stop("Unknown overlay mode: ", mode, call. = FALSE)
  }
  p <- p + coord_equal() + theme_classic(base_size = 11) + labs(x = "Reference PCA 1", y = "Reference PCA 2")
  save_plot(p, outfile)
}

plot_river <- function(df, source_col, target_col, outfile) {
  tab <- as.data.frame(table(df[[source_col]], df[[target_col]]), stringsAsFactors = FALSE)
  colnames(tab) <- c("source", "target", "n")
  tab <- tab[tab$n > 0, , drop = FALSE]
  if (nrow(tab) == 0L) return(invisible(NULL))
  total <- sum(tab$n)
  source_tot <- aggregate(n ~ source, tab, sum)
  target_tot <- aggregate(n ~ target, tab, sum)
  source_tot <- source_tot[order(source_tot$source), , drop = FALSE]
  target_tot <- target_tot[order(target_tot$target), , drop = FALSE]
  source_tot$ymin <- c(0, head(cumsum(source_tot$n), -1)) / total
  source_tot$ymax <- cumsum(source_tot$n) / total
  target_tot$ymin <- c(0, head(cumsum(target_tot$n), -1)) / total
  target_tot$ymax <- cumsum(target_tot$n) / total
  source_cursor <- stats::setNames(source_tot$ymin, source_tot$source)
  target_cursor <- stats::setNames(target_tot$ymin, target_tot$target)
  xs <- seq(0.1, 0.9, length.out = 24)
  smooth <- ((xs - min(xs)) / diff(range(xs)))^2 * (3 - 2 * ((xs - min(xs)) / diff(range(xs))))
  polys <- list()
  for (i in seq_len(nrow(tab))) {
    row <- tab[i, ]
    h <- row$n / total
    y0l <- source_cursor[[row$source]]
    y1l <- y0l + h
    y0r <- target_cursor[[row$target]]
    y1r <- y0r + h
    source_cursor[[row$source]] <- y1l
    target_cursor[[row$target]] <- y1r
    lower <- y0l + (y0r - y0l) * smooth
    upper <- y1l + (y1r - y1l) * smooth
    polys[[i]] <- data.frame(x = c(xs, rev(xs)), y = c(lower, rev(upper)), source = row$source, target = row$target, group = i)
  }
  poly_df <- do.call(rbind, polys)
  p <- ggplot() +
    geom_polygon(data = poly_df, aes(x, y, group = group, fill = source), alpha = 0.62, color = NA) +
    geom_rect(data = source_tot, aes(xmin = 0.02, xmax = 0.08, ymin = ymin, ymax = ymax, fill = source), color = "white", linewidth = 0.2) +
    geom_rect(data = target_tot, aes(xmin = 0.92, xmax = 0.98, ymin = ymin, ymax = ymax), fill = "#4a4a4a", color = "white", linewidth = 0.2) +
    geom_text(data = source_tot, aes(x = 0.01, y = (ymin + ymax) / 2, label = source), hjust = 1, size = 2.6) +
    geom_text(data = target_tot, aes(x = 0.99, y = (ymin + ymax) / 2, label = target), hjust = 0, size = 2.6) +
    scale_fill_manual(values = palette_for(source_tot$source, "Dynamic"), guide = "none") +
    coord_cartesian(xlim = c(-0.28, 1.28), ylim = c(0, 1), clip = "off") +
    theme_void(base_size = 11) +
    theme(plot.margin = margin(10, 90, 10, 120))
  save_plot(p, outfile, width = 11, height = 7.5)
}

main <- function() {
  opt <- parse_args(commandArgs(trailingOnly = TRUE))
  if (is.null(opt$bridge_dir) || is.null(opt$transfer_dir) || is.null(opt$outdir)) {
    stop("--bridge-dir, --transfer-dir, and --outdir are required", call. = FALSE)
  }
  dirs <- file.path(opt$outdir, c("plots", "tables", "seurat"))
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  progress_path <<- file.path(opt$outdir, "pcaproject_capped_visual_progress.tsv")
  if (file.exists(progress_path)) file.remove(progress_path)
  utils::capture.output(sessionInfo(), file = file.path(opt$outdir, "session_info.txt"))

  mark_progress("run", "start", opt$outdir)
  reference <- timed_step("load_reference", load_bridge_object("reference", opt$bridge_dir), opt$bridge_dir)
  query <- timed_step("load_query", load_bridge_object("query", opt$bridge_dir), opt$bridge_dir)
  reference$adult_subtype_label <- as.character(reference[[opt$subtype_label_column]][, 1])
  reference$adult_broad_label <- adult_broad(reference[[opt$broad_label_column]][, 1])
  query$div90_class <- as.character(query[[opt$query_class_col]][, 1])

  if (!is.na(opt$exclude_label) && nzchar(opt$exclude_label) && toupper(opt$exclude_label) != "NONE") {
    reference <- subset(reference, cells = colnames(reference)[reference$adult_subtype_label != opt$exclude_label])
  }
  reference <- subset(reference, cells = stratified_cells(reference[[]], "adult_subtype_label", opt$max_reference_cells, opt$seed))
  query <- subset(query, cells = stratified_cells(query[[]], "div90_class", opt$max_query_cells, opt$seed + 1L))
  mark_progress("cell_scope", "end", paste0("reference=", ncol(reference), "; query=", ncol(query)))

  shared <- intersect(rownames(reference), rownames(query))
  selected <- read_selected_features(opt$transfer_dir, opt$nfeatures)
  if (is.null(selected)) selected <- shared
  features <- intersect(selected, shared)
  features <- features[seq_len(min(length(features), opt$nfeatures))]
  npcs <- min(opt$npcs, length(features) - 1L, ncol(reference) - 1L, ncol(query) - 1L)
  dims_use <- seq_len(min(opt$dims, npcs))
  write_tsv(data.frame(feature = features), file.path(opt$outdir, "seurat", "selected_transfer_features.tsv"))
  mark_progress("feature_scope", "end", paste0("features=", length(features), "; npcs=", npcs, "; dims=1:", max(dims_use)))

  reference <- timed_step("prepare_original_reference", prepare_original(reference, features, npcs, opt$seurat_verbose))
  query <- timed_step("prepare_original_query", prepare_original(query, features, npcs, opt$seurat_verbose))
  stripped_reference <- timed_step("build_stripped_reference", make_stripped_from_data_layer(reference, "RNA", features, npcs, opt$seurat_verbose))
  stripped_query <- timed_step("build_stripped_query", make_stripped_from_data_layer(query, "RNA", features, npcs, opt$seurat_verbose))

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
  mark_progress("find_transfer_anchors", "start", paste0("pcaproject; k.filter=", ifelse(is.na(opt$k_filter), "NA", opt$k_filter)))
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

  query_projected <- timed_step("project_query_pca_for_plots", ProjectDimReduc(
    query = stripped_query,
    reference = stripped_reference,
    mode = "pcaproject",
    reference.reduction = "pca",
    query.assay = "RNA",
    reference.assay = "RNA",
    features = features,
    do.scale = TRUE,
    reduction.name = "ref.pca",
    reduction.key = "refPC_",
    verbose = opt$seurat_verbose
  ))

  ref_pca <- as.data.frame(Embeddings(stripped_reference, "pca")[, dims_use, drop = FALSE])
  colnames(ref_pca) <- paste0("PC_", seq_len(ncol(ref_pca)))
  ref_pca$cell_id <- rownames(ref_pca)
  ref_pca$adult_subtype_label <- stripped_reference$adult_subtype_label
  ref_pca$adult_broad_label <- stripped_reference$adult_broad_label

  q_red <- if ("ref.pca" %in% Reductions(query_projected)) "ref.pca" else "pca"
  query_pca <- as.data.frame(Embeddings(query_projected, q_red)[, dims_use, drop = FALSE])
  colnames(query_pca) <- paste0("PC_", seq_len(ncol(query_pca)))
  query_pca$cell_id <- rownames(query_pca)
  query_pca$div90_class <- stripped_query$div90_class
  query_pca$predicted_broad <- broad_pred$predicted.id
  query_pca$prediction_score_broad <- broad_pred$prediction.score.max
  query_pca$predicted_subtype <- subtype_pred$predicted.id
  query_pca$prediction_score_subtype <- subtype_pred$prediction.score.max

  write_tsv(ref_pca, file.path(opt$outdir, "tables", "reference_pca_coordinates.tsv"))
  write_tsv(query_pca, file.path(opt$outdir, "tables", "query_projected_pca_coordinates_with_predictions.tsv"))
  predictions <- query_pca[, c("cell_id", "div90_class", "predicted_broad", "prediction_score_broad", "predicted_subtype", "prediction_score_subtype")]
  write_tsv(predictions, file.path(opt$outdir, "seurat_pcaproject_capped_per_cell_predictions.tsv"))

  class_counts <- as.data.frame(table(predictions$div90_class, predictions$predicted_broad, predictions$predicted_subtype), stringsAsFactors = FALSE)
  colnames(class_counts) <- c("div90_class", "predicted_broad", "predicted_subtype", "n_cells")
  class_counts <- class_counts[class_counts$n_cells > 0L, , drop = FALSE]
  score_summary <- aggregate(
    cbind(prediction_score_broad, prediction_score_subtype) ~ div90_class + predicted_broad + predicted_subtype,
    predictions,
    function(x) mean(x, na.rm = TRUE)
  )
  class_summary <- merge(class_counts, score_summary, by = c("div90_class", "predicted_broad", "predicted_subtype"), all.x = TRUE)
  write_tsv(class_summary, file.path(opt$outdir, "seurat_pcaproject_capped_prediction_scores_by_class.tsv"))

  anchor_summary <- data.frame(
    method = "capped_stripped_rna_data_pcaproject",
    reference_cells = ncol(stripped_reference),
    query_cells = ncol(stripped_query),
    n_features = length(features),
    npcs = npcs,
    dims = paste0("1:", max(dims_use)),
    k_filter = ifelse(is.na(opt$k_filter), "NA", as.character(opt$k_filter)),
    k_weight = k_weight_used,
    nn_method = opt$nn_method,
    n_anchors = anchor_count,
    find_transfer_anchors_elapsed_sec = round(anchor_elapsed, 3),
    stringsAsFactors = FALSE
  )
  write_tsv(anchor_summary, file.path(opt$outdir, "seurat_pcaproject_capped_anchor_summary.tsv"))

  div90_pal <- palette_for(query_pca$div90_class, "Dynamic")
  plot_pc_scatter(ref_pca, "adult_subtype_label", file.path(opt$outdir, "plots", "PCA_reference_adult_subtypes.png"), "Adult reference PCA: subtype", palette_for(ref_pca$adult_subtype_label, "Dark 3"), alpha = 0.72, size = 0.55)
  plot_pc_scatter(ref_pca, "adult_broad_label", file.path(opt$outdir, "plots", "PCA_reference_adult_broad.png"), "Adult reference PCA: broad", c(Cortical = "#2f6db3", Subpallial = "#d4782a", Unassigned = "#999999"), alpha = 0.72, size = 0.55)
  plot_overlay(ref_pca, query_pca, file.path(opt$outdir, "plots", "PCA_adult_gray_DIV90_classes.png"), "adult_gray_div90_classes", div90_pal)
  plot_overlay(ref_pca, query_pca, file.path(opt$outdir, "plots", "PCA_DIV90_predicted_broad.png"), "predicted_broad", div90_pal)
  plot_overlay(ref_pca, query_pca, file.path(opt$outdir, "plots", "PCA_DIV90_predicted_subtype.png"), "predicted_subtype", div90_pal)

  p_score <- ggplot(predictions, aes(div90_class, prediction_score_subtype, fill = div90_class)) +
    geom_boxplot(outlier.size = 0.25, linewidth = 0.25) +
    scale_fill_manual(values = div90_pal, guide = "none") +
    coord_flip() +
    theme_classic(base_size = 10) +
    labs(x = "DIV90 class", y = "Prediction score max", title = "Subtype prediction scores by DIV90 class")
  save_plot(p_score, file.path(opt$outdir, "plots", "prediction_score_by_DIV90_class.png"), width = 8, height = 7)

  plot_river(predictions, "div90_class", "predicted_broad", file.path(opt$outdir, "plots", "river_DIV90_class_to_predicted_broad.png"))
  plot_river(predictions, "div90_class", "predicted_subtype", file.path(opt$outdir, "plots", "river_DIV90_class_to_predicted_subtype.png"))

  report <- c(
    "# Capped stripped RNA-only pcaproject visual diagnostic",
    "",
    "This is a diagnostic visualization of the PCA/pcaproject approach that worked in the S0 test.",
    "It is capped and stratified, so it is not a final full-reference result.",
    "",
    "## Inputs",
    "",
    paste0("- Bridge: `", opt$bridge_dir, "`"),
    paste0("- Transfer features: `", file.path(opt$transfer_dir, "fast_knn", "selected_transfer_features.tsv"), "`"),
    "",
    "## Parameters",
    "",
    paste0("- Reference cells: ", ncol(stripped_reference)),
    paste0("- Query cells: ", ncol(stripped_query)),
    paste0("- Features: ", length(features)),
    paste0("- Dims: 1:", max(dims_use)),
    paste0("- k.filter: ", ifelse(is.na(opt$k_filter), "NA", opt$k_filter)),
    paste0("- Anchors retained: ", anchor_count),
    paste0("- FindTransferAnchors elapsed sec: ", round(anchor_elapsed, 3)),
    "",
    "## Primary plots",
    "",
    "- `plots/PCA_reference_adult_subtypes.png`",
    "- `plots/PCA_reference_adult_broad.png`",
    "- `plots/PCA_adult_gray_DIV90_classes.png`",
    "- `plots/PCA_DIV90_predicted_broad.png`",
    "- `plots/PCA_DIV90_predicted_subtype.png`",
    "- `plots/river_DIV90_class_to_predicted_broad.png`",
    "- `plots/river_DIV90_class_to_predicted_subtype.png`",
    "- `plots/prediction_score_by_DIV90_class.png`"
  )
  writeLines(report, file.path(opt$outdir, "README_capped_visual_diagnostic.md"))
  mark_progress("run", "end", "completed")
}

main()
