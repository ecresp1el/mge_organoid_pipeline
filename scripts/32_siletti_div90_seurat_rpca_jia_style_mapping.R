#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
})

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_msg <- function(...) message("[R ", timestamp(), "] ", paste0(..., collapse = ""))

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
    k_weight = 50L,
    k_anchor = 5L,
    k_filter = "NA",
    k_score = 30L,
    max_features = 200L,
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
  int_names <- c("nfeatures", "npcs", "dims", "seed", "k_weight", "k_anchor", "k_score", "max_features", "n_trees")
  for (name in int_names) out[[name]] <- as.integer(out[[name]])
  for (name in c("seurat_verbose", "save_rds")) {
    out[[name]] <- tolower(as.character(out[[name]])) %in% c("true", "t", "1", "yes", "y")
  }
  if (toupper(as.character(out$k_filter)) == "NA") out$k_filter <- NA_integer_ else out$k_filter <- as.integer(out$k_filter)
  out
}

read_tsv <- function(path) {
  utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

progress_path <- NULL
mark_progress <- function(step, status, detail = "") {
  row <- data.frame(
    timestamp = timestamp(),
    step = step,
    status = status,
    detail = as.character(detail),
    stringsAsFactors = FALSE
  )
  if (!is.null(progress_path)) {
    con <- file(progress_path, open = if (file.exists(progress_path)) "at" else "wt")
    on.exit(close(con), add = TRUE)
    utils::write.table(row, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = !file.exists(progress_path))
  }
  log_msg(step, " [", status, "]: ", detail)
  flush.console()
  invisible(NULL)
}

timed_step <- function(step, expr, detail = "") {
  mark_progress(step, "start", detail)
  elapsed_start <- proc.time()[["elapsed"]]
  value <- force(expr)
  elapsed <- proc.time()[["elapsed"]] - elapsed_start
  mark_progress(step, "end", paste0("elapsed_sec=", round(elapsed, 2)))
  value
}

sanitize_token <- function(x) {
  y <- gsub("[^A-Za-z0-9]+", "_", trimws(as.character(x)))
  y <- gsub("_+", "_", y)
  y <- gsub("^_|_$", "", y)
  ifelse(nzchar(y), y, "value")
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
  CreateSeuratObject(counts = mat, meta.data = meta)
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

prediction_score_cols <- function(df) {
  setdiff(grep("^prediction\\.score\\.", colnames(df), value = TRUE), "prediction.score.max")
}

choose_transfer_k_weight <- function(anchors, requested = 50L) {
  n_anchors <- nrow(as.data.frame(anchors@anchors))
  if (n_anchors <= requested) max(1L, n_anchors - 1L) else requested
}

make_label_palette <- function(labels) {
  labels <- sort(unique(as.character(labels)))
  cols <- grDevices::hcl.colors(max(3, length(labels)), palette = "Dark 3")
  stats::setNames(cols[seq_along(labels)], labels)
}

div90_palette <- function(labels) {
  labels <- sort(unique(as.character(labels)))
  cols <- grDevices::hcl.colors(max(3, length(labels)), palette = "Dynamic")
  stats::setNames(cols[seq_along(labels)], labels)
}

plot_overlay <- function(ref_df, query_df, outfile, mode, palette = NULL) {
  p <- ggplot()
  if (mode == "adult_gray_div90_classes") {
    pal <- palette
    p <- p +
      geom_point(data = ref_df, aes(umap_1, umap_2), color = "#d0d0d0", alpha = 0.34, size = 0.32) +
      geom_point(data = query_df, aes(umap_1, umap_2, color = div90_class), alpha = 0.85, size = 0.85) +
      scale_color_manual(values = pal, name = "DIV90 class")
  } else if (mode == "adult_broad_div90_overlay") {
    side_pal <- c(Cortical = "#6c8ebf", Subpallial = "#d89555", Unassigned = "#b8b8b8")
    pal <- palette
    p <- p +
      geom_point(data = ref_df, aes(umap_1, umap_2, color = adult_broad_label), alpha = 0.24, size = 0.34) +
      geom_point(data = query_df, aes(umap_1, umap_2, fill = div90_class), shape = 21, color = NA, alpha = 0.86, size = 0.9) +
      scale_color_manual(values = side_pal, name = "Adult broad") +
      scale_fill_manual(values = pal, name = "DIV90 class")
  } else if (mode == "adult_subtypes_div90_black") {
    subtype_pal <- make_label_palette(ref_df$adult_subtype_label)
    p <- p +
      geom_point(data = ref_df, aes(umap_1, umap_2, color = adult_subtype_label), alpha = 0.42, size = 0.35) +
      geom_point(data = query_df, aes(umap_1, umap_2), color = "black", alpha = 0.50, size = 0.72) +
      scale_color_manual(values = subtype_pal, name = "Adult subtype")
  } else if (mode == "div90_predicted_broad") {
    pred_pal <- c(Cortical = "#2f6db3", Subpallial = "#d4782a", Unassigned = "#999999")
    p <- p +
      geom_point(data = ref_df, aes(umap_1, umap_2), color = "#d9d9d9", alpha = 0.16, size = 0.25) +
      geom_point(data = query_df, aes(umap_1, umap_2, color = predicted_broad), alpha = 0.86, size = 0.85) +
      scale_color_manual(values = pred_pal, name = "Predicted broad", na.value = "#999999")
  } else {
    stop("Unknown plot mode: ", mode, call. = FALSE)
  }
  p <- p +
    coord_equal() +
    theme_classic(base_size = 11) +
    labs(x = "Seurat RPCA reference UMAP 1", y = "Seurat RPCA reference UMAP 2")
  ggsave(outfile, p, width = 8.8, height = 7.2, dpi = 300, bg = "white")
}

plot_river <- function(df, source_col, target_col, outfile, fill_title) {
  tab <- as.data.frame(table(df[[source_col]], df[[target_col]]), stringsAsFactors = FALSE)
  colnames(tab) <- c("source", "target", "n")
  tab <- tab[tab$n > 0, , drop = FALSE]
  tab$source <- as.character(tab$source)
  tab$target <- as.character(tab$target)
  tab <- tab[order(tab$source, tab$target), , drop = FALSE]
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
  xs <- seq(0.08, 0.92, length.out = 32)
  t <- (xs - min(xs)) / diff(range(xs))
  smooth <- t * t * (3 - 2 * t)
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
    polys[[i]] <- data.frame(
      flow = i,
      source = row$source,
      target = row$target,
      x = c(xs, rev(xs)),
      y = c(lower, rev(upper)),
      stringsAsFactors = FALSE
    )
  }
  poly_df <- do.call(rbind, polys)
  pal <- make_label_palette(target_tot$target)
  p <- ggplot() +
    geom_polygon(data = poly_df, aes(x, y, group = flow, fill = target), alpha = 0.68, color = NA) +
    geom_rect(data = source_tot, aes(xmin = 0.00, xmax = 0.045, ymin = ymin, ymax = ymax), fill = "#4d4d4d") +
    geom_rect(data = target_tot, aes(xmin = 0.955, xmax = 1.00, ymin = ymin, ymax = ymax, fill = target), color = NA) +
    geom_text(data = source_tot, aes(x = -0.01, y = (ymin + ymax) / 2, label = source), hjust = 1, size = 2.7) +
    geom_text(data = target_tot, aes(x = 1.01, y = (ymin + ymax) / 2, label = target), hjust = 0, size = 2.7) +
    scale_fill_manual(values = pal, name = fill_title) +
    coord_cartesian(xlim = c(-0.42, 1.42), ylim = c(0, 1), clip = "off") +
    theme_void(base_size = 11) +
    theme(legend.position = "none", plot.margin = margin(8, 90, 8, 150))
  ggsave(outfile, p, width = 10.5, height = 6.8, dpi = 300, bg = "white")
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("bridge_dir", "transfer_dir", "outdir")
for (name in required) {
  if (is.null(opt[[name]]) || !nzchar(opt[[name]])) stop("--", gsub("_", "-", name), " is required", call. = FALSE)
}

set.seed(opt$seed)
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(opt$outdir, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(opt$outdir, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(opt$outdir, "seurat"), showWarnings = FALSE, recursive = TRUE)
progress_path <- file.path(opt$outdir, "seurat_rpca_run_audit.tsv")

log_msg("Seurat version: ", as.character(packageVersion("Seurat")))
log_msg("SeuratObject version: ", as.character(packageVersion("SeuratObject")))
log_msg("Bridge dir: ", opt$bridge_dir)
log_msg("Transfer dir: ", opt$transfer_dir)
log_msg("Outdir: ", opt$outdir)
mark_progress("run", "start", paste0("outdir=", opt$outdir))

reference <- timed_step("load_reference", load_bridge_object("reference", opt$bridge_dir), opt$bridge_dir)
query <- timed_step("load_query", load_bridge_object("query", opt$bridge_dir), opt$bridge_dir)
mark_progress("object_sizes", "info", paste0("reference_cells=", ncol(reference), "; query_cells=", ncol(query), "; genes=", nrow(reference)))

if (!(opt$subtype_label_column %in% colnames(reference@meta.data))) {
  stop("Reference metadata missing subtype label column: ", opt$subtype_label_column, call. = FALSE)
}
if (!(opt$broad_label_column %in% colnames(reference@meta.data))) {
  stop("Reference metadata missing broad label column: ", opt$broad_label_column, call. = FALSE)
}
if (!(opt$query_class_col %in% colnames(query@meta.data))) {
  stop("Query metadata missing DIV90 class column: ", opt$query_class_col, call. = FALSE)
}

reference$adult_subtype_label <- as.character(reference@meta.data[[opt$subtype_label_column]])
reference$adult_broad_label <- adult_broad(reference@meta.data[[opt$broad_label_column]])
query$div90_class <- as.character(query@meta.data[[opt$query_class_col]])

if (!is.null(opt$exclude_label) && nzchar(opt$exclude_label) && opt$exclude_label != "NONE") {
  before <- ncol(reference)
  keep <- !is.na(reference$adult_subtype_label) & reference$adult_subtype_label != opt$exclude_label
  reference <- subset(reference, cells = colnames(reference)[keep])
  mark_progress("filter_reference_label", "info", paste0("excluded_label=", opt$exclude_label, "; reference_cells_before=", before, "; reference_cells_after=", ncol(reference)))
}

selected_features <- read_selected_features(opt$transfer_dir, opt$nfeatures)
shared_features <- intersect(rownames(reference), rownames(query))
if (is.null(selected_features)) {
  mark_progress("selected_features", "info", "selected_transfer_features.tsv not found; using Seurat variable shared genes")
} else {
  features <- intersect(selected_features, shared_features)
  if (length(features) < 500L) stop("Only ", length(features), " selected features overlap Seurat objects", call. = FALSE)
  mark_progress("selected_features", "info", paste0("source=existing_fast_knn; n_features=", length(features)))
}

reference <- timed_step("normalize_reference", NormalizeData(reference, normalization.method = "LogNormalize", verbose = opt$seurat_verbose))
query <- timed_step("normalize_query", NormalizeData(query, normalization.method = "LogNormalize", verbose = opt$seurat_verbose))

if (is.null(selected_features)) {
  reference <- timed_step("variable_features_reference", FindVariableFeatures(reference, selection.method = "vst", nfeatures = opt$nfeatures, verbose = opt$seurat_verbose))
  query <- timed_step("variable_features_query", FindVariableFeatures(query, selection.method = "vst", nfeatures = opt$nfeatures, verbose = opt$seurat_verbose))
  features <- unique(c(VariableFeatures(reference), VariableFeatures(query)))
  features <- intersect(features, shared_features)
  features <- features[seq_len(min(length(features), opt$nfeatures))]
}
write_tsv(data.frame(feature = features), file.path(opt$outdir, "seurat", "selected_transfer_features.tsv"))

VariableFeatures(reference) <- features
VariableFeatures(query) <- features
reference <- timed_step("scale_reference", ScaleData(reference, features = features, verbose = opt$seurat_verbose))
query <- timed_step("scale_query", ScaleData(query, features = features, verbose = opt$seurat_verbose))
reference <- timed_step("pca_reference", RunPCA(reference, features = features, npcs = opt$npcs, verbose = opt$seurat_verbose))
query <- timed_step("pca_query", RunPCA(query, features = features, npcs = opt$npcs, verbose = opt$seurat_verbose))

dims_use <- seq_len(min(opt$dims, opt$npcs))
reference <- timed_step(
  "reference_umap_model",
  RunUMAP(reference, reduction = "pca", dims = dims_use, reduction.name = "umap", return.model = TRUE, verbose = opt$seurat_verbose),
  paste0("dims=1:", max(dims_use))
)

anchor_args <- list(
  reference = reference,
  query = query,
  normalization.method = "LogNormalize",
  reference.assay = "RNA",
  query.assay = "RNA",
  reduction = "rpca",
  reference.reduction = "pca",
  features = features,
  dims = dims_use,
  npcs = opt$npcs,
  k.anchor = opt$k_anchor,
  k.filter = opt$k_filter,
  k.score = opt$k_score,
  max.features = opt$max_features,
  nn.method = opt$nn_method,
  n.trees = opt$n_trees,
  verbose = opt$seurat_verbose
)
mark_progress("find_transfer_anchors_args", "info", paste0("reduction=rpca; dims=1:", max(dims_use), "; features=", length(features), "; k.filter=", opt$k_filter))
anchors <- timed_step("find_transfer_anchors", do.call(FindTransferAnchors, anchor_args))
anchor_count <- nrow(as.data.frame(anchors@anchors))
k_weight_used <- choose_transfer_k_weight(anchors, opt$k_weight)
mark_progress("anchor_summary", "info", paste0("n_anchors=", anchor_count, "; transfer_k_weight=", k_weight_used))

broad_pred <- timed_step("transfer_broad", as.data.frame(TransferData(
  anchorset = anchors,
  refdata = reference$adult_broad_label,
  dims = dims_use,
  k.weight = k_weight_used,
  verbose = opt$seurat_verbose
)))
subtype_pred <- timed_step("transfer_subtype", as.data.frame(TransferData(
  anchorset = anchors,
  refdata = reference$adult_subtype_label,
  dims = dims_use,
  k.weight = k_weight_used,
  verbose = opt$seurat_verbose
)))

mark_progress("map_query", "start", "MapQuery with refdata broad/subtype and reduction.model=umap")
query_mapped <- MapQuery(
  anchorset = anchors,
  query = query,
  reference = reference,
  refdata = list(
    predicted_broad = "adult_broad_label",
    predicted_subtype = "adult_subtype_label"
  ),
  new.reduction.name = "integrated.rpca",
  reference.reduction = "pca",
  reference.dims = dims_use,
  query.dims = dims_use,
  reduction.model = "umap",
  transferdata.args = list(k.weight = k_weight_used),
  integrateembeddings.args = list(dims.to.integrate = dims_use, k.weight = k_weight_used),
  projectumap.args = list(k.param = 30, n.trees = opt$n_trees),
  verbose = opt$seurat_verbose
)
mark_progress("map_query", "end", "complete")

query_umap_name <- if ("ref.umap" %in% names(query_mapped@reductions)) "ref.umap" else if ("umap" %in% names(query_mapped@reductions)) "umap" else tail(names(query_mapped@reductions), 1)
ref_umap <- as.data.frame(Embeddings(reference, "umap"))
query_umap <- as.data.frame(Embeddings(query_mapped, query_umap_name))
colnames(ref_umap)[1:2] <- c("umap_1", "umap_2")
colnames(query_umap)[1:2] <- c("umap_1", "umap_2")

ref_df <- cbind(
  data.frame(cell_id = colnames(reference), adult_broad_label = reference$adult_broad_label, adult_subtype_label = reference$adult_subtype_label, stringsAsFactors = FALSE),
  ref_umap[, c("umap_1", "umap_2")]
)
query_df <- cbind(
  data.frame(
    cell_id = colnames(query_mapped),
    div90_class = query_mapped$div90_class,
    predicted_broad = broad_pred$predicted.id,
    predicted_broad_score = broad_pred$prediction.score.max,
    predicted_subtype = subtype_pred$predicted.id,
    predicted_subtype_score = subtype_pred$prediction.score.max,
    stringsAsFactors = FALSE
  ),
  query_umap[, c("umap_1", "umap_2")]
)

per_cell <- cbind(
  query_mapped@meta.data[, setdiff(colnames(query_mapped@meta.data), c("predicted_broad", "predicted_subtype")), drop = FALSE],
  query_df[, c("cell_id", "div90_class", "predicted_broad", "predicted_broad_score", "predicted_subtype", "predicted_subtype_score", "umap_1", "umap_2"), drop = FALSE]
)
per_cell <- per_cell[, !duplicated(colnames(per_cell)), drop = FALSE]
write_tsv(per_cell, file.path(opt$outdir, "seurat_rpca_per_cell_predictions.tsv"))

class_summary <- aggregate(
  cbind(predicted_broad_score, predicted_subtype_score) ~ div90_class,
  query_df,
  function(x) c(mean = mean(x, na.rm = TRUE), median = stats::median(x, na.rm = TRUE))
)
class_summary <- do.call(data.frame, class_summary)
names(class_summary) <- gsub("\\.", "_", names(class_summary))
broad_counts <- as.data.frame.matrix(table(query_df$div90_class, query_df$predicted_broad))
broad_counts$div90_class <- rownames(broad_counts)
subtype_counts <- as.data.frame.matrix(table(query_df$div90_class, query_df$predicted_subtype))
subtype_counts$div90_class <- rownames(subtype_counts)
class_summary <- merge(class_summary, broad_counts, by = "div90_class", all.x = TRUE)
write_tsv(class_summary, file.path(opt$outdir, "seurat_rpca_prediction_scores_by_class.tsv"))

pal <- div90_palette(query_df$div90_class)
plot_overlay(ref_df, query_df, file.path(opt$outdir, "FINAL_seurat_rpca_adult_gray_div90_classes.png"), "adult_gray_div90_classes", pal)
plot_overlay(ref_df, query_df, file.path(opt$outdir, "FINAL_seurat_rpca_adult_broad_div90_overlay.png"), "adult_broad_div90_overlay", pal)
plot_overlay(ref_df, query_df, file.path(opt$outdir, "FINAL_seurat_rpca_adult_subtypes_div90_black.png"), "adult_subtypes_div90_black", pal)
plot_overlay(ref_df, query_df, file.path(opt$outdir, "FINAL_seurat_rpca_div90_predicted_broad.png"), "div90_predicted_broad", pal)
plot_river(query_df, "div90_class", "predicted_broad", file.path(opt$outdir, "FINAL_seurat_rpca_river_div90_class_to_adult_broad.png"), "Adult broad")
plot_river(query_df, "div90_class", "predicted_subtype", file.path(opt$outdir, "FINAL_seurat_rpca_river_div90_class_to_adult_subtype.png"), "Adult subtype")

audit <- data.frame(
  seurat_version = as.character(packageVersion("Seurat")),
  seuratobject_version = as.character(packageVersion("SeuratObject")),
  bridge_dir = opt$bridge_dir,
  transfer_dir = opt$transfer_dir,
  n_reference_cells = ncol(reference),
  n_query_cells = ncol(query_mapped),
  n_features = length(features),
  npcs = opt$npcs,
  dims = paste0("1:", max(dims_use)),
  reduction = "rpca",
  n_anchors = anchor_count,
  k_weight_requested = opt$k_weight,
  k_weight_used = k_weight_used,
  query_umap_reduction = query_umap_name,
  stringsAsFactors = FALSE
)
write_tsv(audit, file.path(opt$outdir, "seurat_rpca_run_audit.tsv"))
if (opt$save_rds) {
  timed_step("save_reference_rds", saveRDS(reference, file.path(opt$outdir, "seurat", "reference_seurat_rpca_mapping.rds")))
  timed_step("save_query_rds", saveRDS(query_mapped, file.path(opt$outdir, "seurat", "query_seurat_rpca_mapped.rds")))
}

mark_progress("run", "end", paste0("outdir=", opt$outdir))
print(audit)
