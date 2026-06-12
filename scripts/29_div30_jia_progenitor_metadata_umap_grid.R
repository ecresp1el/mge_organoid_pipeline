#!/usr/bin/env Rscript

# Plot DIV30 progenitor reclustering UMAPs by recluster, sample, and cell line.
#
# This is a post hoc visualization step. It reads the saved progenitor-only
# Seurat object and reuses the already computed progenitor_umap coordinates.
# It does not rerun clustering or UMAP.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = "")))
  flush.console()
}

parse_args <- function(args) {
  out <- list(
    `seurat-rds` = NULL,
    outdir = NULL,
    `reduction` = "progenitor_umap",
    `cluster-col` = "div30_progenitor_cluster",
    `sample-col` = "",
    `cell-line-col` = "",
    `sample-map` = "",
    `point-size` = "0.32",
    `alpha` = "0.9",
    help = FALSE
  )
  i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (a %in% c("--help", "-h")) {
      out$help <- TRUE
      i <- i + 1L
      next
    }
    if (!startsWith(a, "--")) stop("Unknown argument: ", a, call. = FALSE)
    key <- substring(a, 3L)
    if (!(key %in% names(out))) stop("Unknown argument: ", a, call. = FALSE)
    if (i == length(args)) stop("Missing value for argument: ", a, call. = FALSE)
    out[[key]] <- args[[i + 1L]]
    i <- i + 2L
  }
  out
}

print_usage <- function() {
  cat(paste(
    "Usage:",
    "  Rscript scripts/29_div30_jia_progenitor_metadata_umap_grid.R --seurat-rds <rds> --outdir <dir>",
    "",
    "Outputs:",
    "  tables/progenitor_umap_metadata_columns.tsv",
    "  tables/progenitor_umap_metadata_counts.tsv",
    "  plots/progenitor_umap_grid_cluster_sample_cellline.png",
    "  plots/progenitor_umap_reclusters_by_sample_id_grid.png",
    "  plots/progenitor_umap_reclusters_by_cell_line_grid.png",
    sep = "\n"
  ))
}

as_num <- function(x, name) {
  value <- suppressWarnings(as.numeric(x))
  if (is.na(value)) stop(name, " must be numeric; got ", x, call. = FALSE)
  value
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

save_plot_pair <- function(plot, stem, width, height, dpi = 300) {
  ggsave(paste0(stem, ".png"), plot = plot, width = width, height = height, dpi = dpi, bg = "white")
  ggsave(paste0(stem, ".pdf"), plot = plot, width = width, height = height, bg = "white")
}

first_existing <- function(meta_cols, requested, candidates, label) {
  if (nzchar(requested)) {
    if (!(requested %in% meta_cols)) stop(label, " column not found: ", requested, call. = FALSE)
    return(requested)
  }
  hit <- candidates[candidates %in% meta_cols]
  if (length(hit) > 0) return(hit[[1]])
  ci <- match(tolower(candidates), tolower(meta_cols))
  ci <- ci[!is.na(ci)]
  if (length(ci) > 0) return(meta_cols[[ci[[1]]]])
  ""
}

metadata_summary <- function(meta) {
  rows <- lapply(colnames(meta), function(cc) {
    x <- as.character(meta[[cc]])
    ux <- unique(x[!is.na(x)])
    data.frame(
      column = cc,
      class = paste(class(meta[[cc]]), collapse = ";"),
      n_non_na = sum(!is.na(meta[[cc]])),
      n_unique = length(ux),
      example_values = paste(head(ux, 8), collapse = " | "),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

metadata_counts <- function(df, cols) {
  cols <- cols[nzchar(cols)]
  rows <- list()
  for (cc in cols) {
    tab <- sort(table(df[[cc]]), decreasing = TRUE)
    rows[[cc]] <- data.frame(
      column = cc,
      value = names(tab),
      n_cells = as.integer(tab),
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0) return(data.frame(column = character(), value = character(), n_cells = integer()))
  do.call(rbind, rows)
}

umap_df <- function(obj, reduction, cluster_col, sample_col, cell_line_col) {
  coords <- as.data.frame(Embeddings(obj, reduction))
  if (ncol(coords) < 2) stop("Reduction has fewer than two dimensions: ", reduction, call. = FALSE)
  colnames(coords)[1:2] <- c("UMAP_1", "UMAP_2")
  coords$cell_id <- rownames(coords)
  meta <- obj@meta.data
  meta$cell_id <- rownames(meta)
  keep <- unique(c("cell_id", cluster_col, sample_col, cell_line_col, "div30_parent_cluster", "div30_parent_progenitor_label"))
  keep <- keep[keep %in% colnames(meta)]
  df <- merge(coords[, c("cell_id", "UMAP_1", "UMAP_2")], meta[, keep, drop = FALSE], by = "cell_id", all.x = TRUE)
  df[[cluster_col]] <- factor(as.character(df[[cluster_col]]), levels = sort(unique(as.character(df[[cluster_col]]))))
  if (nzchar(sample_col)) df[[sample_col]] <- factor(as.character(df[[sample_col]]), levels = names(sort(table(as.character(df[[sample_col]])), decreasing = TRUE)))
  if (nzchar(cell_line_col)) df[[cell_line_col]] <- factor(as.character(df[[cell_line_col]]), levels = names(sort(table(as.character(df[[cell_line_col]])), decreasing = TRUE)))
  df
}

cell_line_from_label <- function(x) {
  out <- sub("_rep[0-9]+$", "", x)
  out <- sub("_old$", "", out)
  out
}

sample_key_variants <- function(x) {
  unique(c(
    x,
    sub("^9583-", "9853-", x),
    sub("^9853-", "9583-", x)
  ))
}

attach_sample_map <- function(df, sample_col, sample_map_path) {
  if (!nzchar(sample_col) || !nzchar(sample_map_path) || !file.exists(sample_map_path)) return(df)
  sample_map <- read.delim(sample_map_path, stringsAsFactors = FALSE, check.names = FALSE)
  needed <- c("DIV", "run_sample_id", "biological_label")
  if (!all(needed %in% colnames(sample_map))) return(df)
  sample_map <- sample_map[sample_map$DIV == "DIV30", needed, drop = FALSE]
  sample_map$cell_line_from_sample_map <- cell_line_from_label(sample_map$biological_label)

  expanded <- do.call(rbind, lapply(seq_len(nrow(sample_map)), function(i) {
    keys <- sample_key_variants(sample_map$run_sample_id[[i]])
    data.frame(
      sample_join_key = keys,
      run_sample_id_from_sample_map = sample_map$run_sample_id[[i]],
      biological_label = sample_map$biological_label[[i]],
      cell_line_from_sample_map = sample_map$cell_line_from_sample_map[[i]],
      stringsAsFactors = FALSE
    )
  }))
  expanded <- expanded[!duplicated(expanded$sample_join_key), , drop = FALSE]
  df$sample_join_key <- as.character(df[[sample_col]])
  out <- merge(df, expanded, by = "sample_join_key", all.x = TRUE)
  out$sample_join_key <- NULL
  out
}

plain_umap <- function(df, color_col, title, point_size, alpha, label_clusters = FALSE, cluster_col = NULL) {
  p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = .data[[color_col]])) +
    geom_point(size = point_size, alpha = alpha, stroke = 0) +
    coord_equal() +
    theme_classic(base_size = 9) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_blank(),
      legend.key.height = unit(0.35, "cm")
    ) +
    labs(title = title, x = "PROGUMAP_1", y = "PROGUMAP_2")
  if (isTRUE(label_clusters) && !is.null(cluster_col)) {
    centers <- aggregate(df[, c("UMAP_1", "UMAP_2")], by = list(cluster = df[[cluster_col]]), FUN = median)
    p <- p + geom_text(data = centers, aes(x = UMAP_1, y = UMAP_2, label = cluster), inherit.aes = FALSE, size = 3, fontface = "bold")
  }
  p
}

facet_background <- function(df, facet_col) {
  facet_levels <- levels(df[[facet_col]])
  if (is.null(facet_levels)) facet_levels <- sort(unique(as.character(df[[facet_col]])))
  pieces <- lapply(facet_levels, function(level) {
    out <- df[, c("UMAP_1", "UMAP_2"), drop = FALSE]
    out[[facet_col]] <- factor(level, levels = facet_levels)
    out
  })
  do.call(rbind, pieces)
}

facet_highlight_grid <- function(df, facet_col, cluster_col, title, point_size, alpha) {
  bg <- facet_background(df, facet_col)
  ggplot() +
    geom_point(data = bg, aes(x = UMAP_1, y = UMAP_2), color = "grey86", size = point_size * 0.85, alpha = 0.28, stroke = 0) +
    geom_point(data = df, aes(x = UMAP_1, y = UMAP_2, color = .data[[cluster_col]]), size = point_size, alpha = alpha, stroke = 0) +
    facet_wrap(stats::as.formula(paste("~", facet_col))) +
    coord_equal() +
    theme_classic(base_size = 8) +
    theme(
      strip.text = element_text(face = "bold", size = 7),
      legend.title = element_blank(),
      legend.position = "right",
      axis.text = element_blank(),
      axis.ticks = element_blank()
    ) +
    labs(title = title, x = "PROGUMAP_1", y = "PROGUMAP_2")
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(opt$help)) {
  print_usage()
  quit(save = "no", status = 0)
}
if (is.null(opt$`seurat-rds`) || !nzchar(opt$`seurat-rds`)) stop("--seurat-rds is required", call. = FALSE)
if (is.null(opt$outdir) || !nzchar(opt$outdir)) stop("--outdir is required", call. = FALSE)

cfg <- list(
  seurat_rds = opt$`seurat-rds`,
  outdir = opt$outdir,
  table_dir = file.path(opt$outdir, "tables"),
  plot_dir = file.path(opt$outdir, "plots"),
  reduction = opt$reduction,
  cluster_col = opt$`cluster-col`,
  sample_col_requested = opt$`sample-col`,
  cell_line_col_requested = opt$`cell-line-col`,
  sample_map = opt$`sample-map`,
  point_size = as_num(opt$`point-size`, "point-size"),
  alpha = as_num(opt$alpha, "alpha")
)
dir.create(cfg$table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$plot_dir, recursive = TRUE, showWarnings = FALSE)

log_msg("Reading Seurat object: ", cfg$seurat_rds)
obj <- readRDS(cfg$seurat_rds)
if (!(cfg$reduction %in% Reductions(obj))) stop("Reduction not found: ", cfg$reduction, call. = FALSE)
if (!(cfg$cluster_col %in% colnames(obj@meta.data))) stop("Cluster column not found: ", cfg$cluster_col, call. = FALSE)

meta_cols <- colnames(obj@meta.data)
sample_col <- first_existing(
  meta_cols,
  cfg$sample_col_requested,
  c("sample_id", "sample", "Sample", "orig.ident", "sample_name", "sampleID", "sample.id"),
  "sample"
)
cell_line_col <- first_existing(
  meta_cols,
  cfg$cell_line_col_requested,
  c("cell_line", "cellline", "cell.line", "line", "line_id", "clone", "genotype", "cell_line_id"),
  "cell line"
)
if (!nzchar(sample_col)) log_msg("No sample column auto-detected.")
if (!nzchar(cell_line_col)) log_msg("No cell-line column auto-detected.")
log_msg("Using sample_col=", ifelse(nzchar(sample_col), sample_col, "none"), "; cell_line_col=", ifelse(nzchar(cell_line_col), cell_line_col, "none"))

meta_summary <- metadata_summary(obj@meta.data)
meta_summary$selected_for_plot <- meta_summary$column %in% c(cfg$cluster_col, sample_col, cell_line_col)
write_tsv(meta_summary, file.path(cfg$table_dir, "progenitor_umap_metadata_columns.tsv"))

df <- umap_df(obj, cfg$reduction, cfg$cluster_col, sample_col, cell_line_col)
df <- attach_sample_map(df, sample_col, cfg$sample_map)
if (!nzchar(cell_line_col) && "cell_line_from_sample_map" %in% colnames(df) && any(!is.na(df$cell_line_from_sample_map))) {
  cell_line_col <- "cell_line_from_sample_map"
  df[[cell_line_col]] <- factor(
    as.character(df[[cell_line_col]]),
    levels = names(sort(table(as.character(df[[cell_line_col]])), decreasing = TRUE))
  )
  log_msg("Using sample-map-derived cell_line_col=", cell_line_col)
}
write_tsv(metadata_counts(df, c(cfg$cluster_col, sample_col, cell_line_col)), file.path(cfg$table_dir, "progenitor_umap_metadata_counts.tsv"))
write_tsv(df, file.path(cfg$table_dir, "progenitor_umap_metadata_plot_cells.tsv"))

plots <- list(
  plain_umap(df, cfg$cluster_col, "Reclusters", cfg$point_size, cfg$alpha, TRUE, cfg$cluster_col)
)
if (nzchar(sample_col)) plots[[length(plots) + 1L]] <- plain_umap(df, sample_col, sample_col, cfg$point_size, cfg$alpha)
if (nzchar(cell_line_col)) plots[[length(plots) + 1L]] <- plain_umap(df, cell_line_col, cell_line_col, cfg$point_size, cfg$alpha)

combined <- patchwork::wrap_plots(plots, nrow = 1, guides = "collect") &
  theme(legend.position = "right")
save_plot_pair(combined, file.path(cfg$plot_dir, "progenitor_umap_grid_cluster_sample_cellline"), width = 15, height = 5)

if (nzchar(sample_col)) {
  n_sample <- length(unique(df[[sample_col]]))
  h <- max(5, ceiling(n_sample / 4) * 3)
  p_sample_grid <- facet_highlight_grid(df, sample_col, cfg$cluster_col, "Reclusters within each sample_id", cfg$point_size, cfg$alpha)
  save_plot_pair(p_sample_grid, file.path(cfg$plot_dir, "progenitor_umap_reclusters_by_sample_id_grid"), width = 14, height = h)
}

if (nzchar(cell_line_col)) {
  n_line <- length(unique(df[[cell_line_col]]))
  h <- max(5, ceiling(n_line / 4) * 3)
  p_line_grid <- facet_highlight_grid(df, cell_line_col, cfg$cluster_col, "Reclusters within each cell line", cfg$point_size, cfg$alpha)
  save_plot_pair(p_line_grid, file.path(cfg$plot_dir, "progenitor_umap_reclusters_by_cell_line_grid"), width = 14, height = h)
}

write_tsv(
  data.frame(
    completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    seurat_rds = cfg$seurat_rds,
    reduction = cfg$reduction,
    cluster_col = cfg$cluster_col,
    sample_col = sample_col,
    cell_line_col = cell_line_col,
    sample_map = cfg$sample_map,
    n_cells = nrow(df),
    stringsAsFactors = FALSE
  ),
  file.path(cfg$table_dir, "progenitor_umap_metadata_plot_complete.tsv")
)
log_msg("Metadata UMAP plots complete: ", cfg$outdir)
