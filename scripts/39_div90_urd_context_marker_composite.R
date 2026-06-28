#!/usr/bin/env Rscript

# Composite DIV90 URD figure:
#   top row: cluster UMAP + cluster lineage tree, pseudotime UMAP + pseudotime lineage tree
#   bottom row: six marker-expression lineage overlays with the final grey-floor/blue color logic

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = " ")))
  flush.console()
}

parse_args <- function(args) {
  out <- list(
    `tree-rds` = NULL,
    outdir = NULL,
    genes = "HES1,NKX2-1,LHX6,LHX8,CRABP1,KCNC1",
    `gene-labels` = "Hes1,Nkx2.1,Lhx6,Lhx8,Crabp1,Kcnc1",
    `annotation-col` = "cluster_number_name",
    `pseudotime-name` = "div90_jia_rootscore_root",
    `expression-color-floor` = "1",
    `vmax-quantile` = "0.99",
    `point-size` = "0.26",
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

required <- c("Matrix", "ggplot2", "cowplot", "scales")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Missing required R package(s): ", paste(missing, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
  library(cowplot)
})

as_num <- function(x, name) {
  value <- suppressWarnings(as.numeric(x))
  if (is.na(value)) stop(name, " must be numeric; got ", x, call. = FALSE)
  value
}

split_csv <- function(x) {
  if (!nzchar(x)) return(character())
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

write_tsv <- function(x, path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, open = "wt") else file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.table(x, con, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

save_plot_set <- function(plot, prefix, width, height, dpi = 300) {
  png_path <- paste0(prefix, ".png")
  pdf_path <- paste0(prefix, ".pdf")
  svg_path <- paste0(prefix, ".svg")
  ggplot2::ggsave(png_path, plot, width = width, height = height, dpi = dpi, bg = "white")
  ggplot2::ggsave(pdf_path, plot, width = width, height = height, device = grDevices::cairo_pdf, bg = "white")
  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(svg_path, width = width, height = height)
    print(plot)
    grDevices::dev.off()
  } else {
    unlink(svg_path)
    log_msg("Skipping SVG because svglite is unavailable; editable vector text is available in PDF: ", pdf_path)
  }
}

repair_logupx_dimnames <- function(object) {
  expr <- object@logupx.data
  counts <- object@count.data
  if (is.null(rownames(expr)) && !is.null(rownames(counts)) && nrow(expr) == nrow(counts)) {
    rownames(expr) <- rownames(counts)
  }
  if (is.null(colnames(expr)) && !is.null(colnames(counts)) && ncol(expr) == ncol(counts)) {
    colnames(expr) <- colnames(counts)
  }
  expr
}

palette_for_clusters <- function(n) {
  base <- c(
    "#1f78b4", "#33a02c", "#6a3d9a", "#b15928", "#a6cee3",
    "#fb9a99", "#999999", "#bdbdbd", "#e31a1c", "#cab2d6",
    "#fdbf6f", "#ff7f00", "#ffff99", "#8dd3c7", "#bebada"
  )
  rep(base, length.out = n)
}

cluster_id_from_label <- function(x) {
  suppressWarnings(as.integer(sub("^([0-9]+).*", "\\1", as.character(x))))
}

clean_biology_name <- function(x) {
  sub("^[0-9]+\\s*-\\s*", "", as.character(x))
}

wrapped <- function(x, width = 34) {
  vapply(strwrap(x, width = width, simplify = FALSE), paste, character(1), collapse = "\n")
}

extract_pseudotime <- function(object, pseudotime_name) {
  if (!is.null(object@pseudotime) && ncol(object@pseudotime) > 0) {
    if (nzchar(pseudotime_name) && pseudotime_name %in% colnames(object@pseudotime)) {
      values <- as.numeric(object@pseudotime[, pseudotime_name])
      names(values) <- rownames(object@pseudotime)
      return(values)
    }
    values <- as.numeric(object@pseudotime[, 1])
    names(values) <- rownames(object@pseudotime)
    return(values)
  }
  values <- as.numeric(object@tree$pseudotime)
  names(values) <- names(object@tree$pseudotime)
  values
}

orient_tree_left_to_right <- function(layout, cells) {
  # URD stores pseudotime on the y-axis; rotate to put root/trunk left and tips right.
  layout$x1_plot <- layout$y1
  layout$y1_plot <- -layout$x1
  layout$x2_plot <- layout$y2
  layout$y2_plot <- -layout$x2
  cells$x_plot <- cells$y
  cells$y_plot <- -cells$x
  layout$x1 <- layout$x1_plot
  layout$y1 <- layout$y1_plot
  layout$x2 <- layout$x2_plot
  layout$y2 <- layout$y2_plot
  cells$x <- cells$x_plot
  cells$y <- cells$y_plot
  list(layout = layout, cells = cells)
}

context_data <- function(object, annotation_col, pseudotime_name) {
  meta <- as.data.frame(object@meta, stringsAsFactors = FALSE)
  meta$cell <- rownames(meta)
  if (!all(c("UMAP_1", "UMAP_2") %in% colnames(meta))) {
    stop("URD metadata is missing UMAP_1/UMAP_2.", call. = FALSE)
  }
  if (!(annotation_col %in% colnames(meta))) stop("Missing annotation column: ", annotation_col, call. = FALSE)
  pt <- extract_pseudotime(object, pseudotime_name)
  meta$pseudotime <- pt[meta$cell]
  meta$annotation <- as.character(meta[[annotation_col]])
  meta$cluster_id <- cluster_id_from_label(meta$annotation)
  mapping <- aggregate(cell ~ cluster_id + annotation, meta, length)
  colnames(mapping)[colnames(mapping) == "cell"] <- "n_cells"
  mapping$biology_name <- clean_biology_name(mapping$annotation)
  mapping <- mapping[order(mapping$cluster_id, mapping$annotation), , drop = FALSE]
  mapping$cluster_factor <- factor(mapping$annotation, levels = mapping$annotation)
  colors <- setNames(palette_for_clusters(nrow(mapping)), mapping$annotation)

  layout <- as.data.frame(object@tree$tree.layout, stringsAsFactors = FALSE)
  cells <- as.data.frame(object@tree$cell.layout, stringsAsFactors = FALSE)
  cells$annotation <- meta[cells$cell, annotation_col]
  cells$annotation <- factor(as.character(cells$annotation), levels = mapping$annotation)
  cells$pseudotime <- pt[cells$cell]
  oriented <- orient_tree_left_to_right(layout, cells)
  layout <- oriented$layout
  cells <- oriented$cells

  list(meta = meta, mapping = mapping, colors = colors, layout = layout, cells = cells)
}

base_tree_plot <- function(layout) {
  geom_segment(
    data = layout,
    aes(x = x1, y = y1, xend = x2, yend = y2),
    linewidth = 0.22,
    color = "grey55",
    alpha = 0.85
  )
}

tree_limits <- function(layout, cells, right_extra = 0.04) {
  x_limits <- range(c(layout$x1, layout$x2, cells$x), na.rm = TRUE)
  y_limits <- range(c(layout$y1, layout$y2, cells$y), na.rm = TRUE)
  x_pad <- diff(x_limits) * 0.04
  y_pad <- diff(y_limits) * 0.04
  list(x = x_limits + c(-x_pad, diff(x_limits) * right_extra), y = y_limits + c(-y_pad, y_pad))
}

cluster_umap_plot <- function(df, mapping, colors) {
  df <- df[is.finite(df$UMAP_1) & is.finite(df$UMAP_2) & !is.na(df$annotation), , drop = FALSE]
  df$annotation <- factor(df$annotation, levels = mapping$annotation)
  centers <- aggregate(cbind(UMAP_1, UMAP_2) ~ cluster_id + annotation, df, median)
  centers <- centers[order(centers$cluster_id), , drop = FALSE]
  ggplot(df, aes(UMAP_1, UMAP_2, color = annotation)) +
    geom_point(size = 0.13, alpha = 0.7) +
    geom_label(
      data = centers,
      aes(UMAP_1, UMAP_2, label = cluster_id),
      inherit.aes = FALSE,
      size = 2.7,
      fontface = "bold",
      linewidth = 0.18,
      fill = "white",
      color = "black"
    ) +
    scale_color_manual(values = colors, guide = "none") +
    coord_equal() +
    theme_void(base_size = 8) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 9)) +
    labs(title = "DIV90 UMAP: clusters")
}

cluster_key_plot <- function(mapping) {
  key <- sprintf("%s = %s", mapping$cluster_id, wrapped(mapping$biology_name, 24))
  cowplot::ggdraw() +
    cowplot::draw_label("DIV90 cluster names", x = 0, y = 0.99, hjust = 0, vjust = 1, fontface = "bold", size = 8.5) +
    cowplot::draw_label(paste(key, collapse = "\n"), x = 0, y = 0.93, hjust = 0, vjust = 1, size = 4.6, lineheight = 0.84)
}

cluster_tree_tip_labels <- function(cells, mapping) {
  rows <- lapply(mapping$annotation, function(annotation) {
    df <- cells[as.character(cells$annotation) == annotation & is.finite(cells$x) & is.finite(cells$y), , drop = FALSE]
    if (nrow(df) == 0) return(NULL)
    x_cut <- stats::quantile(df$x, probs = 0.86, na.rm = TRUE, names = FALSE)
    tip_df <- df[df$x >= x_cut, , drop = FALSE]
    if (nrow(tip_df) == 0) tip_df <- df[which.max(df$x), , drop = FALSE]
    m <- mapping[mapping$annotation == annotation, , drop = FALSE]
    data.frame(
      annotation = annotation,
      cluster_id = m$cluster_id[[1]],
      label = paste0(m$cluster_id[[1]], " - ", wrapped(m$biology_name[[1]], 31)),
      x = max(tip_df$x, na.rm = TRUE),
      y = stats::median(tip_df$y, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  labels <- do.call(rbind, rows)
  if (is.null(labels) || nrow(labels) == 0) return(data.frame())
  labels <- labels[order(labels$y), , drop = FALSE]
  y_range <- diff(range(cells$y, na.rm = TRUE))
  min_gap <- y_range / 14
  for (i in seq_len(nrow(labels))[-1]) {
    if (labels$y[[i]] - labels$y[[i - 1]] < min_gap) labels$y[[i]] <- labels$y[[i - 1]] + min_gap
  }
  labels
}

cluster_tree_plot <- function(layout, cells, mapping, colors) {
  lim <- tree_limits(layout, cells, right_extra = 0.72)
  labels <- cluster_tree_tip_labels(cells, mapping)
  x_span <- diff(range(c(layout$x1, layout$x2, cells$x), na.rm = TRUE))
  if (nrow(labels) > 0) labels$x <- labels$x + x_span * 0.025
  ggplot() +
    base_tree_plot(layout) +
    geom_point(data = cells, aes(x = x, y = y, color = annotation), size = 0.34, alpha = 0.88) +
    geom_text(
      data = labels,
      aes(x = x, y = y, label = label, color = annotation),
      hjust = 0,
      size = 1.95,
      lineheight = 0.9,
      fontface = "plain",
      show.legend = FALSE
    ) +
    scale_color_manual(values = colors, guide = "none", na.value = "grey85") +
    coord_cartesian(xlim = lim$x, ylim = lim$y, expand = FALSE) +
    theme_void(base_size = 8) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 9)) +
    labs(title = "DIV90 URD tree: clusters")
}

pseudotime_umap_plot <- function(df, pseudotime_limits) {
  df <- df[is.finite(df$UMAP_1) & is.finite(df$UMAP_2), , drop = FALSE]
  ggplot(df, aes(UMAP_1, UMAP_2, color = pseudotime)) +
    geom_point(size = 0.13, alpha = 0.8) +
    coord_equal() +
    scale_color_viridis_c(name = "Pseudotime", limits = pseudotime_limits, oob = scales::squish, na.value = "grey85") +
    theme_void(base_size = 8) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 9),
      legend.key.height = grid::unit(16, "pt"),
      legend.key.width = grid::unit(5, "pt"),
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 5)
    ) +
    labs(title = "DIV90 UMAP: pseudotime")
}

pseudotime_tree_plot <- function(layout, cells, pseudotime_limits) {
  lim <- tree_limits(layout, cells)
  ggplot() +
    base_tree_plot(layout) +
    geom_point(data = cells, aes(x = x, y = y, color = pseudotime), size = 0.28, alpha = 0.85) +
    scale_color_viridis_c(name = "Pseudotime", limits = pseudotime_limits, oob = scales::squish, na.value = "grey85") +
    coord_cartesian(xlim = lim$x, ylim = lim$y, expand = FALSE) +
    theme_void(base_size = 8) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 9),
      legend.key.height = grid::unit(16, "pt"),
      legend.key.width = grid::unit(5, "pt"),
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 5)
    ) +
    labs(title = "DIV90 URD tree: pseudotime")
}

marker_spec <- function(genes, labels) {
  genes <- split_csv(genes)
  labels <- split_csv(labels)
  if (length(labels) == 0) labels <- genes
  if (length(genes) != length(labels)) {
    stop("--genes and --gene-labels must contain the same number of comma-separated values.", call. = FALSE)
  }
  data.frame(gene = genes, display_label = labels, display_order = seq_along(genes), stringsAsFactors = FALSE)
}

build_marker_table <- function(expr, cells, spec) {
  rows <- lapply(seq_len(nrow(spec)), function(i) {
    gene <- spec$gene[[i]]
    present <- gene %in% rownames(expr)
    values <- if (present) as.numeric(expr[gene, cells]) else rep(NA_real_, length(cells))
    data.frame(
      gene = gene,
      display_label = spec$display_label[[i]],
      display_order = spec$display_order[[i]],
      cell = cells,
      expression_logupx = values,
      gene_present = present,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

marker_summary <- function(marker_df) {
  rows <- lapply(split(marker_df, marker_df$gene), function(df) {
    values <- df$expression_logupx
    has_values <- any(!is.na(values))
    data.frame(
      display_order = df$display_order[[1]],
      gene = df$gene[[1]],
      display_label = df$display_label[[1]],
      gene_present = df$gene_present[[1]],
      n_tree_cells = nrow(df),
      pct_above_floor = if (has_values) mean(values > 1, na.rm = TRUE) else NA_real_,
      mean_logupx = if (has_values) mean(values, na.rm = TRUE) else NA_real_,
      max_logupx = if (has_values) max(values, na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$display_order), , drop = FALSE]
}

expression_color_vmax <- function(values, floor_value, quantile_value) {
  values <- values[is.finite(values)]
  positive <- values[values > 0]
  if (length(positive) == 0) return(floor_value)
  max(as.numeric(stats::quantile(positive, probs = quantile_value, names = FALSE, na.rm = TRUE)), floor_value)
}

marker_tree_plot <- function(layout, cells, marker_df, gene, point_size, color_floor, vmax_quantile) {
  df <- marker_df[marker_df$gene == gene, , drop = FALSE]
  plot_df <- merge(cells[, c("cell", "x", "y")], df[, c("cell", "expression_logupx", "gene_present", "display_label")], by = "cell", all.x = TRUE)
  plot_df <- plot_df[order(plot_df$expression_logupx, na.last = TRUE), , drop = FALSE]
  plot_df$expression_color_value <- ifelse(
    is.finite(plot_df$expression_logupx) & plot_df$expression_logupx > color_floor,
    plot_df$expression_logupx,
    NA_real_
  )
  lim <- tree_limits(layout, plot_df)
  max_expr <- expression_color_vmax(plot_df$expression_logupx, color_floor, vmax_quantile)
  color_vmax <- max(max_expr, color_floor + 1e-6)
  floor_fraction <- max(0, min(1, color_floor / color_vmax))

  p <- ggplot() +
    base_tree_plot(layout) +
    geom_point(data = plot_df, aes(x = x, y = y), size = point_size * 0.72, color = "#d0d0d0", alpha = 0.7) +
    geom_point(
      data = plot_df[is.finite(plot_df$expression_color_value), , drop = FALSE],
      aes(x = x, y = y, color = expression_color_value),
      size = point_size,
      alpha = 0.9
    ) +
    scale_color_gradientn(
      colours = c("#d0d0d0", "#d0d0d0", "#0000ff"),
      values = c(0, floor_fraction, 1),
      limits = c(0, color_vmax),
      oob = scales::squish,
      name = "logUPX",
      breaks = c(0, color_vmax),
      labels = c("0", formatC(max_expr, format = "fg", digits = 2))
    ) +
    coord_cartesian(xlim = lim$x, ylim = lim$y, expand = FALSE) +
    theme_void(base_size = 8) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 9),
      legend.key.height = grid::unit(14, "pt"),
      legend.key.width = grid::unit(5, "pt"),
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 5)
    ) +
    labs(title = df$display_label[[1]])
  if (!isTRUE(df$gene_present[[1]])) {
    p <- p + annotate("text", x = mean(lim$x), y = mean(lim$y), label = "missing", size = 3.5, color = "grey20")
  }
  p
}

print_usage <- function() {
  cat("Usage: Rscript scripts/39_div90_urd_context_marker_composite.R --tree-rds <tree.rds> --outdir <dir>\n")
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(opt$help)) {
  print_usage()
  quit(save = "no", status = 0)
}
if (is.null(opt$`tree-rds`) || !nzchar(opt$`tree-rds`)) stop("--tree-rds is required", call. = FALSE)
if (is.null(opt$outdir) || !nzchar(opt$outdir)) stop("--outdir is required", call. = FALSE)

cfg <- list(
  tree_rds = opt$`tree-rds`,
  outdir = opt$outdir,
  plot_dir = file.path(opt$outdir, "plots"),
  table_dir = file.path(opt$outdir, "tables"),
  annotation_col = opt$`annotation-col`,
  pseudotime_name = opt$`pseudotime-name`,
  expression_color_floor = as_num(opt$`expression-color-floor`, "expression-color-floor"),
  vmax_quantile = as_num(opt$`vmax-quantile`, "vmax-quantile"),
  point_size = as_num(opt$`point-size`, "point-size")
)
if (cfg$vmax_quantile <= 0 || cfg$vmax_quantile > 1) stop("vmax-quantile must be in (0, 1]", call. = FALSE)
dir.create(cfg$plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$table_dir, recursive = TRUE, showWarnings = FALSE)

log_msg("Reading URD tree object: ", cfg$tree_rds)
urd <- readRDS(cfg$tree_rds)
expr <- repair_logupx_dimnames(urd)
ctx <- context_data(urd, cfg$annotation_col, cfg$pseudotime_name)
spec <- marker_spec(opt$genes, opt$`gene-labels`)
marker_df <- build_marker_table(expr, ctx$cells$cell, spec)

write_tsv(spec, file.path(cfg$table_dir, "context_marker_gene_order.tsv"))
write_tsv(ctx$mapping[, c("cluster_id", "annotation", "biology_name", "n_cells")], file.path(cfg$table_dir, "context_cluster_key.tsv"))
write_tsv(marker_summary(marker_df), file.path(cfg$table_dir, "context_marker_expression_summary.tsv"))

missing_genes <- spec$gene[!spec$gene %in% rownames(expr)]
if (length(missing_genes) > 0) {
  log_msg("Missing requested marker(s): ", paste(missing_genes, collapse = ", "))
} else {
  log_msg("All requested marker genes are present.")
}

pseudotime_limits <- range(c(ctx$meta$pseudotime, ctx$cells$pseudotime), finite = TRUE, na.rm = TRUE)

cluster_block <- cowplot::plot_grid(
  cluster_umap_plot(ctx$meta, ctx$mapping, ctx$colors),
  cluster_tree_plot(ctx$layout, ctx$cells, ctx$mapping, ctx$colors),
  nrow = 1,
  rel_widths = c(0.9, 1.85)
)
pseudotime_block <- cowplot::plot_grid(
  pseudotime_umap_plot(ctx$meta, pseudotime_limits),
  pseudotime_tree_plot(ctx$layout, ctx$cells, pseudotime_limits),
  nrow = 1,
  rel_widths = c(1, 1.08)
)
top_row <- cowplot::plot_grid(cluster_block, pseudotime_block, nrow = 1, rel_widths = c(3.05, 2.05))

marker_plots <- lapply(spec$gene, function(gene) {
  marker_tree_plot(ctx$layout, ctx$cells, marker_df, gene, cfg$point_size, cfg$expression_color_floor, cfg$vmax_quantile)
})
bottom_row <- cowplot::plot_grid(plotlist = marker_plots, nrow = 1, rel_widths = rep(1, length(marker_plots)))

title <- cowplot::ggdraw() +
  cowplot::draw_label(
    "DIV90 URD context and marker-expression lineage overlays",
    x = 0,
    y = 0.55,
    hjust = 0,
    vjust = 0.5,
    fontface = "bold",
    size = 14
  )
composite <- cowplot::plot_grid(
  title,
  top_row,
  bottom_row,
  ncol = 1,
  rel_heights = c(0.08, 1.02, 0.98)
)

save_plot_set(
  composite,
  file.path(cfg$plot_dir, "div90_urd_context_marker_composite"),
  width = 22,
  height = 11,
  dpi = 300
)

render_status <- data.frame(
  field = c(
    "rendered_at",
    "tree_rds",
    "annotation_col",
    "pseudotime_name",
    "marker_panel",
    "marker_color_map",
    "expression_color_floor",
    "vmax_rule",
    "cluster_palette",
    "tree_orientation",
    "pseudotime_scale",
    "cluster_tree_tip_labels",
    "editable_vector_outputs"
  ),
  value = c(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    cfg$tree_rds,
    cfg$annotation_col,
    cfg$pseudotime_name,
    paste(spec$display_label, collapse = ","),
    "cross-study-compatible grey floor to #0000ff blue",
    as.character(cfg$expression_color_floor),
    paste0("q", cfg$vmax_quantile, " positive expression per gene"),
    "palette_for_clusters from scripts/26_div90_umap_cluster_label_audit.R",
    "left-to-right rotated tree; tips at right",
    paste0("shared UMAP/tree limits ", paste(signif(pseudotime_limits, 4), collapse = " to ")),
    "cluster number and wrapped cluster name labels drawn at cluster tree tips",
    if (requireNamespace("svglite", quietly = TRUE)) "pdf,svg" else "pdf"
  ),
  stringsAsFactors = FALSE
)
write_tsv(render_status, file.path(cfg$outdir, "context_marker_composite_render_status.tsv"))
log_msg("Wrote composite outputs under: ", cfg$plot_dir)
