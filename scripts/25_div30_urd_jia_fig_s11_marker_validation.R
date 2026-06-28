#!/usr/bin/env Rscript

# Jia Fig. S11-style marker-expression validation on the DIV30 URD tree.
#
# This is a marker-expression figure only. It intentionally does not use
# branchpoint DE genes, lineage program averages, or marker substitutions.
#
# Panel A marker order:
#   HES1 | CACNA1E | DLX2 | DCX
#
# Panel B marker order and labels:
#   LHX8   -> LHX8/ISL1      -> Subpallial cholinergic interneuron
#   NR2F1  -> NR2F1/NR2F2    -> Subpallial GABAergic interneuron
#   EPHA5  -> EPHA5/MEF2C    -> Cortical GABAergic interneuron
#   MEF2C  -> LHX6/NFIA      -> Cortical interneuron lineage
#   CRABP1 -> CRABP1/ANGPT2  -> Subpallial GABAergic interneuron
#
# Expression source:
#   URD count.data transformed at plot time to log1p(CP10K):
#   log1p(gene_count / total_cell_counts * 10000).

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = " ")))
  flush.console()
}

parse_args <- function(args) {
  out <- list(
    `tree-rds` = NULL,
    outdir = NULL,
    genes = "",
    `gene-labels` = "",
    `panel-title` = "",
    `point-size` = "0.28",
    `expression-color-floor` = "1",
    `vmax-quantile` = "0.99",
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

required <- c("Matrix", "ggplot2", "cowplot")
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

write_tsv <- function(x, path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, open = "wt") else file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.table(x, con, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

save_plot_set <- function(plot, png_path, pdf_path, svg_path, width, height, dpi = 300) {
  ggsave(png_path, plot, width = width, height = height, dpi = dpi, bg = "white")
  ggsave(pdf_path, plot, width = width, height = height, device = grDevices::cairo_pdf, bg = "white")
  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(svg_path, width = width, height = height)
    print(plot)
    grDevices::dev.off()
  } else {
    unlink(svg_path)
    log_msg("Skipping SVG because svglite is not available; editable vector text is available in PDF: ", pdf_path)
  }
}

repair_count_dimnames <- function(object) {
  counts <- object@count.data
  expr <- object@logupx.data
  if (is.null(rownames(counts)) && !is.null(rownames(expr)) && nrow(counts) == nrow(expr)) {
    rownames(counts) <- rownames(expr)
  }
  if (is.null(colnames(counts)) && !is.null(colnames(expr)) && ncol(counts) == ncol(expr)) {
    colnames(counts) <- colnames(expr)
  }
  counts
}

marker_values_log1p_cp10k <- function(counts, totals, gene, cells) {
  if (!(gene %in% rownames(counts))) return(rep(NA_real_, length(cells)))
  totals <- totals[cells]
  values <- as.numeric(counts[gene, cells])
  out <- rep(NA_real_, length(cells))
  ok <- is.finite(totals) & totals > 0
  out[ok] <- log1p(values[ok] / totals[ok] * 10000)
  out
}

marker_spec <- function() {
  data.frame(
    panel = c(rep("A_developmental_progression", 4), rep("B_inhibitory_lineages", 5)),
    display_order = c(seq_len(4), seq_len(5)),
    category = c("VZ RGC", "SVZ RGC", "IPC", "Newborn neuron", "LHX8 lineage", "NR2F1 lineage", "EPHA5 lineage", "MEF2C lineage", "CRABP1 lineage"),
    gene = c("HES1", "CACNA1E", "DLX2", "DCX", "LHX8", "NR2F1", "EPHA5", "MEF2C", "CRABP1"),
    display_label = c("HES1", "CACNA1E", "DLX2", "DCX", "LHX8", "NR2F1", "EPHA5", "MEF2C", "CRABP1"),
    lineage_pair = c(NA, NA, NA, NA, "LHX8/ISL1", "NR2F1/NR2F2", "EPHA5/MEF2C", "LHX6/NFIA", "CRABP1/ANGPT2"),
    terminal_identity = c(NA, NA, NA, NA, "Subpallial cholinergic interneuron", "Subpallial GABAergic interneuron", "Cortical GABAergic interneuron", "Cortical interneuron lineage", "Subpallial GABAergic interneuron"),
    stringsAsFactors = FALSE
  )
}

split_csv <- function(x) {
  if (!nzchar(x)) return(character())
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

custom_marker_spec <- function(genes, labels = character()) {
  genes <- split_csv(genes)
  if (length(genes) == 0) return(NULL)
  labels <- split_csv(labels)
  if (length(labels) == 0) labels <- genes
  if (length(labels) != length(genes)) {
    stop("--gene-labels must contain the same number of comma-separated values as --genes", call. = FALSE)
  }
  data.frame(
    panel = rep("custom_marker_panel", length(genes)),
    display_order = seq_along(genes),
    category = labels,
    gene = genes,
    display_label = labels,
    lineage_pair = NA_character_,
    terminal_identity = NA_character_,
    stringsAsFactors = FALSE
  )
}

build_marker_table <- function(counts, cells, spec) {
  totals <- Matrix::colSums(counts)
  rows <- lapply(seq_len(nrow(spec)), function(i) {
    gene <- spec$gene[[i]]
    present <- gene %in% rownames(counts)
    values <- marker_values_log1p_cp10k(counts, totals, gene, cells)
    data.frame(
      panel = spec$panel[[i]],
      display_order = spec$display_order[[i]],
      category = spec$category[[i]],
      lineage_pair = spec$lineage_pair[[i]],
      terminal_identity = spec$terminal_identity[[i]],
      gene = gene,
      display_label = spec$display_label[[i]],
      cell = cells,
      expression_log1p_cp10k = values,
      gene_present = present,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

marker_summary <- function(marker_df) {
  rows <- lapply(split(marker_df, marker_df$gene), function(df) {
    values <- df$expression_log1p_cp10k
    has_values <- any(!is.na(values))
    data.frame(
      panel = df$panel[[1]],
      display_order = df$display_order[[1]],
      category = df$category[[1]],
      lineage_pair = df$lineage_pair[[1]],
      terminal_identity = df$terminal_identity[[1]],
      gene = df$gene[[1]],
      display_label = df$display_label[[1]],
      gene_present = df$gene_present[[1]],
      n_cells = nrow(df),
      pct_expressed = if (has_values) mean(values > 0, na.rm = TRUE) else NA_real_,
      mean_log1p_cp10k = if (has_values) mean(values, na.rm = TRUE) else NA_real_,
      median_log1p_cp10k = if (has_values) median(values, na.rm = TRUE) else NA_real_,
      max_log1p_cp10k = if (has_values) max(values, na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$panel, out$display_order), , drop = FALSE]
}

expression_color_vmax <- function(values, floor_value, quantile_value) {
  values <- values[is.finite(values)]
  positive <- values[values > 0]
  if (length(positive) == 0) return(floor_value)
  max(as.numeric(stats::quantile(positive, probs = quantile_value, names = FALSE, na.rm = TRUE)), floor_value)
}

expression_floor_palette <- function(color_vmax, floor_value, background = "#d0d0d0", high = "#0000ff") {
  color_vmax <- max(color_vmax, floor_value + 1e-6)
  floor_fraction <- max(0, min(1, floor_value / color_vmax))
  scales::gradient_n_pal(
    colours = c(background, background, high),
    values = c(0, floor_fraction, 1)
  )
}

tree_marker_plot <- function(layout, cells, marker_df, gene, panel_title, point_size, color_floor, vmax_quantile, show_legend = FALSE) {
  # Match the cross-study marker-expression logic: values at/below the floor
  # are drawn as background grey, and the blue scale clips at per-gene q99.
  df <- marker_df[marker_df$gene == gene, , drop = FALSE]
  plot_df <- merge(cells, df[, c("cell", "expression_log1p_cp10k", "gene_present")], by = "cell", all.x = TRUE)
  plot_df <- plot_df[order(plot_df$expression_log1p_cp10k, na.last = TRUE), , drop = FALSE]
  plot_df$expression_color_value <- ifelse(
    is.finite(plot_df$expression_log1p_cp10k) & plot_df$expression_log1p_cp10k > color_floor,
    plot_df$expression_log1p_cp10k,
    NA_real_
  )
  x_limits <- range(c(layout$x1, layout$x2, cells$x), na.rm = TRUE)
  y_limits <- range(c(layout$y1, layout$y2, cells$y), na.rm = TRUE)
  x_pad <- diff(x_limits) * 0.04
  y_pad <- diff(y_limits) * 0.04
  x_limits <- x_limits + c(-x_pad, x_pad)
  y_limits <- y_limits + c(-y_pad, y_pad)
  max_expr <- expression_color_vmax(plot_df$expression_log1p_cp10k, color_floor, vmax_quantile)
  color_vmax <- max(max_expr, color_floor + 1e-6)

  p <- ggplot() +
    geom_segment(data = layout, aes(x = x1, y = y1, xend = x2, yend = y2), linewidth = 0.25, color = "grey60") +
    geom_point(data = plot_df, aes(x = x, y = y), size = point_size * 0.72, color = "#d0d0d0", alpha = 0.7) +
    geom_point(
      data = plot_df[is.finite(plot_df$expression_color_value), , drop = FALSE],
      aes(x = x, y = y, color = expression_color_value),
      size = point_size,
      alpha = 0.9
    ) +
    scale_color_gradientn(
      colours = c("#d0d0d0", "#d0d0d0", "#0000ff"),
      values = c(0, color_floor / color_vmax, 1),
      limits = c(0, color_vmax),
      oob = scales::squish,
      name = "log1p(CP10K)",
      breaks = c(0, color_vmax),
      labels = c("0", formatC(max_expr, format = "fg", digits = 2))
    ) +
    coord_cartesian(xlim = x_limits, ylim = y_limits, expand = FALSE) +
    theme_void(base_size = 8) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 9),
      plot.subtitle = element_text(hjust = 0.5, size = 7),
      legend.position = if (show_legend) "right" else "none",
      legend.key.height = grid::unit(16, "pt"),
      legend.key.width = grid::unit(6, "pt"),
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 5)
    ) +
    labs(title = df$display_label[[1]], subtitle = panel_title)
  if (!isTRUE(df$gene_present[[1]])) {
    p <- p + annotate("text", x = mean(x_limits), y = mean(y_limits), label = "missing", size = 4, color = "grey20")
  }
  p
}

panel_plots <- function(layout, cells, marker_df, spec, point_size, color_floor, vmax_quantile) {
  lapply(seq_len(nrow(spec)), function(i) {
    subtitle <- if (spec$panel[[i]] == "custom_marker_panel") {
      ""
    } else if (spec$panel[[i]] == "A_developmental_progression") {
      spec$category[[i]]
    } else {
      paste(spec$lineage_pair[[i]], spec$terminal_identity[[i]], sep = "\n")
    }
    tree_marker_plot(layout, cells, marker_df, spec$gene[[i]], subtitle, point_size, color_floor, vmax_quantile, show_legend = TRUE)
  })
}

safe_gene_filename <- function(gene) {
  gsub("[^A-Za-z0-9]+", "_", gene)
}

print_usage <- function() {
  cat("Usage: Rscript scripts/25_div30_urd_jia_fig_s11_marker_validation.R --tree-rds <tree.rds> --outdir <dir>\n")
  cat("Optional: --genes HES1,NKX2-1,LHX6 --gene-labels Hes1,Nkx2.1,Lhx6 --panel-title <title>\n")
  cat("Optional color controls: --expression-color-floor 1 --vmax-quantile 0.99\n")
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
  table_dir = file.path(opt$outdir, "tables"),
  plot_dir = file.path(opt$outdir, "plots"),
  point_size = as_num(opt$`point-size`, "point-size"),
  expression_color_floor = as_num(opt$`expression-color-floor`, "expression-color-floor"),
  vmax_quantile = as_num(opt$`vmax-quantile`, "vmax-quantile")
)
if (cfg$vmax_quantile <= 0 || cfg$vmax_quantile > 1) stop("vmax-quantile must be in (0, 1]", call. = FALSE)
dir.create(cfg$table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$plot_dir, recursive = TRUE, showWarnings = FALSE)

log_msg("Reading URD tree: ", cfg$tree_rds)
urd <- readRDS(cfg$tree_rds)
counts <- repair_count_dimnames(urd)
layout <- as.data.frame(urd@tree$tree.layout, stringsAsFactors = FALSE)
cells <- as.data.frame(urd@tree$cell.layout, stringsAsFactors = FALSE)
custom_spec <- custom_marker_spec(opt$genes, opt$`gene-labels`)
spec <- if (is.null(custom_spec)) marker_spec() else custom_spec
panel_title <- if (nzchar(opt$`panel-title`)) opt$`panel-title` else "Jia Fig. S11-style URD marker validation"

marker_df <- build_marker_table(counts, cells$cell, spec)
summary_df <- marker_summary(marker_df)
write_tsv(spec, file.path(cfg$table_dir, "jia_fig_s11_marker_order.tsv"))
write_tsv(summary_df, file.path(cfg$table_dir, "jia_fig_s11_marker_expression_summary.tsv"))
write_tsv(marker_df, file.path(cfg$table_dir, "jia_fig_s11_marker_expression_by_cell.tsv.gz"))

missing_genes <- spec$gene[!spec$gene %in% rownames(counts)]
if (length(missing_genes) > 0) {
  log_msg("Missing requested marker(s): ", paste(missing_genes, collapse = ", "))
} else {
  log_msg("All requested Jia Fig. S11 markers are present.")
}

if ("custom_marker_panel" %in% spec$panel) {
  custom_panel <- cowplot::plot_grid(plotlist = panel_plots(layout, cells, marker_df, spec, cfg$point_size, cfg$expression_color_floor, cfg$vmax_quantile), nrow = 1, labels = NULL)
  combined <- cowplot::plot_grid(
    cowplot::ggdraw() + cowplot::draw_label(panel_title, x = 0, hjust = 0, fontface = "bold", size = 11),
    custom_panel,
    ncol = 1,
    rel_heights = c(0.08, 1)
  )
  save_plot_set(
    combined,
    file.path(cfg$plot_dir, "jia_fig_s11_style_urd_marker_validation.png"),
    file.path(cfg$plot_dir, "jia_fig_s11_style_urd_marker_validation.pdf"),
    file.path(cfg$plot_dir, "jia_fig_s11_style_urd_marker_validation.svg"),
    width = 24,
    height = 4.8
  )
} else {
  panel_a_spec <- spec[spec$panel == "A_developmental_progression", , drop = FALSE]
  panel_b_spec <- spec[spec$panel == "B_inhibitory_lineages", , drop = FALSE]
  panel_a <- cowplot::plot_grid(plotlist = panel_plots(layout, cells, marker_df, panel_a_spec, cfg$point_size, cfg$expression_color_floor, cfg$vmax_quantile), nrow = 1, labels = NULL)
  panel_b <- cowplot::plot_grid(plotlist = panel_plots(layout, cells, marker_df, panel_b_spec, cfg$point_size, cfg$expression_color_floor, cfg$vmax_quantile), nrow = 1, labels = NULL)
  panel_a_labeled <- cowplot::plot_grid(
    cowplot::ggdraw() + cowplot::draw_label("A. Developmental progression markers: VZ RGC -> SVZ RGC -> IPC -> newborn neuron", x = 0, hjust = 0, fontface = "bold", size = 11),
    panel_a,
    ncol = 1,
    rel_heights = c(0.08, 1)
  )
  panel_b_labeled <- cowplot::plot_grid(
    cowplot::ggdraw() + cowplot::draw_label("B. Jia inhibitory neuron lineage markers: VZ-RGC -> SVZ-RGC -> Neurogenesis -> lineage branches", x = 0, hjust = 0, fontface = "bold", size = 11),
    panel_b,
    ncol = 1,
    rel_heights = c(0.08, 1)
  )
  combined <- cowplot::plot_grid(panel_a_labeled, panel_b_labeled, ncol = 1, rel_heights = c(1, 1))

  save_plot_set(panel_a_labeled, file.path(cfg$plot_dir, "jia_fig_s11_panel_a_developmental_markers.png"), file.path(cfg$plot_dir, "jia_fig_s11_panel_a_developmental_markers.pdf"), file.path(cfg$plot_dir, "jia_fig_s11_panel_a_developmental_markers.svg"), width = 16, height = 4.2)
  save_plot_set(panel_b_labeled, file.path(cfg$plot_dir, "jia_fig_s11_panel_b_lineage_markers.png"), file.path(cfg$plot_dir, "jia_fig_s11_panel_b_lineage_markers.pdf"), file.path(cfg$plot_dir, "jia_fig_s11_panel_b_lineage_markers.svg"), width = 20, height = 4.2)
  save_plot_set(combined, file.path(cfg$plot_dir, "jia_fig_s11_style_urd_marker_validation.png"), file.path(cfg$plot_dir, "jia_fig_s11_style_urd_marker_validation.pdf"), file.path(cfg$plot_dir, "jia_fig_s11_style_urd_marker_validation.svg"), width = 20, height = 8.8)
}

for (i in seq_len(nrow(spec))) {
  gene <- spec$gene[[i]]
  subtitle <- if (spec$panel[[i]] == "custom_marker_panel") {
    ""
  } else if (spec$panel[[i]] == "A_developmental_progression") {
    spec$category[[i]]
  } else {
    paste(spec$lineage_pair[[i]], spec$terminal_identity[[i]], sep = "\n")
  }
  p_gene <- tree_marker_plot(layout, cells, marker_df, gene, subtitle, cfg$point_size, cfg$expression_color_floor, cfg$vmax_quantile, show_legend = TRUE)
  base <- paste0("jia_fig_s11_marker_tree_overlay_", safe_gene_filename(gene))
  save_plot_set(
    p_gene,
    file.path(cfg$plot_dir, paste0(base, ".png")),
    file.path(cfg$plot_dir, paste0(base, ".pdf")),
    file.path(cfg$plot_dir, paste0(base, ".svg")),
    width = 6.2,
    height = 5.2
  )
}

report <- c(
  "# Jia Fig. S11-Style URD Marker Validation",
  "",
  paste0("- Tree object: `", cfg$tree_rds, "`"),
  "- Expression source: `URD count.data`, transformed at plot time as `log1p(count / total_cell_counts * 10000)`.",
  "- No z-score scaling across genes.",
  "- No branchpoint DE genes.",
  "- No lineage program averages.",
  "- No marker substitutions.",
  paste0("- Expression color map: cross-study marker-expression whiteBlue logic (`#d0d0d0` background floor to `#0000ff`)."),
  paste0("- Colorbar begins at 0; values from 0 to ", cfg$expression_color_floor, " are drawn as background grey."),
  paste0("- Per-gene color maximum: q", cfg$vmax_quantile, " of positive expression, clipped at the maximum color."),
  "",
  "## Marker Panel",
  "",
  paste0("`", paste(spec$display_label, collapse = " | "), "`"),
  "",
  "## Outputs",
  "",
  "- `plots/jia_fig_s11_style_urd_marker_validation.png`",
  "- `plots/jia_fig_s11_style_urd_marker_validation.pdf`",
  "- `plots/jia_fig_s11_style_urd_marker_validation.svg`",
  if ("custom_marker_panel" %in% spec$panel) character() else "- `plots/jia_fig_s11_panel_a_developmental_markers.png`",
  if ("custom_marker_panel" %in% spec$panel) character() else "- `plots/jia_fig_s11_panel_b_lineage_markers.png`",
  "- `plots/jia_fig_s11_marker_tree_overlay_<GENE>.png` and `.pdf` for each individual marker, with visible log1p(CP10K) colorbars",
  "- `tables/jia_fig_s11_marker_order.tsv`",
  "- `tables/jia_fig_s11_marker_expression_summary.tsv`",
  "",
  "## Marker Summary",
  "",
  paste(capture.output(print(summary_df, row.names = FALSE)), collapse = "\n")
)
writeLines(report, file.path(cfg$outdir, "jia_fig_s11_style_marker_validation_report.md"))
log_msg("Done: ", file.path(cfg$outdir, "jia_fig_s11_style_marker_validation_report.md"))
