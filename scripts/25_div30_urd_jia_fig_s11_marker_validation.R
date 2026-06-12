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
#   URD logupx.data. Some saved URD tree objects may have logupx.data without
#   dimnames after serialization; when dimensions match count.data, this script
#   restores logupx.data row/column names from count.data before extracting
#   exact requested genes.

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = " ")))
  flush.console()
}

parse_args <- function(args) {
  out <- list(
    `tree-rds` = NULL,
    outdir = NULL,
    `point-size` = "0.28",
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

save_plot_pair <- function(plot, png_path, pdf_path, width, height, dpi = 300) {
  ggsave(png_path, plot, width = width, height = height, dpi = dpi, bg = "white")
  ggsave(pdf_path, plot, width = width, height = height, bg = "white")
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

marker_spec <- function() {
  data.frame(
    panel = c(rep("A_developmental_progression", 4), rep("B_inhibitory_lineages", 5)),
    display_order = c(seq_len(4), seq_len(5)),
    category = c("VZ RGC", "SVZ RGC", "IPC", "Newborn neuron", "LHX8 lineage", "NR2F1 lineage", "EPHA5 lineage", "MEF2C lineage", "CRABP1 lineage"),
    gene = c("HES1", "CACNA1E", "DLX2", "DCX", "LHX8", "NR2F1", "EPHA5", "MEF2C", "CRABP1"),
    lineage_pair = c(NA, NA, NA, NA, "LHX8/ISL1", "NR2F1/NR2F2", "EPHA5/MEF2C", "LHX6/NFIA", "CRABP1/ANGPT2"),
    terminal_identity = c(NA, NA, NA, NA, "Subpallial cholinergic interneuron", "Subpallial GABAergic interneuron", "Cortical GABAergic interneuron", "Cortical interneuron lineage", "Subpallial GABAergic interneuron"),
    stringsAsFactors = FALSE
  )
}

build_marker_table <- function(expr, cells, spec) {
  rows <- lapply(seq_len(nrow(spec)), function(i) {
    gene <- spec$gene[[i]]
    present <- gene %in% rownames(expr)
    values <- if (present) as.numeric(expr[gene, cells]) else rep(NA_real_, length(cells))
    data.frame(
      panel = spec$panel[[i]],
      display_order = spec$display_order[[i]],
      category = spec$category[[i]],
      lineage_pair = spec$lineage_pair[[i]],
      terminal_identity = spec$terminal_identity[[i]],
      gene = gene,
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
    data.frame(
      panel = df$panel[[1]],
      display_order = df$display_order[[1]],
      category = df$category[[1]],
      lineage_pair = df$lineage_pair[[1]],
      terminal_identity = df$terminal_identity[[1]],
      gene = df$gene[[1]],
      gene_present = df$gene_present[[1]],
      n_cells = nrow(df),
      pct_expressed = mean(df$expression_logupx > 0, na.rm = TRUE),
      mean_logupx = mean(df$expression_logupx, na.rm = TRUE),
      median_logupx = median(df$expression_logupx, na.rm = TRUE),
      max_logupx = max(df$expression_logupx, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$panel, out$display_order), , drop = FALSE]
}

tree_marker_plot <- function(layout, cells, marker_df, gene, panel_title, point_size, show_legend = FALSE) {
  df <- marker_df[marker_df$gene == gene, , drop = FALSE]
  plot_df <- merge(cells, df[, c("cell", "expression_logupx", "gene_present")], by = "cell", all.x = TRUE)
  plot_df <- plot_df[order(plot_df$expression_logupx, na.last = TRUE), , drop = FALSE]
  x_limits <- range(c(layout$x1, layout$x2, cells$x), na.rm = TRUE)
  y_limits <- range(c(layout$y1, layout$y2, cells$y), na.rm = TRUE)
  x_pad <- diff(x_limits) * 0.04
  y_pad <- diff(y_limits) * 0.04
  x_limits <- x_limits + c(-x_pad, x_pad)
  y_limits <- y_limits + c(-y_pad, y_pad)
  max_expr <- max(plot_df$expression_logupx, na.rm = TRUE)
  if (!is.finite(max_expr) || max_expr <= 0) max_expr <- 1

  p <- ggplot() +
    geom_segment(data = layout, aes(x = x1, y = y1, xend = x2, yend = y2), linewidth = 0.25, color = "grey60") +
    geom_point(data = plot_df, aes(x = x, y = y, color = expression_logupx), size = point_size, alpha = 0.85) +
    scale_color_gradient(low = "grey90", high = "#b2182b", limits = c(0, max_expr), na.value = "grey88", name = "logUPX") +
    coord_cartesian(xlim = x_limits, ylim = y_limits, expand = FALSE) +
    theme_void(base_size = 8) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 9),
      plot.subtitle = element_text(hjust = 0.5, size = 7),
      legend.position = if (show_legend) "right" else "none"
    ) +
    labs(title = gene, subtitle = panel_title)
  if (!isTRUE(df$gene_present[[1]])) {
    p <- p + annotate("text", x = mean(x_limits), y = mean(y_limits), label = "missing", size = 4, color = "grey20")
  }
  p
}

panel_plots <- function(layout, cells, marker_df, spec, point_size) {
  lapply(seq_len(nrow(spec)), function(i) {
    subtitle <- if (spec$panel[[i]] == "A_developmental_progression") {
      spec$category[[i]]
    } else {
      paste(spec$lineage_pair[[i]], spec$terminal_identity[[i]], sep = "\n")
    }
    tree_marker_plot(layout, cells, marker_df, spec$gene[[i]], subtitle, point_size)
  })
}

print_usage <- function() {
  cat("Usage: Rscript scripts/25_div30_urd_jia_fig_s11_marker_validation.R --tree-rds <tree.rds> --outdir <dir>\n")
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
  point_size = as_num(opt$`point-size`, "point-size")
)
dir.create(cfg$table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$plot_dir, recursive = TRUE, showWarnings = FALSE)

log_msg("Reading URD tree: ", cfg$tree_rds)
urd <- readRDS(cfg$tree_rds)
expr <- repair_logupx_dimnames(urd)
layout <- as.data.frame(urd@tree$tree.layout, stringsAsFactors = FALSE)
cells <- as.data.frame(urd@tree$cell.layout, stringsAsFactors = FALSE)
spec <- marker_spec()

marker_df <- build_marker_table(expr, cells$cell, spec)
summary_df <- marker_summary(marker_df)
write_tsv(spec, file.path(cfg$table_dir, "jia_fig_s11_marker_order.tsv"))
write_tsv(summary_df, file.path(cfg$table_dir, "jia_fig_s11_marker_expression_summary.tsv"))
write_tsv(marker_df, file.path(cfg$table_dir, "jia_fig_s11_marker_expression_by_cell.tsv.gz"))

missing_genes <- spec$gene[!spec$gene %in% rownames(expr)]
if (length(missing_genes) > 0) {
  log_msg("Missing requested marker(s): ", paste(missing_genes, collapse = ", "))
} else {
  log_msg("All requested Jia Fig. S11 markers are present.")
}

panel_a_spec <- spec[spec$panel == "A_developmental_progression", , drop = FALSE]
panel_b_spec <- spec[spec$panel == "B_inhibitory_lineages", , drop = FALSE]
panel_a <- cowplot::plot_grid(plotlist = panel_plots(layout, cells, marker_df, panel_a_spec, cfg$point_size), nrow = 1, labels = NULL)
panel_b <- cowplot::plot_grid(plotlist = panel_plots(layout, cells, marker_df, panel_b_spec, cfg$point_size), nrow = 1, labels = NULL)
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

save_plot_pair(panel_a_labeled, file.path(cfg$plot_dir, "jia_fig_s11_panel_a_developmental_markers.png"), file.path(cfg$plot_dir, "jia_fig_s11_panel_a_developmental_markers.pdf"), width = 11, height = 3.2)
save_plot_pair(panel_b_labeled, file.path(cfg$plot_dir, "jia_fig_s11_panel_b_lineage_markers.png"), file.path(cfg$plot_dir, "jia_fig_s11_panel_b_lineage_markers.pdf"), width = 13, height = 3.2)
save_plot_pair(combined, file.path(cfg$plot_dir, "jia_fig_s11_style_urd_marker_validation.png"), file.path(cfg$plot_dir, "jia_fig_s11_style_urd_marker_validation.pdf"), width = 13, height = 6.6)

report <- c(
  "# Jia Fig. S11-Style URD Marker Validation",
  "",
  paste0("- Tree object: `", cfg$tree_rds, "`"),
  "- Expression source: `URD logupx.data`; dimnames repaired from `count.data` when absent.",
  "- No z-score scaling across genes.",
  "- No branchpoint DE genes.",
  "- No lineage program averages.",
  "- No marker substitutions.",
  "",
  "## Panel A",
  "",
  "`HES1 | CACNA1E | DLX2 | DCX`",
  "",
  "Interpretation order: `VZ RGC -> SVZ RGC -> IPC -> newborn neuron`",
  "",
  "## Panel B",
  "",
  "`LHX8 | NR2F1 | EPHA5 | MEF2C | CRABP1`",
  "",
  "Lineage labels:",
  "",
  "- `LHX8/ISL1`: subpallial cholinergic interneuron",
  "- `NR2F1/NR2F2`: subpallial GABAergic interneuron",
  "- `EPHA5/MEF2C`: cortical GABAergic interneuron",
  "- `LHX6/NFIA`: cortical interneuron lineage",
  "- `CRABP1/ANGPT2`: subpallial GABAergic interneuron",
  "",
  "## Outputs",
  "",
  "- `plots/jia_fig_s11_style_urd_marker_validation.png`",
  "- `plots/jia_fig_s11_style_urd_marker_validation.pdf`",
  "- `plots/jia_fig_s11_panel_a_developmental_markers.png`",
  "- `plots/jia_fig_s11_panel_b_lineage_markers.png`",
  "- `tables/jia_fig_s11_marker_order.tsv`",
  "- `tables/jia_fig_s11_marker_expression_summary.tsv`",
  "",
  "## Marker Summary",
  "",
  paste(capture.output(print(summary_df, row.names = FALSE)), collapse = "\n")
)
writeLines(report, file.path(cfg$outdir, "jia_fig_s11_style_marker_validation_report.md"))
log_msg("Done: ", file.path(cfg$outdir, "jia_fig_s11_style_marker_validation_report.md"))
