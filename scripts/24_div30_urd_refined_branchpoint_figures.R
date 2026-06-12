#!/usr/bin/env Rscript

# Refined DIV30 URD figure sets for Jia lineage interpretation.
#
# This script keeps three analyses separate:
#   Figure Set 1: Jia lineage localization by URD segment.
#   Figure Set 2: branchpoint decision genes from branchpoint-local DE.
#   Figure Set 3: validation of branch identity with Jia markers found in DE.
#
# It uses the saved URD tree and URD segment assignments directly. No Seurat
# clusters are used to define roots, tips, branches, or Jia lineage scores.

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = " ")))
  flush.console()
}

parse_args <- function(args) {
  out <- list(
    `tree-rds` = NULL,
    outdir = NULL,
    `pseudotime-name` = "",
    window = "0.06",
    `min-cells` = "40",
    genes = "variable",
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

required <- c("Matrix", "ggplot2", "cowplot", "ggrepel", "scales")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Missing required R package(s): ", paste(missing, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
  library(cowplot)
  library(ggrepel)
})

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

as_num <- function(x, name) {
  value <- suppressWarnings(as.numeric(x))
  if (is.na(value)) stop(name, " must be numeric; got ", x, call. = FALSE)
  value
}

as_int <- function(x, name) {
  value <- suppressWarnings(as.integer(x))
  if (is.na(value)) stop(name, " must be an integer; got ", x, call. = FALSE)
  value
}

save_plot_pair <- function(plot, png_path, pdf_path, width, height, dpi = 300) {
  ggsave(png_path, plot, width = width, height = height, dpi = dpi, bg = "white")
  ggsave(pdf_path, plot, width = width, height = height, bg = "white")
}

extract_pseudotime <- function(object, requested) {
  pt <- object@pseudotime
  if (!is.data.frame(pt) && !is.matrix(pt)) stop("URD object does not have a tabular pseudotime slot.", call. = FALSE)
  name <- requested
  if (!nzchar(name)) name <- colnames(pt)[[1]]
  if (!(name %in% colnames(pt))) stop("Pseudotime column not found: ", name, call. = FALSE)
  values <- as.numeric(pt[rownames(object@meta), name])
  names(values) <- rownames(object@meta)
  list(name = name, values = values)
}

segment_label <- function(seg) {
  labels <- c("5" = "5 root/trunk", "4" = "4 PV/MGE trunk", "1" = "1 SST", "2" = "2 PV", "3" = "3 MGE")
  out <- labels[as.character(seg)]
  out[is.na(out)] <- as.character(seg)[is.na(out)]
  unname(out)
}

segment_cell_table <- function(object, pt_values) {
  rows <- lapply(names(object@tree$cells.in.segment), function(seg) {
    cells <- intersect(object@tree$cells.in.segment[[seg]], names(pt_values))
    data.frame(
      cell_id = cells,
      segment = seg,
      segment_label = segment_label(seg),
      pseudotime = as.numeric(pt_values[cells]),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$segment <- factor(out$segment, levels = c("5", "4", "1", "2", "3"))
  out$segment_label <- factor(out$segment_label, levels = segment_label(c("5", "4", "1", "2", "3")))
  rownames(out) <- NULL
  out
}

select_window_cells <- function(seg_df, segment, split_pt, side, window, min_cells) {
  seg <- seg_df[as.character(seg_df$segment) == segment & is.finite(seg_df$pseudotime), , drop = FALSE]
  max_width <- max(abs(range(seg$pseudotime, na.rm = TRUE) - split_pt), window)
  widths <- unique(c(window, seq(window, max_width, length.out = 10), max_width))
  for (w in widths) {
    if (side == "upstream") {
      keep <- seg$pseudotime <= split_pt & seg$pseudotime >= split_pt - w
    } else {
      keep <- seg$pseudotime >= split_pt & seg$pseudotime <= split_pt + w
    }
    selected <- seg[keep, , drop = FALSE]
    if (nrow(selected) >= min_cells || identical(w, tail(widths, 1))) {
      selected$selection_width <- w
      selected$selection_side <- side
      return(selected)
    }
  }
  seg[FALSE, , drop = FALSE]
}

jia_lineage_sets <- function() {
  list(
    LHX8_ISL1 = c("LHX8", "ISL1", "ZIC4", "HMGA1", "ZIC1", "HEY1", "ID4", "JUN"),
    NR2F1_NR2F2 = c("NR2F1", "NR2F2", "SOX3", "E2F1", "FOXP2", "RARB", "EMX2"),
    EPHA5_MEF2C = c("EPHA5", "MEF2C", "ARX", "SOX11", "HMGA2", "ASCL1", "ELF2", "PBX1", "NRF1", "BCL11A", "PBX3", "DLX5", "DLX2", "DLX1"),
    LHX6_NFIA = c("LHX6", "NFIA", "ARX", "SOX4", "TCF4", "NFIB", "NFIC", "POU3F2", "POU3F4", "NFAT5", "TAF1", "ZNF713", "RFX4", "PROX1", "SOX5"),
    CRABP1_ANGPT2 = c("CRABP1", "ANGPT2", "SPOCK1", "DSCAM", "RIPOR2", "RBP1", "BRCA1", "ETV1", "PURA", "ZEB1", "FOXO3", "KLF12")
  )
}

candidate_markers <- function() {
  list(
    PV = c("LHX6", "ARX", "NFIB", "TCF4"),
    SST = c("EPHA5", "ASCL1", "DLX1", "DLX2", "HMGA2"),
    MGE = c("SPOCK1", "LHX8", "ISL1", "ID4")
  )
}

lineage_scores_by_cell <- function(object, sets) {
  rows <- lapply(names(sets), function(lineage) {
    genes <- intersect(sets[[lineage]], rownames(object@logupx.data))
    if (length(genes) == 0) {
      score <- rep(NA_real_, ncol(object@logupx.data))
      names(score) <- colnames(object@logupx.data)
    } else {
      score <- Matrix::colMeans(object@logupx.data[genes, , drop = FALSE])
    }
    data.frame(
      cell_id = names(score),
      lineage = lineage,
      score = as.numeric(score),
      n_genes_scored = length(genes),
      scored_genes = paste(genes, collapse = ","),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

lineage_segment_summary <- function(scores, seg_df) {
  x <- merge(scores, seg_df[, c("cell_id", "segment", "segment_label")], by = "cell_id")
  rows <- lapply(split(x, paste(x$segment_label, x$lineage, sep = "|||")), function(df) {
    data.frame(
      segment = as.character(df$segment[[1]]),
      segment_label = as.character(df$segment_label[[1]]),
      lineage = df$lineage[[1]],
      n_cells = nrow(df),
      mean_score = mean(df$score, na.rm = TRUE),
      median_score = median(df$score, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$segment_label <- factor(out$segment_label, levels = segment_label(c("5", "4", "1", "2", "3")))
  out$lineage <- factor(out$lineage, levels = names(jia_lineage_sets()))
  out <- out[order(out$segment_label, out$lineage), , drop = FALSE]
  rownames(out) <- NULL
  out
}

add_lineage_z <- function(summary_df) {
  summary_df$z_score <- NA_real_
  for (lineage in levels(summary_df$lineage)) {
    idx <- which(summary_df$lineage == lineage)
    z <- as.numeric(scale(summary_df$mean_score[idx]))
    z[!is.finite(z)] <- 0
    summary_df$z_score[idx] <- z
  }
  summary_df
}

lineage_segment_enrichment <- function(scores, seg_df) {
  x <- merge(scores, seg_df[, c("cell_id", "segment", "segment_label")], by = "cell_id")
  rows <- list()
  i <- 1L
  for (lineage in unique(x$lineage)) {
    for (seg in levels(seg_df$segment_label)) {
      in_seg <- x$score[x$lineage == lineage & x$segment_label == seg]
      out_seg <- x$score[x$lineage == lineage & x$segment_label != seg]
      p <- suppressWarnings(stats::wilcox.test(in_seg, out_seg)$p.value)
      rows[[i]] <- data.frame(
        segment_label = seg,
        lineage = lineage,
        n_segment_cells = length(in_seg),
        mean_segment_score = mean(in_seg, na.rm = TRUE),
        mean_other_score = mean(out_seg, na.rm = TRUE),
        delta_mean = mean(in_seg, na.rm = TRUE) - mean(out_seg, na.rm = TRUE),
        p_value = p,
        stringsAsFactors = FALSE
      )
      i <- i + 1L
    }
  }
  out <- do.call(rbind, rows)
  out$adjusted_p_value <- p.adjust(out$p_value, method = "BH")
  out <- out[order(out$adjusted_p_value, -out$delta_mean), , drop = FALSE]
  rownames(out) <- NULL
  out
}

plot_lineage_heatmap <- function(summary_df, value_col, title, fill_label, diverging = FALSE) {
  p <- ggplot(summary_df, aes(segment_label, lineage, fill = .data[[value_col]])) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = sprintf("%.2f", .data[[value_col]])), size = 2.4) +
    theme_bw(base_size = 8) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1)) +
    labs(title = title, x = NULL, y = NULL, fill = fill_label)
  if (diverging) {
    p + scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, oob = scales::squish)
  } else {
    p + scale_fill_gradient(low = "white", high = "#b2182b")
  }
}

safe_wilcox <- function(x, y) {
  if (length(unique(c(x, y))) < 2) return(NA_real_)
  suppressWarnings(stats::wilcox.test(x, y)$p.value)
}

de_between_groups <- function(object, cells_a, cells_b, label_a, label_b, genes) {
  cells_a <- intersect(cells_a, colnames(object@logupx.data))
  cells_b <- intersect(cells_b, colnames(object@logupx.data))
  mat_a <- object@logupx.data[genes, cells_a, drop = FALSE]
  mat_b <- object@logupx.data[genes, cells_b, drop = FALSE]
  mean_a <- Matrix::rowMeans(mat_a)
  mean_b <- Matrix::rowMeans(mat_b)
  pct_a <- Matrix::rowSums(mat_a > 0) / ncol(mat_a)
  pct_b <- Matrix::rowSums(mat_b > 0) / ncol(mat_b)
  p <- vapply(genes, function(g) safe_wilcox(as.numeric(mat_a[g, ]), as.numeric(mat_b[g, ])), numeric(1))
  out <- data.frame(
    gene = genes,
    comparison = paste(label_a, "vs", label_b),
    group_a = label_a,
    group_b = label_b,
    n_a = length(cells_a),
    n_b = length(cells_b),
    mean_logupx_a = as.numeric(mean_a),
    mean_logupx_b = as.numeric(mean_b),
    logFC = as.numeric(mean_a - mean_b),
    p_value = p,
    adjusted_p_value = p.adjust(p, method = "BH"),
    pct_expressed_a = as.numeric(pct_a),
    pct_expressed_b = as.numeric(pct_b),
    pct_expressed_delta = as.numeric(pct_a - pct_b),
    stringsAsFactors = FALSE
  )
  out$rank_score <- abs(out$logFC) * -log10(pmax(out$adjusted_p_value, .Machine$double.xmin))
  out <- out[order(out$adjusted_p_value, -abs(out$logFC), out$gene), , drop = FALSE]
  rownames(out) <- NULL
  out
}

rank_branch <- function(de, branch, direction, n) {
  x <- if (direction == "positive") de[de$logFC > 0, , drop = FALSE] else de[de$logFC < 0, , drop = FALSE]
  x$branch <- branch
  x$branch_logFC <- if (direction == "positive") x$logFC else -x$logFC
  x$pct_expressed_branch <- if (direction == "positive") x$pct_expressed_a else x$pct_expressed_b
  x$pct_expressed_other <- if (direction == "positive") x$pct_expressed_b else x$pct_expressed_a
  x$branch_rank_score <- x$branch_logFC * -log10(pmax(x$adjusted_p_value, .Machine$double.xmin))
  x <- x[order(x$adjusted_p_value, -x$branch_logFC, x$gene), , drop = FALSE]
  head(x, n)
}

mean_expression_by_group <- function(object, genes, groups) {
  genes <- intersect(unique(genes), rownames(object@logupx.data))
  rows <- lapply(names(groups), function(group) {
    cells <- intersect(groups[[group]], colnames(object@logupx.data))
    vals <- Matrix::rowMeans(object@logupx.data[genes, cells, drop = FALSE])
    data.frame(group = group, gene = genes, mean_logupx = as.numeric(vals), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

zscore_rows <- function(mat) {
  out <- t(scale(t(as.matrix(mat))))
  out[!is.finite(out)] <- 0
  out
}

plot_decision_heatmap <- function(object, genes, groups, title) {
  df <- mean_expression_by_group(object, genes, groups)
  wide <- reshape(df, idvar = "gene", timevar = "group", direction = "wide")
  rownames(wide) <- wide$gene
  mat <- as.matrix(wide[, setdiff(colnames(wide), "gene"), drop = FALSE])
  colnames(mat) <- sub("^mean_logupx\\.", "", colnames(mat))
  mat <- mat[intersect(genes, rownames(mat)), names(groups), drop = FALSE]
  z <- zscore_rows(mat)
  long <- data.frame(gene = rep(rownames(z), times = ncol(z)), group = rep(colnames(z), each = nrow(z)), z = as.vector(z), stringsAsFactors = FALSE)
  long$gene <- factor(long$gene, levels = rev(rownames(z)))
  long$group <- factor(long$group, levels = names(groups))
  ggplot(long, aes(group, gene, fill = z)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, limits = c(-2, 2), oob = scales::squish) +
    theme_bw(base_size = 8) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1), axis.text.y = element_text(size = 5.8)) +
    labs(title = title, x = NULL, y = NULL, fill = "Row z")
}

plot_volcano <- function(de, title, left_label, right_label, n_label = 12) {
  df <- de
  df$neg_log10_fdr <- -log10(pmax(df$adjusted_p_value, .Machine$double.xmin))
  df$direction <- ifelse(df$logFC > 0, right_label, left_label)
  df$highlight <- df$adjusted_p_value < 0.05 & abs(df$logFC) >= 0.5
  lab <- head(df[order(df$adjusted_p_value, -abs(df$logFC)), , drop = FALSE], n_label)
  ggplot(df, aes(logFC, neg_log10_fdr)) +
    geom_point(aes(color = direction, alpha = highlight), size = 0.85) +
    ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 2.4, max.overlaps = Inf, min.segment.length = 0) +
    geom_vline(xintercept = 0, linewidth = 0.25, color = "grey40") +
    scale_color_manual(values = setNames(c("#377eb8", "#b2182b"), c(left_label, right_label))) +
    scale_alpha_manual(values = c(`FALSE` = 0.35, `TRUE` = 0.9), guide = "none") +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(), legend.position = "bottom") +
    labs(title = title, x = "logFC", y = "-log10 adjusted p-value", color = NULL)
}

tree_overlay_plot <- function(object, gene) {
  layout <- as.data.frame(object@tree$tree.layout, stringsAsFactors = FALSE)
  cells <- as.data.frame(object@tree$cell.layout, stringsAsFactors = FALSE)
  if (gene %in% rownames(object@logupx.data)) {
    cells$expr <- as.numeric(object@logupx.data[gene, cells$cell])
  } else {
    cells$expr <- NA_real_
  }
  ggplot() +
    geom_segment(data = layout, aes(x = x1, y = y1, xend = x2, yend = y2), linewidth = 0.3, color = "grey50") +
    geom_point(data = cells, aes(x = x, y = y, color = expr), size = 0.25, alpha = 0.8) +
    scale_color_gradient(low = "grey92", high = "#b2182b", na.value = "grey80") +
    theme_void(base_size = 7) +
    theme(legend.position = "none") +
    labs(title = gene)
}

marker_dotplot_data <- function(object, genes, seg_df) {
  genes <- intersect(unique(genes), rownames(object@logupx.data))
  rows <- list()
  i <- 1L
  for (seg in levels(seg_df$segment_label)) {
    cells <- as.character(seg_df$cell_id[seg_df$segment_label == seg])
    cells <- intersect(cells, colnames(object@logupx.data))
    mat <- object@logupx.data[genes, cells, drop = FALSE]
    for (g in genes) {
      rows[[i]] <- data.frame(
        segment_label = seg,
        gene = g,
        mean_logupx = mean(as.numeric(mat[g, ]), na.rm = TRUE),
        pct_expressed = mean(as.numeric(mat[g, ]) > 0, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      i <- i + 1L
    }
  }
  out <- do.call(rbind, rows)
  out$segment_label <- factor(out$segment_label, levels = levels(seg_df$segment_label))
  out$gene <- factor(out$gene, levels = genes)
  out
}

plot_marker_dotplot <- function(df, title) {
  ggplot(df, aes(segment_label, gene)) +
    geom_point(aes(size = pct_expressed, color = mean_logupx)) +
    scale_size(range = c(0.4, 5), labels = scales::percent_format()) +
    scale_color_gradient(low = "grey90", high = "#b2182b") +
    theme_bw(base_size = 8) +
    theme(panel.grid = element_line(linewidth = 0.15, color = "grey90"), axis.text.x = element_text(angle = 35, hjust = 1)) +
    labs(title = title, x = NULL, y = NULL, size = "% expressed", color = "Mean logUPX")
}

plot_marker_heatmap <- function(df, title) {
  wide <- reshape(df[, c("segment_label", "gene", "mean_logupx")], idvar = "gene", timevar = "segment_label", direction = "wide")
  rownames(wide) <- wide$gene
  mat <- as.matrix(wide[, setdiff(colnames(wide), "gene"), drop = FALSE])
  colnames(mat) <- sub("^mean_logupx\\.", "", colnames(mat))
  z <- zscore_rows(mat)
  long <- data.frame(gene = rep(rownames(z), times = ncol(z)), segment_label = rep(colnames(z), each = nrow(z)), z = as.vector(z), stringsAsFactors = FALSE)
  long$gene <- factor(long$gene, levels = rev(rownames(z)))
  long$segment_label <- factor(long$segment_label, levels = colnames(z))
  ggplot(long, aes(segment_label, gene, fill = z)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, limits = c(-2, 2), oob = scales::squish) +
    theme_bw(base_size = 8) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1)) +
    labs(title = title, x = NULL, y = NULL, fill = "Row z")
}

print_usage <- function() {
  cat("Usage: Rscript scripts/24_div30_urd_refined_branchpoint_figures.R --tree-rds <tree.rds> --outdir <dir> --pseudotime-name <name>\n")
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
  pseudotime_name = opt$`pseudotime-name`,
  window = as_num(opt$window, "window"),
  min_cells = as_int(opt$`min-cells`, "min-cells"),
  genes = opt$genes
)
dir.create(cfg$table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$plot_dir, recursive = TRUE, showWarnings = FALSE)

log_msg("Reading URD tree: ", cfg$tree_rds)
urd <- readRDS(cfg$tree_rds)
pt <- extract_pseudotime(urd, cfg$pseudotime_name)
cfg$pseudotime_name <- pt$name
seg_df <- segment_cell_table(urd, pt$values)
write_tsv(seg_df, file.path(cfg$table_dir, "tree_segment_cell_assignments.tsv"))

de_genes <- if (cfg$genes == "variable") intersect(urd@var.genes, rownames(urd@logupx.data)) else rownames(urd@logupx.data)
joins <- urd@tree$segment.joins
split_1 <- as.numeric(joins$pseudotime[joins$parent == "5" & joins$child == "1"][[1]])
split_2 <- as.numeric(joins$pseudotime[joins$parent == "4" & joins$child == "2"][[1]])

bp1_up <- select_window_cells(seg_df, "5", split_1, "upstream", cfg$window, cfg$min_cells)
bp1_sst <- select_window_cells(seg_df, "1", split_1, "downstream", cfg$window, cfg$min_cells)
bp1_pvmge <- select_window_cells(seg_df, "4", split_1, "downstream", cfg$window, cfg$min_cells)
bp2_up <- select_window_cells(seg_df, "4", split_2, "upstream", cfg$window, cfg$min_cells)
bp2_pv <- select_window_cells(seg_df, "2", split_2, "downstream", cfg$window, cfg$min_cells)
bp2_mge <- select_window_cells(seg_df, "3", split_2, "downstream", cfg$window, cfg$min_cells)

selection_summary <- data.frame(
  branchpoint = c("SST_vs_PV_MGE", "SST_vs_PV_MGE", "SST_vs_PV_MGE", "PV_vs_MGE", "PV_vs_MGE", "PV_vs_MGE"),
  region = c("upstream_segment5", "downstream_SST_segment1", "downstream_PV_MGE_segment4", "upstream_segment4", "downstream_PV_segment2", "downstream_MGE_segment3"),
  split_pseudotime = c(split_1, split_1, split_1, split_2, split_2, split_2),
  n_cells = c(nrow(bp1_up), nrow(bp1_sst), nrow(bp1_pvmge), nrow(bp2_up), nrow(bp2_pv), nrow(bp2_mge)),
  median_pseudotime = c(median(bp1_up$pseudotime), median(bp1_sst$pseudotime), median(bp1_pvmge$pseudotime), median(bp2_up$pseudotime), median(bp2_pv$pseudotime), median(bp2_mge$pseudotime)),
  selection_width = c(unique(bp1_up$selection_width), unique(bp1_sst$selection_width), unique(bp1_pvmge$selection_width), unique(bp2_up$selection_width), unique(bp2_pv$selection_width), unique(bp2_mge$selection_width)),
  stringsAsFactors = FALSE
)
write_tsv(selection_summary, file.path(cfg$table_dir, "branchpoint_window_summary.tsv"))

log_msg("Figure Set 1: Jia lineage localization by segment")
sets <- jia_lineage_sets()
scores <- lineage_scores_by_cell(urd, sets)
write_tsv(scores, file.path(cfg$table_dir, "jia_lineage_scores_by_cell.tsv"))
lineage_summary <- add_lineage_z(lineage_segment_summary(scores, seg_df))
write_tsv(lineage_summary, file.path(cfg$table_dir, "jia_lineage_scores_by_segment.tsv"))
lineage_enrichment <- lineage_segment_enrichment(scores, seg_df)
write_tsv(lineage_enrichment, file.path(cfg$table_dir, "jia_lineage_segment_enrichment_summary.tsv"))
score_heatmap <- plot_lineage_heatmap(lineage_summary, "mean_score", "Jia lineage scores by URD segment", "Mean score")
z_heatmap <- plot_lineage_heatmap(lineage_summary, "z_score", "Jia lineage z-scores by URD segment", "Segment z", TRUE)
enrich_plot <- ggplot(lineage_enrichment, aes(segment_label, lineage, fill = -log10(pmax(adjusted_p_value, .Machine$double.xmin)))) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = sprintf("d=%.2f", delta_mean)), size = 2.2) +
  scale_fill_gradient(low = "white", high = "#b2182b") +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(title = "Segment lineage enrichment", x = NULL, y = NULL, fill = "-log10 FDR")
save_plot_pair(score_heatmap, file.path(cfg$plot_dir, "figset1_segment_lineage_score_heatmap.png"), file.path(cfg$plot_dir, "figset1_segment_lineage_score_heatmap.pdf"), 7.4, 4.5)
save_plot_pair(z_heatmap, file.path(cfg$plot_dir, "figset1_segment_lineage_zscore_heatmap.png"), file.path(cfg$plot_dir, "figset1_segment_lineage_zscore_heatmap.pdf"), 7.4, 4.5)
save_plot_pair(enrich_plot, file.path(cfg$plot_dir, "figset1_segment_lineage_enrichment_summary.png"), file.path(cfg$plot_dir, "figset1_segment_lineage_enrichment_summary.pdf"), 7.4, 4.5)

log_msg("Figure Set 2: Branchpoint decision genes")
bp1_de <- de_between_groups(urd, bp1_sst$cell_id, bp1_pvmge$cell_id, "SST", "PV_MGE", de_genes)
bp2_de <- de_between_groups(urd, bp2_pv$cell_id, bp2_mge$cell_id, "PV", "MGE", de_genes)
write_tsv(bp1_de, file.path(cfg$table_dir, "branchpoint1_sst_vs_pv_mge_ranked_marker_table.tsv"))
write_tsv(bp2_de, file.path(cfg$table_dir, "branchpoint2_pv_vs_mge_ranked_marker_table.tsv"))
for (n in c(20, 50)) {
  write_tsv(head(bp1_de, n), file.path(cfg$table_dir, paste0("branchpoint1_top", n, "_decision_genes.tsv")))
  write_tsv(head(bp2_de, n), file.path(cfg$table_dir, paste0("branchpoint2_top", n, "_decision_genes.tsv")))
  write_tsv(rank_branch(bp1_de, "SST", "positive", n), file.path(cfg$table_dir, paste0("branchpoint1_sst_enriched_top", n, ".tsv")))
  write_tsv(rank_branch(bp1_de, "PV_MGE", "negative", n), file.path(cfg$table_dir, paste0("branchpoint1_pv_mge_enriched_top", n, ".tsv")))
  write_tsv(rank_branch(bp2_de, "PV", "positive", n), file.path(cfg$table_dir, paste0("branchpoint2_pv_enriched_top", n, ".tsv")))
  write_tsv(rank_branch(bp2_de, "MGE", "negative", n), file.path(cfg$table_dir, paste0("branchpoint2_mge_enriched_top", n, ".tsv")))
}
bp1_volcano <- plot_volcano(bp1_de, "Branchpoint 1: SST vs PV/MGE", "PV/MGE higher", "SST higher")
bp2_volcano <- plot_volcano(bp2_de, "Branchpoint 2: PV vs MGE", "MGE higher", "PV higher")
bp1_heat20 <- plot_decision_heatmap(urd, head(bp1_de$gene, 20), list(upstream = bp1_up$cell_id, SST = bp1_sst$cell_id, PV_MGE = bp1_pvmge$cell_id), "BP1 top 20 decision genes")
bp1_heat50 <- plot_decision_heatmap(urd, head(bp1_de$gene, 50), list(upstream = bp1_up$cell_id, SST = bp1_sst$cell_id, PV_MGE = bp1_pvmge$cell_id), "BP1 top 50 decision genes")
bp2_heat20 <- plot_decision_heatmap(urd, head(bp2_de$gene, 20), list(upstream = bp2_up$cell_id, PV = bp2_pv$cell_id, MGE = bp2_mge$cell_id), "BP2 top 20 decision genes")
bp2_heat50 <- plot_decision_heatmap(urd, head(bp2_de$gene, 50), list(upstream = bp2_up$cell_id, PV = bp2_pv$cell_id, MGE = bp2_mge$cell_id), "BP2 top 50 decision genes")
save_plot_pair(bp1_volcano, file.path(cfg$plot_dir, "figset2_branchpoint1_volcano.png"), file.path(cfg$plot_dir, "figset2_branchpoint1_volcano.pdf"), 6.5, 5.3)
save_plot_pair(bp2_volcano, file.path(cfg$plot_dir, "figset2_branchpoint2_volcano.png"), file.path(cfg$plot_dir, "figset2_branchpoint2_volcano.pdf"), 6.5, 5.3)
save_plot_pair(bp1_heat20, file.path(cfg$plot_dir, "figset2_branchpoint1_top20_heatmap.png"), file.path(cfg$plot_dir, "figset2_branchpoint1_top20_heatmap.pdf"), 5.6, 6.0)
save_plot_pair(bp1_heat50, file.path(cfg$plot_dir, "figset2_branchpoint1_top50_heatmap.png"), file.path(cfg$plot_dir, "figset2_branchpoint1_top50_heatmap.pdf"), 6.0, 10.0)
save_plot_pair(bp2_heat20, file.path(cfg$plot_dir, "figset2_branchpoint2_top20_heatmap.png"), file.path(cfg$plot_dir, "figset2_branchpoint2_top20_heatmap.pdf"), 5.6, 6.0)
save_plot_pair(bp2_heat50, file.path(cfg$plot_dir, "figset2_branchpoint2_top50_heatmap.png"), file.path(cfg$plot_dir, "figset2_branchpoint2_top50_heatmap.pdf"), 6.0, 10.0)

log_msg("Figure Set 3: Jia marker validation")
markers <- candidate_markers()
marker_genes <- unique(unlist(markers, use.names = FALSE))
marker_inventory <- data.frame(
  marker_group = rep(names(markers), lengths(markers)),
  gene = marker_genes,
  in_expression_matrix = marker_genes %in% rownames(urd@logupx.data),
  in_de_universe = marker_genes %in% de_genes,
  in_bp1_de = marker_genes %in% bp1_de$gene,
  in_bp2_de = marker_genes %in% bp2_de$gene,
  stringsAsFactors = FALSE
)
write_tsv(marker_inventory, file.path(cfg$table_dir, "jia_marker_validation_inventory.tsv"))
marker_df <- marker_dotplot_data(urd, marker_genes, seg_df)
write_tsv(marker_df, file.path(cfg$table_dir, "jia_marker_validation_by_segment.tsv"))
dotplot <- plot_marker_dotplot(marker_df, "Jia marker validation by URD segment")
marker_heatmap <- plot_marker_heatmap(marker_df, "Jia marker validation heatmap by URD segment")
save_plot_pair(dotplot, file.path(cfg$plot_dir, "figset3_jia_marker_dotplot_by_segment.png"), file.path(cfg$plot_dir, "figset3_jia_marker_dotplot_by_segment.pdf"), 7.8, 5.4)
save_plot_pair(marker_heatmap, file.path(cfg$plot_dir, "figset3_jia_marker_heatmap_by_segment.png"), file.path(cfg$plot_dir, "figset3_jia_marker_heatmap_by_segment.pdf"), 7.8, 5.4)

overlay_genes <- intersect(marker_genes, rownames(urd@logupx.data))
overlay_plots <- lapply(overlay_genes, function(g) tree_overlay_plot(urd, g))
overlay_panel <- cowplot::plot_grid(plotlist = overlay_plots, ncol = 4)
save_plot_pair(overlay_panel, file.path(cfg$plot_dir, "figset3_jia_marker_tree_overlays.png"), file.path(cfg$plot_dir, "figset3_jia_marker_tree_overlays.pdf"), 10.5, 8.5)

figure_set_panel <- cowplot::plot_grid(
  z_heatmap,
  cowplot::plot_grid(bp1_volcano, bp2_volcano, ncol = 2),
  marker_heatmap,
  ncol = 1,
  labels = c("A", "B", "C"),
  rel_heights = c(1, 1.3, 1)
)
save_plot_pair(figure_set_panel, file.path(cfg$plot_dir, "refined_urd_branchpoint_summary_panel.png"), file.path(cfg$plot_dir, "refined_urd_branchpoint_summary_panel.pdf"), 12, 13)

top_lineage_hits <- lineage_enrichment[lineage_enrichment$delta_mean > 0, , drop = FALSE]
top_lineage_hits <- top_lineage_hits[order(top_lineage_hits$adjusted_p_value, -top_lineage_hits$delta_mean), , drop = FALSE]

report <- c(
  "# Refined URD Branchpoint Figure Sets",
  "",
  paste0("- Tree object: `", cfg$tree_rds, "`"),
  paste0("- Pseudotime: `", cfg$pseudotime_name, "`"),
  paste0("- Branchpoint window: ", cfg$window),
  paste0("- DE universe: ", cfg$genes, " genes (n=", length(de_genes), ")"),
  "",
  "## Figure Set 1: Jia Lineage Localization",
  "",
  "- `plots/figset1_segment_lineage_score_heatmap.png`",
  "- `plots/figset1_segment_lineage_zscore_heatmap.png`",
  "- `plots/figset1_segment_lineage_enrichment_summary.png`",
  "- `tables/jia_lineage_segment_enrichment_summary.tsv`",
  "",
  "Top positive segment-lineage enrichments:",
  "",
  paste(capture.output(print(head(top_lineage_hits, 12), row.names = FALSE)), collapse = "\n"),
  "",
  "## Figure Set 2: Branchpoint Decision Genes",
  "",
  "- `plots/figset2_branchpoint1_volcano.png`",
  "- `plots/figset2_branchpoint2_volcano.png`",
  "- `plots/figset2_branchpoint1_top20_heatmap.png`",
  "- `plots/figset2_branchpoint1_top50_heatmap.png`",
  "- `plots/figset2_branchpoint2_top20_heatmap.png`",
  "- `plots/figset2_branchpoint2_top50_heatmap.png`",
  "- `tables/branchpoint1_sst_vs_pv_mge_ranked_marker_table.tsv`",
  "- `tables/branchpoint2_pv_vs_mge_ranked_marker_table.tsv`",
  "",
  "## Figure Set 3: Jia Marker Validation",
  "",
  "- `plots/figset3_jia_marker_tree_overlays.png`",
  "- `plots/figset3_jia_marker_dotplot_by_segment.png`",
  "- `plots/figset3_jia_marker_heatmap_by_segment.png`",
  "- `tables/jia_marker_validation_by_segment.tsv`",
  "- `tables/jia_marker_validation_inventory.tsv`",
  "",
  "## Summary Panel",
  "",
  "- `plots/refined_urd_branchpoint_summary_panel.png`",
  "- `plots/refined_urd_branchpoint_summary_panel.pdf`"
)
writeLines(report, file.path(cfg$outdir, "refined_urd_branchpoint_figure_sets_report.md"))
log_msg("Done: ", file.path(cfg$outdir, "refined_urd_branchpoint_figure_sets_report.md"))
