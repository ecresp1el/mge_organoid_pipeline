#!/usr/bin/env Rscript

# Publication-oriented Jia lineage interpretation for a saved DIV30 URD tree.
#
# This script starts from the built URD tree object, not from Seurat clusters.
# It uses the tree segments and branchpoint pseudotimes already encoded by URD:
#   segment 5: root/trunk before the SST vs PV/MGE split
#   segment 1: SST+ cIN branch
#   segment 4: shared PV/MGE branch before the PV vs MGE split
#   segment 2: PV neuron precursor branch
#   segment 3: MGE subpallial neuron branch
#
# Outputs are split into two conceptual layers:
#   1. Decision genes: branchpoint-local DE immediately downstream of splits.
#   2. Jia interpretation: overlap and score maps for Jia lineage programs.

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = " ")))
  flush.console()
}

parse_args <- function(args) {
  out <- list(
    `tree-rds` = NULL,
    outdir = NULL,
    `pseudotime-name` = "",
    `rgc1-col` = "jia_score_RGC1",
    `rgc2-col` = "jia_score_RGC2",
    `ipc-col` = "jia_score_IPC",
    window = "0.06",
    `min-cells` = "40",
    `top-ns` = "50,100,250",
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

required <- c("Matrix", "ggplot2", "cowplot", "gridExtra", "scales")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Missing required R package(s): ", paste(missing, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
  library(cowplot)
  library(gridExtra)
})

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

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

parse_top_ns <- function(x) {
  out <- suppressWarnings(as.integer(strsplit(x, ",", fixed = TRUE)[[1]]))
  out <- sort(unique(out[is.finite(out) & out > 0]))
  if (length(out) == 0) stop("--top-ns must contain at least one positive integer", call. = FALSE)
  out
}

extract_pseudotime <- function(object, requested) {
  pt <- object@pseudotime
  if (!is.data.frame(pt) && !is.matrix(pt)) stop("URD object does not have a tabular pseudotime slot.", call. = FALSE)
  if (ncol(pt) == 0) stop("URD object has no pseudotime columns.", call. = FALSE)
  name <- requested
  if (!nzchar(name)) name <- colnames(pt)[[1]]
  if (!(name %in% colnames(pt))) {
    stop("Pseudotime column not found: ", name, ". Available: ", paste(colnames(pt), collapse = ", "), call. = FALSE)
  }
  values <- as.numeric(pt[rownames(object@meta), name])
  names(values) <- rownames(object@meta)
  list(name = name, values = values)
}

segment_cell_table <- function(object, pt_values) {
  rows <- lapply(names(object@tree$cells.in.segment), function(seg) {
    cells <- intersect(object@tree$cells.in.segment[[seg]], names(pt_values))
    data.frame(
      cell_id = cells,
      segment = seg,
      pseudotime = as.numeric(pt_values[cells]),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

segment_label <- function(seg) {
  labels <- c("5" = "5 root/trunk", "4" = "4 PV/MGE trunk", "1" = "1 SST", "2" = "2 PV", "3" = "3 MGE")
  out <- labels[as.character(seg)]
  out[is.na(out)] <- as.character(seg)[is.na(out)]
  unname(out)
}

select_window_cells <- function(seg_df, segment, split_pt, side, window, min_cells) {
  seg <- seg_df[seg_df$segment == segment & is.finite(seg_df$pseudotime), , drop = FALSE]
  if (nrow(seg) == 0) return(seg)
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
  out <- out[order(out$adjusted_p_value, -abs(out$logFC), out$gene), , drop = FALSE]
  rownames(out) <- NULL
  out
}

rank_branch <- function(de, branch, direction, n) {
  if (direction == "positive") {
    x <- de[de$logFC > 0, , drop = FALSE]
  } else {
    x <- de[de$logFC < 0, , drop = FALSE]
  }
  x$branch <- branch
  x$branch_logFC <- if (direction == "positive") x$logFC else -x$logFC
  x$pct_expressed_branch <- if (direction == "positive") x$pct_expressed_a else x$pct_expressed_b
  x$pct_expressed_other <- if (direction == "positive") x$pct_expressed_b else x$pct_expressed_a
  x <- x[order(x$adjusted_p_value, -x$branch_logFC, x$gene), , drop = FALSE]
  head(x, n)
}

zscore_rows <- function(mat) {
  out <- t(scale(t(as.matrix(mat))))
  out[!is.finite(out)] <- 0
  out
}

mean_expression_by_group <- function(object, genes, groups) {
  genes <- intersect(unique(genes), rownames(object@logupx.data))
  if (length(genes) == 0) return(data.frame())
  rows <- lapply(names(groups), function(group) {
    cells <- intersect(groups[[group]], colnames(object@logupx.data))
    vals <- Matrix::rowMeans(object@logupx.data[genes, cells, drop = FALSE])
    data.frame(group = group, gene = genes, mean_logupx = as.numeric(vals), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

plot_decision_heatmap <- function(object, genes, groups, title) {
  df <- mean_expression_by_group(object, genes, groups)
  if (nrow(df) == 0) return(ggplot() + theme_void() + labs(title = title))
  wide <- reshape(df, idvar = "gene", timevar = "group", direction = "wide")
  rownames(wide) <- wide$gene
  mat <- as.matrix(wide[, setdiff(colnames(wide), "gene"), drop = FALSE])
  colnames(mat) <- sub("^mean_logupx\\.", "", colnames(mat))
  mat <- mat[intersect(genes, rownames(mat)), names(groups), drop = FALSE]
  z <- zscore_rows(mat)
  long <- data.frame(
    gene = rep(rownames(z), times = ncol(z)),
    group = rep(colnames(z), each = nrow(z)),
    z = as.vector(z),
    mean_logupx = as.vector(mat),
    stringsAsFactors = FALSE
  )
  long$gene <- factor(long$gene, levels = rev(rownames(z)))
  long$group <- factor(long$group, levels = names(groups))
  ggplot(long, aes(group, gene, fill = z)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, limits = c(-2, 2), oob = scales::squish) +
    theme_bw(base_size = 8) +
    theme(panel.grid = element_blank(), axis.text.y = element_text(size = 5.8), axis.text.x = element_text(angle = 30, hjust = 1)) +
    labs(title = title, x = NULL, y = NULL, fill = "Row z")
}

jia_segment_summary <- function(object, seg_df, cfg) {
  needed <- c(cfg$rgc1_col, cfg$rgc2_col, cfg$ipc_col)
  missing_cols <- setdiff(needed, colnames(object@meta))
  if (length(missing_cols) > 0) stop("Missing Jia column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  rows <- lapply(split(seg_df, seg_df$segment), function(x) {
    meta <- object@meta[x$cell_id, , drop = FALSE]
    data.frame(
      segment = unique(x$segment),
      segment_label = segment_label(unique(x$segment)),
      n_cells = nrow(x),
      pseudotime_min = min(x$pseudotime, na.rm = TRUE),
      pseudotime_median = median(x$pseudotime, na.rm = TRUE),
      pseudotime_max = max(x$pseudotime, na.rm = TRUE),
      mean_RGC1 = mean(meta[[cfg$rgc1_col]], na.rm = TRUE),
      mean_RGC2 = mean(meta[[cfg$rgc2_col]], na.rm = TRUE),
      mean_IPC = mean(meta[[cfg$ipc_col]], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$segment <- as.character(out$segment)
  out <- out[match(c("5", "4", "1", "2", "3"), out$segment, nomatch = 0), , drop = FALSE]
  rownames(out) <- NULL
  out
}

plot_jia_segment_heatmap <- function(summary_df, title = "Jia progenitor programs by URD segment") {
  long <- rbind(
    data.frame(segment_label = summary_df$segment_label, program = "RGC1", mean_score = summary_df$mean_RGC1, stringsAsFactors = FALSE),
    data.frame(segment_label = summary_df$segment_label, program = "RGC2", mean_score = summary_df$mean_RGC2, stringsAsFactors = FALSE),
    data.frame(segment_label = summary_df$segment_label, program = "IPC", mean_score = summary_df$mean_IPC, stringsAsFactors = FALSE)
  )
  long$segment_label <- factor(long$segment_label, levels = summary_df$segment_label)
  long$program <- factor(long$program, levels = c("RGC1", "RGC2", "IPC"))
  ggplot(long, aes(segment_label, program, fill = mean_score)) +
    geom_tile(color = "white", linewidth = 0.35) +
    geom_text(aes(label = sprintf("%.3f", mean_score)), size = 2.4) +
    scale_fill_gradient(low = "white", high = "#2b8cbe") +
    theme_bw(base_size = 8) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1)) +
    labs(title = title, x = NULL, y = NULL, fill = "Mean score")
}

tree_branchpoint_plot <- function(object, seg_df, split_1, split_2) {
  layout <- as.data.frame(object@tree$tree.layout, stringsAsFactors = FALSE)
  cells <- as.data.frame(object@tree$cell.layout, stringsAsFactors = FALSE)
  cells <- merge(cells, seg_df[, c("cell_id", "segment")], by.x = "cell", by.y = "cell_id", all.x = TRUE)
  cells$segment_label <- segment_label(cells$segment)
  cells$segment_label <- factor(cells$segment_label, levels = segment_label(c("5", "4", "1", "2", "3")))
  branchpoints <- data.frame(
    label = c("BP1\nSST vs PV/MGE", "BP2\nPV vs MGE"),
    x = c(1.25, 0.50),
    y = c(split_1, split_2),
    stringsAsFactors = FALSE
  )
  ggplot() +
    geom_segment(data = layout, aes(x = x1, y = y1, xend = x2, yend = y2), linewidth = 0.55, color = "grey20") +
    geom_point(data = cells, aes(x = x, y = y, color = segment_label), size = 0.33, alpha = 0.75) +
    geom_point(data = branchpoints, aes(x = x, y = y), size = 2.2, shape = 21, fill = "white", color = "black", stroke = 0.6) +
    geom_label(data = branchpoints, aes(x = x, y = y, label = label), size = 2.5, linewidth = 0.2, fill = "white") +
    scale_color_manual(values = c("5 root/trunk" = "#4d4d4d", "4 PV/MGE trunk" = "#756bb1", "1 SST" = "#1b9e77", "2 PV" = "#d95f02", "3 MGE" = "#377eb8"), drop = FALSE) +
    theme_void(base_size = 9) +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    labs(title = "URD tree with branchpoints")
}

jia_lineage_sets <- function() {
  list(
    LHX8_ISL1 = list(
      interpretation = "Early VZ-derived subpallial cholinergic inhibitory neuron lineage",
      core = c("LHX8", "ISL1"),
      additional = character(0),
      tf = c("ZIC4", "HMGA1", "ZIC1", "HEY1", "ID4", "JUN")
    ),
    NR2F1_NR2F2 = list(
      interpretation = "Early VZ-derived subpallial GABAergic inhibitory neuron lineage",
      core = c("NR2F1", "NR2F2"),
      additional = character(0),
      tf = c("SOX3", "E2F1", "FOXP2", "RARB", "EMX2")
    ),
    EPHA5_MEF2C = list(
      interpretation = "SVZ-derived cortex-bound inhibitory neuron lineage",
      core = c("EPHA5", "MEF2C"),
      additional = character(0),
      tf = c("ARX", "SOX11", "HMGA2", "ASCL1", "ELF2", "PBX1", "NRF1", "BCL11A", "PBX3", "DLX5", "DLX2", "DLX1")
    ),
    LHX6_NFIA = list(
      interpretation = "SVZ-derived cortical inhibitory neuron lineage, including chandelier-cell related populations",
      core = c("LHX6", "NFIA"),
      additional = character(0),
      tf = c("ARX", "SOX4", "TCF4", "NFIB", "NFIC", "POU3F2", "POU3F4", "NFAT5", "TAF1", "ZNF713", "RFX4", "PROX1", "SOX5")
    ),
    CRABP1_ANGPT2 = list(
      interpretation = "SVZ-derived subpallial GABAergic inhibitory neuron lineage",
      core = c("CRABP1", "ANGPT2"),
      additional = c("SPOCK1", "DSCAM", "RIPOR2", "RBP1"),
      tf = c("BRCA1", "ETV1", "PURA", "ZEB1", "FOXO3", "KLF12")
    )
  )
}

lineage_set_table <- function(sets, expr_genes, de_genes) {
  rows <- lapply(names(sets), function(name) {
    set <- sets[[name]]
    genes <- unique(c(set$core, set$additional, set$tf))
    data.frame(
      lineage = name,
      category = c(rep("core", length(set$core)), rep("additional", length(set$additional)), rep("tf_program", length(set$tf))),
      gene = c(set$core, set$additional, set$tf),
      in_expression_matrix = c(set$core, set$additional, set$tf) %in% expr_genes,
      in_de_universe = c(set$core, set$additional, set$tf) %in% de_genes,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

branch_enrichment <- function(branch_tables, sets, universe) {
  rows <- list()
  universe <- unique(universe)
  i <- 1L
  for (branch_name in names(branch_tables)) {
    genes <- unique(branch_tables[[branch_name]]$gene)
    genes <- intersect(genes, universe)
    parts <- strsplit(branch_name, "__", fixed = TRUE)[[1]]
    for (lineage in names(sets)) {
      categories <- list(
        core = sets[[lineage]]$core,
        tf_program = sets[[lineage]]$tf,
        all_lineage_genes = unique(c(sets[[lineage]]$core, sets[[lineage]]$additional, sets[[lineage]]$tf))
      )
      for (category in names(categories)) {
        set_genes <- intersect(unique(categories[[category]]), universe)
        overlap <- intersect(genes, set_genes)
        p <- if (length(set_genes) == 0 || length(genes) == 0) NA_real_ else {
          phyper(length(overlap) - 1L, length(set_genes), length(universe) - length(set_genes), length(genes), lower.tail = FALSE)
        }
        rows[[i]] <- data.frame(
          branchpoint = parts[[1]],
          branch = parts[[2]],
          top_n = parts[[3]],
          lineage = lineage,
          category = category,
          universe_size = length(universe),
          branch_gene_count = length(genes),
          lineage_gene_count_in_universe = length(set_genes),
          overlap_count = length(overlap),
          overlap_genes = paste(overlap, collapse = ","),
          p_value = p,
          interpretation = sets[[lineage]]$interpretation,
          stringsAsFactors = FALSE
        )
        i <- i + 1L
      }
    }
  }
  out <- do.call(rbind, rows)
  out$adjusted_p_value <- p.adjust(out$p_value, method = "BH")
  out <- out[order(out$branchpoint, out$branch, as.integer(out$top_n), out$adjusted_p_value, -out$overlap_count), , drop = FALSE]
  rownames(out) <- NULL
  out
}

program_scores <- function(object, sets) {
  rows <- lapply(names(sets), function(name) {
    genes <- intersect(unique(c(sets[[name]]$core, sets[[name]]$additional, sets[[name]]$tf)), rownames(object@logupx.data))
    if (length(genes) == 0) {
      score <- rep(NA_real_, ncol(object@logupx.data))
      names(score) <- colnames(object@logupx.data)
    } else {
      score <- Matrix::colMeans(object@logupx.data[genes, , drop = FALSE])
    }
    data.frame(cell_id = names(score), lineage = name, score = as.numeric(score), n_genes_scored = length(genes), scored_genes = paste(genes, collapse = ","), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

plot_lineage_score_tree <- function(object, scores, lineage) {
  layout <- as.data.frame(object@tree$tree.layout, stringsAsFactors = FALSE)
  cells <- as.data.frame(object@tree$cell.layout, stringsAsFactors = FALSE)
  score <- scores[scores$lineage == lineage, c("cell_id", "score"), drop = FALSE]
  cells <- merge(cells, score, by.x = "cell", by.y = "cell_id", all.x = TRUE)
  ggplot() +
    geom_segment(data = layout, aes(x = x1, y = y1, xend = x2, yend = y2), linewidth = 0.35, color = "grey35") +
    geom_point(data = cells, aes(x = x, y = y, color = score), size = 0.32, alpha = 0.85) +
    scale_color_gradient(low = "grey92", high = "#b2182b", na.value = "grey85") +
    theme_void(base_size = 8) +
    theme(legend.position = "right") +
    labs(title = lineage, color = "Mean\nlogUPX")
}

segment_expression_heatmap <- function(object, genes, seg_df, title) {
  genes <- intersect(unique(genes), rownames(object@logupx.data))
  segment_order <- c("5", "4", "1", "2", "3")
  groups <- lapply(segment_order, function(seg) seg_df$cell_id[seg_df$segment == seg])
  names(groups) <- segment_label(segment_order)
  df <- mean_expression_by_group(object, genes, groups)
  wide <- reshape(df, idvar = "gene", timevar = "group", direction = "wide")
  rownames(wide) <- wide$gene
  mat <- as.matrix(wide[, setdiff(colnames(wide), "gene"), drop = FALSE])
  colnames(mat) <- sub("^mean_logupx\\.", "", colnames(mat))
  mat <- mat[genes, names(groups), drop = FALSE]
  z <- zscore_rows(mat)
  long <- data.frame(
    gene = rep(rownames(z), times = ncol(z)),
    group = rep(colnames(z), each = nrow(z)),
    z = as.vector(z),
    stringsAsFactors = FALSE
  )
  long$gene <- factor(long$gene, levels = rev(rownames(z)))
  long$group <- factor(long$group, levels = names(groups))
  ggplot(long, aes(group, gene, fill = z)) +
    geom_tile(color = "white", linewidth = 0.15) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, limits = c(-2, 2), oob = scales::squish) +
    theme_bw(base_size = 8) +
    theme(panel.grid = element_blank(), axis.text.y = element_text(size = 5.4), axis.text.x = element_text(angle = 35, hjust = 1)) +
    labs(title = title, x = NULL, y = NULL, fill = "Row z")
}

save_plot_pair <- function(plot, png_path, pdf_path, width, height, dpi = 300) {
  ggsave(png_path, plot, width = width, height = height, dpi = dpi, bg = "white")
  ggsave(pdf_path, plot, width = width, height = height, bg = "white")
}

print_usage <- function() {
  cat("Usage: Rscript scripts/23_div30_urd_jia_lineage_branchpoint_figure.R --tree-rds <tree.rds> --outdir <dir> --pseudotime-name <name>\n")
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
  rgc1_col = opt$`rgc1-col`,
  rgc2_col = opt$`rgc2-col`,
  ipc_col = opt$`ipc-col`,
  window = as_num(opt$window, "window"),
  min_cells = as_int(opt$`min-cells`, "min-cells"),
  top_ns = parse_top_ns(opt$`top-ns`),
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

if (cfg$genes == "variable") {
  de_genes <- intersect(urd@var.genes, rownames(urd@logupx.data))
} else {
  de_genes <- rownames(urd@logupx.data)
}
if (length(de_genes) == 0) stop("No genes available for DE.", call. = FALSE)

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
  pseudotime_median = c(median(bp1_up$pseudotime), median(bp1_sst$pseudotime), median(bp1_pvmge$pseudotime), median(bp2_up$pseudotime), median(bp2_pv$pseudotime), median(bp2_mge$pseudotime)),
  selection_width = c(unique(bp1_up$selection_width), unique(bp1_sst$selection_width), unique(bp1_pvmge$selection_width), unique(bp2_up$selection_width), unique(bp2_pv$selection_width), unique(bp2_mge$selection_width)),
  stringsAsFactors = FALSE
)
write_tsv(selection_summary, file.path(cfg$table_dir, "branchpoint_cell_selection_summary.tsv"))

log_msg("Running branchpoint DE over ", length(de_genes), " genes")
bp1_de <- de_between_groups(urd, bp1_sst$cell_id, bp1_pvmge$cell_id, "SST", "PV_MGE", de_genes)
bp2_de <- de_between_groups(urd, bp2_pv$cell_id, bp2_mge$cell_id, "PV", "MGE", de_genes)
write_tsv(bp1_de, file.path(cfg$table_dir, "branchpoint1_all_de_sst_vs_pv_mge.tsv"))
write_tsv(bp2_de, file.path(cfg$table_dir, "branchpoint2_all_de_pv_vs_mge.tsv"))

branch_tables <- list()
for (n in cfg$top_ns) {
  branch_tables[[paste("SST_vs_PV_MGE", "SST", n, sep = "__")]] <- rank_branch(bp1_de, "SST", "positive", n)
  branch_tables[[paste("SST_vs_PV_MGE", "PV_MGE", n, sep = "__")]] <- rank_branch(bp1_de, "PV_MGE", "negative", n)
  branch_tables[[paste("PV_vs_MGE", "PV", n, sep = "__")]] <- rank_branch(bp2_de, "PV", "positive", n)
  branch_tables[[paste("PV_vs_MGE", "MGE", n, sep = "__")]] <- rank_branch(bp2_de, "MGE", "negative", n)
}
for (name in names(branch_tables)) {
  parts <- strsplit(name, "__", fixed = TRUE)[[1]]
  file <- paste0(tolower(parts[[1]]), "_", tolower(parts[[2]]), "_enriched_top", parts[[3]], ".tsv")
  write_tsv(branch_tables[[name]], file.path(cfg$table_dir, file))
}

log_msg("Computing Jia lineage overlaps and scores")
sets <- jia_lineage_sets()
set_tbl <- lineage_set_table(sets, rownames(urd@logupx.data), de_genes)
write_tsv(set_tbl, file.path(cfg$table_dir, "jia_lineage_marker_inventory.tsv"))

overlap <- branch_enrichment(branch_tables, sets, de_genes)
write_tsv(overlap, file.path(cfg$table_dir, "jia_lineage_overlap_enrichment.tsv"))

all_lineage_overlap <- overlap[overlap$category == "all_lineage_genes", , drop = FALSE]
best <- do.call(rbind, lapply(split(all_lineage_overlap, paste(all_lineage_overlap$branchpoint, all_lineage_overlap$branch, all_lineage_overlap$top_n, sep = "__")), function(x) {
  x <- x[order(x$adjusted_p_value, -x$overlap_count), , drop = FALSE]
  head(x, 1)
}))
interpretation <- data.frame(
  branchpoint = best$branchpoint,
  branch = best$branch,
  top_n = best$top_n,
  best_matching_jia_lineage = ifelse(best$overlap_count > 0, best$lineage, "No Jia marker overlap"),
  overlap_count = best$overlap_count,
  supporting_genes = best$overlap_genes,
  adjusted_p_value = ifelse(best$overlap_count > 0, best$adjusted_p_value, NA_real_),
  interpretation = ifelse(best$overlap_count > 0, best$interpretation, "No Jia marker genes from this set were present in this branch-enriched list."),
  stringsAsFactors = FALSE
)
write_tsv(interpretation, file.path(cfg$table_dir, "branchpoint_jia_lineage_interpretation.tsv"))

jia_summary <- jia_segment_summary(urd, seg_df, cfg)
write_tsv(jia_summary, file.path(cfg$table_dir, "jia_rgc_ipc_program_by_tree_segment.tsv"))

scores <- program_scores(urd, sets)
write_tsv(scores, file.path(cfg$table_dir, "jia_lineage_program_scores_by_cell.tsv"))

all_lineage_genes <- unique(unlist(lapply(sets, function(x) c(x$core, x$additional, x$tf)), use.names = FALSE))
all_lineage_genes <- intersect(all_lineage_genes, rownames(urd@logupx.data))
seg_groups <- split(seg_df$cell_id, segment_label(seg_df$segment))
segment_expr <- mean_expression_by_group(urd, all_lineage_genes, seg_groups)
write_tsv(segment_expr, file.path(cfg$table_dir, "jia_lineage_marker_expression_by_segment.tsv"))

bp_groups <- list(
  BP1_upstream = bp1_up$cell_id,
  BP1_SST = bp1_sst$cell_id,
  BP1_PV_MGE = bp1_pvmge$cell_id,
  BP2_upstream = bp2_up$cell_id,
  BP2_PV = bp2_pv$cell_id,
  BP2_MGE = bp2_mge$cell_id
)
branchpoint_expr <- mean_expression_by_group(urd, all_lineage_genes, bp_groups)
write_tsv(branchpoint_expr, file.path(cfg$table_dir, "jia_lineage_marker_expression_by_branchpoint_region.tsv"))

log_msg("Generating publication figures")
panel_a <- tree_branchpoint_plot(urd, seg_df, split_1, split_2)
panel_b <- plot_decision_heatmap(
  urd,
  head(bp1_de$gene, 20),
  list(upstream = bp1_up$cell_id, SST = bp1_sst$cell_id, PV_MGE = bp1_pvmge$cell_id),
  "BP1 decision genes: SST vs PV/MGE"
)
panel_c <- plot_decision_heatmap(
  urd,
  head(bp2_de$gene, 20),
  list(upstream = bp2_up$cell_id, PV = bp2_pv$cell_id, MGE = bp2_mge$cell_id),
  "BP2 decision genes: PV vs MGE"
)
panel_d <- plot_jia_segment_heatmap(jia_summary)
publication <- cowplot::plot_grid(panel_a, panel_b, panel_c, panel_d, labels = c("A", "B", "C", "D"), ncol = 2, align = "hv")
save_plot_pair(publication, file.path(cfg$plot_dir, "urd_branchpoint_jia_publication_panel.png"), file.path(cfg$plot_dir, "urd_branchpoint_jia_publication_panel.pdf"), width = 12, height = 10)
save_plot_pair(panel_a, file.path(cfg$plot_dir, "panel_a_urd_tree_branchpoints.png"), file.path(cfg$plot_dir, "panel_a_urd_tree_branchpoints.pdf"), width = 7, height = 6)
save_plot_pair(panel_b, file.path(cfg$plot_dir, "panel_b_branchpoint1_top20_heatmap.png"), file.path(cfg$plot_dir, "panel_b_branchpoint1_top20_heatmap.pdf"), width = 5.5, height = 6.2)
save_plot_pair(panel_c, file.path(cfg$plot_dir, "panel_c_branchpoint2_top20_heatmap.png"), file.path(cfg$plot_dir, "panel_c_branchpoint2_top20_heatmap.pdf"), width = 5.5, height = 6.2)
save_plot_pair(panel_d, file.path(cfg$plot_dir, "panel_d_jia_rgc_ipc_segment_heatmap.png"), file.path(cfg$plot_dir, "panel_d_jia_rgc_ipc_segment_heatmap.pdf"), width = 6.5, height = 3.4)

score_plots <- lapply(names(sets), function(lineage) plot_lineage_score_tree(urd, scores, lineage))
score_panel <- cowplot::plot_grid(plotlist = score_plots, labels = names(sets), ncol = 3)
save_plot_pair(score_panel, file.path(cfg$plot_dir, "jia_lineage_program_scores_on_tree.png"), file.path(cfg$plot_dir, "jia_lineage_program_scores_on_tree.pdf"), width = 12, height = 8)

marker_segment_plot <- segment_expression_heatmap(urd, all_lineage_genes, seg_df, "Jia lineage marker expression by URD segment")
save_plot_pair(marker_segment_plot, file.path(cfg$plot_dir, "jia_lineage_marker_segment_heatmap.png"), file.path(cfg$plot_dir, "jia_lineage_marker_segment_heatmap.pdf"), width = 8, height = 9)

bp_gene_union <- unique(c(head(bp1_de$gene, 20), head(bp2_de$gene, 20)))
write_tsv(data.frame(panel = c(rep("B_BP1", 20), rep("C_BP2", 20)), gene = c(head(bp1_de$gene, 20), head(bp2_de$gene, 20)), stringsAsFactors = FALSE), file.path(cfg$table_dir, "publication_panel_top20_genes.tsv"))

report <- c(
  "# Jia Lineage Branchpoint Interpretation",
  "",
  paste0("- Tree object: `", cfg$tree_rds, "`"),
  paste0("- Pseudotime: `", cfg$pseudotime_name, "`"),
  paste0("- Branchpoint window: ", cfg$window),
  paste0("- DE universe: ", cfg$genes, " genes (n=", length(de_genes), ")"),
  "",
  "## Core Figure",
  "",
  "- `plots/urd_branchpoint_jia_publication_panel.png`",
  "- `plots/urd_branchpoint_jia_publication_panel.pdf`",
  "",
  "## Branchpoint Selection",
  "",
  paste(capture.output(print(selection_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Jia RGC/RGC2/IPC Segment Localization",
  "",
  paste(capture.output(print(jia_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Best Jia Lineage Matches",
  "",
  paste(capture.output(print(interpretation, row.names = FALSE)), collapse = "\n"),
  "",
  "## Main Tables",
  "",
  "- `tables/branchpoint1_all_de_sst_vs_pv_mge.tsv`",
  "- `tables/branchpoint2_all_de_pv_vs_mge.tsv`",
  "- `tables/*_enriched_top50.tsv`, `*_top100.tsv`, `*_top250.tsv`",
  "- `tables/jia_lineage_overlap_enrichment.tsv`",
  "- `tables/branchpoint_jia_lineage_interpretation.tsv`",
  "- `tables/jia_lineage_marker_expression_by_segment.tsv`",
  "- `tables/jia_lineage_program_scores_by_cell.tsv`"
)
writeLines(report, file.path(cfg$outdir, "jia_lineage_branchpoint_interpretation_report.md"))

log_msg("Done: ", file.path(cfg$outdir, "jia_lineage_branchpoint_interpretation_report.md"))
