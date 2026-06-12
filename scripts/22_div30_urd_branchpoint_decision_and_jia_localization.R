#!/usr/bin/env Rscript

# Branchpoint decision-gene analysis plus Jia progenitor-program localization
# for a saved DIV30 URD tree object.
#
# These are intentionally separate analyses:
#   A. Jia localization asks where RGC1/RGC2/IPC programs sit on URD segments.
#   B. Branchpoint DE asks which genes distinguish daughter branches immediately
#      downstream of branch splits.
#
# Current tree convention from scripts/16/17:
#   tip 1 = SST+ cIN
#   tip 2 = PV neuron precursor
#   tip 3 = MGE subpallial neurons
#   segment 4 = PV/MGE shared trunk
#   segment 5 = root/trunk before SST vs PV/MGE split

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
    `window` = "0.06",
    `min-cells` = "40",
    `top-n` = "50",
    `genes` = "variable",
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

required <- c("Matrix", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Missing required R package(s): ", paste(missing, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
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
    cells <- object@tree$cells.in.segment[[seg]]
    cells <- intersect(cells, names(pt_values))
    data.frame(
      cell_id = cells,
      segment = seg,
      pseudotime = as.numeric(pt_values[cells]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

jia_segment_summary <- function(object, seg_df, cfg) {
  needed <- c(cfg$rgc1_col, cfg$rgc2_col, cfg$ipc_col)
  missing_cols <- setdiff(needed, colnames(object@meta))
  if (length(missing_cols) > 0) stop("Missing Jia column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  rows <- lapply(split(seg_df, seg_df$segment), function(x) {
    meta <- object@meta[x$cell_id, , drop = FALSE]
    data.frame(
      segment = unique(x$segment),
      n_cells = nrow(x),
      pseudotime_min = min(x$pseudotime, na.rm = TRUE),
      pseudotime_median = median(x$pseudotime, na.rm = TRUE),
      pseudotime_max = max(x$pseudotime, na.rm = TRUE),
      mean_RGC1 = mean(meta[[cfg$rgc1_col]], na.rm = TRUE),
      mean_RGC2 = mean(meta[[cfg$rgc2_col]], na.rm = TRUE),
      mean_IPC = mean(meta[[cfg$ipc_col]], na.rm = TRUE),
      median_RGC1 = median(meta[[cfg$rgc1_col]], na.rm = TRUE),
      median_RGC2 = median(meta[[cfg$rgc2_col]], na.rm = TRUE),
      median_IPC = median(meta[[cfg$ipc_col]], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(as.numeric(out$segment)), , drop = FALSE]
}

jia_long <- function(summary_df) {
  rbind(
    data.frame(segment = summary_df$segment, n_cells = summary_df$n_cells, pseudotime_median = summary_df$pseudotime_median, program = "RGC1", mean_score = summary_df$mean_RGC1, stringsAsFactors = FALSE),
    data.frame(segment = summary_df$segment, n_cells = summary_df$n_cells, pseudotime_median = summary_df$pseudotime_median, program = "RGC2", mean_score = summary_df$mean_RGC2, stringsAsFactors = FALSE),
    data.frame(segment = summary_df$segment, n_cells = summary_df$n_cells, pseudotime_median = summary_df$pseudotime_median, program = "IPC", mean_score = summary_df$mean_IPC, stringsAsFactors = FALSE)
  )
}

plot_jia_segment_heatmap <- function(summary_df, path) {
  long <- jia_long(summary_df)
  p <- ggplot(long, aes(segment, program, fill = mean_score)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.3f", mean_score)), size = 3) +
    scale_fill_gradient(low = "white", high = "#2b8cbe") +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank()) +
    labs(title = "Jia progenitor programs by URD segment", x = "URD segment", y = NULL, fill = "Mean score")
  ggsave(path, p, width = 7, height = 4, dpi = 240, bg = "white")
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

de_between_groups <- function(object, cells_a, cells_b, label_a, label_b, genes, top_n) {
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
    pct_expressed = as.numeric(pct_a - pct_b),
    stringsAsFactors = FALSE
  )
  out <- out[order(out$adjusted_p_value, -abs(out$logFC), out$gene), , drop = FALSE]
  head(out, top_n)
}

plot_decision_heatmap <- function(object, genes, groups, path, title) {
  genes <- unique(genes)
  rows <- lapply(names(groups), function(group) {
    cells <- intersect(groups[[group]], colnames(object@logupx.data))
    vals <- Matrix::rowMeans(object@logupx.data[genes, cells, drop = FALSE])
    data.frame(group = group, gene = genes, mean_logupx = as.numeric(vals), stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows)
  df$gene <- factor(df$gene, levels = rev(unique(genes)))
  p <- ggplot(df, aes(group, gene, fill = mean_logupx)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "#b2182b") +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(), axis.text.y = element_text(size = 6)) +
    labs(title = title, x = NULL, y = NULL, fill = "Mean\nlogUPX")
  ggsave(path, p, width = 7, height = max(5, length(genes) * 0.13), dpi = 240, bg = "white")
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(opt$help)) {
  cat("Usage: Rscript scripts/22_div30_urd_branchpoint_decision_and_jia_localization.R --tree-rds <tree.rds> --outdir <dir>\n")
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
  top_n = as_int(opt$`top-n`, "top-n"),
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

log_msg("Computing Jia program localization by segment")
jia_summary <- jia_segment_summary(urd, seg_df, cfg)
write_tsv(jia_summary, file.path(cfg$table_dir, "jia_program_by_tree_segment.tsv"))
plot_jia_segment_heatmap(jia_summary, file.path(cfg$plot_dir, "jia_program_by_tree_segment_heatmap.png"))

genes <- if (cfg$genes == "variable") intersect(urd@var.genes, rownames(urd@logupx.data)) else rownames(urd@logupx.data)
if (length(genes) == 0) stop("No genes available for DE.", call. = FALSE)

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
  pseudotime_min = c(min(bp1_up$pseudotime), min(bp1_sst$pseudotime), min(bp1_pvmge$pseudotime), min(bp2_up$pseudotime), min(bp2_pv$pseudotime), min(bp2_mge$pseudotime)),
  pseudotime_median = c(median(bp1_up$pseudotime), median(bp1_sst$pseudotime), median(bp1_pvmge$pseudotime), median(bp2_up$pseudotime), median(bp2_pv$pseudotime), median(bp2_mge$pseudotime)),
  pseudotime_max = c(max(bp1_up$pseudotime), max(bp1_sst$pseudotime), max(bp1_pvmge$pseudotime), max(bp2_up$pseudotime), max(bp2_pv$pseudotime), max(bp2_mge$pseudotime)),
  selection_width = c(unique(bp1_up$selection_width), unique(bp1_sst$selection_width), unique(bp1_pvmge$selection_width), unique(bp2_up$selection_width), unique(bp2_pv$selection_width), unique(bp2_mge$selection_width)),
  stringsAsFactors = FALSE
)
write_tsv(selection_summary, file.path(cfg$table_dir, "branchpoint_cell_selection_summary.tsv"))
write_tsv(data.frame(branchpoint = "SST_vs_PV_MGE", region = "upstream", cell_id = bp1_up$cell_id), file.path(cfg$table_dir, "branchpoint1_upstream_cells.tsv"))
write_tsv(rbind(
  data.frame(branchpoint = "SST_vs_PV_MGE", branch = "SST", cell_id = bp1_sst$cell_id),
  data.frame(branchpoint = "SST_vs_PV_MGE", branch = "PV_MGE", cell_id = bp1_pvmge$cell_id)
), file.path(cfg$table_dir, "branchpoint1_downstream_cells.tsv"))
write_tsv(data.frame(branchpoint = "PV_vs_MGE", region = "upstream", cell_id = bp2_up$cell_id), file.path(cfg$table_dir, "branchpoint2_upstream_cells.tsv"))
write_tsv(rbind(
  data.frame(branchpoint = "PV_vs_MGE", branch = "PV", cell_id = bp2_pv$cell_id),
  data.frame(branchpoint = "PV_vs_MGE", branch = "MGE", cell_id = bp2_mge$cell_id)
), file.path(cfg$table_dir, "branchpoint2_downstream_cells.tsv"))

log_msg("Running branchpoint DE")
bp1_de <- de_between_groups(urd, bp1_sst$cell_id, bp1_pvmge$cell_id, "SST_downstream", "PV_MGE_downstream", genes, cfg$top_n)
bp2_de <- de_between_groups(urd, bp2_pv$cell_id, bp2_mge$cell_id, "PV_downstream", "MGE_downstream", genes, cfg$top_n)
write_tsv(bp1_de, file.path(cfg$table_dir, "branchpoint1_decision_genes.tsv"))
write_tsv(bp2_de, file.path(cfg$table_dir, "branchpoint2_decision_genes.tsv"))

plot_decision_heatmap(
  urd,
  head(bp1_de$gene, cfg$top_n),
  list(upstream = bp1_up$cell_id, SST = bp1_sst$cell_id, PV_MGE = bp1_pvmge$cell_id),
  file.path(cfg$plot_dir, "branchpoint1_decision_gene_heatmap.png"),
  "Branchpoint 1 decision genes: SST vs PV/MGE"
)
plot_decision_heatmap(
  urd,
  head(bp2_de$gene, cfg$top_n),
  list(upstream = bp2_up$cell_id, PV = bp2_pv$cell_id, MGE = bp2_mge$cell_id),
  file.path(cfg$plot_dir, "branchpoint2_decision_gene_heatmap.png"),
  "Branchpoint 2 decision genes: PV vs MGE"
)

report <- c(
  "# URD Branchpoint Decision Genes And Jia Localization",
  "",
  paste0("- Tree object: `", cfg$tree_rds, "`"),
  paste0("- Pseudotime: `", cfg$pseudotime_name, "`"),
  paste0("- Initial branchpoint window: ", cfg$window),
  paste0("- Minimum cells per selected region: ", cfg$min_cells),
  "",
  "## Jia Program Localization",
  "",
  paste(capture.output(print(jia_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Branchpoint Cell Selection",
  "",
  paste(capture.output(print(selection_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Decision Gene Tables",
  "",
  "- `tables/branchpoint1_decision_genes.tsv`",
  "- `tables/branchpoint2_decision_genes.tsv`"
)
writeLines(report, file.path(cfg$outdir, "branchpoint_decision_and_jia_localization_report.md"))
log_msg("Done: ", file.path(cfg$outdir, "branchpoint_decision_and_jia_localization_report.md"))
