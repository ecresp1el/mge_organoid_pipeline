#!/usr/bin/env Rscript

# Project candidate DIV90 lineage markers onto a completed URD tree.
#
# This script is specific to the corrected DIV90 v2 smoke logic:
#   - tips are LHX8/ISL1-like, LHX6/NFIA-like, and CRABP1/ANGPT2-like;
#   - clusters 3 and 11 are retained in the manifold but not supplied as tips;
#   - after tree construction, exact marker expression is used to ask whether
#     clusters 3/11 look closer to CRABP1/ANGPT2, LHX6/NFIA, or a separate
#     lineage state.
#
# The script does not rebuild URD, change roots, change tips, or run DE.

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = " ")))
  flush.console()
}

parse_args <- function(args) {
  out <- list(
    `tree-rds` = NULL,
    outdir = NULL,
    `annotation-col` = "cluster_number_name",
    `cluster-col` = "cluster_id_numeric",
    `tip-group-col` = "div90_jia_tip_group",
    `candidate-clusters` = "3,11",
    `marker-genes` = "MEF2C,EPHA5,LHX6,CRABP1,LHX8,NR2F1,NR2F2",
    `point-size` = "0.25",
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

split_csv <- function(x) {
  parts <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  parts[nzchar(parts)]
}

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

segment_membership <- function(object) {
  segs <- object@tree$cells.in.segment
  if (is.null(segs) || length(segs) == 0) {
    return(data.frame(cell = character(), tree_segment = character(), stringsAsFactors = FALSE))
  }
  rows <- lapply(seq_along(segs), function(i) {
    cells <- as.character(segs[[i]])
    if (length(cells) == 0) return(NULL)
    seg_name <- names(segs)[[i]]
    if (is.null(seg_name) || is.na(seg_name) || !nzchar(seg_name)) seg_name <- as.character(i)
    data.frame(cell = cells, tree_segment = seg_name, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

build_marker_table <- function(expr, cells, marker_genes) {
  rows <- lapply(marker_genes, function(gene) {
    present <- gene %in% rownames(expr)
    values <- if (present) as.numeric(expr[gene, cells]) else rep(NA_real_, length(cells))
    data.frame(
      gene = gene,
      cell = cells,
      expression_logupx = values,
      gene_present = present,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

summarize_long <- function(df, group_cols) {
  split_terms <- c(
    lapply(df[, group_cols, drop = FALSE], function(x) {
      x <- as.character(x)
      x[is.na(x) | !nzchar(x)] <- "NA"
      x
    }),
    list(gene = as.character(df$gene))
  )
  groups <- split(df, do.call(interaction, c(split_terms, list(drop = TRUE, sep = "||"))))
  rows <- lapply(groups, function(x) {
    group_values <- x[1, group_cols, drop = FALSE]
    data.frame(
      group_values,
      gene = x$gene[[1]],
      gene_present = if ("gene_present" %in% colnames(x)) x$gene_present[[1]] else TRUE,
      n_cells = length(unique(x$cell)),
      pct_expressed = mean(x$expression_logupx > 0, na.rm = TRUE),
      mean_logupx = mean(x$expression_logupx, na.rm = TRUE),
      median_logupx = median(x$expression_logupx, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

safe_row_mean <- function(mat) {
  if (is.null(mat) || ncol(mat) == 0) return(rep(NA_real_, nrow(mat)))
  out <- rowMeans(mat, na.rm = TRUE)
  out[is.nan(out)] <- NA_real_
  out
}

wide_marker_matrix <- function(marker_df, cells) {
  genes <- unique(marker_df$gene)
  out <- matrix(NA_real_, nrow = length(cells), ncol = length(genes), dimnames = list(cells, genes))
  for (gene in genes) {
    df <- marker_df[marker_df$gene == gene, c("cell", "expression_logupx"), drop = FALSE]
    out[df$cell, gene] <- df$expression_logupx
  }
  out
}

lineage_scores <- function(marker_mat) {
  gene_or_na <- function(gene) {
    if (gene %in% colnames(marker_mat)) marker_mat[, gene] else rep(NA_real_, nrow(marker_mat))
  }
  data.frame(
    cell = rownames(marker_mat),
    LHX8_ISL1_proxy = gene_or_na("LHX8"),
    NR2F1_NR2F2_proxy = safe_row_mean(cbind(gene_or_na("NR2F1"), gene_or_na("NR2F2"))),
    EPHA5_MEF2C_proxy = safe_row_mean(cbind(gene_or_na("EPHA5"), gene_or_na("MEF2C"))),
    LHX6_NFIA_proxy = gene_or_na("LHX6"),
    CRABP1_ANGPT2_proxy = gene_or_na("CRABP1"),
    stringsAsFactors = FALSE
  )
}

tree_marker_plot <- function(layout, cells, marker_df, gene, point_size) {
  df <- marker_df[marker_df$gene == gene, , drop = FALSE]
  plot_df <- merge(cells, df[, c("cell", "expression_logupx", "gene_present")], by = "cell", all.x = TRUE)
  plot_df <- plot_df[order(plot_df$expression_logupx, na.last = TRUE), , drop = FALSE]
  has_segments <- all(c("x1", "y1", "x2", "y2") %in% colnames(layout))
  x_vals <- if (has_segments) c(layout$x1, layout$x2, cells$x) else cells$x
  y_vals <- if (has_segments) c(layout$y1, layout$y2, cells$y) else cells$y
  x_limits <- range(x_vals, na.rm = TRUE)
  y_limits <- range(y_vals, na.rm = TRUE)
  x_pad <- diff(x_limits) * 0.04
  y_pad <- diff(y_limits) * 0.04
  max_expr <- max(plot_df$expression_logupx, na.rm = TRUE)
  if (!is.finite(max_expr) || max_expr <= 0) max_expr <- 1

  p <- ggplot()
  if (has_segments) {
    p <- p + geom_segment(data = layout, aes(x = x1, y = y1, xend = x2, yend = y2), linewidth = 0.25, color = "grey60")
  }
  p <- p +
    geom_point(data = plot_df, aes(x = x, y = y, color = expression_logupx), size = point_size, alpha = 0.86) +
    scale_color_gradient(low = "grey90", high = "#b2182b", limits = c(0, max_expr), na.value = "grey88", name = "logUPX") +
    coord_cartesian(xlim = x_limits + c(-x_pad, x_pad), ylim = y_limits + c(-y_pad, y_pad), expand = FALSE) +
    theme_void(base_size = 8) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 9), legend.position = "none") +
    labs(title = gene)
  if (!isTRUE(df$gene_present[[1]])) {
    p <- p + annotate("text", x = mean(x_limits), y = mean(y_limits), label = "missing", size = 4, color = "grey20")
  }
  p
}

print_usage <- function() {
  cat("Usage: Rscript scripts/27_div90_urd_project_candidate_lineage_markers.R --tree-rds <tree.rds> --outdir <dir>\n")
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
  annotation_col = opt$`annotation-col`,
  cluster_col = opt$`cluster-col`,
  tip_group_col = opt$`tip-group-col`,
  candidate_clusters = split_csv(opt$`candidate-clusters`),
  marker_genes = split_csv(opt$`marker-genes`),
  point_size = as_num(opt$`point-size`, "point-size")
)
dir.create(cfg$table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$plot_dir, recursive = TRUE, showWarnings = FALSE)

log_msg("Reading URD tree: ", cfg$tree_rds)
urd <- readRDS(cfg$tree_rds)
expr <- repair_logupx_dimnames(urd)
meta <- urd@meta
meta$cell <- rownames(meta)

for (col in c(cfg$annotation_col, cfg$cluster_col, cfg$tip_group_col)) {
  if (!col %in% colnames(meta)) stop("URD metadata is missing column: ", col, call. = FALSE)
}

layout <- as.data.frame(urd@tree$tree.layout, stringsAsFactors = FALSE)
cells <- as.data.frame(urd@tree$cell.layout, stringsAsFactors = FALSE)
if (!"cell" %in% colnames(cells)) stop("URD tree cell.layout is missing a `cell` column", call. = FALSE)

tree_meta <- merge(cells, meta, by = "cell", all.x = TRUE)
seg_map <- segment_membership(urd)
tree_meta <- merge(tree_meta, seg_map, by = "cell", all.x = TRUE)
tree_meta$candidate_cluster_status <- ifelse(
  as.character(tree_meta[[cfg$cluster_col]]) %in% cfg$candidate_clusters,
  "candidate_unassigned_cluster_3_11",
  "other_retained_tree_cell"
)

marker_df <- build_marker_table(expr, cells$cell, cfg$marker_genes)
marker_with_meta <- merge(
  marker_df,
  tree_meta[, c("cell", cfg$annotation_col, cfg$cluster_col, cfg$tip_group_col, "tree_segment", "candidate_cluster_status"), drop = FALSE],
  by = "cell",
  all.x = TRUE
)

write_tsv(marker_df, file.path(cfg$table_dir, "div90_candidate_marker_expression_by_cell.tsv.gz"))
write_tsv(marker_with_meta, file.path(cfg$table_dir, "div90_candidate_marker_expression_by_cell_with_metadata.tsv.gz"))
write_tsv(summarize_long(marker_with_meta, c(cfg$annotation_col, cfg$cluster_col)), file.path(cfg$table_dir, "div90_marker_expression_summary_by_cluster.tsv"))
write_tsv(summarize_long(marker_with_meta, c("tree_segment")), file.path(cfg$table_dir, "div90_marker_expression_summary_by_tree_segment.tsv"))

candidate_expr <- marker_with_meta[as.character(marker_with_meta[[cfg$cluster_col]]) %in% cfg$candidate_clusters, , drop = FALSE]
write_tsv(summarize_long(candidate_expr, c(cfg$annotation_col, cfg$cluster_col, "tree_segment")), file.path(cfg$table_dir, "div90_candidate_marker_expression_by_cluster_segment.tsv"))

marker_mat <- wide_marker_matrix(marker_df, cells$cell)
score_df <- merge(lineage_scores(marker_mat), tree_meta[, c("cell", cfg$annotation_col, cfg$cluster_col, cfg$tip_group_col, "tree_segment", "candidate_cluster_status"), drop = FALSE], by = "cell", all.x = TRUE)
write_tsv(score_df, file.path(cfg$table_dir, "div90_candidate_lineage_proxy_scores_by_cell.tsv.gz"))

score_cols <- grep("_proxy$", colnames(score_df), value = TRUE)
score_summary <- summarize_long(
  reshape(
    score_df[, c("cell", cfg$annotation_col, cfg$cluster_col, "tree_segment", score_cols), drop = FALSE],
    varying = score_cols,
    v.names = "expression_logupx",
    timevar = "gene",
    times = score_cols,
    direction = "long"
  ),
  c(cfg$annotation_col, cfg$cluster_col, "tree_segment")
)
colnames(score_summary)[colnames(score_summary) == "gene"] <- "lineage_proxy"
write_tsv(score_summary, file.path(cfg$table_dir, "div90_candidate_lineage_proxy_scores_by_cluster_segment.tsv"))

group_labels <- rep(NA_character_, nrow(tree_meta))
candidate_hit <- as.character(tree_meta[[cfg$cluster_col]]) %in% cfg$candidate_clusters
group_labels[candidate_hit] <- paste0("candidate_cluster_", tree_meta[[cfg$cluster_col]][candidate_hit])
tip_hit <- !is.na(tree_meta[[cfg$tip_group_col]]) & nzchar(as.character(tree_meta[[cfg$tip_group_col]]))
tip_only_hit <- tip_hit & is.na(group_labels)
group_labels[tip_only_hit] <- as.character(tree_meta[[cfg$tip_group_col]][tip_only_hit])
profile_meta <- tree_meta[!is.na(group_labels), , drop = FALSE]
profile_meta$profile_group <- group_labels[!is.na(group_labels)]
profile_mat <- marker_mat[profile_meta$cell, cfg$marker_genes, drop = FALSE]
profile_rows <- lapply(split(seq_len(nrow(profile_meta)), profile_meta$profile_group), function(idx) {
  data.frame(
    profile_group = profile_meta$profile_group[idx][[1]],
    n_cells = length(idx),
    t(colMeans(profile_mat[idx, , drop = FALSE], na.rm = TRUE)),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
})
profile_df <- do.call(rbind, profile_rows)
rownames(profile_df) <- NULL
write_tsv(profile_df, file.path(cfg$table_dir, "div90_candidate_and_tip_marker_profiles.tsv"))

candidate_groups <- paste0("candidate_cluster_", cfg$candidate_clusters)
tip_groups <- setdiff(profile_df$profile_group, candidate_groups)
cor_rows <- list()
for (candidate in intersect(candidate_groups, profile_df$profile_group)) {
  candidate_vec <- as.numeric(profile_df[profile_df$profile_group == candidate, cfg$marker_genes, drop = TRUE])
  for (tip in tip_groups) {
    tip_vec <- as.numeric(profile_df[profile_df$profile_group == tip, cfg$marker_genes, drop = TRUE])
    ok <- is.finite(candidate_vec) & is.finite(tip_vec)
    cor_rows[[length(cor_rows) + 1L]] <- data.frame(
      candidate_cluster = sub("^candidate_cluster_", "", candidate),
      candidate_profile_group = candidate,
      tip_profile_group = tip,
      n_marker_pairs = sum(ok),
      pearson_correlation = if (sum(ok) >= 2) cor(candidate_vec[ok], tip_vec[ok], method = "pearson") else NA_real_,
      stringsAsFactors = FALSE
    )
  }
}
cor_df <- if (length(cor_rows) > 0) do.call(rbind, cor_rows) else data.frame()
write_tsv(cor_df, file.path(cfg$table_dir, "div90_candidate_cluster_to_tip_profile_correlations.tsv"))

assignment_rows <- lapply(cfg$candidate_clusters, function(cluster_id) {
  group <- paste0("candidate_cluster_", cluster_id)
  group_cells <- score_df[as.character(score_df[[cfg$cluster_col]]) == cluster_id, , drop = FALSE]
  proxy_means <- colMeans(group_cells[, score_cols, drop = FALSE], na.rm = TRUE)
  best_proxy <- names(proxy_means)[which.max(proxy_means)]
  cor_sub <- cor_df[cor_df$candidate_cluster == cluster_id, , drop = FALSE]
  if (nrow(cor_sub) > 0 && any(is.finite(cor_sub$pearson_correlation))) {
    best_cor <- cor_sub[which.max(cor_sub$pearson_correlation), , drop = FALSE]
    best_tip <- best_cor$tip_profile_group[[1]]
    best_cor_value <- best_cor$pearson_correlation[[1]]
  } else {
    best_tip <- NA_character_
    best_cor_value <- NA_real_
  }
  data.frame(
    candidate_cluster = cluster_id,
    n_tree_cells = nrow(group_cells),
    best_marker_profile_match = best_tip,
    best_marker_profile_pearson = best_cor_value,
    highest_lineage_proxy = best_proxy,
    highest_lineage_proxy_mean = proxy_means[[best_proxy]],
    interpretation_guardrail = "candidate only; inspect tree position, marker overlays, and profile correlation before assigning as a tip",
    stringsAsFactors = FALSE
  )
})
assignment_df <- do.call(rbind, assignment_rows)
write_tsv(assignment_df, file.path(cfg$table_dir, "div90_candidate_cluster_lineage_assignment.tsv"))

overlay_plots <- lapply(cfg$marker_genes, function(gene) tree_marker_plot(layout, cells, marker_df, gene, cfg$point_size))
overlay_grid <- cowplot::plot_grid(plotlist = overlay_plots, nrow = 1)
save_plot_pair(
  overlay_grid,
  file.path(cfg$plot_dir, "div90_candidate_marker_tree_overlays.png"),
  file.path(cfg$plot_dir, "div90_candidate_marker_tree_overlays.pdf"),
  width = 15,
  height = 3.2
)

profile_long <- reshape(
  profile_df,
  varying = cfg$marker_genes,
  v.names = "mean_logupx",
  timevar = "gene",
  times = cfg$marker_genes,
  direction = "long"
)
profile_long$profile_group <- factor(profile_long$profile_group, levels = rev(unique(profile_df$profile_group)))
profile_long$gene <- factor(profile_long$gene, levels = cfg$marker_genes)
profile_heatmap <- ggplot(profile_long, aes(x = gene, y = profile_group, fill = mean_logupx)) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_fill_gradient(low = "grey95", high = "#b2182b", na.value = "grey80", name = "mean logUPX") +
  theme_minimal(base_size = 9) +
  theme(axis.title = element_blank(), panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "DIV90 candidate clusters 3/11 versus v2 tip marker profiles")
save_plot_pair(
  profile_heatmap,
  file.path(cfg$plot_dir, "div90_candidate_and_tip_marker_profile_heatmap.png"),
  file.path(cfg$plot_dir, "div90_candidate_and_tip_marker_profile_heatmap.pdf"),
  width = 8.5,
  height = 4.5
)

report <- c(
  "# DIV90 Candidate PV Marker Projection",
  "",
  paste0("- Tree object: `", cfg$tree_rds, "`"),
  paste0("- Candidate clusters retained but not used as tips: `", paste(cfg$candidate_clusters, collapse = ", "), "`"),
  paste0("- Marker genes projected: `", paste(cfg$marker_genes, collapse = ", "), "`"),
  "- Expression source: `URD logupx.data`; no z-score scaling.",
  "- This report does not rebuild URD and does not perform branchpoint DE.",
  "",
  "## Assignment Summary",
  "",
  paste(capture.output(print(assignment_df, row.names = FALSE)), collapse = "\n"),
  "",
  "## Outputs",
  "",
  "- `plots/div90_candidate_marker_tree_overlays.png`",
  "- `plots/div90_candidate_and_tip_marker_profile_heatmap.png`",
  "- `tables/div90_candidate_cluster_lineage_assignment.tsv`",
  "- `tables/div90_candidate_cluster_to_tip_profile_correlations.tsv`",
  "- `tables/div90_candidate_and_tip_marker_profiles.tsv`",
  "- `tables/div90_marker_expression_summary_by_cluster.tsv`",
  "- `tables/div90_marker_expression_summary_by_tree_segment.tsv`"
)
writeLines(report, file.path(cfg$outdir, "div90_candidate_marker_projection_report.md"))
log_msg("Done: ", file.path(cfg$outdir, "div90_candidate_marker_projection_report.md"))
