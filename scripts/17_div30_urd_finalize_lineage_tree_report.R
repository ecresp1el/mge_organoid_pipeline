#!/usr/bin/env Rscript

# Finalize tables and lightweight plots from a saved DIV30 URD lineage tree.
#
# This script deliberately does not rerun random walks or `buildTree()`. It is
# for the post-tree reporting layer only: tree status, segment joins, tip
# composition, branch-specific genes, and simple tree-layout PNGs.

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = " ")))
  flush.console()
}

parse_args <- function(args) {
  out <- list(
    `tree-rds` = NULL,
    outdir = NULL,
    `annotation-col` = "paper_cluster_annotation",
    `pseudotime-name` = "paper_radial_glia_root",
    `tip-mapping` = "",
    `dataset-label` = "DIV30",
    `root-label` = "Radial glia",
    `top-n-branch-genes` = "50",
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
    "  Rscript scripts/17_div30_urd_finalize_lineage_tree_report.R --tree-rds <div30_urd_lineage_tree_object.rds> --outdir <tree_dir>",
    "",
    "Outputs:",
    "  tables/tree_status.tsv",
    "  tables/tree_segment_joins.tsv",
    "  tables/tree_tip_composition.tsv",
    "  tables/branch_specific_genes.tsv",
    "  plots/urd_tree_annotation.png",
    "  plots/urd_tree_pseudotime.png",
    "  urd_lineage_tree_report.md",
    sep = "\n"
  ))
}

required <- c("Matrix", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Missing required R package(s): ", paste(missing, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
})

as_int <- function(x, name) {
  value <- suppressWarnings(as.integer(x))
  if (is.na(value)) stop(name, " must be an integer; got ", x, call. = FALSE)
  value
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

read_tip_mapping <- function(path, object) {
  if (nzchar(path) && file.exists(path)) {
    mapping <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
    if (all(c("tip_id", "paper_cluster_annotation") %in% colnames(mapping))) return(mapping)
  }
  tip_ids <- sort(setdiff(unique(as.character(object@meta$paper_tree_tip_id)), NA))
  labels <- vapply(tip_ids, function(id) {
    cells <- rownames(object@meta)[object@meta$paper_tree_tip_id %in% id]
    tab <- sort(table(object@meta[cells, "paper_cluster_annotation"]), decreasing = TRUE)
    names(tab)[[1]]
  }, character(1))
  data.frame(tip_id = tip_ids, paper_cluster_annotation = labels, stringsAsFactors = FALSE)
}

tree_status <- function(object, tip_mapping) {
  joins <- object@tree$segment.joins
  n_joins <- if (is.null(joins)) 0L else nrow(joins)
  segments <- object@tree$segments
  data.frame(
    tree_slot_length = length(object@tree),
    n_requested_tips = nrow(tip_mapping),
    n_segment_joins = n_joins,
    n_segments = if (is.null(segments)) 0L else length(segments),
    has_distinct_branching = n_joins > 0,
    tips = paste(tip_mapping$tip_id, tip_mapping$paper_cluster_annotation, sep = "=", collapse = "; "),
    stringsAsFactors = FALSE
  )
}

segment_join_table <- function(object) {
  joins <- object@tree$segment.joins
  if (is.null(joins) || nrow(joins) == 0) {
    return(data.frame(parent = character(), child = character(), pseudotime = numeric()))
  }
  as.data.frame(joins, stringsAsFactors = FALSE)
}

tip_composition <- function(object, annotation_col, tip_mapping) {
  rows <- lapply(seq_len(nrow(tip_mapping)), function(i) {
    id <- tip_mapping$tip_id[[i]]
    cells <- rownames(object@meta)[object@meta$paper_tree_tip_id %in% id]
    tab <- table(object@meta[cells, annotation_col])
    data.frame(
      tip_id = id,
      expected_annotation = tip_mapping$paper_cluster_annotation[[i]],
      n_tip_cells = length(cells),
      median_pseudotime = median(object@tree$pseudotime[cells], na.rm = TRUE),
      annotation_composition = paste(names(tab), as.integer(tab), sep = "=", collapse = "; "),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

branch_genes <- function(object, tip_mapping, top_n) {
  genes <- intersect(object@var.genes, rownames(object@logupx.data))
  if (length(genes) == 0) genes <- rownames(object@logupx.data)
  rows <- lapply(seq_len(nrow(tip_mapping)), function(i) {
    tip_id <- tip_mapping$tip_id[[i]]
    cells_1 <- intersect(rownames(object@meta)[object@meta$paper_tree_tip_id == tip_id], colnames(object@logupx.data))
    cells_2 <- intersect(
      rownames(object@meta)[!is.na(object@meta$paper_tree_tip_id) & object@meta$paper_tree_tip_id != tip_id],
      colnames(object@logupx.data)
    )
    if (length(cells_1) == 0 || length(cells_2) == 0) return(NULL)
    mat_1 <- object@logupx.data[genes, cells_1, drop = FALSE]
    mat_2 <- object@logupx.data[genes, cells_2, drop = FALSE]
    mean_1 <- Matrix::rowMeans(mat_1)
    mean_2 <- Matrix::rowMeans(mat_2)
    frac_1 <- Matrix::rowSums(mat_1 > 0) / ncol(mat_1)
    frac_2 <- Matrix::rowSums(mat_2 > 0) / ncol(mat_2)
    log2_fc <- log2((2^mean_1 - 1 + 1e-6) / (2^mean_2 - 1 + 1e-6))
    out <- data.frame(
      tip_id = tip_id,
      tip_annotation = tip_mapping$paper_cluster_annotation[[i]],
      gene = genes,
      mean_logupx_tip = as.numeric(mean_1),
      mean_logupx_other_tips = as.numeric(mean_2),
      frac_expressed_tip = as.numeric(frac_1),
      frac_expressed_other_tips = as.numeric(frac_2),
      log2_fold_change_tip_vs_other_tips = as.numeric(log2_fc),
      expression_fraction_delta = as.numeric(frac_1 - frac_2),
      n_tip_cells = length(cells_1),
      n_other_tip_cells = length(cells_2),
      stringsAsFactors = FALSE
    )
    out$branch_specificity_score <- out$log2_fold_change_tip_vs_other_tips * pmax(out$expression_fraction_delta, 0)
    out <- out[is.finite(out$branch_specificity_score), , drop = FALSE]
    out <- out[order(out$branch_specificity_score, out$log2_fold_change_tip_vs_other_tips, decreasing = TRUE), , drop = FALSE]
    head(out, top_n)
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) {
    return(data.frame(status = "not_available", note = "No branch marker genes could be computed.", stringsAsFactors = FALSE))
  }
  rownames(out) <- NULL
  out
}

plot_tree_layout <- function(object, label_values, label_name, path) {
  layout <- as.data.frame(object@tree$tree.layout, stringsAsFactors = FALSE)
  cells <- as.data.frame(object@tree$cell.layout, stringsAsFactors = FALSE)
  cells[[label_name]] <- label_values[cells$cell]
  p <- ggplot() +
    geom_segment(
      data = layout,
      aes(x = x1, y = y1, xend = x2, yend = y2),
      linewidth = 0.55,
      color = "grey20",
      alpha = 0.9
    ) +
    geom_point(data = cells, aes(x = x, y = y, color = .data[[label_name]]), size = 0.45, alpha = 0.8) +
    theme_void(base_size = 11) +
    theme(plot.background = element_rect(fill = "white", color = NA), legend.position = "right") +
    labs(title = paste("URD tree colored by", label_name), color = label_name)
  if (is.numeric(cells[[label_name]])) {
    p <- p + scale_color_viridis_c(na.value = "grey85")
  }
  ggsave(path, p, width = 8, height = 7, dpi = 240, bg = "white")
}

write_report <- function(path, tree_rds, dataset_label, root_label, status, tip_mapping, tip_comp, joins, branch_gene_path) {
  lines <- c(
    paste0("# ", dataset_label, " URD Lineage Tree Report"),
    "",
    paste0("- Tree object: `", tree_rds, "`"),
    paste0("- Root: `", root_label, "`"),
    "- Tip group column: `paper_tree_tip_id`",
    "",
    "## Tip Mapping",
    "",
    paste(capture.output(print(tip_mapping, row.names = FALSE)), collapse = "\n"),
    "",
    "## Tree Status",
    "",
    paste(capture.output(print(status, row.names = FALSE)), collapse = "\n"),
    "",
    "## Tip Composition",
    "",
    paste(capture.output(print(tip_comp, row.names = FALSE)), collapse = "\n"),
    "",
    "## Segment Joins",
    "",
    if (nrow(joins) > 0) paste(capture.output(print(joins, row.names = FALSE)), collapse = "\n") else "No segment joins were produced.",
    "",
    "## Branch Genes",
    "",
    paste0("Branch-specific genes: `", branch_gene_path, "`")
  )
  writeLines(lines, path)
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(opt$help)) {
  print_usage()
  quit(save = "no", status = 0)
}
if (is.null(opt$`tree-rds`) || !nzchar(opt$`tree-rds`)) stop("--tree-rds is required", call. = FALSE)
if (is.null(opt$outdir) || !nzchar(opt$outdir)) stop("--outdir is required", call. = FALSE)

table_dir <- file.path(opt$outdir, "tables")
plot_dir <- file.path(opt$outdir, "plots")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

log_msg("Reading saved URD tree: ", opt$`tree-rds`)
urd <- readRDS(opt$`tree-rds`)
tip_mapping <- read_tip_mapping(opt$`tip-mapping`, urd)

status <- tree_status(urd, tip_mapping)
joins <- segment_join_table(urd)
tip_comp <- tip_composition(urd, opt$`annotation-col`, tip_mapping)
branches <- branch_genes(urd, tip_mapping, as_int(opt$`top-n-branch-genes`, "top-n-branch-genes"))

write_tsv(status, file.path(table_dir, "tree_status.tsv"))
write_tsv(joins, file.path(table_dir, "tree_segment_joins.tsv"))
write_tsv(tip_comp, file.path(table_dir, "tree_tip_composition.tsv"))
branch_gene_path <- file.path(table_dir, "branch_specific_genes.tsv")
write_tsv(branches, branch_gene_path)

annotation <- as.character(urd@meta[, opt$`annotation-col`])
names(annotation) <- rownames(urd@meta)
pt <- as.numeric(urd@pseudotime[, opt$`pseudotime-name`])
names(pt) <- rownames(urd@pseudotime)
plot_tree_layout(urd, annotation, "paper_cluster_annotation", file.path(plot_dir, "urd_tree_annotation.png"))
plot_tree_layout(urd, pt, "pseudotime", file.path(plot_dir, "urd_tree_pseudotime.png"))

write_report(
  file.path(opt$outdir, "urd_lineage_tree_report.md"),
  opt$`tree-rds`,
  opt$`dataset-label`,
  opt$`root-label`,
  status,
  tip_mapping,
  tip_comp,
  joins,
  branch_gene_path
)

log_msg("Done. Final tree report: ", file.path(opt$outdir, "urd_lineage_tree_report.md"))
