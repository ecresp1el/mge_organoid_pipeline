#!/usr/bin/env Rscript

# Generate a lineage decision-tree report from a DIV30 URD object.
#
# This report is intentionally post hoc: it does not recompute URD. It reads an
# existing `div30_first_urd_object.rds`, extracts the root/pseudotime metadata,
# ranks annotations by pseudotime, computes gene-pseudotime correlations, and
# records whether a true URD branch tree exists. If no branch tree exists, the
# decision-tree figure is a linear pseudotime ordering rather than a branch
# decision tree.

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = " ")))
  flush.console()
}

parse_args <- function(args) {
  out <- list(
    `urd-rds` = NULL,
    outdir = NULL,
    `annotation-col` = "paper_cluster_annotation",
    `pseudotime-name` = "",
    `root-col` = "urd_root_candidate",
    `top-n` = "50",
    `correlation-genes` = "variable",
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
    "  Rscript scripts/15_div30_urd_lineage_decision_report.R --urd-rds <div30_first_urd_object.rds> --outdir <report_dir>",
    "",
    "Outputs:",
    "  tables/root_annotation_composition.tsv",
    "  tables/pseudotime_ordering_by_annotation.tsv",
    "  tables/top_negative_pseudotime_genes.tsv",
    "  tables/top_positive_pseudotime_genes.tsv",
    "  tables/branch_structure_status.tsv",
    "  tables/decision_genes_between_branches.tsv",
    "  tables/flood_stability_summary.tsv",
    "  tables/gene_cascade_heatmap_matrix.tsv",
    "  plots/umap_pseudotime.png",
    "  plots/diffusion_map_pseudotime.png",
    "  plots/diffusion_map_annotation.png",
    "  plots/flood_stability.png",
    "  plots/tree_visualization.png",
    "  plots/gene_cascade_heatmap.png",
    "  plots/lineage_decision_tree.png",
    "  lineage_decision_tree_report.md",
    sep = "\n"
  ))
}

required <- c("Matrix", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing required R package(s): ", paste(missing, collapse = ", "), call. = FALSE)
}
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

extract_pseudotime <- function(object, pseudotime_name) {
  pt <- object@pseudotime
  if ((is.data.frame(pt) || is.matrix(pt)) && ncol(pt) > 0) {
    name <- pseudotime_name
    if (!nzchar(name)) name <- colnames(pt)[[1]]
    if (!(name %in% colnames(pt))) {
      stop("Pseudotime name not found: ", name, ". Available: ", paste(colnames(pt), collapse = ", "), call. = FALSE)
    }
    return(list(name = name, values = as.numeric(pt[rownames(object@meta), name])))
  }
  if (is.list(pt) && length(pt) > 0) {
    name <- pseudotime_name
    if (!nzchar(name)) name <- names(pt)[[1]]
    if (!(name %in% names(pt))) {
      stop("Pseudotime name not found: ", name, ". Available: ", paste(names(pt), collapse = ", "), call. = FALSE)
    }
    return(list(name = name, values = as.numeric(pt[[name]][rownames(object@meta)])))
  }
  stop("URD object does not contain extractable pseudotime.", call. = FALSE)
}

root_annotation_composition <- function(meta, annotation_col, root_col) {
  if (!(annotation_col %in% colnames(meta))) stop("Missing annotation column: ", annotation_col, call. = FALSE)
  if (!(root_col %in% colnames(meta))) stop("Missing root column: ", root_col, call. = FALSE)
  root <- as.logical(meta[[root_col]])
  root[is.na(root)] <- FALSE
  tab <- as.data.frame(table(meta[[annotation_col]][root]), stringsAsFactors = FALSE)
  colnames(tab) <- c(annotation_col, "n_root_cells")
  tab <- tab[tab$n_root_cells > 0, , drop = FALSE]
  tab$fraction_root_cells <- tab$n_root_cells / sum(tab$n_root_cells)
  tab[order(tab$n_root_cells, decreasing = TRUE), , drop = FALSE]
}

pseudotime_ordering <- function(meta, pt_values, annotation_col) {
  df <- data.frame(annotation = meta[[annotation_col]], pseudotime = pt_values, stringsAsFactors = FALSE)
  df <- df[is.finite(df$pseudotime) & !is.na(df$annotation) & nzchar(df$annotation), , drop = FALSE]
  split_pt <- split(df$pseudotime, df$annotation)
  out <- data.frame(
    annotation = names(split_pt),
    n_cells = vapply(split_pt, length, integer(1)),
    median_pseudotime = vapply(split_pt, median, numeric(1), na.rm = TRUE),
    mean_pseudotime = vapply(split_pt, mean, numeric(1), na.rm = TRUE),
    q25_pseudotime = vapply(split_pt, quantile, numeric(1), probs = 0.25, na.rm = TRUE),
    q75_pseudotime = vapply(split_pt, quantile, numeric(1), probs = 0.75, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  out$rank_by_median_pseudotime <- rank(out$median_pseudotime, ties.method = "first")
  out[order(out$median_pseudotime), , drop = FALSE]
}

gene_correlations <- function(object, pt_values, mode = "variable") {
  genes <- rownames(object@logupx.data)
  if (mode == "variable") {
    genes <- intersect(object@var.genes, genes)
  } else if (mode != "all") {
    stop("--correlation-genes must be 'variable' or 'all'", call. = FALSE)
  }
  if (length(genes) == 0) stop("No genes selected for correlation.", call. = FALSE)

  mat <- object@logupx.data[genes, rownames(object@meta), drop = FALSE]
  ok_cells <- is.finite(pt_values)
  mat <- mat[, ok_cells, drop = FALSE]
  pt <- pt_values[ok_cells]

  # Spearman correlation is Pearson correlation of ranks. This loop is clearer
  # and memory-friendlier than densifying all genes for larger follow-up runs.
  pt_rank <- rank(pt, ties.method = "average")
  cors <- vapply(
    seq_along(genes),
    function(i) {
      x <- as.numeric(mat[i, ])
      if (stats::sd(x) == 0) return(NA_real_)
      suppressWarnings(stats::cor(rank(x, ties.method = "average"), pt_rank, method = "pearson"))
    },
    numeric(1)
  )
  data.frame(
    gene = genes,
    spearman_with_pseudotime = cors,
    mean_logupx = Matrix::rowMeans(mat),
    pct_cells_expressed = Matrix::rowSums(mat > 0) / ncol(mat),
    stringsAsFactors = FALSE
  )
}

branch_status <- function(object) {
  has_tree_slot <- "tree" %in% slotNames(object)
  tree_len <- if (has_tree_slot) length(object@tree) else 0L
  meta_cols <- colnames(object@meta)
  candidate_cols <- meta_cols[grepl("branch|segment|tip|tree", meta_cols, ignore.case = TRUE)]
  data.frame(
    has_tree_slot = has_tree_slot,
    tree_length = tree_len,
    has_branch_structure = has_tree_slot && tree_len > 0,
    branch_like_metadata_columns = if (length(candidate_cols)) paste(candidate_cols, collapse = ",") else "",
    note = if (has_tree_slot && tree_len > 0) {
      "URD tree slot is populated; branch decision-gene extraction can be added for this object."
    } else {
      "No populated URD tree slot was found. This object contains flood pseudotime, not a reconstructed branch tree."
    },
    stringsAsFactors = FALSE
  )
}

make_decision_figure <- function(ordering, neg_genes, pos_genes, status, plot_path) {
  ordering <- ordering[order(ordering$median_pseudotime), , drop = FALSE]
  ordering$x <- ordering$median_pseudotime
  ordering$y <- 0
  ordering$label <- paste0(ordering$annotation, "\nmedian=", sprintf("%.3f", ordering$median_pseudotime))
  ordering$label_y <- 0.12 + ((seq_len(nrow(ordering)) - 1L) %% 3L) * 0.095

  edges <- data.frame(
    x = head(ordering$x, -1),
    xend = tail(ordering$x, -1),
    y = 0,
    yend = 0
  )
  root_text <- paste(head(neg_genes$gene, 8), collapse = ", ")
  terminal_text <- paste(head(pos_genes$gene, 8), collapse = ", ")
  subtitle <- if (isTRUE(status$has_branch_structure[[1]])) {
    "URD branch tree detected"
  } else {
    "No populated URD branch tree; figure shows pseudotime-ordered annotation states"
  }

  p <- ggplot() +
    geom_segment(data = edges, aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.8, arrow = arrow(length = unit(0.12, "inches"))) +
    geom_segment(data = ordering, aes(x = x, xend = x, y = y + 0.015, yend = label_y - 0.035), linewidth = 0.25, color = "grey55") +
    geom_point(data = ordering, aes(x = x, y = y, size = n_cells, color = annotation), alpha = 0.9) +
    geom_label(data = ordering, aes(x = x, y = label_y, label = label), size = 3.0, lineheight = 0.9, linewidth = 0, fill = "white", alpha = 0.95) +
    annotate("text", x = min(ordering$x), y = -0.18, hjust = 0, size = 3.2, label = paste("Root-associated genes:", root_text)) +
    annotate("text", x = max(ordering$x), y = -0.28, hjust = 1, size = 3.2, label = paste("Terminal-associated genes:", terminal_text)) +
    scale_size_continuous(range = c(3, 10)) +
    coord_cartesian(ylim = c(-0.36, 0.42), clip = "off") +
    theme_void(base_size = 11) +
    theme(
      legend.position = "none",
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(color = "black", face = "bold"),
      plot.subtitle = element_text(color = "black"),
      plot.margin = margin(20, 30, 50, 30)
    ) +
    labs(title = "DIV30 URD Lineage Decision-Tree Report", subtitle = subtitle)
  ggsave(plot_path, p, width = 12, height = 5, dpi = 240, bg = "white")
}

plot_umap_pseudotime <- function(meta, pt_values, annotation_col, plot_path) {
  if (!all(c("UMAP_1", "UMAP_2") %in% colnames(meta))) return(FALSE)
  df <- data.frame(
    UMAP_1 = as.numeric(meta$UMAP_1),
    UMAP_2 = as.numeric(meta$UMAP_2),
    pseudotime = pt_values,
    annotation = meta[[annotation_col]],
    stringsAsFactors = FALSE
  )
  df <- df[is.finite(df$UMAP_1) & is.finite(df$UMAP_2), , drop = FALSE]
  p <- ggplot(df, aes(UMAP_1, UMAP_2, color = pseudotime)) +
    geom_point(size = 0.35, alpha = 0.85) +
    coord_equal() +
    scale_color_viridis_c(na.value = "grey85") +
    theme_void(base_size = 11) +
    theme(plot.background = element_rect(fill = "white", color = NA)) +
    labs(title = "URD pseudotime on existing Seurat UMAP", color = "Pseudotime")
  ggsave(plot_path, p, width = 7, height = 6, dpi = 240, bg = "white")
  TRUE
}

diffusion_map_dataframe <- function(object, pt_values, annotation_col) {
  eig <- as.data.frame(object@dm@eigenvectors[, seq_len(min(4, ncol(object@dm@eigenvectors))), drop = FALSE])
  eig$cell_id <- rownames(object@dm@eigenvectors)
  eig$pseudotime <- pt_values[eig$cell_id]
  eig$annotation <- object@meta[eig$cell_id, annotation_col]
  eig
}

plot_diffusion_maps <- function(dm_df, pt_path, annotation_path) {
  p_pt <- ggplot(dm_df, aes(DC1, DC2, color = pseudotime)) +
    geom_point(size = 0.45, alpha = 0.85) +
    scale_color_viridis_c(na.value = "grey85") +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank()) +
    labs(title = "Diffusion map colored by URD pseudotime", color = "Pseudotime")
  ggsave(pt_path, p_pt, width = 7, height = 6, dpi = 240, bg = "white")

  p_anno <- ggplot(dm_df, aes(DC1, DC2, color = annotation)) +
    geom_point(size = 0.45, alpha = 0.85) +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(), legend.position = "right") +
    labs(title = "Diffusion map colored by annotation", color = NULL)
  ggsave(annotation_path, p_anno, width = 8, height = 6, dpi = 240, bg = "white")
}

flood_stability_summary <- function(object, pt_values) {
  stab <- object@pseudotime.stability
  if (!is.list(stab) || !all(c("pseudotime", "walks.per.cell") %in% names(stab))) {
    return(data.frame())
  }
  pt_stab <- stab$pseudotime
  walks <- stab$walks.per.cell
  common <- intersect(rownames(pt_stab), names(pt_values))
  pt_final <- pt_values[common]
  out <- lapply(colnames(pt_stab), function(col) {
    x <- as.numeric(pt_stab[common, col])
    ok <- is.finite(x) & is.finite(pt_final)
    delta <- abs(x[ok] - pt_final[ok])
    data.frame(
      stability_column = col,
      median_walks_per_cell = median(as.numeric(walks[common, col]), na.rm = TRUE),
      n_cells = sum(ok),
      spearman_to_final_pseudotime = suppressWarnings(cor(x[ok], pt_final[ok], method = "spearman")),
      median_abs_delta = median(delta, na.rm = TRUE),
      q95_abs_delta = as.numeric(quantile(delta, 0.95, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

plot_flood_stability <- function(stability_df, plot_path) {
  if (nrow(stability_df) == 0) return(FALSE)
  p <- ggplot(stability_df, aes(median_walks_per_cell)) +
    geom_line(aes(y = spearman_to_final_pseudotime), color = "#1f78b4", linewidth = 0.8) +
    geom_point(aes(y = spearman_to_final_pseudotime), color = "#1f78b4", size = 2) +
    geom_line(aes(y = 1 - median_abs_delta), color = "#e31a1c", linewidth = 0.8) +
    geom_point(aes(y = 1 - median_abs_delta), color = "#e31a1c", size = 2) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank()) +
    labs(
      title = "Flood pseudotime stability",
      subtitle = "Blue: Spearman correlation to final pseudotime; red: 1 - median absolute delta",
      x = "Walks per cell",
      y = "Stability metric"
    )
  ggsave(plot_path, p, width = 7, height = 5, dpi = 240, bg = "white")
  TRUE
}

plot_tree_status <- function(status, plot_path) {
  label <- if (isTRUE(status$has_branch_structure[[1]])) {
    "URD tree slot is populated.\nBranch-specific tree visualization can be added."
  } else {
    "No populated URD branch tree exists yet.\nThis object contains flood pseudotime, not a branch tree.\nBranch genes should wait for URD tree reconstruction."
  }
  p <- ggplot(data.frame(x = 0, y = 0, label = label), aes(x, y)) +
    geom_label(aes(label = label), size = 5, lineheight = 1.05, linewidth = 0.35, fill = "white") +
    xlim(-1, 1) +
    ylim(-1, 1) +
    theme_void(base_size = 12) +
    theme(plot.background = element_rect(fill = "white", color = NA)) +
    labs(title = "URD tree status")
  ggsave(plot_path, p, width = 7, height = 4, dpi = 240, bg = "white")
}

gene_cascade_matrix <- function(object, pt_values, neg_genes, pos_genes, n_each = 20L, n_bins = 20L) {
  genes <- unique(c(head(neg_genes$gene, n_each), head(pos_genes$gene, n_each)))
  genes <- intersect(genes, rownames(object@logupx.data))
  ok <- is.finite(pt_values)
  pt <- pt_values[ok]
  cells <- names(pt)
  breaks <- unique(stats::quantile(pt, probs = seq(0, 1, length.out = n_bins + 1L), na.rm = TRUE))
  if (length(breaks) < 3) return(data.frame())
  bins <- cut(pt, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  mat <- object@logupx.data[genes, cells, drop = FALSE]
  rows <- lapply(genes, function(gene) {
    x <- as.numeric(mat[gene, ])
    vals <- tapply(x, bins, mean, na.rm = TRUE)
    vals <- as.numeric(vals)
    if (stats::sd(vals, na.rm = TRUE) > 0) {
      z <- as.numeric(scale(vals))
    } else {
      z <- rep(0, length(vals))
    }
    data.frame(
      gene = gene,
      bin = seq_along(vals),
      mean_logupx = vals,
      z_mean_logupx = z,
      pseudotime_bin_midpoint = tapply(pt, bins, median, na.rm = TRUE),
      direction = if (gene %in% neg_genes$gene) "root_negative" else "terminal_positive",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

plot_gene_cascade_heatmap <- function(cascade_df, plot_path) {
  if (nrow(cascade_df) == 0) return(FALSE)
  gene_order <- unique(cascade_df$gene[order(cascade_df$direction, cascade_df$gene)])
  cascade_df$gene <- factor(cascade_df$gene, levels = rev(gene_order))
  p <- ggplot(cascade_df, aes(bin, gene, fill = z_mean_logupx)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, na.value = "grey90") +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(), axis.text.y = element_text(size = 7)) +
    labs(title = "Gene cascade across URD pseudotime", x = "Pseudotime bin", y = NULL, fill = "Scaled\nmean")
  ggsave(plot_path, p, width = 8, height = 9, dpi = 240, bg = "white")
  TRUE
}

write_markdown_report <- function(path, cfg, status, root_comp, ordering, neg, pos, outputs) {
  lines <- c(
    "# DIV30 URD Lineage Decision-Tree Report",
    "",
    paste0("- URD object: `", cfg$urd_rds, "`"),
    paste0("- Annotation column: `", cfg$annotation_col, "`"),
    paste0("- Pseudotime: `", cfg$pseudotime_name, "`"),
    paste0("- Correlation gene set: `", cfg$correlation_genes, "`"),
    "",
    "## Branch Status",
    "",
    if (isTRUE(status$has_branch_structure[[1]])) {
      "A populated URD tree slot was detected."
    } else {
      "No populated URD tree slot was detected. This object currently supports root-to-terminal pseudotime ordering, not branch decision-gene analysis."
    },
    "",
    "## Root Annotation Composition",
    "",
    paste(capture.output(print(root_comp, row.names = FALSE)), collapse = "\n"),
    "",
    "## Pseudotime Ordering By Annotation",
    "",
    paste(capture.output(print(ordering, row.names = FALSE)), collapse = "\n"),
    "",
    "## Top Root Genes",
    "",
    paste(capture.output(print(head(neg, 15), row.names = FALSE)), collapse = "\n"),
    "",
    "## Top Terminal Genes",
    "",
    paste(capture.output(print(head(pos, 15), row.names = FALSE)), collapse = "\n"),
    "",
    "## URD Validation Layers",
    "",
    "- Pseudotime UMAP: existing Seurat UMAP colored by URD pseudotime.",
    "- Diffusion map plots: DC1/DC2 colored by pseudotime and annotation.",
    "- Flood stability: convergence of intermediate flood pseudotime estimates toward final pseudotime.",
    "- Tree visualization: reports whether a populated URD branch tree exists.",
    "- Branch genes: only available after URD tree reconstruction; not inferred for a pseudotime-only object.",
    "- Gene cascade heatmap: top negative and positive pseudotime genes across pseudotime bins.",
    "",
    "## Output Files",
    "",
    paste0("- `", outputs, "`")
  )
  writeLines(lines, path)
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(opt$help)) {
  print_usage()
  quit(save = "no", status = 0)
}
if (is.null(opt$`urd-rds`) || !nzchar(opt$`urd-rds`)) stop("--urd-rds is required", call. = FALSE)
if (is.null(opt$outdir) || !nzchar(opt$outdir)) stop("--outdir is required", call. = FALSE)

cfg <- list(
  urd_rds = opt$`urd-rds`,
  outdir = opt$outdir,
  annotation_col = opt$`annotation-col`,
  pseudotime_name = opt$`pseudotime-name`,
  root_col = opt$`root-col`,
  top_n = as_int(opt$`top-n`, "top-n"),
  correlation_genes = opt$`correlation-genes`
)
table_dir <- file.path(cfg$outdir, "tables")
plot_dir <- file.path(cfg$outdir, "plots")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

log_msg("Reading URD object: ", cfg$urd_rds)
urd <- readRDS(cfg$urd_rds)
if (!inherits(urd, "URD")) stop("Input is not an URD object: ", cfg$urd_rds, call. = FALSE)

pt <- extract_pseudotime(urd, cfg$pseudotime_name)
cfg$pseudotime_name <- pt$name
meta <- urd@meta[rownames(urd@meta), , drop = FALSE]
names(pt$values) <- rownames(meta)

log_msg("Computing root annotation composition")
root_comp <- root_annotation_composition(meta, cfg$annotation_col, cfg$root_col)
log_msg("Computing pseudotime ordering")
ordering <- pseudotime_ordering(meta, pt$values, cfg$annotation_col)
log_msg("Computing gene-pseudotime correlations on ", cfg$correlation_genes, " genes")
gene_cor <- gene_correlations(urd, pt$values, cfg$correlation_genes)
gene_cor <- gene_cor[is.finite(gene_cor$spearman_with_pseudotime), , drop = FALSE]
neg <- head(gene_cor[order(gene_cor$spearman_with_pseudotime), , drop = FALSE], cfg$top_n)
pos <- head(gene_cor[order(gene_cor$spearman_with_pseudotime, decreasing = TRUE), , drop = FALSE], cfg$top_n)
status <- branch_status(urd)

decision_genes <- data.frame(
  status = if (isTRUE(status$has_branch_structure[[1]])) "branch_tree_detected_not_implemented" else "not_available_no_branch_tree",
  comparison = "",
  gene = "",
  statistic = NA_real_,
  note = status$note,
  stringsAsFactors = FALSE
)

paths <- c(
  root_comp = file.path(table_dir, "root_annotation_composition.tsv"),
  ordering = file.path(table_dir, "pseudotime_ordering_by_annotation.tsv"),
  negative = file.path(table_dir, "top_negative_pseudotime_genes.tsv"),
  positive = file.path(table_dir, "top_positive_pseudotime_genes.tsv"),
  branch_status = file.path(table_dir, "branch_structure_status.tsv"),
  decision_genes = file.path(table_dir, "decision_genes_between_branches.tsv"),
  flood_stability = file.path(table_dir, "flood_stability_summary.tsv"),
  cascade = file.path(table_dir, "gene_cascade_heatmap_matrix.tsv"),
  umap = file.path(plot_dir, "umap_pseudotime.png"),
  dm_pseudotime = file.path(plot_dir, "diffusion_map_pseudotime.png"),
  dm_annotation = file.path(plot_dir, "diffusion_map_annotation.png"),
  flood_stability_plot = file.path(plot_dir, "flood_stability.png"),
  tree_visualization = file.path(plot_dir, "tree_visualization.png"),
  cascade_plot = file.path(plot_dir, "gene_cascade_heatmap.png"),
  figure = file.path(plot_dir, "lineage_decision_tree.png"),
  report = file.path(cfg$outdir, "lineage_decision_tree_report.md")
)

stability_df <- flood_stability_summary(urd, pt$values)
cascade_df <- gene_cascade_matrix(urd, pt$values, neg, pos)
dm_df <- diffusion_map_dataframe(urd, pt$values, cfg$annotation_col)

write_tsv(root_comp, paths[["root_comp"]])
write_tsv(ordering, paths[["ordering"]])
write_tsv(neg, paths[["negative"]])
write_tsv(pos, paths[["positive"]])
write_tsv(status, paths[["branch_status"]])
write_tsv(decision_genes, paths[["decision_genes"]])
write_tsv(stability_df, paths[["flood_stability"]])
write_tsv(cascade_df, paths[["cascade"]])
invisible(plot_umap_pseudotime(meta, pt$values, cfg$annotation_col, paths[["umap"]]))
invisible(plot_diffusion_maps(dm_df, paths[["dm_pseudotime"]], paths[["dm_annotation"]]))
invisible(plot_flood_stability(stability_df, paths[["flood_stability_plot"]]))
invisible(plot_tree_status(status, paths[["tree_visualization"]]))
invisible(plot_gene_cascade_heatmap(cascade_df, paths[["cascade_plot"]]))
make_decision_figure(ordering, neg, pos, status, paths[["figure"]])
write_markdown_report(paths[["report"]], cfg, status, root_comp, ordering, neg, pos, unname(paths))

log_msg("Wrote report: ", paths[["report"]])
log_msg("Wrote decision-tree figure: ", paths[["figure"]])
