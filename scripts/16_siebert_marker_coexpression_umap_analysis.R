#!/usr/bin/env Rscript

# Non-final Siebert 2026 marker/coexpression UMAP analysis.
# Extracts requested genes from the Seurat object as log1p(CP10K) from counts,
# plots cluster/marker UMAPs, then overlays sequential marker-positive states.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

PNG_DPI <- 450

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = "")))
  flush.console()
}

get_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  opt <- list(
    project_root = get_env("PROJECT_ROOT", "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder"),
    seurat_rds = get_env("SIEBERT_SEURAT_RDS", ""),
    out_dir = get_env("SIEBERT_MARKER_COEXPR_OUT_DIR", ""),
    study_id = get_env("MARKER_COEXPR_STUDY_ID", "siebert_2026"),
    study_label = get_env("MARKER_COEXPR_STUDY_LABEL", "Siebert 2026"),
    sample_col = get_env("MARKER_COEXPR_SAMPLE_COL", "sample"),
    assay = get_env("SIEBERT_MARKER_COEXPR_ASSAY", "RNA"),
    reduction = get_env("SIEBERT_MARKER_COEXPR_REDUCTION", "umap"),
    cluster_col = get_env("SIEBERT_MARKER_COEXPR_CLUSTER_COL", "seurat_clusters"),
    output_prefix = get_env("MARKER_COEXPR_OUTPUT_PREFIX", "siebert_marker"),
    threshold = as.numeric(get_env("SIEBERT_MARKER_COEXPR_THRESHOLD", "0.5"))
  )

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    value <- if (i < length(args)) args[[i + 1L]] else ""
    if (key == "--project-root") {
      opt$project_root <- value; i <- i + 2L
    } else if (key == "--seurat-rds") {
      opt$seurat_rds <- value; i <- i + 2L
    } else if (key == "--out-dir") {
      opt$out_dir <- value; i <- i + 2L
    } else if (key == "--study-id") {
      opt$study_id <- value; i <- i + 2L
    } else if (key == "--study-label") {
      opt$study_label <- value; i <- i + 2L
    } else if (key == "--sample-col") {
      opt$sample_col <- value; i <- i + 2L
    } else if (key == "--assay") {
      opt$assay <- value; i <- i + 2L
    } else if (key == "--reduction") {
      opt$reduction <- value; i <- i + 2L
    } else if (key == "--cluster-col") {
      opt$cluster_col <- value; i <- i + 2L
    } else if (key == "--output-prefix") {
      opt$output_prefix <- value; i <- i + 2L
    } else if (key == "--threshold") {
      opt$threshold <- as.numeric(value); i <- i + 2L
    } else {
      stop("Unknown argument: ", key, call. = FALSE)
    }
  }

  if (!nzchar(opt$seurat_rds)) {
    opt$seurat_rds <- file.path(opt$project_root, "results/siebert_2026/siebert_2026_seurat.rds")
  }
  if (!nzchar(opt$out_dir)) {
    opt$out_dir <- file.path(
      opt$project_root,
      "results/siebert_2026/analysis/siebert_marker_coexpression_umap_v3_sample_breakdown"
    )
  }
  if (!is.finite(opt$threshold)) stop("--threshold must be numeric", call. = FALSE)
  opt
}

write_tsv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

write_tsv_gz <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

theme_umap <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0),
      legend.title = element_text(size = base_size - 1),
      legend.text = element_text(size = base_size - 1),
      plot.margin = margin(3, 3, 3, 3)
    )
}

wrap_title <- function(x, width = 24) {
  paste(strwrap(x, width = width), collapse = "\n")
}

mixed_cluster_levels <- function(values) {
  values <- unique(as.character(values))
  numeric_like <- grepl("^[0-9]+$", values)
  c(sort(as.integer(values[numeric_like])), sort(values[!numeric_like]))
}

gene_aliases <- list(
  DCX = c("DCX"),
  PCDH19 = c("PCDH19"),
  VIM = c("VIM"),
  MKI67 = c("MKI67", "KI67", "KI67M", "Ki67", "Ki67m"),
  SOX2 = c("SOX2")
)

display_gene <- c(
  DCX = "DCX",
  PCDH19 = "PCDH19",
  VIM = "VIM",
  MKI67 = "Ki67 (MKI67)",
  SOX2 = "SOX2"
)

match_gene <- function(requested_gene, features) {
  aliases <- gene_aliases[[requested_gene]]
  if (is.null(aliases)) aliases <- requested_gene
  feature_upper <- toupper(features)
  for (alias in aliases) {
    exact_idx <- match(alias, features)
    if (!is.na(exact_idx)) {
      return(list(gene = requested_gene, matched_feature = features[[exact_idx]], match_type = "exact"))
    }
    upper_idx <- match(toupper(alias), feature_upper)
    if (!is.na(upper_idx)) {
      return(list(gene = requested_gene, matched_feature = features[[upper_idx]], match_type = "case_insensitive"))
    }
  }
  list(gene = requested_gene, matched_feature = "", match_type = "missing")
}

get_assay_layers <- function(obj, assay, layer) {
  if (!assay %in% Assays(obj)) stop("Assay not found: ", assay, call. = FALSE)
  assay_obj <- obj[[assay]]
  if (inherits(assay_obj, "Assay5")) {
    layers <- Layers(assay_obj)
    hits <- layers[layers == layer]
    if (length(hits) == 0L) hits <- layers[grepl(paste0("^", layer, "(\\.|$)"), layers)]
    if (length(hits) == 0L) {
      stop("No ", assay, "/", layer, " layer found", call. = FALSE)
    }
    out <- lapply(hits, function(hit) LayerData(assay_obj, layer = hit))
    names(out) <- hits
    return(out)
  }
  out <- tryCatch(
    GetAssayData(obj, assay = assay, layer = layer),
    error = function(err) GetAssayData(obj, assay = assay, slot = layer)
  )
  list(layer = out)
}

count_layer_cells <- function(count_layers) {
  unique(unlist(lapply(count_layers, colnames), use.names = FALSE))
}

count_layer_features <- function(count_layers) {
  unique(unlist(lapply(count_layers, rownames), use.names = FALSE))
}

extract_log1p_cp10k <- function(count_layers, matched_features, cells) {
  genes <- names(matched_features)
  totals <- rep(NA_real_, length(cells))
  names(totals) <- cells
  marker_counts_dense <- matrix(
    0,
    nrow = length(cells),
    ncol = length(genes),
    dimnames = list(cells, genes)
  )

  for (layer_name in names(count_layers)) {
    counts <- count_layers[[layer_name]]
    layer_cells <- intersect(colnames(counts), cells)
    if (length(layer_cells) == 0L) next
    totals[layer_cells] <- as.numeric(Matrix::colSums(counts[, layer_cells, drop = FALSE]))

    present_genes <- genes[matched_features[genes] %in% rownames(counts)]
    if (length(present_genes) == 0L) next
    sub <- counts[unname(matched_features[present_genes]), layer_cells, drop = FALSE]
    marker_counts_dense[layer_cells, present_genes] <- as.matrix(t(sub))
  }

  totals[!is.finite(totals) | totals <= 0] <- NA_real_
  expr <- t(t(marker_counts_dense) / totals) * 10000
  expr[is.na(expr)] <- 0
  expr <- log1p(expr)
  list(
    expr = as.data.frame(expr, check.names = FALSE),
    raw_marker_counts = marker_counts_dense,
    total_counts = totals
  )
}

plot_cluster_umap <- function(df) {
  label_df <- aggregate(cbind(umap_1, umap_2) ~ cluster, data = df, FUN = median)
  ggplot(df, aes(umap_1, umap_2, color = cluster)) +
    geom_point(size = 0.11, alpha = 0.86, stroke = 0) +
    geom_text(
      data = label_df,
      aes(x = umap_1, y = umap_2, label = cluster),
      color = "black",
      size = 2.1,
      fontface = "bold",
      inherit.aes = FALSE
    ) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1), ncol = 2)) +
    labs(title = "Clusters", color = "Cluster") +
    theme_umap()
}

sample_levels <- function(values) {
  preferred <- c(
    "Young_1", "Young_2", "Old_1", "Old_2",
    "9583-MW-1", "9583-MW-2", "9583-MW-3", "9583-MW-4", "9583-MW-5", "9583-MW-6"
  )
  values <- unique(as.character(values))
  c(intersect(preferred, values), sort(setdiff(values, preferred)))
}

sample_palette <- function(values) {
  levels <- sample_levels(values)
  preset <- c(
    "Young_1" = "#1f78b4",
    "Young_2" = "#33a02c",
    "Old_1" = "#e31a1c",
    "Old_2" = "#ff7f00",
    "9583-MW-1" = "#1f78b4",
    "9583-MW-2" = "#33a02c",
    "9583-MW-3" = "#e31a1c",
    "9583-MW-4" = "#ff7f00",
    "9583-MW-5" = "#6a3d9a",
    "9583-MW-6" = "#b15928"
  )
  missing <- setdiff(levels, names(preset))
  if (length(missing) > 0L) {
    fallback <- hue_pal()(length(missing))
    names(fallback) <- missing
    preset <- c(preset, fallback)
  }
  preset[levels]
}

plot_sample_umap <- function(df) {
  plot_df <- df
  plot_df$sample <- factor(plot_df$sample, levels = sample_levels(plot_df$sample))
  label_df <- aggregate(cbind(umap_1, umap_2) ~ sample, data = plot_df, FUN = median)
  ggplot(plot_df, aes(umap_1, umap_2, color = sample)) +
    geom_point(size = 0.11, alpha = 0.82, stroke = 0) +
    geom_text(
      data = label_df,
      aes(x = umap_1, y = umap_2, label = sample),
      color = "black",
      size = 2.1,
      fontface = "bold",
      inherit.aes = FALSE
    ) +
    scale_color_manual(values = sample_palette(plot_df$sample), na.value = "#777777") +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    labs(title = "Samples", color = "Sample") +
    theme_umap()
}

plot_expression_umap <- function(df, gene) {
  bg <- df
  plot_df <- df[df[[gene]] > 0, , drop = FALSE]
  plot_df <- plot_df[order(plot_df[[gene]], na.last = TRUE), , drop = FALSE]
  vmax <- as.numeric(stats::quantile(plot_df[[gene]], probs = 0.99, na.rm = TRUE))
  if (!is.finite(vmax) || vmax <= 0) vmax <- max(df[[gene]], na.rm = TRUE)
  if (!is.finite(vmax) || vmax <= 0) vmax <- 1
  plot_df$plot_value <- pmin(plot_df[[gene]], vmax)

  ggplot() +
    geom_point(data = bg, aes(umap_1, umap_2), color = "#c9c9c9", size = 0.055, alpha = 0.52, stroke = 0) +
    geom_point(data = plot_df, aes(umap_1, umap_2, color = plot_value), size = 0.13, alpha = 0.98, stroke = 0) +
    scale_color_gradientn(
      colors = c("#1b2a6b", "#1f78b4", "#2fb47c", "#fee08b", "#f46d43", "#b2182b"),
      limits = c(0, vmax),
      oob = scales::squish,
      name = "log1p(CP10K)"
    ) +
    labs(title = unname(display_gene[[gene]])) +
    theme_umap()
}

plot_state_umap <- function(df, state_col, title, color) {
  bg <- df[!df[[state_col]], , drop = FALSE]
  fg <- df[df[[state_col]], , drop = FALSE]
  ggplot() +
    geom_point(data = bg, aes(umap_1, umap_2), color = "#c7c7c7", size = 0.055, alpha = 0.42, stroke = 0) +
    geom_point(data = fg, aes(umap_1, umap_2), color = color, size = 0.15, alpha = 0.98, stroke = 0) +
    labs(title = wrap_title(title, width = 22)) +
    theme_umap()
}

plot_summary_bars <- function(summary_df) {
  bar_colors <- c(
    "DCX+ Ki67+" = "#0f9b8e",
    "DCX+ Ki67+ SOX2+" = "#e08d2d",
    "DCX+ Ki67+ SOX2+ PCDH19+" = "#b13f8a"
  )
  ggplot(summary_df, aes(x = state_label, y = percent_cells, fill = state_label)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.2) +
    geom_text(aes(label = sprintf("%.2f%%", percent_cells)), hjust = -0.08, size = 2.3) +
    coord_flip(clip = "off") +
    scale_y_continuous(limits = c(0, max(5, max(summary_df$percent_cells, na.rm = TRUE) * 1.22)), expand = c(0, 0)) +
    scale_fill_manual(values = bar_colors) +
    labs(title = "Sequential percentage", x = NULL, y = "% cells") +
    theme_classic(base_size = 8) +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 7),
      axis.text.x = element_text(size = 7),
      axis.title.x = element_text(size = 8),
      plot.title = element_text(face = "bold", size = 9),
      plot.margin = margin(4, 12, 4, 3)
    )
}

plot_all_four_callout <- function(summary_df) {
  all_row <- summary_df[summary_df$state_id == "dcx_mki67_sox2_pcdh19", , drop = FALSE]
  label <- sprintf(
    "All four coexpressed\n%s / %s cells\n%.2f%%",
    comma(all_row$n_positive),
    comma(all_row$n_total),
    all_row$percent_cells
  )
  ggplot(data.frame(x = 0, y = 0), aes(x, y)) +
    annotate("text", x = 0, y = 0.14, label = label, size = 3.4, fontface = "bold", lineheight = 1.05) +
    annotate("text", x = 0, y = -0.36, label = "cutoff: expression >= 0.5", size = 2.3, color = "#555555") +
    xlim(-1, 1) +
    ylim(-1, 1) +
    labs(title = "Final coexpression") +
    theme_void(base_size = 8) +
    theme(
      plot.title = element_text(face = "bold", size = 9, hjust = 0),
      plot.margin = margin(4, 3, 4, 3)
    )
}

plot_sample_breakdown <- function(sample_summary) {
  plot_df <- sample_summary
  plot_df$sample <- factor(plot_df$sample, levels = rev(sample_levels(plot_df$sample)))
  plot_df$state_label <- factor(
    plot_df$state_label,
    levels = c("DCX+ Ki67+", "DCX+ Ki67+ SOX2+", "DCX+ Ki67+ SOX2+ PCDH19+")
  )
  ggplot(plot_df, aes(x = state_label, y = sample, fill = percent_cells)) +
    geom_tile(color = "white", linewidth = 0.45) +
    geom_text(aes(label = sprintf("%.1f%%", percent_cells)), size = 2.2, fontface = "bold") +
    scale_fill_gradientn(
      colors = c("#f2f2f2", "#a6cee3", "#1f78b4", "#08306b"),
      name = "% cells"
    ) +
    labs(title = "Sample breakdown", x = NULL, y = NULL) +
    theme_classic(base_size = 8) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, size = 6.3),
      axis.text.y = element_text(size = 7),
      axis.ticks = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6.5),
      plot.title = element_text(face = "bold", size = 9),
      plot.margin = margin(4, 3, 4, 3)
    )
}

save_figure <- function(plot, out_base, width, height) {
  ggsave(paste0(out_base, ".png"), plot = plot, width = width, height = height, units = "in", dpi = PNG_DPI, bg = "white")
  ggsave(paste0(out_base, ".pdf"), plot = plot, width = width, height = height, units = "in", bg = "white")
  svg_path <- paste0(out_base, ".svg")
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(svg_path, plot = plot, width = width, height = height, units = "in", bg = "white")
  } else {
    grDevices::svg(svg_path, width = width, height = height, bg = "white", onefile = FALSE)
    print(plot)
    grDevices::dev.off()
  }
}

main <- function() {
  opt <- parse_args()
  dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
  table_dir <- file.path(opt$out_dir, "tables")
  plot_dir <- file.path(opt$out_dir, "plots")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  prefix <- opt$output_prefix
  log_msg("Reading ", opt$study_label, " Seurat object: ", opt$seurat_rds)
  obj <- readRDS(opt$seurat_rds)
  if (!opt$reduction %in% Reductions(obj)) stop("Reduction not found: ", opt$reduction, call. = FALSE)

  log_msg("Extracting UMAP, metadata, and counts layer")
  emb <- Embeddings(obj, reduction = opt$reduction)
  meta <- obj@meta.data
  count_layers <- get_assay_layers(obj, assay = opt$assay, layer = "counts")
  features <- count_layer_features(count_layers)

  requested_genes <- c("DCX", "PCDH19", "VIM", "MKI67", "SOX2")
  matches <- lapply(requested_genes, match_gene, features = features)
  match_table <- do.call(
    rbind,
    lapply(matches, function(x) {
      data.frame(
        requested_gene = x$gene,
        display_gene = unname(display_gene[[x$gene]]),
        matched_feature = x$matched_feature,
        match_type = x$match_type,
        stringsAsFactors = FALSE
      )
    })
  )
  write_tsv(match_table, file.path(table_dir, paste0(prefix, "_gene_matches.tsv")))
  if (any(!nzchar(match_table$matched_feature))) {
    stop("Missing requested marker(s): ", paste(match_table$requested_gene[!nzchar(match_table$matched_feature)], collapse = ","), call. = FALSE)
  }

  assay_layers <- if (inherits(obj[[opt$assay]], "Assay5")) paste(Layers(obj[[opt$assay]]), collapse = ",") else "counts,data,scale.data"
  object_qc <- data.frame(
    key = c("study_id", "study_label", "seurat_rds", "assays", "assay_used", "assay_used_layers", "reduction_used", "sample_col", "cluster_col", "n_cells", "n_features"),
    value = c(
      opt$study_id,
      opt$study_label,
      opt$seurat_rds,
      paste(Assays(obj), collapse = ","),
      opt$assay,
      assay_layers,
      opt$reduction,
      opt$sample_col,
      opt$cluster_col,
      as.character(ncol(obj)),
      as.character(nrow(obj))
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(object_qc, file.path(table_dir, paste0(prefix, "_seurat_object_qc.tsv")))

  cells <- Reduce(intersect, list(count_layer_cells(count_layers), rownames(emb), rownames(meta)))
  if (length(cells) == 0L) stop("No common cells across counts, UMAP, and metadata", call. = FALSE)
  emb <- emb[cells, , drop = FALSE]
  meta <- meta[cells, , drop = FALSE]
  extracted <- extract_log1p_cp10k(
    count_layers = count_layers,
    matched_features = setNames(match_table$matched_feature, match_table$requested_gene),
    cells = cells
  )
  expr <- extracted$expr
  colnames(expr) <- requested_genes

  marker_qc <- do.call(
    rbind,
    lapply(requested_genes, function(gene) {
      matched_feature <- match_table$matched_feature[match_table$requested_gene == gene][[1]]
      raw_counts <- extracted$raw_marker_counts[, gene]
      raw_nonzero <- sum(raw_counts > 0)
      vals <- expr[[gene]]
      qs <- stats::quantile(vals, probs = c(0, 0.5, 0.9, 0.99, 1), na.rm = TRUE)
      data.frame(
        requested_gene = gene,
        matched_feature = matched_feature,
        raw_total_counts = as.numeric(sum(raw_counts, na.rm = TRUE)),
        raw_nonzero_cells = raw_nonzero,
        raw_pct_nonzero_cells = 100 * raw_nonzero / length(cells),
        log1p_cp10k_min = unname(qs[[1]]),
        log1p_cp10k_median = unname(qs[[2]]),
        log1p_cp10k_q90 = unname(qs[[3]]),
        log1p_cp10k_q99 = unname(qs[[4]]),
        log1p_cp10k_max = unname(qs[[5]]),
        pct_cells_ge_threshold = 100 * sum(vals >= opt$threshold, na.rm = TRUE) / length(vals),
        stringsAsFactors = FALSE
      )
    })
  )
  write_tsv(marker_qc, file.path(table_dir, paste0(prefix, "_object_expression_qc.tsv")))

  cluster <- if (opt$cluster_col %in% colnames(meta)) as.character(meta[[opt$cluster_col]]) else as.character(Idents(obj)[cells])
  cluster[is.na(cluster) | !nzchar(cluster)] <- "NA"
  cluster <- factor(cluster, levels = mixed_cluster_levels(cluster))

  df <- data.frame(
    cell_id = cells,
    study_id = opt$study_id,
    study_label = opt$study_label,
    cluster = cluster,
    sample = if (opt$sample_col %in% colnames(meta)) as.character(meta[[opt$sample_col]]) else "",
    orig_ident = if ("orig.ident" %in% colnames(meta)) as.character(meta$orig.ident) else "",
    predictions = if ("predictions" %in% colnames(meta)) as.character(meta$predictions) else "",
    predictions_class = if ("predictions_class" %in% colnames(meta)) as.character(meta$predictions_class) else "",
    umap_1 = emb[, 1],
    umap_2 = emb[, 2],
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  df <- cbind(df, expr)

  threshold <- opt$threshold
  df$dcx_mki67 <- df$DCX >= threshold & df$MKI67 >= threshold
  df$dcx_mki67_sox2 <- df$dcx_mki67 & df$SOX2 >= threshold
  df$dcx_mki67_sox2_pcdh19 <- df$dcx_mki67_sox2 & df$PCDH19 >= threshold

  state_specs <- data.frame(
    state_id = c("dcx_mki67", "dcx_mki67_sox2", "dcx_mki67_sox2_pcdh19"),
    state_label = c("DCX+ Ki67+", "DCX+ Ki67+ SOX2+", "DCX+ Ki67+ SOX2+ PCDH19+"),
    genes = c("DCX,MKI67", "DCX,MKI67,SOX2", "DCX,MKI67,SOX2,PCDH19"),
    color = c("#0f9b8e", "#e08d2d", "#b13f8a"),
    stringsAsFactors = FALSE
  )
  summary_df <- do.call(
    rbind,
    lapply(seq_len(nrow(state_specs)), function(i) {
      state_id <- state_specs$state_id[[i]]
      n_pos <- sum(df[[state_id]], na.rm = TRUE)
      data.frame(
        state_id = state_id,
        state_label = state_specs$state_label[[i]],
        genes = state_specs$genes[[i]],
        threshold = threshold,
        n_positive = n_pos,
        n_total = nrow(df),
        percent_cells = 100 * n_pos / nrow(df),
        stringsAsFactors = FALSE
      )
    })
  )
  summary_df$state_label <- factor(summary_df$state_label, levels = rev(state_specs$state_label))
  write_tsv(summary_df, file.path(table_dir, paste0(prefix, "_coexpression_summary.tsv")))

  cluster_summary <- do.call(
    rbind,
    lapply(levels(df$cluster), function(cluster_id) {
      sub <- df[df$cluster == cluster_id, , drop = FALSE]
      do.call(
        rbind,
        lapply(seq_len(nrow(state_specs)), function(i) {
          state_id <- state_specs$state_id[[i]]
          n_pos <- sum(sub[[state_id]], na.rm = TRUE)
          data.frame(
            cluster = cluster_id,
            state_id = state_id,
            state_label = state_specs$state_label[[i]],
            threshold = threshold,
            n_positive = n_pos,
            n_total = nrow(sub),
            percent_cells = if (nrow(sub) > 0L) 100 * n_pos / nrow(sub) else NA_real_,
            stringsAsFactors = FALSE
          )
        })
      )
    })
  )
  write_tsv(cluster_summary, file.path(table_dir, paste0(prefix, "_coexpression_by_cluster.tsv")))

  sample_summary <- do.call(
    rbind,
    lapply(sample_levels(df$sample), function(sample_id) {
      sub <- df[df$sample == sample_id, , drop = FALSE]
      do.call(
        rbind,
        lapply(seq_len(nrow(state_specs)), function(i) {
          state_id <- state_specs$state_id[[i]]
          n_pos <- sum(sub[[state_id]], na.rm = TRUE)
          data.frame(
            sample = sample_id,
            state_id = state_id,
            state_label = state_specs$state_label[[i]],
            threshold = threshold,
            n_positive = n_pos,
            n_total = nrow(sub),
            percent_cells = if (nrow(sub) > 0L) 100 * n_pos / nrow(sub) else NA_real_,
            stringsAsFactors = FALSE
          )
        })
      )
    })
  )
  write_tsv(sample_summary, file.path(table_dir, paste0(prefix, "_coexpression_by_sample.tsv")))
  write_tsv_gz(df, file.path(table_dir, paste0(prefix, "_expression_and_coexpression_per_cell.tsv.gz")))

  log_msg("Rendering UMAP panels")
  top_plots <- list(
    plot_cluster_umap(df),
    plot_expression_umap(df, "DCX"),
    plot_expression_umap(df, "PCDH19"),
    plot_expression_umap(df, "VIM"),
    plot_expression_umap(df, "MKI67"),
    plot_sample_umap(df)
  )
  bottom_plots <- list(
    plot_state_umap(df, "dcx_mki67", "DCX+ Ki67+", "#0f9b8e"),
    plot_state_umap(df, "dcx_mki67_sox2", "DCX+ Ki67+ SOX2+", "#e08d2d"),
    plot_state_umap(df, "dcx_mki67_sox2_pcdh19", "DCX+ Ki67+ SOX2+ PCDH19+", "#b13f8a"),
    plot_summary_bars(summary_df),
    plot_all_four_callout(summary_df),
    plot_sample_breakdown(sample_summary)
  )
  combined <- wrap_plots(c(top_plots, bottom_plots), ncol = 6, guides = "collect") +
    plot_annotation(
      title = paste0(opt$study_label, " marker expression and sequential coexpression on UMAP"),
      subtitle = "Expression values are log1p(CP10K) from RNA counts; positive cutoff is expression >= 0.5"
    ) &
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, color = "#4b4b4b"),
      legend.position = "right"
    )
  save_figure(combined, file.path(plot_dir, paste0(prefix, "_coexpression_umap_grid")), width = 20, height = 7.3)

  manifest <- data.frame(
    key = c(
      "seurat_rds", "out_dir", "assay", "reduction", "cluster_col", "threshold",
      "combined_png", "combined_pdf", "combined_svg", "per_cell_table", "summary_table", "sample_summary_table"
    ),
    value = c(
      opt$seurat_rds, opt$out_dir, opt$assay, opt$reduction, opt$cluster_col, as.character(threshold),
      file.path(plot_dir, paste0(prefix, "_coexpression_umap_grid.png")),
      file.path(plot_dir, paste0(prefix, "_coexpression_umap_grid.pdf")),
      file.path(plot_dir, paste0(prefix, "_coexpression_umap_grid.svg")),
      file.path(table_dir, paste0(prefix, "_expression_and_coexpression_per_cell.tsv.gz")),
      file.path(table_dir, paste0(prefix, "_coexpression_summary.tsv")),
      file.path(table_dir, paste0(prefix, "_coexpression_by_sample.tsv"))
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(manifest, file.path(opt$out_dir, paste0(prefix, "_coexpression_manifest.tsv")))
  log_msg("Done: ", opt$out_dir)
}

main()
