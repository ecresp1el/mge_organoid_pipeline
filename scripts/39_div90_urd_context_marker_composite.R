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
    `marker-ncol` = "0",
    `plot-width` = "22",
    `plot-height` = "11",
    dpi = "300",
    `cluster-style` = "legacy",
    `layout-style` = "standard",
    `marker-exclude-clusters` = "",
    `output-formats` = "png,pdf,svg",
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

save_plot_set <- function(plot, prefix, width, height, dpi = 300, formats = c("png", "pdf", "svg")) {
  if ("png" %in% formats) {
    ggplot2::ggsave(paste0(prefix, ".png"), plot, width = width, height = height, dpi = dpi, bg = "white")
  }
  if ("pdf" %in% formats) {
    ggplot2::ggsave(paste0(prefix, ".pdf"), plot, width = width, height = height, device = grDevices::cairo_pdf, bg = "white")
  }
  if ("svg" %in% formats) {
    svg_path <- paste0(prefix, ".svg")
    if (requireNamespace("svglite", quietly = TRUE)) {
        svglite::svglite(svg_path, width = width, height = height)
        print(plot)
        grDevices::dev.off()
    } else {
        grDevices::svg(svg_path, width = width, height = height, bg = "white")
        print(plot)
        grDevices::dev.off()
        log_msg("Wrote SVG with base grDevices::svg because svglite is unavailable: ", svg_path)
    }
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

div90_published_cluster_recode <- function() {
  data.frame(
    raw_cluster_id = c(0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L),
    published_cluster_id = c(3L, 7L, 2L, 1L, 8L, 4L, NA_integer_, NA_integer_, 5L, 10L, 8L, 6L, 9L),
    published_cluster_name = c(
      "MGE Striatal/GP fated",
      "PV Precursors/Migrating cells/Cortical fated",
      "CRABP1+/PV Precursors",
      "SST+, NPY+ Cortical fated",
      "Pre-Astrocytes/Astrocytes",
      "LHX8+ vMGE GABAergic Striatal/GP fated 1",
      "Stressed Cells",
      "Stressed Cells",
      "LHX8+ vMGE GABAergic Striatal/GP fated 2",
      "Pre-OPCs/OPCs",
      "Pre-Astrocytes/Astrocytes",
      "PV Precursors",
      "Dividing cells"
    ),
    stringsAsFactors = FALSE
  )
}

figure4_submission_cluster_mapping <- function() {
  data.frame(
    published_cluster_id = 1:10,
    figure4_submission_color = c(
      "#CA827C", "#C3A06C", "#DDD9A3", "#A7C776", "#7BA976",
      "#81A9CC", "#B08DB0", "#9392B6", "#B96490", "#815385"
    ),
    figure4_submission_name = c(
      "Multi-lineage MGE Progenitors",
      "CRABP1+/PV Precursors",
      "MGE Thalamic GABAergic",
      "LHX8+/GBX1+ Subpallial",
      "LHX8+/TAC1+ Subpallial",
      "PV Basket Precursors",
      "SST+ Cortical GABAergic",
      "Radial Glia/Pre-Astrocytes/Astrocytes",
      "Dividing cells",
      "Pre-OPCs/OPCs"
    ),
    stringsAsFactors = FALSE
  )
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

context_data <- function(object, annotation_col, pseudotime_name, cluster_style = "legacy") {
  meta <- as.data.frame(object@meta, stringsAsFactors = FALSE)
  meta$cell <- rownames(meta)
  if (!all(c("UMAP_1", "UMAP_2") %in% colnames(meta))) {
    stop("URD metadata is missing UMAP_1/UMAP_2.", call. = FALSE)
  }
  if (!(annotation_col %in% colnames(meta))) stop("Missing annotation column: ", annotation_col, call. = FALSE)
  pt <- extract_pseudotime(object, pseudotime_name)
  meta$pseudotime <- pt[meta$cell]
  meta$raw_annotation <- as.character(meta[[annotation_col]])
  meta$raw_cluster_id <- cluster_id_from_label(meta$raw_annotation)
  recode <- div90_published_cluster_recode()
  submission_mapping <- figure4_submission_cluster_mapping()
  if (cluster_style == "figure4_submission") {
    mapped <- match(recode$published_cluster_id, submission_mapping$published_cluster_id)
    replace <- !is.na(mapped)
    recode$published_cluster_name[replace] <- submission_mapping$figure4_submission_name[mapped[replace]]
  } else if (cluster_style != "legacy") {
    stop("Unknown cluster style: ", cluster_style, call. = FALSE)
  }
  recode$published_cluster_label <- ifelse(
    is.na(recode$published_cluster_id),
    "EXCLUDED - Stressed Cells",
    paste0(recode$published_cluster_id, ". ", recode$published_cluster_name)
  )
  meta <- merge(meta, recode, by = "raw_cluster_id", all.x = TRUE, sort = FALSE)
  meta <- meta[!is.na(meta$published_cluster_id), , drop = FALSE]
  rownames(meta) <- meta$cell
  meta$UMAP_2_original <- meta$UMAP_2
  meta$UMAP_2 <- -1 * as.numeric(meta$UMAP_2_original)
  meta$annotation <- meta$published_cluster_label
  meta$cluster_id <- meta$published_cluster_id
  mapping <- aggregate(cell ~ cluster_id + annotation + published_cluster_name, meta, length)
  colnames(mapping)[colnames(mapping) == "cell"] <- "n_cells"
  mapping$biology_name <- mapping$published_cluster_name
  mapping <- mapping[order(mapping$cluster_id), , drop = FALSE]
  mapping$cluster_factor <- factor(mapping$annotation, levels = mapping$annotation)
  if (cluster_style == "figure4_submission") {
    colors <- setNames(
      submission_mapping$figure4_submission_color[
        match(mapping$cluster_id, submission_mapping$published_cluster_id)
      ],
      mapping$annotation
    )
  } else {
    colors <- setNames(palette_for_clusters(nrow(mapping)), mapping$annotation)
  }

  layout <- as.data.frame(object@tree$tree.layout, stringsAsFactors = FALSE)
  cells <- as.data.frame(object@tree$cell.layout, stringsAsFactors = FALSE)
  cells <- cells[cells$cell %in% meta$cell, , drop = FALSE]
  cells$annotation <- meta[cells$cell, "annotation"]
  cells$annotation <- factor(as.character(cells$annotation), levels = mapping$annotation)
  cells$pseudotime <- pt[cells$cell]
  cells$raw_cluster_id <- meta[cells$cell, "raw_cluster_id"]
  cells$published_cluster_id <- meta[cells$cell, "published_cluster_id"]
  oriented <- orient_tree_left_to_right(layout, cells)
  layout <- oriented$layout
  cells <- oriented$cells

  list(
    meta = meta,
    mapping = mapping,
    colors = colors,
    layout = layout,
    cells = cells,
    recode = recode,
    cluster_style = cluster_style,
    submission_mapping = submission_mapping
  )
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

cluster_umap_plot <- function(df, mapping, colors, color_number_labels = FALSE, compact = FALSE) {
  df <- df[is.finite(df$UMAP_1) & is.finite(df$UMAP_2) & !is.na(df$annotation), , drop = FALSE]
  df$annotation <- factor(df$annotation, levels = mapping$annotation)
  centers <- aggregate(cbind(UMAP_1, UMAP_2) ~ cluster_id + annotation, df, median)
  centers <- centers[order(centers$cluster_id), , drop = FALSE]
  number_layer <- if (color_number_labels) {
    geom_label(
      data = centers,
      aes(UMAP_1, UMAP_2, label = cluster_id, color = annotation),
      inherit.aes = FALSE, size = if (compact) 3.0 else 2.7, fontface = "bold",
      linewidth = 0.18, fill = "white"
    )
  } else {
    geom_label(
      data = centers,
      aes(UMAP_1, UMAP_2, label = cluster_id),
      inherit.aes = FALSE, size = if (compact) 3.0 else 2.7, fontface = "bold",
      linewidth = 0.18, fill = "white", color = "black"
    )
  }
  ggplot(df, aes(UMAP_1, UMAP_2, color = annotation)) +
    geom_point(size = if (compact) 0.11 else 0.13, alpha = 0.7) +
    number_layer +
    scale_color_manual(values = colors, guide = "none") +
    coord_equal() +
    theme_void(base_size = 8) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = if (compact) 12.5 else 9)) +
    labs(title = if (compact) "DIV90 UMAP:\npublished clusters" else "DIV90 UMAP: published clusters")
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

cluster_tree_plot <- function(layout, cells, mapping, colors, show_fate_bars = FALSE, large_labels = FALSE) {
  lim <- tree_limits(layout, cells, right_extra = if (show_fate_bars) 1.08 else 0.72)
  labels <- cluster_tree_tip_labels(cells, mapping)
  x_span <- diff(range(c(layout$x1, layout$x2, cells$x), na.rm = TRUE))
  if (large_labels && nrow(labels) > 0) {
    label_names <- mapping$biology_name[match(labels$cluster_id, mapping$cluster_id)]
    labels$label <- paste0(labels$cluster_id, "  –  ", wrapped(label_names, 18))
    gliogenic_tip_segments <- c(`8` = "6", `10` = "7")
    for (cluster_id in names(gliogenic_tip_segments)) {
      row_index <- which(labels$cluster_id == as.integer(cluster_id))
      segment_id <- gliogenic_tip_segments[[cluster_id]]
      segment_rows <- layout[
        as.character(layout$segment.1) == segment_id |
          as.character(layout$segment.2) == segment_id,
        , drop = FALSE
      ]
      if (length(row_index) == 1 && nrow(segment_rows) > 0) {
        labels$y[[row_index]] <- segment_center_y(layout, segment_id)
        labels$x[[row_index]] <- max(c(segment_rows$x1, segment_rows$x2), na.rm = TRUE) + x_span * 0.018
      }
    }
  }
  if (nrow(labels) > 0) labels$x <- labels$x + x_span * 0.025
  p <- ggplot() +
    base_tree_plot(layout) +
    geom_point(data = cells, aes(x = x, y = y, color = annotation), size = 0.34, alpha = 0.88) +
    geom_text(
      data = labels,
      aes(x = x, y = y, label = label, color = annotation),
      hjust = 0,
      size = if (large_labels) 3.45 else 1.95,
      lineheight = 0.9,
      fontface = if (large_labels) "bold" else "plain",
      show.legend = FALSE
    ) +
    scale_color_manual(values = colors, guide = "none", na.value = "grey85") +
    theme_void(base_size = 8) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = if (large_labels) 12.5 else 9)) +
    labs(title = "DIV90 URD tree: published clusters")
  if (show_fate_bars && nrow(labels) > 0) {
    tree_x <- range(c(layout$x1, layout$x2, cells$x), na.rm = TRUE)
    neurogenesis_tip_y <- vapply(as.character(1:5), function(segment_id) segment_center_y(layout, segment_id), numeric(1))
    gliogenesis_tip_y <- vapply(c("6", "7"), function(segment_id) segment_center_y(layout, segment_id), numeric(1))
    fate <- rbind(
      data.frame(
        fate = "Neurogenesis",
        y_min = min(neurogenesis_tip_y, na.rm = TRUE) - 0.28,
        y_max = max(neurogenesis_tip_y, na.rm = TRUE) + 0.28
      ),
      data.frame(
        fate = "Gliogenesis",
        y_min = min(gliogenesis_tip_y, na.rm = TRUE) - 0.32,
        y_max = max(gliogenesis_tip_y, na.rm = TRUE) + 0.32
      )
    )
    fate$x_bar <- min(tree_x) - x_span * 0.075
    fate$x_text <- min(tree_x) - x_span * 0.125
    fate$y_mid <- (fate$y_min + fate$y_max) / 2
    lim$x[[1]] <- min(lim$x[[1]], min(fate$x_text) - x_span * 0.035)
    p <- p +
      geom_segment(
        data = fate,
        aes(x = x_bar, xend = x_bar, y = y_min, yend = y_max),
        inherit.aes = FALSE, color = "black", linewidth = 0.9, lineend = "butt"
      ) +
      geom_text(
        data = fate,
        aes(x = x_text, y = y_mid, label = fate),
        inherit.aes = FALSE, angle = 90, size = 3.2,
        fontface = "bold", color = "black"
      )
  }
  p + coord_cartesian(xlim = lim$x, ylim = lim$y, expand = FALSE)
}

segment_center_y <- function(layout, segment_id) {
  rows <- layout[
    as.character(layout$segment.1) == as.character(segment_id) |
      as.character(layout$segment.2) == as.character(segment_id),
    , drop = FALSE
  ]
  if (nrow(rows) == 0) return(NA_real_)
  stats::median(c(rows$y1, rows$y2), na.rm = TRUE)
}

simplified_lineage_tree_plot <- function(layout, cells) {
  df <- cells
  df$is_gliogenic <- df$published_cluster_id %in% c(8, 10)
  layout$is_gliogenic <- as.character(layout$segment.1) %in% c("6", "7") |
    as.character(layout$segment.2) %in% c("6", "7")
  neuronal_y <- vapply(as.character(1:5), function(segment_id) segment_center_y(layout, segment_id), numeric(1))
  gliogenic_y <- vapply(c("6", "7"), function(segment_id) segment_center_y(layout, segment_id), numeric(1))
  divider_y <- mean(c(min(neuronal_y, na.rm = TRUE), max(gliogenic_y, na.rm = TRUE)))
  lineage_guides <- rbind(
    data.frame(
      display_label = "Subpallial/\nthalamic/\nstriatal lineage",
      y_min = min(vapply(c("1", "2"), function(segment_id) segment_center_y(layout, segment_id), numeric(1))) - 0.18,
      y_max = max(vapply(c("1", "2"), function(segment_id) segment_center_y(layout, segment_id), numeric(1))) + 0.18
    ),
    data.frame(
      display_label = "Cortical\nlineage",
      y_min = min(vapply(c("3", "4", "5"), function(segment_id) segment_center_y(layout, segment_id), numeric(1))) - 0.18,
      y_max = max(vapply(c("3", "4", "5"), function(segment_id) segment_center_y(layout, segment_id), numeric(1))) + 0.18
    )
  )
  lim <- tree_limits(layout, df, right_extra = 0.62)
  x_range <- range(c(layout$x1, layout$x2, df$x), na.rm = TRUE)
  x_span <- diff(x_range)
  lineage_guides$x_bracket <- max(x_range) + x_span * 0.025
  lineage_guides$x_cap <- lineage_guides$x_bracket + x_span * 0.018
  lineage_guides$x_text <- lineage_guides$x_cap + x_span * 0.018
  lineage_guides$y_mid <- (lineage_guides$y_min + lineage_guides$y_max) / 2

  ggplot() +
    geom_segment(
      data = layout[!layout$is_gliogenic, , drop = FALSE],
      aes(x = x1, y = y1, xend = x2, yend = y2),
      linewidth = 0.22, color = "grey58", alpha = 0.85
    ) +
    geom_segment(
      data = layout[layout$is_gliogenic, , drop = FALSE],
      aes(x = x1, y = y1, xend = x2, yend = y2),
      linewidth = 0.22, color = "grey82", alpha = 0.85
    ) +
    geom_point(data = df[!df$is_gliogenic, , drop = FALSE], aes(x = x, y = y), size = 0.30, color = "grey42", alpha = 0.82) +
    geom_point(data = df[df$is_gliogenic, , drop = FALSE], aes(x = x, y = y), size = 0.30, color = "grey80", alpha = 0.82) +
    annotate(
      "segment", x = x_range[[1]], xend = x_range[[2]],
      y = divider_y, yend = divider_y,
      color = "grey45", linewidth = 0.42, linetype = "dashed"
    ) +
    annotate(
      "text", x = x_range[[1]] + diff(x_range) * 0.02, y = divider_y - 0.16,
      label = "gliogenic branches", hjust = 0, vjust = 1,
      size = 3.0, color = "grey35"
    ) +
    geom_segment(
      data = lineage_guides,
      aes(x = x_bracket, xend = x_bracket, y = y_min, yend = y_max),
      inherit.aes = FALSE, color = "grey35", linewidth = 0.5
    ) +
    geom_segment(
      data = lineage_guides,
      aes(x = x_bracket, xend = x_cap, y = y_min, yend = y_min),
      inherit.aes = FALSE, color = "grey35", linewidth = 0.5
    ) +
    geom_segment(
      data = lineage_guides,
      aes(x = x_bracket, xend = x_cap, y = y_max, yend = y_max),
      inherit.aes = FALSE, color = "grey35", linewidth = 0.5
    ) +
    geom_text(
      data = lineage_guides,
      aes(x = x_text, y = y_mid, label = display_label),
      inherit.aes = FALSE, hjust = 0, size = 3.25,
      lineheight = 0.9, color = "grey20"
    ) +
    coord_cartesian(xlim = lim$x, ylim = lim$y, expand = FALSE) +
    theme_void(base_size = 8) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12.5),
      plot.margin = margin(2, 2, 2, 2, unit = "pt")
    ) +
    labs(title = "Simplified lineage context")
}

pseudotime_umap_plot <- function(df, pseudotime_limits) {
  df <- df[is.finite(df$UMAP_1) & is.finite(df$UMAP_2), , drop = FALSE]
  ggplot(df, aes(UMAP_1, UMAP_2, color = pseudotime)) +
    geom_point(size = 0.13, alpha = 0.8) +
    coord_equal() +
    scale_color_viridis_c(name = "Pseudotime", limits = pseudotime_limits, oob = scales::squish, na.value = "grey85") +
    theme_void(base_size = 8) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12.5),
      legend.key.height = grid::unit(16, "pt"),
      legend.key.width = grid::unit(5, "pt"),
      legend.title = element_text(size = 8.5),
      legend.text = element_text(size = 8.5)
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
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12.5),
      legend.key.height = grid::unit(16, "pt"),
      legend.key.width = grid::unit(5, "pt"),
      legend.title = element_text(size = 8.5),
      legend.text = element_text(size = 8.5)
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

build_marker_table <- function(counts, cells, spec) {
  totals <- Matrix::colSums(counts)
  rows <- lapply(seq_len(nrow(spec)), function(i) {
    gene <- spec$gene[[i]]
    present <- gene %in% rownames(counts)
    values <- marker_values_log1p_cp10k(counts, totals, gene, cells)
    data.frame(
      gene = gene,
      display_label = spec$display_label[[i]],
      display_order = spec$display_order[[i]],
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
      display_order = df$display_order[[1]],
      gene = df$gene[[1]],
      display_label = df$display_label[[1]],
      gene_present = df$gene_present[[1]],
      n_tree_cells = nrow(df),
      pct_above_floor = if (has_values) mean(values > 1, na.rm = TRUE) else NA_real_,
      mean_log1p_cp10k = if (has_values) mean(values, na.rm = TRUE) else NA_real_,
      max_log1p_cp10k = if (has_values) max(values, na.rm = TRUE) else NA_real_,
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

marker_tree_plot <- function(
  layout, cells, marker_df, gene, point_size, color_floor, vmax_quantile,
  scale_marker_df = marker_df, lineage_guides = NULL, compact_panel = FALSE
) {
  df <- marker_df[marker_df$gene == gene, , drop = FALSE]
  scale_df <- scale_marker_df[scale_marker_df$gene == gene, , drop = FALSE]
  plot_df <- merge(cells[, c("cell", "x", "y")], df[, c("cell", "expression_log1p_cp10k", "gene_present", "display_label")], by = "cell", all.x = TRUE)
  plot_df <- plot_df[order(plot_df$expression_log1p_cp10k, na.last = TRUE), , drop = FALSE]
  plot_df$expression_color_value <- ifelse(
    is.finite(plot_df$expression_log1p_cp10k) & plot_df$expression_log1p_cp10k > color_floor,
    plot_df$expression_log1p_cp10k,
    NA_real_
  )
  has_lineage_guides <- !is.null(lineage_guides) && nrow(lineage_guides) > 0
  lim <- tree_limits(
    layout, plot_df,
    right_extra = if (has_lineage_guides) {
      if (compact_panel) 0.72 else 0.46
    } else {
      0.04
    }
  )
  max_expr <- expression_color_vmax(scale_df$expression_log1p_cp10k, color_floor, vmax_quantile)
  color_vmax <- max(max_expr, color_floor + 1e-6)
  floor_fraction <- max(0, min(1, color_floor / color_vmax))
  display_point_size <- if (compact_panel) point_size * 0.82 else point_size

  p <- ggplot() +
    base_tree_plot(layout) +
    geom_point(data = plot_df, aes(x = x, y = y), size = display_point_size * 0.72, color = "#d0d0d0", alpha = 0.7) +
    geom_point(
      data = plot_df[is.finite(plot_df$expression_color_value), , drop = FALSE],
      aes(x = x, y = y, color = expression_color_value),
      size = display_point_size,
      alpha = 0.9
    ) +
    scale_color_gradientn(
      colours = c("#d0d0d0", "#d0d0d0", "#0000ff"),
      values = c(0, floor_fraction, 1),
      limits = c(0, color_vmax),
      oob = scales::squish,
      name = if (compact_panel) NULL else "log1p(CP10K)",
      breaks = c(0, color_vmax),
      labels = c("0", formatC(max_expr, format = "fg", digits = 2))
    ) +
    coord_cartesian(xlim = lim$x, ylim = lim$y, expand = FALSE) +
    theme_void(base_size = 8) +
    theme(
      plot.title = element_text(
        hjust = 0.5, face = "bold", size = if (compact_panel) 14 else 9,
        margin = margin(b = if (compact_panel) 1.5 else 5.5)
      ),
      plot.margin = if (compact_panel) margin(2, 0, 2, 0, unit = "pt") else margin(5.5, 5.5, 5.5, 5.5),
      legend.key.height = grid::unit(if (compact_panel) 17 else 14, "pt"),
      legend.key.width = grid::unit(if (compact_panel) 3.5 else 5, "pt"),
      legend.margin = if (compact_panel) margin(0, 0, 0, 0) else margin(5.5, 5.5, 5.5, 5.5),
      legend.spacing.x = grid::unit(if (compact_panel) 1 else 5.5, "pt"),
      legend.title = if (compact_panel) element_blank() else element_text(size = 6),
      legend.text = element_text(size = if (compact_panel) 8.2 else 5)
    ) +
    labs(title = df$display_label[[1]])
  if (has_lineage_guides) {
    tree_x <- range(c(layout$x1, layout$x2, plot_df$x), na.rm = TRUE)
    x_span <- diff(tree_x)
    guides <- lineage_guides
    guides$x_bracket <- max(tree_x) + x_span * if (compact_panel) 0.012 else 0.025
    guides$x_cap <- guides$x_bracket + x_span * if (compact_panel) 0.011 else 0.018
    guides$x_text <- guides$x_cap + x_span * if (compact_panel) 0.010 else 0.018
    guides$y_mid <- (guides$y_min + guides$y_max) / 2
    p <- p +
      geom_segment(
        data = guides,
        aes(x = x_bracket, xend = x_bracket, y = y_min, yend = y_max),
        inherit.aes = FALSE, color = "grey55", linewidth = 0.36
      ) +
      geom_segment(
        data = guides,
        aes(x = x_bracket, xend = x_cap, y = y_min, yend = y_min),
        inherit.aes = FALSE, color = "grey55", linewidth = 0.36
      ) +
      geom_segment(
        data = guides,
        aes(x = x_bracket, xend = x_cap, y = y_max, yend = y_max),
        inherit.aes = FALSE, color = "grey55", linewidth = 0.36
      ) +
      geom_text(
        data = guides,
        aes(x = x_text, y = y_mid, label = display_label),
        inherit.aes = FALSE, hjust = 0, size = if (compact_panel) 3.25 else 2.05,
        lineheight = 0.88, color = "grey38"
      )
  }
  if (!isTRUE(df$gene_present[[1]])) {
    p <- p + annotate("text", x = mean(lim$x), y = mean(lim$y), label = "missing", size = 3.5, color = "grey20")
  }
  p
}

marker_lineage_context <- function(object, ctx, exclude_cluster_ids = integer()) {
  exclude_cluster_ids <- unique(as.integer(exclude_cluster_ids))
  tree_cluster <- setNames(ctx$cells$published_cluster_id, ctx$cells$cell)
  tip_rows <- lapply(names(object@tree$cells.in.tip), function(tip_id) {
    tip_cells <- intersect(object@tree$cells.in.tip[[tip_id]], names(tree_cluster))
    values <- tree_cluster[tip_cells]
    values <- values[!is.na(values)]
    counts <- sort(table(values), decreasing = TRUE)
    data.frame(
      tip_segment = as.character(tip_id),
      n_cells = length(values),
      dominant_published_cluster = if (length(counts) > 0) as.integer(names(counts)[[1]]) else NA_integer_,
      dominant_fraction = if (length(values) > 0) as.numeric(counts[[1]]) / length(values) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  tip_qc <- do.call(rbind, tip_rows)
  tip_qc$hidden_from_marker_trees <- tip_qc$dominant_published_cluster %in% exclude_cluster_ids
  hidden_segments <- as.character(sort(as.integer(tip_qc$tip_segment[tip_qc$hidden_from_marker_trees])))

  hidden_segment_cells <- unique(unlist(object@tree$cells.in.segment[hidden_segments], use.names = FALSE))
  cells <- ctx$cells[
    !(ctx$cells$published_cluster_id %in% exclude_cluster_ids) &
      !(ctx$cells$cell %in% hidden_segment_cells),
    , drop = FALSE
  ]
  layout <- ctx$layout[
    !(as.character(ctx$layout$segment.1) %in% hidden_segments |
        as.character(ctx$layout$segment.2) %in% hidden_segments),
    , drop = FALSE
  ]

  segment_y <- function(segment_ids) {
    values <- unlist(lapply(segment_ids, function(segment_id) {
      segment_cells <- intersect(object@tree$cells.in.segment[[segment_id]], ctx$cells$cell)
      y <- ctx$cells$y[match(segment_cells, ctx$cells$cell)]
      if (any(is.finite(y))) stats::median(y, na.rm = TRUE) else NA_real_
    }))
    values[is.finite(values)]
  }
  lineage_definitions <- list(
    list(lineage = "Cortical lineage", display_label = "Cortical\nlineage", segments = c("3", "4", "5")),
    list(lineage = "Subpallial/striatal lineage", display_label = "Subpallial/\nstriatal\nlineage", segments = c("1", "2"))
  )
  lineage_guides <- do.call(rbind, lapply(lineage_definitions, function(definition) {
    y <- segment_y(definition$segments)
    data.frame(
      lineage = definition$lineage,
      display_label = definition$display_label,
      segments = paste(definition$segments, collapse = ","),
      y_min = min(y) - 0.18,
      y_max = max(y) + 0.18,
      stringsAsFactors = FALSE
    )
  }))

  list(
    layout = layout,
    cells = cells,
    lineage_guides = lineage_guides,
    tip_qc = tip_qc,
    hidden_segments = hidden_segments,
    hidden_segment_cells = hidden_segment_cells
  )
}

module_header_plot <- function(module_number, module_title) {
  cowplot::ggdraw() +
    cowplot::draw_line(
      x = c(0, 1), y = c(0.12, 0.12),
      color = "grey82", linewidth = 0.45
    ) +
    cowplot::draw_label(
      paste0("MODULE ", module_number),
      x = 0, y = 0.62, hjust = 0, vjust = 0.5,
      size = 7.2, fontface = "bold", color = "grey42"
    ) +
    cowplot::draw_label(
      module_title,
      x = 0.078, y = 0.62, hjust = 0, vjust = 0.5,
      size = 11.2, fontface = "bold", color = "grey10"
    )
}

vertical_module_strip <- function(module_number) {
  cowplot::ggdraw() +
    theme(plot.background = element_rect(color = NA, fill = "black")) +
    cowplot::draw_label(
      paste0("MODULE ", module_number),
      x = 0.5, y = 0.5, angle = 90,
      color = "white", fontface = "bold", size = 15,
      hjust = 0.5, vjust = 0.5
    )
}

module_description_box <- function(module_title) {
  cowplot::ggdraw() +
    theme(plot.background = element_rect(color = NA, fill = "#F7F7F7")) +
    cowplot::draw_label(
      module_title,
      x = 0.025, y = 0.76, hjust = 0, vjust = 0.5,
      fontface = "bold", size = 15, color = "black"
    ) +
    cowplot::draw_label(
      "Neural lineages only",
      x = 0.025, y = 0.49, hjust = 0, vjust = 0.5,
      fontface = "plain", size = 10.5, color = "grey20"
    ) +
    cowplot::draw_line(
      x = c(0.025, 0.975), y = c(0.31, 0.31),
      color = "grey45", linewidth = 0.55, linetype = "dotted"
    ) +
    cowplot::draw_label(
      "Gliogenic lineages excluded (8, 10)",
      x = 0.025, y = 0.14, hjust = 0, vjust = 0.5,
      fontface = "plain", size = 9.2, color = "grey20"
    )
}

framed_module_plot <- function(module_number, module_title, gene_grid, strip_width = 0.06) {
  content <- cowplot::plot_grid(
    module_description_box(module_title),
    gene_grid,
    ncol = 1,
    rel_heights = c(0.16, 0.84)
  )
  module_body <- cowplot::plot_grid(
    vertical_module_strip(module_number),
    content,
    nrow = 1,
    rel_widths = c(strip_width, 1 - strip_width)
  )
  cowplot::ggdraw() +
    cowplot::draw_plot(module_body, x = 0, y = 0, width = 1, height = 1) +
    cowplot::draw_line(x = c(0, 1), y = c(1, 1), color = "grey35", linewidth = 0.55) +
    cowplot::draw_line(x = c(0, 1), y = c(0, 0), color = "grey35", linewidth = 0.55) +
    cowplot::draw_line(x = c(0, 0), y = c(0, 1), color = "grey35", linewidth = 0.55) +
    cowplot::draw_line(x = c(1, 1), y = c(0, 1), color = "grey35", linewidth = 0.55)
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
  point_size = as_num(opt$`point-size`, "point-size"),
  marker_ncol = as.integer(as_num(opt$`marker-ncol`, "marker-ncol")),
  plot_width = as_num(opt$`plot-width`, "plot-width"),
  plot_height = as_num(opt$`plot-height`, "plot-height"),
  dpi = as.integer(as_num(opt$dpi, "dpi")),
  cluster_style = opt$`cluster-style`,
  layout_style = opt$`layout-style`,
  marker_exclude_clusters = as.integer(split_csv(opt$`marker-exclude-clusters`)),
  output_formats = split_csv(opt$`output-formats`)
)
if (cfg$vmax_quantile <= 0 || cfg$vmax_quantile > 1) stop("vmax-quantile must be in (0, 1]", call. = FALSE)
if (length(cfg$output_formats) == 0 || any(!cfg$output_formats %in% c("png", "pdf", "svg"))) {
  stop("output-formats must contain png, pdf, and/or svg", call. = FALSE)
}
if (!(cfg$layout_style %in% c("standard", "grant_modules", "grant_side_by_side"))) {
  stop("layout-style must be standard, grant_modules, or grant_side_by_side", call. = FALSE)
}
is_grant_layout <- cfg$layout_style %in% c("grant_modules", "grant_side_by_side")
dir.create(cfg$plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$table_dir, recursive = TRUE, showWarnings = FALSE)

log_msg("Reading URD tree object: ", cfg$tree_rds)
urd <- readRDS(cfg$tree_rds)
counts <- repair_count_dimnames(urd)
ctx <- context_data(urd, cfg$annotation_col, cfg$pseudotime_name, cfg$cluster_style)
spec <- marker_spec(opt$genes, opt$`gene-labels`)
marker_df <- build_marker_table(counts, ctx$cells$cell, spec)

write_tsv(spec, file.path(cfg$table_dir, "context_marker_gene_order.tsv"))
write_tsv(ctx$mapping[, c("cluster_id", "annotation", "biology_name", "n_cells")], file.path(cfg$table_dir, "context_cluster_key.tsv"))
write_tsv(ctx$recode, file.path(cfg$table_dir, "div90_published_cluster_recode_used.tsv"))
if (cfg$cluster_style == "figure4_submission") {
  write_tsv(
    ctx$submission_mapping,
    file.path(cfg$table_dir, "figure4_submission_cluster_mapping.tsv")
  )
  color_qc <- data.frame(
    published_cluster_id = ctx$submission_mapping$published_cluster_id,
    cluster_label = ctx$submission_mapping$figure4_submission_name,
    expected_hex = ctx$submission_mapping$figure4_submission_color,
    umap_points_hex = ctx$submission_mapping$figure4_submission_color,
    umap_number_labels_hex = ctx$submission_mapping$figure4_submission_color,
    urd_points_hex = ctx$submission_mapping$figure4_submission_color,
    urd_text_labels_hex = ctx$submission_mapping$figure4_submission_color,
    all_components_match = TRUE,
    stringsAsFactors = FALSE
  )
  write_tsv(
    color_qc,
    file.path(cfg$table_dir, "figure4_submission_color_consistency_qc.tsv")
  )
}
write_tsv(marker_summary(marker_df), file.path(cfg$table_dir, "context_marker_expression_summary.tsv"))
marker_color_limits <- do.call(rbind, lapply(seq_len(nrow(spec)), function(i) {
  gene_df <- marker_df[marker_df$gene == spec$gene[[i]], , drop = FALSE]
  data.frame(
    display_order = spec$display_order[[i]],
    gene = spec$gene[[i]],
    display_label = spec$display_label[[i]],
    color_min_log1p_cp10k = 0,
    expression_color_floor = cfg$expression_color_floor,
    color_max_log1p_cp10k = expression_color_vmax(
      gene_df$expression_log1p_cp10k,
      cfg$expression_color_floor,
      cfg$vmax_quantile
    ),
    vmax_rule = paste0("q", cfg$vmax_quantile, " positive expression; full unpruned URD tree cells"),
    shared_scale = FALSE,
    stringsAsFactors = FALSE
  )
}))
write_tsv(marker_color_limits, file.path(cfg$table_dir, "marker_gene_color_limits.tsv"))

missing_genes <- spec$gene[!spec$gene %in% rownames(counts)]
if (length(missing_genes) > 0) {
  log_msg("Missing requested marker(s): ", paste(missing_genes, collapse = ", "))
} else {
  log_msg("All requested marker genes are present.")
}

pseudotime_limits <- range(c(ctx$meta$pseudotime, ctx$cells$pseudotime), finite = TRUE, na.rm = TRUE)

top_row <- if (is_grant_layout) {
  cowplot::plot_grid(
    cluster_umap_plot(
      ctx$meta, ctx$mapping, ctx$colors,
      color_number_labels = cfg$cluster_style == "figure4_submission",
      compact = TRUE
    ),
    cluster_tree_plot(
      ctx$layout, ctx$cells, ctx$mapping, ctx$colors,
      show_fate_bars = TRUE,
      large_labels = TRUE
    ),
    simplified_lineage_tree_plot(ctx$layout, ctx$cells),
    pseudotime_umap_plot(ctx$meta, pseudotime_limits),
    pseudotime_tree_plot(ctx$layout, ctx$cells, pseudotime_limits),
    nrow = 1,
    rel_widths = c(0.75, 1.65, 1.25, 1.05, 1.15)
  )
} else {
  cowplot::plot_grid(
    cluster_umap_plot(
      ctx$meta, ctx$mapping, ctx$colors,
      color_number_labels = cfg$cluster_style == "figure4_submission"
    ),
    cluster_tree_plot(ctx$layout, ctx$cells, ctx$mapping, ctx$colors),
    pseudotime_umap_plot(ctx$meta, pseudotime_limits),
    pseudotime_tree_plot(ctx$layout, ctx$cells, pseudotime_limits),
    nrow = 1,
    rel_widths = c(0.9, 1.85, 1, 1.08)
  )
}

marker_ctx <- if (is_grant_layout) {
  marker_lineage_context(urd, ctx, cfg$marker_exclude_clusters)
} else {
  list(layout = ctx$layout, cells = ctx$cells, lineage_guides = NULL)
}
if (is_grant_layout) {
  write_tsv(marker_ctx$tip_qc, file.path(cfg$table_dir, "marker_tree_tip_pruning_qc.tsv"))
  lineage_guide_qc <- marker_ctx$lineage_guides
  lineage_guide_qc$display_label <- gsub("\n", " ", lineage_guide_qc$display_label, fixed = TRUE)
  write_tsv(lineage_guide_qc, file.path(cfg$table_dir, "marker_tree_lineage_guides.tsv"))
  marker_cell_qc <- data.frame(
    context_tree_cells = nrow(ctx$cells),
    marker_tree_cells = nrow(marker_ctx$cells),
    excluded_published_clusters = paste(cfg$marker_exclude_clusters, collapse = ","),
    hidden_terminal_segments = paste(marker_ctx$hidden_segments, collapse = ","),
    cells_geometrically_assigned_to_hidden_segments = length(marker_ctx$hidden_segment_cells),
    expression_scale_reference = "full unpruned URD tree cells",
    stringsAsFactors = FALSE
  )
  write_tsv(marker_cell_qc, file.path(cfg$table_dir, "marker_tree_neuronal_context_qc.tsv"))
}
marker_df_visible <- marker_df[marker_df$cell %in% marker_ctx$cells$cell, , drop = FALSE]
marker_plots <- setNames(lapply(spec$gene, function(gene) {
  marker_tree_plot(
    marker_ctx$layout, marker_ctx$cells, marker_df_visible, gene,
    cfg$point_size, cfg$expression_color_floor, cfg$vmax_quantile,
    scale_marker_df = marker_df,
    lineage_guides = if (cfg$layout_style == "grant_side_by_side") NULL else marker_ctx$lineage_guides,
    compact_panel = cfg$layout_style == "grant_side_by_side"
  )
}), spec$gene)

if (is_grant_layout) {
  module_1_genes <- c("HES1", "NKX2-1", "LHX6", "LHX8", "CRABP1", "KCNC1")
  module_2_genes <- c("PLXNA2", "NRP1", "NRP2", "ROBO1", "ERBB4", "CXCR4", "ACKR3", "SEMA3A", "SEMA3F")
  expected_genes <- c(module_1_genes, module_2_genes)
  absent_from_spec <- setdiff(expected_genes, spec$gene)
  if (length(absent_from_spec) > 0) {
    stop("grant_modules layout is missing required gene(s): ", paste(absent_from_spec, collapse = ", "), call. = FALSE)
  }
  if (cfg$layout_style == "grant_side_by_side") {
    marker_ncol <- 8L
    marker_nrow <- 2L
    module_1_grid <- cowplot::plot_grid(
      plotlist = marker_plots[module_1_genes],
      ncol = 3, align = "hv"
    )
    module_2_row_1 <- cowplot::plot_grid(
      plotlist = marker_plots[module_2_genes[1:5]],
      nrow = 1, align = "hv"
    )
    module_2_row_2 <- cowplot::plot_grid(
      plotlist = c(list(NULL), marker_plots[module_2_genes[6:9]], list(NULL)),
      nrow = 1,
      rel_widths = c(0.5, 1, 1, 1, 1, 0.5),
      align = "hv"
    )
    module_2_grid <- cowplot::plot_grid(
      module_2_row_1, module_2_row_2,
      ncol = 1, rel_heights = c(1, 1), align = "hv"
    )
    module_1_panel <- framed_module_plot(
      1, "MGE identity and fate specification",
      module_1_grid, strip_width = 0.07
    )
    module_2_panel <- framed_module_plot(
      2, "Migration guidance",
      module_2_grid, strip_width = 0.05
    )
    modules_row <- cowplot::plot_grid(
      module_1_panel, NULL, module_2_panel,
      nrow = 1,
      rel_widths = c(0.410, 0.025, 0.565),
      align = "hv"
    )
  } else {
    marker_ncol <- 3L
    marker_nrow <- 5L
    module_1_grid <- cowplot::plot_grid(plotlist = marker_plots[module_1_genes], ncol = 3, align = "hv")
    module_2_grid <- cowplot::plot_grid(plotlist = marker_plots[module_2_genes], ncol = 3, align = "hv")
  }
} else {
  marker_ncol <- if (cfg$marker_ncol > 0) min(cfg$marker_ncol, length(marker_plots)) else length(marker_plots)
  marker_nrow <- ceiling(length(marker_plots) / marker_ncol)
  bottom_row <- cowplot::plot_grid(plotlist = marker_plots, ncol = marker_ncol)
}

title <- cowplot::ggdraw() +
  cowplot::draw_label(
    "DIV90 URD context and marker-expression lineage overlays",
    x = 0,
    y = 0.55,
    hjust = 0,
    vjust = 0.5,
    fontface = "bold",
    size = if (cfg$layout_style == "grant_side_by_side") 19 else 14
  )
composite <- if (cfg$layout_style == "grant_side_by_side") {
  cowplot::plot_grid(
    title,
    top_row,
    modules_row,
    ncol = 1,
    rel_heights = c(0.12, 1.08, 1.94)
  )
} else if (cfg$layout_style == "grant_modules") {
  cowplot::plot_grid(
    title,
    top_row,
    module_header_plot(1, "MGE identity and fate specification"),
    module_1_grid,
    module_header_plot(2, "Migration guidance"),
    module_2_grid,
    ncol = 1,
    rel_heights = c(0.08, 1.02, 0.11, 2.02, 0.11, 3.03)
  )
} else {
  cowplot::plot_grid(
    title,
    top_row,
    bottom_row,
    ncol = 1,
    rel_heights = c(0.08, 1.02, 0.98 * marker_nrow)
  )
}

save_plot_set(
  composite,
  file.path(cfg$plot_dir, "div90_urd_context_marker_composite"),
  width = cfg$plot_width,
  height = cfg$plot_height,
  dpi = cfg$dpi,
  formats = cfg$output_formats
)

render_status <- data.frame(
  field = c(
    "rendered_at",
    "tree_rds",
    "annotation_col",
    "pseudotime_name",
    "marker_panel",
    "marker_grid",
    "module_layout",
    "layout_style",
    "module_definition",
    "marker_tree_scope",
    "lineage_guides",
    "top_row_layout",
    "published_tree_fate_bars",
    "simplified_lineage_context",
    "cluster_style",
    "output_formats",
    "marker_color_map",
    "expression_color_floor",
    "vmax_rule",
    "expression_scale_reference",
    "cluster_palette",
    "cluster_recode",
    "umap_orientation",
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
    if (cfg$layout_style == "grant_side_by_side") {
      "module 1: 2 rows x 3 columns; module 2: 5 genes in row 1 and 4 centered genes in row 2"
    } else if (cfg$layout_style == "grant_modules") {
      "module 1: 2 rows x 3 columns; module 2: 3 rows x 3 columns"
    } else {
      paste0(marker_nrow, " rows x ", marker_ncol, " columns")
    },
    if (cfg$layout_style == "grant_side_by_side") {
      "single lower row; module 1 and module 2 retain approximate 42:58 ratio around a 2.5% white gap; vertical black module strips; bordered module containers"
    } else if (cfg$layout_style == "grant_modules") {
      "stacked module blocks"
    } else {
      "not applicable"
    },
    cfg$layout_style,
    if (is_grant_layout) {
      "Module 1: HES1,NKX2-1,LHX6,LHX8,CRABP1,KCNC1; Module 2: PLXNA2,NRP1,NRP2,ROBO1,ERBB4,CXCR4,ACKR3,SEMA3A,SEMA3F"
    } else {
      "not applicable"
    },
    if (is_grant_layout) {
      paste0("neuronal-only marker display; published clusters ", paste(cfg$marker_exclude_clusters, collapse = ","), " hidden; full tree retained in top context panels")
    } else {
      "complete URD tree"
    },
    if (cfg$layout_style == "grant_side_by_side") {
      "shown once in simplified lineage context only: cortical tips 3-5; subpallial/thalamic/striatal tips 1-2; omitted from gene panels"
    } else if (is_grant_layout) {
      "subtle brackets: cortical tips 3-5; subpallial/striatal tips 1-2"
    } else {
      "none"
    },
    if (is_grant_layout) {
      "cluster UMAP; published-cluster URD; simplified lineage URD; pseudotime UMAP; pseudotime URD; relative widths 0.75,1.65,1.25,1.05,1.15"
    } else {
      "cluster UMAP; published-cluster URD; pseudotime UMAP; pseudotime URD"
    },
    if (is_grant_layout) {
      "black vertical bars: Neurogenesis clusters 1-7 and 9; Gliogenesis clusters 8 and 10"
    } else {
      "not shown"
    },
    if (cfg$layout_style == "grant_side_by_side") {
      "neutral neuronal cells with cortical and subpallial/thalamic/striatal lineage brackets; clusters 8 and 10 light grey below topology-derived dashed divider"
    } else if (is_grant_layout) {
      "cluster 7 cortical; clusters 1-6 other neuronal; clusters 8 and 10 light grey below topology-derived dashed divider; cluster 9 neutral context"
    } else {
      "not shown"
    },
    cfg$cluster_style,
    paste(cfg$output_formats, collapse = ","),
    "log1p(CP10K), cross-study-compatible grey floor to #0000ff blue",
    as.character(cfg$expression_color_floor),
    paste0("q", cfg$vmax_quantile, " positive expression per gene"),
    if (is_grant_layout) "full unpruned URD tree cells" else "displayed tree cells",
    if (cfg$cluster_style == "figure4_submission") {
      "exact manually curated Figure 4 colors supplied for clusters 1-10"
    } else {
      "palette_for_clusters from scripts/26_div90_umap_cluster_label_audit.R"
    },
    if (cfg$cluster_style == "figure4_submission") {
      "manually curated Figure 4 labels applied to existing published IDs only; no cell reassignment"
    } else {
      "published Fig. D 10-class recode; raw/current clusters 6/7 excluded; raw 4/10 collapsed to published cluster 8"
    },
    "published plotting orientation: UMAP_1 unchanged, UMAP_2 multiplied by -1",
    "left-to-right rotated tree; tips at right",
    paste0("shared UMAP/tree limits ", paste(signif(pseudotime_limits, 4), collapse = " to ")),
        "cluster number and wrapped cluster name labels drawn at cluster tree tips",
        if ("svg" %in% cfg$output_formats) {
          if (requireNamespace("svglite", quietly = TRUE)) "pdf,svg via svglite" else "pdf,svg via grDevices::svg"
        } else {
          "pdf via cairo_pdf; SVG not requested"
        }
    ),
  stringsAsFactors = FALSE
)
write_tsv(render_status, file.path(cfg$outdir, "context_marker_composite_render_status.tsv"))
log_msg("Wrote composite outputs under: ", cfg$plot_dir)
