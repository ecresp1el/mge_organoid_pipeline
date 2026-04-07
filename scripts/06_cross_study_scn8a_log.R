#!/usr/bin/env Rscript

# Build a cross-study SCN8A-only UMAP panel from existing Seurat objects.
# The plotted expression is always log-transformed via log1p(), regardless
# of whether values are fetched from data/counts slots.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

TARGET_GENE <- "SCN8A"
PNG_DPI <- 600

log_msg <- function(...) {
  msg <- paste0(..., collapse = " ")
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
  flush.console()
}

fmt_count <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x)) return("NA")
  format(round(as.numeric(x), digits = 0), trim = TRUE, scientific = FALSE, big.mark = ",")
}

trim_trailing_slash <- function(x) sub("/+$", "", x)

normalize_abs <- function(path, must_work = TRUE) {
  normalizePath(path, winslash = "/", mustWork = must_work)
}

parse_args <- function(args) {
  out <- list(
    config = NULL,
    `project-root` = NULL,
    `run-label` = NULL,
    help = FALSE
  )
  i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (a == "--help" || a == "-h") {
      out$help <- TRUE
      i <- i + 1L
      next
    }
    if (!startsWith(a, "--")) stop("Unknown argument: ", a, call. = FALSE)
    key <- substring(a, 3L)
    if (!(key %in% names(out))) stop("Unknown argument: ", a, call. = FALSE)
    if (identical(key, "help")) {
      out$help <- TRUE
      i <- i + 1L
      next
    }
    if (i == length(args)) stop("Missing value for argument: ", a, call. = FALSE)
    out[[key]] <- args[[i + 1L]]
    i <- i + 2L
  }
  out
}

print_usage <- function() {
  cat(
    paste(
      "Usage:",
      "  Rscript scripts/06_cross_study_scn8a_log.R --config <config.R> [--project-root <path>] [--run-label <label>]",
      "",
      "Outputs:",
      "  PROJECT_ROOT/results/<run_label>/plots/scn8a_cross_study_log.png",
      "  PROJECT_ROOT/results/<run_label>/plots/scn8a_cross_study_study_status.tsv",
      sep = "\n"
    )
  )
}

read_config <- function(path) {
  cfg <- dget(path)
  if (!is.list(cfg)) stop("Config must be a list readable by dget().", call. = FALSE)
  if (is.null(cfg$studies) || !is.data.frame(cfg$studies)) {
    stop("Config must contain studies as a data.frame.", call. = FALSE)
  }
  req <- c("study_id", "study_label", "object_path")
  missing_cols <- setdiff(req, colnames(cfg$studies))
  if (length(missing_cols) > 0) {
    stop("studies is missing required column(s): ", paste(missing_cols, collapse = ","), call. = FALSE)
  }

  s <- cfg$studies
  if (!("reduction" %in% colnames(s))) s$reduction <- "umap"
  if (!("assay" %in% colnames(s))) s$assay <- "RNA"
  if (!("expression_slot" %in% colnames(s))) s$expression_slot <- "data"
  s$reduction[is.na(s$reduction) | !nzchar(s$reduction)] <- "umap"
  s$assay[is.na(s$assay) | !nzchar(s$assay)] <- "RNA"
  s$expression_slot[is.na(s$expression_slot) | !nzchar(s$expression_slot)] <- "data"
  cfg$studies <- s
  cfg
}

plot_study_order <- function(studies_info) {
  labels <- vapply(studies_info, function(x) x$study_label, character(1))
  ids <- tolower(vapply(studies_info, function(x) x$study_id, character(1)))
  labels_norm <- tolower(labels)
  varela_idx <- which(grepl("varela", ids, fixed = TRUE) | grepl("varela", labels_norm, fixed = TRUE))
  if (length(varela_idx) == 0) return(labels)
  varela_idx <- unique(varela_idx)
  c(labels[varela_idx], labels[setdiff(seq_along(labels), varela_idx)])
}

resolve_path <- function(path, project_root) {
  if (startsWith(path, "/")) return(normalize_abs(path, must_work = FALSE))
  normalize_abs(file.path(project_root, path), must_work = FALSE)
}

fetch_gene_expression_log <- function(obj, assay, gene, cells, preferred_slot = "data") {
  old_assay <- tryCatch(as.character(DefaultAssay(obj)), error = function(e) "")
  on.exit({
    if (nzchar(old_assay)) {
      tryCatch(DefaultAssay(obj) <- old_assay, error = function(e) NULL)
    }
  }, add = TRUE)
  tryCatch(DefaultAssay(obj) <- assay, error = function(e) NULL)

  slots <- unique(c(preferred_slot, "data", "counts"))
  for (slot_name in slots) {
    fetched <- tryCatch(
      Seurat::FetchData(obj, vars = gene, cells = cells, slot = slot_name),
      error = function(e) NULL
    )
    if (is.null(fetched) || nrow(fetched) == 0 || !(gene %in% colnames(fetched))) next

    common_cells <- intersect(cells, rownames(fetched))
    if (length(common_cells) == 0) next

    vals <- as.numeric(fetched[common_cells, gene, drop = TRUE])
    vals <- log1p(pmax(vals, 0))
    names(vals) <- common_cells
    return(list(status = "ok", expression_slot = slot_name, expr = vals))
  }

  list(status = "missing_gene", expression_slot = NA_character_, expr = numeric(0))
}

prepare_study <- function(study_row, project_root) {
  info <- list(
    study_id = study_row$study_id[[1]],
    study_label = study_row$study_label[[1]],
    status = "ok",
    reason = "",
    object_path = study_row$object_path[[1]],
    object_path_resolved = "",
    reduction = study_row$reduction[[1]],
    assay = study_row$assay[[1]],
    preferred_slot = study_row$expression_slot[[1]],
    expression_slot = NA_character_,
    n_cells = NA_integer_,
    coords = NULL,
    expr = NULL
  )

  obj_path <- resolve_path(info$object_path, project_root)
  info$object_path_resolved <- obj_path
  if (!file.exists(obj_path)) {
    info$status <- "missing_object"
    info$reason <- "Missing object"
    return(info)
  }

  obj <- tryCatch(readRDS(obj_path), error = function(e) e)
  if (inherits(obj, "error")) {
    info$status <- "read_error"
    info$reason <- conditionMessage(obj)
    return(info)
  }

  reds <- tryCatch(Reductions(obj), error = function(e) character(0))
  if (!(info$reduction %in% reds)) {
    info$status <- "missing_umap"
    info$reason <- paste0("Missing reduction: ", info$reduction)
    return(info)
  }

  coords_raw <- tryCatch(Embeddings(obj, reduction = info$reduction), error = function(e) e)
  if (inherits(coords_raw, "error") || is.null(coords_raw) || nrow(coords_raw) == 0 || ncol(coords_raw) < 2) {
    info$status <- "missing_umap"
    info$reason <- "Invalid reduction coordinates"
    return(info)
  }

  coords <- as.data.frame(coords_raw[, 1:2, drop = FALSE], stringsAsFactors = FALSE)
  colnames(coords) <- c("UMAP_1", "UMAP_2")
  coords$cell_id <- rownames(coords)
  cells <- coords$cell_id
  info$n_cells <- length(cells)

  expr_res <- fetch_gene_expression_log(
    obj = obj,
    assay = info$assay,
    gene = TARGET_GENE,
    cells = cells,
    preferred_slot = info$preferred_slot
  )
  info$expression_slot <- expr_res$expression_slot
  if (!identical(expr_res$status, "ok")) {
    info$status <- "missing_gene"
    info$reason <- paste0("Gene not found: ", TARGET_GENE)
    return(info)
  }

  common_cells <- intersect(cells, names(expr_res$expr))
  if (length(common_cells) == 0) {
    info$status <- "missing_assay_data"
    info$reason <- "No overlapping expression cells"
    return(info)
  }

  coords <- coords[match(common_cells, coords$cell_id), c("UMAP_1", "UMAP_2", "cell_id"), drop = FALSE]
  rownames(coords) <- coords$cell_id
  info$coords <- coords
  info$expr <- expr_res$expr[common_cells]
  info$status <- "ok"
  info$reason <- ""
  info
}

build_status_table <- function(studies_info) {
  data.frame(
    study_id = vapply(studies_info, function(x) x$study_id, character(1)),
    study_label = vapply(studies_info, function(x) x$study_label, character(1)),
    status = vapply(studies_info, function(x) x$status, character(1)),
    reason = vapply(studies_info, function(x) x$reason, character(1)),
    object_path = vapply(studies_info, function(x) x$object_path_resolved, character(1)),
    reduction = vapply(studies_info, function(x) x$reduction, character(1)),
    assay = vapply(studies_info, function(x) x$assay, character(1)),
    expression_slot_requested = vapply(studies_info, function(x) x$preferred_slot, character(1)),
    expression_slot_selected = vapply(studies_info, function(x) ifelse(is.na(x$expression_slot), "", x$expression_slot), character(1)),
    n_cells = vapply(studies_info, function(x) ifelse(is.na(x$n_cells), NA_integer_, as.integer(x$n_cells)), integer(1)),
    stringsAsFactors = FALSE
  )
}

build_plot <- function(studies_info, ordered_labels) {
  point_chunks <- list()
  placeholder_chunks <- list()
  pooled_expr <- numeric()

  for (study_info in studies_info) {
    if (identical(study_info$status, "ok")) {
      chunk <- data.frame(
        study_label = study_info$study_label,
        UMAP_1 = study_info$coords$UMAP_1,
        UMAP_2 = study_info$coords$UMAP_2,
        expr = as.numeric(study_info$expr),
        stringsAsFactors = FALSE
      )
      point_chunks[[length(point_chunks) + 1]] <- chunk
      pooled_expr <- c(pooled_expr, chunk$expr[is.finite(chunk$expr)])
    } else {
      placeholder_chunks[[length(placeholder_chunks) + 1]] <- data.frame(
        study_label = study_info$study_label,
        UMAP_1 = 0,
        UMAP_2 = 0,
        reason = study_info$reason,
        stringsAsFactors = FALSE
      )
    }
  }

  point_df <- if (length(point_chunks) > 0) do.call(rbind, point_chunks) else data.frame(
    study_label = character(), UMAP_1 = numeric(), UMAP_2 = numeric(), expr = numeric(), stringsAsFactors = FALSE
  )
  placeholder_df <- if (length(placeholder_chunks) > 0) do.call(rbind, placeholder_chunks) else data.frame(
    study_label = character(), UMAP_1 = numeric(), UMAP_2 = numeric(), reason = character(), stringsAsFactors = FALSE
  )

  point_df$study_label <- factor(point_df$study_label, levels = ordered_labels)
  placeholder_df$study_label <- factor(placeholder_df$study_label, levels = ordered_labels)

  limits <- if (length(pooled_expr) > 0) range(pooled_expr, na.rm = TRUE, finite = TRUE) else c(0, 1)
  if (!all(is.finite(limits))) limits <- c(0, 1)
  if (limits[1] == limits[2]) limits[2] <- limits[2] + 1e-6

  facets_seed <- data.frame(
    study_label = factor(ordered_labels, levels = ordered_labels),
    UMAP_1 = 0, UMAP_2 = 0, stringsAsFactors = FALSE
  )
  legend_seed <- data.frame(
    study_label = factor(ordered_labels[[1]], levels = ordered_labels),
    UMAP_1 = c(0, 0), UMAP_2 = c(0, 0), expr = limits, stringsAsFactors = FALSE
  )

  p <- ggplot() +
    geom_blank(data = facets_seed, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(
      data = legend_seed,
      aes(x = UMAP_1, y = UMAP_2, color = expr),
      inherit.aes = FALSE,
      alpha = 0, size = 0.01, show.legend = TRUE
    )

  if (nrow(point_df) > 0) {
    p <- p + geom_point(
      data = point_df,
      aes(x = UMAP_1, y = UMAP_2, color = expr),
      size = 0.20, alpha = 0.85, stroke = 0, show.legend = FALSE
    )
  }

  if (nrow(placeholder_df) > 0) {
    p <- p + geom_text(
      data = placeholder_df,
      aes(x = UMAP_1, y = UMAP_2, label = reason),
      inherit.aes = FALSE,
      color = "grey35", size = 3.0, lineheight = 0.9
    )
  }

  p +
    facet_wrap(~study_label, nrow = 1, drop = FALSE, scales = "fixed") +
    coord_equal() +
    scale_color_gradientn(
      colours = c("#d9d9d9", "#6baed6", "#08306b"),
      limits = limits,
      oob = scales::squish,
      name = paste0(TARGET_GENE, " (log1p)")
    ) +
    guides(
      color = guide_colorbar(
        title.position = "left", title.hjust = 0.5, direction = "horizontal",
        barwidth = grid::unit(88, "pt"), barheight = grid::unit(8, "pt")
      )
    ) +
    labs(
      title = paste0(TARGET_GENE, " expression across studies (always log1p)"),
      subtitle = "Single-gene cross-study UMAP panel",
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.3),
      strip.background = element_rect(fill = "grey96", color = "grey80", linewidth = 0.3),
      strip.text = element_text(size = 9, face = "bold"),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(size = 12, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 9, hjust = 0),
      panel.spacing = grid::unit(5, "mm"),
      legend.position = "bottom"
    )
}

main <- function(cli_args = commandArgs(trailingOnly = TRUE)) {
  args <- parse_args(cli_args)
  if (isTRUE(args$help)) {
    print_usage()
    quit(save = "no", status = 0)
  }
  if (is.null(args$config) || !nzchar(args$config)) {
    print_usage()
    stop("--config is required.", call. = FALSE)
  }

  cfg <- read_config(args$config)
  project_root <- args[["project-root"]]
  if (is.null(project_root) || !nzchar(project_root)) project_root <- cfg$project_root
  if (!nzchar(project_root)) stop("PROJECT_ROOT not set.", call. = FALSE)
  project_root <- normalize_abs(trim_trailing_slash(project_root), must_work = FALSE)
  if (!dir.exists(project_root)) stop("PROJECT_ROOT does not exist: ", project_root, call. = FALSE)

  run_label <- args[["run-label"]]
  if (is.null(run_label) || !nzchar(run_label)) run_label <- "panel_b_scn8a_log_v1"

  log_msg("Run label: ", run_label)
  log_msg("Project root: ", project_root)
  log_msg("Target gene: ", TARGET_GENE)

  studies_info <- vector("list", nrow(cfg$studies))
  for (i in seq_len(nrow(cfg$studies))) {
    row <- cfg$studies[i, , drop = FALSE]
    log_msg("Study ", i, "/", nrow(cfg$studies), ": ", row$study_label[[1]], " (", row$study_id[[1]], ")")
    studies_info[[i]] <- prepare_study(row, project_root = project_root)
    log_msg("Study ", row$study_label[[1]], " status=", studies_info[[i]]$status,
            ifelse(nzchar(studies_info[[i]]$reason), paste0(" (", studies_info[[i]]$reason, ")"), ""))
  }

  ordered_labels <- plot_study_order(studies_info)
  ordered_info <- studies_info[match(ordered_labels, vapply(studies_info, function(x) x$study_label, character(1)))]
  p <- build_plot(ordered_info, ordered_labels = ordered_labels)

  out_dir <- normalize_abs(file.path(project_root, "results", run_label, "plots"), must_work = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  png_path <- file.path(out_dir, "scn8a_cross_study_log.png")
  status_path <- file.path(out_dir, "scn8a_cross_study_study_status.tsv")

  status_tbl <- build_status_table(ordered_info)
  utils::write.table(status_tbl, file = status_path, sep = "\t", quote = FALSE, row.names = FALSE)
  log_msg("Wrote study status table: ", status_path)

  fig_width <- max(16, length(ordered_labels) * 2.6 + 2.0)
  fig_height <- 4.8
  log_msg("Writing PNG: ", png_path)
  ggsave(
    filename = png_path, plot = p, width = fig_width, height = fig_height,
    units = "in", dpi = PNG_DPI, bg = "white", limitsize = FALSE
  )
  log_msg("Done.")
}

if (sys.nframe() == 0) {
  main()
}

