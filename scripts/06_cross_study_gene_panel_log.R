#!/usr/bin/env Rscript

# Build a cross-study multi-gene UMAP panel from existing Seurat objects.
# Expression values are always log1p-transformed.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

PNG_DPI <- 600
PALETTE_GREY_PURPLE <- c("#d9d9d9", "#b39ddb", "#6a1b9a")

log_msg <- function(...) {
  msg <- paste0(..., collapse = " ")
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
  flush.console()
}

trim_trailing_slash <- function(x) sub("/+$", "", x)
normalize_abs <- function(path, must_work = TRUE) normalizePath(path, winslash = "/", mustWork = must_work)

canonical_gene <- function(gene) {
  toupper(gsub("\\.", "-", trimws(as.character(gene))))
}

gene_aliases <- function(gene) {
  g <- trimws(as.character(gene))
  if (!nzchar(g)) return(character(0))
  g <- toupper(g)
  unique(c(g, gsub("\\.", "-", g), gsub("-", ".", g)))
}

parse_args <- function(args) {
  out <- list(
    config = NULL,
    `project-root` = NULL,
    `run-label` = NULL,
    genes = NULL,
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
      "  Rscript scripts/06_cross_study_gene_panel_log.R --config <config.R> [--project-root <path>] [--run-label <label>] [--genes <g1,g2,...>]",
      "",
      "Config fields:",
      "  project_root, run_label, optional genes, studies(data.frame)",
      "",
      "Outputs:",
      "  PROJECT_ROOT/results/<run_label>/plots/gene_panel_cross_study_log.png",
      "  PROJECT_ROOT/results/<run_label>/plots/gene_panel_cross_study_status.tsv",
      sep = "\n"
    )
  )
}

parse_gene_vector <- function(value) {
  if (is.null(value)) return(character(0))
  if (length(value) == 1 && is.character(value) && grepl(",", value, fixed = TRUE)) {
    value <- strsplit(value, ",", fixed = TRUE)[[1]]
  }
  value <- trimws(as.character(value))
  value <- value[nzchar(value)]
  unique(vapply(value, canonical_gene, character(1)))
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

  cfg$genes <- parse_gene_vector(cfg$genes)
  cfg
}

resolve_path <- function(path, project_root) {
  if (startsWith(path, "/")) return(normalize_abs(path, must_work = FALSE))
  normalize_abs(file.path(project_root, path), must_work = FALSE)
}

resolve_gene_in_assay <- function(obj, assay, gene_request) {
  features <- tryCatch(rownames(obj[[assay]]), error = function(e) character(0))
  if (length(features) == 0) return(list(found = FALSE, resolved = NA_character_))
  features_up <- toupper(features)
  aliases <- gene_aliases(gene_request)
  for (a in aliases) {
    idx <- match(a, features_up)
    if (!is.na(idx)) {
      return(list(found = TRUE, resolved = features[[idx]], alias = a))
    }
  }
  list(found = FALSE, resolved = NA_character_)
}

fetch_gene_expression_log <- function(obj, assay, gene_request, cells, preferred_slot = "data") {
  old_assay <- tryCatch(as.character(DefaultAssay(obj)), error = function(e) "")
  on.exit({
    if (nzchar(old_assay)) {
      tryCatch(DefaultAssay(obj) <- old_assay, error = function(e) NULL)
    }
  }, add = TRUE)
  tryCatch(DefaultAssay(obj) <- assay, error = function(e) NULL)

  resolved <- resolve_gene_in_assay(obj = obj, assay = assay, gene_request = gene_request)
  if (!isTRUE(resolved$found)) {
    return(list(status = "missing_gene", resolved_gene = NA_character_, expression_slot = NA_character_, expr = numeric(0)))
  }

  slots <- unique(c(preferred_slot, "data", "counts"))
  for (slot_name in slots) {
    fetched <- tryCatch(
      Seurat::FetchData(obj, vars = resolved$resolved, cells = cells, slot = slot_name),
      error = function(e) NULL
    )
    if (is.null(fetched) || nrow(fetched) == 0 || !(resolved$resolved %in% colnames(fetched))) next

    common_cells <- intersect(cells, rownames(fetched))
    if (length(common_cells) == 0) next

    vals <- as.numeric(fetched[common_cells, resolved$resolved, drop = TRUE])
    vals <- log1p(pmax(vals, 0))
    names(vals) <- common_cells
    return(list(status = "ok", resolved_gene = resolved$resolved, expression_slot = slot_name, expr = vals))
  }

  list(status = "missing_assay_data", resolved_gene = resolved$resolved, expression_slot = NA_character_, expr = numeric(0))
}

prepare_study_base <- function(study_row, project_root) {
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
    n_cells = NA_integer_,
    coords = NULL,
    obj = NULL
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
  info$n_cells <- nrow(coords)
  info$coords <- coords
  info$obj <- obj
  info
}

build_plot_data <- function(studies_info, genes) {
  point_chunks <- list()
  placeholder_chunks <- list()
  status_chunks <- list()
  pooled_expr <- numeric()

  for (study_info in studies_info) {
    if (!identical(study_info$status, "ok")) {
      for (g in genes) {
        status_chunks[[length(status_chunks) + 1]] <- data.frame(
          study_id = study_info$study_id,
          study_label = study_info$study_label,
          gene_requested = g,
          gene_resolved = "",
          status = study_info$status,
          reason = study_info$reason,
          object_path = study_info$object_path_resolved,
          reduction = study_info$reduction,
          assay = study_info$assay,
          expression_slot_requested = study_info$preferred_slot,
          expression_slot_selected = "",
          n_cells = study_info$n_cells,
          stringsAsFactors = FALSE
        )
        placeholder_chunks[[length(placeholder_chunks) + 1]] <- data.frame(
          study_label = study_info$study_label,
          gene_label = g,
          UMAP_1 = 0,
          UMAP_2 = 0,
          reason = study_info$reason,
          stringsAsFactors = FALSE
        )
      }
      next
    }

    coords <- study_info$coords
    cells <- coords$cell_id
    for (g in genes) {
      expr_res <- fetch_gene_expression_log(
        obj = study_info$obj,
        assay = study_info$assay,
        gene_request = g,
        cells = cells,
        preferred_slot = study_info$preferred_slot
      )

      st <- if (identical(expr_res$status, "ok")) "ok" else expr_res$status
      rsn <- if (identical(st, "ok")) "" else if (identical(st, "missing_gene")) paste0("Gene not found: ", g) else "No overlapping expression cells"

      status_chunks[[length(status_chunks) + 1]] <- data.frame(
        study_id = study_info$study_id,
        study_label = study_info$study_label,
        gene_requested = g,
        gene_resolved = ifelse(is.na(expr_res$resolved_gene), "", expr_res$resolved_gene),
        status = st,
        reason = rsn,
        object_path = study_info$object_path_resolved,
        reduction = study_info$reduction,
        assay = study_info$assay,
        expression_slot_requested = study_info$preferred_slot,
        expression_slot_selected = ifelse(is.na(expr_res$expression_slot), "", expr_res$expression_slot),
        n_cells = study_info$n_cells,
        stringsAsFactors = FALSE
      )

      if (!identical(st, "ok")) {
        placeholder_chunks[[length(placeholder_chunks) + 1]] <- data.frame(
          study_label = study_info$study_label,
          gene_label = g,
          UMAP_1 = 0,
          UMAP_2 = 0,
          reason = rsn,
          stringsAsFactors = FALSE
        )
        next
      }

      common_cells <- intersect(cells, names(expr_res$expr))
      if (length(common_cells) == 0) {
        placeholder_chunks[[length(placeholder_chunks) + 1]] <- data.frame(
          study_label = study_info$study_label,
          gene_label = g,
          UMAP_1 = 0,
          UMAP_2 = 0,
          reason = "No overlapping expression cells",
          stringsAsFactors = FALSE
        )
        next
      }

      coords_sub <- coords[match(common_cells, coords$cell_id), , drop = FALSE]
      expr_vals <- as.numeric(expr_res$expr[common_cells])
      chunk <- data.frame(
        study_label = study_info$study_label,
        gene_label = g,
        UMAP_1 = coords_sub$UMAP_1,
        UMAP_2 = coords_sub$UMAP_2,
        expr = expr_vals,
        stringsAsFactors = FALSE
      )
      point_chunks[[length(point_chunks) + 1]] <- chunk
      pooled_expr <- c(pooled_expr, expr_vals[is.finite(expr_vals)])
    }
  }

  point_df <- if (length(point_chunks) > 0) do.call(rbind, point_chunks) else data.frame(
    study_label = character(), gene_label = character(), UMAP_1 = numeric(), UMAP_2 = numeric(), expr = numeric(), stringsAsFactors = FALSE
  )
  placeholder_df <- if (length(placeholder_chunks) > 0) do.call(rbind, placeholder_chunks) else data.frame(
    study_label = character(), gene_label = character(), UMAP_1 = numeric(), UMAP_2 = numeric(), reason = character(), stringsAsFactors = FALSE
  )
  status_df <- if (length(status_chunks) > 0) do.call(rbind, status_chunks) else data.frame(
    study_id = character(), study_label = character(), gene_requested = character(), gene_resolved = character(),
    status = character(), reason = character(), object_path = character(), reduction = character(), assay = character(),
    expression_slot_requested = character(), expression_slot_selected = character(), n_cells = integer(), stringsAsFactors = FALSE
  )

  list(point_df = point_df, placeholder_df = placeholder_df, status_df = status_df, pooled_expr = pooled_expr)
}

build_plot <- function(plot_data, study_labels, genes) {
  point_df <- plot_data$point_df
  placeholder_df <- plot_data$placeholder_df
  pooled_expr <- plot_data$pooled_expr

  # Rescale coordinates per study for plotting so each study uses panel space.
  # This keeps UMAP aspect ratio within each study while avoiding cross-study squishing.
  if (nrow(point_df) > 0) {
    for (study in unique(as.character(point_df$study_label))) {
      idx <- which(as.character(point_df$study_label) == study)
      if (length(idx) == 0) next
      xr <- range(point_df$UMAP_1[idx], finite = TRUE, na.rm = TRUE)
      yr <- range(point_df$UMAP_2[idx], finite = TRUE, na.rm = TRUE)
      span <- max(diff(xr), diff(yr))
      if (!is.finite(span) || span <= 0) span <- 1
      xmid <- mean(xr)
      ymid <- mean(yr)
      point_df$UMAP_1[idx] <- (point_df$UMAP_1[idx] - xmid) / span
      point_df$UMAP_2[idx] <- (point_df$UMAP_2[idx] - ymid) / span
    }
  }

  point_df$study_label <- factor(point_df$study_label, levels = study_labels)
  point_df$gene_label <- factor(point_df$gene_label, levels = genes)
  placeholder_df$study_label <- factor(placeholder_df$study_label, levels = study_labels)
  placeholder_df$gene_label <- factor(placeholder_df$gene_label, levels = genes)

  limits <- if (length(pooled_expr) > 0) range(pooled_expr, na.rm = TRUE, finite = TRUE) else c(0, 1)
  if (!all(is.finite(limits))) limits <- c(0, 1)
  if (limits[1] == limits[2]) limits[2] <- limits[2] + 1e-6

  facet_seed <- expand.grid(
    study_label = factor(study_labels, levels = study_labels),
    gene_label = factor(genes, levels = genes),
    stringsAsFactors = FALSE
  )
  facet_seed$UMAP_1 <- 0
  facet_seed$UMAP_2 <- 0

  legend_seed <- data.frame(
    study_label = factor(study_labels[[1]], levels = study_labels),
    gene_label = factor(genes[[1]], levels = genes),
    UMAP_1 = c(0, 0), UMAP_2 = c(0, 0), expr = limits,
    stringsAsFactors = FALSE
  )

  p <- ggplot() +
    geom_blank(data = facet_seed, aes(x = UMAP_1, y = UMAP_2)) +
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
      size = 0.40, alpha = 0.85, stroke = 0, show.legend = FALSE
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
    facet_grid(rows = vars(study_label), cols = vars(gene_label), drop = FALSE, scales = "fixed") +
    coord_equal() +
    scale_color_gradientn(
      colours = PALETTE_GREY_PURPLE,
      limits = limits,
      oob = scales::squish,
      name = "Expression (log1p)"
    ) +
    guides(
      color = guide_colorbar(
        title.position = "left", title.hjust = 0.5, direction = "horizontal",
        barwidth = grid::unit(88, "pt"), barheight = grid::unit(8, "pt")
      )
    ) +
    labs(
      title = "Cross-study marker expression (log1p)",
      subtitle = "Color scale: grey -> purple; per-study UMAP display scaling",
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
      panel.spacing = grid::unit(4, "mm"),
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
  if (is.null(run_label) || !nzchar(run_label)) run_label <- cfg$run_label
  if (is.null(run_label) || !nzchar(run_label)) run_label <- "gene_panel_cross_study_v1"

  genes <- parse_gene_vector(args$genes)
  if (length(genes) == 0) genes <- cfg$genes
  if (length(genes) == 0) stop("No genes provided. Set config$genes or --genes.", call. = FALSE)

  log_msg("Run label: ", run_label)
  log_msg("Project root: ", project_root)
  log_msg("Genes: ", paste(genes, collapse = ","))

  studies_info <- vector("list", nrow(cfg$studies))
  for (i in seq_len(nrow(cfg$studies))) {
    row <- cfg$studies[i, , drop = FALSE]
    log_msg("Study ", i, "/", nrow(cfg$studies), ": ", row$study_label[[1]], " (", row$study_id[[1]], ")")
    studies_info[[i]] <- prepare_study_base(row, project_root = project_root)
    log_msg("Study ", row$study_label[[1]], " status=", studies_info[[i]]$status,
            ifelse(nzchar(studies_info[[i]]$reason), paste0(" (", studies_info[[i]]$reason, ")"), ""))
  }

  study_labels <- vapply(studies_info, function(x) x$study_label, character(1))
  plot_data <- build_plot_data(studies_info = studies_info, genes = genes)
  p <- build_plot(plot_data = plot_data, study_labels = study_labels, genes = genes)

  out_dir <- normalize_abs(file.path(project_root, "results", run_label, "plots"), must_work = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  png_path <- file.path(out_dir, "gene_panel_cross_study_log.png")
  status_path <- file.path(out_dir, "gene_panel_cross_study_status.tsv")

  utils::write.table(plot_data$status_df, file = status_path, sep = "\t", quote = FALSE, row.names = FALSE)
  log_msg("Wrote status table: ", status_path)

  fig_width <- max(10, length(genes) * 3.2 + 1.8)
  fig_height <- max(4.8, length(study_labels) * 2.8 + 1.2)
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
