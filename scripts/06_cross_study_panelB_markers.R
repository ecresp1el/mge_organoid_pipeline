#!/usr/bin/env Rscript

# Assemble a single cross-study marker-expression figure (Panel B only)
# from precomputed Seurat objects. This script performs no recomputation.
#
# Scope and assumptions:
# - No analysis is performed (no normalization/scaling/PCA/UMAP/clustering/integration).
# - Study identities, labels, ordering, object paths, assays, and reductions come from config.
# - Fixed marker order is enforced to match the reference panel.
# - Every study object path must be explicitly provided and must reside under PROJECT_ROOT.
# - Missing object/reduction/assay/gene conditions are rendered as same-size placeholders.
#
# Inputs:
# - Config file (`--config`) containing:
#   - run_label
#   - optional project_root
#   - studies data.frame with explicit paths and plotting metadata
# - Optional CLI overrides: --project-root and --run-label
#
# Outputs:
# - PROJECT_ROOT/results/<run_label>/plots/panel_b_cross_study_markers.pdf
# - PROJECT_ROOT/results/<run_label>/plots/panel_b_cross_study_markers.svg
# - stdout audit of missing studies/genes/components

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

ON_TARGET_GENES <- c("DCX", "GAD2", "DLX5", "LHX6", "MAF", "SST", "LHX8", "SP8")
OFF_TARGET_GENES <- c("PAX6", "NEUROD2", "ISL1", "ACHF")
GENE_ORDER <- c(ON_TARGET_GENES, OFF_TARGET_GENES)

log_msg <- function(...) {
  msg <- paste0(..., collapse = " ")
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
  flush.console()
}

print_usage <- function() {
  cat(
    paste(
      "Usage:",
      "  Rscript scripts/06_cross_study_panelB_markers.R --config <config.R> [--project-root <path>] [--run-label <label>] [--retain-seurat <true|false>] [--export-global <true|false>]",
      "",
      "Config file format:",
      "  A single R object (list) readable by dget(), with fields:",
      "    project_root (optional if --project-root or PROJECT_ROOT is set)",
      "    run_label",
      "    studies (data.frame with columns: study_id, study_label, object_path, reduction, assay,",
      "             and optional expression_slot)",
      "",
      "Optional flags:",
      "  --retain-seurat true|false    keep full Seurat objects in returned study list (interactive debugging)",
      "  --export-global true|false    export run objects to .GlobalEnv (panel_b_result, panel_b_studies, panel_b_rows, panel_b_issues)",
      "  --show-progress true|false    print each gene-row plot during assembly (interactive use)",
      "",
      "Outputs:",
      "  PROJECT_ROOT/results/<run_label>/plots/panel_b_cross_study_markers.{pdf,svg}",
      sep = "\n"
    )
  )
}

parse_args <- function(args) {
  # Minimal flag parser to avoid external dependencies.
  out <- list()
  i <- 1
  while (i <= length(args)) {
    token <- args[[i]]
    if (token %in% c("-h", "--help")) {
      out$help <- TRUE
      i <- i + 1
      next
    }
    if (!startsWith(token, "--")) {
      stop("Unexpected argument: ", token, call. = FALSE)
    }
    token <- sub("^--", "", token)
    if (grepl("=", token, fixed = TRUE)) {
      kv <- strsplit(token, "=", fixed = TRUE)[[1]]
      key <- kv[[1]]
      value <- if (length(kv) > 1) paste(kv[-1], collapse = "=") else ""
      out[[key]] <- value
      i <- i + 1
      next
    }
    key <- token
    if (i == length(args)) {
      stop("Missing value for --", key, call. = FALSE)
    }
    out[[key]] <- args[[i + 1]]
    i <- i + 2
  }
  out
}

parse_bool_flag <- function(x, default = FALSE) {
  if (is.null(x) || !nzchar(as.character(x))) return(default)
  val <- tolower(as.character(x))
  if (val %in% c("1", "true", "t", "yes", "y")) return(TRUE)
  if (val %in% c("0", "false", "f", "no", "n")) return(FALSE)
  stop("Invalid boolean flag value: ", x, call. = FALSE)
}

trim_trailing_slash <- function(x) {
  sub("/+$", "", x)
}

normalize_abs <- function(path, must_work = FALSE) {
  normalizePath(path, winslash = "/", mustWork = must_work)
}

is_subpath <- function(path, root) {
  path_norm <- normalize_abs(path, must_work = FALSE)
  root_norm <- normalize_abs(root, must_work = FALSE)
  if (identical(path_norm, root_norm)) return(TRUE)
  startsWith(path_norm, paste0(root_norm, "/"))
}

resolve_under_project_root <- function(path_value, project_root) {
  # Enforce explicit, non-globbed study paths under PROJECT_ROOT only.
  if (!is.character(path_value) || length(path_value) != 1 || !nzchar(path_value)) {
    stop("Invalid path in config: must be non-empty string", call. = FALSE)
  }
  if (grepl("[*?\\[]", path_value)) {
    stop("Path globbing is not allowed in config path: ", path_value, call. = FALSE)
  }
  candidate <- if (startsWith(path_value, "/")) {
    path_value
  } else {
    file.path(project_root, path_value)
  }
  candidate <- normalize_abs(candidate, must_work = FALSE)
  if (!is_subpath(candidate, project_root)) {
    stop(
      "All study object paths must be under PROJECT_ROOT. Invalid path: ",
      candidate,
      call. = FALSE
    )
  }
  candidate
}

read_config <- function(config_path) {
  # Strict schema checks keep study ordering/labels/config behavior deterministic.
  cfg <- tryCatch(
    dget(config_path),
    error = function(e) stop("Failed to read config via dget(): ", conditionMessage(e), call. = FALSE)
  )
  if (!is.list(cfg)) stop("Config must be a list.", call. = FALSE)
  if (is.null(cfg$run_label) || !nzchar(as.character(cfg$run_label))) {
    stop("Config missing run_label.", call. = FALSE)
  }
  if (is.null(cfg$studies) || !is.data.frame(cfg$studies)) {
    stop("Config missing studies data.frame.", call. = FALSE)
  }

  required_cols <- c("study_id", "study_label", "object_path", "reduction", "assay")
  missing_cols <- setdiff(required_cols, colnames(cfg$studies))
  if (length(missing_cols) > 0) {
    stop("studies is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  studies <- cfg$studies
  for (col in required_cols) {
    studies[[col]] <- as.character(studies[[col]])
    if (any(!nzchar(studies[[col]]))) {
      stop("studies$", col, " contains empty values.", call. = FALSE)
    }
  }

  if ("expression_slot" %in% colnames(studies)) {
    studies$expression_slot <- as.character(studies$expression_slot)
    studies$expression_slot[is.na(studies$expression_slot) | !nzchar(studies$expression_slot)] <- "data"
  } else {
    studies$expression_slot <- "data"
  }

  if (anyDuplicated(studies$study_id) > 0) {
    stop("studies$study_id must be unique.", call. = FALSE)
  }
  if (nrow(studies) == 0) {
    stop("studies must include at least one row.", call. = FALSE)
  }
  if (anyDuplicated(studies$study_label) > 0) {
    stop("studies$study_label must be unique for fixed panel layout.", call. = FALSE)
  }

  cfg$run_label <- as.character(cfg$run_label)
  if (!is.null(cfg$project_root)) cfg$project_root <- as.character(cfg$project_root)
  cfg$studies <- studies
  cfg
}

read_rds_any <- function(path) {
  # Robust reader for plain .rds and several gzip wrapping variants seen in GEO files.
  if (!grepl("\\.gz$", path, ignore.case = TRUE)) {
    return(readRDS(path))
  }

  attempts <- list(
    function() {
      con <- gzfile(path, open = "rb")
      on.exit(close(con), add = TRUE)
      readRDS(con)
    },
    function() {
      con <- gzcon(gzfile(path, open = "rb"))
      on.exit(close(con), add = TRUE)
      readRDS(con)
    },
    function() {
      con <- pipe(sprintf("gunzip -c %s", shQuote(path)), open = "rb")
      on.exit(close(con), add = TRUE)
      readRDS(con)
    },
    function() {
      con <- pipe(sprintf("gunzip -c %s | gunzip -c", shQuote(path)), open = "rb")
      on.exit(close(con), add = TRUE)
      readRDS(con)
    }
  )

  last_error <- NULL
  for (fn in attempts) {
    candidate <- tryCatch(
      fn(),
      error = function(e) {
        last_error <<- e
        NULL
      }
    )
    if (!is.null(candidate)) return(candidate)
  }
  stop("Unable to read gzipped RDS: ", path, " (", conditionMessage(last_error), ")", call. = FALSE)
}

choose_expression_slot <- function(obj, assay, preferred_slot) {
  # Prefer configured slot, then safe fallbacks commonly present in Seurat objects.
  candidates <- unique(c(preferred_slot, "data", "counts"))
  for (slot_name in candidates) {
    ok <- tryCatch(
      {
        mat <- suppressWarnings(GetAssayData(obj, assay = assay, slot = slot_name))
        !is.null(mat) && nrow(mat) > 0 && ncol(mat) > 0
      },
      error = function(e) FALSE
    )
    if (ok) return(slot_name)
  }
  NA_character_
}

prepare_study <- function(study_row, project_root, retain_seurat = FALSE) {
  # Per-study preprocessing:
  # - validate object/reduction/assay availability
  # - extract UMAP coords
  # - extract only target genes to reduce memory footprint
  info <- list(
    study_id = study_row$study_id,
    study_label = study_row$study_label,
    object_path = resolve_under_project_root(study_row$object_path, project_root),
    reduction = study_row$reduction,
    assay = study_row$assay,
    preferred_slot = study_row$expression_slot,
    expression_slot = NA_character_,
    status = "ok",
    reason = NA_character_,
    detail = NA_character_,
    seurat_obj = NULL,
    expr_sub = NULL,
    coords = NULL
  )

  if (!file.exists(info$object_path)) {
    info$status <- "missing_object"
    info$reason <- "Missing object"
    return(info)
  }

  obj <- tryCatch(read_rds_any(info$object_path), error = function(e) e)
  if (inherits(obj, "error")) {
    info$status <- "unreadable_object"
    info$reason <- "Unreadable object"
    info$detail <- conditionMessage(obj)
    return(info)
  }

  if (!inherits(obj, "Seurat")) {
    info$status <- "invalid_object"
    info$reason <- "Invalid object"
    info$detail <- paste("Class:", paste(class(obj), collapse = ","))
    return(info)
  }
  if (retain_seurat) {
    info$seurat_obj <- obj
  }

  if (!(info$assay %in% names(obj@assays))) {
    info$status <- "missing_assay"
    info$reason <- "Missing assay"
    info$detail <- paste0("assay=", info$assay)
    return(info)
  }

  if (!(info$reduction %in% names(obj@reductions))) {
    info$status <- "missing_umap"
    info$reason <- "Missing UMAP"
    info$detail <- paste0("reduction=", info$reduction)
    return(info)
  }

  coords <- tryCatch(Embeddings(obj, reduction = info$reduction), error = function(e) e)
  if (inherits(coords, "error") || is.null(dim(coords)) || ncol(coords) < 2 || nrow(coords) == 0) {
    info$status <- "missing_umap"
    info$reason <- "Missing UMAP"
    info$detail <- paste0("reduction=", info$reduction)
    return(info)
  }

  slot_name <- choose_expression_slot(obj, info$assay, info$preferred_slot)
  if (is.na(slot_name)) {
    info$status <- "missing_assay_data"
    info$reason <- "Missing assay data"
    info$detail <- paste0("assay=", info$assay)
    return(info)
  }

  coords <- as.data.frame(coords[, 1:2, drop = FALSE], stringsAsFactors = FALSE)
  colnames(coords) <- c("UMAP_1", "UMAP_2")
  coords$cell_id <- rownames(coords)

  mat <- tryCatch(
    suppressWarnings(GetAssayData(obj, assay = info$assay, slot = slot_name)),
    error = function(e) e
  )
  if (inherits(mat, "error") || is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
    info$status <- "missing_assay_data"
    info$reason <- "Missing assay data"
    info$detail <- paste0("assay=", info$assay, "; slot=", slot_name)
    return(info)
  }

  cells <- coords$cell_id
  if (!all(cells %in% colnames(mat))) {
    info$status <- "missing_assay_data"
    info$reason <- "Cell mismatch"
    info$detail <- "UMAP cells not all present in assay matrix."
    return(info)
  }

  gene_hits <- GENE_ORDER[GENE_ORDER %in% rownames(mat)]
  expr_sub <- if (length(gene_hits) > 0) {
    mat[gene_hits, cells, drop = FALSE]
  } else {
    NULL
  }

  if (retain_seurat) {
    rm(mat)
  } else {
    rm(mat, obj)
  }
  invisible(gc(verbose = FALSE))

  info$expression_slot <- slot_name
  info$expr_sub <- expr_sub
  info$coords <- coords
  info
}

extract_gene_values <- function(study_info, gene) {
  # Pull one gene vector from the pre-extracted matrix for the row builder.
  if (is.null(study_info$expr_sub) || !(gene %in% rownames(study_info$expr_sub))) {
    return(list(status = "missing_gene", expr = NULL, reason = "Gene not found"))
  }
  vals <- as.numeric(study_info$expr_sub[gene, , drop = TRUE])
  list(status = "ok", expr = vals, reason = NA_character_)
}

build_study_status_table <- function(studies_info) {
  chunks <- lapply(studies_info, function(study_info) {
    data.frame(
      study_id = study_info$study_id,
      study_label = study_info$study_label,
      status = study_info$status,
      reason = ifelse(is.na(study_info$reason), "", study_info$reason),
      detail = ifelse(is.na(study_info$detail), "", study_info$detail),
      reduction = study_info$reduction,
      assay = study_info$assay,
      expression_slot = ifelse(is.na(study_info$expression_slot), "", study_info$expression_slot),
      n_cells = if (!is.null(study_info$coords)) nrow(study_info$coords) else NA_integer_,
      n_marker_genes_available = if (!is.null(study_info$expr_sub)) nrow(study_info$expr_sub) else NA_integer_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, chunks)
}

build_gene_row <- function(gene, gene_group, studies_info, ordered_labels) {
  # Build one full row (one gene across all studies):
  # - collect expression points where available
  # - emit placeholder panels where unavailable
  # - compute shared row-wise limits for cross-study comparability
  point_chunks <- list()
  placeholder_chunks <- list()
  issue_chunks <- list()
  pooled_expr <- numeric()

  for (study_info in studies_info) {
    if (study_info$status != "ok") {
      placeholder_chunks[[length(placeholder_chunks) + 1]] <- data.frame(
        study_label = study_info$study_label,
        UMAP_1 = 0,
        UMAP_2 = 0,
        reason = study_info$reason,
        stringsAsFactors = FALSE
      )
      issue_chunks[[length(issue_chunks) + 1]] <- data.frame(
        study_id = study_info$study_id,
        study_label = study_info$study_label,
        gene = gene,
        reason = study_info$reason,
        scope = "study",
        stringsAsFactors = FALSE
      )
      next
    }

    expr_res <- tryCatch(
      extract_gene_values(study_info, gene),
      error = function(e) list(status = "error", expr = NULL, reason = "Expression read error")
    )

    if (expr_res$status != "ok") {
      reason <- if (expr_res$status == "missing_gene") "Gene not found" else expr_res$reason
      placeholder_chunks[[length(placeholder_chunks) + 1]] <- data.frame(
        study_label = study_info$study_label,
        UMAP_1 = 0,
        UMAP_2 = 0,
        reason = reason,
        stringsAsFactors = FALSE
      )
      issue_chunks[[length(issue_chunks) + 1]] <- data.frame(
        study_id = study_info$study_id,
        study_label = study_info$study_label,
        gene = gene,
        reason = reason,
        scope = if (reason == "Gene not found") "gene" else "study",
        stringsAsFactors = FALSE
      )
      next
    }

    chunk <- study_info$coords
    chunk$expr <- expr_res$expr
    chunk$study_label <- study_info$study_label
    chunk <- chunk[, c("study_label", "UMAP_1", "UMAP_2", "expr"), drop = FALSE]
    point_chunks[[length(point_chunks) + 1]] <- chunk
    pooled_expr <- c(pooled_expr, chunk$expr[is.finite(chunk$expr)])
  }

  point_df <- if (length(point_chunks) > 0) do.call(rbind, point_chunks) else data.frame(
    study_label = character(),
    UMAP_1 = numeric(),
    UMAP_2 = numeric(),
    expr = numeric(),
    stringsAsFactors = FALSE
  )
  placeholder_df <- if (length(placeholder_chunks) > 0) do.call(rbind, placeholder_chunks) else data.frame(
    study_label = character(),
    UMAP_1 = numeric(),
    UMAP_2 = numeric(),
    reason = character(),
    stringsAsFactors = FALSE
  )
  issue_df <- if (length(issue_chunks) > 0) do.call(rbind, issue_chunks) else data.frame(
    study_id = character(),
    study_label = character(),
    gene = character(),
    reason = character(),
    scope = character(),
    stringsAsFactors = FALSE
  )

  point_df$study_label <- factor(point_df$study_label, levels = ordered_labels)
  placeholder_df$study_label <- factor(placeholder_df$study_label, levels = ordered_labels)

  limits <- if (length(pooled_expr) > 0) range(pooled_expr, na.rm = TRUE, finite = TRUE) else c(0, 1)
  if (!all(is.finite(limits))) limits <- c(0, 1)
  if (limits[1] == limits[2]) limits[2] <- limits[2] + 1e-6

  facets_seed <- data.frame(
    study_label = factor(ordered_labels, levels = ordered_labels),
    UMAP_1 = 0,
    UMAP_2 = 0,
    stringsAsFactors = FALSE
  )
  legend_seed <- data.frame(
    study_label = factor(ordered_labels[[1]], levels = ordered_labels),
    UMAP_1 = c(0, 0),
    UMAP_2 = c(0, 0),
    expr = limits,
    stringsAsFactors = FALSE
  )

  p <- ggplot() +
    geom_blank(data = facets_seed, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(
      data = legend_seed,
      aes(x = UMAP_1, y = UMAP_2, color = expr),
      inherit.aes = FALSE,
      alpha = 0,
      size = 0.01,
      show.legend = TRUE
    )

  if (nrow(point_df) > 0) {
    p <- p + geom_point(
      data = point_df,
      aes(x = UMAP_1, y = UMAP_2, color = expr),
      size = 0.08,
      alpha = 0.85,
      stroke = 0,
      show.legend = FALSE
    )
  }

  if (nrow(placeholder_df) > 0) {
    p <- p + geom_text(
      data = placeholder_df,
      aes(x = UMAP_1, y = UMAP_2, label = reason),
      inherit.aes = FALSE,
      color = "grey35",
      size = 2.8,
      lineheight = 0.9
    )
  }

  p <- p +
    facet_wrap(~study_label, nrow = 1, drop = FALSE, scales = "fixed") +
    coord_equal() +
    scale_color_gradientn(
      # Use a visible low-expression gray tone to avoid "blank-looking" panels.
      colours = c("#d9d9d9", "#6baed6", "#08306b"),
      limits = limits,
      oob = scales::squish,
      name = gene
    ) +
    guides(
      color = guide_colorbar(
        title.position = "top",
        barheight = grid::unit(22, "pt"),
        frame.colour = "grey65"
      )
    ) +
    labs(title = paste0(gene, " (", gene_group, ")"), x = NULL, y = NULL) +
    theme_minimal(base_size = 8) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "grey80", fill = NA, size = 0.3),
      strip.background = element_rect(fill = "grey95", color = "grey80"),
      strip.text = element_text(size = 8, face = "bold"),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(size = 9, face = "bold", hjust = 0),
      panel.spacing = grid::unit(2, "mm"),
      legend.position = "right",
      legend.title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 7),
      plot.margin = margin(2, 4, 2, 4)
    )

  row_summary <- data.frame(
    gene = gene,
    gene_group = gene_group,
    scale_min = limits[[1]],
    scale_max = limits[[2]],
    n_points = nrow(point_df),
    n_placeholder_panels = nrow(placeholder_df),
    n_studies_with_expression = length(unique(as.character(point_df$study_label))),
    stringsAsFactors = FALSE
  )

  list(plot = p, issues = issue_df, summary = row_summary)
}

print_audit <- function(issue_df) {
  # Compact audit report required for deterministic missing-data tracking.
  cat("\n=== Panel B audit ===\n")
  if (nrow(issue_df) == 0) {
    cat("No missing studies, missing UMAPs, or missing genes.\n")
    return(invisible(NULL))
  }

  study_level <- unique(issue_df[issue_df$scope == "study", c("study_label", "reason"), drop = FALSE])
  if (nrow(study_level) > 0) {
    cat("Study-level placeholders:\n")
    for (i in seq_len(nrow(study_level))) {
      cat(" - ", study_level$study_label[[i]], ": ", study_level$reason[[i]], "\n", sep = "")
    }
  } else {
    cat("Study-level placeholders: none\n")
  }

  gene_level <- issue_df[issue_df$reason == "Gene not found", c("study_label", "gene"), drop = FALSE]
  if (nrow(gene_level) > 0) {
    cat("Gene-not-found placeholders:\n")
    split_genes <- split(gene_level$gene, gene_level$study_label)
    for (study in names(split_genes)) {
      cat(" - ", study, ": ", paste(unique(split_genes[[study]]), collapse = ", "), "\n", sep = "")
    }
  } else {
    cat("Gene-not-found placeholders: none\n")
  }
}

main <- function(cli_args = commandArgs(trailingOnly = TRUE)) {
  # Entry point:
  # - parse/validate args + config
  # - prepare per-study extracts
  # - build all gene rows
  # - write PDF + SVG under PROJECT_ROOT/results/<run_label>/plots
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
  if (is.null(project_root) || !nzchar(project_root)) project_root <- Sys.getenv("PROJECT_ROOT", "")
  if (!nzchar(project_root)) {
    stop("PROJECT_ROOT not set. Use --project-root, config$project_root, or env PROJECT_ROOT.", call. = FALSE)
  }
  project_root <- trim_trailing_slash(project_root)
  project_root <- normalize_abs(project_root, must_work = FALSE)
  if (!dir.exists(project_root)) {
    stop("PROJECT_ROOT does not exist: ", project_root, call. = FALSE)
  }

  run_label <- args[["run-label"]]
  if (is.null(run_label) || !nzchar(run_label)) run_label <- cfg$run_label
  if (!nzchar(run_label)) stop("run_label is empty.", call. = FALSE)
  retain_seurat <- parse_bool_flag(args[["retain-seurat"]], default = FALSE)
  export_global <- parse_bool_flag(args[["export-global"]], default = FALSE)
  show_progress <- parse_bool_flag(args[["show-progress"]], default = FALSE)

  log_msg("Preparing studies from config: ", args$config)
  studies_info <- lapply(seq_len(nrow(cfg$studies)), function(i) {
    prepare_study(cfg$studies[i, , drop = FALSE], project_root, retain_seurat = retain_seurat)
  })
  ordered_labels <- vapply(studies_info, function(x) x$study_label, character(1))
  study_status <- build_study_status_table(studies_info)

  out_dir <- normalize_abs(file.path(project_root, "results", run_label, "plots"), must_work = FALSE)
  if (!is_subpath(out_dir, project_root)) {
    stop("Output path must remain under PROJECT_ROOT.", call. = FALSE)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  log_msg("Building Panel B rows for ", length(GENE_ORDER), " genes across ", length(studies_info), " studies.")
  row_plots <- vector("list", length(GENE_ORDER))
  names(row_plots) <- GENE_ORDER
  issue_chunks <- list()
  row_summary_chunks <- list()
  for (i in seq_along(GENE_ORDER)) {
    gene <- GENE_ORDER[[i]]
    gene_group <- if (gene %in% ON_TARGET_GENES) "ON-target" else "OFF-target"
    row_res <- build_gene_row(gene, gene_group, studies_info, ordered_labels)
    row_plots[[i]] <- row_res$plot
    if (show_progress) {
      print(row_res$plot)
    }
    if (nrow(row_res$issues) > 0) issue_chunks[[length(issue_chunks) + 1]] <- row_res$issues
    row_summary_chunks[[length(row_summary_chunks) + 1]] <- row_res$summary
  }
  row_summary <- do.call(rbind, row_summary_chunks)
  issues <- if (length(issue_chunks) > 0) do.call(rbind, issue_chunks) else data.frame(
    study_id = character(),
    study_label = character(),
    gene = character(),
    reason = character(),
    scope = character(),
    stringsAsFactors = FALSE
  )

  fig <- wrap_plots(row_plots, ncol = 1) +
    plot_layout(heights = rep(1, length(row_plots))) +
    plot_annotation(
      title = "Panel B: Cross-study marker expression on existing UMAPs",
      theme = theme(
        plot.title = element_text(size = 12, face = "bold", hjust = 0)
      )
    )

  fig_width <- max(10, length(studies_info) * 2.15 + 1.8)
  fig_height <- max(12, length(GENE_ORDER) * 1.65 + 1.4)
  pdf_path <- file.path(out_dir, "panel_b_cross_study_markers.pdf")
  svg_path <- file.path(out_dir, "panel_b_cross_study_markers.svg")

  log_msg("Writing PDF: ", pdf_path)
  ggsave(filename = pdf_path, plot = fig, width = fig_width, height = fig_height, units = "in")
  log_msg("Writing SVG: ", svg_path)
  ggsave(filename = svg_path, plot = fig, width = fig_width, height = fig_height, units = "in", device = grDevices::svg)

  log_msg("Done.")
  print_audit(issues)

  result <- list(
    project_root = project_root,
    run_label = run_label,
    genes = GENE_ORDER,
    study_status = study_status,
    row_summary = row_summary,
    issues = issues,
    output_paths = list(pdf = pdf_path, svg = svg_path),
    row_plots = row_plots,
    final_plot = fig,
    studies_info = studies_info
  )

  if (export_global) {
    assign("panel_b_result", result, envir = .GlobalEnv)
    assign("panel_b_studies", study_status, envir = .GlobalEnv)
    assign("panel_b_rows", row_summary, envir = .GlobalEnv)
    assign("panel_b_issues", issues, envir = .GlobalEnv)
    assign("panel_b_row_plots", row_plots, envir = .GlobalEnv)
    assign("panel_b_final_plot", fig, envir = .GlobalEnv)
  }

  invisible(result)
}

run_panel_b_local <- function(config_path,
                              project_root = NULL,
                              run_label = NULL,
                              retain_seurat = FALSE,
                              export_global = TRUE,
                              show_progress_plots = interactive()) {
  # Interactive helper:
  # - returns result list
  # - optionally exports tables/objects to .GlobalEnv for inspection
  args <- c("--config", config_path)
  if (!is.null(project_root)) args <- c(args, "--project-root", project_root)
  if (!is.null(run_label)) args <- c(args, "--run-label", run_label)
  args <- c(
    args,
    "--retain-seurat", if (isTRUE(retain_seurat)) "true" else "false",
    "--export-global", if (isTRUE(export_global)) "true" else "false",
    "--show-progress", if (isTRUE(show_progress_plots)) "true" else "false"
  )
  main(args)
}

if (sys.nframe() == 0) {
  main()
}
