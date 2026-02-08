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
# - PROJECT_ROOT/results/<run_label>/plots/panel_b_study_status.tsv
# - PROJECT_ROOT/results/<run_label>/plots/panel_b_marker_presence.tsv
# - PROJECT_ROOT/results/<run_label>/plots/panel_b_row_summary.tsv
# - PROJECT_ROOT/results/<run_label>/plots/panel_b_issues.tsv
# - PROJECT_ROOT/results/<run_label>/plots/panel_b_assay_slot_summary.tsv
# - PROJECT_ROOT/results/<run_label>/plots/panel_b_reduction_summary.tsv
# - PROJECT_ROOT/results/<run_label>/plots/panel_b_metadata_summary.tsv
# - PROJECT_ROOT/results/<run_label>/plots/panel_b_metadata_columns.tsv
# - PROJECT_ROOT/results/<run_label>/plots/panel_b_ident_counts.tsv
# - PROJECT_ROOT/results/<run_label>/plots/panel_b_feature_space.tsv
# - stdout audit of missing studies/genes/components

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

ON_TARGET_GENES <- c("DCX", "GAD2", "DLX5", "LHX6", "MAF", "SST", "LHX8", "SP8")
OFF_TARGET_GENES <- c("PAX6", "NEUROD2", "ISL1", "ACHF")
GENE_ORDER <- c(ON_TARGET_GENES, OFF_TARGET_GENES)
DETAILED_LOG <- TRUE

log_msg <- function(...) {
  msg <- paste0(..., collapse = " ")
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
  flush.console()
}

set_detailed_log <- function(flag) {
  DETAILED_LOG <<- isTRUE(flag)
}

log_detail <- function(...) {
  if (isTRUE(DETAILED_LOG)) {
    log_msg(...)
  }
}

fmt_count <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x)) return("NA")
  format(round(as.numeric(x), digits = 0), trim = TRUE, scientific = FALSE, big.mark = ",")
}

fmt_num <- function(x, digits = 5) {
  if (is.null(x) || length(x) == 0 || is.na(x) || !is.finite(x)) return("NA")
  format(round(as.numeric(x), digits = digits), trim = TRUE, scientific = FALSE)
}

fmt_bool <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x)) return("NA")
  if (isTRUE(x)) "TRUE" else "FALSE"
}

collapse_csv <- function(x) {
  if (is.null(x) || length(x) == 0) return("")
  paste(as.character(x), collapse = ",")
}

fmt_bytes <- function(bytes) {
  if (is.null(bytes) || length(bytes) == 0 || is.na(bytes) || !is.finite(bytes)) return("NA")
  b <- as.numeric(bytes)
  if (b < 1024) return(paste0(fmt_num(b, digits = 0), " B"))
  if (b < 1024^2) return(paste0(fmt_num(b / 1024, digits = 2), " KB"))
  if (b < 1024^3) return(paste0(fmt_num(b / 1024^2, digits = 2), " MB"))
  paste0(fmt_num(b / 1024^3, digits = 2), " GB")
}

safe_dim <- function(x) {
  if (is.null(x)) return(c(NA_integer_, NA_integer_))
  nr <- tryCatch(nrow(x), error = function(e) NA)
  nc <- tryCatch(ncol(x), error = function(e) NA)
  if (length(nr) == 0 || is.na(nr[[1]]) || !is.finite(as.numeric(nr[[1]]))) nr <- NA_integer_ else nr <- as.integer(nr[[1]])
  if (length(nc) == 0 || is.na(nc[[1]]) || !is.finite(as.numeric(nc[[1]]))) nc <- NA_integer_ else nc <- as.integer(nc[[1]])
  c(nr, nc)
}

safe_range <- function(x) {
  out <- tryCatch(range(x, na.rm = TRUE, finite = TRUE), error = function(e) c(NA_real_, NA_real_))
  if (length(out) != 2) out <- c(NA_real_, NA_real_)
  out
}

safe_ncol <- function(x) {
  nc <- tryCatch(ncol(x), error = function(e) NA)
  if (length(nc) == 0 || is.na(nc[[1]]) || !is.finite(as.numeric(nc[[1]]))) return(NA_integer_)
  as.integer(nc[[1]])
}

safe_nrow <- function(x) {
  nr <- tryCatch(nrow(x), error = function(e) NA)
  if (length(nr) == 0 || is.na(nr[[1]]) || !is.finite(as.numeric(nr[[1]]))) return(NA_integer_)
  as.integer(nr[[1]])
}

detect_feature_id_type <- function(features) {
  if (is.null(features) || length(features) == 0) return("unavailable")
  x <- as.character(features)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) return("unavailable")
  probe <- utils::head(x, 5000)
  n <- length(probe)
  is_ensembl <- grepl("^ENSG[0-9]+(\\.[0-9]+)?$", probe)
  is_symbol <- grepl("^[A-Za-z][A-Za-z0-9_.-]*$", probe) & !is_ensembl
  frac_ensembl <- sum(is_ensembl) / n
  frac_symbol <- sum(is_symbol) / n
  if (frac_ensembl >= 0.8) return("ensembl_id")
  if (frac_symbol >= 0.8) return("symbol_like")
  "mixed_or_other"
}

infer_missing_gene_reason <- function(study_info) {
  if (!is.null(study_info$feature_id_type) && identical(study_info$feature_id_type, "ensembl_id")) {
    return("Gene not found (feature IDs are Ensembl-like)")
  }
  "Gene not found"
}

is_assay_runtime_incompatible <- function(detail, assay_class = "") {
  txt <- tolower(paste(detail, assay_class, collapse = " "))
  grepl("not an assay", txt, fixed = TRUE) ||
    grepl("assay5", txt, fixed = TRUE) ||
    grepl("layer", txt, fixed = TRUE)
}

get_assay_features <- function(obj, assay) {
  f1 <- tryCatch(rownames(obj), error = function(e) character(0))
  if (length(f1) > 0) return(unique(as.character(f1)))
  f2 <- tryCatch(Seurat::Features(obj, assay = assay), error = function(e) character(0))
  if (length(f2) > 0) return(unique(as.character(f2)))
  f3 <- tryCatch(rownames(obj[[assay]]), error = function(e) character(0))
  if (length(f3) > 0) return(unique(as.character(f3)))
  character(0)
}

fetch_marker_matrix <- function(obj, assay, preferred_slot, cells, genes) {
  # Fallback extractor for assay classes where GetAssayData/slots are unavailable.
  old_assay <- tryCatch(as.character(DefaultAssay(obj)), error = function(e) "")
  on.exit(
    {
      if (nzchar(old_assay)) {
        tryCatch(DefaultAssay(obj) <- old_assay, error = function(e) NULL)
      }
    },
    add = TRUE
  )
  tryCatch(DefaultAssay(obj) <- assay, error = function(e) NULL)

  features <- get_assay_features(obj, assay)
  features_known <- length(features) > 0
  genes_present <- if (features_known) genes[genes %in% features] else character(0)
  genes_missing <- setdiff(genes, genes_present)
  vars_to_fetch <- if (features_known) genes_present else genes
  slot_candidates <- unique(c(preferred_slot, "data", "counts"))

  if (features_known && length(genes_present) == 0) {
    return(list(
      status = "ok",
      expression_slot = preferred_slot,
      expr_sub = NULL,
      genes_present = genes_present,
      genes_missing = genes_missing,
      n_features_assay = ifelse(length(features) > 0, length(features), NA_integer_),
      n_cells_assay = length(cells),
      n_cells_common = length(cells),
      n_cells_umap_not_in_assay = 0L,
      n_cells_assay_not_in_umap = 0L,
      detail = "No requested marker genes were found in assay features."
    ))
  }

  last_error <- NULL
  for (slot_name in slot_candidates) {
    fetched <- tryCatch(
      Seurat::FetchData(obj, vars = vars_to_fetch, cells = cells, slot = slot_name),
      error = function(e) {
        last_error <<- e
        NULL
      }
    )
    if (is.null(fetched) || nrow(fetched) == 0 || ncol(fetched) == 0) next

    common_cells <- intersect(cells, rownames(fetched))
    if (length(common_cells) == 0) next
    fetched <- fetched[common_cells, , drop = FALSE]
    cols <- intersect(colnames(fetched), genes)
    if (length(cols) == 0) next

    mat <- t(as.matrix(fetched[, cols, drop = FALSE]))
    rownames(mat) <- cols
    colnames(mat) <- rownames(fetched)

    return(list(
      status = "ok",
      expression_slot = slot_name,
      expr_sub = mat,
      genes_present = rownames(mat),
      genes_missing = setdiff(genes, rownames(mat)),
      n_features_assay = ifelse(length(features) > 0, length(features), NA_integer_),
      n_cells_assay = nrow(fetched),
      n_cells_common = length(common_cells),
      n_cells_umap_not_in_assay = length(setdiff(cells, rownames(fetched))),
      n_cells_assay_not_in_umap = length(setdiff(rownames(fetched), cells)),
      detail = ""
    ))
  }

  list(
    status = "error",
    expression_slot = NA_character_,
    expr_sub = NULL,
    genes_present = character(0),
    genes_missing = genes,
    n_features_assay = ifelse(length(features) > 0, length(features), NA_integer_),
    n_cells_assay = NA_integer_,
    n_cells_common = NA_integer_,
    n_cells_umap_not_in_assay = NA_integer_,
    n_cells_assay_not_in_umap = NA_integer_,
    detail = if (!is.null(last_error)) conditionMessage(last_error) else "FetchData fallback returned no usable data."
  )
}

summarize_assay_slots <- function(obj, study_id, study_label) {
  assays <- names(obj@assays)
  slots <- c("counts", "data", "scale.data")
  chunks <- list()
  k <- 0
  for (assay_name in assays) {
    for (slot_name in slots) {
      k <- k + 1
      mat <- tryCatch(
        suppressWarnings(GetAssayData(obj, assay = assay_name, slot = slot_name)),
        error = function(e) e
      )
      if (inherits(mat, "error") || is.null(mat)) {
        chunks[[k]] <- data.frame(
          study_id = study_id,
          study_label = study_label,
          assay = assay_name,
          slot = slot_name,
          status = "missing_or_error",
          matrix_class = "",
          n_features = NA_integer_,
          n_cells = NA_integer_,
          detail = if (inherits(mat, "error")) conditionMessage(mat) else "slot missing",
          stringsAsFactors = FALSE
        )
      } else {
        d <- safe_dim(mat)
        chunks[[k]] <- data.frame(
          study_id = study_id,
          study_label = study_label,
          assay = assay_name,
          slot = slot_name,
          status = ifelse(!is.na(d[[1]]) && !is.na(d[[2]]) && d[[1]] > 0 && d[[2]] > 0, "ok", "empty"),
          matrix_class = collapse_csv(class(mat)),
          n_features = d[[1]],
          n_cells = d[[2]],
          detail = "",
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(chunks) == 0) {
    return(data.frame(
      study_id = character(),
      study_label = character(),
      assay = character(),
      slot = character(),
      status = character(),
      matrix_class = character(),
      n_features = integer(),
      n_cells = integer(),
      detail = character(),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, chunks)
}

summarize_reductions <- function(obj, study_id, study_label) {
  reductions <- names(obj@reductions)
  chunks <- list()
  k <- 0
  for (reduction_name in reductions) {
    k <- k + 1
    emb <- tryCatch(Embeddings(obj, reduction = reduction_name), error = function(e) e)
    if (inherits(emb, "error") || is.null(emb)) {
      chunks[[k]] <- data.frame(
        study_id = study_id,
        study_label = study_label,
        reduction = reduction_name,
        status = "missing_or_error",
        n_cells = NA_integer_,
        n_dims = NA_integer_,
        dim1_min = NA_real_,
        dim1_max = NA_real_,
        dim2_min = NA_real_,
        dim2_max = NA_real_,
        detail = if (inherits(emb, "error")) conditionMessage(emb) else "reduction missing",
        stringsAsFactors = FALSE
      )
    } else {
      d <- safe_dim(emb)
      r1 <- if (!is.na(d[[2]]) && d[[2]] >= 1) safe_range(emb[, 1]) else c(NA_real_, NA_real_)
      r2 <- if (!is.na(d[[2]]) && d[[2]] >= 2) safe_range(emb[, 2]) else c(NA_real_, NA_real_)
      chunks[[k]] <- data.frame(
        study_id = study_id,
        study_label = study_label,
        reduction = reduction_name,
        status = "ok",
        n_cells = d[[1]],
        n_dims = d[[2]],
        dim1_min = r1[[1]],
        dim1_max = r1[[2]],
        dim2_min = r2[[1]],
        dim2_max = r2[[2]],
        detail = "",
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(chunks) == 0) {
    return(data.frame(
      study_id = character(),
      study_label = character(),
      reduction = character(),
      status = character(),
      n_cells = integer(),
      n_dims = integer(),
      dim1_min = numeric(),
      dim1_max = numeric(),
      dim2_min = numeric(),
      dim2_max = numeric(),
      detail = character(),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, chunks)
}

summarize_metadata <- function(obj, study_id, study_label) {
  md <- tryCatch(obj@meta.data, error = function(e) NULL)
  if (is.null(md)) {
    return(list(
      summary = data.frame(
        study_id = study_id,
        study_label = study_label,
        status = "error",
        n_meta_rows = NA_integer_,
        n_meta_cols = NA_integer_,
        meta_rows_match_object_cells = NA,
        has_seurat_clusters = NA,
        n_seurat_clusters = NA_integer_,
        has_sample_id = NA,
        n_sample_id = NA_integer_,
        has_orig_ident = NA,
        n_orig_ident = NA_integer_,
        has_domain = NA,
        n_domain = NA_integer_,
        metadata_columns_preview = "",
        stringsAsFactors = FALSE
      ),
      columns = data.frame(
        study_id = character(),
        study_label = character(),
        column_name = character(),
        column_class = character(),
        n_non_na = integer(),
        n_unique_non_na = integer(),
        example_values = character(),
        stringsAsFactors = FALSE
      )
    ))
  }

  cols <- colnames(md)
  has_sc <- "seurat_clusters" %in% cols
  has_sample <- "sample_id" %in% cols
  has_orig <- "orig.ident" %in% cols
  has_domain <- "domain" %in% cols

  safe_unique <- function(x) {
    tryCatch(length(unique(as.character(x[!is.na(x)]))), error = function(e) NA_integer_)
  }
  preview_values <- function(x, n = 3) {
    vals <- tryCatch(unique(as.character(x[!is.na(x)])), error = function(e) character(0))
    vals <- vals[nzchar(vals)]
    collapse_csv(utils::head(vals, n))
  }

  summary_row <- data.frame(
    study_id = study_id,
    study_label = study_label,
    status = "ok",
    n_meta_rows = tryCatch(as.integer(nrow(md)), error = function(e) NA_integer_),
    n_meta_cols = tryCatch(as.integer(ncol(md)), error = function(e) NA_integer_),
    meta_rows_match_object_cells = tryCatch(identical(rownames(md), colnames(obj)), error = function(e) NA),
    has_seurat_clusters = has_sc,
    n_seurat_clusters = if (has_sc) safe_unique(md$seurat_clusters) else NA_integer_,
    has_sample_id = has_sample,
    n_sample_id = if (has_sample) safe_unique(md$sample_id) else NA_integer_,
    has_orig_ident = has_orig,
    n_orig_ident = if (has_orig) safe_unique(md$orig.ident) else NA_integer_,
    has_domain = has_domain,
    n_domain = if (has_domain) safe_unique(md$domain) else NA_integer_,
    metadata_columns_preview = collapse_csv(utils::head(cols, 20)),
    stringsAsFactors = FALSE
  )

  column_rows <- list()
  for (i in seq_along(cols)) {
    col_name <- cols[[i]]
    v <- md[[col_name]]
    column_rows[[i]] <- data.frame(
      study_id = study_id,
      study_label = study_label,
      column_name = col_name,
      column_class = collapse_csv(class(v)),
      n_non_na = tryCatch(sum(!is.na(v)), error = function(e) NA_integer_),
      n_unique_non_na = safe_unique(v),
      example_values = preview_values(v),
      stringsAsFactors = FALSE
    )
  }
  columns_df <- if (length(column_rows) > 0) do.call(rbind, column_rows) else data.frame(
    study_id = character(),
    study_label = character(),
    column_name = character(),
    column_class = character(),
    n_non_na = integer(),
    n_unique_non_na = integer(),
    example_values = character(),
    stringsAsFactors = FALSE
  )

  list(summary = summary_row, columns = columns_df)
}

summarize_idents <- function(obj, study_id, study_label) {
  ids <- tryCatch(Idents(obj), error = function(e) NULL)
  if (is.null(ids)) {
    return(data.frame(
      study_id = character(),
      study_label = character(),
      ident_value = character(),
      n_cells = integer(),
      stringsAsFactors = FALSE
    ))
  }
  tb <- sort(table(as.character(ids)), decreasing = TRUE)
  data.frame(
    study_id = study_id,
    study_label = study_label,
    ident_value = names(tb),
    n_cells = as.integer(tb),
    stringsAsFactors = FALSE
  )
}

combine_table_field <- function(studies_info, field_name, empty_table) {
  chunks <- lapply(studies_info, function(study_info) study_info[[field_name]])
  chunks <- chunks[!vapply(chunks, is.null, logical(1))]
  chunks <- chunks[vapply(chunks, function(x) is.data.frame(x) && nrow(x) > 0, logical(1))]
  if (length(chunks) == 0) return(empty_table)
  out <- do.call(rbind, chunks)
  rownames(out) <- NULL
  out
}

print_usage <- function() {
  cat(
    paste(
      "Usage:",
      "  Rscript scripts/06_cross_study_panelB_markers.R --config <config.R> [--project-root <path>] [--run-label <label>] [--retain-seurat <true|false>] [--export-global <true|false>] [--detailed-log <true|false>]",
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
      "  --detailed-log true|false     print per-study diagnostics (object structure, assay/reduction dims, cell matching)",
      "",
      "Outputs:",
      "  PROJECT_ROOT/results/<run_label>/plots/panel_b_cross_study_markers.{pdf,svg}",
      "  PROJECT_ROOT/results/<run_label>/plots/panel_b_study_status.tsv",
      "  PROJECT_ROOT/results/<run_label>/plots/panel_b_marker_presence.tsv",
      "  PROJECT_ROOT/results/<run_label>/plots/panel_b_row_summary.tsv",
      "  PROJECT_ROOT/results/<run_label>/plots/panel_b_issues.tsv",
      "  PROJECT_ROOT/results/<run_label>/plots/panel_b_assay_slot_summary.tsv",
      "  PROJECT_ROOT/results/<run_label>/plots/panel_b_reduction_summary.tsv",
      "  PROJECT_ROOT/results/<run_label>/plots/panel_b_metadata_summary.tsv",
      "  PROJECT_ROOT/results/<run_label>/plots/panel_b_metadata_columns.tsv",
      "  PROJECT_ROOT/results/<run_label>/plots/panel_b_ident_counts.tsv",
      "  PROJECT_ROOT/results/<run_label>/plots/panel_b_feature_space.tsv",
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
    assay_slot_requested = study_row$expression_slot,
    expression_slot = NA_character_,
    status = "ok",
    reason = NA_character_,
    detail = NA_character_,
    seurat_obj = NULL,
    expr_sub = NULL,
    coords = NULL,
    object_bytes = NA_real_,
    object_class = NA_character_,
    default_assay = NA_character_,
    assays_available = NA_character_,
    assay_class = NA_character_,
    assay_slot_names = NA_character_,
    reductions_available = NA_character_,
    n_cells_object = NA_integer_,
    n_features_assay = NA_integer_,
    n_features_detected = NA_integer_,
    feature_id_type = NA_character_,
    feature_id_examples = NA_character_,
    data_access_mode = NA_character_,
    runtime_note = NA_character_,
    n_cells_assay = NA_integer_,
    n_cells_umap = NA_integer_,
    n_dims_umap = NA_integer_,
    n_cells_common = NA_integer_,
    n_cells_umap_not_in_assay = NA_integer_,
    n_cells_assay_not_in_umap = NA_integer_,
    n_marker_genes_requested = length(GENE_ORDER),
    n_marker_genes_present = NA_integer_,
    marker_genes_present = NA_character_,
    marker_genes_missing = NA_character_,
    assay_slot_dims = NA_character_,
    umap_dim1_min = NA_real_,
    umap_dim1_max = NA_real_,
    umap_dim2_min = NA_real_,
    umap_dim2_max = NA_real_,
    n_meta_rows = NA_integer_,
    n_meta_cols = NA_integer_,
    meta_rows_match_object_cells = NA,
    has_seurat_clusters = NA,
    n_seurat_clusters = NA_integer_,
    has_sample_id = NA,
    n_sample_id = NA_integer_,
    has_orig_ident = NA,
    n_orig_ident = NA_integer_,
    has_domain = NA,
    n_domain = NA_integer_,
    metadata_columns_preview = NA_character_,
    n_ident_levels = NA_integer_,
    n_ident_cells = NA_integer_,
    assay_slot_summary = NULL,
    reduction_summary = NULL,
    metadata_summary = NULL,
    metadata_columns = NULL,
    ident_counts = NULL,
    load_seconds = NA_real_
  )
  start_ts <- Sys.time()

  log_detail("---- Study ", info$study_label, " (", info$study_id, ") ----")
  log_detail("Configured object_path: ", study_row$object_path)
  log_detail("Resolved object path: ", info$object_path)

  if (!file.exists(info$object_path)) {
    info$status <- "missing_object"
    info$reason <- "Missing object"
    info$detail <- "Object file does not exist."
    log_detail("Status: missing_object | reason: Missing object")
    return(info)
  }

  fi <- tryCatch(file.info(info$object_path), error = function(e) NULL)
  if (!is.null(fi) && nrow(fi) == 1 && !is.na(fi$size[[1]])) {
    info$object_bytes <- as.numeric(fi$size[[1]])
    log_detail("Object size on disk: ", fmt_bytes(info$object_bytes), " (", fmt_count(info$object_bytes), " bytes)")
  }

  obj <- tryCatch(read_rds_any(info$object_path), error = function(e) e)
  if (inherits(obj, "error")) {
    info$status <- "unreadable_object"
    info$reason <- "Unreadable object"
    info$detail <- conditionMessage(obj)
    log_detail("Status: unreadable_object | ", info$detail)
    return(info)
  }

  info$object_class <- collapse_csv(class(obj))
  log_detail("Loaded object class: ", info$object_class)

  if (!inherits(obj, "Seurat")) {
    info$status <- "invalid_object"
    info$reason <- "Invalid object"
    info$detail <- paste("Class:", info$object_class)
    log_detail("Status: invalid_object | ", info$detail)
    return(info)
  }
  if (retain_seurat) {
    info$seurat_obj <- obj
  }

  info$assay_slot_summary <- summarize_assay_slots(obj, info$study_id, info$study_label)
  info$reduction_summary <- summarize_reductions(obj, info$study_id, info$study_label)
  md_res <- summarize_metadata(obj, info$study_id, info$study_label)
  info$metadata_summary <- md_res$summary
  info$metadata_columns <- md_res$columns
  info$ident_counts <- summarize_idents(obj, info$study_id, info$study_label)

  if (nrow(info$metadata_summary) > 0) {
    ms <- info$metadata_summary[1, , drop = FALSE]
    info$n_meta_rows <- ms$n_meta_rows[[1]]
    info$n_meta_cols <- ms$n_meta_cols[[1]]
    info$meta_rows_match_object_cells <- ms$meta_rows_match_object_cells[[1]]
    info$has_seurat_clusters <- ms$has_seurat_clusters[[1]]
    info$n_seurat_clusters <- ms$n_seurat_clusters[[1]]
    info$has_sample_id <- ms$has_sample_id[[1]]
    info$n_sample_id <- ms$n_sample_id[[1]]
    info$has_orig_ident <- ms$has_orig_ident[[1]]
    info$n_orig_ident <- ms$n_orig_ident[[1]]
    info$has_domain <- ms$has_domain[[1]]
    info$n_domain <- ms$n_domain[[1]]
    info$metadata_columns_preview <- ms$metadata_columns_preview[[1]]
  }
  if (!is.null(info$ident_counts) && nrow(info$ident_counts) > 0) {
    info$n_ident_levels <- nrow(info$ident_counts)
    info$n_ident_cells <- sum(info$ident_counts$n_cells, na.rm = TRUE)
  }

  info$default_assay <- tryCatch(as.character(DefaultAssay(obj)), error = function(e) "")
  info$assays_available <- collapse_csv(names(obj@assays))
  info$reductions_available <- collapse_csv(names(obj@reductions))
  info$n_cells_object <- safe_ncol(obj)
  if (is.na(info$n_cells_object) && !is.na(info$n_meta_rows) && info$n_meta_rows > 0) {
    info$n_cells_object <- as.integer(info$n_meta_rows)
  }
  log_detail("Default assay: ", ifelse(nzchar(info$default_assay), info$default_assay, "<none>"))
  log_detail("Available assays: ", ifelse(nzchar(info$assays_available), info$assays_available, "<none>"))
  log_detail("Available reductions: ", ifelse(nzchar(info$reductions_available), info$reductions_available, "<none>"))
  log_detail("Object cells (ncol): ", fmt_count(info$n_cells_object))
  log_detail("Configured assay/reduction/slot: ", info$assay, " / ", info$reduction, " / ", info$preferred_slot)
  log_detail(
    "Metadata rows x cols: ", fmt_count(info$n_meta_rows), " x ", fmt_count(info$n_meta_cols),
    "; rows_match_object_cells=", fmt_bool(info$meta_rows_match_object_cells)
  )
  log_detail(
    "Metadata key columns: seurat_clusters=", fmt_bool(info$has_seurat_clusters), " (n=", fmt_count(info$n_seurat_clusters),
    "), sample_id=", fmt_bool(info$has_sample_id), " (n=", fmt_count(info$n_sample_id),
    "), orig.ident=", fmt_bool(info$has_orig_ident), " (n=", fmt_count(info$n_orig_ident),
    "), domain=", fmt_bool(info$has_domain), " (n=", fmt_count(info$n_domain), ")"
  )
  log_detail("Metadata column preview: ", ifelse(nzchar(info$metadata_columns_preview), info$metadata_columns_preview, "<none>"))
  log_detail("Idents: levels=", fmt_count(info$n_ident_levels), "; cells=", fmt_count(info$n_ident_cells))

  if (!is.null(info$assay_slot_summary) && nrow(info$assay_slot_summary) > 0) {
    log_detail("Assay-slot inventory:")
    for (k in seq_len(nrow(info$assay_slot_summary))) {
      rr <- info$assay_slot_summary[k, , drop = FALSE]
      log_detail(
        "  - ", rr$assay[[1]], "::", rr$slot[[1]], " -> ",
        fmt_count(rr$n_features[[1]]), " x ", fmt_count(rr$n_cells[[1]]),
        " [", rr$status[[1]], "] ", rr$matrix_class[[1]]
      )
    }
  }
  if (!is.null(info$reduction_summary) && nrow(info$reduction_summary) > 0) {
    log_detail("Reduction inventory:")
    for (k in seq_len(nrow(info$reduction_summary))) {
      rr <- info$reduction_summary[k, , drop = FALSE]
      log_detail(
        "  - ", rr$reduction[[1]], " -> ",
        fmt_count(rr$n_cells[[1]]), " x ", fmt_count(rr$n_dims[[1]]),
        " [", rr$status[[1]], "]",
        "; dim1=[", fmt_num(rr$dim1_min[[1]]), ",", fmt_num(rr$dim1_max[[1]]), "]",
        "; dim2=[", fmt_num(rr$dim2_min[[1]]), ",", fmt_num(rr$dim2_max[[1]]), "]"
      )
    }
  }
  if (!is.null(info$ident_counts) && nrow(info$ident_counts) > 0) {
    top_ids <- utils::head(info$ident_counts, 10)
    log_detail("Top identity levels:")
    for (k in seq_len(nrow(top_ids))) {
      rr <- top_ids[k, , drop = FALSE]
      log_detail("  - ", rr$ident_value[[1]], ": ", fmt_count(rr$n_cells[[1]]), " cells")
    }
  }

  if (!(info$assay %in% names(obj@assays))) {
    info$status <- "missing_assay"
    info$reason <- "Missing assay"
    info$detail <- paste0(
      "Configured assay=", info$assay,
      "; available assays=", ifelse(nzchar(info$assays_available), info$assays_available, "<none>")
    )
    log_detail("Status: missing_assay | ", info$detail)
    return(info)
  }

  assay_obj <- tryCatch(obj[[info$assay]], error = function(e) NULL)
  assay_dims <- safe_dim(assay_obj)
  info$assay_class <- if (!is.null(assay_obj)) collapse_csv(class(assay_obj)) else ""
  info$assay_slot_names <- if (!is.null(assay_obj)) {
    collapse_csv(tryCatch(slotNames(assay_obj), error = function(e) character(0)))
  } else {
    ""
  }
  assay_features <- get_assay_features(obj, info$assay)
  if (length(assay_features) > 0) {
    info$n_features_detected <- as.integer(length(unique(assay_features)))
    info$feature_id_type <- detect_feature_id_type(assay_features)
    info$feature_id_examples <- collapse_csv(utils::head(assay_features, 8))
  } else {
    info$n_features_detected <- NA_integer_
    info$feature_id_type <- "unavailable"
    info$feature_id_examples <- ""
  }
  info$n_features_assay <- assay_dims[[1]]
  if (is.na(info$n_features_assay) && !is.na(info$n_features_detected)) {
    info$n_features_assay <- info$n_features_detected
  }
  if (!is.na(assay_dims[[2]])) {
    info$n_cells_assay <- assay_dims[[2]]
  }
  log_detail(
    "Assay object dims (features x cells): ",
    fmt_count(info$n_features_assay), " x ", fmt_count(assay_dims[[2]])
  )
  log_detail(
    "Assay class/slots: ",
    ifelse(nzchar(info$assay_class), info$assay_class, "<unknown>"),
    " | slots=", ifelse(nzchar(info$assay_slot_names), info$assay_slot_names, "<none>")
  )
  log_detail(
    "Feature namespace: type=", ifelse(nzchar(info$feature_id_type), info$feature_id_type, "<unknown>"),
    "; detected_features=", fmt_count(info$n_features_detected),
    ifelse(nzchar(info$feature_id_examples), paste0("; examples=", info$feature_id_examples), "")
  )
  if (identical(info$feature_id_type, "ensembl_id")) {
    info$runtime_note <- "Feature IDs are Ensembl-like; gene symbols require mapping for marker lookup."
    log_detail("Runtime note: ", info$runtime_note)
  }
  if (!is.na(info$assay_class) && grepl("Assay5", info$assay_class, fixed = TRUE)) {
    note <- "Assay5 detected; full expression access may require Seurat v5 runtime."
    info$runtime_note <- if (is.na(info$runtime_note) || !nzchar(info$runtime_note)) {
      note
    } else {
      paste(info$runtime_note, note, sep = " | ")
    }
    log_detail("Runtime note: ", note)
  }

  if (!(info$reduction %in% names(obj@reductions))) {
    info$status <- "missing_umap"
    info$reason <- "Missing UMAP"
    info$detail <- paste0(
      "Configured reduction=", info$reduction,
      "; available reductions=", ifelse(nzchar(info$reductions_available), info$reductions_available, "<none>")
    )
    log_detail("Status: missing_umap | ", info$detail)
    return(info)
  }

  coords_raw <- tryCatch(Embeddings(obj, reduction = info$reduction), error = function(e) e)
  if (inherits(coords_raw, "error")) {
    info$status <- "missing_umap"
    info$reason <- "Missing UMAP"
    info$detail <- paste0("Failed to read reduction=", info$reduction, "; ", conditionMessage(coords_raw))
    log_detail("Status: missing_umap | ", info$detail)
    return(info)
  }

  coords_dims <- safe_dim(coords_raw)
  info$n_cells_umap <- coords_dims[[1]]
  info$n_dims_umap <- coords_dims[[2]]
  if (is.na(info$n_cells_umap) || is.na(info$n_dims_umap) || info$n_cells_umap == 0 || info$n_dims_umap < 2) {
    info$status <- "missing_umap"
    info$reason <- "Missing UMAP"
    info$detail <- paste0(
      "Reduction has invalid dims (cells x dims)=",
      fmt_count(info$n_cells_umap), " x ", fmt_count(info$n_dims_umap)
    )
    log_detail("Status: missing_umap | ", info$detail)
    return(info)
  }

  d1 <- safe_range(coords_raw[, 1])
  d2 <- safe_range(coords_raw[, 2])
  info$umap_dim1_min <- d1[[1]]
  info$umap_dim1_max <- d1[[2]]
  info$umap_dim2_min <- d2[[1]]
  info$umap_dim2_max <- d2[[2]]
  log_detail(
    "Reduction dims (cells x dims): ",
    fmt_count(info$n_cells_umap), " x ", fmt_count(info$n_dims_umap)
  )
  log_detail(
    "Reduction ranges: UMAP_1=[", fmt_num(info$umap_dim1_min), ", ", fmt_num(info$umap_dim1_max),
    "], UMAP_2=[", fmt_num(info$umap_dim2_min), ", ", fmt_num(info$umap_dim2_max), "]"
  )

  coords <- as.data.frame(coords_raw[, 1:2, drop = FALSE], stringsAsFactors = FALSE)
  colnames(coords) <- c("UMAP_1", "UMAP_2")
  coords$cell_id <- rownames(coords)
  cells <- coords$cell_id

  mat <- NULL
  expr_sub <- NULL
  gene_hits <- character(0)
  gene_missing <- GENE_ORDER
  info$data_access_mode <- "unresolved"
  slot_name <- choose_expression_slot(obj, info$assay, info$preferred_slot)
  direct_assay_error <- NA_character_

  if (!is.na(slot_name)) {
    info$expression_slot <- slot_name
    log_detail("Expression slot selected: ", slot_name, " (requested ", info$preferred_slot, ")")
    mat_try <- tryCatch(
      suppressWarnings(GetAssayData(obj, assay = info$assay, slot = slot_name)),
      error = function(e) e
    )
    mat_try_nrow <- if (!inherits(mat_try, "error") && !is.null(mat_try)) safe_nrow(mat_try) else NA_integer_
    mat_try_ncol <- if (!inherits(mat_try, "error") && !is.null(mat_try)) safe_ncol(mat_try) else NA_integer_
    mat_try_invalid <- inherits(mat_try, "error") ||
      is.null(mat_try) ||
      is.na(mat_try_nrow) || is.na(mat_try_ncol) ||
      mat_try_nrow == 0 || mat_try_ncol == 0
    if (mat_try_invalid) {
      direct_assay_error <- if (inherits(mat_try, "error")) {
        paste0("assay=", info$assay, "; slot=", slot_name, "; ", conditionMessage(mat_try))
      } else {
        paste0("assay=", info$assay, "; slot=", slot_name, " has empty matrix")
      }
      log_detail("GetAssayData unavailable/empty; switching to FetchData fallback: ", direct_assay_error)
    } else {
      mat <- mat_try
      info$data_access_mode <- "direct_getassaydata"
    }
  } else {
    direct_assay_error <- paste0(
      "No non-empty slot among GetAssayData candidates: ",
      collapse_csv(unique(c(info$preferred_slot, "data", "counts")))
    )
    log_detail(direct_assay_error, "; switching to FetchData fallback.")
  }

  if (!is.null(mat)) {
    mat_dims <- safe_dim(mat)
    info$n_features_assay <- mat_dims[[1]]
    info$n_cells_assay <- mat_dims[[2]]
    info$assay_slot_dims <- paste0(fmt_count(mat_dims[[1]]), "x", fmt_count(mat_dims[[2]]))
    log_detail(
      "Assay slot dims (features x cells): ",
      fmt_count(mat_dims[[1]]), " x ", fmt_count(mat_dims[[2]])
    )

    umap_cells <- unique(cells)
    assay_cells <- unique(colnames(mat))
    common_cells <- intersect(umap_cells, assay_cells)
    missing_in_assay <- setdiff(umap_cells, assay_cells)
    extra_in_assay <- setdiff(assay_cells, umap_cells)
    info$n_cells_common <- length(common_cells)
    info$n_cells_umap_not_in_assay <- length(missing_in_assay)
    info$n_cells_assay_not_in_umap <- length(extra_in_assay)
    log_detail(
      "Cell overlap: common=", fmt_count(info$n_cells_common),
      "; UMAP-only=", fmt_count(info$n_cells_umap_not_in_assay),
      "; assay-only=", fmt_count(info$n_cells_assay_not_in_umap)
    )

    if (length(missing_in_assay) > 0) {
      info$status <- "missing_assay_data"
      info$reason <- "Cell mismatch"
      info$detail <- paste0(
        "UMAP cells not in assay matrix: ", fmt_count(length(missing_in_assay)),
        if (length(missing_in_assay) > 0) {
          paste0(" (examples: ", paste(utils::head(missing_in_assay, 5), collapse = ", "), ")")
        } else {
          ""
        }
      )
      log_detail("Status: missing_assay_data | ", info$detail)
      return(info)
    }

    gene_hits <- GENE_ORDER[GENE_ORDER %in% rownames(mat)]
    gene_missing <- setdiff(GENE_ORDER, gene_hits)
    expr_sub <- if (length(gene_hits) > 0) {
      mat[gene_hits, cells, drop = FALSE]
    } else {
      NULL
    }
  } else {
    fallback_res <- fetch_marker_matrix(
      obj = obj,
      assay = info$assay,
      preferred_slot = info$preferred_slot,
      cells = cells,
      genes = GENE_ORDER
    )
    if (!identical(fallback_res$status, "ok")) {
      info$data_access_mode <- "assay_data_unavailable"
      detail_parts <- character(0)
      if (!is.na(direct_assay_error) && nzchar(direct_assay_error)) {
        detail_parts <- c(detail_parts, direct_assay_error)
      }
      if (!is.null(fallback_res$detail) && nzchar(fallback_res$detail)) {
        detail_parts <- c(detail_parts, paste0("FetchData fallback: ", fallback_res$detail))
      }
      if (length(detail_parts) == 0) {
        detail_parts <- "Unable to retrieve marker expression by GetAssayData or FetchData."
      }
      info$detail <- paste(detail_parts, collapse = " | ")
      if (is_assay_runtime_incompatible(info$detail, info$assay_class)) {
        info$status <- "assay_runtime_incompatible"
        info$reason <- "Assay runtime incompatible"
        info$data_access_mode <- "assay_runtime_incompatible"
        note <- paste0(
          "Assay=", info$assay,
          "; class=", ifelse(nzchar(info$assay_class), info$assay_class, "<unknown>"),
          "; runtime cannot read this assay type."
        )
        info$runtime_note <- if (is.na(info$runtime_note) || !nzchar(info$runtime_note)) note else paste(info$runtime_note, note, sep = " | ")
        log_detail("Status: assay_runtime_incompatible | ", info$detail)
      } else {
        info$status <- "missing_assay_data"
        info$reason <- "Missing assay data"
        log_detail("Status: missing_assay_data | ", info$detail)
      }
      return(info)
    }

    info$data_access_mode <- "fetchdata_fallback"
    info$expression_slot <- fallback_res$expression_slot
    if (!is.na(fallback_res$n_features_assay)) {
      info$n_features_assay <- as.integer(fallback_res$n_features_assay)
    }
    if (!is.na(fallback_res$n_cells_assay)) {
      info$n_cells_assay <- as.integer(fallback_res$n_cells_assay)
    }
    info$n_cells_common <- as.integer(fallback_res$n_cells_common)
    info$n_cells_umap_not_in_assay <- as.integer(fallback_res$n_cells_umap_not_in_assay)
    info$n_cells_assay_not_in_umap <- as.integer(fallback_res$n_cells_assay_not_in_umap)
    expr_sub <- fallback_res$expr_sub
    gene_hits <- fallback_res$genes_present
    gene_missing <- setdiff(GENE_ORDER, gene_hits)

    if (!is.null(expr_sub)) {
      common_cells <- intersect(cells, colnames(expr_sub))
      if (length(common_cells) > 0) {
        expr_sub <- expr_sub[, common_cells, drop = FALSE]
        coords <- coords[match(common_cells, coords$cell_id), , drop = FALSE]
        rownames(coords) <- coords$cell_id
      } else {
        expr_sub <- NULL
      }
    }

    if (!is.null(expr_sub)) {
      expr_dims <- safe_dim(expr_sub)
      info$assay_slot_dims <- paste0(fmt_count(expr_dims[[1]]), "x", fmt_count(expr_dims[[2]]))
    } else {
      info$assay_slot_dims <- "0x0"
    }

    log_detail(
      "Expression data resolved via FetchData fallback using slot ",
      info$expression_slot,
      " (requested ", info$preferred_slot, ")"
    )
    log_detail(
      "Cell overlap: common=", fmt_count(info$n_cells_common),
      "; UMAP-only=", fmt_count(info$n_cells_umap_not_in_assay),
      "; assay-only=", fmt_count(info$n_cells_assay_not_in_umap)
    )

    if (is.finite(info$n_cells_umap_not_in_assay) && info$n_cells_umap_not_in_assay > 0) {
      info$status <- "missing_assay_data"
      info$reason <- "Cell mismatch"
      info$detail <- paste0(
        "UMAP cells not in assay matrix: ",
        fmt_count(info$n_cells_umap_not_in_assay)
      )
      log_detail("Status: missing_assay_data | ", info$detail)
      return(info)
    }
  }

  info$n_marker_genes_present <- length(gene_hits)
  info$marker_genes_present <- collapse_csv(gene_hits)
  info$marker_genes_missing <- collapse_csv(gene_missing)
  if (info$n_marker_genes_present == 0 && identical(info$feature_id_type, "ensembl_id")) {
    note <- "None of the symbol-based markers were found because assay features are Ensembl-like."
    info$runtime_note <- if (is.na(info$runtime_note) || !nzchar(info$runtime_note)) note else paste(info$runtime_note, note, sep = " | ")
  }
  log_detail(
    "Marker genes: requested=", fmt_count(info$n_marker_genes_requested),
    "; present=", fmt_count(info$n_marker_genes_present),
    "; missing=", fmt_count(length(gene_missing))
  )
  if (length(gene_missing) > 0) {
    log_detail("Missing marker genes: ", collapse_csv(gene_missing))
  }
  if (!is.null(expr_sub)) {
    expr_dims <- safe_dim(expr_sub)
    log_detail(
      "Extracted marker matrix dims (features x cells): ",
      fmt_count(expr_dims[[1]]), " x ", fmt_count(expr_dims[[2]])
    )
  } else {
    log_detail("Extracted marker matrix dims: 0 x 0")
  }

  if (retain_seurat) {
    rm(mat)
  } else {
    rm(mat, obj)
  }
  invisible(gc(verbose = FALSE))

  info$expr_sub <- expr_sub
  info$coords <- coords
  info$load_seconds <- as.numeric(difftime(Sys.time(), start_ts, units = "secs"))
  log_detail("Load/inspect time (sec): ", fmt_num(info$load_seconds, digits = 2))
  log_detail("Status: ok")
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
      object_path = study_info$object_path,
      object_bytes = ifelse(is.na(study_info$object_bytes), "", as.character(as.numeric(study_info$object_bytes))),
      object_class = ifelse(is.na(study_info$object_class), "", study_info$object_class),
      default_assay = ifelse(is.na(study_info$default_assay), "", study_info$default_assay),
      assays_available = ifelse(is.na(study_info$assays_available), "", study_info$assays_available),
      assay_class = ifelse(is.na(study_info$assay_class), "", study_info$assay_class),
      assay_slot_names = ifelse(is.na(study_info$assay_slot_names), "", study_info$assay_slot_names),
      reductions_available = ifelse(is.na(study_info$reductions_available), "", study_info$reductions_available),
      reduction = study_info$reduction,
      assay = study_info$assay,
      assay_slot_requested = ifelse(is.na(study_info$assay_slot_requested), "", study_info$assay_slot_requested),
      expression_slot = ifelse(is.na(study_info$expression_slot), "", study_info$expression_slot),
      data_access_mode = ifelse(is.na(study_info$data_access_mode), "", study_info$data_access_mode),
      assay_slot_dims = ifelse(is.na(study_info$assay_slot_dims), "", study_info$assay_slot_dims),
      n_cells_object = study_info$n_cells_object,
      n_features_assay = study_info$n_features_assay,
      n_features_detected = study_info$n_features_detected,
      n_cells_assay = study_info$n_cells_assay,
      n_cells_umap = study_info$n_cells_umap,
      n_dims_umap = study_info$n_dims_umap,
      n_cells_common = study_info$n_cells_common,
      n_cells_umap_not_in_assay = study_info$n_cells_umap_not_in_assay,
      n_cells_assay_not_in_umap = study_info$n_cells_assay_not_in_umap,
      n_meta_rows = study_info$n_meta_rows,
      n_meta_cols = study_info$n_meta_cols,
      meta_rows_match_object_cells = study_info$meta_rows_match_object_cells,
      has_seurat_clusters = study_info$has_seurat_clusters,
      n_seurat_clusters = study_info$n_seurat_clusters,
      has_sample_id = study_info$has_sample_id,
      n_sample_id = study_info$n_sample_id,
      has_orig_ident = study_info$has_orig_ident,
      n_orig_ident = study_info$n_orig_ident,
      has_domain = study_info$has_domain,
      n_domain = study_info$n_domain,
      metadata_columns_preview = ifelse(is.na(study_info$metadata_columns_preview), "", study_info$metadata_columns_preview),
      n_ident_levels = study_info$n_ident_levels,
      n_ident_cells = study_info$n_ident_cells,
      umap_dim1_min = study_info$umap_dim1_min,
      umap_dim1_max = study_info$umap_dim1_max,
      umap_dim2_min = study_info$umap_dim2_min,
      umap_dim2_max = study_info$umap_dim2_max,
      n_marker_genes_requested = study_info$n_marker_genes_requested,
      n_marker_genes_present = study_info$n_marker_genes_present,
      marker_genes_present = ifelse(is.na(study_info$marker_genes_present), "", study_info$marker_genes_present),
      marker_genes_missing = ifelse(is.na(study_info$marker_genes_missing), "", study_info$marker_genes_missing),
      feature_id_type = ifelse(is.na(study_info$feature_id_type), "", study_info$feature_id_type),
      feature_id_examples = ifelse(is.na(study_info$feature_id_examples), "", study_info$feature_id_examples),
      runtime_note = ifelse(is.na(study_info$runtime_note), "", study_info$runtime_note),
      load_seconds = study_info$load_seconds,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, chunks)
}

build_marker_presence_table <- function(studies_info) {
  chunks <- list()
  k <- 0
  for (study_info in studies_info) {
    study_ok <- identical(study_info$status, "ok")
    genes_present <- if (!is.null(study_info$expr_sub)) rownames(study_info$expr_sub) else character(0)
    missing_gene_reason <- infer_missing_gene_reason(study_info)
    for (gene in GENE_ORDER) {
      k <- k + 1
      present <- study_ok && (gene %in% genes_present)
      reason <- if (!study_ok) {
        ifelse(is.na(study_info$reason), "Study unavailable", study_info$reason)
      } else if (!present) {
        missing_gene_reason
      } else {
        ""
      }
      chunks[[k]] <- data.frame(
        study_id = study_info$study_id,
        study_label = study_info$study_label,
        gene = gene,
        present = present,
        reason = reason,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, chunks)
}

build_assay_slot_summary_table <- function(studies_info) {
  empty <- data.frame(
    study_id = character(),
    study_label = character(),
    assay = character(),
    slot = character(),
    status = character(),
    matrix_class = character(),
    n_features = integer(),
    n_cells = integer(),
    detail = character(),
    stringsAsFactors = FALSE
  )
  out <- combine_table_field(studies_info, "assay_slot_summary", empty)
  if (nrow(out) == 0) {
    fallback <- lapply(studies_info, function(study_info) {
      data.frame(
        study_id = study_info$study_id,
        study_label = study_info$study_label,
        assay = study_info$assay,
        slot = study_info$preferred_slot,
        status = study_info$status,
        matrix_class = "",
        n_features = NA_integer_,
        n_cells = NA_integer_,
        detail = ifelse(is.na(study_info$reason), "", study_info$reason),
        stringsAsFactors = FALSE
      )
    })
    out <- do.call(rbind, fallback)
  }
  out
}

build_reduction_summary_table <- function(studies_info) {
  empty <- data.frame(
    study_id = character(),
    study_label = character(),
    reduction = character(),
    status = character(),
    n_cells = integer(),
    n_dims = integer(),
    dim1_min = numeric(),
    dim1_max = numeric(),
    dim2_min = numeric(),
    dim2_max = numeric(),
    detail = character(),
    stringsAsFactors = FALSE
  )
  out <- combine_table_field(studies_info, "reduction_summary", empty)
  if (nrow(out) == 0) {
    fallback <- lapply(studies_info, function(study_info) {
      data.frame(
        study_id = study_info$study_id,
        study_label = study_info$study_label,
        reduction = study_info$reduction,
        status = study_info$status,
        n_cells = NA_integer_,
        n_dims = NA_integer_,
        dim1_min = NA_real_,
        dim1_max = NA_real_,
        dim2_min = NA_real_,
        dim2_max = NA_real_,
        detail = ifelse(is.na(study_info$reason), "", study_info$reason),
        stringsAsFactors = FALSE
      )
    })
    out <- do.call(rbind, fallback)
  }
  out
}

build_metadata_summary_table <- function(studies_info) {
  rows <- lapply(studies_info, function(study_info) {
    if (!is.null(study_info$metadata_summary) && nrow(study_info$metadata_summary) > 0) {
      row <- study_info$metadata_summary[1, , drop = FALSE]
      row$status <- study_info$status
      return(row)
    }
    data.frame(
      study_id = study_info$study_id,
      study_label = study_info$study_label,
      status = study_info$status,
      n_meta_rows = NA_integer_,
      n_meta_cols = NA_integer_,
      meta_rows_match_object_cells = NA,
      has_seurat_clusters = NA,
      n_seurat_clusters = NA_integer_,
      has_sample_id = NA,
      n_sample_id = NA_integer_,
      has_orig_ident = NA,
      n_orig_ident = NA_integer_,
      has_domain = NA,
      n_domain = NA_integer_,
      metadata_columns_preview = "",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

build_metadata_columns_table <- function(studies_info) {
  empty <- data.frame(
    study_id = character(),
    study_label = character(),
    column_name = character(),
    column_class = character(),
    n_non_na = integer(),
    n_unique_non_na = integer(),
    example_values = character(),
    stringsAsFactors = FALSE
  )
  combine_table_field(studies_info, "metadata_columns", empty)
}

build_ident_counts_table <- function(studies_info) {
  empty <- data.frame(
    study_id = character(),
    study_label = character(),
    ident_value = character(),
    n_cells = integer(),
    stringsAsFactors = FALSE
  )
  combine_table_field(studies_info, "ident_counts", empty)
}

build_feature_space_table <- function(studies_info) {
  chunks <- lapply(studies_info, function(study_info) {
    data.frame(
      study_id = study_info$study_id,
      study_label = study_info$study_label,
      status = study_info$status,
      assay = study_info$assay,
      assay_class = ifelse(is.na(study_info$assay_class), "", study_info$assay_class),
      assay_slot_names = ifelse(is.na(study_info$assay_slot_names), "", study_info$assay_slot_names),
      n_features_detected = study_info$n_features_detected,
      feature_id_type = ifelse(is.na(study_info$feature_id_type), "", study_info$feature_id_type),
      feature_id_examples = ifelse(is.na(study_info$feature_id_examples), "", study_info$feature_id_examples),
      data_access_mode = ifelse(is.na(study_info$data_access_mode), "", study_info$data_access_mode),
      runtime_note = ifelse(is.na(study_info$runtime_note), "", study_info$runtime_note),
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
      reason <- if (expr_res$status == "missing_gene") infer_missing_gene_reason(study_info) else expr_res$reason
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
        scope = if (startsWith(reason, "Gene not found")) "gene" else "study",
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

  gene_level <- issue_df[grepl("^Gene not found", issue_df$reason), c("study_label", "gene"), drop = FALSE]
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

print_study_diagnostics <- function(study_status) {
  cat("\n=== Study diagnostics ===\n")
  for (i in seq_len(nrow(study_status))) {
    row <- study_status[i, , drop = FALSE]
    cat(
      paste0(
        "[", row$study_label[[1]], "] ",
        "status=", row$status[[1]],
        "; reason=", ifelse(nzchar(row$reason[[1]]), row$reason[[1]], "<none>"),
        "; file=", row$object_path[[1]], "\n"
      )
    )
    cat(
      paste0(
        "  class=", ifelse(nzchar(row$object_class[[1]]), row$object_class[[1]], "<NA>"),
        "; bytes=", ifelse(nzchar(row$object_bytes[[1]]), row$object_bytes[[1]], "NA"),
        "; default_assay=", ifelse(nzchar(row$default_assay[[1]]), row$default_assay[[1]], "<NA>"),
        "; load_seconds=", fmt_num(row$load_seconds[[1]], digits = 2), "\n"
      )
    )
    cat(
      paste0(
        "  assays=", ifelse(nzchar(row$assays_available[[1]]), row$assays_available[[1]], "<none>"),
        "; reductions=", ifelse(nzchar(row$reductions_available[[1]]), row$reductions_available[[1]], "<none>"), "\n"
      )
    )
    cat(
      paste0(
        "  requested assay/reduction/slot=", row$assay[[1]], "/", row$reduction[[1]], "/", row$assay_slot_requested[[1]],
        "; selected slot=", ifelse(nzchar(row$expression_slot[[1]]), row$expression_slot[[1]], "<NA>"),
        "; access_mode=", ifelse(nzchar(row$data_access_mode[[1]]), row$data_access_mode[[1]], "<NA>"), "\n"
      )
    )
    cat(
      paste0(
        "  assay class=", ifelse(nzchar(row$assay_class[[1]]), row$assay_class[[1]], "<NA>"),
        "; assay slots=", ifelse(nzchar(row$assay_slot_names[[1]]), row$assay_slot_names[[1]], "<NA>"),
        "; feature_id_type=", ifelse(nzchar(row$feature_id_type[[1]]), row$feature_id_type[[1]], "<NA>"), "\n"
      )
    )
    cat(
      paste0(
        "  cells: object=", fmt_count(row$n_cells_object[[1]]),
        ", assay=", fmt_count(row$n_cells_assay[[1]]),
        ", umap=", fmt_count(row$n_cells_umap[[1]]),
        ", metadata_rows=", fmt_count(row$n_meta_rows[[1]]),
        ", common=", fmt_count(row$n_cells_common[[1]]),
        ", umap_only=", fmt_count(row$n_cells_umap_not_in_assay[[1]]),
        ", assay_only=", fmt_count(row$n_cells_assay_not_in_umap[[1]]), "\n"
      )
    )
    cat(
      paste0(
        "  metadata: cols=", fmt_count(row$n_meta_cols[[1]]),
        ", rows_match_obj=", fmt_bool(row$meta_rows_match_object_cells[[1]]),
        ", seurat_clusters=", fmt_bool(row$has_seurat_clusters[[1]]), " (n=", fmt_count(row$n_seurat_clusters[[1]]), ")",
        ", sample_id=", fmt_bool(row$has_sample_id[[1]]), " (n=", fmt_count(row$n_sample_id[[1]]), ")",
        ", orig.ident=", fmt_bool(row$has_orig_ident[[1]]), " (n=", fmt_count(row$n_orig_ident[[1]]), ")",
        ", domain=", fmt_bool(row$has_domain[[1]]), " (n=", fmt_count(row$n_domain[[1]]), ")\n"
      )
    )
    if (nzchar(row$metadata_columns_preview[[1]])) {
      cat(paste0("  metadata columns preview=", row$metadata_columns_preview[[1]], "\n"))
    }
    cat(
      paste0(
        "  features: assay=", fmt_count(row$n_features_assay[[1]]),
        "; markers present=", fmt_count(row$n_marker_genes_present[[1]]),
        "/", fmt_count(row$n_marker_genes_requested[[1]]),
        "; ident_levels=", fmt_count(row$n_ident_levels[[1]]), "\n"
      )
    )
    if (nzchar(row$marker_genes_missing[[1]])) {
      cat(paste0("  missing markers=", row$marker_genes_missing[[1]], "\n"))
    }
    if (nzchar(row$detail[[1]])) {
      cat(paste0("  detail=", row$detail[[1]], "\n"))
    }
    if (nzchar(row$runtime_note[[1]])) {
      cat(paste0("  runtime_note=", row$runtime_note[[1]], "\n"))
    }
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
  detailed_log <- parse_bool_flag(args[["detailed-log"]], default = TRUE)
  set_detailed_log(detailed_log)
  options(warn = 1)

  log_msg("Run label: ", run_label)
  log_msg("Project root: ", project_root)
  log_msg("Detailed logging: ", ifelse(detailed_log, "enabled", "disabled"))
  log_msg("R version: ", R.version.string)
  log_msg("Seurat version: ", as.character(utils::packageVersion("Seurat")))
  log_msg("ggplot2 version: ", as.character(utils::packageVersion("ggplot2")))
  log_msg("patchwork version: ", as.character(utils::packageVersion("patchwork")))
  log_msg("Warnings: immediate (options(warn=1))")
  log_msg("No analysis recomputation is performed in this script (load existing objects only).")
  log_msg("Markers requested (", length(GENE_ORDER), "): ", collapse_csv(GENE_ORDER))
  log_msg("Studies in config (ordered): ", collapse_csv(cfg$studies$study_id))

  log_msg("Preparing studies from config: ", args$config)
  studies_info <- vector("list", nrow(cfg$studies))
  for (i in seq_len(nrow(cfg$studies))) {
    row <- cfg$studies[i, , drop = FALSE]
    log_msg("Study ", i, "/", nrow(cfg$studies), ": ", row$study_label[[1]], " (", row$study_id[[1]], ")")
    studies_info[[i]] <- prepare_study(row, project_root, retain_seurat = retain_seurat)
    log_msg(
      "Study ", row$study_label[[1]], " complete: status=", studies_info[[i]]$status,
      if (!is.na(studies_info[[i]]$reason)) paste0(" (", studies_info[[i]]$reason, ")") else ""
    )
  }
  ordered_labels <- vapply(studies_info, function(x) x$study_label, character(1))
  study_status <- build_study_status_table(studies_info)
  marker_presence <- build_marker_presence_table(studies_info)
  assay_slot_summary <- build_assay_slot_summary_table(studies_info)
  reduction_summary <- build_reduction_summary_table(studies_info)
  metadata_summary <- build_metadata_summary_table(studies_info)
  metadata_columns <- build_metadata_columns_table(studies_info)
  ident_counts <- build_ident_counts_table(studies_info)
  feature_space <- build_feature_space_table(studies_info)

  out_dir <- normalize_abs(file.path(project_root, "results", run_label, "plots"), must_work = FALSE)
  if (!is_subpath(out_dir, project_root)) {
    stop("Output path must remain under PROJECT_ROOT.", call. = FALSE)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  study_status_path <- file.path(out_dir, "panel_b_study_status.tsv")
  marker_presence_path <- file.path(out_dir, "panel_b_marker_presence.tsv")
  assay_slot_summary_path <- file.path(out_dir, "panel_b_assay_slot_summary.tsv")
  reduction_summary_path <- file.path(out_dir, "panel_b_reduction_summary.tsv")
  metadata_summary_path <- file.path(out_dir, "panel_b_metadata_summary.tsv")
  metadata_columns_path <- file.path(out_dir, "panel_b_metadata_columns.tsv")
  ident_counts_path <- file.path(out_dir, "panel_b_ident_counts.tsv")
  feature_space_path <- file.path(out_dir, "panel_b_feature_space.tsv")

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
    log_detail(
      "Row ", i, "/", length(GENE_ORDER), " ", gene,
      " | points=", fmt_count(row_res$summary$n_points[[1]]),
      " | placeholders=", fmt_count(row_res$summary$n_placeholder_panels[[1]]),
      " | studies_with_expression=", fmt_count(row_res$summary$n_studies_with_expression[[1]]),
      " | scale=[", fmt_num(row_res$summary$scale_min[[1]]), ", ", fmt_num(row_res$summary$scale_max[[1]]), "]"
    )
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
  row_summary_path <- file.path(out_dir, "panel_b_row_summary.tsv")
  issues_path <- file.path(out_dir, "panel_b_issues.tsv")

  utils::write.table(study_status, file = study_status_path, sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(marker_presence, file = marker_presence_path, sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(assay_slot_summary, file = assay_slot_summary_path, sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(reduction_summary, file = reduction_summary_path, sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(metadata_summary, file = metadata_summary_path, sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(metadata_columns, file = metadata_columns_path, sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(ident_counts, file = ident_counts_path, sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(feature_space, file = feature_space_path, sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(row_summary, file = row_summary_path, sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(issues, file = issues_path, sep = "\t", quote = FALSE, row.names = FALSE)
  log_msg("Wrote study diagnostics table: ", study_status_path)
  log_msg("Wrote marker presence table: ", marker_presence_path)
  log_msg("Wrote assay-slot summary table: ", assay_slot_summary_path)
  log_msg("Wrote reduction summary table: ", reduction_summary_path)
  log_msg("Wrote metadata summary table: ", metadata_summary_path)
  log_msg("Wrote metadata columns table: ", metadata_columns_path)
  log_msg("Wrote identity counts table: ", ident_counts_path)
  log_msg("Wrote feature-space table: ", feature_space_path)
  log_msg("Wrote row summary table: ", row_summary_path)
  log_msg("Wrote issue table: ", issues_path)

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
  print_study_diagnostics(study_status)
  print_audit(issues)

  result <- list(
    project_root = project_root,
    run_label = run_label,
    genes = GENE_ORDER,
    study_status = study_status,
    marker_presence = marker_presence,
    assay_slot_summary = assay_slot_summary,
    reduction_summary = reduction_summary,
    metadata_summary = metadata_summary,
    metadata_columns = metadata_columns,
    ident_counts = ident_counts,
    feature_space = feature_space,
    row_summary = row_summary,
    issues = issues,
    output_paths = list(
      pdf = pdf_path,
      svg = svg_path,
      study_status = study_status_path,
      marker_presence = marker_presence_path,
      assay_slot_summary = assay_slot_summary_path,
      reduction_summary = reduction_summary_path,
      metadata_summary = metadata_summary_path,
      metadata_columns = metadata_columns_path,
      ident_counts = ident_counts_path,
      feature_space = feature_space_path,
      row_summary = row_summary_path,
      issues = issues_path
    ),
    row_plots = row_plots,
    final_plot = fig,
    studies_info = studies_info
  )

  if (export_global) {
    assign("panel_b_result", result, envir = .GlobalEnv)
    assign("panel_b_studies", study_status, envir = .GlobalEnv)
    assign("panel_b_marker_presence", marker_presence, envir = .GlobalEnv)
    assign("panel_b_assay_slot_summary", assay_slot_summary, envir = .GlobalEnv)
    assign("panel_b_reduction_summary", reduction_summary, envir = .GlobalEnv)
    assign("panel_b_metadata_summary", metadata_summary, envir = .GlobalEnv)
    assign("panel_b_metadata_columns", metadata_columns, envir = .GlobalEnv)
    assign("panel_b_ident_counts", ident_counts, envir = .GlobalEnv)
    assign("panel_b_feature_space", feature_space, envir = .GlobalEnv)
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
                              show_progress_plots = interactive(),
                              detailed_log = TRUE) {
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
    "--show-progress", if (isTRUE(show_progress_plots)) "true" else "false",
    "--detailed-log", if (isTRUE(detailed_log)) "true" else "false"
  )
  main(args)
}

if (sys.nframe() == 0) {
  main()
}
