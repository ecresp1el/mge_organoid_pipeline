#!/usr/bin/env Rscript

# Stage 07 placeholder for cross-study comparison against human developing GE.
#
# This script does not yet run the biological comparison. It prepares and audits
# validated per-study Seurat object paths exported by Stage 06 so downstream
# development can consume a stable manifest.
#
# Inputs:
# - Config file (`--config`) with:
#   - run_label
#   - optional project_root
#   - prepared_object_manifest_path
#   - optional required_studies
#   - optional human_ge_reference_path (placeholder)
#
# Outputs:
# - PROJECT_ROOT/results/<run_label>/human_ge_comparison_tbd/human_ge_input_manifest.tsv
# - PROJECT_ROOT/results/<run_label>/human_ge_comparison_tbd/human_ge_missing_required_studies.tsv
# - PROJECT_ROOT/results/<run_label>/human_ge_comparison_tbd/human_ge_comparison_tbd_notes.txt

log_msg <- function(...) {
  msg <- paste0(..., collapse = " ")
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
  flush.console()
}

print_usage <- function() {
  cat(
    paste(
      "Usage:",
      "  Rscript scripts/07_compare_human_developing_ge_tbd.R --config <config.R> [--project-root <path>] [--run-label <label>] [--manifest-path <path>]",
      "",
      "Config file format:",
      "  list(",
      "    project_root = <optional>,",
      "    run_label = <string>,",
      "    prepared_object_manifest_path = <string>,",
      "    prepared_object_studies_dir = <optional string>,",
      "    required_studies = <optional character vector>,",
      "    human_ge_reference_path = <optional string>",
      "  )",
      "",
      "Outputs:",
      "  PROJECT_ROOT/results/<run_label>/human_ge_comparison_tbd/human_ge_input_manifest.tsv",
      "  PROJECT_ROOT/results/<run_label>/human_ge_comparison_tbd/human_ge_missing_required_studies.tsv",
      "  PROJECT_ROOT/results/<run_label>/human_ge_comparison_tbd/human_ge_comparison_tbd_notes.txt",
      sep = "\n"
    )
  )
}

parse_args <- function(args) {
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
      out[[kv[[1]]]] <- if (length(kv) > 1) paste(kv[-1], collapse = "=") else ""
      i <- i + 1
      next
    }
    if (i == length(args)) stop("Missing value for --", token, call. = FALSE)
    out[[token]] <- args[[i + 1]]
    i <- i + 2
  }
  out
}

trim_trailing_slash <- function(x) sub("/+$", "", x)
normalize_abs <- function(path, must_work = FALSE) normalizePath(path, winslash = "/", mustWork = must_work)

is_subpath <- function(path, root) {
  path_norm <- normalize_abs(path, must_work = FALSE)
  root_norm <- normalize_abs(root, must_work = FALSE)
  if (identical(path_norm, root_norm)) return(TRUE)
  startsWith(path_norm, paste0(root_norm, "/"))
}

resolve_under_project_root <- function(path_value, project_root) {
  if (!is.character(path_value) || length(path_value) != 1 || !nzchar(path_value)) {
    stop("Invalid path value: must be non-empty string", call. = FALSE)
  }
  candidate <- if (startsWith(path_value, "/")) path_value else file.path(project_root, path_value)
  candidate <- normalize_abs(candidate, must_work = FALSE)
  if (!is_subpath(candidate, project_root)) {
    stop("Path must remain under PROJECT_ROOT: ", candidate, call. = FALSE)
  }
  candidate
}

read_config <- function(config_path) {
  cfg <- tryCatch(
    dget(config_path),
    error = function(e) stop("Failed to read config via dget(): ", conditionMessage(e), call. = FALSE)
  )
  if (!is.list(cfg)) stop("Config must be a list.", call. = FALSE)
  if (is.null(cfg$run_label) || !nzchar(as.character(cfg$run_label))) {
    stop("Config missing run_label.", call. = FALSE)
  }
  if (is.null(cfg$prepared_object_manifest_path) || !nzchar(as.character(cfg$prepared_object_manifest_path))) {
    stop("Config missing prepared_object_manifest_path.", call. = FALSE)
  }

  cfg$run_label <- as.character(cfg$run_label)
  if (!is.null(cfg$project_root)) cfg$project_root <- as.character(cfg$project_root)
  cfg$prepared_object_manifest_path <- as.character(cfg$prepared_object_manifest_path)
  if (!is.null(cfg$prepared_object_studies_dir)) {
    cfg$prepared_object_studies_dir <- as.character(cfg$prepared_object_studies_dir)
  } else {
    cfg$prepared_object_studies_dir <- "results/panel_b_prepared_objects/studies"
  }
  if (!is.null(cfg$required_studies)) {
    cfg$required_studies <- unique(as.character(cfg$required_studies))
  } else {
    cfg$required_studies <- character(0)
  }
  if (!is.null(cfg$human_ge_reference_path)) {
    cfg$human_ge_reference_path <- as.character(cfg$human_ge_reference_path)
  } else {
    cfg$human_ge_reference_path <- ""
  }
  cfg
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

  manifest_path <- args[["manifest-path"]]
  if (is.null(manifest_path) || !nzchar(manifest_path)) manifest_path <- cfg$prepared_object_manifest_path
  manifest_path <- resolve_under_project_root(manifest_path, project_root)
  if (!file.exists(manifest_path)) {
    stop("Prepared object manifest not found: ", manifest_path, call. = FALSE)
  }

  studies_dir <- resolve_under_project_root(cfg$prepared_object_studies_dir, project_root)
  out_dir <- normalize_abs(file.path(project_root, "results", run_label, "human_ge_comparison_tbd"), must_work = FALSE)
  if (!is_subpath(out_dir, project_root)) {
    stop("Output path must remain under PROJECT_ROOT.", call. = FALSE)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  log_msg("Project root: ", project_root)
  log_msg("Run label: ", run_label)
  log_msg("Prepared object manifest: ", manifest_path)
  log_msg("Prepared object studies dir: ", studies_dir)

  manifest <- utils::read.table(manifest_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  required_cols <- c("study_id", "study_label", "status", "publish_status", "prepared_object_path")
  missing_cols <- setdiff(required_cols, colnames(manifest))
  if (length(missing_cols) > 0) {
    stop("Prepared manifest missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  manifest$prepared_object_path <- as.character(manifest$prepared_object_path)
  manifest$file_exists <- nzchar(manifest$prepared_object_path) & file.exists(manifest$prepared_object_path)

  eligible <- manifest[
    manifest$status == "ok" &
      manifest$publish_status == "published" &
      manifest$file_exists,
    ,
    drop = FALSE
  ]
  if (nrow(eligible) == 0) {
    stop("No eligible prepared Seurat objects found in manifest.", call. = FALSE)
  }

  required_ids <- cfg$required_studies
  missing_required <- character(0)
  if (length(required_ids) > 0) {
    missing_required <- setdiff(required_ids, eligible$study_id)
  }
  missing_df <- if (length(missing_required) > 0) {
    data.frame(
      missing_study_id = missing_required,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      missing_study_id = character(),
      stringsAsFactors = FALSE
    )
  }

  eligible_path <- file.path(out_dir, "human_ge_input_manifest.tsv")
  missing_path <- file.path(out_dir, "human_ge_missing_required_studies.tsv")
  notes_path <- file.path(out_dir, "human_ge_comparison_tbd_notes.txt")

  utils::write.table(eligible, file = eligible_path, sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(missing_df, file = missing_path, sep = "\t", quote = FALSE, row.names = FALSE)

  notes_lines <- c(
    "Human Developing GE Comparison: TBD Placeholder",
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("Source manifest: ", manifest_path),
    paste0("Eligible prepared studies: ", nrow(eligible)),
    paste0("Required studies: ", if (length(required_ids) > 0) paste(required_ids, collapse = ",") else "<none specified>"),
    paste0("Missing required studies: ", if (length(missing_required) > 0) paste(missing_required, collapse = ",") else "<none>"),
    paste0("Human GE reference path (placeholder): ", if (nzchar(cfg$human_ge_reference_path)) cfg$human_ge_reference_path else "<unset>"),
    "",
    "TBD Next Development Tasks",
    "1. Load prepared study objects from human_ge_input_manifest.tsv.",
    "2. Load/validate human developing GE reference data.",
    "3. Harmonize gene namespace across studies and reference (symbols/Ensembl).",
    "4. Define and implement comparison metrics and visualization outputs.",
    "5. Emit methods-locked outputs under this run_label for reproducibility."
  )
  writeLines(notes_lines, con = notes_path, useBytes = TRUE)

  log_msg("Eligible prepared study objects: ", nrow(eligible))
  if (length(missing_required) > 0) {
    log_msg("Missing required studies: ", paste(missing_required, collapse = ","))
  } else {
    log_msg("All required studies are present.")
  }
  log_msg("Wrote human GE input manifest: ", eligible_path)
  log_msg("Wrote missing-required table: ", missing_path)
  log_msg("Wrote placeholder notes: ", notes_path)

  invisible(list(
    project_root = project_root,
    run_label = run_label,
    manifest_path = manifest_path,
    eligible_path = eligible_path,
    missing_path = missing_path,
    notes_path = notes_path,
    n_eligible = nrow(eligible),
    missing_required = missing_required
  ))
}

if (sys.nframe() == 0) {
  main()
}
