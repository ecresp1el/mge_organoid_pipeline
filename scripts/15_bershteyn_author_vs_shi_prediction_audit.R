#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_msg <- function(...) message("[", timestamp(), "] ", paste0(..., collapse = ""))

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    project_root = Sys.getenv("PROJECT_ROOT", "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder"),
    run_label = Sys.getenv("CROSS_STUDY_SHI_RUN_LABEL", "cross_study_shi_seurat_label_transfer_v2_ge_only_age"),
    seurat = "",
    obs = "",
    outdir = ""
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    value <- if (i < length(args)) args[[i + 1L]] else ""
    if (key == "--project-root") {
      out$project_root <- value; i <- i + 2L
    } else if (key == "--run-label") {
      out$run_label <- value; i <- i + 2L
    } else if (key == "--seurat") {
      out$seurat <- value; i <- i + 2L
    } else if (key == "--obs") {
      out$obs <- value; i <- i + 2L
    } else if (key == "--outdir") {
      out$outdir <- value; i <- i + 2L
    } else {
      stop("Unknown argument: ", key, call. = FALSE)
    }
  }
  run_root <- file.path(out$project_root, "results/cross_study_shi_seurat_label_transfer", out$run_label)
  if (!nzchar(out$seurat)) {
    out$seurat <- file.path(out$project_root, "results/bershteyn_2025/bershteyn_2025_seurat.rds")
  }
  if (!nzchar(out$obs)) {
    out$obs <- file.path(run_root, "tables/per_study/bershteyn_2025_shi_seurat_label_transfer_obs.tsv.gz")
  }
  if (!nzchar(out$outdir)) {
    out$outdir <- file.path(run_root, "tables/author_comparison")
  }
  out
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

write_tsv_gz <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

score_inventory_rows <- function(location, names, class, length_value) {
  if (length(names) == 0L) {
    return(data.frame(location = character(), name = character(), class = character(), length = integer()))
  }
  data.frame(
    location = location,
    name = names,
    class = class,
    length = length_value,
    stringsAsFactors = FALSE
  )
}

read_tsv_maybe_gz <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  utils::read.delim(con, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

normalize_gw_label <- function(x) {
  text <- as.character(x)
  out <- rep(NA_character_, length(text))
  hit <- regexpr("GW\\s*0*([0-9]+)", text, ignore.case = TRUE, perl = TRUE)
  has_hit <- hit > 0
  if (any(has_hit)) {
    matched <- regmatches(text, hit)
    nums <- sub(".*GW\\s*0*([0-9]+).*", "\\1", matched, ignore.case = TRUE, perl = TRUE)
    out[has_hit] <- sprintf("GW%02d", as.integer(nums))
  }
  out[!has_hit & nzchar(text)] <- text[!has_hit & nzchar(text)]
  out
}

parse_gw_numeric_safe <- function(x) {
  label <- normalize_gw_label(x)
  out <- rep(NA_real_, length(label))
  hit <- regexpr("[0-9]+", label, perl = TRUE)
  has_hit <- hit > 0
  out[has_hit] <- as.numeric(regmatches(label, hit)[has_hit])
  out
}

expected_gw_from_score_cols <- function(df, cols) {
  labels <- sub("^.*_prediction_score_", "", cols)
  gw_numeric <- parse_gw_numeric_safe(labels)
  if (any(!is.finite(gw_numeric))) {
    stop("Could not parse GW labels from score columns: ", paste(cols[!is.finite(gw_numeric)], collapse = ", "), call. = FALSE)
  }
  mat <- as.matrix(df[, cols, drop = FALSE])
  storage.mode(mat) <- "double"
  denom <- rowSums(mat, na.rm = TRUE)
  expected <- as.numeric(mat %*% gw_numeric) / denom
  expected[!is.finite(expected)] <- NA_real_
  expected
}

count_fraction <- function(values, source, total) {
  values <- as.character(values)
  values[is.na(values)] <- "<NA>"
  tab <- sort(table(values), decreasing = TRUE)
  data.frame(
    source = source,
    label = names(tab),
    n_cells = as.integer(tab),
    fraction = as.numeric(tab) / total,
    stringsAsFactors = FALSE
  )
}

confusion_table <- function(author, ours, table_name) {
  author <- as.character(author)
  ours <- as.character(ours)
  author[is.na(author)] <- "<NA>"
  ours[is.na(ours)] <- "<NA>"
  tab <- as.data.frame(table(author_label = author, our_label = ours), stringsAsFactors = FALSE)
  tab <- tab[tab$Freq > 0L, , drop = FALSE]
  colnames(tab)[3] <- "n_cells"
  totals <- stats::setNames(as.integer(table(author)), names(table(author)))
  tab$author_total <- totals[tab$author_label]
  tab$fraction_within_author_label <- tab$n_cells / tab$author_total
  tab$table_name <- table_name
  tab <- tab[order(tab$author_label, -tab$n_cells, tab$our_label), c("table_name", "author_label", "our_label", "n_cells", "author_total", "fraction_within_author_label")]
  rownames(tab) <- NULL
  tab
}

label_count <- function(df, label) {
  hit <- df[df$label == label, , drop = FALSE]
  if (nrow(hit) == 0L) return(list(n = 0L, fraction = 0))
  list(n = hit$n_cells[[1]], fraction = hit$fraction[[1]])
}

main <- function() {
  opt <- parse_args()
  if (!file.exists(opt$seurat)) stop("Missing Seurat object: ", opt$seurat, call. = FALSE)
  if (!file.exists(opt$obs)) stop("Missing Shi transfer obs table: ", opt$obs, call. = FALSE)
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

  log_msg("Reading Bershteyn Seurat object: ", opt$seurat)
  obj <- readRDS(opt$seurat)
  if (!inherits(obj, "Seurat")) stop("Object is not a Seurat object: ", opt$seurat, call. = FALSE)
  meta <- obj@meta.data
  meta$cell_id <- rownames(meta)
  author_cols <- c(
    "cell_id", "sample", "seurat_clusters", "predicted.GEtype",
    "predicted.GEcluster", "predicted.GEgws", "predicted.macaclass",
    "predicted.musclass", "type"
  )
  keep <- intersect(author_cols, colnames(meta))
  author_meta <- meta[, keep, drop = FALSE]
  log_msg("Author metadata rows: ", nrow(author_meta), "; columns: ", paste(keep, collapse = ","))

  score_patterns <- "(score|prediction\\.score|predicted|GEgws|GEtype|macaclass|musclass)"
  score_like_rows <- list(
    score_inventory_rows(
      location = "meta.data",
      names = grep(score_patterns, colnames(meta), value = TRUE, ignore.case = TRUE),
      class = "metadata_column",
      length_value = nrow(meta)
    )
  )
  for (slot_name in c("misc", "tools", "commands")) {
    value <- tryCatch(slot(obj, slot_name), error = function(e) NULL)
    if (is.null(value) || length(value) == 0L) next
    value_names <- names(value)
    if (is.null(value_names)) value_names <- rep("", length(value))
    hits <- grep(score_patterns, value_names, value = TRUE, ignore.case = TRUE)
    if (length(hits) > 0L) {
      score_like_rows[[length(score_like_rows) + 1L]] <- data.frame(
        location = slot_name,
        name = hits,
        class = vapply(value[hits], function(x) paste(class(x), collapse = ","), character(1)),
        length = vapply(value[hits], length, integer(1)),
        stringsAsFactors = FALSE
      )
    }
  }
  reductions <- names(obj@reductions)
  assays <- names(obj@assays)
  score_like_rows[[length(score_like_rows) + 1L]] <- score_inventory_rows(
    location = "reductions",
    names = grep(score_patterns, reductions, value = TRUE, ignore.case = TRUE),
    class = "reduction_name",
    length_value = NA_integer_
  )
  score_like_rows[[length(score_like_rows) + 1L]] <- score_inventory_rows(
    location = "assays",
    names = grep(score_patterns, assays, value = TRUE, ignore.case = TRUE),
    class = "assay_name",
    length_value = NA_integer_
  )
  score_like_inventory <- do.call(rbind, score_like_rows)
  score_like_inventory <- score_like_inventory[nzchar(score_like_inventory$name), , drop = FALSE]

  log_msg("Reading our Shi transfer table: ", opt$obs)
  obs <- read_tsv_maybe_gz(opt$obs)
  obs$cell_id <- as.character(obs$cell_id)
  author_meta$cell_id <- as.character(author_meta$cell_id)

  log_msg("Joining author metadata to our per-cell transfer table")
  joined <- merge(
    obs,
    author_meta,
    by = "cell_id",
    all.x = TRUE,
    suffixes = c(".our", ".author"),
    sort = FALSE
  )
  if (nrow(joined) != nrow(obs)) {
    stop("Join changed row count: obs=", nrow(obs), " joined=", nrow(joined), call. = FALSE)
  }
  n_cells <- nrow(joined)
  n_author_missing <- sum(is.na(joined$predicted.GEtype))

  joined$author_predicted_GEgws_normalized <- normalize_gw_label(joined$predicted.GEgws)
  joined$our_full_week_normalized <- normalize_gw_label(joined$shi_seurat_full_predicted_shi_week_label)
  joined$our_ge_only_week_normalized <- normalize_gw_label(joined$shi_seurat_ge_only_predicted_shi_week_label)

  full_week_cols <- grep("^shi_seurat_full_week_prediction_score_GW", colnames(joined), value = TRUE)
  ge_week_cols <- grep("^shi_seurat_ge_only_week_prediction_score_GW", colnames(joined), value = TRUE)
  joined$shi_seurat_full_expected_shi_gw_numeric_corrected <- expected_gw_from_score_cols(joined, full_week_cols)
  joined$shi_seurat_ge_only_expected_shi_gw_numeric_corrected <- expected_gw_from_score_cols(joined, ge_week_cols)

  author_ge <- count_fraction(joined$predicted.GEtype, "author_predicted_GEtype", n_cells)
  author_gw <- count_fraction(joined$author_predicted_GEgws_normalized, "author_predicted_GEgws_normalized", n_cells)
  author_type <- count_fraction(joined$type, "author_type", n_cells)
  our_major <- count_fraction(joined$shi_seurat_full_predicted_shi_label, "our_shi_full_major_label", n_cells)
  our_full_week <- count_fraction(joined$our_full_week_normalized, "our_shi_full_week_label_normalized", n_cells)
  our_ge_week <- count_fraction(joined$our_ge_only_week_normalized, "our_shi_ge_only_week_label_normalized", n_cells)
  distributions <- rbind(author_ge, author_gw, author_type, our_major, our_full_week, our_ge_week)

  get_n <- function(dist, label) label_count(dist, label)$n
  get_frac <- function(dist, label) label_count(dist, label)$fraction
  coarse <- rbind(
    do.call(rbind, lapply(c("MGE", "LGE", "CGE"), function(label) {
      data.frame(
        comparison = label,
        author_column = "predicted.GEtype",
        author_n = get_n(author_ge, label),
        author_fraction = get_frac(author_ge, label),
        our_column = "shi_seurat_full_predicted_shi_label",
        our_n = get_n(our_major, label),
        our_fraction = get_frac(our_major, label),
        fraction_delta_our_minus_author = get_frac(our_major, label) - get_frac(author_ge, label),
        stringsAsFactors = FALSE
      )
    })),
    do.call(rbind, lapply(c("GW09", "GW12", "GW13", "GW16", "GW18"), function(label) {
      data.frame(
        comparison = label,
        author_column = "predicted.GEgws",
        author_n = get_n(author_gw, label),
        author_fraction = get_frac(author_gw, label),
        our_column = "normalized shi_seurat_full_predicted_shi_week_label",
        our_n = get_n(our_full_week, label),
        our_fraction = get_frac(our_full_week, label),
        fraction_delta_our_minus_author = get_frac(our_full_week, label) - get_frac(author_gw, label),
        stringsAsFactors = FALSE
      )
    }))
  )

  confusion <- rbind(
    confusion_table(joined$predicted.GEtype, joined$shi_seurat_full_predicted_shi_label, "author_GEtype_vs_our_full_major"),
    confusion_table(joined$author_predicted_GEgws_normalized, joined$our_full_week_normalized, "author_GEgws_vs_our_full_week"),
    confusion_table(joined$author_predicted_GEgws_normalized, joined$our_ge_only_week_normalized, "author_GEgws_vs_our_ge_only_week"),
    confusion_table(joined$type, joined$shi_seurat_full_predicted_shi_label, "author_type_vs_our_full_major")
  )

  confidence_major <- aggregate(
    shi_seurat_full_prediction_score ~ shi_seurat_full_predicted_shi_label,
    data = joined,
    FUN = function(x) c(mean = mean(x, na.rm = TRUE), median = stats::median(x, na.rm = TRUE))
  )
  confidence_major <- data.frame(
    label = confidence_major$shi_seurat_full_predicted_shi_label,
    mean_max_score = confidence_major$shi_seurat_full_prediction_score[, "mean"],
    median_max_score = confidence_major$shi_seurat_full_prediction_score[, "median"],
    stringsAsFactors = FALSE
  )
  confidence_major <- merge(
    our_major[, c("label", "n_cells", "fraction")],
    confidence_major,
    by = "label",
    all.x = TRUE,
    sort = FALSE
  )

  confidence_week <- aggregate(
    cbind(shi_seurat_full_week_prediction_score, shi_seurat_full_expected_shi_gw_numeric_corrected) ~ our_full_week_normalized,
    data = joined,
    FUN = function(x) c(mean = mean(x, na.rm = TRUE), median = stats::median(x, na.rm = TRUE))
  )
  confidence_week <- data.frame(
    label = confidence_week$our_full_week_normalized,
    mean_max_score = confidence_week$shi_seurat_full_week_prediction_score[, "mean"],
    median_max_score = confidence_week$shi_seurat_full_week_prediction_score[, "median"],
    mean_expected_gw_corrected = confidence_week$shi_seurat_full_expected_shi_gw_numeric_corrected[, "mean"],
    median_expected_gw_corrected = confidence_week$shi_seurat_full_expected_shi_gw_numeric_corrected[, "median"],
    stringsAsFactors = FALSE
  )
  confidence_week <- merge(
    our_full_week[, c("label", "n_cells", "fraction")],
    confidence_week,
    by = "label",
    all.x = TRUE,
    sort = FALSE
  )

  confidence_ge_week <- aggregate(
    cbind(shi_seurat_ge_only_week_prediction_score, shi_seurat_ge_only_expected_shi_gw_numeric_corrected) ~ our_ge_only_week_normalized,
    data = joined,
    FUN = function(x) c(mean = mean(x, na.rm = TRUE), median = stats::median(x, na.rm = TRUE))
  )
  confidence_ge_week <- data.frame(
    label = confidence_ge_week$our_ge_only_week_normalized,
    mean_max_score = confidence_ge_week$shi_seurat_ge_only_week_prediction_score[, "mean"],
    median_max_score = confidence_ge_week$shi_seurat_ge_only_week_prediction_score[, "median"],
    mean_expected_gw_corrected = confidence_ge_week$shi_seurat_ge_only_expected_shi_gw_numeric_corrected[, "mean"],
    median_expected_gw_corrected = confidence_ge_week$shi_seurat_ge_only_expected_shi_gw_numeric_corrected[, "median"],
    stringsAsFactors = FALSE
  )
  confidence_ge_week <- merge(
    our_ge_week[, c("label", "n_cells", "fraction")],
    confidence_ge_week,
    by = "label",
    all.x = TRUE,
    sort = FALSE
  )

  sample_week <- as.data.frame(table(
    sample = joined$sample.our,
    shi_seurat_full_predicted_shi_week_label_normalized = joined$our_full_week_normalized
  ), stringsAsFactors = FALSE)
  colnames(sample_week)[3] <- "n_cells"
  sample_week <- sample_week[sample_week$n_cells > 0L, , drop = FALSE]
  sample_totals <- stats::setNames(as.integer(table(joined$sample.our)), names(table(joined$sample.our)))
  sample_week$sample_total <- sample_totals[sample_week$sample]
  sample_week$fraction <- sample_week$n_cells / sample_week$sample_total

  summary <- data.frame(
    metric = c(
      "n_cells",
      "n_cells_missing_author_metadata_after_join",
      "author_GEtype_MGE_fraction",
      "our_full_major_MGE_fraction",
      "author_GEgws_GW18_fraction",
      "our_full_week_GW18_fraction_normalized",
      "our_ge_only_week_GW18_fraction_normalized",
      "our_full_expected_gw_mean_corrected_from_scores",
      "our_ge_only_expected_gw_mean_corrected_from_scores",
      "per_study_obs_full_expected_gw_mean_stale_buggy",
      "per_study_obs_ge_only_expected_gw_mean_stale_buggy",
      "median_major_prediction_score",
      "fraction_major_score_ge_0_75",
      "median_ge_only_week_prediction_score"
    ),
    value = c(
      n_cells,
      n_author_missing,
      get_frac(author_ge, "MGE"),
      get_frac(our_major, "MGE"),
      get_frac(author_gw, "GW18"),
      get_frac(our_full_week, "GW18"),
      get_frac(our_ge_week, "GW18"),
      mean(joined$shi_seurat_full_expected_shi_gw_numeric_corrected, na.rm = TRUE),
      mean(joined$shi_seurat_ge_only_expected_shi_gw_numeric_corrected, na.rm = TRUE),
      mean(joined$shi_seurat_full_expected_shi_gw_numeric, na.rm = TRUE),
      mean(joined$shi_seurat_ge_only_expected_shi_gw_numeric, na.rm = TRUE),
      stats::median(joined$shi_seurat_full_prediction_score, na.rm = TRUE),
      mean(joined$shi_seurat_full_prediction_score >= 0.75, na.rm = TRUE),
      stats::median(joined$shi_seurat_ge_only_week_prediction_score, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )

  write_tsv(summary, file.path(opt$outdir, "bershteyn_2025_author_vs_our_shi_prediction_summary.tsv"))
  write_tsv(score_like_inventory, file.path(opt$outdir, "bershteyn_2025_author_score_like_object_inventory.tsv"))
  write_tsv(distributions, file.path(opt$outdir, "bershteyn_2025_author_vs_our_shi_prediction_distributions.tsv"))
  write_tsv(coarse, file.path(opt$outdir, "bershteyn_2025_author_vs_our_shi_prediction_coarse_comparison.tsv"))
  write_tsv(confusion, file.path(opt$outdir, "bershteyn_2025_author_vs_our_shi_prediction_cell_level_confusion.tsv"))
  write_tsv(confidence_major, file.path(opt$outdir, "bershteyn_2025_our_shi_major_confidence_by_label.tsv"))
  write_tsv(confidence_week, file.path(opt$outdir, "bershteyn_2025_our_shi_full_week_confidence_by_label.tsv"))
  write_tsv(confidence_ge_week, file.path(opt$outdir, "bershteyn_2025_our_shi_ge_only_week_confidence_by_label.tsv"))
  write_tsv(sample_week, file.path(opt$outdir, "bershteyn_2025_our_shi_full_week_by_sample.tsv"))

  joined_export_cols <- c(
    "cell_id", "sample.our", "cluster", "predicted.GEtype", "predicted.GEcluster",
    "author_predicted_GEgws_normalized", "predicted.macaclass", "predicted.musclass",
    "type", "shi_seurat_full_predicted_shi_label", "shi_seurat_full_prediction_score",
    "our_full_week_normalized", "shi_seurat_full_week_prediction_score",
    "shi_seurat_full_expected_shi_gw_numeric_corrected", "our_ge_only_week_normalized",
    "shi_seurat_ge_only_week_prediction_score", "shi_seurat_ge_only_expected_shi_gw_numeric_corrected"
  )
  joined_export_cols <- intersect(joined_export_cols, colnames(joined))
  write_tsv_gz(joined[, joined_export_cols, drop = FALSE], file.path(opt$outdir, "bershteyn_2025_author_vs_our_shi_prediction_joined_cells.tsv.gz"))

  write_tsv(
    data.frame(
      note = c(
        "Author columns are read from the downloaded Bershteyn Seurat object metadata.",
        "Our columns are read from the cross-study Shi Seurat TransferData v2 per-study obs table.",
        "Prediction inputs for our transfer are RNA expression only; author labels are used here only for post-hoc audit.",
        "The per-study obs table contains stale/buggy expected-GW numeric columns from the R transfer step because duplicate labels such as GW12_01 were parsed by averaging all numbers. Corrected expected-GW values are recomputed here from score columns using the GW prefix only."
      ),
      stringsAsFactors = FALSE
    ),
    file.path(opt$outdir, "bershteyn_2025_author_vs_our_shi_prediction_notes.tsv")
  )

  log_msg("Wrote author-vs-our audit tables to: ", opt$outdir)
  print(summary)
}

main()
