#!/usr/bin/env Rscript

# Identify Jia-like proliferative VZ-RGC root candidates from a score-defined
# root pool.
#
# Per-cell score:
#   RootScore = z(jia_score_RGC1) + z(proliferation_score) - z(jia_score_IPC)
#
# The proliferation score is the mean logUPX expression of available
# proliferation genes. Candidate roots are the top 1%, 2%, 5%, and 10% cells by
# RootScore, either across all cells or within an explicit root pool such as
# paper/manual Radial glia. This script reports candidate marker expression and
# candidate placement on existing UMAP and URD pseudotime, then writes a scored
# URD object with logical root-candidate columns.

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = " ")))
  flush.console()
}

parse_args <- function(args) {
  out <- list(
    `urd-rds` = NULL,
    outdir = NULL,
    `rgc1-col` = "jia_score_RGC1",
    `rgc2-col` = "jia_score_RGC2",
    `ipc-col` = "jia_score_IPC",
    `pseudotime-name` = "",
    `proliferation-genes` = "MKI67,TOP2A,CENPF,UBE2C,PCNA,MCM2,MCM3,MCM4,MCM5,MCM6,MCM7,STMN1,HMGB2,AURKB,CDK1,CCNB1,CCNB2",
    `marker-genes` = "HES1,FGFR2,NOTCH1,NOTCH2,VIM,NES,DACH1,DLX1,DLX2,ASCL1",
    `candidate-pcts` = "1,2,5,10",
    `selected-pct` = "2",
    `pool-col` = "",
    `pool-value` = "",
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
    "  Rscript scripts/19_div30_jia_rootscore_candidates.R --urd-rds <object.rds> --outdir <dir>",
    "",
    "Outputs:",
    "  div30_urd_jia_rootscore_object.rds",
    "  tables/root_score_all_cells.tsv",
    "  tables/root_score_root_pool_cells.tsv",
    "  tables/root_score_distribution.tsv",
    "  tables/root_score_candidate_counts.tsv",
    "  tables/root_score_program_marker_summary.tsv",
    "  tables/root_score_candidate_cells_top{pct}.tsv",
    "  tables/root_score_marker_expression_by_candidate_set.tsv",
    "  tables/root_score_pseudotime_by_candidate_set.tsv",
    "  plots/root_score_distribution.png",
    "  plots/root_score_umap.png",
    "  plots/root_score_candidate_umap.png",
    "  plots/root_score_by_pseudotime.png",
    "  jia_rootscore_candidate_report.md",
    sep = "\n"
  ))
}

required <- c("Matrix", "URD", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Missing required R package(s): ", paste(missing, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages({
  library(Matrix)
  library(URD)
  library(ggplot2)
})

split_csv <- function(x) {
  vals <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  vals[nzchar(vals)]
}

as_num <- function(x, name) {
  value <- suppressWarnings(as.numeric(x))
  if (any(is.na(value))) stop(name, " must be numeric; got ", paste(x, collapse = ","), call. = FALSE)
  value
}

zscore <- function(x) {
  x <- as.numeric(x)
  mu <- mean(x, na.rm = TRUE)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mu) / s
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

find_genes_case_insensitive <- function(requested, available) {
  idx <- match(toupper(requested), toupper(available))
  found <- available[idx[!is.na(idx)]]
  names(found) <- requested[!is.na(idx)]
  found
}

extract_pseudotime <- function(object, requested) {
  pt <- object@pseudotime
  if (!is.data.frame(pt) && !is.matrix(pt)) return(list(name = "", values = rep(NA_real_, nrow(object@meta))))
  if (ncol(pt) == 0) return(list(name = "", values = rep(NA_real_, nrow(object@meta))))
  name <- requested
  if (!nzchar(name)) name <- colnames(pt)[[1]]
  if (!(name %in% colnames(pt))) {
    stop("Pseudotime column not found: ", name, ". Available: ", paste(colnames(pt), collapse = ", "), call. = FALSE)
  }
  values <- as.numeric(pt[rownames(object@meta), name])
  names(values) <- rownames(object@meta)
  list(name = name, values = values)
}

expression_matrix <- function(object, genes, cells) {
  object@logupx.data[genes, cells, drop = FALSE]
}

compute_scores <- function(object, cfg) {
  meta <- object@meta[rownames(object@meta), , drop = FALSE]
  needed <- c(cfg$rgc1_col, cfg$rgc2_col, cfg$ipc_col)
  missing_cols <- setdiff(needed, colnames(meta))
  if (length(missing_cols) > 0) stop("Missing metadata column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)

  requested_prolif <- split_csv(cfg$proliferation_genes)
  requested_markers <- split_csv(cfg$marker_genes)
  available_genes <- rownames(object@logupx.data)
  proliferation_genes <- find_genes_case_insensitive(requested_prolif, available_genes)
  marker_genes <- find_genes_case_insensitive(requested_markers, available_genes)
  if (length(proliferation_genes) == 0) stop("No proliferation genes were found in URD logupx.data.", call. = FALSE)

  cells <- rownames(meta)
  proliferation_score <- Matrix::colMeans(expression_matrix(object, unname(proliferation_genes), cells))
  names(proliferation_score) <- cells

  pt <- extract_pseudotime(object, cfg$pseudotime_name)
  score_df <- data.frame(
    cell_id = cells,
    jia_score_RGC1 = as.numeric(meta[[cfg$rgc1_col]]),
    jia_score_RGC2 = as.numeric(meta[[cfg$rgc2_col]]),
    jia_score_IPC = as.numeric(meta[[cfg$ipc_col]]),
    proliferation_score = as.numeric(proliferation_score[cells]),
    urd_pseudotime = as.numeric(pt$values[cells]),
    stringsAsFactors = FALSE
  )
  if (all(c("UMAP_1", "UMAP_2") %in% colnames(meta))) {
    score_df$UMAP_1 <- as.numeric(meta$UMAP_1)
    score_df$UMAP_2 <- as.numeric(meta$UMAP_2)
  }
  if (nzchar(cfg$pool_col)) {
    if (!(cfg$pool_col %in% colnames(meta))) {
      stop("Pool column not found in URD metadata: ", cfg$pool_col, call. = FALSE)
    }
    score_df[[cfg$pool_col]] <- meta[[cfg$pool_col]]
  }
  summary_genes <- c("HES1", "VIM", "NES", "DLX1", "DLX2", "ASCL1")
  found_summary_genes <- find_genes_case_insensitive(summary_genes, available_genes)
  for (gene in summary_genes) {
    col <- paste0("logupx_", gene)
    if (gene %in% names(found_summary_genes)) {
      score_df[[col]] <- as.numeric(expression_matrix(object, unname(found_summary_genes[[gene]]), cells))
    } else {
      score_df[[col]] <- NA_real_
    }
  }
  score_df$z_RGC1 <- zscore(score_df$jia_score_RGC1)
  score_df$z_IPC <- zscore(score_df$jia_score_IPC)
  score_df$z_proliferation <- zscore(score_df$proliferation_score)
  score_df$RootScore <- score_df$z_RGC1 + score_df$z_proliferation - score_df$z_IPC
  score_df <- score_df[order(score_df$RootScore, decreasing = TRUE), , drop = FALSE]
  score_df$RootScore_rank <- seq_len(nrow(score_df))
  score_df$RootScore_percentile_top <- 100 * score_df$RootScore_rank / nrow(score_df)

  list(
    score_df = score_df,
    proliferation_genes = proliferation_genes,
    marker_genes = marker_genes,
    pseudotime_name = pt$name
  )
}

select_root_pool <- function(score_df, cfg) {
  if (!nzchar(cfg$pool_col)) {
    out <- score_df
    out$RootScore_pool_rank <- seq_len(nrow(out))
    out$RootScore_pool_percentile_top <- 100 * out$RootScore_pool_rank / nrow(out)
    return(out)
  }
  if (!nzchar(cfg$pool_value)) {
    stop("--pool-value is required when --pool-col is set", call. = FALSE)
  }
  pool <- score_df[as.character(score_df[[cfg$pool_col]]) == cfg$pool_value, , drop = FALSE]
  if (nrow(pool) == 0) {
    stop("Root pool selected zero cells: ", cfg$pool_col, " == ", cfg$pool_value, call. = FALSE)
  }
  pool <- pool[order(pool$RootScore, decreasing = TRUE), , drop = FALSE]
  pool$RootScore_pool_rank <- seq_len(nrow(pool))
  pool$RootScore_pool_percentile_top <- 100 * pool$RootScore_pool_rank / nrow(pool)
  pool
}

distribution_table <- function(score_df) {
  probs <- c(0, 0.01, 0.02, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.98, 0.99, 1)
  q <- stats::quantile(score_df$RootScore, probs = probs, na.rm = TRUE)
  data.frame(
    metric = c("n_cells", paste0("quantile_", names(q)), "mean", "sd", "min", "max"),
    value = c(
      nrow(score_df),
      as.numeric(q),
      mean(score_df$RootScore, na.rm = TRUE),
      stats::sd(score_df$RootScore, na.rm = TRUE),
      min(score_df$RootScore, na.rm = TRUE),
      max(score_df$RootScore, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
}

candidate_counts <- function(score_df, candidate_pcts, selected_pct) {
  rows <- lapply(candidate_pcts, function(pct) {
    n <- ceiling(nrow(score_df) * pct / 100)
    cells <- head(score_df, n)
    data.frame(
      candidate_set = paste0("top", pct, "pct"),
      top_percent = pct,
      n_cells = nrow(cells),
      rootscore_min = min(cells$RootScore, na.rm = TRUE),
      rootscore_median = median(cells$RootScore, na.rm = TRUE),
      rootscore_max = max(cells$RootScore, na.rm = TRUE),
      median_RGC1 = median(cells$jia_score_RGC1, na.rm = TRUE),
      median_IPC = median(cells$jia_score_IPC, na.rm = TRUE),
      median_proliferation = median(cells$proliferation_score, na.rm = TRUE),
      median_urd_pseudotime = median(cells$urd_pseudotime, na.rm = TRUE),
      selected_for_alternative_root = pct == selected_pct,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

root_program_marker_summary <- function(score_df, root_pool_df, selected_pct, cfg) {
  selected_n <- ceiling(nrow(root_pool_df) * selected_pct / 100)
  selected_df <- head(root_pool_df, selected_n)
  metric_cols <- c(
    "jia_score_RGC1",
    "jia_score_RGC2",
    "jia_score_IPC",
    "logupx_HES1",
    "logupx_VIM",
    "logupx_NES",
    "logupx_DLX1",
    "logupx_DLX2",
    "logupx_ASCL1",
    "RootScore",
    "proliferation_score"
  )
  summarize_one <- function(df, group, pool_col = cfg$pool_col, pool_value = cfg$pool_value) {
    values <- vapply(metric_cols, function(col) mean(as.numeric(df[[col]]), na.rm = TRUE), numeric(1))
    values[!is.finite(values)] <- NA_real_
    out <- data.frame(
      comparison_group = group,
      root_pool_col = pool_col,
      root_pool_value = pool_value,
      selected_top_percent = selected_pct,
      n_cells = nrow(df),
      percent_of_all_scored_cells = 100 * nrow(df) / nrow(score_df),
      percent_of_root_pool = 100 * nrow(df) / nrow(root_pool_df),
      stringsAsFactors = FALSE
    )
    for (col in metric_cols) out[[paste0("mean_", col)]] <- unname(values[[col]])
    out
  }

  rows <- list(
    selected_roots = summarize_one(selected_df, paste0("selected_roots_top", selected_pct, "pct")),
    root_pool = summarize_one(root_pool_df, "root_pool_all"),
    all_scored = summarize_one(score_df, "all_scored_cells", "", "")
  )
  if (nzchar(cfg$pool_col) && cfg$pool_col %in% colnames(score_df)) {
    pool_values <- sort(unique(as.character(score_df[[cfg$pool_col]])))
    for (value in pool_values) {
      df <- score_df[as.character(score_df[[cfg$pool_col]]) == value, , drop = FALSE]
      key <- paste0("pool_level_", make.names(value))
      rows[[key]] <- summarize_one(df, paste0(cfg$pool_col, "=", value), cfg$pool_col, value)
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

marker_summary <- function(object, score_df, marker_genes, candidate_pcts) {
  if (length(marker_genes) == 0) {
    return(data.frame(status = "no_requested_marker_genes_found", stringsAsFactors = FALSE))
  }
  rows <- list()
  all_cells <- score_df$cell_id
  for (pct in candidate_pcts) {
    n <- ceiling(nrow(score_df) * pct / 100)
    candidate_cells <- head(score_df$cell_id, n)
    other_cells <- setdiff(all_cells, candidate_cells)
    mat_candidate <- expression_matrix(object, unname(marker_genes), candidate_cells)
    mat_other <- expression_matrix(object, unname(marker_genes), other_cells)
    rows[[paste0("top", pct)]] <- data.frame(
      candidate_set = paste0("top", pct, "pct"),
      top_percent = pct,
      marker = names(marker_genes),
      gene = unname(marker_genes),
      mean_logupx_candidate = as.numeric(Matrix::rowMeans(mat_candidate)),
      mean_logupx_other = as.numeric(Matrix::rowMeans(mat_other)),
      median_logupx_candidate = apply(mat_candidate, 1, median),
      pct_expressed_candidate = as.numeric(Matrix::rowSums(mat_candidate > 0) / ncol(mat_candidate)),
      pct_expressed_other = as.numeric(Matrix::rowSums(mat_other > 0) / ncol(mat_other)),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

pseudotime_summary <- function(score_df, candidate_pcts) {
  rows <- lapply(candidate_pcts, function(pct) {
    n <- ceiling(nrow(score_df) * pct / 100)
    x <- head(score_df$urd_pseudotime, n)
    data.frame(
      candidate_set = paste0("top", pct, "pct"),
      top_percent = pct,
      n_cells = n,
      median_urd_pseudotime = median(x, na.rm = TRUE),
      mean_urd_pseudotime = mean(x, na.rm = TRUE),
      q25_urd_pseudotime = as.numeric(stats::quantile(x, 0.25, na.rm = TRUE)),
      q75_urd_pseudotime = as.numeric(stats::quantile(x, 0.75, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

add_candidate_columns <- function(object, score_df, root_pool_df, candidate_pcts, selected_pct) {
  object@meta$jia_rootscore <- NA_real_
  object@meta$jia_rootscore_z_RGC1 <- NA_real_
  object@meta$jia_rootscore_z_IPC <- NA_real_
  object@meta$jia_rootscore_z_proliferation <- NA_real_
  object@meta$jia_rootscore_proliferation_score <- NA_real_
  object@meta$jia_rootscore_root_pool <- FALSE
  object@meta$jia_rootscore_selected_root <- FALSE
  object@meta$jia_rootscore_selected_top_percent <- NA_real_
  object@meta[score_df$cell_id, "jia_rootscore"] <- score_df$RootScore
  object@meta[score_df$cell_id, "jia_rootscore_z_RGC1"] <- score_df$z_RGC1
  object@meta[score_df$cell_id, "jia_rootscore_z_IPC"] <- score_df$z_IPC
  object@meta[score_df$cell_id, "jia_rootscore_z_proliferation"] <- score_df$z_proliferation
  object@meta[score_df$cell_id, "jia_rootscore_proliferation_score"] <- score_df$proliferation_score
  object@meta[root_pool_df$cell_id, "jia_rootscore_root_pool"] <- TRUE
  for (pct in candidate_pcts) {
    col <- paste0("jia_rootscore_top", gsub("\\.", "p", as.character(pct)), "pct_candidate")
    object@meta[[col]] <- FALSE
    cells <- head(root_pool_df$cell_id, ceiling(nrow(root_pool_df) * pct / 100))
    object@meta[cells, col] <- TRUE
  }
  selected_cells <- head(root_pool_df$cell_id, ceiling(nrow(root_pool_df) * selected_pct / 100))
  object@meta[selected_cells, "jia_rootscore_selected_root"] <- TRUE
  object@meta[selected_cells, "jia_rootscore_selected_top_percent"] <- selected_pct
  object
}

plot_distribution <- function(score_df, path) {
  p <- ggplot(score_df, aes(RootScore)) +
    geom_histogram(bins = 80, fill = "#4d9221", color = "white", linewidth = 0.1) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank()) +
    labs(title = "Jia proliferative VZ-RGC RootScore distribution", x = "RootScore", y = "Cells")
  ggsave(path, p, width = 7, height = 5, dpi = 240, bg = "white")
}

plot_umap <- function(score_df, candidate_pcts, path_score, path_candidates) {
  if (!all(c("UMAP_1", "UMAP_2") %in% colnames(score_df))) return(invisible(FALSE))
  rank_col <- if ("RootScore_pool_rank" %in% colnames(score_df)) "RootScore_pool_rank" else "RootScore_rank"
  p1 <- ggplot(score_df, aes(UMAP_1, UMAP_2, color = RootScore)) +
    geom_point(size = 0.35, alpha = 0.85) +
    coord_equal() +
    scale_color_viridis_c() +
    theme_void(base_size = 11) +
    theme(plot.background = element_rect(fill = "white", color = NA)) +
    labs(title = "RootScore on UMAP", color = "RootScore")
  ggsave(path_score, p1, width = 7, height = 6, dpi = 240, bg = "white")

  rows <- lapply(candidate_pcts, function(pct) {
    df <- score_df
    df$candidate_set <- paste0("top", pct, "pct")
    df$is_candidate <- df[[rank_col]] <= ceiling(nrow(df) * pct / 100)
    df
  })
  cand <- do.call(rbind, rows)
  p2 <- ggplot(cand, aes(UMAP_1, UMAP_2)) +
    geom_point(color = "grey82", size = 0.25, alpha = 0.55) +
    geom_point(data = cand[cand$is_candidate, , drop = FALSE], color = "#d73027", size = 0.42, alpha = 0.9) +
    coord_equal() +
    facet_wrap(~candidate_set) +
    theme_void(base_size = 10) +
    theme(plot.background = element_rect(fill = "white", color = NA), strip.text = element_text(face = "bold")) +
    labs(title = "Top RootScore candidate roots on UMAP")
  ggsave(path_candidates, p2, width = 10, height = 8, dpi = 240, bg = "white")
  invisible(TRUE)
}

plot_pseudotime <- function(score_df, candidate_pcts, path) {
  rank_col <- if ("RootScore_pool_rank" %in% colnames(score_df)) "RootScore_pool_rank" else "RootScore_rank"
  rows <- lapply(candidate_pcts, function(pct) {
    df <- score_df
    df$candidate_set <- paste0("top", pct, "pct")
    df$is_candidate <- df[[rank_col]] <= ceiling(nrow(df) * pct / 100)
    df
  })
  df <- do.call(rbind, rows)
  p <- ggplot(df, aes(urd_pseudotime, RootScore)) +
    geom_point(color = "grey80", size = 0.25, alpha = 0.4) +
    geom_point(data = df[df$is_candidate, , drop = FALSE], color = "#d73027", size = 0.45, alpha = 0.9) +
    facet_wrap(~candidate_set) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank()) +
    labs(title = "Top RootScore candidates over URD pseudotime", x = "Existing URD pseudotime", y = "RootScore")
  ggsave(path, p, width = 10, height = 8, dpi = 240, bg = "white")
}

write_report <- function(path, cfg, genes, counts, dist, root_pool_df, selected_pct) {
  pool_text <- if (nzchar(cfg$pool_col)) {
    paste0("`", cfg$pool_col, " == ", cfg$pool_value, "`")
  } else {
    "all scored cells"
  }
  lines <- c(
    "# Jia Proliferative VZ-RGC RootScore Candidate Report",
    "",
    paste0("- Input URD object: `", cfg$urd_rds, "`"),
    paste0("- Root pool: ", pool_text),
    paste0("- Root pool cells: ", nrow(root_pool_df)),
    paste0("- Selected alternative root set: top ", selected_pct, "% by RootScore within the root pool"),
    "- Scoring uses the formula below; candidate percentage is applied after optional root-pool filtering.",
    "",
    "## Formula",
    "",
    "`RootScore = z(RGC1_score) + z(proliferation_score) - z(IPC_score)`",
    "",
    "## Genes Used",
    "",
    paste0("- Proliferation genes found: `", paste(unname(genes$proliferation_genes), collapse = ", "), "`"),
    paste0("- Marker genes found: `", paste(unname(genes$marker_genes), collapse = ", "), "`"),
    "",
    "## Candidate Counts",
    "",
    paste(capture.output(print(counts, row.names = FALSE)), collapse = "\n"),
    "",
    "## RootScore Distribution",
    "",
    paste(capture.output(print(dist, row.names = FALSE)), collapse = "\n")
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
  table_dir = file.path(opt$outdir, "tables"),
  plot_dir = file.path(opt$outdir, "plots"),
  rgc1_col = opt$`rgc1-col`,
  rgc2_col = opt$`rgc2-col`,
  ipc_col = opt$`ipc-col`,
  pseudotime_name = opt$`pseudotime-name`,
  proliferation_genes = opt$`proliferation-genes`,
  marker_genes = opt$`marker-genes`,
  candidate_pcts = as_num(split_csv(opt$`candidate-pcts`), "candidate-pcts"),
  selected_pct = as_num(opt$`selected-pct`, "selected-pct"),
  pool_col = opt$`pool-col`,
  pool_value = opt$`pool-value`
)
dir.create(cfg$table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$plot_dir, recursive = TRUE, showWarnings = FALSE)

log_msg("Reading URD object: ", cfg$urd_rds)
urd <- readRDS(cfg$urd_rds)
if (!inherits(urd, "URD")) stop("Input is not an URD object.", call. = FALSE)

computed <- compute_scores(urd, cfg)
score_df <- computed$score_df
root_pool_df <- select_root_pool(score_df, cfg)
dist <- distribution_table(score_df)
counts <- candidate_counts(root_pool_df, cfg$candidate_pcts, cfg$selected_pct)
program_marker_summary <- root_program_marker_summary(score_df, root_pool_df, cfg$selected_pct, cfg)
markers <- marker_summary(urd, root_pool_df, computed$marker_genes, cfg$candidate_pcts)
pt_summary <- pseudotime_summary(root_pool_df, cfg$candidate_pcts)
urd <- add_candidate_columns(urd, score_df, root_pool_df, cfg$candidate_pcts, cfg$selected_pct)

write_tsv(score_df, file.path(cfg$table_dir, "root_score_all_cells.tsv"))
write_tsv(root_pool_df, file.path(cfg$table_dir, "root_score_root_pool_cells.tsv"))
write_tsv(dist, file.path(cfg$table_dir, "root_score_distribution.tsv"))
write_tsv(counts, file.path(cfg$table_dir, "root_score_candidate_counts.tsv"))
write_tsv(program_marker_summary, file.path(cfg$table_dir, "root_score_program_marker_summary.tsv"))
write_tsv(markers, file.path(cfg$table_dir, "root_score_marker_expression_by_candidate_set.tsv"))
write_tsv(pt_summary, file.path(cfg$table_dir, "root_score_pseudotime_by_candidate_set.tsv"))
write_tsv(
  data.frame(requested = names(computed$proliferation_genes), gene = unname(computed$proliferation_genes), role = "proliferation", stringsAsFactors = FALSE),
  file.path(cfg$table_dir, "root_score_proliferation_genes_used.tsv")
)
write_tsv(
  data.frame(requested = names(computed$marker_genes), gene = unname(computed$marker_genes), role = "marker", stringsAsFactors = FALSE),
  file.path(cfg$table_dir, "root_score_marker_genes_used.tsv")
)
for (pct in cfg$candidate_pcts) {
  n <- ceiling(nrow(root_pool_df) * pct / 100)
  write_tsv(head(root_pool_df, n), file.path(cfg$table_dir, paste0("root_score_candidate_cells_top", pct, "pct.tsv")))
}

plot_distribution(score_df, file.path(cfg$plot_dir, "root_score_distribution.png"))
plot_umap(root_pool_df, cfg$candidate_pcts, file.path(cfg$plot_dir, "root_score_umap.png"), file.path(cfg$plot_dir, "root_score_candidate_umap.png"))
plot_pseudotime(root_pool_df, cfg$candidate_pcts, file.path(cfg$plot_dir, "root_score_by_pseudotime.png"))

out_rds <- file.path(cfg$outdir, "div30_urd_jia_rootscore_object.rds")
saveRDS(urd, out_rds)
write_report(
  file.path(cfg$outdir, "jia_rootscore_candidate_report.md"),
  cfg,
  computed,
  counts,
  dist,
  root_pool_df,
  cfg$selected_pct
)
log_msg("Saved scored URD object: ", out_rds)
log_msg("Wrote RootScore report: ", file.path(cfg$outdir, "jia_rootscore_candidate_report.md"))
