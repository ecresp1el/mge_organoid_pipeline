#!/usr/bin/env Rscript

# QC reconstruction sweep for Wang et al. 2025 hnbMO scRNA-seq (GSE286235).
# GEO provides raw_feature_bc_matrix files only, so the "Cell Ranger-style"
# stage here is an explicit cell-like droplet proxy, not the authors' exact
# Cell Ranger filtered_feature_bc_matrix output.

suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
})

parse_args <- function(args) {
  out <- list(
    `project-root` = Sys.getenv("PROJECT_ROOT", unset = ""),
    `tenx-dir` = "",
    outdir = "",
    `cell-call-min-features` = "200,500,1000,1200,1360,1500,1800",
    `mt-thresholds` = "5,10,15,20",
    `target-cells` = "14245",
    `plot-sample-size` = "200000",
    seed = "11",
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
    if (a == "-p") a <- "--project-root"
    if (a == "-o") a <- "--outdir"
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
  cat(
    paste(
      "Usage:",
      "  Rscript scripts/05l_liu_2025_hnbmo_qc_threshold_sweep.R --project-root <PROJECT_ROOT> [options]",
      "",
      "Options:",
      "  --tenx-dir <path>                 Converted raw 10x directories by sample_id",
      "  --outdir <path>                   Output directory",
      "  --cell-call-min-features <csv>    Candidate cell-like droplet nFeature gates",
      "  --mt-thresholds <csv>             Mitochondrial percentage thresholds",
      "  --target-cells <int>              Reported healthy-cell target (default: 14245)",
      "  --plot-sample-size <int>          Max raw/proxy points per sample for distribution plots",
      "  --seed <int>                      Random seed for downsampled plotting",
      sep = "\n"
    )
  )
}

parse_number_csv <- function(x) {
  as.numeric(strsplit(x, ",", fixed = TRUE)[[1L]])
}

write_tsv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(df, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
}

read_tsv_gz <- function(path, header = FALSE) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  read.delim(con, header = header, stringsAsFactors = FALSE)
}

sample_info <- data.frame(
  sample_id = c("BF_H9_D36", "BF_H9_D63", "BFCO_IMR_D63"),
  figure1c_sample = c("H9 D36", "H9 D63", "IMR90-4 D63"),
  cell_line = c("H9", "H9", "IMR90-4"),
  day = c("D36", "D63", "D63"),
  stringsAsFactors = FALSE
)

read_sample_metrics <- function(tenx_dir, sample_id, sample_label) {
  sample_dir <- file.path(tenx_dir, sample_id)
  matrix_path <- file.path(sample_dir, "matrix.mtx.gz")
  features_path <- file.path(sample_dir, "features.tsv.gz")
  barcodes_path <- file.path(sample_dir, "barcodes.tsv.gz")

  if (!file.exists(matrix_path)) stop("Missing matrix: ", matrix_path)
  if (!file.exists(features_path)) stop("Missing features: ", features_path)
  if (!file.exists(barcodes_path)) stop("Missing barcodes: ", barcodes_path)

  message("Reading ", sample_id, " from ", sample_dir)
  mat <- readMM(matrix_path)
  features <- read_tsv_gz(features_path, header = FALSE)
  barcodes <- read_tsv_gz(barcodes_path, header = FALSE)

  if (nrow(features) != nrow(mat)) stop("Feature count mismatch for ", sample_id)
  if (nrow(barcodes) != ncol(mat)) stop("Barcode count mismatch for ", sample_id)

  gene_names <- make.unique(as.character(features[[2L]]))
  mt_idx <- grepl("^MT-", gene_names)

  n_count <- as.numeric(Matrix::colSums(mat))
  n_feature <- as.integer(Matrix::colSums(mat > 0))
  mt_count <- if (any(mt_idx)) as.numeric(Matrix::colSums(mat[mt_idx, , drop = FALSE])) else rep(0, ncol(mat))
  percent_mt <- ifelse(n_count > 0, 100 * mt_count / n_count, NA_real_)

  data.frame(
    sample_id = sample_id,
    figure1c_sample = sample_label,
    barcode = as.character(barcodes[[1L]]),
    nCount_RNA = n_count,
    nFeature_RNA = n_feature,
    percent.mt = percent_mt,
    stringsAsFactors = FALSE
  )
}

summarize_metric <- function(df, stage, min_features, mt_threshold) {
  if (nrow(df) == 0L) {
    return(data.frame(
      stage = stage,
      cell_call_min_features = min_features,
      mt_threshold = mt_threshold,
      sample_id = NA_character_,
      figure1c_sample = NA_character_,
      n_barcodes = 0L,
      median_nFeature_RNA = NA_real_,
      median_nCount_RNA = NA_real_,
      q05_nFeature_RNA = NA_real_,
      q10_nFeature_RNA = NA_real_,
      q25_nFeature_RNA = NA_real_,
      q75_nFeature_RNA = NA_real_,
      q95_nFeature_RNA = NA_real_,
      q05_nCount_RNA = NA_real_,
      q10_nCount_RNA = NA_real_,
      q25_nCount_RNA = NA_real_,
      q75_nCount_RNA = NA_real_,
      q95_nCount_RNA = NA_real_
    ))
  }

  split_df <- split(df, df$figure1c_sample, drop = TRUE)
  rows <- lapply(names(split_df), function(label) {
    x <- split_df[[label]]
    data.frame(
      stage = stage,
      cell_call_min_features = min_features,
      mt_threshold = mt_threshold,
      sample_id = x$sample_id[[1L]],
      figure1c_sample = label,
      n_barcodes = nrow(x),
      median_nFeature_RNA = median(x$nFeature_RNA, na.rm = TRUE),
      median_nCount_RNA = median(x$nCount_RNA, na.rm = TRUE),
      q05_nFeature_RNA = as.numeric(quantile(x$nFeature_RNA, 0.05, na.rm = TRUE)),
      q10_nFeature_RNA = as.numeric(quantile(x$nFeature_RNA, 0.10, na.rm = TRUE)),
      q25_nFeature_RNA = as.numeric(quantile(x$nFeature_RNA, 0.25, na.rm = TRUE)),
      q75_nFeature_RNA = as.numeric(quantile(x$nFeature_RNA, 0.75, na.rm = TRUE)),
      q95_nFeature_RNA = as.numeric(quantile(x$nFeature_RNA, 0.95, na.rm = TRUE)),
      q05_nCount_RNA = as.numeric(quantile(x$nCount_RNA, 0.05, na.rm = TRUE)),
      q10_nCount_RNA = as.numeric(quantile(x$nCount_RNA, 0.10, na.rm = TRUE)),
      q25_nCount_RNA = as.numeric(quantile(x$nCount_RNA, 0.25, na.rm = TRUE)),
      q75_nCount_RNA = as.numeric(quantile(x$nCount_RNA, 0.75, na.rm = TRUE)),
      q95_nCount_RNA = as.numeric(quantile(x$nCount_RNA, 0.95, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

downsample_by_sample <- function(df, n_per_sample, seed) {
  set.seed(seed)
  split_df <- split(df, df$figure1c_sample, drop = TRUE)
  rows <- lapply(split_df, function(x) {
    if (nrow(x) <= n_per_sample) return(x)
    x[sample.int(nrow(x), n_per_sample), , drop = FALSE]
  })
  do.call(rbind, rows)
}

make_long_metrics <- function(df) {
  rbind(
    data.frame(
      figure1c_sample = df$figure1c_sample,
      metric = "nFeature_RNA",
      value = df$nFeature_RNA,
      stringsAsFactors = FALSE
    ),
    data.frame(
      figure1c_sample = df$figure1c_sample,
      metric = "nCount_RNA",
      value = df$nCount_RNA,
      stringsAsFactors = FALSE
    )
  )
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(opt$help)) {
  print_usage()
  quit(save = "no", status = 0)
}

if (!nzchar(opt$`project-root`)) stop("PROJECT_ROOT or --project-root is required")
project_root <- sub("/+$", "", opt$`project-root`)
tenx_dir <- if (nzchar(opt$`tenx-dir`)) opt$`tenx-dir` else file.path(project_root, "data/raw/liu_2025_hnbmo_geo_files/suppl/10x")
outdir <- if (nzchar(opt$outdir)) opt$outdir else file.path(project_root, "results/liu_2025_hnbmo_qc_sweep")
plot_dir <- file.path(outdir, "plots")
table_dir <- file.path(outdir, "tables")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

cell_call_min_features <- as.integer(parse_number_csv(opt$`cell-call-min-features`))
mt_thresholds <- parse_number_csv(opt$`mt-thresholds`)
target_cells <- as.integer(opt$`target-cells`)
plot_sample_size <- as.integer(opt$`plot-sample-size`)
seed <- as.integer(opt$seed)

message("Project root: ", project_root)
message("10x dir: ", tenx_dir)
message("Output dir: ", outdir)
message("Cell-like nFeature gates: ", paste(cell_call_min_features, collapse = ", "))
message("percent.mt thresholds: ", paste(mt_thresholds, collapse = ", "))

metrics_list <- lapply(seq_len(nrow(sample_info)), function(i) {
  read_sample_metrics(tenx_dir, sample_info$sample_id[[i]], sample_info$figure1c_sample[[i]])
})
metrics <- do.call(rbind, metrics_list)
metrics$figure1c_sample <- factor(metrics$figure1c_sample, levels = sample_info$figure1c_sample)

raw_summary <- summarize_metric(metrics, "raw_observed_barcodes", NA_integer_, NA_real_)
write_tsv(raw_summary, file.path(table_dir, "raw_observed_barcode_qc_summary.tsv"))

rank_rows <- lapply(split(metrics, metrics$figure1c_sample, drop = TRUE), function(x) {
  counts <- sort(x$nCount_RNA, decreasing = TRUE)
  data.frame(
    figure1c_sample = x$figure1c_sample[[1L]],
    rank = seq_along(counts),
    nCount_RNA = counts,
    stringsAsFactors = FALSE
  )
})
rank_df <- do.call(rbind, rank_rows)

p_rank <- ggplot(rank_df, aes(x = rank, y = nCount_RNA, color = figure1c_sample)) +
  geom_line(linewidth = 0.35) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Barcode rank", y = "nCount_RNA", color = NULL, title = "Raw barcode-rank curves") +
  theme_classic(base_size = 13)
ggsave(file.path(plot_dir, "raw_barcode_rank_nCount_loglog.png"), p_rank, width = 7.5, height = 5, dpi = 300)
ggsave(file.path(plot_dir, "raw_barcode_rank_nCount_loglog.pdf"), p_rank, width = 7.5, height = 5)

raw_plot_df <- downsample_by_sample(metrics, plot_sample_size, seed)
raw_long <- make_long_metrics(raw_plot_df)
p_raw <- ggplot(raw_long, aes(x = figure1c_sample, y = value, fill = figure1c_sample)) +
  geom_violin(scale = "width", trim = TRUE, color = "gray25", linewidth = 0.2) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  scale_y_log10() +
  labs(x = NULL, y = NULL, title = "Raw observed barcode distributions before QC") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(plot_dir, "raw_observed_barcodes_nFeature_nCount_violin_log10.png"), p_raw, width = 9, height = 5, dpi = 300)
ggsave(file.path(plot_dir, "raw_observed_barcodes_nFeature_nCount_violin_log10.pdf"), p_raw, width = 9, height = 5)

summary_rows <- list()
selected_metric_rows <- list()
best_violin_inputs <- list()

for (min_feature in cell_call_min_features) {
  cell_like <- metrics[metrics$nFeature_RNA >= min_feature, , drop = FALSE]
  summary_rows[[paste0("cell_like_", min_feature)]] <- summarize_metric(
    cell_like,
    "cellranger_style_proxy_before_mt",
    min_feature,
    NA_real_
  )

  cell_like_plot <- downsample_by_sample(cell_like, plot_sample_size, seed)
  cell_like_long <- make_long_metrics(cell_like_plot)
  p_cell_like <- ggplot(cell_like_long, aes(x = figure1c_sample, y = value, fill = figure1c_sample)) +
    geom_violin(scale = "width", trim = TRUE, color = "gray25", linewidth = 0.2) +
    geom_boxplot(width = 0.12, outlier.shape = NA, color = "gray20", fill = "white", alpha = 0.75) +
    facet_wrap(~metric, scales = "free_y", nrow = 1) +
    labs(x = NULL, y = NULL, title = paste0("Cell Ranger-style proxy before percent.mt: nFeature_RNA >= ", min_feature)) +
    theme_classic(base_size = 13) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))
  ggsave(
    file.path(plot_dir, paste0("cellranger_style_proxy_before_mt_minFeature_", min_feature, "_violins.png")),
    p_cell_like,
    width = 9,
    height = 5,
    dpi = 300
  )

  mt_long_rows <- list()
  for (mt_threshold in mt_thresholds) {
    post_qc <- cell_like[!is.na(cell_like$percent.mt) & cell_like$percent.mt < mt_threshold, , drop = FALSE]
    key <- paste0("min", min_feature, "_mt", mt_threshold)
    summary_rows[[key]] <- summarize_metric(post_qc, "post_percent_mt_qc", min_feature, mt_threshold)

    sampled <- downsample_by_sample(post_qc, plot_sample_size, seed)
    sampled_long <- make_long_metrics(sampled)
    sampled_long$mt_threshold <- paste0("percent.mt < ", mt_threshold)
    sampled_long$cell_call_gate <- paste0("nFeature >= ", min_feature)
    mt_long_rows[[as.character(mt_threshold)]] <- sampled_long

    if (min_feature %in% c(1000L, 1200L, 1360L, 1500L)) {
      selected_metric_rows[[key]] <- sampled_long
    }
    if (min_feature == 1360L && mt_threshold == 20) {
      best_violin_inputs[[key]] <- sampled_long
    }
  }

  mt_long <- do.call(rbind, mt_long_rows)
  mt_long$mt_threshold <- factor(mt_long$mt_threshold, levels = paste0("percent.mt < ", mt_thresholds))
  p_mt <- ggplot(mt_long, aes(x = figure1c_sample, y = value, fill = figure1c_sample)) +
    geom_violin(scale = "width", trim = TRUE, color = "gray25", linewidth = 0.2) +
    geom_boxplot(width = 0.12, outlier.shape = NA, color = "gray20", fill = "white", alpha = 0.75) +
    facet_grid(metric ~ mt_threshold, scales = "free_y") +
    labs(x = NULL, y = NULL, title = paste0("Post-QC violins: nFeature_RNA >= ", min_feature)) +
    theme_classic(base_size = 12) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))
  ggsave(
    file.path(plot_dir, paste0("post_qc_mt_sweep_minFeature_", min_feature, "_violins.png")),
    p_mt,
    width = 13,
    height = 7,
    dpi = 300
  )
  ggsave(
    file.path(plot_dir, paste0("post_qc_mt_sweep_minFeature_", min_feature, "_violins.pdf")),
    p_mt,
    width = 13,
    height = 7
  )
}

all_summary <- do.call(rbind, summary_rows)
write_tsv(all_summary, file.path(table_dir, "qc_sweep_summary_by_sample.tsv"))

post_summary <- all_summary[all_summary$stage == "post_percent_mt_qc", , drop = FALSE]
total_summary <- aggregate(
  n_barcodes ~ cell_call_min_features + mt_threshold,
  data = post_summary,
  FUN = sum
)
names(total_summary)[names(total_summary) == "n_barcodes"] <- "total_cells"

median_summary <- aggregate(
  cbind(median_nFeature_RNA, median_nCount_RNA, q10_nFeature_RNA, q10_nCount_RNA) ~ cell_call_min_features + mt_threshold,
  data = post_summary,
  FUN = median
)

ranking <- merge(total_summary, median_summary, by = c("cell_call_min_features", "mt_threshold"), all = TRUE)
ranking$reported_target_cells <- target_cells
ranking$abs_delta_from_reported_14245 <- abs(ranking$total_cells - target_cells)
ranking <- ranking[order(ranking$abs_delta_from_reported_14245, ranking$cell_call_min_features, ranking$mt_threshold), ]
write_tsv(ranking, file.path(table_dir, "qc_sweep_ranking_against_reported_cell_count.tsv"))

selected_long <- do.call(rbind, selected_metric_rows)
if (!is.null(selected_long) && nrow(selected_long) > 0L) {
  selected_long$cell_call_gate <- factor(
    selected_long$cell_call_gate,
    levels = paste0("nFeature >= ", c(1000L, 1200L, 1360L, 1500L))
  )
  selected_long$mt_threshold <- factor(selected_long$mt_threshold, levels = paste0("percent.mt < ", mt_thresholds))
  p_selected <- ggplot(selected_long, aes(x = figure1c_sample, y = value, fill = figure1c_sample)) +
    geom_violin(scale = "width", trim = TRUE, color = "gray25", linewidth = 0.18) +
    geom_boxplot(width = 0.11, outlier.shape = NA, color = "gray20", fill = "white", alpha = 0.75) +
    facet_grid(metric + cell_call_gate ~ mt_threshold, scales = "free_y") +
    labs(x = NULL, y = NULL, title = "Selected post-QC violin panels for Figure S1M visual comparison") +
    theme_classic(base_size = 12) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))
  ggsave(file.path(plot_dir, "selected_post_qc_violin_panels_for_figS1M_comparison.png"), p_selected, width = 13, height = 12, dpi = 300)
  ggsave(file.path(plot_dir, "selected_post_qc_violin_panels_for_figS1M_comparison.pdf"), p_selected, width = 13, height = 12)
}

best_min_feature <- ranking$cell_call_min_features[[1L]]
best_mt <- ranking$mt_threshold[[1L]]
best_df <- metrics[
  metrics$nFeature_RNA >= best_min_feature &
    !is.na(metrics$percent.mt) &
    metrics$percent.mt < best_mt,
  ,
  drop = FALSE
]
best_long <- make_long_metrics(best_df)
p_best <- ggplot(best_long, aes(x = figure1c_sample, y = value, fill = figure1c_sample)) +
  geom_violin(scale = "width", trim = TRUE, color = "gray25", linewidth = 0.25) +
  geom_boxplot(width = 0.12, outlier.shape = NA, color = "gray20", fill = "white", alpha = 0.75) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  labs(
    x = NULL,
    y = NULL,
    title = paste0(
      "Best count match: nFeature_RNA >= ", best_min_feature,
      ", percent.mt < ", best_mt,
      " (n = ", ranking$total_cells[[1L]], ")"
    )
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(plot_dir, "best_count_match_nCount_nFeature_violin.png"), p_best, width = 9, height = 5, dpi = 300)
ggsave(file.path(plot_dir, "best_count_match_nCount_nFeature_violin.pdf"), p_best, width = 9, height = 5)

assumptions <- data.frame(
  key = c(
    "study",
    "reported_healthy_cells",
    "raw_data_limitation",
    "cellranger_style_proxy",
    "mt_thresholds_tested",
    "ranking_rule"
  ),
  value = c(
    "Wang et al. 2025 GSE286235 hnbMO healthy samples: H9 D36, H9 D63, IMR90-4 D63",
    as.character(target_cells),
    "GEO provides raw_feature_bc_matrix files, not filtered_feature_bc_matrix files or author QC code.",
    "Candidate nFeature_RNA gates on raw observed barcodes; this is a reconstruction proxy, not exact Cell Ranger v3.1 cell calling.",
    paste(mt_thresholds, collapse = ","),
    "Primary automatic ranking uses absolute difference from the reported 14,245 healthy-cell count; median/lower-tail quantities are reported for visual comparison to Figure S1M."
  ),
  stringsAsFactors = FALSE
)
write_tsv(assumptions, file.path(table_dir, "qc_sweep_assumptions.tsv"))

message("Best count-match combination: nFeature_RNA >= ", best_min_feature,
        ", percent.mt < ", best_mt,
        "; retained cells = ", ranking$total_cells[[1L]])
message("Wrote tables to: ", table_dir)
message("Wrote plots to: ", plot_dir)
