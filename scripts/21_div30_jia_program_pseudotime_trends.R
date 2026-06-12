#!/usr/bin/env Rscript

# Plot Jia RGC1/RGC2/IPC program trends over an existing URD pseudotime.
#
# Programmatic flow:
#   1. Read a URD object with Jia program score metadata.
#   2. Extract the requested pseudotime vector.
#   3. Bin cells by pseudotime.
#   4. Compute mean RGC1, RGC2, and IPC score per bin.
#   5. Plot individual score-vs-pseudotime panels and a combined binned trend.
#   6. Report the pseudotime bin where each program reaches its maximum.

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = " ")))
  flush.console()
}

parse_args <- function(args) {
  out <- list(
    `urd-rds` = NULL,
    outdir = NULL,
    `pseudotime-name` = "",
    `rgc1-col` = "jia_score_RGC1",
    `rgc2-col` = "jia_score_RGC2",
    `ipc-col` = "jia_score_IPC",
    `n-bins` = "30",
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

required <- c("ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Missing required R package(s): ", paste(missing, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(library(ggplot2))

as_int <- function(x, name) {
  value <- suppressWarnings(as.integer(x))
  if (is.na(value)) stop(name, " must be an integer; got ", x, call. = FALSE)
  value
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

extract_pseudotime <- function(object, requested) {
  pt <- object@pseudotime
  if (!is.data.frame(pt) && !is.matrix(pt)) stop("URD object does not have a tabular pseudotime slot.", call. = FALSE)
  if (ncol(pt) == 0) stop("URD object has no pseudotime columns.", call. = FALSE)
  name <- requested
  if (!nzchar(name)) name <- colnames(pt)[[1]]
  if (!(name %in% colnames(pt))) {
    stop("Pseudotime column not found: ", name, ". Available: ", paste(colnames(pt), collapse = ", "), call. = FALSE)
  }
  values <- as.numeric(pt[rownames(object@meta), name])
  names(values) <- rownames(object@meta)
  list(name = name, values = values)
}

make_long <- function(df) {
  rbind(
    data.frame(cell_id = df$cell_id, pseudotime = df$pseudotime, program = "RGC1", score = df$RGC1, stringsAsFactors = FALSE),
    data.frame(cell_id = df$cell_id, pseudotime = df$pseudotime, program = "RGC2", score = df$RGC2, stringsAsFactors = FALSE),
    data.frame(cell_id = df$cell_id, pseudotime = df$pseudotime, program = "IPC", score = df$IPC, stringsAsFactors = FALSE)
  )
}

bin_programs <- function(df, n_bins) {
  rng <- range(df$pseudotime, na.rm = TRUE)
  breaks <- seq(rng[[1]], rng[[2]], length.out = n_bins + 1L)
  breaks <- unique(breaks)
  if (length(breaks) < 3) stop("Pseudotime range is too small to bin.", call. = FALSE)
  df$bin <- cut(df$pseudotime, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  rows <- lapply(sort(unique(df$bin)), function(bin_id) {
    x <- df[df$bin == bin_id, , drop = FALSE]
    data.frame(
      bin = bin_id,
      n_cells = nrow(x),
      pseudotime_min = min(x$pseudotime, na.rm = TRUE),
      pseudotime_mid = median(x$pseudotime, na.rm = TRUE),
      pseudotime_max = max(x$pseudotime, na.rm = TRUE),
      mean_RGC1 = mean(x$RGC1, na.rm = TRUE),
      mean_RGC2 = mean(x$RGC2, na.rm = TRUE),
      mean_IPC = mean(x$IPC, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

bin_long <- function(binned) {
  rbind(
    data.frame(bin = binned$bin, n_cells = binned$n_cells, pseudotime_mid = binned$pseudotime_mid, program = "RGC1", mean_score = binned$mean_RGC1, stringsAsFactors = FALSE),
    data.frame(bin = binned$bin, n_cells = binned$n_cells, pseudotime_mid = binned$pseudotime_mid, program = "RGC2", mean_score = binned$mean_RGC2, stringsAsFactors = FALSE),
    data.frame(bin = binned$bin, n_cells = binned$n_cells, pseudotime_mid = binned$pseudotime_mid, program = "IPC", mean_score = binned$mean_IPC, stringsAsFactors = FALSE)
  )
}

peak_table <- function(binned) {
  rows <- lapply(
    c(RGC1 = "mean_RGC1", RGC2 = "mean_RGC2", IPC = "mean_IPC"),
    function(col) {
      idx <- which.max(binned[[col]])
      data.frame(
        program = sub("mean_", "", col),
        peak_bin = binned$bin[[idx]],
        peak_pseudotime_mid = binned$pseudotime_mid[[idx]],
        peak_pseudotime_min = binned$pseudotime_min[[idx]],
        peak_pseudotime_max = binned$pseudotime_max[[idx]],
        peak_mean_score = binned[[col]][[idx]],
        n_cells_in_peak_bin = binned$n_cells[[idx]],
        stringsAsFactors = FALSE
      )
    }
  )
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

interpret_order <- function(peaks, binned) {
  p <- setNames(peaks$peak_pseudotime_mid, peaks$program)
  bin_width <- median(diff(sort(unique(binned$pseudotime_mid))), na.rm = TRUE)
  if (!is.finite(bin_width)) bin_width <- 0
  rg_overlap <- abs(p[["RGC1"]] - p[["RGC2"]]) <= bin_width
  if (p[["RGC1"]] < p[["RGC2"]] && p[["RGC2"]] < p[["IPC"]] && !rg_overlap) {
    call <- "RGC1 -> RGC2 -> IPC"
  } else if (rg_overlap && max(p[["RGC1"]], p[["RGC2"]]) < p[["IPC"]]) {
    call <- "RGC1/RGC2 overlap -> IPC"
  } else {
    call <- "alternative structure"
  }
  data.frame(
    interpretation = call,
    RGC1_peak = p[["RGC1"]],
    RGC2_peak = p[["RGC2"]],
    IPC_peak = p[["IPC"]],
    median_bin_spacing = bin_width,
    note = "Peaks are based on binned mean scores over all cells ordered by URD pseudotime.",
    stringsAsFactors = FALSE
  )
}

plot_individual <- function(long_df, path) {
  p <- ggplot(long_df, aes(pseudotime, score)) +
    geom_point(size = 0.25, alpha = 0.25, color = "grey35") +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, span = 0.5, color = "#1b7837", linewidth = 0.9) +
    facet_wrap(~program, scales = "free_y", ncol = 1) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank()) +
    labs(title = "Jia program scores vs URD pseudotime", x = "URD pseudotime", y = "Jia program score")
  ggsave(path, p, width = 8, height = 10, dpi = 240, bg = "white")
}

plot_combined <- function(binned_long, peaks, path) {
  p <- ggplot(binned_long, aes(pseudotime_mid, mean_score, color = program)) +
    geom_line(linewidth = 0.9) +
    geom_point(aes(size = n_cells), alpha = 0.8) +
    geom_vline(data = peaks, aes(xintercept = peak_pseudotime_mid, color = program), linetype = "dashed", linewidth = 0.45, show.legend = FALSE) +
    scale_size_continuous(range = c(1.2, 4)) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank()) +
    labs(title = "Binned Jia program trends over URD pseudotime", x = "URD pseudotime bin median", y = "Mean Jia score", color = "Program", size = "Cells/bin")
  ggsave(path, p, width = 8, height = 5.5, dpi = 240, bg = "white")
}

write_report <- function(path, cfg, peaks, interp) {
  lines <- c(
    "# Jia Program Pseudotime Trends",
    "",
    paste0("- URD object: `", cfg$urd_rds, "`"),
    paste0("- Pseudotime: `", cfg$pseudotime_name, "`"),
    paste0("- Bins: ", cfg$n_bins),
    "",
    "## Peak Pseudotimes",
    "",
    paste(capture.output(print(peaks, row.names = FALSE)), collapse = "\n"),
    "",
    "## Interpretation",
    "",
    paste(capture.output(print(interp, row.names = FALSE)), collapse = "\n")
  )
  writeLines(lines, path)
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(opt$help)) {
  cat("Usage: Rscript scripts/21_div30_jia_program_pseudotime_trends.R --urd-rds <object.rds> --outdir <dir> --pseudotime-name <name>\n")
  quit(save = "no", status = 0)
}
if (is.null(opt$`urd-rds`) || !nzchar(opt$`urd-rds`)) stop("--urd-rds is required", call. = FALSE)
if (is.null(opt$outdir) || !nzchar(opt$outdir)) stop("--outdir is required", call. = FALSE)

cfg <- list(
  urd_rds = opt$`urd-rds`,
  outdir = opt$outdir,
  table_dir = file.path(opt$outdir, "tables"),
  plot_dir = file.path(opt$outdir, "plots"),
  pseudotime_name = opt$`pseudotime-name`,
  rgc1_col = opt$`rgc1-col`,
  rgc2_col = opt$`rgc2-col`,
  ipc_col = opt$`ipc-col`,
  n_bins = as_int(opt$`n-bins`, "n-bins")
)
dir.create(cfg$table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$plot_dir, recursive = TRUE, showWarnings = FALSE)

log_msg("Reading URD object: ", cfg$urd_rds)
urd <- readRDS(cfg$urd_rds)
needed <- c(cfg$rgc1_col, cfg$rgc2_col, cfg$ipc_col)
missing_cols <- setdiff(needed, colnames(urd@meta))
if (length(missing_cols) > 0) stop("Missing metadata column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
pt <- extract_pseudotime(urd, cfg$pseudotime_name)
cfg$pseudotime_name <- pt$name

df <- data.frame(
  cell_id = rownames(urd@meta),
  pseudotime = as.numeric(pt$values[rownames(urd@meta)]),
  RGC1 = as.numeric(urd@meta[[cfg$rgc1_col]]),
  RGC2 = as.numeric(urd@meta[[cfg$rgc2_col]]),
  IPC = as.numeric(urd@meta[[cfg$ipc_col]]),
  stringsAsFactors = FALSE
)
df <- df[is.finite(df$pseudotime), , drop = FALSE]
long_df <- make_long(df)
binned <- bin_programs(df, cfg$n_bins)
binned_long <- bin_long(binned)
peaks <- peak_table(binned)
interp <- interpret_order(peaks, binned)

write_tsv(df, file.path(cfg$table_dir, "jia_program_scores_by_cell.tsv"))
write_tsv(binned, file.path(cfg$table_dir, "jia_program_binned_means.tsv"))
write_tsv(peaks, file.path(cfg$table_dir, "jia_program_peak_pseudotimes.tsv"))
write_tsv(interp, file.path(cfg$table_dir, "jia_program_ordering_interpretation.tsv"))
plot_individual(long_df, file.path(cfg$plot_dir, "jia_program_scores_vs_pseudotime.png"))
plot_combined(binned_long, peaks, file.path(cfg$plot_dir, "jia_program_binned_trends.png"))
write_report(file.path(cfg$outdir, "jia_program_pseudotime_trends_report.md"), cfg, peaks, interp)

log_msg("Wrote Jia program trend report: ", file.path(cfg$outdir, "jia_program_pseudotime_trends_report.md"))
