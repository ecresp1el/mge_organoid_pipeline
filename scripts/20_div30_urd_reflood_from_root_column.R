#!/usr/bin/env Rscript

# Re-run URD flood pseudotime from an existing logical root-candidate metadata
# column. This is for score-defined roots that are selected after the diffusion
# geometry already exists, such as `jia_rootscore_selected_root`.
#
# It does not recompute PCA or the diffusion map.

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = " ")))
  flush.console()
}

parse_args <- function(args) {
  out <- list(
    `urd-rds` = NULL,
    outdir = NULL,
    `root-col` = NULL,
    `pseudotime-name` = NULL,
    seed = "7",
    `n-floods` = "20",
    `flood-minimum-cells` = "2",
    `flood-max-frac-na` = "0.4",
    `flood-stability-div` = "10",
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

required <- c("URD")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Missing required R package(s): ", paste(missing, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(library(URD))

as_int <- function(x, name) {
  value <- suppressWarnings(as.integer(x))
  if (is.na(value)) stop(name, " must be an integer; got ", x, call. = FALSE)
  value
}

as_num <- function(x, name) {
  value <- suppressWarnings(as.numeric(x))
  if (is.na(value)) stop(name, " must be numeric; got ", x, call. = FALSE)
  value
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(opt$help)) {
  cat("Usage: Rscript scripts/20_div30_urd_reflood_from_root_column.R --urd-rds <object.rds> --outdir <dir> --root-col <logical_col> --pseudotime-name <name>\n")
  quit(save = "no", status = 0)
}
if (is.null(opt$`urd-rds`) || !nzchar(opt$`urd-rds`)) stop("--urd-rds is required", call. = FALSE)
if (is.null(opt$outdir) || !nzchar(opt$outdir)) stop("--outdir is required", call. = FALSE)
if (is.null(opt$`root-col`) || !nzchar(opt$`root-col`)) stop("--root-col is required", call. = FALSE)
if (is.null(opt$`pseudotime-name`) || !nzchar(opt$`pseudotime-name`)) stop("--pseudotime-name is required", call. = FALSE)

cfg <- list(
  urd_rds = opt$`urd-rds`,
  outdir = opt$outdir,
  table_dir = file.path(opt$outdir, "tables"),
  root_col = opt$`root-col`,
  pseudotime_name = opt$`pseudotime-name`,
  seed = as_int(opt$seed, "seed"),
  n_floods = as_int(opt$`n-floods`, "n-floods"),
  flood_minimum_cells = as_int(opt$`flood-minimum-cells`, "flood-minimum-cells"),
  flood_max_frac_na = as_num(opt$`flood-max-frac-na`, "flood-max-frac-na"),
  flood_stability_div = as_int(opt$`flood-stability-div`, "flood-stability-div")
)
dir.create(cfg$table_dir, recursive = TRUE, showWarnings = FALSE)

log_msg("Reading URD object: ", cfg$urd_rds)
urd <- readRDS(cfg$urd_rds)
if (!(cfg$root_col %in% colnames(urd@meta))) stop("Missing root column: ", cfg$root_col, call. = FALSE)
root_flag <- as.logical(urd@meta[[cfg$root_col]])
root_flag[is.na(root_flag)] <- FALSE
root_cells <- rownames(urd@meta)[root_flag]
if (length(root_cells) == 0) stop("Root column selected zero cells.", call. = FALSE)

set.seed(cfg$seed)
log_msg("Running floodPseudotime from ", length(root_cells), " root cells")
floods <- floodPseudotime(
  urd,
  root.cells = root_cells,
  n = cfg$n_floods,
  minimum.cells.flooded = cfg$flood_minimum_cells,
  verbose = TRUE
)
urd <- floodPseudotimeProcess(
  urd,
  floods,
  floods.name = cfg$pseudotime_name,
  max.frac.NA = cfg$flood_max_frac_na,
  pseudotime.fun = mean,
  stability.div = cfg$flood_stability_div
)

summary <- data.frame(
  key = c("input_rds", "root_col", "pseudotime_name", "n_root_cells", "n_floods", "flood_minimum_cells", "flood_max_frac_na", "flood_stability_div"),
  value = c(cfg$urd_rds, cfg$root_col, cfg$pseudotime_name, length(root_cells), cfg$n_floods, cfg$flood_minimum_cells, cfg$flood_max_frac_na, cfg$flood_stability_div),
  stringsAsFactors = FALSE
)
write_tsv(summary, file.path(cfg$table_dir, "reflood_root_summary.tsv"))
write_tsv(data.frame(cell_id = root_cells, stringsAsFactors = FALSE), file.path(cfg$table_dir, "reflood_root_cells.tsv"))

out_rds <- file.path(cfg$outdir, "div30_urd_reflood_object.rds")
saveRDS(urd, out_rds)
log_msg("Saved reflooded URD object: ", out_rds)
