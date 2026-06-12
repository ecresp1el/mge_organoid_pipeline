#!/usr/bin/env Rscript

# Recompute DIV30 URD flood pseudotime using a Jia hierarchy-defined
# VZ-RGC-like alternative root.
#
# Programmatic root definition:
#   1. Start within cells annotated as Radial glia.
#   2. Keep high RGC1, low RGC2, and low IPC cells using configurable
#      quantile thresholds within Radial glia.
#   3. If proliferation marker genes are available, optionally restrict to the
#      most proliferative subset by a mean logUPX proliferation score.
#   4. Re-run only flood pseudotime from this root on the existing URD diffusion
#      geometry. PCA and diffusion map are not recomputed.

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = " ")))
  flush.console()
}

parse_args <- function(args) {
  out <- list(
    `urd-rds` = NULL,
    outdir = NULL,
    `annotation-col` = "paper_cluster_annotation",
    `radial-glia-label` = "Radial glia",
    `pseudotime-name` = "jia_vz_rgc_root",
    `rgc1-col` = "jia_score_RGC1",
    `rgc2-col` = "jia_score_RGC2",
    `ipc-col` = "jia_score_IPC",
    `rgc1-high-quantile` = "0.75",
    `rgc2-low-quantile` = "0.50",
    `ipc-low-quantile` = "0.50",
    `proliferation-genes` = "MKI67,TOP2A,PCNA,MCM2,MCM3,MCM4,MCM5,MCM6,MCM7,STMN1,HMGB2,UBE2C,AURKB,CDK1,CCNB1,CCNB2",
    `proliferation-high-quantile` = "0.50",
    `require-proliferation` = "false",
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

required <- c("Matrix", "URD", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Missing required R package(s): ", paste(missing, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages({
  library(Matrix)
  library(URD)
  library(ggplot2)
})

as_num <- function(x, name) {
  value <- suppressWarnings(as.numeric(x))
  if (is.na(value)) stop(name, " must be numeric; got ", x, call. = FALSE)
  value
}

as_int <- function(x, name) {
  value <- suppressWarnings(as.integer(x))
  if (is.na(value)) stop(name, " must be an integer; got ", x, call. = FALSE)
  value
}

as_bool <- function(x) tolower(x) %in% c("true", "t", "1", "yes", "y")

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

split_csv <- function(x) {
  vals <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  vals[nzchar(vals)]
}

find_genes_case_insensitive <- function(requested, available) {
  idx <- match(toupper(requested), toupper(available))
  available[idx[!is.na(idx)]]
}

select_jia_vz_root <- function(object, cfg) {
  meta <- object@meta
  needed <- c(cfg$annotation_col, cfg$rgc1_col, cfg$rgc2_col, cfg$ipc_col)
  missing_cols <- setdiff(needed, colnames(meta))
  if (length(missing_cols) > 0) stop("Missing metadata column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)

  rg_cells <- rownames(meta)[meta[[cfg$annotation_col]] == cfg$radial_glia_label]
  if (length(rg_cells) == 0) stop("No Radial glia cells found.", call. = FALSE)
  rg <- meta[rg_cells, , drop = FALSE]

  rgc1_threshold <- unname(quantile(rg[[cfg$rgc1_col]], cfg$rgc1_high_quantile, na.rm = TRUE))
  rgc2_threshold <- unname(quantile(rg[[cfg$rgc2_col]], cfg$rgc2_low_quantile, na.rm = TRUE))
  ipc_threshold <- unname(quantile(rg[[cfg$ipc_col]], cfg$ipc_low_quantile, na.rm = TRUE))

  jia_flag <- rg[[cfg$rgc1_col]] >= rgc1_threshold &
    rg[[cfg$rgc2_col]] <= rgc2_threshold &
    rg[[cfg$ipc_col]] <= ipc_threshold
  jia_flag[is.na(jia_flag)] <- FALSE

  requested_prolif <- split_csv(cfg$proliferation_genes)
  found_prolif <- find_genes_case_insensitive(requested_prolif, rownames(object@logupx.data))
  proliferation_score <- rep(NA_real_, nrow(meta))
  names(proliferation_score) <- rownames(meta)
  used_proliferation <- FALSE
  proliferation_threshold <- NA_real_
  final_flag <- rep(FALSE, nrow(rg))
  names(final_flag) <- rg_cells
  final_flag[jia_flag] <- TRUE

  if (length(found_prolif) > 0) {
    score <- Matrix::colMeans(object@logupx.data[found_prolif, rownames(meta), drop = FALSE])
    proliferation_score[names(score)] <- as.numeric(score)
    candidate_cells <- rg_cells[jia_flag]
    proliferation_threshold <- unname(quantile(proliferation_score[candidate_cells], cfg$proliferation_high_quantile, na.rm = TRUE))
    pro_flag <- proliferation_score[rg_cells] >= proliferation_threshold
    pro_flag[is.na(pro_flag)] <- FALSE
    final_flag <- jia_flag & pro_flag
    used_proliferation <- TRUE
  } else if (cfg$require_proliferation) {
    stop("No proliferation marker genes were found, but --require-proliferation is true.", call. = FALSE)
  }

  root_cells <- rg_cells[final_flag]
  if (length(root_cells) == 0) stop("Jia VZ root selection produced zero cells.", call. = FALSE)

  object@meta$jia_vz_rgc_root_candidate <- rownames(meta) %in% root_cells
  object@meta$jia_vz_rgc_jia_hierarchy_candidate <- FALSE
  object@meta[rg_cells, "jia_vz_rgc_jia_hierarchy_candidate"] <- jia_flag
  object@meta$jia_vz_rgc_proliferation_score <- proliferation_score[rownames(meta)]

  thresholds <- data.frame(
    criterion = c("Radial glia label", "RGC1 high", "RGC2 low", "IPC low", "Proliferation refinement"),
    column_or_genes = c(
      cfg$annotation_col,
      cfg$rgc1_col,
      cfg$rgc2_col,
      cfg$ipc_col,
      paste(found_prolif, collapse = ",")
    ),
    quantile = c(NA, cfg$rgc1_high_quantile, cfg$rgc2_low_quantile, cfg$ipc_low_quantile, if (used_proliferation) cfg$proliferation_high_quantile else NA),
    threshold = c(NA, rgc1_threshold, rgc2_threshold, ipc_threshold, proliferation_threshold),
    used = c(TRUE, TRUE, TRUE, TRUE, used_proliferation),
    note = c(
      cfg$radial_glia_label,
      "Keep cells at or above this threshold within Radial glia.",
      "Keep cells at or below this threshold within Radial glia.",
      "Keep cells at or below this threshold within Radial glia.",
      if (used_proliferation) "Proliferation markers found and used." else "No proliferation marker genes found; skipped."
    ),
    stringsAsFactors = FALSE
  )

  summary <- data.frame(
    key = c(
      "n_total_cells",
      "n_radial_glia_cells",
      "n_jia_hierarchy_candidates",
      "n_final_jia_vz_root_cells",
      "used_proliferation_refinement",
      "found_proliferation_genes"
    ),
    value = c(
      nrow(meta),
      length(rg_cells),
      sum(jia_flag),
      length(root_cells),
      used_proliferation,
      paste(found_prolif, collapse = ",")
    ),
    stringsAsFactors = FALSE
  )

  list(object = object, root_cells = root_cells, thresholds = thresholds, summary = summary)
}

run_flood <- function(object, root_cells, cfg) {
  set.seed(cfg$seed)
  floods <- floodPseudotime(
    object,
    root.cells = root_cells,
    n = cfg$n_floods,
    minimum.cells.flooded = cfg$flood_minimum_cells,
    verbose = TRUE
  )
  floodPseudotimeProcess(
    object,
    floods,
    floods.name = cfg$pseudotime_name,
    max.frac.NA = cfg$flood_max_frac_na,
    pseudotime.fun = mean,
    stability.div = cfg$flood_stability_div
  )
}

write_root_umap <- function(object, cfg, path) {
  if (!all(c("UMAP_1", "UMAP_2") %in% colnames(object@meta))) return(invisible(FALSE))
  df <- object@meta
  df$cell_id <- rownames(df)
  df$root_status <- ifelse(df$jia_vz_rgc_root_candidate, "Jia VZ-RGC root", as.character(df[[cfg$annotation_col]]))
  p <- ggplot(df, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = root_status == "Jia VZ-RGC root"), size = 0.3, alpha = 0.75) +
    coord_equal() +
    scale_color_manual(values = c("TRUE" = "#d73027", "FALSE" = "grey75"), labels = c("Other", "Jia VZ-RGC root")) +
    theme_void(base_size = 11) +
    theme(plot.background = element_rect(fill = "white", color = NA)) +
    labs(title = "Jia VZ-RGC-like alternative root", color = NULL)
  ggsave(path, p, width = 7, height = 6, dpi = 240, bg = "white")
  invisible(TRUE)
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(opt$help)) {
  cat("Usage: Rscript scripts/18_div30_urd_jia_vz_root_pseudotime.R --urd-rds <object.rds> --outdir <dir>\n")
  quit(save = "no", status = 0)
}
if (is.null(opt$`urd-rds`) || !nzchar(opt$`urd-rds`)) stop("--urd-rds is required", call. = FALSE)
if (is.null(opt$outdir) || !nzchar(opt$outdir)) stop("--outdir is required", call. = FALSE)

cfg <- list(
  urd_rds = opt$`urd-rds`,
  outdir = opt$outdir,
  table_dir = file.path(opt$outdir, "tables"),
  plot_dir = file.path(opt$outdir, "plots"),
  annotation_col = opt$`annotation-col`,
  radial_glia_label = opt$`radial-glia-label`,
  pseudotime_name = opt$`pseudotime-name`,
  rgc1_col = opt$`rgc1-col`,
  rgc2_col = opt$`rgc2-col`,
  ipc_col = opt$`ipc-col`,
  rgc1_high_quantile = as_num(opt$`rgc1-high-quantile`, "rgc1-high-quantile"),
  rgc2_low_quantile = as_num(opt$`rgc2-low-quantile`, "rgc2-low-quantile"),
  ipc_low_quantile = as_num(opt$`ipc-low-quantile`, "ipc-low-quantile"),
  proliferation_genes = opt$`proliferation-genes`,
  proliferation_high_quantile = as_num(opt$`proliferation-high-quantile`, "proliferation-high-quantile"),
  require_proliferation = as_bool(opt$`require-proliferation`),
  seed = as_int(opt$seed, "seed"),
  n_floods = as_int(opt$`n-floods`, "n-floods"),
  flood_minimum_cells = as_int(opt$`flood-minimum-cells`, "flood-minimum-cells"),
  flood_max_frac_na = as_num(opt$`flood-max-frac-na`, "flood-max-frac-na"),
  flood_stability_div = as_int(opt$`flood-stability-div`, "flood-stability-div")
)
dir.create(cfg$table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$plot_dir, recursive = TRUE, showWarnings = FALSE)

log_msg("Reading URD object: ", cfg$urd_rds)
urd <- readRDS(cfg$urd_rds)
selected <- select_jia_vz_root(urd, cfg)
urd <- selected$object

log_msg("Selected Jia VZ-RGC-like root cells: ", length(selected$root_cells))
log_msg("Recomputing flood pseudotime: ", cfg$pseudotime_name)
urd <- run_flood(urd, selected$root_cells, cfg)

out_rds <- file.path(cfg$outdir, "div30_urd_jia_vz_root_object.rds")
saveRDS(urd, out_rds)
write_tsv(selected$thresholds, file.path(cfg$table_dir, "jia_vz_root_thresholds.tsv"))
write_tsv(selected$summary, file.path(cfg$table_dir, "jia_vz_root_summary.tsv"))
write_tsv(data.frame(cell_id = selected$root_cells, stringsAsFactors = FALSE), file.path(cfg$table_dir, "jia_vz_root_cells.tsv"))
write_root_umap(urd, cfg, file.path(cfg$plot_dir, "jia_vz_root_umap.png"))

log_msg("Saved Jia VZ-root URD object: ", out_rds)
