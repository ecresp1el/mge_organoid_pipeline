#!/usr/bin/env Rscript
#' Run approved Step 03 scDblFinder detection in native R.
#'
#' The script treats every retained cell as part of one GEX_1 capture, uses
#' clusters=TRUE and dbr.sd=1, leaves all model parameters at defaults, retains
#' the internal PCA only for diagnostics, repeats classification under a second
#' seed, and never removes a cell.

suppressPackageStartupMessages({
  library(data.table)
  library(DropletUtils)
  library(R6)
  library(scDblFinder)
  library(SingleCellExperiment)
})

#' Parse named command-line options supplied as --name value pairs.
#'
#' @return A named character list containing every provided option.
parse_options <- function(arguments) {
  if (length(arguments) %% 2L != 0L || any(!startsWith(arguments[seq(1L, length(arguments), 2L)], "--"))) {
    stop("Options must be supplied as --name value pairs")
  }
  keys <- sub("^--", "", arguments[seq(1L, length(arguments), 2L)])
  values <- arguments[seq(2L, length(arguments), 2L)]
  stats::setNames(as.list(values), keys)
}

#' Require a named command-line option.
#'
#' @return The non-empty option value.
require_option <- function(options, name) {
  value <- options[[name]]
  if (is.null(value) || !nzchar(value)) stop("Missing required option: --", name)
  value
}

#' Load and validate the one-capture native-R input.
Step03InputLoader <- R6Class(
  "Step03InputLoader",
  public = list(
    bridge_h5 = NULL,
    metadata_tsv = NULL,
    expected_cells = NULL,
    expected_genes = NULL,
    capture_id = NULL,

    #' Construct a loader for one exact bridge and ordered metadata table.
    initialize = function(bridge_h5, metadata_tsv, expected_cells, expected_genes, capture_id) {
      self$bridge_h5 <- bridge_h5
      self$metadata_tsv <- metadata_tsv
      self$expected_cells <- as.integer(expected_cells)
      self$expected_genes <- as.integer(expected_genes)
      self$capture_id <- capture_id
    },

    #' Read counts and attach the approved capture and reporting metadata.
    load = function() {
      sce <- DropletUtils::read10xCounts(self$bridge_h5, type = "HDF5", col.names = TRUE)
      metadata <- data.table::fread(self$metadata_tsv, sep = "\t", data.table = FALSE)
      stopifnot(identical(dim(sce), c(self$expected_genes, self$expected_cells)))
      stopifnot(nrow(metadata) == self$expected_cells)
      stopifnot(!anyDuplicated(metadata$cell_id))
      stopifnot(identical(as.character(colnames(sce)), as.character(metadata$cell_id)))
      stopifnot(identical(unique(as.character(metadata$capture_id)), self$capture_id))
      stopifnot(all(metadata$capture_id == self$capture_id))
      for (field in colnames(metadata)) {
        if (field != "cell_id") colData(sce)[[field]] <- metadata[[field]]
      }
      sce
    }
  )
)

#' Execute the approved primary and reproducibility scDblFinder passes.
Step03ScDblFinderRunner <- R6Class(
  "Step03ScDblFinderRunner",
  public = list(
    primary_seed = NULL,
    reproducibility_seed = NULL,

    #' Store the two predeclared deterministic random seeds.
    initialize = function(primary_seed, reproducibility_seed) {
      self$primary_seed <- as.integer(primary_seed)
      self$reproducibility_seed <- as.integer(reproducibility_seed)
    },

    #' Run one model with only the approved scientific arguments changed.
    run_once = function(sce, seed, return_type) {
      set.seed(seed)
      scDblFinder::scDblFinder(
        sce,
        samples = "capture_id",
        clusters = TRUE,
        dbr.sd = 1,
        returnType = return_type
      )
    },

    #' Run the primary full diagnostic return and an independent-seed repeat.
    run = function(sce) {
      primary <- self$run_once(sce, self$primary_seed, "full")
      real <- which(as.character(primary$src) == "real")
      stopifnot(length(real) == ncol(sce))
      stopifnot(identical(as.character(colnames(primary)[real]), as.character(colnames(sce))))
      primary_table <- as.data.frame(colData(primary)[real, , drop = FALSE])
      primary_table$cell_id <- colnames(primary)[real]
      primary_pca <- reducedDim(primary, "PCA")[real, , drop = FALSE]
      rownames(primary_pca) <- primary_table$cell_id
      primary_stats <- metadata(primary)$scDblFinder.stats
      rm(primary)
      invisible(gc())

      replicate <- self$run_once(sce, self$reproducibility_seed, "scores")
      replicate_table <- as.data.frame(replicate)
      replicate_table$cell_id <- rownames(replicate)
      stopifnot(identical(as.character(replicate_table$cell_id), as.character(colnames(sce))))
      list(
        primary_table = primary_table,
        primary_pca = primary_pca,
        primary_stats = primary_stats,
        replicate_table = replicate_table
      )
    }
  )
)

#' Persist complete scDblFinder results and exact package/runtime provenance.
Step03ResultWriter <- R6Class(
  "Step03ResultWriter",
  public = list(
    output_dir = NULL,

    #' Construct a writer that refuses to replace an existing R output folder.
    initialize = function(output_dir) {
      self$output_dir <- output_dir
      if (dir.exists(output_dir) || file.exists(output_dir)) stop("Refusing existing R output directory: ", output_dir)
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    },

    #' Normalize one computed field name for a stable cross-language schema.
    normalize_name = function(name) {
      name <- sub("^scDblFinder\\.", "", name)
      gsub("[^A-Za-z0-9]+", "_", name)
    },

    #' Write aligned primary/replicate fields, PCA, model stats, and versions.
    write = function(result) {
      primary <- result$primary_table
      replicate <- result$replicate_table
      primary_names <- vapply(colnames(primary), self$normalize_name, character(1L))
      colnames(primary) <- paste0("primary_detail_", primary_names)
      colnames(primary)[colnames(primary) == "primary_detail_cell_id"] <- "cell_id"
      replicate <- replicate[, c("cell_id", "score", "class"), drop = FALSE]
      colnames(replicate) <- c("cell_id", "replicate_score", "replicate_class")
      output <- merge(primary, replicate, by = "cell_id", sort = FALSE)
      output <- output[match(rownames(result$primary_pca), output$cell_id), , drop = FALSE]
      required <- c("primary_detail_score", "primary_detail_class", "primary_detail_cluster")
      stopifnot(all(required %in% colnames(output)))
      setnames(output, required, c("primary_score", "primary_class", "primary_cluster"))
      data.table::fwrite(output, file.path(self$output_dir, "scdblfinder_per_cell_results.tsv.gz"), sep = "\t", compress = "gzip")
      self$write_pca(result$primary_pca)
      saveRDS(result$primary_stats, file.path(self$output_dir, "scdblfinder_primary_stats.rds"), compress = "xz")
      self$write_versions()
      self$write_method_contract()
    },

    #' Write the real-cell internal PCA as aligned portable and native objects.
    write_pca = function(pca) {
      pca_table <- data.frame(cell_id = rownames(pca), as.matrix(pca), check.names = FALSE)
      data.table::fwrite(pca_table, file.path(self$output_dir, "scdblfinder_internal_pca.tsv.gz"), sep = "\t", compress = "gzip")
      saveRDS(pca, file.path(self$output_dir, "scdblfinder_internal_pca.rds"), compress = "gzip")
    },

    #' Write exact R and package versions used by the native detector.
    write_versions = function() {
      packages <- c("scDblFinder", "SingleCellExperiment", "DropletUtils", "data.table", "R6", "rhdf5", "BiocParallel")
      versions <- data.frame(
        component = c("R", packages),
        version = c(R.version.string, vapply(packages, function(x) as.character(utils::packageVersion(x)), character(1L)))
      )
      data.table::fwrite(versions, file.path(self$output_dir, "r_software_versions.tsv"), sep = "\t")
      capture.output(sessionInfo(), file = file.path(self$output_dir, "r_session_info.txt"))
    },

    #' Write the exact scientific and output-only arguments used by both runs.
    write_method_contract = function() {
      contract <- data.frame(
        setting = c("samples", "clusters", "dbr.sd", "dbr", "other_model_parameters", "primary_returnType", "replicate_returnType", "cell_removal"),
        value = c("capture_id (constant GEX_1)", "TRUE", "1", "not supplied", "package defaults", "full", "scores", "none"),
        role = c("approved scientific", "approved scientific", "approved scientific", "approved scientific", "approved scientific", "output-only: retain internal PCA and diagnostics", "output-only: retain repeat scores/calls", "review boundary")
      )
      data.table::fwrite(contract, file.path(self$output_dir, "scdblfinder_method_contract.tsv"), sep = "\t")
    }
  )
)

#' Coordinate native-R loading, classification, validation, and persistence.
Step03RWorkflow <- R6Class(
  "Step03RWorkflow",
  public = list(
    options = NULL,

    #' Store parsed command-line options for one frozen run.
    initialize = function(options) {
      self$options <- options
    },

    #' Execute two non-filtering scDblFinder passes and persist review evidence.
    run = function() {
      loader <- Step03InputLoader$new(
        require_option(self$options, "bridge-h5"),
        require_option(self$options, "metadata-tsv"),
        require_option(self$options, "expected-cells"),
        require_option(self$options, "expected-genes"),
        require_option(self$options, "capture-id")
      )
      sce <- loader$load()
      runner <- Step03ScDblFinderRunner$new(
        require_option(self$options, "primary-seed"),
        require_option(self$options, "reproducibility-seed")
      )
      result <- runner$run(sce)
      writer <- Step03ResultWriter$new(require_option(self$options, "output-dir"))
      writer$write(result)
      message("Native R Step 03 results written; no cells were removed.")
    }
  )
)

options <- parse_options(commandArgs(trailingOnly = TRUE))
Step03RWorkflow$new(options)$run()
