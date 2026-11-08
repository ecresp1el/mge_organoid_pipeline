# URD pseudotime and lineage reconstruction functions.
# Counts are genes x cells; metadata row names must match count columns.

# Write analysis progress to the R message stream.
# Inputs: text fragments; fragments are concatenated without adding separators.
# Output: a message; no URD slots or analysis parameters are modified.
log_msg <- function(...) message(paste0(..., collapse = ""))

# Merge named parameter overrides into a default configuration.
# Inputs: default values and a named list of overrides.
# Returns: the merged list. Unknown names stop execution to catch misspelled settings.
check_parameters <- function(defaults, overrides) {
  unknown <- setdiff(names(overrides), names(defaults))
  if (length(unknown)) stop("Unknown parameters: ", paste(unknown, collapse = ", "))
  utils::modifyList(defaults, overrides)
}

# Create and normalize an URD object, then identify surviving root cells.
# Inputs: bundle$counts (genes x cells), bundle$metadata, and cfg filtering thresholds.
# Returns: a list containing the filtered URD object and retained root cell IDs.
# Root flags come from metadata; an empty root set after filtering stops execution.
create_filtered_urd <- function(bundle, cfg) {
  log_msg(
    "Creating URD object with min.genes=", cfg$min_genes,
    " min.cells=", cfg$min_cells,
    " min.counts=", cfg$min_counts
  )
  urd <- createURD(
    count.data = bundle$counts,
    meta = bundle$metadata,
    min.cells = cfg$min_cells,
    min.genes = cfg$min_genes,
    min.counts = cfg$min_counts,
    verbose = TRUE
  )

  # Derive roots from retained cells so filtered-out barcodes cannot seed floods.
  cells <- rownames(urd@meta)
  root_cells <- cells[as.logical(urd@meta[cells, "urd_root_candidate"])]
  if (length(root_cells) == 0) {
    stop("No root candidate cells remained after URD filtering", call. = FALSE)
  }
  log_msg("Root cells retained after filtering: ", length(root_cells))
  list(urd = urd, root_cells = root_cells)
}

# Select genes with the greatest variance in normalized log expression.
# Inputs: an URD object with sparse logupx.data and cfg$num_variable_genes.
# Returns: the object with ordered gene names in var.genes; expression is unchanged.
# Variance is E[X^2] - E[X]^2, without binning genes by mean expression.
select_variable_genes <- function(urd, cfg) {
  log_msg("Selecting top variable genes from URD logupx.data")
  log_data <- urd@logupx.data
  gene_means <- Matrix::rowMeans(log_data)
  # Squaring sparse nonzero entries preserves implicit zeros without densifying.
  log_data_sq <- log_data
  log_data_sq@x <- log_data_sq@x^2
  gene_vars <- Matrix::rowMeans(log_data_sq) - gene_means^2
  # Rank only finite variances, then retain at most the requested number of genes.
  gene_vars <- gene_vars[is.finite(gene_vars)]
  gene_vars <- sort(gene_vars, decreasing = TRUE)
  urd@var.genes <- names(gene_vars)[seq_len(min(cfg$num_variable_genes, length(gene_vars)))]
  log_msg("Variable genes stored in URD object: ", length(urd@var.genes))
  urd
}

# Calculate PCA followed by the configured graph construction method.
# Inputs: an URD object with variable genes and the geometry configuration.
# Returns: the object with PCA scores and graph information; sets cfg$seed.
run_urd_geometry <- function(urd, cfg) {
  set.seed(cfg$seed)
  log_msg("Running calcPCA with mp.factor=", cfg$pca_mp_factor)
  urd <- calcPCA(urd, mp.factor = cfg$pca_mp_factor)
  run_urd_diffusion_map(urd, cfg)
}

# Calculate URD principal components from the selected variable genes.
# Inputs: an URD object, cfg$seed, and the Marchenko-Pastur factor cfg$pca_mp_factor.
# Returns: the object with PCA results. The random seed is reset for this stage.
run_urd_pca <- function(urd, cfg) {
  set.seed(cfg$seed)
  log_msg("Running calcPCA with mp.factor=", cfg$pca_mp_factor)
  calcPCA(urd, mp.factor = cfg$pca_mp_factor)
}

# Return the largest index represented by R's standard integer type.
# Output: a numeric scalar used to check the n_cells^2 size of pairwise geometry.
matrix_index_ceiling <- function() {
  as.numeric(.Machine$integer.max)
}

# Resolve an explicit graph method or the internal automatic selection rule.
# Inputs: the number of cells in urd@meta and cfg$dm_method.
# Returns: urd_expression or pca_knn. In auto mode the pairwise index ceiling
# selects the method; the public urd_pseudotime() interface requires an explicit choice.
resolve_dm_method <- function(urd, cfg) {
  n_cells <- nrow(urd@meta)
  pairwise_entries <- as.numeric(n_cells) * as.numeric(n_cells)
  if (cfg$dm_method == "auto" && pairwise_entries > matrix_index_ceiling()) {
    return("pca_knn")
  }
  if (cfg$dm_method == "auto") return("urd_expression")
  cfg$dm_method
}

# Check the pairwise index limit before expression-based graph construction.
# Inputs: an URD object and its geometry configuration.
# Returns: invisible TRUE, or stops when n_cells^2 exceeds the integer index limit.
# This check addresses indexing capacity; it does not estimate available memory.
guard_urd_expression_dm <- function(urd, cfg) {
  n_cells <- nrow(urd@meta)
  pairwise_entries <- as.numeric(n_cells) * as.numeric(n_cells)
  if (pairwise_entries > matrix_index_ceiling()) {
    stop(
      "URD calcDM cannot process ", n_cells, " cells because n_cells^2=",
      format(pairwise_entries, scientific = FALSE),
      " exceeds R Matrix's 2^31-1 sparse-index ceiling. Use --dm-method pca_knn ",
      "or downsample below ~46k cells before running expression-based calcDM.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# Return cell-aligned PC coordinates for nearest-neighbor search.
# Inputs: populated urd@pca.scores and cfg$dm_pca_dims (NA selects all stored PCs).
# Returns: a dense cells x PCs matrix ordered like urd@meta, with at least two PCs.
extract_pca_matrix <- function(urd, cfg) {
  if (nrow(urd@pca.scores) == 0 || ncol(urd@pca.scores) == 0) {
    stop("pca_knn diffusion requires populated urd@pca.scores. Run calcPCA first.", call. = FALSE)
  }
  pca <- as.matrix(urd@pca.scores)
  if (is.null(rownames(pca))) rownames(pca) <- rownames(urd@meta)
  # Graph rows and columns must follow the metadata cell order exactly.
  pca <- pca[rownames(urd@meta), , drop = FALSE]
  n_dims <- if (is.na(cfg$dm_pca_dims)) ncol(pca) else min(cfg$dm_pca_dims, ncol(pca))
  if (n_dims < 2) stop("pca_knn diffusion requires at least two PC dimensions.", call. = FALSE)
  pca[, seq_len(n_dims), drop = FALSE]
}

# Build a sparse, symmetric Gaussian affinity graph in PC space.
# Inputs: a named cells x PCs matrix, neighbor count, and kernel-width settings.
# Returns: a cell-named sparse affinity matrix, before pseudotime weighting.
# Local widths use each cell's kth-neighbor distance; numeric or global widths
# use a shared bandwidth. Directed neighbor weights are averaged with their transpose.
pca_knn_transition_matrix <- function(coords, cfg) {
  if (!requireNamespace("RANN", quietly = TRUE)) {
    stop("pca_knn diffusion requires R package RANN.", call. = FALSE)
  }
  n_cells <- nrow(coords)
  k <- min(cfg$knn, n_cells - 1L)
  if (k < 2) stop("pca_knn diffusion requires at least three cells.", call. = FALSE)
  log_msg("Building PCA-kNN graph with cells=", n_cells, " dims=", ncol(coords), " k=", k)
  # Request one extra neighbor because the first match is the query cell itself.
  nn <- RANN::nn2(coords, k = k + 1L)
  idx <- nn$nn.idx[, -1L, drop = FALSE]
  dist <- nn$nn.dists[, -1L, drop = FALSE]

  sigma_text <- cfg$sigma_text
  if (identical(tolower(sigma_text), "local")) {
    # Local bandwidths adapt to density; each edge uses sigma_i * sigma_j.
    sigma <- dist[, k]
    positive <- sigma[is.finite(sigma) & sigma > 0]
    if (length(positive) == 0) stop("Could not estimate local sigma from PCA-kNN distances.", call. = FALSE)
    # Replace zero/invalid widths to keep the Gaussian denominator positive.
    sigma[!is.finite(sigma) | sigma <= 0] <- stats::median(positive)
    weights <- exp(-(dist^2) / pmax(sigma[row(dist)] * sigma[idx], .Machine$double.eps))
  } else if (is.numeric(cfg$sigma)) {
    weights <- exp(-(dist^2) / pmax(cfg$sigma^2, .Machine$double.eps))
  } else {
    positive <- dist[is.finite(dist) & dist > 0]
    global_sigma <- stats::median(positive)
    log_msg("Using median PCA-kNN distance as global sigma=", signif(global_sigma, 4))
    weights <- exp(-(dist^2) / pmax(global_sigma^2, .Machine$double.eps))
  }

  i <- rep(seq_len(n_cells), each = k)
  tm <- Matrix::sparseMatrix(
    i = i,
    j = as.vector(t(idx)),
    x = as.vector(t(weights)),
    dims = c(n_cells, n_cells)
  )
  # Average reciprocal edges; this is an affinity matrix, not a row-normalized
  # probability matrix. URD consumes it before pseudotime-directed weighting.
  tm <- Matrix::forceSymmetric((tm + Matrix::t(tm)) / 2, uplo = "U")
  rownames(tm) <- rownames(coords)
  colnames(tm) <- rownames(coords)
  Matrix::drop0(tm)
}

# Import a PC-space affinity graph into URD's graph representation.
# Inputs: an object with PCA scores and the neighbor/kernel configuration.
# Returns: an object with the affinity graph imported through a DiffusionMap container.
# The container stores scaled PCs as coordinates and NA eigenvalues: this path
# does not calculate diffusion eigenvectors or a diffusion eigendecomposition.
run_pca_knn_diffusion_map <- function(urd, cfg) {
  coords <- extract_pca_matrix(urd, cfg)
  tm <- pca_knn_transition_matrix(coords, cfg)
  # Store scaled PC coordinates in the container; no eigensolver runs here.
  eig <- scale(coords)
  eig[!is.finite(eig)] <- 0
  n_eigs <- min(200L, ncol(eig))
  eig <- eig[, seq_len(n_eigs), drop = FALSE]
  colnames(eig) <- paste0("DC", seq_len(ncol(eig)))
  rownames(eig) <- rownames(coords)

  # Supply the slots expected by the DiffusionMap class used by importDM().
  sigmas <- methods::new(
    "Sigmas",
    log_sigmas = NULL,
    dim_norms = NULL,
    optimal_sigma = NULL,
    optimal_idx = NULL,
    avrd_norms = NULL
  )
  dm <- methods::new(
    "DiffusionMap",
    eigenvalues = rep(NA_real_, ncol(eig)),
    eigenvectors = eig,
    sigmas = sigmas,
    data_env = new.env(parent = emptyenv()),
    eigenvec0 = rep(1, nrow(eig)),
    transitions = tm,
    d = as.numeric(Matrix::rowSums(tm)),
    d_norm = as.numeric(Matrix::rowSums(tm)),
    k = cfg$knn,
    n_local = min(cfg$knn, 7L),
    density_norm = FALSE,
    rotate = FALSE,
    distance = "euclidean",
    censor_val = NULL,
    censor_range = NULL,
    missing_range = NULL,
    vars = colnames(coords)
  )
  log_msg("Importing PCA-kNN transition matrix into URD object; nnz=", length(dm@transitions@x))
  importDM(urd, dm)
}

# Dispatch graph construction to the selected expression or PC-space method.
# Inputs: an URD object and cfg geometry settings, including the random seed.
# Returns: the updated object. Expression mode calls calcDM with knn and sigma.use;
# PC-space mode imports the affinity graph built by pca_knn_transition_matrix().
run_urd_diffusion_map <- function(urd, cfg) {
  set.seed(cfg$seed)
  method <- resolve_dm_method(urd, cfg)
  log_msg("Resolved diffusion map method: ", method, " (requested=", cfg$dm_method, ")")
  if (method == "pca_knn") {
    return(run_pca_knn_diffusion_map(urd, cfg))
  }
  guard_urd_expression_dm(urd, cfg)
  log_msg("Running URD calcDM with knn=", cfg$knn, " sigma=", cfg$sigma_text)
  calcDM(urd, knn = cfg$knn, sigma.use = cfg$sigma)
}

# Estimate root-oriented pseudotime by repeated stochastic flooding.
# Inputs: an URD graph, named root cells, flood count, and processing thresholds.
# Returns: the object with pseudotime stored under cfg$pseudotime_name.
# Flood replicates are combined with the mean; max.frac.NA controls missingness
# and stability.div controls the flood-processing stability evaluation.
run_flood_pseudotime <- function(urd, root_cells, cfg) {
  set.seed(cfg$seed)
  log_msg(
    "Running floodPseudotime n=", cfg$n_floods,
    " minimum.cells.flooded=", cfg$flood_minimum_cells
  )
  # Independent floods sample graph reachability from the same root set.
  floods <- floodPseudotime(
    urd,
    root.cells = root_cells,
    n = cfg$n_floods,
    minimum.cells.flooded = cfg$flood_minimum_cells,
    verbose = TRUE
  )
  log_msg(
    "Processing flood pseudotime max.frac.NA=", cfg$flood_max_frac_na,
    " stability.div=", cfg$flood_stability_div
  )
  # Consolidate repeated floods into the named pseudotime column.
  floodPseudotimeProcess(
    urd,
    floods,
    floods.name = cfg$pseudotime_name,
    max.frac.NA = cfg$flood_max_frac_na,
    pseudotime.fun = mean,
    stability.div = cfg$flood_stability_div
  )
}

# Convert a comma-separated sequence of terminal labels to numeric tip IDs.
# Input: tip_labels in the required order; whitespace and empty entries are removed.
# Returns: a table linking character IDs 1, 2, ... to terminal annotations.
# At least two labels are required; tip IDs follow input order, not alphabetic order.
make_tip_mapping <- function(tip_labels) {
  labels <- trimws(strsplit(tip_labels, ",", fixed = TRUE)[[1]])
  labels <- labels[nzchar(labels)]
  if (length(labels) < 2) stop("At least two tip labels are required to build a tree.", call. = FALSE)
  data.frame(
    tip_id = as.character(seq_along(labels)),
    paper_cluster_annotation = labels,
    stringsAsFactors = FALSE
  )
}

# Attach terminal tip IDs and root flags to the URD cell metadata.
# Inputs: annotation column, root flag column or root label, and a tip mapping table.
# Returns: the updated object, root cell IDs, and requested numeric tip IDs.
# Cells outside the supplied terminal groups retain NA tip IDs and remain in the graph.
add_tree_groups <- function(object, annotation_col, root_label, root_col, tip_mapping) {
  if (!(annotation_col %in% colnames(object@meta))) stop("Missing annotation column: ", annotation_col, call. = FALSE)
  meta <- object@meta
  if (nzchar(root_col)) {
    if (!(root_col %in% colnames(meta))) stop("Missing root column: ", root_col, call. = FALSE)
    root_flag <- as.logical(meta[[root_col]])
    root_flag[is.na(root_flag)] <- FALSE
    root_cells <- rownames(meta)[root_flag]
  } else {
    root_cells <- rownames(meta)[meta[[annotation_col]] == root_label]
  }
  if (length(root_cells) == 0) stop("No root cells found for label: ", root_label, call. = FALSE)

  # Nonterminal cells receive NA tip IDs but are retained as intermediate states.
  tip_id_by_label <- setNames(tip_mapping$tip_id, tip_mapping$paper_cluster_annotation)
  tip_ids <- unname(tip_id_by_label[as.character(meta[[annotation_col]])])
  tip_ids[is.na(tip_ids)] <- NA_character_
  names(tip_ids) <- rownames(meta)

  if (!all(tip_mapping$paper_cluster_annotation %in% meta[[annotation_col]])) {
    missing_labels <- setdiff(tip_mapping$paper_cluster_annotation, unique(meta[[annotation_col]]))
    stop("No cells found for tip label(s): ", paste(missing_labels, collapse = ", "), call. = FALSE)
  }

  object@group.ids[rownames(meta), "paper_tree_tip_id"] <- tip_ids[rownames(meta)]
  object@meta[rownames(meta), "paper_tree_tip_id"] <- tip_ids[rownames(meta)]
  object@meta[rownames(meta), "paper_tree_root_candidate"] <- rownames(meta) %in% root_cells
  list(object = object, root_cells = root_cells, tip_ids = tip_mapping$tip_id)
}


# Root and tip selection are supplied explicitly in metadata, before filtering.
# Choose a graph method explicitly; the PCA-kNN method is not a diffusion eigendecomposition.
# Run filtering, variable-gene selection, geometry, and flood pseudotime.
# Inputs: named raw counts (genes x cells), matching metadata, an explicit dm_method,
# and a named parameters list overriding the defaults below. Metadata must contain
# logical urd_root_candidate values without NA; supplied roots are checked after filtering.
# Returns: list(object, parameters, root_cells). The function does not write files.
# The object contains the selected genes, graph, and named flood pseudotime;
# parameters records the resolved settings for downstream analysis.
urd_pseudotime <- function(counts, metadata, dm_method, parameters = list()) {
  if (!requireNamespace("URD", quietly = TRUE) || !requireNamespace("Matrix", quietly = TRUE)) {
    stop("Install URD and Matrix before running this analysis.")
  }
  suppressPackageStartupMessages(library(URD))
  if (missing(dm_method) || !dm_method %in% c("urd_expression", "pca_knn")) {
    stop("Choose dm_method = 'urd_expression' or 'pca_knn' explicitly.")
  }
  if (is.null(rownames(counts)) || is.null(colnames(counts)) ||
      anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts))) {
    stop("Counts require unique gene and cell names.")
  }
  if (anyDuplicated(rownames(metadata)) || !setequal(colnames(counts), rownames(metadata))) {
    stop("Metadata row names must uniquely match count column names.")
  }
  # Align metadata by cell name before passing both inputs to createURD().
  metadata <- metadata[colnames(counts), , drop = FALSE]
  if (!"urd_root_candidate" %in% names(metadata)) stop("Metadata needs urd_root_candidate.")
  if (!is.logical(metadata$urd_root_candidate) || anyNA(metadata$urd_root_candidate)) {
    stop("urd_root_candidate must contain TRUE/FALSE without missing values.")
  }
  # Named overrides replace these analytical defaults; unknown keys are rejected.
  cfg <- check_parameters(list(
    seed = 7L, min_genes = 500L, min_cells = 3L, min_counts = 10L,
    num_variable_genes = 3000L, pca_mp_factor = 2, dm_pca_dims = NA_integer_,
    knn = 100L, sigma_text = "local", n_floods = 20L,
    flood_minimum_cells = 2L, flood_max_frac_na = 0.4,
    flood_stability_div = NULL, pseudotime_name = "root_pseudotime"
  ), parameters)
  cfg$dm_method <- dm_method
  cfg$sigma <- if (toupper(cfg$sigma_text) == "NULL") NULL else {
    value <- suppressWarnings(as.numeric(cfg$sigma_text))
    if (is.na(value)) cfg$sigma_text else value
  }
  if (is.null(cfg$flood_stability_div)) cfg$flood_stability_div <- min(10L, cfg$n_floods)
  # Apply cell/gene filters before feature selection and graph construction.
  created <- create_filtered_urd(list(counts = counts, metadata = metadata), cfg)
  object <- select_variable_genes(created$urd, cfg)
  object <- run_urd_pca(object, cfg)
  object <- run_urd_diffusion_map(object, cfg)
  object <- run_flood_pseudotime(object, created$root_cells, cfg)
  list(object = object, parameters = cfg, root_cells = created$root_cells)
}

# tip_labels order determines numeric tip IDs and must match the intended analysis.
# Reconstruct lineage branching from a graph with calculated pseudotime.
# Inputs: an URD object, terminal annotation column, ordered tip_labels, pseudotime
# column name, logical root column, and optional named tree-parameter overrides.
# Returns: list(object, parameters, tip_mapping, logistic); the object contains the tree.
# Walks start in terminal groups and are biased toward lower pseudotime. Their
# visitation profiles inform the divergence tests used to merge terminal branches.
# Root and terminal definitions are supplied annotations, not inferred cell identities.
urd_lineage <- function(object, annotation_col, tip_labels, pseudotime_name,
                        root_col = "urd_root_candidate", parameters = list()) {
  suppressPackageStartupMessages(library(URD))
  if (anyDuplicated(tip_labels) || length(tip_labels) < 2L) stop("Supply at least two distinct tip labels.")
  if (!pseudotime_name %in% colnames(object@pseudotime)) stop("Pseudotime column not found.")
  if (!root_col %in% names(object@meta) || !is.logical(object@meta[[root_col]]) ||
      anyNA(object@meta[[root_col]])) stop("Root flags must be logical without missing values.")
  # Named overrides replace these analytical defaults; unknown keys are rejected.
  cfg <- check_parameters(list(
    seed = 7L, optimal_cells_forward = 100L, max_cells_back = 40L,
    pseudotime_direction = "<", n_per_tip = 5000L, root_visits = 1L,
    max_steps = ncol(object@logupx.data), process_n_subsample = 10L,
    cells_per_pseudotime_bin = 80L, bins_per_pseudotime_window = 5L,
    minimum_visits = 10L, visit_threshold = 0.7, divergence_method = "ks",
    p_thresh = 0.01, min_cells_per_segment = 1L,
    min_pseudotime_per_segment = 0.01, dendro_node_size = 100L
  ), parameters)
  set.seed(cfg$seed)
  # Keep the supplied endpoint order when assigning URD numeric tip IDs.
  mapping <- make_tip_mapping(paste(tip_labels, collapse = ","))
  grouped <- add_tree_groups(object, annotation_col, "", root_col, mapping)
  object <- loadTipCells(grouped$object, tips = "paper_tree_tip_id")
  # Fit the allowed pseudotime movement using forward/backward neighborhood sizes.
  logistic <- pseudotimeDetermineLogistic(
    object, pseudotime = pseudotime_name,
    optimal.cells.forward = cfg$optimal_cells_forward, max.cells.back = cfg$max_cells_back,
    pseudotime.direction = cfg$pseudotime_direction, do.plot = FALSE, print.values = FALSE
  )
  # Bias graph transitions toward the root using the fitted logistic curve.
  transition <- pseudotimeWeightTransitionMatrix(
    object, pseudotime = pseudotime_name, logistic.params = logistic,
    pseudotime.direction = cfg$pseudotime_direction, verbose = TRUE
  )
  # Walk from each terminal group until root visits or the step limit is reached.
  walks <- simulateRandomWalksFromTips(
    object, tip.group.id = "paper_tree_tip_id", root.cells = grouped$root_cells,
    transition.matrix = transition, n.per.tip = cfg$n_per_tip,
    root.visits = cfg$root_visits, max.steps = cfg$max_steps, verbose = TRUE
  )
  # Convert simulated paths to tip-specific visitation profiles for each cell.
  object <- processRandomWalksFromTips(object, walks.list = walks,
                                     n.subsample = cfg$process_n_subsample, verbose = TRUE)
  # Compare visitation profiles across pseudotime windows to join branches.
  # Segment size, visit support, and divergence thresholds come from cfg.
  object <- buildTree(
    object, pseudotime = pseudotime_name, tips.use = grouped$tip_ids,
    divergence.method = cfg$divergence_method, weighted.fusion = TRUE,
    use.only.original.tips = TRUE, cells.per.pseudotime.bin = cfg$cells_per_pseudotime_bin,
    bins.per.pseudotime.window = cfg$bins_per_pseudotime_window,
    minimum.visits = cfg$minimum_visits, visit.threshold = cfg$visit_threshold,
    save.breakpoint.plots = NULL, save.all.breakpoint.info = TRUE,
    p.thresh = cfg$p_thresh, min.cells.per.segment = cfg$min_cells_per_segment,
    min.pseudotime.per.segment = cfg$min_pseudotime_per_segment,
    dendro.node.size = cfg$dendro_node_size, verbose = TRUE
  )
  list(object = object, parameters = cfg, tip_mapping = mapping, logistic = logistic)
}
