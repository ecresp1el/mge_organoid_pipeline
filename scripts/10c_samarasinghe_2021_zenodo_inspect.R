#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

get_arg <- function(args, key, default = "") {
  hit <- which(args == key)
  if (length(hit) > 0 && hit[1] < length(args)) return(args[hit[1] + 1])
  prefix <- paste0(key, "=")
  hit <- grep(paste0("^", prefix), args)
  if (length(hit) > 0) return(sub(prefix, "", args[hit[1]]))
  default
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

args <- commandArgs(trailingOnly = TRUE)
project_root <- get_arg(args, "--project-root", Sys.getenv("PROJECT_ROOT", unset = ""))
if (!nzchar(project_root)) {
  project_root <- "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder"
}
project_root <- sub("/+$", "", project_root)

rdata <- get_arg(
  args,
  "--rdata",
  file.path(project_root, "data/raw/samarasinghe_2021_zenodo/Samarasinghe_2021_seurat_object.RData")
)
out_dir <- get_arg(
  args,
  "--out-dir",
  file.path(project_root, "results/samarasinghe_2021_zenodo_inspection")
)
stopifnot(file.exists(rdata))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading Zenodo RData: ", rdata)
loaded <- load(rdata)
message("Loaded object names: ", paste(loaded, collapse = ", "))

object_inventory <- do.call(
  rbind,
  lapply(loaded, function(nm) {
    x <- get(nm)
    data.frame(
      object_name = nm,
      class = paste(class(x), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
)
write_tsv(object_inventory, file.path(out_dir, "loaded_objects.tsv"))

seurat_names <- loaded[vapply(loaded, function(nm) inherits(get(nm), "Seurat"), logical(1))]
if (length(seurat_names) == 0) stop("No Seurat object found in RData")

obj_name <- seurat_names[[1]]
obj <- get(obj_name)
message("Inspecting Seurat object: ", obj_name)

meta <- obj[[]]
assays <- Assays(obj)
reductions <- Reductions(obj)

inventory <- data.frame(
  object_name = obj_name,
  n_cells = ncol(obj),
  n_features_default_assay = nrow(obj),
  default_assay = DefaultAssay(obj),
  assays = paste(assays, collapse = ","),
  reductions = paste(reductions, collapse = ","),
  n_metadata_columns = ncol(meta),
  stringsAsFactors = FALSE
)
write_tsv(inventory, file.path(out_dir, "seurat_inventory.tsv"))

metadata_columns <- data.frame(metadata_column = colnames(meta), stringsAsFactors = FALSE)
write_tsv(metadata_columns, file.path(out_dir, "metadata_columns.tsv"))

metadata_summary <- do.call(
  rbind,
  lapply(colnames(meta), function(col) {
    vals <- meta[[col]]
    data.frame(
      metadata_column = col,
      class = paste(class(vals), collapse = ";"),
      n_unique = length(unique(vals)),
      n_missing = sum(is.na(vals)),
      example_values = paste(head(unique(as.character(vals)), 12), collapse = "|"),
      stringsAsFactors = FALSE
    )
  })
)
write_tsv(metadata_summary, file.path(out_dir, "metadata_summary.tsv"))

for (col in colnames(meta)) {
  vals <- meta[[col]]
  if (length(unique(vals)) <= 100) {
    tab <- as.data.frame(table(vals, useNA = "ifany"), stringsAsFactors = FALSE)
    colnames(tab) <- c(col, "n_cells")
    write_tsv(tab, file.path(out_dir, paste0("metadata_counts__", make.names(col), ".tsv")))
  }
}

reduction_dims <- if (length(reductions) > 0) {
  do.call(
    rbind,
    lapply(reductions, function(red) {
      emb <- Embeddings(obj, red)
      data.frame(
        reduction = red,
        n_cells = nrow(emb),
        n_dims = ncol(emb),
        dim_names = paste(colnames(emb), collapse = ","),
        stringsAsFactors = FALSE
      )
    })
  )
} else {
  data.frame(reduction = character(), n_cells = integer(), n_dims = integer(), dim_names = character())
}
write_tsv(reduction_dims, file.path(out_dir, "reduction_dimensions.tsv"))

assay_dims <- do.call(
  rbind,
  lapply(assays, function(assay_name) {
    assay_obj <- obj[[assay_name]]
    slot_dims <- lapply(c("counts", "data", "scale.data"), function(slot_name) {
      if (slot_name %in% slotNames(assay_obj)) {
        dim(slot(assay_obj, slot_name))
      } else {
        c(NA_integer_, NA_integer_)
      }
    })
    names(slot_dims) <- c("counts", "data", "scale.data")
    data.frame(
      assay = assay_name,
      counts_features = slot_dims$counts[[1]],
      counts_cells = slot_dims$counts[[2]],
      data_features = slot_dims$data[[1]],
      data_cells = slot_dims$data[[2]],
      scale_features = slot_dims$scale.data[[1]],
      scale_cells = slot_dims$scale.data[[2]],
      stringsAsFactors = FALSE
    )
  })
)
write_tsv(assay_dims, file.path(out_dir, "assay_dimensions.tsv"))

message("Inspection complete: ", out_dir)
print(inventory)
print(reduction_dims)
