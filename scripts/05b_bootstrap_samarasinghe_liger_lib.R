#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(utils)
})

parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]
    if (startsWith(arg, "--")) {
      keyval <- sub("^--", "", arg)
      if (grepl("=", keyval, fixed = TRUE)) {
        kv <- strsplit(keyval, "=", fixed = TRUE)[[1]]
        out[[kv[[1]]]] <- if (length(kv) > 1) kv[[2]] else ""
      } else {
        out[[keyval]] <- if (i < length(args) && !startsWith(args[[i + 1]], "--")) args[[i + 1]] else ""
        if (i < length(args) && !startsWith(args[[i + 1]], "--")) i <- i + 1
      }
    }
    i <- i + 1
  }
  out
}

arg_value <- function(args, name, default = "") {
  value <- args[[name]]
  if (!is.null(value) && nzchar(value)) value else default
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
project_root <- arg_value(args, "project-root", Sys.getenv("PROJECT_ROOT", unset = ""))
if (!nzchar(project_root)) stop("PROJECT_ROOT or --project-root is required")
project_root <- sub("/+$", "", project_root)

lib_path <- arg_value(
  args,
  "lib-path",
  file.path(project_root, "software/Rlibs/samarasinghe_liger_R4.4")
)
dir.create(lib_path, recursive = TRUE, showWarnings = FALSE)
unlink(Sys.glob(file.path(lib_path, "00LOCK*")), recursive = TRUE, force = TRUE)
.libPaths(c(lib_path, .libPaths()))

options(repos = c(CRAN = "https://cloud.r-project.org"))
message("Bootstrap library: ", lib_path)
message("R version: ", paste(R.version$major, R.version$minor, sep = "."))
message(".libPaths():\n  ", paste(.libPaths(), collapse = "\n  "))

install_cran_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing CRAN package: ", pkg)
    install.packages(pkg, lib = lib_path, dependencies = c("Depends", "Imports", "LinkingTo"))
  } else {
    message("Already available: ", pkg, " ", as.character(packageVersion(pkg)))
  }
}

install_cran_if_missing("remotes")
install_cran_if_missing("BiocManager")
install_cran_if_missing("R.utils")

bioc_pkgs <- c("S4Vectors", "DelayedArray", "HDF5Array")
missing_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc) > 0) {
  message("Installing Bioconductor package(s): ", paste(missing_bioc, collapse = ", "))
  BiocManager::install(
    missing_bioc,
    lib = lib_path,
    ask = FALSE,
    update = FALSE,
    dependencies = c("Depends", "Imports", "LinkingTo")
  )
}

if (!requireNamespace("rliger", quietly = TRUE)) {
  message("Installing rliger into isolated library")
  install.packages("rliger", lib = lib_path, dependencies = c("Depends", "Imports", "LinkingTo"))
} else {
  message("Already available: rliger ", as.character(packageVersion("rliger")))
}

if (!requireNamespace("SeuratWrappers", quietly = TRUE)) {
  message("Installing SeuratWrappers from satijalab/seurat-wrappers into isolated library")
  remotes::install_github(
    "satijalab/seurat-wrappers",
    lib = lib_path,
    dependencies = c("Depends", "Imports", "LinkingTo"),
    upgrade = "never",
    build_vignettes = FALSE
  )
} else {
  message("Already available: SeuratWrappers ", as.character(packageVersion("SeuratWrappers")))
}

required <- c("SeuratWrappers", "rliger")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Bootstrap incomplete; missing: ", paste(missing, collapse = ", "))
}

message("Bootstrap complete.")
for (pkg in required) {
  message(pkg, " ", as.character(packageVersion(pkg)))
}
