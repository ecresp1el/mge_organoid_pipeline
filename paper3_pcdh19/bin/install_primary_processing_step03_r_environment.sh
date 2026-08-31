#!/usr/bin/env bash
# PURPOSE
#   Recreate or verify the project-scoped native-R environment for Step 03.
#
# OUTPUT
#   Installs the declared Bioconductor 3.22 packages under the Step 03 turbo
#   library without updating unrelated packages, then writes a full inventory.
set -Eeuo pipefail

BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
# shellcheck disable=SC1090
source "${BUNDLE_DIR}/config/primary_processing_step03.env"

module load "${PRIMARY_PROCESSING_STEP03_R_MODULE}"
export R_LIBS_USER="${PRIMARY_PROCESSING_STEP03_R_LIBRARY}"
mkdir -p "${R_LIBS_USER}"

Rscript - "${R_LIBS_USER}" <<'RS'
arguments <- commandArgs(trailingOnly = TRUE)
library_path <- arguments[[1]]
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org", lib = library_path)
}
BiocManager::install(version = "3.22", ask = FALSE, update = FALSE)
BiocManager::install(
  c("scDblFinder", "DropletUtils", "HDF5Array", "zellkonverter", "ggplot2", "patchwork", "data.table"),
  lib = library_path,
  ask = FALSE,
  update = FALSE,
  Ncpus = 8
)
required <- c("scDblFinder", "SingleCellExperiment", "DropletUtils", "data.table", "R6", "rhdf5", "BiocParallel")
missing <- required[!vapply(required, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing)) stop("Missing required packages after installation: ", paste(missing, collapse = ", "))
inventory <- as.data.frame(installed.packages()[, c("Package", "Version", "LibPath")])
write.table(inventory, file.path(dirname(library_path), "r_package_inventory.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
RS

echo "Verified Step 03 R environment: ${R_LIBS_USER}"
