#!/usr/bin/env bash
set -euo pipefail

# Render the DIV90 cluster/pseudotime context and two marker modules.
# Arguments: a URD tree RDS file and the output folder. Rscript must be on PATH.

if [[ $# -ne 2 ]]; then
  printf 'Usage: bash render_div90_figure.sh <tree.rds> <output_folder>\n' >&2
  exit 1
fi

# Locate the renderer relative to this script, independent of the working directory.
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# Use log1p(CP10K) marker colours with a floor of 1 and gene-specific 99th-percentile
# upper limits. The neuronal marker panels exclude annotated clusters 8 and 10;
# the top context panels retain the full displayed tree. Export one 300 dpi PNG.
Rscript --vanilla "${script_dir}/plot_urd.R" \
  --tree-rds "$1" \
  --outdir "$2" \
  --genes "HES1,NKX2-1,LHX6,LHX8,CRABP1,KCNC1,PLXNA2,NRP1,NRP2,ROBO1,ERBB4,CXCR4,ACKR3,SEMA3A,SEMA3F" \
  --gene-labels "HES1,NKX2-1,LHX6,LHX8,CRABP1,KCNC1,PLXNA2,NRP1,NRP2,ROBO1,ERBB4,CXCR4,CXCR7 (ACKR3),SEMA3A,SEMA3F" \
  --annotation-col cluster_number_name \
  --pseudotime-name div90_jia_rootscore_root \
  --expression-color-floor 1 \
  --vmax-quantile 0.99 \
  --point-size 0.26 \
  --marker-ncol 3 \
  --plot-width 16 \
  --plot-height 13 \
  --dpi 300 \
  --cluster-style annotated \
  --layout-style side_by_side \
  --marker-exclude-clusters "8,10" \
  --output-formats png
