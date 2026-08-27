#!/usr/bin/env bash

set -Eeuo pipefail

: "${PAPER2_GENE_RUN_DIR:?run directory is required}"
: "${CONDA_ENV_BIN:?Python environment bin directory is required}"
: "${MAKE_SVG:=false}"
: "${FIGURE_DPI:=300}"

failed() {
  printf 'FAILED\njob_id=%s\nline=%s\nexit_status=%s\ntime=%s\n' \
    "${SLURM_JOB_ID:-unknown}" "$1" "$2" "$(date --iso-8601=seconds)" \
    > "${PAPER2_GENE_RUN_DIR}/FIGURES_FAILED.txt"
}
trap 'failed "${LINENO}" "$?"' ERR

echo "Step 02 plot-only overlap rendering started: $(date --iso-8601=seconds)"
"${CONDA_ENV_BIN}/python" \
  "${PAPER2_GENE_RUN_DIR}/code/02_render_gene_overlap_figures.py" \
  --run-dir "${PAPER2_GENE_RUN_DIR}" \
  --dpi "${FIGURE_DPI}" \
  --svg "${MAKE_SVG}"

[[ -s "${PAPER2_GENE_RUN_DIR}/FIGURES_SUCCESS.txt" ]]
[[ "$(find "${PAPER2_GENE_RUN_DIR}/figures/png" -type f -name '*.png' | wc -l)" -eq 4 ]]
[[ "$(find "${PAPER2_GENE_RUN_DIR}/figures/pdf" -type f -name '*.pdf' | wc -l)" -eq 4 ]]
rm -f "${PAPER2_GENE_RUN_DIR}/FIGURES_FAILED.txt"
echo "Step 02 plot-only overlap rendering completed: $(date --iso-8601=seconds)"
