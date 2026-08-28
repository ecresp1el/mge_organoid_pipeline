#!/usr/bin/env bash
# PRELIMINARY Step 09: filtered-cell gene counts -> sample/state pseudobulk -> edgeR.
set -Eeuo pipefail

BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
source "${BUNDLE_DIR}/config/greatlakes.env"
: "${PAPER3_ROOT:?PAPER3_ROOT is required}"

PYTHON_BIN="${PAPER3_PYTHON_BIN:-python3}"
PYTHON_SITE="${PAPER3_ROOT}/software/gse94641_e15_5_label_transfer/python3.6-site-packages"
PYTHONPATH="${PYTHON_SITE}" "${PYTHON_BIN}" -c 'import h5py,numpy,scipy; assert h5py.__version__=="3.1.0"; assert numpy.__version__=="1.19.5"; assert scipy.__version__=="1.5.4"'

module load R/4.4.3
R_SCRIPT_BIN="$(command -v Rscript)"
R_LIBRARY="${PAPER3_ROOT}/software/step_09_pcdh19_preliminary_pseudobulk_differential_expression/R-4.4-library"
export R_LIBS_USER="${R_LIBRARY}"
"${R_SCRIPT_BIN}" -e 'stopifnot(as.character(packageVersion("edgeR"))=="4.4.2", as.character(packageVersion("limma"))=="3.62.2")'

LOCK="${BUNDLE_DIR}/config/step_09_pcdh19_preliminary_pseudobulk_de.lock.json"
REQUIREMENTS="${BUNDLE_DIR}/config/step_09_pcdh19_preliminary_pseudobulk_de.requirements.txt"
SCRIPT="${BUNDLE_DIR}/scripts/Step_09_PCDH19_Preliminary_Pseudobulk_Differential_Expression.py"
R_MODULE="${BUNDLE_DIR}/scripts/step_09_preliminary_edger_pseudobulk.R"
LABEL_ROOT="${PAPER3_ROOT}/results/mge_reference_mapping_gse94641/label_transfer_e15_5"
OUTPUT_ROOT="${PAPER3_ROOT}/results/step_09_pcdh19_preliminary_pseudobulk_differential_expression"
LOG_ROOT="${PAPER3_ROOT}/logs/step_09_pcdh19_preliminary_pseudobulk_differential_expression"
mkdir -p "${LOG_ROOT}"
RUN_STAMP="$(date -u +%Y%m%dT%H%M%SZ)"; LOG_PATH="${LOG_ROOT}/run_${RUN_STAMP}_$$.log"
export PYTHONHASHSEED=0 LC_ALL=C TZ=UTC

{
  echo "analysis_status=PRELIMINARY"
  echo "started_utc=$(date -u --iso-8601=seconds)"
  echo "host=$(hostname)"
  echo "step=09_pcdh19_preliminary_pseudobulk_differential_expression"
  echo "output_root=${OUTPUT_ROOT}"
  PYTHONPATH="${PYTHON_SITE}" "${PYTHON_BIN}" "${SCRIPT}" \
    --lock "${LOCK}" --requirements "${REQUIREMENTS}" \
    --sample-key "${BUNDLE_DIR}/config/sample_key.csv" \
    --label-root "${LABEL_ROOT}" --rscript "${R_SCRIPT_BIN}" \
    --r-module "${R_MODULE}" --output-root "${OUTPUT_ROOT}" "$@"
  echo "completed_utc=$(date -u --iso-8601=seconds)"
} 2>&1 | tee "${LOG_PATH}"

