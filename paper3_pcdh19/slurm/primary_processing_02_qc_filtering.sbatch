#!/usr/bin/env bash
# PURPOSE
#   Execute frozen Step 02 per-sample 5-MAD filtering on Great Lakes.
#
# BOUNDARY
#   Remove only the approved low-count/low-gene/high-mitochondrial union.
#   Stop before any scDblFinder calculation or doublet-based filtering.
set -Eeuo pipefail

: "${PAPER3_PRIMARY_RUN_DIR:?Submit through submit_primary_processing_step_02.sh}"
: "${PAPER3_PRIMARY_WORKFLOW_ROOT:?Workflow root is required}"
: "${PAPER3_PRIMARY_PYTHON_BIN:?Python interpreter is required}"
RUN_DIR="${PAPER3_PRIMARY_RUN_DIR}"

failed() {
  local line="$1" status="$2"
  printf 'COMPUTATION_FAILED\nstep=02_qc_filtering\njob_id=%s\nline=%s\nexit_status=%s\ntime=%s\n' "${SLURM_JOB_ID:-unknown}" "${line}" "${status}" "$(date --iso-8601=seconds)" > "${RUN_DIR}/COMPUTATION_FAILED.txt"
}
trap 'failed "${LINENO}" "$?"' ERR

# shellcheck disable=SC1090
source "${RUN_DIR}/config/resolved.env"
export PYTHONPATH="${RUN_DIR}/code"
export PYTHONHASHSEED=0
export LC_ALL=C
export TZ=UTC
export MPLCONFIGDIR="${RUN_DIR}/provenance/matplotlib-cache"
mkdir -p "${MPLCONFIGDIR}"
{
  echo "started=$(date --iso-8601=seconds)"
  echo "host=$(hostname)"
  echo "slurm_job_id=${SLURM_JOB_ID:-unknown}"
  echo "slurm_cpus=${SLURM_CPUS_PER_TASK:-unknown}"
  echo "python=$(${PAPER3_PRIMARY_PYTHON_BIN} --version 2>&1)"
} > "${RUN_DIR}/provenance/runtime_environment.txt"

"${PAPER3_PRIMARY_PYTHON_BIN}" -m primary_processing.step02_cli \
  --input-h5ad "${INPUT_H5AD}" \
  --input-step01-status "${RUN_DIR}/config/input_step01_STEP_STATUS.tsv" \
  --input-step01-manifest "${RUN_DIR}/config/input_step01_output_manifest.tsv" \
  --input-flags "${INPUT_FLAGS}" \
  --input-step01a-status "${RUN_DIR}/config/input_step01a_STEP_STATUS.tsv" \
  --input-step01a-manifest "${RUN_DIR}/config/input_step01a_output_manifest.tsv" \
  --approval-ledger-snapshot "${RUN_DIR}/config/approval_ledger_at_submission.tsv" \
  --run-dir "${RUN_DIR}" \
  --workflow-root "${WORKFLOW_ROOT}" \
  --frozen-code-dir "${RUN_DIR}/code/primary_processing" \
  --expected-step01-run-id "${PRIMARY_PROCESSING_STEP02_STEP01_RUN_ID}" \
  --expected-step01-bytes "${PRIMARY_PROCESSING_STEP02_STEP01_BYTES}" \
  --expected-step01-sha256 "${PRIMARY_PROCESSING_STEP02_STEP01_SHA256}" \
  --expected-step01a-run-id "${PRIMARY_PROCESSING_STEP02_STEP01A_RUN_ID}" \
  --expected-flags-bytes "${PRIMARY_PROCESSING_STEP02_FLAGS_BYTES}" \
  --expected-flags-sha256 "${PRIMARY_PROCESSING_STEP02_FLAGS_SHA256}" \
  --expected-cells-before "${PRIMARY_PROCESSING_STEP02_EXPECTED_CELLS_BEFORE}" \
  --expected-cells-after "${PRIMARY_PROCESSING_STEP02_EXPECTED_CELLS_AFTER}" \
  --expected-removed "${PRIMARY_PROCESSING_STEP02_EXPECTED_REMOVED}" \
  --expected-genes "${PRIMARY_PROCESSING_STEP02_EXPECTED_GENES}" \
  --selected-stringency "${PRIMARY_PROCESSING_STEP02_SELECTED_STRINGENCY}" \
  --sample-field "${PRIMARY_PROCESSING_STEP02_SAMPLE_FIELD}" \
  --design-field "${PRIMARY_PROCESSING_STEP02_DESIGN_FIELD}" \
  --plot-dpi "${PRIMARY_PROCESSING_STEP02_PLOT_DPI}" \
  --compression "${PRIMARY_PROCESSING_STEP02_COMPRESSION}"

printf 'COMPUTATION_SUCCESS\nstep=02_qc_filtering\nstatus=IN_REVIEW\nfilter=per_sample_5mad_union\nscDblFinder_run=no\njob_id=%s\ncompleted=%s\n' "${SLURM_JOB_ID:-unknown}" "$(date --iso-8601=seconds)" > "${RUN_DIR}/COMPUTATION_SUCCESS.txt"
rm -f "${RUN_DIR}/COMPUTATION_FAILED.txt"
