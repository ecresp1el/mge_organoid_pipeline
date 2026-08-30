#!/usr/bin/env bash
# PURPOSE
#   Execute frozen Scanpy Step 01 QC metrics on Great Lakes.
#
# BOUNDARY
#   Add descriptive QC metadata and figures only. Never filter, threshold,
#   normalize, reduce, cluster, or annotate.
set -Eeuo pipefail

: "${PAPER3_PRIMARY_RUN_DIR:?Submit through submit_primary_processing_step_01.sh}"
: "${PAPER3_PRIMARY_WORKFLOW_ROOT:?Workflow root is required}"
: "${PAPER3_PRIMARY_PYTHON_BIN:?Python interpreter is required}"
RUN_DIR="${PAPER3_PRIMARY_RUN_DIR}"

failed() {
  local line="$1" status="$2"
  printf 'COMPUTATION_FAILED\nstep=01_qc_metrics\njob_id=%s\nline=%s\nexit_status=%s\ntime=%s\n' "${SLURM_JOB_ID:-unknown}" "${line}" "${status}" "$(date --iso-8601=seconds)" > "${RUN_DIR}/COMPUTATION_FAILED.txt"
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

"${PAPER3_PRIMARY_PYTHON_BIN}" -m primary_processing.step01_cli \
  --input-h5ad "${INPUT_H5AD}" \
  --input-step-status "${RUN_DIR}/config/approved_step00_STEP_STATUS.tsv" \
  --input-output-manifest "${RUN_DIR}/config/approved_step00_output_manifest.tsv" \
  --approval-ledger-snapshot "${RUN_DIR}/config/approval_ledger_at_submission.tsv" \
  --run-dir "${RUN_DIR}" \
  --workflow-root "${WORKFLOW_ROOT}" \
  --frozen-code-dir "${RUN_DIR}/code/primary_processing" \
  --expected-input-run-id "${PRIMARY_PROCESSING_STEP01_INPUT_RUN_ID}" \
  --expected-input-bytes "${PRIMARY_PROCESSING_STEP01_INPUT_BYTES}" \
  --expected-input-sha256 "${PRIMARY_PROCESSING_STEP01_INPUT_SHA256}" \
  --expected-cells "${PRIMARY_PROCESSING_STEP01_EXPECTED_CELLS}" \
  --expected-genes "${PRIMARY_PROCESSING_STEP01_EXPECTED_GENES}" \
  --mitochondrial-prefix "${PRIMARY_PROCESSING_STEP01_MITO_PREFIX}" \
  --ribosomal-prefixes "${PRIMARY_PROCESSING_STEP01_RIBO_PREFIXES}" \
  --percent-top "${PRIMARY_PROCESSING_STEP01_PERCENT_TOP}" \
  --histogram-bins "${PRIMARY_PROCESSING_STEP01_HISTOGRAM_BINS}" \
  --plot-dpi "${PRIMARY_PROCESSING_STEP01_PLOT_DPI}" \
  --compression "${PRIMARY_PROCESSING_STEP01_COMPRESSION}"

printf 'COMPUTATION_SUCCESS\nstep=01_qc_metrics\nstatus=IN_REVIEW\njob_id=%s\ncompleted=%s\n' "${SLURM_JOB_ID:-unknown}" "$(date --iso-8601=seconds)" > "${RUN_DIR}/COMPUTATION_SUCCESS.txt"
rm -f "${RUN_DIR}/COMPUTATION_FAILED.txt"
