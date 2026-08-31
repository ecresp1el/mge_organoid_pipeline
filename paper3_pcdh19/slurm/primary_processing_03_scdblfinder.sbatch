#!/usr/bin/env bash
# PURPOSE
#   Execute frozen Step 03 scDblFinder detection on Great Lakes.
#
# BOUNDARY
#   One GEX_1 capture; clusters=TRUE; dbr.sd=1; no manual dbr; otherwise model
#   defaults. Preserve calls and diagnostics; remove no cells and STOP in review.
set -Eeuo pipefail

: "${PAPER3_PRIMARY_RUN_DIR:?Submit through submit_primary_processing_step_03.sh}"
: "${PAPER3_PRIMARY_WORKFLOW_ROOT:?Workflow root is required}"
: "${PAPER3_PRIMARY_PYTHON_BIN:?Python interpreter is required}"
RUN_DIR="${PAPER3_PRIMARY_RUN_DIR}"

failed() {
  local line="$1" status="$2"
  printf 'COMPUTATION_FAILED\nstep=03_scdblfinder\njob_id=%s\nline=%s\nexit_status=%s\ntime=%s\n' "${SLURM_JOB_ID:-unknown}" "${line}" "${status}" "$(date --iso-8601=seconds)" > "${RUN_DIR}/COMPUTATION_FAILED.txt"
}
trap 'failed "${LINENO}" "$?"' ERR

# shellcheck disable=SC1090
source "${RUN_DIR}/config/resolved.env"
module load "${PRIMARY_PROCESSING_STEP03_R_MODULE}"
export R_LIBS_USER="${PRIMARY_PROCESSING_STEP03_R_LIBRARY}"
export PYTHONPATH="${RUN_DIR}/code"
export PYTHONHASHSEED=0
export LC_ALL=C
export TZ=UTC
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MPLCONFIGDIR="${RUN_DIR}/provenance/matplotlib-cache"
mkdir -p "${MPLCONFIGDIR}"

{
  echo "started=$(date --iso-8601=seconds)"
  echo "host=$(hostname)"
  echo "slurm_job_id=${SLURM_JOB_ID:-unknown}"
  echo "slurm_cpus=${SLURM_CPUS_PER_TASK:-unknown}"
  echo "python=$(${PAPER3_PRIMARY_PYTHON_BIN} --version 2>&1)"
  echo "r=$(Rscript --version 2>&1)"
  echo "r_libs_user=${R_LIBS_USER}"
} > "${RUN_DIR}/provenance/runtime_environment.txt"

COMMON=(
  --input-h5ad "${INPUT_H5AD}"
  --input-status "${RUN_DIR}/config/input_step02_STEP_STATUS.tsv"
  --input-manifest "${RUN_DIR}/config/input_step02_output_manifest.tsv"
  --approval-ledger-snapshot "${RUN_DIR}/config/approval_ledger_at_submission.tsv"
  --run-dir "${RUN_DIR}"
  --workflow-root "${WORKFLOW_ROOT}"
  --frozen-code-dir "${RUN_DIR}/code/primary_processing"
  --expected-step02-run-id "${PRIMARY_PROCESSING_STEP03_STEP02_RUN_ID}"
  --expected-input-bytes "${PRIMARY_PROCESSING_STEP03_INPUT_BYTES}"
  --expected-input-sha256 "${PRIMARY_PROCESSING_STEP03_INPUT_SHA256}"
  --expected-cells "${PRIMARY_PROCESSING_STEP03_EXPECTED_CELLS}"
  --expected-genes "${PRIMARY_PROCESSING_STEP03_EXPECTED_GENES}"
  --capture-id "${PRIMARY_PROCESSING_STEP03_CAPTURE_ID}"
  --sample-field "${PRIMARY_PROCESSING_STEP03_SAMPLE_FIELD}"
  --design-field "${PRIMARY_PROCESSING_STEP03_DESIGN_FIELD}"
  --primary-seed "${PRIMARY_PROCESSING_STEP03_PRIMARY_SEED}"
  --reproducibility-seed "${PRIMARY_PROCESSING_STEP03_REPRODUCIBILITY_SEED}"
  --plot-dpi "${PRIMARY_PROCESSING_STEP03_PLOT_DPI}"
  --compression "${PRIMARY_PROCESSING_STEP03_COMPRESSION}"
)

"${PAPER3_PRIMARY_PYTHON_BIN}" -m primary_processing.step03_cli prepare "${COMMON[@]}"

Rscript "${RUN_DIR}/code/primary_processing/step03_scdblfinder.R" \
  --bridge-h5 "${RUN_DIR}/intermediate/step02_counts_for_scdblfinder.h5" \
  --metadata-tsv "${RUN_DIR}/intermediate/step03_cell_metadata.tsv" \
  --output-dir "${RUN_DIR}/intermediate/r_output" \
  --expected-cells "${PRIMARY_PROCESSING_STEP03_EXPECTED_CELLS}" \
  --expected-genes "${PRIMARY_PROCESSING_STEP03_EXPECTED_GENES}" \
  --capture-id "${PRIMARY_PROCESSING_STEP03_CAPTURE_ID}" \
  --primary-seed "${PRIMARY_PROCESSING_STEP03_PRIMARY_SEED}" \
  --reproducibility-seed "${PRIMARY_PROCESSING_STEP03_REPRODUCIBILITY_SEED}"

"${PAPER3_PRIMARY_PYTHON_BIN}" -m primary_processing.step03_cli finalize "${COMMON[@]}"

printf 'COMPUTATION_SUCCESS\nstep=03_scdblfinder\nstatus=IN_REVIEW\ncapture=GEX_1\nclusters=TRUE\ndbr.sd=1\ndbr=not_supplied\ncells_removed=0\njob_id=%s\ncompleted=%s\n' "${SLURM_JOB_ID:-unknown}" "$(date --iso-8601=seconds)" > "${RUN_DIR}/COMPUTATION_SUCCESS.txt"
rm -f "${RUN_DIR}/COMPUTATION_FAILED.txt"
