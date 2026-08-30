#!/usr/bin/env bash
# PURPOSE
#   Execute the frozen primary-processing Step 00 Python package on Great Lakes.
#
# SCIENTIFIC BOUNDARY
#   Construct and structurally validate the canonical raw-count AnnData only.
#   Do not calculate QC metrics, filter, normalize, reduce, cluster, or annotate.
#
# RESOURCES AND LOGS
#   The submission wrapper assigns resources and scheduler log paths so the
#   exact same frozen script supports versioned and guarded replacement runs.
set -Eeuo pipefail

: "${PAPER3_PRIMARY_RUN_DIR:?Submit through submit_primary_processing_step_00.sh}"
: "${PAPER3_PRIMARY_WORKFLOW_ROOT:?Primary-processing workflow root is required}"
: "${PAPER3_PRIMARY_CELLRANGER_ROOT:?Cell Ranger root is required}"
: "${PAPER3_PRIMARY_PYTHON_BIN:?Python interpreter is required}"

RUN_DIR="${PAPER3_PRIMARY_RUN_DIR}"
failed() {
  local line="$1" status="$2"
  printf 'COMPUTATION_FAILED\nstep=00_input_validation_and_canonical_anndata\njob_id=%s\nline=%s\nexit_status=%s\ntime=%s\n' \
    "${SLURM_JOB_ID:-unknown}" "${line}" "${status}" "$(date --iso-8601=seconds)" \
    > "${RUN_DIR}/COMPUTATION_FAILED.txt"
}
trap 'failed "${LINENO}" "$?"' ERR

# shellcheck disable=SC1090
source "${RUN_DIR}/config/resolved.env"
export PYTHONPATH="${RUN_DIR}/code"
export PYTHONHASHSEED=0
export LC_ALL=C
export TZ=UTC

{
  echo "started=$(date --iso-8601=seconds)"
  echo "host=$(hostname)"
  echo "slurm_job_id=${SLURM_JOB_ID:-unknown}"
  echo "slurm_cpus=${SLURM_CPUS_PER_TASK:-unknown}"
  echo "python=$(${PAPER3_PRIMARY_PYTHON_BIN} --version 2>&1)"
  echo "python_executable=${PAPER3_PRIMARY_PYTHON_BIN}"
} > "${RUN_DIR}/provenance/runtime_environment.txt"

"${PAPER3_PRIMARY_PYTHON_BIN}" -m primary_processing.cli \
  --cellranger-root "${PAPER3_PRIMARY_CELLRANGER_ROOT}" \
  --sample-key "${RUN_DIR}/config/sample_key.csv" \
  --technical-manifest "${RUN_DIR}/config/sample_manifest.tsv" \
  --run-dir "${RUN_DIR}" \
  --workflow-root "${PAPER3_PRIMARY_WORKFLOW_ROOT}" \
  --frozen-code-dir "${RUN_DIR}/code/primary_processing" \
  --expected-samples "${PRIMARY_PROCESSING_EXPECTED_SAMPLES}" \
  --expected-cells "${PRIMARY_PROCESSING_EXPECTED_CELLS}" \
  --expected-genes "${PRIMARY_PROCESSING_EXPECTED_GENES}" \
  --expected-genome "${PRIMARY_PROCESSING_EXPECTED_GENOME}" \
  --feature-type "${PRIMARY_PROCESSING_FEATURE_TYPE}" \
  --compression "${PRIMARY_PROCESSING_STEP00_COMPRESSION}"

printf 'COMPUTATION_SUCCESS\nstep=00_input_validation_and_canonical_anndata\nstatus=IN_REVIEW\njob_id=%s\ncompleted=%s\n' \
  "${SLURM_JOB_ID:-unknown}" "$(date --iso-8601=seconds)" \
  > "${RUN_DIR}/COMPUTATION_SUCCESS.txt"
rm -f "${RUN_DIR}/COMPUTATION_FAILED.txt"
