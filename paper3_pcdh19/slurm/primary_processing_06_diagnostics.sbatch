#!/usr/bin/env bash
# PURPOSE
#   Execute the frozen, unintegrated Step 06 A-L diagnostic workflow.
#
# BOUNDARY
#   Consume approved Step 02 directly; use no Step 03-05 output; perform no
#   integration, correction, final annotation, or cell/gene removal.
set -Eeuo pipefail

: "${PAPER3_PRIMARY_RUN_DIR:?Submit through submit_primary_processing_step_06.sh}"
: "${PAPER3_PRIMARY_WORKFLOW_ROOT:?Workflow root is required}"
: "${PAPER3_PRIMARY_PYTHON_BIN:?Python interpreter is required}"
RUN_DIR="${PAPER3_PRIMARY_RUN_DIR}"

# failed(line, status): ERR trap writes a run-local failure marker with job ID/time; no retry.
failed() {
  local line="$1" status="$2"
  printf 'COMPUTATION_FAILED\nstep=06_technical_sample_batch_diagnostics\njob_id=%s\nline=%s\nexit_status=%s\ntime=%s\n' "${SLURM_JOB_ID:-unknown}" "${line}" "${status}" "$(date --iso-8601=seconds)" > "${RUN_DIR}/COMPUTATION_FAILED.txt"
}
trap 'failed "${LINENO}" "$?"' ERR

# shellcheck disable=SC1090
# Frozen configuration -> explicit CLI options below -> dataclass settings. Repository edits have no effect.
source "${RUN_DIR}/config/resolved.env"
export PYTHONPATH="${RUN_DIR}/code"
export PYTHONHASHSEED=0
export LC_ALL=C
export TZ=UTC
# Execution controls only: thread limits follow the scheduler allocation, separate from scientific knobs.
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export NUMBA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MPLCONFIGDIR="${RUN_DIR}/provenance/matplotlib-cache"
mkdir -p "${MPLCONFIGDIR}"

{
  echo "started=$(date --iso-8601=seconds)"
  echo "host=$(hostname)"
  echo "slurm_job_id=${SLURM_JOB_ID:-unknown}"
  echo "slurm_cpus=${SLURM_CPUS_PER_TASK:-unknown}"
  echo "python=$("${PAPER3_PRIMARY_PYTHON_BIN}" --version 2>&1)"
} > "${RUN_DIR}/provenance/runtime_environment.txt"

# This invocation runs the analysis. --help on the module lists flags without starting it.
"${PAPER3_PRIMARY_PYTHON_BIN}" -m primary_processing.step06_cli \
  --input-h5ad "${INPUT_H5AD}" \
  --input-status "${RUN_DIR}/config/input_step02_STEP_STATUS.tsv" \
  --input-manifest "${RUN_DIR}/config/input_step02_output_manifest.tsv" \
  --approval-ledger-snapshot "${RUN_DIR}/config/approval_ledger_at_submission.tsv" \
  --bypass-decision "${RUN_DIR}/config/step04_step05_bypass_decision.md" \
  --run-dir "${RUN_DIR}" \
  --workflow-root "${WORKFLOW_ROOT}" \
  --frozen-code-dir "${RUN_DIR}/code/primary_processing" \
  --expected-step02-run-id "${PRIMARY_PROCESSING_STEP06_STEP02_RUN_ID}" \
  --expected-input-bytes "${PRIMARY_PROCESSING_STEP06_INPUT_BYTES}" \
  --expected-input-sha256 "${PRIMARY_PROCESSING_STEP06_INPUT_SHA256}" \
  --expected-cells "${PRIMARY_PROCESSING_STEP06_EXPECTED_CELLS}" \
  --expected-genes "${PRIMARY_PROCESSING_STEP06_EXPECTED_GENES}" \
  --sample-field "${PRIMARY_PROCESSING_STEP06_SAMPLE_FIELD}" \
  --genotype-field "${PRIMARY_PROCESSING_STEP06_GENOTYPE_FIELD}" \
  --sex-field "${PRIMARY_PROCESSING_STEP06_SEX_FIELD}" \
  --design-field "${PRIMARY_PROCESSING_STEP06_DESIGN_FIELD}" \
  --target-sum "${PRIMARY_PROCESSING_STEP06_TARGET_SUM}" \
  --hvg-flavor "${PRIMARY_PROCESSING_STEP06_HVG_FLAVOR}" \
  --n-top-genes "${PRIMARY_PROCESSING_STEP06_N_TOP_GENES}" \
  --pca-components "${PRIMARY_PROCESSING_STEP06_PCA_COMPONENTS}" \
  --centroid-components "${PRIMARY_PROCESSING_STEP06_CENTROID_COMPONENTS}" \
  --n-neighbors "${PRIMARY_PROCESSING_STEP06_N_NEIGHBORS}" \
  --n-jobs "${PRIMARY_PROCESSING_STEP06_N_JOBS}" \
  --neighbor-metric "${PRIMARY_PROCESSING_STEP06_NEIGHBOR_METRIC}" \
  --umap-min-dist "${PRIMARY_PROCESSING_STEP06_UMAP_MIN_DIST}" \
  --umap-spread "${PRIMARY_PROCESSING_STEP06_UMAP_SPREAD}" \
  --leiden-resolution "${PRIMARY_PROCESSING_STEP06_LEIDEN_RESOLUTION}" \
  --random-seed "${PRIMARY_PROCESSING_STEP06_RANDOM_SEED}" \
  --render-max-cells "${PRIMARY_PROCESSING_STEP06_RENDER_MAX_CELLS}" \
  --plot-dpi "${PRIMARY_PROCESSING_STEP06_PLOT_DPI}" \
  --compression "${PRIMARY_PROCESSING_STEP06_COMPRESSION}"

# This marker means computation finished; the scientific review state remains IN_REVIEW.
printf 'COMPUTATION_SUCCESS\nstep=06_technical_sample_batch_diagnostics\nstatus=IN_REVIEW\ninput=approved_step02\nsteps_03_05=bypassed\nintegration=none\ncells_removed=0\njob_id=%s\ncompleted=%s\n' "${SLURM_JOB_ID:-unknown}" "$(date --iso-8601=seconds)" > "${RUN_DIR}/COMPUTATION_SUCCESS.txt"
rm -f "${RUN_DIR}/COMPUTATION_FAILED.txt"
