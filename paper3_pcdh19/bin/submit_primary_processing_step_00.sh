#!/usr/bin/env bash
# PURPOSE
#   Freeze, submit, and track primary-processing Step 00 on Great Lakes.
#
# INPUTS
#   The original per-sample Cell Ranger matrices plus the registered biological
#   sample key and technical manifest. Scientific logic lives only in the
#   frozen object-oriented Python package.
#
# OUTPUTS
#   A versioned run below results/primary_processing/00_...; an exact submitted
#   scheduler script in PAPER3_ROOT/jobs; and logs below
#   PAPER3_ROOT/logs/primary_processing.
set -Eeuo pipefail

usage() {
  cat >&2 <<'EOF'
Usage: submit_primary_processing_step_00.sh [--dry-run] [--replace-run RUN_ID]

Default behavior creates a new versioned run. --replace-run intentionally
regenerates one existing inactive run inside Step 00 only. It cannot replace
an active job or escape the Step 00 result directory.
EOF
}

ORIGINAL_ARGS=("$@")
DRY_RUN=false
REPLACE_RUN_ID=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run)
      DRY_RUN=true
      shift
      ;;
    --replace-run)
      [[ $# -ge 2 ]] || { usage; exit 2; }
      REPLACE_RUN_ID="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 2
      ;;
  esac
done

BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
GREATLAKES_CONFIG="${BUNDLE_DIR}/config/greatlakes.env"
STEP_CONFIG="${BUNDLE_DIR}/config/primary_processing_step00.env"
REQUIREMENTS="${BUNDLE_DIR}/config/primary_processing_step00.requirements.txt"
SAMPLE_KEY="${BUNDLE_DIR}/config/sample_key.csv"
TECHNICAL_MANIFEST="${BUNDLE_DIR}/config/sample_manifest_draft.tsv"
PYTHON_PACKAGE="${BUNDLE_DIR}/scripts/primary_processing"
SBATCH_SOURCE="${BUNDLE_DIR}/slurm/primary_processing_00_input_validation_and_canonical_anndata.sbatch"
PACKAGE_README="${BUNDLE_DIR}/templates/PRIMARY_PROCESSING_STEP00_OUTPUT_PACKAGE_README.md"
HANDOFF="${BUNDLE_DIR}/PCDH19_PRIMARY_PROCESSING_HANDOFF.md"

# shellcheck disable=SC1090
source "${GREATLAKES_CONFIG}"
# shellcheck disable=SC1090
source "${STEP_CONFIG}"
: "${REPO_ROOT:?REPO_ROOT is required}"
: "${PAPER3_ROOT:?PAPER3_ROOT is required}"
: "${PAPER3_CELLRANGER_ROOT:?PAPER3_CELLRANGER_ROOT is required}"
: "${ACCOUNT:?ACCOUNT is required}"
: "${PARTITION:?PARTITION is required}"
: "${PRIMARY_PROCESSING_PYTHON_BIN:?PRIMARY_PROCESSING_PYTHON_BIN is required}"

for required in \
  "${GREATLAKES_CONFIG}" "${STEP_CONFIG}" "${REQUIREMENTS}" \
  "${SAMPLE_KEY}" "${TECHNICAL_MANIFEST}" "${SBATCH_SOURCE}" \
  "${PACKAGE_README}" "${HANDOFF}" "${PRIMARY_PROCESSING_PYTHON_BIN}"; do
  [[ -f "${required}" ]] || { echo "Missing required file: ${required}" >&2; exit 2; }
done
for module in __init__.py cli.py loaders.py models.py publishing.py validation.py workflow.py; do
  [[ -f "${PYTHON_PACKAGE}/${module}" ]] || { echo "Missing Python module: ${module}" >&2; exit 2; }
done

STEP="00_input_validation_and_canonical_anndata"
WORKFLOW_ROOT="${PAPER3_ROOT}/results/primary_processing"
STEP_ROOT="${WORKFLOW_ROOT}/${STEP}"
LOG_ROOT="${PAPER3_ROOT}/logs/primary_processing"
STAMP="$(date +%Y%m%d_%H%M%S)"
GIT_SHORT="$(git -C "${REPO_ROOT}" rev-parse --short HEAD)"
OUTPUT_MODE=versioned
RUN_ID="${STEP}_${STAMP}_${GIT_SHORT}"
if [[ -n "${REPLACE_RUN_ID}" ]]; then
  OUTPUT_MODE=replace
  RUN_ID="${REPLACE_RUN_ID}"
  if [[ "${RUN_ID}" == *[!A-Za-z0-9._-]* || "${RUN_ID}" != "${STEP}"_* ]]; then
    echo "Unsafe replacement run ID: ${RUN_ID}" >&2
    exit 2
  fi
fi
RUN_DIR="${STEP_ROOT}/${RUN_ID}"
if [[ "$(dirname "$(realpath -m "${RUN_DIR}")")" != "$(realpath -m "${STEP_ROOT}")" ]]; then
  echo "Resolved run directory escaped the Step 00 root: ${RUN_DIR}" >&2
  exit 2
fi

if [[ "${OUTPUT_MODE}" == replace && -f "${RUN_DIR}/provenance/job_ids.tsv" ]]; then
  while IFS= read -r prior_job; do
    [[ -n "${prior_job}" ]] || continue
    prior_state="$(squeue -h -j "${prior_job}" -o '%T' | head -n 1)"
    if [[ -n "${prior_state}" ]]; then
      echo "Refusing to replace ${RUN_ID}; job ${prior_job} is ${prior_state}." >&2
      exit 2
    fi
  done < <(awk -F '\t' 'NR > 1 && $2 != "" {print $2}' "${RUN_DIR}/provenance/job_ids.tsv")
fi

echo "Step: ${STEP}"
echo "Output mode: ${OUTPUT_MODE}"
echo "Run directory: ${RUN_DIR}"
echo "Cell Ranger root: ${PAPER3_CELLRANGER_ROOT}"
echo "Python: ${PRIMARY_PROCESSING_PYTHON_BIN}"
echo "Resources: ${PRIMARY_PROCESSING_STEP00_CPUS} CPU, ${PRIMARY_PROCESSING_STEP00_MEMORY}, ${PRIMARY_PROCESSING_STEP00_WALLTIME}"

"${PRIMARY_PROCESSING_PYTHON_BIN}" - <<'PY'
from importlib import metadata
required = {
    "anndata": "0.12.14",
    "h5py": "3.16.0",
    "numpy": "2.4.4",
    "pandas": "2.3.3",
    "scanpy": "1.11.5",
    "scipy": "1.17.1",
}
observed = {name: metadata.version(name) for name in required}
if observed != required:
    raise SystemExit(f"Primary-processing environment mismatch: {observed!r}")
PY
PYTHONPATH="${BUNDLE_DIR}/scripts" "${PRIMARY_PROCESSING_PYTHON_BIN}" -m compileall -q "${PYTHON_PACKAGE}"
PYTHONPATH="${BUNDLE_DIR}/scripts" "${PRIMARY_PROCESSING_PYTHON_BIN}" -m primary_processing.cli --help >/dev/null
for sample_number in $(seq 1 "${PRIMARY_PROCESSING_EXPECTED_SAMPLES}"); do
  sample_dir="${PAPER3_CELLRANGER_ROOT}/per_sample_outs/15662-JZ-${sample_number}"
  for filename in sample_filtered_feature_bc_matrix.h5 sample_raw_feature_bc_matrix.h5 sample_filtered_barcodes.csv metrics_summary.csv; do
    [[ -r "${sample_dir}/${filename}" ]] || { echo "Unreadable input: ${sample_dir}/${filename}" >&2; exit 2; }
  done
done

if [[ "${DRY_RUN}" == true ]]; then
  echo "Dry run passed; no directory was created and no job was submitted."
  exit 0
fi

"${BUNDLE_DIR}/bin/initialize_turbo.sh" >/dev/null
mkdir -p "${STEP_ROOT}" "${LOG_ROOT}"
JOB_FILE="${PAPER3_ROOT}/jobs/${RUN_ID}.sbatch"
if [[ "${OUTPUT_MODE}" == versioned ]]; then
  for target in "${RUN_DIR}" "${JOB_FILE}"; do
    [[ ! -e "${target}" ]] || { echo "Refusing existing versioned target: ${target}" >&2; exit 2; }
  done
elif [[ -d "${RUN_DIR}" ]]; then
  echo "Intentionally replacing inactive Step 00 run: ${RUN_DIR}"
  find "${RUN_DIR}" -mindepth 1 -depth -delete
elif [[ -e "${RUN_DIR}" ]]; then
  echo "Replacement target is not a directory: ${RUN_DIR}" >&2
  exit 2
fi

mkdir -p "${RUN_DIR}/code/primary_processing" "${RUN_DIR}/config" "${RUN_DIR}/logs" "${RUN_DIR}/provenance"
cp -p "${PYTHON_PACKAGE}"/*.py "${RUN_DIR}/code/primary_processing/"
cp -p "${BASH_SOURCE[0]}" "${SBATCH_SOURCE}" "${RUN_DIR}/code/"
cp -p "${GREATLAKES_CONFIG}" "${RUN_DIR}/config/submitted_greatlakes.env"
cp -p "${STEP_CONFIG}" "${RUN_DIR}/config/submitted_step00.env"
cp -p "${REQUIREMENTS}" "${RUN_DIR}/config/requirements.txt"
cp -p "${SAMPLE_KEY}" "${RUN_DIR}/config/sample_key.csv"
cp -p "${TECHNICAL_MANIFEST}" "${RUN_DIR}/config/sample_manifest.tsv"
cp -p "${HANDOFF}" "${RUN_DIR}/config/authoritative_handoff.md"
cp -p "${PACKAGE_README}" "${RUN_DIR}/README.md"
cp -p "${RUN_DIR}/code/$(basename "${SBATCH_SOURCE}")" "${JOB_FILE}"

{
  echo "REPO_ROOT=${REPO_ROOT}"
  echo "PAPER3_ROOT=${PAPER3_ROOT}"
  echo "PAPER3_CELLRANGER_ROOT=${PAPER3_CELLRANGER_ROOT}"
  echo "STEP=${STEP}"
  echo "RUN_ID=${RUN_ID}"
  echo "RUN_DIR=${RUN_DIR}"
  echo "WORKFLOW_ROOT=${WORKFLOW_ROOT}"
  echo "OUTPUT_MODE=${OUTPUT_MODE}"
  echo "ACCOUNT=${ACCOUNT}"
  echo "PARTITION=${PARTITION}"
  echo "PRIMARY_PROCESSING_PYTHON_BIN=${PRIMARY_PROCESSING_PYTHON_BIN}"
  echo "PRIMARY_PROCESSING_STEP00_COMPRESSION=${PRIMARY_PROCESSING_STEP00_COMPRESSION}"
  echo "PRIMARY_PROCESSING_EXPECTED_SAMPLES=${PRIMARY_PROCESSING_EXPECTED_SAMPLES}"
  echo "PRIMARY_PROCESSING_EXPECTED_CELLS=${PRIMARY_PROCESSING_EXPECTED_CELLS}"
  echo "PRIMARY_PROCESSING_EXPECTED_GENES=${PRIMARY_PROCESSING_EXPECTED_GENES}"
  echo "PRIMARY_PROCESSING_EXPECTED_GENOME=${PRIMARY_PROCESSING_EXPECTED_GENOME}"
  printf 'PRIMARY_PROCESSING_FEATURE_TYPE=%q\n' "${PRIMARY_PROCESSING_FEATURE_TYPE}"
} > "${RUN_DIR}/config/resolved.env"

{
  echo "submitted_at=$(date --iso-8601=seconds)"
  echo "submit_host=$(hostname)"
  echo "submit_user=${USER:-unknown}"
  printf 'command='
  printf '%q ' "$0" "${ORIGINAL_ARGS[@]}"
  printf '\n'
} > "${RUN_DIR}/provenance/submission.txt"
{
  echo "git_commit=$(git -C "${REPO_ROOT}" rev-parse HEAD)"
  echo "git_short=${GIT_SHORT}"
  echo "git_status_begin"
  git -C "${REPO_ROOT}" status --short
  echo "git_status_end"
} > "${RUN_DIR}/provenance/repository_state.txt"
"${PRIMARY_PROCESSING_PYTHON_BIN}" -m pip freeze > "${RUN_DIR}/config/python_pip_freeze.txt"

EXPORTS="ALL,PAPER3_PRIMARY_RUN_DIR=${RUN_DIR},PAPER3_PRIMARY_WORKFLOW_ROOT=${WORKFLOW_ROOT},PAPER3_PRIMARY_CELLRANGER_ROOT=${PAPER3_CELLRANGER_ROOT},PAPER3_PRIMARY_PYTHON_BIN=${PRIMARY_PROCESSING_PYTHON_BIN}"
JOB_ID="$(sbatch --parsable \
  --job-name=pcdh19-primary-00 --account="${ACCOUNT}" --partition="${PARTITION}" \
  --nodes=1 --ntasks=1 --cpus-per-task="${PRIMARY_PROCESSING_STEP00_CPUS}" \
  --mem="${PRIMARY_PROCESSING_STEP00_MEMORY}" --time="${PRIMARY_PROCESSING_STEP00_WALLTIME}" \
  --output="${LOG_ROOT}/${RUN_ID}-%j.out" --error="${LOG_ROOT}/${RUN_ID}-%j.err" \
  --export="${EXPORTS}" "${JOB_FILE}")"
printf 'stage\tjob_id\nstep00\t%s\n' "${JOB_ID}" > "${RUN_DIR}/provenance/job_ids.tsv"
ln -s "${LOG_ROOT}/${RUN_ID}-${JOB_ID}.out" "${RUN_DIR}/logs/scheduler.out"
ln -s "${LOG_ROOT}/${RUN_ID}-${JOB_ID}.err" "${RUN_DIR}/logs/scheduler.err"

echo "Submitted Step 00 job ${JOB_ID}"
echo "Run ID: ${RUN_ID}"
echo "Run directory: ${RUN_DIR}"
echo "Status after successful computation will be IN_REVIEW."
