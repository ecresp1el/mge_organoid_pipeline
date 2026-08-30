#!/usr/bin/env bash
# PURPOSE
#   Freeze and submit Scanpy QC-metric Step 01 without filtering.
#
# INPUT
#   The exact user-approved Step 00 canonical raw-count H5AD.
#
# OUTPUT
#   A versioned 01_qc_metrics run, frozen scheduler job, scheduler logs,
#   unchanged raw-count checkpoint with QC metadata, review tables, and plots.
set -Eeuo pipefail

usage() {
  cat >&2 <<'EOF'
Usage: submit_primary_processing_step_01.sh [--dry-run] [--replace-run RUN_ID]

Default behavior creates a new versioned run. --replace-run regenerates only
one named inactive 01_qc_metrics run and refuses active jobs or unsafe paths.
EOF
}

ORIGINAL_ARGS=("$@")
DRY_RUN=false
REPLACE_RUN_ID=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run) DRY_RUN=true; shift ;;
    --replace-run) [[ $# -ge 2 ]] || { usage; exit 2; }; REPLACE_RUN_ID="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done

BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
GREATLAKES_CONFIG="${BUNDLE_DIR}/config/greatlakes.env"
STEP_CONFIG="${BUNDLE_DIR}/config/primary_processing_step01.env"
REQUIREMENTS="${BUNDLE_DIR}/config/primary_processing_step01.requirements.txt"
PYTHON_PACKAGE="${BUNDLE_DIR}/scripts/primary_processing"
SBATCH_SOURCE="${BUNDLE_DIR}/slurm/primary_processing_01_qc_metrics.sbatch"
PACKAGE_README="${BUNDLE_DIR}/templates/PRIMARY_PROCESSING_STEP01_OUTPUT_PACKAGE_README.md"
HANDOFF="${BUNDLE_DIR}/PCDH19_PRIMARY_PROCESSING_HANDOFF.md"

# shellcheck disable=SC1090
source "${GREATLAKES_CONFIG}"
# shellcheck disable=SC1090
source "${STEP_CONFIG}"
: "${REPO_ROOT:?}" "${PAPER3_ROOT:?}" "${ACCOUNT:?}" "${PARTITION:?}"
: "${PRIMARY_PROCESSING_PYTHON_BIN:?Set by primary_processing_step00.env or environment}"

STEP00_RUN_DIR="${PAPER3_ROOT}/results/primary_processing/00_input_validation_and_canonical_anndata/${PRIMARY_PROCESSING_STEP01_INPUT_RUN_ID}"
INPUT_H5AD="${STEP00_RUN_DIR}/${PRIMARY_PROCESSING_STEP01_INPUT_CHECKPOINT}"
INPUT_STATUS="${STEP00_RUN_DIR}/STEP_STATUS.tsv"
INPUT_MANIFEST="${STEP00_RUN_DIR}/tables/output_manifest.tsv"
APPROVAL_LEDGER="${PAPER3_ROOT}/results/primary_processing/APPROVAL_LEDGER.tsv"
for required in "${GREATLAKES_CONFIG}" "${STEP_CONFIG}" "${REQUIREMENTS}" "${SBATCH_SOURCE}" "${PACKAGE_README}" "${HANDOFF}" "${INPUT_H5AD}" "${INPUT_STATUS}" "${INPUT_MANIFEST}" "${APPROVAL_LEDGER}" "${PRIMARY_PROCESSING_PYTHON_BIN}"; do
  [[ -f "${required}" ]] || { echo "Missing required Step 01 asset: ${required}" >&2; exit 2; }
done
for module in step01_cli.py step01_metrics.py step01_models.py step01_plots.py step01_publishing.py step01_validation.py step01_workflow.py; do
  [[ -f "${PYTHON_PACKAGE}/${module}" ]] || { echo "Missing Step 01 module: ${module}" >&2; exit 2; }
done

awk -F '\t' -v run="${PRIMARY_PROCESSING_STEP01_INPUT_RUN_ID}" 'NR==1 {for(i=1;i<=NF;i++) h[$i]=i; next} $(h["run_id"])==run && $(h["status"])=="APPROVED" && $(h["approved_run_id"])==run {ok=1} END {exit(ok?0:1)}' "${APPROVAL_LEDGER}" || { echo "Exact Step 00 run is not APPROVED in the ledger" >&2; exit 2; }
[[ "$(stat -c '%s' "${INPUT_H5AD}")" == "${PRIMARY_PROCESSING_STEP01_INPUT_BYTES}" ]] || { echo "Approved input size mismatch" >&2; exit 2; }

"${PRIMARY_PROCESSING_PYTHON_BIN}" - <<'PY'
from importlib import metadata
required = {"anndata":"0.12.14","h5py":"3.16.0","matplotlib":"3.10.9","numpy":"2.4.4","pandas":"2.3.3","scanpy":"1.11.5","scipy":"1.17.1"}
observed = {name: metadata.version(name) for name in required}
if observed != required:
    raise SystemExit(f"Primary-processing Step 01 environment mismatch: {observed!r}")
PY
PYTHONPATH="${BUNDLE_DIR}/scripts" "${PRIMARY_PROCESSING_PYTHON_BIN}" -m compileall -q "${PYTHON_PACKAGE}"
PYTHONPATH="${BUNDLE_DIR}/scripts" "${PRIMARY_PROCESSING_PYTHON_BIN}" -m primary_processing.step01_cli --help >/dev/null

STEP="01_qc_metrics"
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
[[ "$(dirname "$(realpath -m "${RUN_DIR}")")" == "$(realpath -m "${STEP_ROOT}")" ]] || { echo "Run escaped Step 01 root" >&2; exit 2; }

PREVIOUS_JOB_IDS=""
if [[ "${OUTPUT_MODE}" == replace && -f "${RUN_DIR}/provenance/job_ids.tsv" ]]; then
  PREVIOUS_JOB_IDS="$(awk -F '\t' 'NR>1 && $2!="" {print $2}' "${RUN_DIR}/provenance/job_ids.tsv" | paste -sd, -)"
  while IFS= read -r prior_job; do
    [[ -n "${prior_job}" ]] || continue
    prior_state="$(squeue -h -j "${prior_job}" -o '%T' | head -n 1)"
    [[ -z "${prior_state}" ]] || { echo "Refusing active replacement: ${prior_job} ${prior_state}" >&2; exit 2; }
  done < <(tr ',' '\n' <<< "${PREVIOUS_JOB_IDS}")
fi

echo "Step: ${STEP}"
echo "Approved input: ${PRIMARY_PROCESSING_STEP01_INPUT_RUN_ID}"
echo "Output mode: ${OUTPUT_MODE}"
echo "Run directory: ${RUN_DIR}"
echo "Resources: ${PRIMARY_PROCESSING_STEP01_CPUS} CPU, ${PRIMARY_PROCESSING_STEP01_MEMORY}, ${PRIMARY_PROCESSING_STEP01_WALLTIME}"
echo "QC: Scanpy calculate_qc_metrics; percent_top=${PRIMARY_PROCESSING_STEP01_PERCENT_TOP}; no filtering"
if [[ "${DRY_RUN}" == true ]]; then
  echo "Dry run passed; no directory was created and no job was submitted."
  exit 0
fi

"${BUNDLE_DIR}/bin/initialize_turbo.sh" >/dev/null
mkdir -p "${STEP_ROOT}" "${LOG_ROOT}"
JOB_FILE="${PAPER3_ROOT}/jobs/${RUN_ID}.sbatch"
if [[ "${OUTPUT_MODE}" == versioned ]]; then
  [[ ! -e "${RUN_DIR}" && ! -e "${JOB_FILE}" ]] || { echo "Refusing existing versioned target" >&2; exit 2; }
elif [[ -d "${RUN_DIR}" ]]; then
  echo "Intentionally replacing inactive Step 01 run: ${RUN_DIR}"
  find "${RUN_DIR}" -mindepth 1 -depth -delete
elif [[ -e "${RUN_DIR}" ]]; then
  echo "Replacement target is not a directory" >&2
  exit 2
fi

mkdir -p "${RUN_DIR}/code/primary_processing" "${RUN_DIR}/config" "${RUN_DIR}/logs" "${RUN_DIR}/provenance"
cp -p "${PYTHON_PACKAGE}"/*.py "${RUN_DIR}/code/primary_processing/"
cp -p "${BASH_SOURCE[0]}" "${SBATCH_SOURCE}" "${RUN_DIR}/code/"
cp -p "${GREATLAKES_CONFIG}" "${RUN_DIR}/config/submitted_greatlakes.env"
cp -p "${STEP_CONFIG}" "${RUN_DIR}/config/submitted_step01.env"
cp -p "${REQUIREMENTS}" "${RUN_DIR}/config/requirements.txt"
cp -p "${HANDOFF}" "${RUN_DIR}/config/authoritative_handoff.md"
cp -p "${INPUT_STATUS}" "${RUN_DIR}/config/approved_step00_STEP_STATUS.tsv"
cp -p "${INPUT_MANIFEST}" "${RUN_DIR}/config/approved_step00_output_manifest.tsv"
cp -p "${APPROVAL_LEDGER}" "${RUN_DIR}/config/approval_ledger_at_submission.tsv"
cp -p "${PACKAGE_README}" "${RUN_DIR}/README.md"
cp -p "${RUN_DIR}/code/$(basename "${SBATCH_SOURCE}")" "${JOB_FILE}"

{
  printf 'REPO_ROOT=%q\n' "${REPO_ROOT}"
  printf 'PAPER3_ROOT=%q\n' "${PAPER3_ROOT}"
  printf 'INPUT_H5AD=%q\n' "${INPUT_H5AD}"
  echo "STEP=${STEP}"
  echo "RUN_ID=${RUN_ID}"
  printf 'RUN_DIR=%q\n' "${RUN_DIR}"
  printf 'WORKFLOW_ROOT=%q\n' "${WORKFLOW_ROOT}"
  echo "OUTPUT_MODE=${OUTPUT_MODE}"
  echo "REPLACED_PREVIOUS_JOB_IDS=${PREVIOUS_JOB_IDS}"
  echo "PRIMARY_PROCESSING_PYTHON_BIN=${PRIMARY_PROCESSING_PYTHON_BIN}"
  grep '^PRIMARY_PROCESSING_STEP01_' "${STEP_CONFIG}"
} > "${RUN_DIR}/config/resolved.env"
{
  echo "submitted_at=$(date --iso-8601=seconds)"
  echo "submit_host=$(hostname)"
  echo "submit_user=${USER:-unknown}"
  printf 'command='; printf '%q ' "$0" "${ORIGINAL_ARGS[@]}"; printf '\n'
} > "${RUN_DIR}/provenance/submission.txt"
{
  echo "git_commit=$(git -C "${REPO_ROOT}" rev-parse HEAD)"
  echo "git_short=${GIT_SHORT}"
  echo "git_status_begin"; git -C "${REPO_ROOT}" status --short; echo "git_status_end"
} > "${RUN_DIR}/provenance/repository_state.txt"
"${PRIMARY_PROCESSING_PYTHON_BIN}" -m pip freeze > "${RUN_DIR}/config/python_pip_freeze.txt"
if [[ -n "${PREVIOUS_JOB_IDS}" ]]; then
  { echo previous_job_id; tr ',' '\n' <<< "${PREVIOUS_JOB_IDS}"; } > "${RUN_DIR}/provenance/replaced_previous_job_ids.txt"
fi

EXPORTS="ALL,PAPER3_PRIMARY_RUN_DIR=${RUN_DIR},PAPER3_PRIMARY_WORKFLOW_ROOT=${WORKFLOW_ROOT},PAPER3_PRIMARY_PYTHON_BIN=${PRIMARY_PROCESSING_PYTHON_BIN}"
JOB_ID="$(sbatch --parsable --job-name=pcdh19-primary-01 --account="${ACCOUNT}" --partition="${PARTITION}" --nodes=1 --ntasks=1 --cpus-per-task="${PRIMARY_PROCESSING_STEP01_CPUS}" --mem="${PRIMARY_PROCESSING_STEP01_MEMORY}" --time="${PRIMARY_PROCESSING_STEP01_WALLTIME}" --output="${LOG_ROOT}/${RUN_ID}-%j.out" --error="${LOG_ROOT}/${RUN_ID}-%j.err" --export="${EXPORTS}" "${JOB_FILE}")"
printf 'stage\tjob_id\nstep01\t%s\n' "${JOB_ID}" > "${RUN_DIR}/provenance/job_ids.tsv"
ln -s "${LOG_ROOT}/${RUN_ID}-${JOB_ID}.out" "${RUN_DIR}/logs/scheduler.out"
ln -s "${LOG_ROOT}/${RUN_ID}-${JOB_ID}.err" "${RUN_DIR}/logs/scheduler.err"
echo "Submitted Step 01 job ${JOB_ID}"
echo "Run ID: ${RUN_ID}"
echo "Run directory: ${RUN_DIR}"
echo "Successful computation will remain IN_REVIEW."
