#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat >&2 <<'EOF'
Usage: submit_developing_mouse_mge_reference_curation.sh [--dry-run] [--replace-run RUN_ID]

Default behavior creates a new versioned run. --replace-run intentionally
clears and regenerates one existing inactive run inside this exact step only; active
SLURM jobs are never replaced.
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
CURATION_CONFIG="${BUNDLE_DIR}/config/developing_mouse_mge_reference_curation.env"
REGISTRY="${BUNDLE_DIR}/config/developing_mouse_mge_source_registry.tsv"
REQUIREMENTS="${BUNDLE_DIR}/config/developing_mouse_mge_reference_curation.requirements.txt"
PY_SCRIPT="${BUNDLE_DIR}/scripts/reference_curation/developing_mouse_mge_reference_curation.py"
R_SCRIPT="${BUNDLE_DIR}/scripts/reference_curation/inspect_bandler_rds.R"
REPORT_SCRIPT="${BUNDLE_DIR}/scripts/reference_curation/generate_reference_curation_report.py"
SOURCE_SBATCH="${BUNDLE_DIR}/slurm/00a_developing_mouse_mge_source_audit.sbatch"
P0_SBATCH="${BUNDLE_DIR}/slurm/00b_developing_mouse_mge_p0_inspection.sbatch"
CHECKPOINT_SBATCH="${BUNDLE_DIR}/slurm/00c_developing_mouse_mge_checkpoint.sbatch"
PACKAGE_README="${BUNDLE_DIR}/templates/REFERENCE_CURATION_OUTPUT_PACKAGE_README.md"
HANDOFF="${BUNDLE_DIR}/PCDH19_DEVELOPING_MOUSE_MGE_REFERENCE_CURATION_HANDOFF.md"

# shellcheck disable=SC1090
source "${GREATLAKES_CONFIG}"
# shellcheck disable=SC1090
source "${CURATION_CONFIG}"
: "${REPO_ROOT:?REPO_ROOT is required}"
: "${PAPER3_ROOT:?PAPER3_ROOT is required}"
: "${ACCOUNT:?ACCOUNT is required}"
: "${PARTITION:?PARTITION is required}"
: "${CURATION_PYTHON_BIN:?CURATION_PYTHON_BIN is required}"
: "${CURATION_R_MODULE:?CURATION_R_MODULE is required}"

for required in \
  "${GREATLAKES_CONFIG}" "${CURATION_CONFIG}" "${REGISTRY}" "${REQUIREMENTS}" \
  "${PY_SCRIPT}" "${R_SCRIPT}" "${REPORT_SCRIPT}" "${SOURCE_SBATCH}" "${P0_SBATCH}" \
  "${CHECKPOINT_SBATCH}" "${PACKAGE_README}" "${HANDOFF}" "${CURATION_PYTHON_BIN}"; do
  [[ -f "${required}" ]] || { echo "Missing required file: ${required}" >&2; exit 2; }
done

STEP="00_developing_mouse_mge_reference_curation"
STEP_ROOT="${PAPER3_ROOT}/results/${STEP}"
SOURCE_ROOT="${PAPER3_ROOT}/inputs/developing_mouse_mge"
STAMP="$(date +%Y%m%d_%H%M%S)"
GIT_SHORT="$(git -C "${REPO_ROOT}" rev-parse --short HEAD)"
OUTPUT_MODE=versioned
RUN_ID="${STEP}_${STAMP}_${GIT_SHORT}"
if [[ -n "${REPLACE_RUN_ID}" ]]; then
  OUTPUT_MODE=replace
  RUN_ID="${REPLACE_RUN_ID}"
  if [[ "${RUN_ID}" == *[!A-Za-z0-9._-]* || "${RUN_ID}" != "${STEP}"_* ]]; then
    echo "Unsafe replacement run ID: ${RUN_ID}" >&2
    echo "It must be a basename beginning with ${STEP}_ and contain only letters, numbers, dot, underscore, or hyphen." >&2
    exit 2
  fi
fi
RUN_DIR="${STEP_ROOT}/${RUN_ID}"
EXPECTED_PARENT="$(realpath -m "${STEP_ROOT}")"
if [[ "$(dirname "$(realpath -m "${RUN_DIR}")")" != "${EXPECTED_PARENT}" ]]; then
  echo "Resolved run directory escaped the permitted step root: ${RUN_DIR}" >&2
  exit 2
fi

SOURCE_JOB_FILE="${PAPER3_ROOT}/jobs/${RUN_ID}.source_audit.sbatch"
P0_JOB_FILE="${PAPER3_ROOT}/jobs/${RUN_ID}.p0_inspection.sbatch"
CHECKPOINT_JOB_FILE="${PAPER3_ROOT}/jobs/${RUN_ID}.checkpoint.sbatch"
PREVIOUS_JOB_IDS=""
if [[ "${OUTPUT_MODE}" == replace && -f "${RUN_DIR}/provenance/job_ids.tsv" ]]; then
  PREVIOUS_JOB_IDS="$(awk -F '\t' 'NR > 1 && $2 != "" {print $2}' "${RUN_DIR}/provenance/job_ids.tsv" | paste -sd, -)"
  while IFS= read -r prior_job; do
    [[ -n "${prior_job}" ]] || continue
    prior_state="$(squeue -h -j "${prior_job}" -o '%T' | head -n 1)"
    if [[ -n "${prior_state}" ]]; then
      echo "Refusing to replace ${RUN_ID}; job ${prior_job} is still ${prior_state}." >&2
      exit 2
    fi
  done < <(tr ',' '\n' <<< "${PREVIOUS_JOB_IDS}")
fi

echo "Step: ${STEP}"
echo "Output mode: ${OUTPUT_MODE}"
echo "Prospective run directory: ${RUN_DIR}"
echo "Reusable P0 source cache: ${SOURCE_ROOT}"
echo "Python: ${CURATION_PYTHON_BIN}"
echo "R module: ${CURATION_R_MODULE}"
echo "Source audit: ${CURATION_SOURCE_AUDIT_CPUS} CPU, ${CURATION_SOURCE_AUDIT_MEMORY}, ${CURATION_SOURCE_AUDIT_WALLTIME}"
echo "P0 array: 3 tasks, max ${CURATION_P0_ARRAY_CONCURRENCY} active; ${CURATION_P0_CPUS} CPU, ${CURATION_P0_MEMORY}, ${CURATION_P0_WALLTIME}"
echo "Checkpoint: ${CURATION_CHECKPOINT_CPUS} CPU, ${CURATION_CHECKPOINT_MEMORY}, ${CURATION_CHECKPOINT_WALLTIME}"
if [[ -n "${PREVIOUS_JOB_IDS}" ]]; then
  echo "Completed prior jobs to replace: ${PREVIOUS_JOB_IDS}"
fi

"${CURATION_PYTHON_BIN}" "${PY_SCRIPT}" --help >/dev/null
if [[ "${DRY_RUN}" == true ]]; then
  echo "Dry run passed; nothing was created, downloaded, deleted, or submitted."
  exit 0
fi

"${BUNDLE_DIR}/bin/initialize_turbo.sh" >/dev/null
if [[ "${OUTPUT_MODE}" == versioned ]]; then
  for target in "${RUN_DIR}" "${SOURCE_JOB_FILE}" "${P0_JOB_FILE}" "${CHECKPOINT_JOB_FILE}"; do
    [[ ! -e "${target}" ]] || { echo "Refusing to reuse existing versioned target: ${target}" >&2; exit 2; }
  done
elif [[ -d "${RUN_DIR}" ]]; then
  echo "Intentionally replacing inactive run within ${STEP}: ${RUN_DIR}"
  find "${RUN_DIR}" -mindepth 1 -depth -delete
elif [[ -e "${RUN_DIR}" ]]; then
  echo "Replacement target is not a directory: ${RUN_DIR}" >&2
  exit 2
fi

mkdir -p \
  "${RUN_DIR}/code" "${RUN_DIR}/config" "${RUN_DIR}/tables" "${RUN_DIR}/logs" "${RUN_DIR}/provenance" \
  "${RUN_DIR}/LaManno2021/metadata" "${RUN_DIR}/LaManno2021/figures" "${RUN_DIR}/LaManno2021/audit" \
  "${RUN_DIR}/Bandler2022/metadata" "${RUN_DIR}/Bandler2022/figures" "${RUN_DIR}/Bandler2022/audit" \
  "${RUN_DIR}/Mayer2018/metadata" "${RUN_DIR}/Mayer2018/figures" "${RUN_DIR}/Mayer2018/audit" \
  "${SOURCE_ROOT}"

cp -p "${PY_SCRIPT}" "${R_SCRIPT}" "${REPORT_SCRIPT}" "${SOURCE_SBATCH}" "${P0_SBATCH}" "${CHECKPOINT_SBATCH}" \
  "${BASH_SOURCE[0]}" "${RUN_DIR}/code/"
cp -p "${GREATLAKES_CONFIG}" "${RUN_DIR}/config/submitted_greatlakes.env"
cp -p "${CURATION_CONFIG}" "${RUN_DIR}/config/submitted_curation.env"
cp -p "${REGISTRY}" "${RUN_DIR}/config/source_registry.tsv"
cp -p "${REQUIREMENTS}" "${RUN_DIR}/config/requirements.txt"
cp -p "${HANDOFF}" "${RUN_DIR}/config/authoritative_handoff.md"
cp -p "${PACKAGE_README}" "${RUN_DIR}/README.md"
cp -p "${RUN_DIR}/code/00a_developing_mouse_mge_source_audit.sbatch" "${SOURCE_JOB_FILE}"
cp -p "${RUN_DIR}/code/00b_developing_mouse_mge_p0_inspection.sbatch" "${P0_JOB_FILE}"
cp -p "${RUN_DIR}/code/00c_developing_mouse_mge_checkpoint.sbatch" "${CHECKPOINT_JOB_FILE}"

{
  echo "REPO_ROOT=${REPO_ROOT}"
  echo "PAPER3_ROOT=${PAPER3_ROOT}"
  echo "STEP=${STEP}"
  echo "RUN_ID=${RUN_ID}"
  echo "RUN_DIR=${RUN_DIR}"
  echo "SOURCE_ROOT=${SOURCE_ROOT}"
  echo "OUTPUT_MODE=${OUTPUT_MODE}"
  echo "REPLACED_PREVIOUS_JOB_IDS=${PREVIOUS_JOB_IDS}"
  echo "ACCOUNT=${ACCOUNT}"
  echo "PARTITION=${PARTITION}"
  echo "CURATION_PYTHON_BIN=${CURATION_PYTHON_BIN}"
  echo "CURATION_R_MODULE=${CURATION_R_MODULE}"
} > "${RUN_DIR}/config/resolved.env"

{
  echo "submitted_at=$(date --iso-8601=seconds)"
  echo "submit_host=$(hostname)"
  echo "submit_user=${USER:-unknown}"
  printf 'command='
  printf '%q ' "$0" "${ORIGINAL_ARGS[@]}"
  printf '\n'
} > "${RUN_DIR}/provenance/submission.txt"

EXPORTS="ALL,PAPER3_CURATION_RUN_DIR=${RUN_DIR},PAPER3_CURATION_SOURCE_ROOT=${SOURCE_ROOT},PAPER3_CURATION_REGISTRY=${RUN_DIR}/config/source_registry.tsv,PAPER3_CURATION_PYTHON_BIN=${CURATION_PYTHON_BIN},PAPER3_CURATION_R_MODULE=${CURATION_R_MODULE}"
printf 'stage\tjob_id\n' > "${RUN_DIR}/provenance/job_ids.tsv"
SOURCE_JOB_ID="$(sbatch --parsable \
  --job-name=mge-ref-source --account="${ACCOUNT}" --partition="${PARTITION}" \
  --nodes=1 --ntasks=1 --cpus-per-task="${CURATION_SOURCE_AUDIT_CPUS}" \
  --mem="${CURATION_SOURCE_AUDIT_MEMORY}" --time="${CURATION_SOURCE_AUDIT_WALLTIME}" \
  --output="${RUN_DIR}/logs/source-audit-%j.out" --error="${RUN_DIR}/logs/source-audit-%j.err" \
  --export="${EXPORTS}" "${SOURCE_JOB_FILE}")"
printf 'source_audit\t%s\n' "${SOURCE_JOB_ID}" >> "${RUN_DIR}/provenance/job_ids.tsv"

P0_JOB_ID="$(sbatch --parsable \
  --job-name=mge-ref-p0 --dependency="afterok:${SOURCE_JOB_ID}" \
  --account="${ACCOUNT}" --partition="${PARTITION}" --nodes=1 --ntasks=1 \
  --cpus-per-task="${CURATION_P0_CPUS}" --mem="${CURATION_P0_MEMORY}" --time="${CURATION_P0_WALLTIME}" \
  --array="0-2%${CURATION_P0_ARRAY_CONCURRENCY}" \
  --output="${RUN_DIR}/logs/p0-%A_%a.out" --error="${RUN_DIR}/logs/p0-%A_%a.err" \
  --export="${EXPORTS}" "${P0_JOB_FILE}")"
printf 'p0_inspection_array\t%s\n' "${P0_JOB_ID}" >> "${RUN_DIR}/provenance/job_ids.tsv"

CHECKPOINT_JOB_ID="$(sbatch --parsable \
  --job-name=mge-ref-check --dependency="afterok:${P0_JOB_ID}" \
  --account="${ACCOUNT}" --partition="${PARTITION}" --nodes=1 --ntasks=1 \
  --cpus-per-task="${CURATION_CHECKPOINT_CPUS}" --mem="${CURATION_CHECKPOINT_MEMORY}" --time="${CURATION_CHECKPOINT_WALLTIME}" \
  --output="${RUN_DIR}/logs/checkpoint-%j.out" --error="${RUN_DIR}/logs/checkpoint-%j.err" \
  --export="${EXPORTS}" "${CHECKPOINT_JOB_FILE}")"
printf 'checkpoint\t%s\n' "${CHECKPOINT_JOB_ID}" >> "${RUN_DIR}/provenance/job_ids.tsv"

echo "Submitted source audit job ${SOURCE_JOB_ID}"
echo "Submitted P0 inspection array ${P0_JOB_ID}"
echo "Submitted checkpoint job ${CHECKPOINT_JOB_ID}"
echo "Run directory: ${RUN_DIR}"
