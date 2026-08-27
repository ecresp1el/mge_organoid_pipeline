#!/usr/bin/env bash

set -Eeuo pipefail

usage() {
  echo "Usage: $0 [--dry-run]" >&2
}

DRY_RUN=false
ORIGINAL_ARGS=("$@")
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run)
      DRY_RUN=true
      shift
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 2
      ;;
  esac
done

BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
CONFIG="${BUNDLE_DIR}/config/greatlakes.env"
MANIFEST="${BUNDLE_DIR}/config/canonical_inputs.tsv"
INPUT_REGISTRY="${BUNDLE_DIR}/config/input_objects.tsv"
R_SCRIPT="${BUNDLE_DIR}/scripts/01_build_canonical_rds.R"
PY_SCRIPT="${BUNDLE_DIR}/scripts/01_write_validate_canonical_h5ad.py"
FINALIZER="${BUNDLE_DIR}/scripts/01_finalize_canonical_inputs.py"
ARRAY_SBATCH="${BUNDLE_DIR}/slurm/01_build_canonical_inputs_array.sbatch"
FINALIZE_SBATCH="${BUNDLE_DIR}/slurm/01_finalize_canonical_inputs.sbatch"

# shellcheck disable=SC1090
source "${CONFIG}"
: "${REPO_ROOT:?REPO_ROOT is required}"
: "${PAPER2_ROOT:?PAPER2_ROOT is required}"
: "${ACCOUNT:?ACCOUNT is required}"
: "${PARTITION:?PARTITION is required}"
: "${CANONICAL_CPUS:?CANONICAL_CPUS is required}"
: "${CANONICAL_MEMORY:?CANONICAL_MEMORY is required}"
: "${CANONICAL_WALLTIME:?CANONICAL_WALLTIME is required}"
: "${CANONICAL_ARRAY_CONCURRENCY:?CANONICAL_ARRAY_CONCURRENCY is required}"
: "${CANONICAL_FINALIZE_CPUS:?CANONICAL_FINALIZE_CPUS is required}"
: "${CANONICAL_FINALIZE_MEMORY:?CANONICAL_FINALIZE_MEMORY is required}"
: "${CANONICAL_FINALIZE_WALLTIME:?CANONICAL_FINALIZE_WALLTIME is required}"

for required in \
  "${CONFIG}" "${MANIFEST}" "${INPUT_REGISTRY}" "${R_SCRIPT}" "${PY_SCRIPT}" \
  "${FINALIZER}" "${ARRAY_SBATCH}" "${FINALIZE_SBATCH}"; do
  [[ -f "${required}" ]] || { echo "Missing required file: ${required}" >&2; exit 2; }
done

[[ "$(awk 'END {print NR - 1}' "${MANIFEST}")" == "6" ]] || {
  echo "Canonical manifest must contain exactly six data rows." >&2
  exit 2
}

INPUTS_ROOT="${PAPER2_ROOT}/inputs"
CANONICAL_DIR="${INPUTS_ROOT}/canonical"
if [[ -e "${CANONICAL_DIR}" ]]; then
  echo "Refusing to overwrite frozen canonical inputs: ${CANONICAL_DIR}" >&2
  exit 2
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
GIT_SHORT="$(git -C "${REPO_ROOT}" rev-parse --short HEAD)"
BUILD_ID="canonical_build_${STAMP}_${GIT_SHORT}"
STAGING_DIR="${INPUTS_ROOT}/.${BUILD_ID}"
ARRAY_JOB_FILE="${PAPER2_ROOT}/jobs/${BUILD_ID}_array.sbatch"
FINALIZE_JOB_FILE="${PAPER2_ROOT}/jobs/${BUILD_ID}_finalize.sbatch"

echo "Step: 01_canonical_inputs"
echo "Repository: ${REPO_ROOT}"
echo "Source manifest: ${MANIFEST}"
echo "Staging directory: ${STAGING_DIR}"
echo "Frozen destination: ${CANONICAL_DIR}"
echo "Array resources: 6 tasks, max ${CANONICAL_ARRAY_CONCURRENCY} concurrent; ${CANONICAL_CPUS} CPUs, ${CANONICAL_MEMORY}, ${CANONICAL_WALLTIME} each"
echo "Finalizer resources: ${CANONICAL_FINALIZE_CPUS} CPUs, ${CANONICAL_FINALIZE_MEMORY}, ${CANONICAL_FINALIZE_WALLTIME}"
if [[ "${DRY_RUN}" == "true" ]]; then
  echo "Dry run passed; nothing was created or submitted."
  exit 0
fi

"${BUNDLE_DIR}/bin/initialize_turbo.sh" >/dev/null
mkdir -p \
  "${STAGING_DIR}/code" \
  "${STAGING_DIR}/config" \
  "${STAGING_DIR}/logs" \
  "${STAGING_DIR}/provenance"

cp -p "${R_SCRIPT}" "${PY_SCRIPT}" "${FINALIZER}" "${ARRAY_SBATCH}" "${FINALIZE_SBATCH}" "${STAGING_DIR}/code/"
cp -p "${CONFIG}" "${STAGING_DIR}/config/submitted_greatlakes.env"
cp -p "${MANIFEST}" "${STAGING_DIR}/config/canonical_inputs.tsv"
cp -p "${INPUT_REGISTRY}" "${STAGING_DIR}/config/input_objects.tsv"
cp -p "${ARRAY_SBATCH}" "${ARRAY_JOB_FILE}"
cp -p "${FINALIZE_SBATCH}" "${FINALIZE_JOB_FILE}"

{
  echo "REPO_ROOT=${REPO_ROOT}"
  echo "PAPER2_ROOT=${PAPER2_ROOT}"
  echo "STEP=01_canonical_inputs"
  echo "BUILD_ID=${BUILD_ID}"
  echo "STAGING_DIR=${STAGING_DIR}"
  echo "CANONICAL_DIR=${CANONICAL_DIR}"
  echo "ACCOUNT=${ACCOUNT}"
  echo "PARTITION=${PARTITION}"
  echo "SEURAT5_MODULE=${SEURAT5_MODULE}"
  echo "CONDA_ENV_BIN=${CONDA_ENV_BIN}"
  echo "CANONICAL_CPUS=${CANONICAL_CPUS}"
  echo "CANONICAL_MEMORY=${CANONICAL_MEMORY}"
  echo "CANONICAL_WALLTIME=${CANONICAL_WALLTIME}"
  echo "CANONICAL_ARRAY_CONCURRENCY=${CANONICAL_ARRAY_CONCURRENCY}"
  echo "CANONICAL_FINALIZE_CPUS=${CANONICAL_FINALIZE_CPUS}"
  echo "CANONICAL_FINALIZE_MEMORY=${CANONICAL_FINALIZE_MEMORY}"
  echo "CANONICAL_FINALIZE_WALLTIME=${CANONICAL_FINALIZE_WALLTIME}"
} > "${STAGING_DIR}/config/resolved.env"

{
  echo "submitted_at=$(date --iso-8601=seconds)"
  echo "submit_host=$(hostname)"
  echo "submit_user=${USER:-unknown}"
  printf 'command='
  printf '%q ' "$0" "${ORIGINAL_ARGS[@]}"
  printf '\n'
} > "${STAGING_DIR}/provenance/submission.txt"
git -C "${REPO_ROOT}" rev-parse HEAD > "${STAGING_DIR}/provenance/git_commit.txt"
git -C "${REPO_ROOT}" status --short > "${STAGING_DIR}/provenance/git_status.txt"
git -C "${REPO_ROOT}" diff --binary > "${STAGING_DIR}/provenance/git_working_tree.patch"

STAGED_MANIFEST="${STAGING_DIR}/config/canonical_inputs.tsv"
ARRAY_JOB_ID="$(sbatch --parsable \
  --job-name="p2-canonical" \
  --account="${ACCOUNT}" \
  --partition="${PARTITION}" \
  --nodes=1 \
  --ntasks=1 \
  --cpus-per-task="${CANONICAL_CPUS}" \
  --mem="${CANONICAL_MEMORY}" \
  --time="${CANONICAL_WALLTIME}" \
  --array="0-5%${CANONICAL_ARRAY_CONCURRENCY}" \
  --output="${STAGING_DIR}/logs/build-%A_%a.out" \
  --error="${STAGING_DIR}/logs/build-%A_%a.err" \
  --export="ALL,PAPER2_CANONICAL_STAGING=${STAGING_DIR},PAPER2_CANONICAL_MANIFEST=${STAGED_MANIFEST},SEURAT5_MODULE=${SEURAT5_MODULE},CONDA_ENV_BIN=${CONDA_ENV_BIN}" \
  "${ARRAY_JOB_FILE}")"
printf '%s\n' "${ARRAY_JOB_ID}" > "${STAGING_DIR}/provenance/array_job_id.txt"

FINALIZE_JOB_ID="$(sbatch --parsable \
  --job-name="p2-canonical-freeze" \
  --dependency="afterok:${ARRAY_JOB_ID}" \
  --account="${ACCOUNT}" \
  --partition="${PARTITION}" \
  --nodes=1 \
  --ntasks=1 \
  --cpus-per-task="${CANONICAL_FINALIZE_CPUS}" \
  --mem="${CANONICAL_FINALIZE_MEMORY}" \
  --time="${CANONICAL_FINALIZE_WALLTIME}" \
  --output="${STAGING_DIR}/logs/finalize-%j.out" \
  --error="${STAGING_DIR}/logs/finalize-%j.err" \
  --export="ALL,PAPER2_CANONICAL_STAGING=${STAGING_DIR},PAPER2_CANONICAL_DIR=${CANONICAL_DIR},PAPER2_CANONICAL_MANIFEST=${STAGED_MANIFEST},CONDA_ENV_BIN=${CONDA_ENV_BIN}" \
  "${FINALIZE_JOB_FILE}")"
printf '%s\n' "${FINALIZE_JOB_ID}" > "${STAGING_DIR}/provenance/finalize_job_id.txt"

echo "Submitted canonical build array job ${ARRAY_JOB_ID}"
echo "Submitted dependent freeze/finalizer job ${FINALIZE_JOB_ID}"
echo "Staging directory: ${STAGING_DIR}"
echo "On complete success it will atomically become: ${CANONICAL_DIR}"
echo "Monitor: squeue -j ${ARRAY_JOB_ID},${FINALIZE_JOB_ID}"
