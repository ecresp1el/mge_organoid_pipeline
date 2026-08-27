#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  echo "Usage: $0 00_input_audit [--dry-run]" >&2
}

if [[ $# -lt 1 || $# -gt 2 ]]; then
  usage
  exit 2
fi

STEP="$1"
DRY_RUN=false
if [[ $# -eq 2 ]]; then
  [[ "$2" == "--dry-run" ]] || { usage; exit 2; }
  DRY_RUN=true
fi

BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
CONFIG="${BUNDLE_DIR}/config/greatlakes.env"
INPUT_MANIFEST="${BUNDLE_DIR}/config/input_objects.tsv"

# shellcheck disable=SC1090
source "${CONFIG}"
: "${REPO_ROOT:?REPO_ROOT is required}"
: "${PAPER2_ROOT:?PAPER2_ROOT is required}"
: "${ACCOUNT:?ACCOUNT is required}"
: "${PARTITION:?PARTITION is required}"

case "${STEP}" in
  00_input_audit)
    STEP_GROUP="00_input_audit"
    SCRIPT="${BUNDLE_DIR}/scripts/00_audit_input_objects.R"
    SBATCH="${BUNDLE_DIR}/slurm/00_audit_input_objects.sbatch"
    CPUS="${AUDIT_CPUS}"
    MEMORY="${AUDIT_MEMORY}"
    WALLTIME="${AUDIT_WALLTIME}"
    ;;
  *)
    echo "Unknown or not-yet-implemented step: ${STEP}" >&2
    usage
    exit 2
    ;;
esac

for required in "${CONFIG}" "${INPUT_MANIFEST}" "${SCRIPT}" "${SBATCH}"; do
  [[ -f "${required}" ]] || { echo "Missing required file: ${required}" >&2; exit 2; }
done

STAMP="$(date +%Y%m%d_%H%M%S)"
GIT_SHORT="$(git -C "${REPO_ROOT}" rev-parse --short HEAD)"
RUN_ID="${STEP}_${STAMP}_${GIT_SHORT}"
RUN_DIR="${PAPER2_ROOT}/results/${STEP_GROUP}/${RUN_ID}"
JOB_FILE="${PAPER2_ROOT}/jobs/${RUN_ID}.sbatch"

echo "Step: ${STEP}"
echo "Repository: ${REPO_ROOT}"
echo "Paper 2 root: ${PAPER2_ROOT}"
echo "Input registry: ${INPUT_MANIFEST}"
echo "Prospective run directory: ${RUN_DIR}"
echo "Resources: account=${ACCOUNT} partition=${PARTITION} cpus=${CPUS} mem=${MEMORY} time=${WALLTIME}"

if [[ "${DRY_RUN}" == "true" ]]; then
  echo "Dry run passed; nothing was created or submitted."
  exit 0
fi

"${BUNDLE_DIR}/bin/initialize_turbo.sh" >/dev/null

if [[ -e "${RUN_DIR}" || -e "${JOB_FILE}" ]]; then
  echo "Refusing to reuse an existing run or job path." >&2
  exit 2
fi

mkdir -p \
  "${RUN_DIR}/code" \
  "${RUN_DIR}/config" \
  "${RUN_DIR}/figures/png" \
  "${RUN_DIR}/figures/pdf" \
  "${RUN_DIR}/figures/svg" \
  "${RUN_DIR}/tables" \
  "${RUN_DIR}/logs" \
  "${RUN_DIR}/provenance"

cp -p "${SCRIPT}" "${RUN_DIR}/code/"
cp -p "${SBATCH}" "${RUN_DIR}/code/"
cp -p "${SBATCH}" "${JOB_FILE}"
cp -p "${CONFIG}" "${RUN_DIR}/config/submitted_greatlakes.env"
cp -p "${INPUT_MANIFEST}" "${RUN_DIR}/config/input_objects.tsv"
cp -p "${BUNDLE_DIR}/templates/OUTPUT_PACKAGE_README.md" "${RUN_DIR}/README.md"

{
  echo "REPO_ROOT=${REPO_ROOT}"
  echo "PAPER2_ROOT=${PAPER2_ROOT}"
  echo "STEP=${STEP}"
  echo "RUN_ID=${RUN_ID}"
  echo "RUN_DIR=${RUN_DIR}"
  echo "ACCOUNT=${ACCOUNT}"
  echo "PARTITION=${PARTITION}"
  echo "CPUS=${CPUS}"
  echo "MEMORY=${MEMORY}"
  echo "WALLTIME=${WALLTIME}"
  echo "SEURAT5_MODULE=${SEURAT5_MODULE}"
} > "${RUN_DIR}/config/resolved.env"

{
  echo "submitted_at=$(date --iso-8601=seconds)"
  echo "submit_host=$(hostname)"
  echo "submit_user=${USER:-unknown}"
  printf 'command='
  printf '%q ' "$0" "$@"
  printf '\n'
} > "${RUN_DIR}/provenance/submission.txt"

git -C "${REPO_ROOT}" rev-parse HEAD > "${RUN_DIR}/provenance/git_commit.txt"
git -C "${REPO_ROOT}" status --short > "${RUN_DIR}/provenance/git_status.txt"
git -C "${REPO_ROOT}" diff --binary > "${RUN_DIR}/provenance/git_working_tree.patch"

JOB_ID="$(sbatch --parsable \
  --job-name="p2-input-audit" \
  --account="${ACCOUNT}" \
  --partition="${PARTITION}" \
  --nodes=1 \
  --ntasks=1 \
  --cpus-per-task="${CPUS}" \
  --mem="${MEMORY}" \
  --time="${WALLTIME}" \
  --output="${RUN_DIR}/logs/slurm-%j.out" \
  --error="${RUN_DIR}/logs/slurm-%j.err" \
  --export="ALL,PAPER2_RUN_DIR=${RUN_DIR},PAPER2_INPUT_MANIFEST=${RUN_DIR}/config/input_objects.tsv,SEURAT5_MODULE=${SEURAT5_MODULE}" \
  "${JOB_FILE}")"

printf '%s\n' "${JOB_ID}" > "${RUN_DIR}/provenance/job_id.txt"
ln -s "${RUN_DIR}/logs/slurm-${JOB_ID}.out" "${PAPER2_ROOT}/logs/${RUN_ID}.out"
ln -s "${RUN_DIR}/logs/slurm-${JOB_ID}.err" "${PAPER2_ROOT}/logs/${RUN_ID}.err"

echo "Submitted SLURM job ${JOB_ID}"
echo "Run directory: ${RUN_DIR}"
echo "Monitor: ${BUNDLE_DIR}/bin/status.sh ${RUN_DIR}"
