#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  echo "Usage: $0 00_input_audit [--dry-run] [--replace-run RUN_ID] [--svg true|false]" >&2
}

if [[ $# -lt 1 ]]; then
  usage
  exit 2
fi

ORIGINAL_ARGS=("$@")
STEP="$1"
shift
DRY_RUN=false
REPLACE_RUN_ID=""
SVG_OVERRIDE=""
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
    --svg)
      [[ $# -ge 2 ]] || { usage; exit 2; }
      SVG_OVERRIDE="$2"
      shift 2
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
INPUT_MANIFEST="${BUNDLE_DIR}/config/input_objects.tsv"
CLUSTER_OVERRIDES="${BUNDLE_DIR}/config/cluster_label_overrides.tsv"

# shellcheck disable=SC1090
source "${CONFIG}"
: "${REPO_ROOT:?REPO_ROOT is required}"
: "${PAPER2_ROOT:?PAPER2_ROOT is required}"
: "${ACCOUNT:?ACCOUNT is required}"
: "${PARTITION:?PARTITION is required}"

MAKE_SVG="${SVG_OVERRIDE:-${MAKE_SVG:-false}}"
[[ "${MAKE_SVG}" == "true" || "${MAKE_SVG}" == "false" ]] || {
  echo "--svg must be true or false; got: ${MAKE_SVG}" >&2
  exit 2
}
: "${PNG_DPI:?PNG_DPI is required}"
: "${PDF_DPI:?PDF_DPI is required}"
: "${SVG_DPI:?SVG_DPI is required}"
for dpi_name in PNG_DPI PDF_DPI SVG_DPI; do
  [[ "${!dpi_name}" == "300" ]] || {
    echo "${dpi_name} must remain 300; got ${!dpi_name}" >&2
    exit 2
  }
done

case "${STEP}" in
  00_input_audit)
    STEP_GROUP="00_input_audit"
    SCRIPT="${BUNDLE_DIR}/scripts/00_audit_input_objects.R"
    PLOT_SCRIPT="${BUNDLE_DIR}/scripts/00_render_input_umap_inventory.py"
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

for required in "${CONFIG}" "${INPUT_MANIFEST}" "${CLUSTER_OVERRIDES}" "${SCRIPT}" "${PLOT_SCRIPT}" "${SBATCH}"; do
  [[ -f "${required}" ]] || { echo "Missing required file: ${required}" >&2; exit 2; }
done

STAMP="$(date +%Y%m%d_%H%M%S)"
GIT_SHORT="$(git -C "${REPO_ROOT}" rev-parse --short HEAD)"
OUTPUT_MODE="versioned"
RUN_ID="${STEP}_${STAMP}_${GIT_SHORT}"
if [[ -n "${REPLACE_RUN_ID}" ]]; then
  OUTPUT_MODE="replace"
  RUN_ID="${REPLACE_RUN_ID}"
  if [[ "${RUN_ID}" == *[!A-Za-z0-9._-]* || "${RUN_ID}" != "${STEP}"_* ]]; then
    echo "Unsafe replacement run ID: ${RUN_ID}" >&2
    echo "It must be a basename beginning with ${STEP}_ and contain only letters, numbers, dot, underscore, or hyphen." >&2
    exit 2
  fi
fi
RUN_DIR="${PAPER2_ROOT}/results/${STEP_GROUP}/${RUN_ID}"
JOB_FILE="${PAPER2_ROOT}/jobs/${RUN_ID}.sbatch"
EXPECTED_PARENT="${PAPER2_ROOT}/results/${STEP_GROUP}"
if [[ "$(dirname "$(realpath -m "${RUN_DIR}")")" != "$(realpath -m "${EXPECTED_PARENT}")" ]]; then
  echo "Resolved run directory escaped its permitted parent: ${RUN_DIR}" >&2
  exit 2
fi

PREVIOUS_JOB_ID=""
PREVIOUS_JOB_STATE=""
if [[ "${OUTPUT_MODE}" == "replace" && -f "${RUN_DIR}/provenance/job_id.txt" ]]; then
  PREVIOUS_JOB_ID="$(<"${RUN_DIR}/provenance/job_id.txt")"
  PREVIOUS_JOB_STATE="$(squeue -h -j "${PREVIOUS_JOB_ID}" -o '%T' | head -n 1)"
  if [[ -n "${PREVIOUS_JOB_STATE}" ]]; then
    echo "Refusing to replace run ${RUN_ID}; job ${PREVIOUS_JOB_ID} is still ${PREVIOUS_JOB_STATE}." >&2
    exit 2
  fi
fi

echo "Step: ${STEP}"
echo "Repository: ${REPO_ROOT}"
echo "Paper 2 root: ${PAPER2_ROOT}"
echo "Input registry: ${INPUT_MANIFEST}"
echo "Prospective run directory: ${RUN_DIR}"
echo "Output mode: ${OUTPUT_MODE}"
echo "SVG enabled: ${MAKE_SVG} (SVG rasterized layers: ${SVG_DPI} dpi)"
echo "PNG/PDF: always enabled (${PNG_DPI}/${PDF_DPI} dpi)"
if [[ -n "${PREVIOUS_JOB_ID}" ]]; then
  echo "Completed prior job to replace: ${PREVIOUS_JOB_ID}"
fi
echo "Resources: account=${ACCOUNT} partition=${PARTITION} cpus=${CPUS} mem=${MEMORY} time=${WALLTIME}"

if [[ "${DRY_RUN}" == "true" ]]; then
  echo "Dry run passed; nothing was created or submitted."
  exit 0
fi

"${BUNDLE_DIR}/bin/initialize_turbo.sh" >/dev/null

if [[ "${OUTPUT_MODE}" == "versioned" ]]; then
  if [[ -e "${RUN_DIR}" || -e "${JOB_FILE}" ]]; then
    echo "Refusing to reuse an existing versioned run or job path." >&2
    exit 2
  fi
elif [[ -d "${RUN_DIR}" ]]; then
  echo "Replacing contents of completed working run: ${RUN_DIR}"
  find "${RUN_DIR}" -mindepth 1 -depth -delete
elif [[ -e "${RUN_DIR}" ]]; then
  echo "Replacement target exists but is not a directory: ${RUN_DIR}" >&2
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
cp -p "${PLOT_SCRIPT}" "${RUN_DIR}/code/"
cp -p "${SBATCH}" "${RUN_DIR}/code/"
cp -p "${SBATCH}" "${JOB_FILE}"
cp -p "${CONFIG}" "${RUN_DIR}/config/submitted_greatlakes.env"
cp -p "${INPUT_MANIFEST}" "${RUN_DIR}/config/input_objects.tsv"
cp -p "${CLUSTER_OVERRIDES}" "${RUN_DIR}/config/cluster_label_overrides.tsv"
cp -p "${BUNDLE_DIR}/templates/OUTPUT_PACKAGE_README.md" "${RUN_DIR}/README.md"

{
  echo "REPO_ROOT=${REPO_ROOT}"
  echo "PAPER2_ROOT=${PAPER2_ROOT}"
  echo "STEP=${STEP}"
  echo "RUN_ID=${RUN_ID}"
  echo "RUN_DIR=${RUN_DIR}"
  echo "OUTPUT_MODE=${OUTPUT_MODE}"
  echo "REPLACED_PREVIOUS_JOB_ID=${PREVIOUS_JOB_ID}"
  echo "ACCOUNT=${ACCOUNT}"
  echo "PARTITION=${PARTITION}"
  echo "CPUS=${CPUS}"
  echo "MEMORY=${MEMORY}"
  echo "WALLTIME=${WALLTIME}"
  echo "SEURAT5_MODULE=${SEURAT5_MODULE}"
  echo "CONDA_ENV_BIN=${CONDA_ENV_BIN}"
  echo "MAKE_SVG=${MAKE_SVG}"
  echo "PNG_DPI=${PNG_DPI}"
  echo "PDF_DPI=${PDF_DPI}"
  echo "SVG_DPI=${SVG_DPI}"
} > "${RUN_DIR}/config/resolved.env"

{
  echo "submitted_at=$(date --iso-8601=seconds)"
  echo "submit_host=$(hostname)"
  echo "submit_user=${USER:-unknown}"
  printf 'command='
  printf '%q ' "$0" "${ORIGINAL_ARGS[@]}"
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
  --export="ALL,PAPER2_RUN_DIR=${RUN_DIR},PAPER2_INPUT_MANIFEST=${RUN_DIR}/config/input_objects.tsv,PAPER2_CLUSTER_OVERRIDES=${RUN_DIR}/config/cluster_label_overrides.tsv,SEURAT5_MODULE=${SEURAT5_MODULE},CONDA_ENV_BIN=${CONDA_ENV_BIN},MAKE_SVG=${MAKE_SVG},PNG_DPI=${PNG_DPI},PDF_DPI=${PDF_DPI},SVG_DPI=${SVG_DPI}" \
  "${JOB_FILE}")"

printf '%s\n' "${JOB_ID}" > "${RUN_DIR}/provenance/job_id.txt"
for extension in out err; do
  convenience_log="${PAPER2_ROOT}/logs/${RUN_ID}.${extension}"
  if [[ -e "${convenience_log}" && ! -L "${convenience_log}" ]]; then
    echo "Refusing to replace non-symlink convenience log: ${convenience_log}" >&2
    exit 2
  fi
done
ln -sfn "${RUN_DIR}/logs/slurm-${JOB_ID}.out" "${PAPER2_ROOT}/logs/${RUN_ID}.out"
ln -sfn "${RUN_DIR}/logs/slurm-${JOB_ID}.err" "${PAPER2_ROOT}/logs/${RUN_ID}.err"

echo "Submitted SLURM job ${JOB_ID}"
echo "Run directory: ${RUN_DIR}"
echo "Monitor: ${BUNDLE_DIR}/bin/status.sh ${RUN_DIR}"
