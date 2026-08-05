#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  echo "Usage: $0 CONFIG.env [--dry-run]" >&2
}

if [[ $# -lt 1 || $# -gt 2 ]]; then
  usage
  exit 2
fi

BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
CONFIG_INPUT="$1"
DRY_RUN=false
if [[ $# -eq 2 ]]; then
  if [[ "$2" != "--dry-run" ]]; then
    usage
    exit 2
  fi
  DRY_RUN=true
fi
if [[ ! -f "${CONFIG_INPUT}" ]]; then
  echo "Configuration not found: ${CONFIG_INPUT}" >&2
  exit 2
fi
CONFIG_INPUT="$(cd "$(dirname "${CONFIG_INPUT}")" && pwd -P)/$(basename "${CONFIG_INPUT}")"

# shellcheck disable=SC1090
source "${CONFIG_INPUT}"
: "${MODE:?set MODE to production or smoke}"
: "${RUN_PARENT:?set RUN_PARENT}"
: "${RUN_NAME:?set RUN_NAME}"
: "${OUTPUT_PREFIX:?set OUTPUT_PREFIX}"
: "${ACCOUNT:?set ACCOUNT to a Great Lakes allocation}"
: "${PARTITION:?set PARTITION}"
: "${CPUS:?set CPUS}"
: "${MEMORY:?set MEMORY}"
: "${WALLTIME:?set WALLTIME}"
: "${JAVA_HEAP:?set JAVA_HEAP}"
: "${CHANNELS:?set CHANNELS}"
CHECKSUM_OUTPUTS="${CHECKSUM_OUTPUTS:-true}"

if [[ "${MODE}" != "production" && "${MODE}" != "smoke" ]]; then
  echo "MODE must be production or smoke, got ${MODE}" >&2
  exit 2
fi
if [[ "${MODE}" == "production" ]]; then
  : "${INPUT_DIR:?production mode requires INPUT_DIR}"
  : "${LAYOUT_FILE:?production mode requires LAYOUT_FILE}"
  if [[ "${LAYOUT_FILE}" != /* ]]; then
    LAYOUT_FILE="${BUNDLE_DIR}/${LAYOUT_FILE}"
  fi
  [[ -d "${INPUT_DIR}" ]] || { echo "Input directory not found: ${INPUT_DIR}" >&2; exit 2; }
  [[ -f "${LAYOUT_FILE}" ]] || { echo "Layout not found: ${LAYOUT_FILE}" >&2; exit 2; }
else
  INPUT_DIR=""
  LAYOUT_FILE=""
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
RUN_DIR="${RUN_PARENT%/}/${RUN_NAME}_${STAMP}"
if [[ -e "${RUN_DIR}" ]]; then
  echo "Refusing to reuse existing run directory: ${RUN_DIR}" >&2
  exit 2
fi

echo "Bundle: ${BUNDLE_DIR}"
echo "Configuration: ${CONFIG_INPUT}"
echo "New isolated run: ${RUN_DIR}"
echo "Mode: ${MODE}"
if [[ "${DRY_RUN}" == "true" ]]; then
  echo "Dry run passed; nothing was created or submitted."
  exit 0
fi

mkdir -p "${RUN_DIR}/code" "${RUN_DIR}/config" "${RUN_DIR}/inputs" "${RUN_DIR}/logs"
cp -a "${BUNDLE_DIR}/." "${RUN_DIR}/code/"
cp -a "${CONFIG_INPUT}" "${RUN_DIR}/config/submitted.env"

write_setting() {
  printf '%s=%q\n' "$1" "${!1}"
}
SETTINGS=(
  MODE INPUT_DIR LAYOUT_FILE RUN_PARENT RUN_NAME OUTPUT_PREFIX ACCOUNT PARTITION
  CPUS MEMORY WALLTIME JAVA_HEAP CHANNELS CHECKSUM_OUTPUTS
)
for setting in "${SETTINGS[@]}"; do
  write_setting "${setting}"
done > "${RUN_DIR}/config/resolved.env"

module purge
module load python3.11-anaconda/2024.02
if [[ "${MODE}" == "production" ]]; then
  python "${BUNDLE_DIR}/src/pipeline_io.py" stage-inputs \
    --source-dir "${INPUT_DIR}" --layout "${LAYOUT_FILE}" --destination "${RUN_DIR}/inputs"
fi

printf 'submitted_at=%s\nsubmit_host=%s\nsubmit_user=%s\ncommand=' \
  "$(date --iso-8601=seconds)" "$(hostname)" "${USER:-unknown}" \
  > "${RUN_DIR}/config/submission.txt"
printf '%q ' "$0" "$@" >> "${RUN_DIR}/config/submission.txt"
printf '\n' >> "${RUN_DIR}/config/submission.txt"

JOB_ID="$(sbatch --parsable \
  --job-name="fiji-${MODE}" \
  --account="${ACCOUNT}" --partition="${PARTITION}" \
  --nodes=1 --ntasks=1 --cpus-per-task="${CPUS}" \
  --mem="${MEMORY}" --time="${WALLTIME}" \
  --output="${RUN_DIR}/logs/slurm-%j.out" \
  --error="${RUN_DIR}/logs/slurm-%j.err" \
  --export="ALL,PIPELINE_RUN_DIR=${RUN_DIR}" \
  "${RUN_DIR}/code/pipeline.sbatch")"
printf '%s\n' "${JOB_ID}" > "${RUN_DIR}/config/job_id.txt"

echo "Submitted Great Lakes job ${JOB_ID}"
echo "Run directory: ${RUN_DIR}"
echo "Monitor: squeue -j ${JOB_ID}"
echo "After completion: cat ${RUN_DIR}/SUCCESS.txt ${RUN_DIR}/validation.json"
