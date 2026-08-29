#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat >&2 <<'EOF'
Usage: submit_bandler_e15_barcode_recovery.sh [--dry-run] RUN_DIRECTORY

Capture public expression vectors on the submission host, freeze all executable
files, and submit the deposited-barcode recovery as an overwritable follow-up.
EOF
}

DRY_RUN=false
if [[ "${1:-}" == --dry-run ]]; then DRY_RUN=true; shift; fi
[[ $# -eq 1 ]] || { usage; exit 2; }

BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
RUN_DIR="$(realpath -m "$1")"
GREATLAKES_CONFIG="${BUNDLE_DIR}/config/greatlakes.env"
CURATION_CONFIG="${BUNDLE_DIR}/config/developing_mouse_mge_reference_curation.env"
PY_SCRIPT="${BUNDLE_DIR}/scripts/reference_curation/recover_bandler_e15_barcodes.py"
R_SCRIPT="${BUNDLE_DIR}/scripts/reference_curation/export_bandler_e15_fingerprints.R"
REPORT_SCRIPT="${BUNDLE_DIR}/scripts/reference_curation/generate_reference_curation_report.py"
SBATCH_SCRIPT="${BUNDLE_DIR}/slurm/00g_recover_bandler_e15_barcodes.sbatch"

# shellcheck disable=SC1090
source "${GREATLAKES_CONFIG}"
# shellcheck disable=SC1090
source "${CURATION_CONFIG}"
: "${PAPER3_ROOT:?PAPER3_ROOT is required}"
: "${ACCOUNT:?ACCOUNT is required}"
: "${PARTITION:?PARTITION is required}"
: "${CURATION_PYTHON_BIN:?CURATION_PYTHON_BIN is required}"
: "${CURATION_R_MODULE:?CURATION_R_MODULE is required}"
: "${CURATION_BARCODE_CPUS:?CURATION_BARCODE_CPUS is required}"
: "${CURATION_BARCODE_MEMORY:?CURATION_BARCODE_MEMORY is required}"
: "${CURATION_BARCODE_WALLTIME:?CURATION_BARCODE_WALLTIME is required}"

STEP_ROOT="$(realpath -m "${PAPER3_ROOT}/results/00_developing_mouse_mge_reference_curation")"
[[ "$(dirname "${RUN_DIR}")" == "${STEP_ROOT}" ]] || { echo "Run is outside the permitted curation step: ${RUN_DIR}" >&2; exit 2; }
[[ -d "${RUN_DIR}" && -f "${RUN_DIR}/SUCCESS.txt" ]] || { echo "Completed curation run not found: ${RUN_DIR}" >&2; exit 2; }
for required in "${PY_SCRIPT}" "${R_SCRIPT}" "${REPORT_SCRIPT}" "${SBATCH_SCRIPT}" "${CURATION_PYTHON_BIN}"; do
  [[ -f "${required}" ]] || { echo "Missing required file: ${required}" >&2; exit 2; }
done

if [[ -f "${RUN_DIR}/provenance/job_ids.tsv" ]]; then
  while IFS=$'\t' read -r stage job_id; do
    [[ "${stage}" == bandler_e15_barcode_recovery ]] || continue
    state="$(squeue -h -j "${job_id}" -o '%T' | head -n 1)"
    [[ -z "${state}" ]] || { echo "Barcode recovery job ${job_id} is still ${state}; refusing a concurrent overwrite." >&2; exit 2; }
  done < "${RUN_DIR}/provenance/job_ids.tsv"
fi

SOURCE_ROOT="${PAPER3_ROOT}/inputs/developing_mouse_mge/Bandler2022/source"
OUTPUT_DIR="${RUN_DIR}/Bandler2022/interactive_atlas/barcode_recovery"
echo "Run directory: ${RUN_DIR}"
echo "Immutable source cache: ${SOURCE_ROOT}"
echo "Output package: ${OUTPUT_DIR}"
echo "Rerun behavior: overwrite only this barcode_recovery package and refresh the main report"
echo "Resources: ${CURATION_BARCODE_CPUS} CPU, ${CURATION_BARCODE_MEMORY}, ${CURATION_BARCODE_WALLTIME}"
"${CURATION_PYTHON_BIN}" -m py_compile "${PY_SCRIPT}" "${REPORT_SCRIPT}"
bash -n "${SBATCH_SCRIPT}" "${BASH_SOURCE[0]}"
if [[ "${DRY_RUN}" == true ]]; then
  echo "Dry run passed; nothing was copied, overwritten, downloaded, or submitted."
  exit 0
fi

mkdir -p "${RUN_DIR}/code" "${RUN_DIR}/logs" "${RUN_DIR}/provenance" "${SOURCE_ROOT}"
if [[ -d "${OUTPUT_DIR}" ]]; then
  find "${OUTPUT_DIR}" -mindepth 1 -depth -delete
fi
mkdir -p "${OUTPUT_DIR}/audit" "${OUTPUT_DIR}/metadata" "${OUTPUT_DIR}/figures" "${OUTPUT_DIR}/public_features"
cp -p "${PY_SCRIPT}" "${R_SCRIPT}" "${REPORT_SCRIPT}" "${SBATCH_SCRIPT}" "${BASH_SOURCE[0]}" "${RUN_DIR}/code/"

echo "Acquiring deposited E15 matrices and intended public feature PDFs on the submission host"
"${CURATION_PYTHON_BIN}" "${RUN_DIR}/code/recover_bandler_e15_barcodes.py" acquire \
  --run-dir "${RUN_DIR}" --source-root "${SOURCE_ROOT}"

EXPORTS="ALL,PAPER3_CURATION_RUN_DIR=${RUN_DIR},PAPER3_CURATION_SOURCE_ROOT=${SOURCE_ROOT},PAPER3_CURATION_PYTHON_BIN=${CURATION_PYTHON_BIN},PAPER3_CURATION_R_MODULE=${CURATION_R_MODULE}"
JOB_ID="$(sbatch --parsable \
  --job-name=mge-bandler-bc --account="${ACCOUNT}" --partition="${PARTITION}" \
  --nodes=1 --ntasks=1 --cpus-per-task="${CURATION_BARCODE_CPUS}" \
  --mem="${CURATION_BARCODE_MEMORY}" --time="${CURATION_BARCODE_WALLTIME}" \
  --output="${RUN_DIR}/logs/bandler-barcode-%j.out" \
  --error="${RUN_DIR}/logs/bandler-barcode-%j.err" \
  --export="${EXPORTS}" "${RUN_DIR}/code/00g_recover_bandler_e15_barcodes.sbatch")"
printf 'bandler_e15_barcode_recovery\t%s\n' "${JOB_ID}" >> "${RUN_DIR}/provenance/job_ids.tsv"
echo "Submitted Bandler E15 barcode recovery job ${JOB_ID}"

