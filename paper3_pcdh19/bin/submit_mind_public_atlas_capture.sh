#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat >&2 <<'EOF'
Usage: submit_mind_public_atlas_capture.sh [--dry-run] RUN_DIRECTORY

Freeze and submit the public MIND-atlas capture as an overwritable follow-up
inside one completed developing-mouse MGE reference-curation run.
EOF
}

DRY_RUN=false
if [[ "${1:-}" == --dry-run ]]; then
  DRY_RUN=true
  shift
fi
[[ $# -eq 1 ]] || { usage; exit 2; }

BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
RUN_DIR="$(realpath -m "$1")"
GREATLAKES_CONFIG="${BUNDLE_DIR}/config/greatlakes.env"
CURATION_CONFIG="${BUNDLE_DIR}/config/developing_mouse_mge_reference_curation.env"
CAPTURE_SCRIPT="${BUNDLE_DIR}/scripts/reference_curation/capture_mind_shiny_public_outputs.py"
REPORT_SCRIPT="${BUNDLE_DIR}/scripts/reference_curation/generate_reference_curation_report.py"
SBATCH_SCRIPT="${BUNDLE_DIR}/slurm/00f_capture_mind_shiny_public_outputs.sbatch"

# shellcheck disable=SC1090
source "${GREATLAKES_CONFIG}"
# shellcheck disable=SC1090
source "${CURATION_CONFIG}"
: "${PAPER3_ROOT:?PAPER3_ROOT is required}"
: "${ACCOUNT:?ACCOUNT is required}"
: "${PARTITION:?PARTITION is required}"
: "${CURATION_PYTHON_BIN:?CURATION_PYTHON_BIN is required}"
: "${CURATION_ATLAS_CPUS:?CURATION_ATLAS_CPUS is required}"
: "${CURATION_ATLAS_MEMORY:?CURATION_ATLAS_MEMORY is required}"
: "${CURATION_ATLAS_WALLTIME:?CURATION_ATLAS_WALLTIME is required}"

STEP_ROOT="$(realpath -m "${PAPER3_ROOT}/results/00_developing_mouse_mge_reference_curation")"
[[ "$(dirname "${RUN_DIR}")" == "${STEP_ROOT}" ]] || {
  echo "Run is outside the permitted curation step: ${RUN_DIR}" >&2
  exit 2
}
[[ -d "${RUN_DIR}" && -f "${RUN_DIR}/SUCCESS.txt" ]] || {
  echo "Completed curation run not found: ${RUN_DIR}" >&2
  exit 2
}
for required in "${CAPTURE_SCRIPT}" "${REPORT_SCRIPT}" "${SBATCH_SCRIPT}" "${CURATION_PYTHON_BIN}"; do
  [[ -f "${required}" ]] || { echo "Missing required file: ${required}" >&2; exit 2; }
done

if [[ -f "${RUN_DIR}/provenance/job_ids.tsv" ]]; then
  while IFS=$'\t' read -r stage job_id; do
    [[ "${stage}" == mind_public_atlas_capture ]] || continue
    state="$(squeue -h -j "${job_id}" -o '%T' | head -n 1)"
    [[ -z "${state}" ]] || {
      echo "Atlas capture job ${job_id} is still ${state}; refusing a concurrent overwrite." >&2
      exit 2
    }
  done < "${RUN_DIR}/provenance/job_ids.tsv"
fi

echo "Run directory: ${RUN_DIR}"
echo "Output package: ${RUN_DIR}/Bandler2022/interactive_atlas"
echo "Rerun behavior: atomically overwrite this follow-up package and refresh the main report"
echo "Resources: ${CURATION_ATLAS_CPUS} CPU, ${CURATION_ATLAS_MEMORY}, ${CURATION_ATLAS_WALLTIME}"
"${CURATION_PYTHON_BIN}" -m py_compile "${CAPTURE_SCRIPT}" "${REPORT_SCRIPT}"
if [[ "${DRY_RUN}" == true ]]; then
  echo "Dry run passed; nothing was copied, overwritten, or submitted."
  exit 0
fi

mkdir -p "${RUN_DIR}/code" "${RUN_DIR}/logs" "${RUN_DIR}/provenance"
cp -p "${CAPTURE_SCRIPT}" "${REPORT_SCRIPT}" "${SBATCH_SCRIPT}" "${BASH_SOURCE[0]}" "${RUN_DIR}/code/"

echo "Capturing intended public Shiny PDFs from the submission host before SLURM validation"
"${CURATION_PYTHON_BIN}" "${RUN_DIR}/code/capture_mind_shiny_public_outputs.py" --run-dir "${RUN_DIR}"

EXPORTS="ALL,PAPER3_CURATION_RUN_DIR=${RUN_DIR},PAPER3_CURATION_PYTHON_BIN=${CURATION_PYTHON_BIN}"
JOB_ID="$(sbatch --parsable \
  --job-name=mge-mind-atlas --account="${ACCOUNT}" --partition="${PARTITION}" \
  --nodes=1 --ntasks=1 --cpus-per-task="${CURATION_ATLAS_CPUS}" \
  --mem="${CURATION_ATLAS_MEMORY}" --time="${CURATION_ATLAS_WALLTIME}" \
  --output="${RUN_DIR}/logs/mind-public-atlas-%j.out" \
  --error="${RUN_DIR}/logs/mind-public-atlas-%j.err" \
  --export="${EXPORTS}" "${RUN_DIR}/code/00f_capture_mind_shiny_public_outputs.sbatch")"
printf 'mind_public_atlas_capture\t%s\n' "${JOB_ID}" >> "${RUN_DIR}/provenance/job_ids.tsv"
echo "Submitted public MIND-atlas capture job ${JOB_ID}"
