#!/usr/bin/env bash
set -Eeuo pipefail

[[ $# -eq 1 ]] || { echo "Usage: republish_e15_mapping_diagnostic_report.sh RUN_ID" >&2; exit 2; }
BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
# shellcheck disable=SC1090
source "${BUNDLE_DIR}/config/greatlakes.env"
# shellcheck disable=SC1090
source "${BUNDLE_DIR}/config/e15_mge_mapping_diagnostic.env"
RUN_ID="$1"; STEP=01_e15_mge_mapping_diagnostic
[[ "${RUN_ID}" != *[!A-Za-z0-9._-]* && "${RUN_ID}" == "${STEP}"_* ]] || { echo "Unsafe run ID" >&2; exit 2; }
RUN_DIR="${PAPER3_ROOT}/results/${STEP}/${RUN_ID}"
[[ -f "${RUN_DIR}/SUCCESS.txt" ]] || { echo "Completed run not found" >&2; exit 2; }
if [[ -n "$(awk -F '\t' '$1=="report_republish" {print $2}' "${RUN_DIR}/provenance/job_ids.tsv")" ]]; then
  echo "A report republish is already recorded" >&2; exit 2
fi
cp -p "${RUN_DIR}/code/build_e15_mapping_diagnostic_report.py" "${RUN_DIR}/code/executed_59182287_build_e15_mapping_diagnostic_report.py"
cp -p "${BUNDLE_DIR}/scripts/reference_mapping/build_e15_mapping_diagnostic_report.py" "${RUN_DIR}/code/build_e15_mapping_diagnostic_report.py"
cp -p "${BUNDLE_DIR}/slurm/01e_republish_e15_mapping_diagnostic_report.sbatch" "${RUN_DIR}/code/01e_republish_e15_mapping_diagnostic_report.sbatch"
cp -p "${BASH_SOURCE[0]}" "${RUN_DIR}/code/"
EXPORTS="ALL,PAPER3_MAPPING_RUN_DIR=${RUN_DIR},PAPER3_MAPMYCELLS_PYTHON=${MAPMYCELLS_PYTHON_BIN},PAPER3_MIND_THRESHOLD=${MIND_UNASSIGNED_THRESHOLD}"
JOB_ID="$(sbatch --parsable --job-name=p3-map-report-pub --account="${ACCOUNT}" --partition="${PARTITION}" \
  --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=64G --time=03:00:00 \
  --output="${RUN_DIR}/logs/report-republish-%j.out" --error="${RUN_DIR}/logs/report-republish-%j.err" \
  --export="${EXPORTS}" "${RUN_DIR}/code/01e_republish_e15_mapping_diagnostic_report.sbatch")"
printf 'report_republish\t%s\n' "${JOB_ID}" >> "${RUN_DIR}/provenance/job_ids.tsv"
echo "Submitted diagnostic report republish ${JOB_ID}"
