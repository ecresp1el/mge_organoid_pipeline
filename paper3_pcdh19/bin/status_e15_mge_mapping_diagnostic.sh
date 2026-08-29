#!/usr/bin/env bash
set -Eeuo pipefail

BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
# shellcheck disable=SC1090
source "${BUNDLE_DIR}/config/greatlakes.env"
STEP_ROOT="${PAPER3_ROOT}/results/01_e15_mge_mapping_diagnostic"
if [[ $# -eq 1 ]]; then
  RUN_DIR="${STEP_ROOT}/$1"
else
  RUN_DIR="$(find "${STEP_ROOT}" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort | tail -n 1)"
  RUN_DIR="${STEP_ROOT}/${RUN_DIR}"
fi
[[ -f "${RUN_DIR}/provenance/job_ids.tsv" ]] || { echo "Mapping run not found: ${RUN_DIR}" >&2; exit 2; }
echo "Run: ${RUN_DIR}"
printf '%-20s %-12s %-12s %-12s %s\n' STAGE JOB_ID STATE EXIT ELAPSED
while IFS=$'\t' read -r stage job_id; do
  [[ "${stage}" == stage ]] && continue
  record="$(sacct -n -X -j "${job_id}" --format=State,ExitCode,Elapsed -P | head -n 1)"
  IFS='|' read -r state exit_code elapsed _ <<< "${record}"
  printf '%-20s %-12s %-12s %-12s %s\n' "${stage}" "${job_id}" "${state:-UNKNOWN}" "${exit_code:--}" "${elapsed:--}"
done < "${RUN_DIR}/provenance/job_ids.tsv"

for marker in PREPARE_SUCCESS.txt BANDLER_MIND_SUCCESS.txt MAPMYCELLS_SUCCESS.txt REPORT_SUCCESS.txt SUCCESS.txt; do
  if [[ -f "${RUN_DIR}/${marker}" ]]; then
    echo "present: ${marker}"
  fi
done
