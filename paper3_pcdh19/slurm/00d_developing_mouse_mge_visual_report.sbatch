#!/usr/bin/env bash
# Scheduler resources and logs are assigned by the submission command.
set -Eeuo pipefail

: "${PAPER3_CURATION_RUN_DIR:?completed curation run directory is required}"
: "${PAPER3_CURATION_PYTHON_BIN:?Python interpreter is required}"

RUN_DIR="${PAPER3_CURATION_RUN_DIR}"
PLOT_SCRIPT="${RUN_DIR}/code/plot_reference_annotation_hierarchy.py"
REPORT_SCRIPT="${RUN_DIR}/code/generate_reference_curation_report.py"

failed() {
  local line="$1" status="$2"
  printf 'FAILED\nstage=hierarchy_visual_report\njob_id=%s\nline=%s\nexit_status=%s\ntime=%s\n' \
    "${SLURM_JOB_ID:-unknown}" "${line}" "${status}" "$(date --iso-8601=seconds)" \
    > "${RUN_DIR}/VISUAL_REPORT_FAILED.txt"
}
trap 'failed "${LINENO}" "$?"' ERR

"${PAPER3_CURATION_PYTHON_BIN}" "${PLOT_SCRIPT}" --run-dir "${RUN_DIR}"
"${PAPER3_CURATION_PYTHON_BIN}" "${REPORT_SCRIPT}" --run-dir "${RUN_DIR}"

rm -f "${RUN_DIR}/VISUAL_REPORT_FAILED.txt"
printf 'PASS\ncompleted=%s\njob_id=%s\n' "$(date --iso-8601=seconds)" "${SLURM_JOB_ID:-unknown}" \
  > "${RUN_DIR}/VISUAL_REPORT_SUCCESS.txt"
