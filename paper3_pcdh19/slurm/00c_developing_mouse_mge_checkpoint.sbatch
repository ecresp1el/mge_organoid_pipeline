#!/usr/bin/env bash
# Scheduler resources, dependency, and logs are assigned by the submission wrapper.
set -Eeuo pipefail

: "${PAPER3_CURATION_RUN_DIR:?submit through bin/submit_developing_mouse_mge_reference_curation.sh}"
: "${PAPER3_CURATION_SOURCE_ROOT:?source cache root is required}"
: "${PAPER3_CURATION_REGISTRY:?submitted source registry is required}"
: "${PAPER3_CURATION_PYTHON_BIN:?Python interpreter is required}"

RUN_DIR="${PAPER3_CURATION_RUN_DIR}"
SCRIPT="${RUN_DIR}/code/developing_mouse_mge_reference_curation.py"
REPORT_SCRIPT="${RUN_DIR}/code/generate_reference_curation_report.py"
PLOT_SCRIPT="${RUN_DIR}/code/plot_reference_annotation_hierarchy.py"

failed() {
  local line="$1" status="$2"
  printf 'FAILED\nstage=checkpoint\njob_id=%s\nline=%s\nexit_status=%s\ntime=%s\n' \
    "${SLURM_JOB_ID:-unknown}" "${line}" "${status}" "$(date --iso-8601=seconds)" > "${RUN_DIR}/FAILED.txt"
}
trap 'failed "${LINENO}" "$?"' ERR

"${PAPER3_CURATION_PYTHON_BIN}" "${SCRIPT}" checkpoint \
  --registry "${PAPER3_CURATION_REGISTRY}" \
  --run-dir "${RUN_DIR}" \
  --source-root "${PAPER3_CURATION_SOURCE_ROOT}"

"${PAPER3_CURATION_PYTHON_BIN}" "${PLOT_SCRIPT}" --run-dir "${RUN_DIR}"
"${PAPER3_CURATION_PYTHON_BIN}" "${REPORT_SCRIPT}" --run-dir "${RUN_DIR}"

echo "Early processed-object checkpoint and observed-object report complete: ${RUN_DIR}"
