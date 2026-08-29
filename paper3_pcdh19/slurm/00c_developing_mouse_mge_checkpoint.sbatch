#!/usr/bin/env bash
# Scheduler resources, dependency, and logs are assigned by the submission wrapper.
set -Eeuo pipefail

: "${PAPER3_CURATION_RUN_DIR:?submit through bin/submit_developing_mouse_mge_reference_curation.sh}"
: "${PAPER3_CURATION_SOURCE_ROOT:?source cache root is required}"
: "${PAPER3_CURATION_PYTHON_BIN:?Python interpreter is required}"

RUN_DIR="${PAPER3_CURATION_RUN_DIR}"
SCRIPT="${RUN_DIR}/code/developing_mouse_mge_reference_curation.py"

failed() {
  local line="$1" status="$2"
  printf 'FAILED\nstage=checkpoint\njob_id=%s\nline=%s\nexit_status=%s\ntime=%s\n' \
    "${SLURM_JOB_ID:-unknown}" "${line}" "${status}" "$(date --iso-8601=seconds)" > "${RUN_DIR}/FAILED.txt"
}
trap 'failed "${LINENO}" "$?"' ERR

"${PAPER3_CURATION_PYTHON_BIN}" "${SCRIPT}" checkpoint \
  --run-dir "${RUN_DIR}" \
  --source-root "${PAPER3_CURATION_SOURCE_ROOT}"

echo "Early processed-object checkpoint complete: ${RUN_DIR}"
