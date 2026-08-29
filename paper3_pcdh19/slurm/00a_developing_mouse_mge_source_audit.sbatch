#!/usr/bin/env bash
# Scheduler resources and logs are assigned by the submission wrapper.
set -Eeuo pipefail

: "${PAPER3_CURATION_RUN_DIR:?submit through bin/submit_developing_mouse_mge_reference_curation.sh}"
: "${PAPER3_CURATION_SOURCE_ROOT:?source cache root is required}"
: "${PAPER3_CURATION_REGISTRY:?submitted source registry is required}"
: "${PAPER3_CURATION_PYTHON_BIN:?Python interpreter is required}"

RUN_DIR="${PAPER3_CURATION_RUN_DIR}"
SCRIPT="${RUN_DIR}/code/developing_mouse_mge_reference_curation.py"

failed() {
  local line="$1" status="$2"
  printf 'FAILED\nstage=source_audit\njob_id=%s\nline=%s\nexit_status=%s\ntime=%s\n' \
    "${SLURM_JOB_ID:-unknown}" "${line}" "${status}" "$(date --iso-8601=seconds)" > "${RUN_DIR}/FAILED.txt"
}
trap 'failed "${LINENO}" "$?"' ERR

{
  echo "started=$(date --iso-8601=seconds)"
  echo "host=$(hostname)"
  echo "slurm_job_id=${SLURM_JOB_ID:-unknown}"
  echo "python=$(${PAPER3_CURATION_PYTHON_BIN} --version 2>&1)"
} > "${RUN_DIR}/provenance/source_audit_environment.txt"

"${PAPER3_CURATION_PYTHON_BIN}" "${SCRIPT}" source-audit \
  --registry "${PAPER3_CURATION_REGISTRY}" \
  --run-dir "${RUN_DIR}" \
  --source-root "${PAPER3_CURATION_SOURCE_ROOT}"

printf 'PASS\ncompleted=%s\njob_id=%s\n' "$(date --iso-8601=seconds)" "${SLURM_JOB_ID:-unknown}" \
  > "${RUN_DIR}/provenance/source_audit_SUCCESS.txt"
rm -f "${RUN_DIR}/FAILED.txt"
