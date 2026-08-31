#!/usr/bin/env bash
# Instrumented, serial, scores-only scDblFinder diagnostic.
# This job reads production inputs but writes only to its isolated run folder.
set -Eeuo pipefail

: "${PAPER3_DIAGNOSTIC_RUN_DIR:?Diagnostic run directory is required}"
: "${PAPER3_DIAGNOSTIC_R_SCRIPT:?Frozen diagnostic R script is required}"
: "${PAPER3_DIAGNOSTIC_BRIDGE_H5:?Read-only production bridge is required}"
: "${PAPER3_DIAGNOSTIC_METADATA_TSV:?Read-only production metadata is required}"
: "${PAPER3_DIAGNOSTIC_PREPARE_CHECKS:?Read-only production validation table is required}"
: "${PAPER3_DIAGNOSTIC_INPUT_H5AD:?Approved Step 02 input is required}"
: "${PAPER3_DIAGNOSTIC_R_MODULE:?R module is required}"
: "${PAPER3_DIAGNOSTIC_R_LIBRARY:?R library is required}"

RUN_DIR="${PAPER3_DIAGNOSTIC_RUN_DIR}"
R_STDOUT="${RUN_DIR}/logs/r.stdout.log"
R_STDERR="${RUN_DIR}/logs/r.stderr.log"
HEARTBEAT="${RUN_DIR}/logs/resource_heartbeat.tsv"

failed() {
  local line="$1" status="$2"
  printf 'DIAGNOSTIC_FAILED\njob_id=%s\nline=%s\nexit_status=%s\ntime=%s\ncells_removed=0\npublished_as_step03=FALSE\n' \
    "${SLURM_JOB_ID:-unknown}" "${line}" "${status}" "$(date --iso-8601=seconds)" > "${RUN_DIR}/DIAGNOSTIC_FAILED.txt"
}
trap 'failed "${LINENO}" "$?"' ERR

module load "${PAPER3_DIAGNOSTIC_R_MODULE}"
export R_LIBS_USER="${PAPER3_DIAGNOSTIC_R_LIBRARY}"
export LC_ALL=C
export TZ=UTC
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
mkdir -p "${RUN_DIR}/logs"
: > "${R_STDOUT}"
: > "${R_STDERR}"

printf 'timestamp\tjob_id\telapsed\tCPUTime\tAveCPU\tRSS_kb\tMaxRSS\tVMSize_kb\tMaxVMSize\tstate\tnode\tpid\n' > "${HEARTBEAT}"

heartbeat() {
  local r_pid="$1"
  while kill -0 "${r_pid}" 2>/dev/null; do
    local timestamp elapsed state node stats ave_cpu max_rss max_vm_size ps_values cpu_time rss_kb vm_size_kb
    timestamp="$(date --iso-8601=seconds)"
    elapsed="$(squeue -h -j "${SLURM_JOB_ID}" -o '%M' | head -n 1)"
    state="$(squeue -h -j "${SLURM_JOB_ID}" -o '%T' | head -n 1)"
    node="$(squeue -h -j "${SLURM_JOB_ID}" -o '%N' | head -n 1)"
    stats="$(sstat -j "${SLURM_JOB_ID}.batch" --noheader --parsable2 --format=AveCPU,MaxRSS,MaxVMSize 2>/dev/null | head -n 1 || true)"
    IFS='|' read -r ave_cpu max_rss max_vm_size <<< "${stats}"
    ps_values="$(ps -p "${r_pid}" -o time=,rss=,vsz= 2>/dev/null | awk '{$1=$1; print}' || true)"
    read -r cpu_time rss_kb vm_size_kb <<< "${ps_values}"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${timestamp}" "${SLURM_JOB_ID}" "${elapsed:-NA}" "${cpu_time:-NA}" "${ave_cpu:-NA}" \
      "${rss_kb:-NA}" "${max_rss:-NA}" "${vm_size_kb:-NA}" "${max_vm_size:-NA}" \
      "${state:-NA}" "${node:-NA}" "${r_pid}" >> "${HEARTBEAT}"
    sleep 60
  done
}

stdbuf -oL -eL Rscript --vanilla "${PAPER3_DIAGNOSTIC_R_SCRIPT}" \
  --bridge-h5 "${PAPER3_DIAGNOSTIC_BRIDGE_H5}" \
  --metadata-tsv "${PAPER3_DIAGNOSTIC_METADATA_TSV}" \
  --prepare-checks-tsv "${PAPER3_DIAGNOSTIC_PREPARE_CHECKS}" \
  --input-h5ad "${PAPER3_DIAGNOSTIC_INPUT_H5AD}" \
  --expected-input-sha256 "${PAPER3_DIAGNOSTIC_EXPECTED_INPUT_SHA256}" \
  --output-dir "${RUN_DIR}/results" \
  --expected-cells "${PAPER3_DIAGNOSTIC_EXPECTED_CELLS}" \
  --expected-genes "${PAPER3_DIAGNOSTIC_EXPECTED_GENES}" \
  --capture-id "${PAPER3_DIAGNOSTIC_CAPTURE_ID}" \
  --primary-seed "${PAPER3_DIAGNOSTIC_PRIMARY_SEED}" \
  > >(stdbuf -oL tee -a "${R_STDOUT}") \
  2> >(stdbuf -oL tee -a "${R_STDERR}" >&2) &
R_PID=$!

heartbeat "${R_PID}" &
HEARTBEAT_PID=$!

set +e
wait "${R_PID}"
R_STATUS=$?
wait "${HEARTBEAT_PID}"
HEARTBEAT_STATUS=$?
set -e

if [[ "${R_STATUS}" -ne 0 ]]; then
  echo "Diagnostic R process failed with status ${R_STATUS}" >&2
  exit "${R_STATUS}"
fi
if [[ "${HEARTBEAT_STATUS}" -ne 0 ]]; then
  echo "Diagnostic heartbeat failed with status ${HEARTBEAT_STATUS}" >&2
  exit "${HEARTBEAT_STATUS}"
fi

printf 'DIAGNOSTIC_SUCCESS\njob_id=%s\ncells_removed=0\npublished_as_step03=FALSE\ncompleted=%s\n' \
  "${SLURM_JOB_ID}" "$(date --iso-8601=seconds)" > "${RUN_DIR}/DIAGNOSTIC_SUCCESS.txt"
rm -f "${RUN_DIR}/DIAGNOSTIC_FAILED.txt"
