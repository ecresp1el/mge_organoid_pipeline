#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 RUN_DIRECTORY" >&2
  exit 2
fi

RUN_DIR="${1%/}"
JOB_ID_FILE="${RUN_DIR}/provenance/job_id.txt"
[[ -f "${JOB_ID_FILE}" ]] || { echo "Missing job ID: ${JOB_ID_FILE}" >&2; exit 2; }
JOB_ID="$(<"${JOB_ID_FILE}")"

echo "Run directory: ${RUN_DIR}"
sacct -j "${JOB_ID}" \
  --format=JobID,JobName%24,State,ExitCode,Elapsed,MaxRSS,AllocCPUS,ReqMem,Start,End

if [[ -f "${RUN_DIR}/SUCCESS.txt" ]]; then
  echo
  cat "${RUN_DIR}/SUCCESS.txt"
elif [[ -f "${RUN_DIR}/FAILED.txt" ]]; then
  echo
  cat "${RUN_DIR}/FAILED.txt"
else
  echo
  echo "Run is pending or active; no terminal marker exists yet."
fi

for log_file in "${RUN_DIR}/logs/slurm-${JOB_ID}.out" "${RUN_DIR}/logs/slurm-${JOB_ID}.err"; do
  if [[ -f "${log_file}" ]]; then
    echo
    echo "Last 30 lines: ${log_file}"
    tail -n 30 "${log_file}"
  fi
done

