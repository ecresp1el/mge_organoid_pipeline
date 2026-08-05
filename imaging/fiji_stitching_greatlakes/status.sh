#!/usr/bin/env bash
set -euo pipefail
if [[ $# -ne 1 ]]; then
  echo "Usage: $0 RUN_DIRECTORY" >&2
  exit 2
fi
RUN_DIR="$1"
JOB_ID="$(<"${RUN_DIR}/config/job_id.txt")"
sacct -j "${JOB_ID}" --format=JobID,JobName%20,State,ExitCode,Elapsed,MaxRSS
if [[ -f "${RUN_DIR}/SUCCESS.txt" ]]; then
  cat "${RUN_DIR}/SUCCESS.txt"
  cat "${RUN_DIR}/validation.json"
elif [[ -f "${RUN_DIR}/FAILED.txt" ]]; then
  cat "${RUN_DIR}/FAILED.txt"
else
  echo "Run is pending or active."
fi
