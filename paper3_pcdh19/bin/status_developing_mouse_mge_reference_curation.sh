#!/usr/bin/env bash
set -Eeuo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 RUN_DIRECTORY" >&2
  exit 2
fi

RUN_DIR="${1%/}"
JOB_TABLE="${RUN_DIR}/provenance/job_ids.tsv"
[[ -f "${JOB_TABLE}" ]] || { echo "Missing job table: ${JOB_TABLE}" >&2; exit 2; }

echo "Run directory: ${RUN_DIR}"
while IFS=$'\t' read -r stage job_id; do
  [[ "${stage}" != stage ]] || continue
  echo
  echo "Stage: ${stage}; job: ${job_id}"
  sacct -j "${job_id}" --format=JobID,JobName%22,State,ExitCode,Elapsed,MaxRSS,AllocCPUS,ReqMem,Start,End
done < "${JOB_TABLE}"

echo
if [[ -f "${RUN_DIR}/SUCCESS.txt" ]]; then
  cat "${RUN_DIR}/SUCCESS.txt"
elif [[ -f "${RUN_DIR}/FAILED.txt" ]]; then
  cat "${RUN_DIR}/FAILED.txt"
else
  echo "Run is pending, active, or dependency-blocked; no terminal checkpoint marker exists yet."
fi

find "${RUN_DIR}/logs" -maxdepth 1 -type f -printf '%T@\t%p\n' 2>/dev/null \
  | sort -n | tail -n 6 | cut -f2- | while IFS= read -r log_file; do
      echo
      echo "Last 20 lines: ${log_file}"
      tail -n 20 "${log_file}"
    done
