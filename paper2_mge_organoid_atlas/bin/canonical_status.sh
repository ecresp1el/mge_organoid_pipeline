#!/usr/bin/env bash

set -Eeuo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 CANONICAL_STAGING_OR_FROZEN_DIRECTORY" >&2
  exit 2
fi

TARGET="${1%/}"
[[ -d "${TARGET}" ]] || { echo "Directory does not exist: ${TARGET}" >&2; exit 2; }
ARRAY_FILE="${TARGET}/provenance/array_job_id.txt"
FINALIZE_FILE="${TARGET}/provenance/finalize_job_id.txt"
[[ -f "${ARRAY_FILE}" ]] || { echo "Missing array job ID: ${ARRAY_FILE}" >&2; exit 2; }
[[ -f "${FINALIZE_FILE}" ]] || { echo "Missing finalizer job ID: ${FINALIZE_FILE}" >&2; exit 2; }
ARRAY_JOB_ID="$(<"${ARRAY_FILE}")"
FINALIZE_JOB_ID="$(<"${FINALIZE_FILE}")"

echo "Canonical build directory: ${TARGET}"
echo "Array job: ${ARRAY_JOB_ID}"
echo "Finalizer job: ${FINALIZE_JOB_ID}"
sacct -j "${ARRAY_JOB_ID},${FINALIZE_JOB_ID}" \
  --format=JobID,JobName%24,State,ExitCode,Elapsed,MaxRSS,AllocCPUS,ReqMem,Start,End

echo
echo "Per-study terminal markers:"
for study_id in varela_div30 varela_div90 walsh bershteyn_2025 bershteyn_2023 siebert_2026; do
  if [[ -s "${TARGET}/${study_id}/SUCCESS.txt" ]]; then
    echo "${study_id}: PASS"
  elif [[ -s "${TARGET}/${study_id}/FAILED.txt" ]]; then
    echo "${study_id}: FAILED"
    sed -n '1,12p' "${TARGET}/${study_id}/FAILED.txt"
  else
    echo "${study_id}: pending/running"
  fi
done

if [[ -s "${TARGET}/FROZEN.txt" ]]; then
  echo
  cat "${TARGET}/FROZEN.txt"
fi

for log_file in "${TARGET}"/logs/*.err; do
  [[ -s "${log_file}" ]] || continue
  echo
  echo "Nonempty stderr tail: ${log_file}"
  tail -n 20 "${log_file}"
done
