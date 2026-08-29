#!/usr/bin/env bash
set -Eeuo pipefail

: "${PAPER3_MAPPING_RUN_DIR:?}"
: "${PAPER3_MAPMYCELLS_PYTHON:?}"
: "${PAPER3_MIND_THRESHOLD:?}"
RUN_DIR="${PAPER3_MAPPING_RUN_DIR}"
trap 'printf "FAILED\nstage=report_republish\njob_id=%s\nline=%s\n" "${SLURM_JOB_ID:-unknown}" "${LINENO}" > "${RUN_DIR}/REPORT_REPUBLISH_FAILED.txt"' ERR
[[ -f "${RUN_DIR}/SUCCESS.txt" ]]
[[ -f "${RUN_DIR}/Bandler_MIND/COMBINED_SEURAT_SUCCESS.txt" ]]
"${PAPER3_MAPMYCELLS_PYTHON}" "${RUN_DIR}/code/build_e15_mapping_diagnostic_report.py" \
  --mind "${RUN_DIR}/Bandler_MIND/mind_label_transfer_per_cell.tsv.gz" \
  --mapmycells "${RUN_DIR}/MapMyCells/mapmycells_per_cell.tsv.gz" \
  --output-dir "${RUN_DIR}" \
  --threshold "${PAPER3_MIND_THRESHOLD}"
rm -f "${RUN_DIR}/REPORT_REPUBLISH_FAILED.txt"
printf 'PASS\njob_id=%s\ncompleted=%s\n' "${SLURM_JOB_ID:-unknown}" "$(date --iso-8601=seconds)" > "${RUN_DIR}/REPORT_REPUBLISH_SUCCESS.txt"
