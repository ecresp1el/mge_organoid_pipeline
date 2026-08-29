#!/usr/bin/env bash
set -Eeuo pipefail

: "${PAPER3_MAPPING_RUN_DIR:?}"
: "${PAPER3_MAPMYCELLS_PYTHON:?}"
: "${PAPER3_MIND_THRESHOLD:?}"
: "${PAPER3_MAPPING_R_MODULE:?}"

RUN_DIR="${PAPER3_MAPPING_RUN_DIR}"
trap 'printf "FAILED\nstage=report\njob_id=%s\nline=%s\n" "${SLURM_JOB_ID:-unknown}" "${LINENO}" > "${RUN_DIR}/REPORT_FAILED.txt"' ERR

[[ -f "${RUN_DIR}/Bandler_MIND/SUCCESS.txt" ]]
[[ -f "${RUN_DIR}/MapMyCells/SUCCESS.txt" ]]
module load Bioinformatics
module load "${PAPER3_MAPPING_R_MODULE}"
Rscript "${RUN_DIR}/code/attach_mapmycells_to_seurat.R" \
  --seurat "${RUN_DIR}/Bandler_MIND/paper3_query_with_mind_labels.rds" \
  --mapmycells "${RUN_DIR}/MapMyCells/mapmycells_per_cell.tsv.gz" \
  --output "${RUN_DIR}/Bandler_MIND/paper3_query_with_mind_and_mapmycells_labels.rds"
"${PAPER3_MAPMYCELLS_PYTHON}" "${RUN_DIR}/code/build_e15_mapping_diagnostic_report.py" \
  --mind "${RUN_DIR}/Bandler_MIND/mind_label_transfer_per_cell.tsv.gz" \
  --mapmycells "${RUN_DIR}/MapMyCells/mapmycells_per_cell.tsv.gz" \
  --output-dir "${RUN_DIR}" \
  --threshold "${PAPER3_MIND_THRESHOLD}"

rm -f "${RUN_DIR}/REPORT_FAILED.txt"
printf 'PASS\njob_id=%s\ncompleted=%s\n' "${SLURM_JOB_ID:-unknown}" "$(date --iso-8601=seconds)" > "${RUN_DIR}/REPORT_SUCCESS.txt"
