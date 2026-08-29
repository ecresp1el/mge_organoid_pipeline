#!/usr/bin/env bash
set -Eeuo pipefail

: "${PAPER3_MAPPING_RUN_DIR:?}"
: "${PAPER3_CELLRANGER_ROOT:?}"
: "${PAPER3_MAPMYCELLS_PYTHON:?}"

RUN_DIR="${PAPER3_MAPPING_RUN_DIR}"
trap 'printf "FAILED\nstage=prepare_query\njob_id=%s\nline=%s\n" "${SLURM_JOB_ID:-unknown}" "${LINENO}" > "${RUN_DIR}/PREPARE_FAILED.txt"' ERR

"${PAPER3_MAPMYCELLS_PYTHON}" "${RUN_DIR}/code/prepare_e15_mge_mapping_query.py" \
  --cellranger-root "${PAPER3_CELLRANGER_ROOT}" \
  --sample-key "${RUN_DIR}/config/sample_key.csv" \
  --output-dir "${RUN_DIR}/query"

rm -f "${RUN_DIR}/PREPARE_FAILED.txt"
printf 'PASS\njob_id=%s\ncompleted=%s\n' "${SLURM_JOB_ID:-unknown}" "$(date --iso-8601=seconds)" > "${RUN_DIR}/PREPARE_SUCCESS.txt"
