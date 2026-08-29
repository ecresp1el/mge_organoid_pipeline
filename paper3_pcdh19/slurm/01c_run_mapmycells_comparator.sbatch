#!/usr/bin/env bash
set -Eeuo pipefail

: "${PAPER3_MAPPING_RUN_DIR:?}"
: "${PAPER3_MAPMYCELLS_PYTHON:?}"
: "${PAPER3_MAPMYCELLS_MARKERS:?}"
: "${PAPER3_MAPMYCELLS_STATS:?}"
: "${PAPER3_MAPMYCELLS_DROP_LEVEL:?}"
: "${PAPER3_MAPPING_SEED:?}"

RUN_DIR="${PAPER3_MAPPING_RUN_DIR}"
trap 'printf "FAILED\nstage=mapmycells\njob_id=%s\nline=%s\n" "${SLURM_JOB_ID:-unknown}" "${LINENO}" > "${RUN_DIR}/MAPMYCELLS_FAILED.txt"' ERR

export OMP_NUM_THREADS=1
"${PAPER3_MAPMYCELLS_PYTHON}" "${RUN_DIR}/code/run_mapmycells_comparator.py" \
  --query-h5ad "${RUN_DIR}/query/query_raw_counts_mapmycells.h5ad" \
  --query-metadata "${RUN_DIR}/query/query_cell_metadata.tsv.gz" \
  --marker-lookup "${PAPER3_MAPMYCELLS_MARKERS}" \
  --precomputed-stats "${PAPER3_MAPMYCELLS_STATS}" \
  --drop-level "${PAPER3_MAPMYCELLS_DROP_LEVEL}" \
  --output-dir "${RUN_DIR}/MapMyCells" \
  --cpus "${SLURM_CPUS_PER_TASK:-1}" \
  --seed "${PAPER3_MAPPING_SEED}" \
  --bootstrap-iteration 100 \
  --bootstrap-factor 0.5

rm -f "${RUN_DIR}/MAPMYCELLS_FAILED.txt"
printf 'PASS\njob_id=%s\ncompleted=%s\n' "${SLURM_JOB_ID:-unknown}" "$(date --iso-8601=seconds)" > "${RUN_DIR}/MAPMYCELLS_SUCCESS.txt"
