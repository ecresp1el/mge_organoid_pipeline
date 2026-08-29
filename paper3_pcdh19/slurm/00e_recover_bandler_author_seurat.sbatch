#!/usr/bin/env bash
# Scheduler resources and logs are assigned by the submission command.
set -Eeuo pipefail

: "${PAPER3_CURATION_RUN_DIR:?completed curation run directory is required}"
: "${PAPER3_CURATION_SOURCE_ROOT:?reference source cache root is required}"
: "${PAPER3_CURATION_PYTHON_BIN:?Python interpreter is required}"
: "${PAPER3_CURATION_R_MODULE:?R/Seurat module is required}"

RUN_DIR="${PAPER3_CURATION_RUN_DIR}"
AUTHOR_OBJECT="${PAPER3_CURATION_SOURCE_ROOT}/Bandler2022/author_objects/STICR.seuratobject.RDS"
CA301_OBJECT="${PAPER3_CURATION_SOURCE_ROOT}/Bandler2022/source/GSM5684876_CA301_filtered_RNA_counts.RDS.gz"
OUTPUT_DIR="${RUN_DIR}/Bandler2022/author_object_audit"
PYTHON_SCRIPT="${RUN_DIR}/code/recover_bandler_author_object.py"
R_SCRIPT="${RUN_DIR}/code/inspect_bandler_author_seurat.R"

failed() {
  local line="$1" status="$2"
  printf 'FAILED\nstage=bandler_author_object_recovery\njob_id=%s\nline=%s\nexit_status=%s\ntime=%s\n' \
    "${SLURM_JOB_ID:-unknown}" "${line}" "${status}" "$(date --iso-8601=seconds)" \
    > "${RUN_DIR}/BANDLER_AUTHOR_OBJECT_FAILED.txt"
}
trap 'failed "${LINENO}" "$?"' ERR

"${PAPER3_CURATION_PYTHON_BIN}" "${PYTHON_SCRIPT}" \
  --run-dir "${RUN_DIR}" \
  --destination "${AUTHOR_OBJECT}"

module purge
module load Bioinformatics
module load "${PAPER3_CURATION_R_MODULE}"
Rscript "${R_SCRIPT}" \
  --object "${AUTHOR_OBJECT}" \
  --ca301 "${CA301_OBJECT}" \
  --output-dir "${OUTPUT_DIR}"

rm -f "${RUN_DIR}/BANDLER_AUTHOR_OBJECT_FAILED.txt"
printf 'PASS\ncompleted=%s\njob_id=%s\n' "$(date --iso-8601=seconds)" "${SLURM_JOB_ID:-unknown}" \
  > "${RUN_DIR}/BANDLER_AUTHOR_OBJECT_SUCCESS.txt"
