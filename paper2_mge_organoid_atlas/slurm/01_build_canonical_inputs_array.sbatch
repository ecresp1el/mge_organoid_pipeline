#!/usr/bin/env bash

set -Eeuo pipefail

: "${PAPER2_CANONICAL_STAGING:?canonical staging directory is required}"
: "${PAPER2_CANONICAL_MANIFEST:?canonical manifest is required}"
: "${CONDA_ENV_BIN:?Python environment bin directory is required}"
: "${SLURM_ARRAY_TASK_ID:?array task ID is required}"

STUDY_ID="$(awk -F '\t' -v row="$((SLURM_ARRAY_TASK_ID + 2))" 'NR == row {print $1}' "${PAPER2_CANONICAL_MANIFEST}")"
[[ -n "${STUDY_ID}" ]] || { echo "No study for array task ${SLURM_ARRAY_TASK_ID}" >&2; exit 2; }
STUDY_DIR="${PAPER2_CANONICAL_STAGING}/${STUDY_ID}"
R_SCRIPT="${PAPER2_CANONICAL_STAGING}/code/01_build_canonical_rds.R"
PY_SCRIPT="${PAPER2_CANONICAL_STAGING}/code/01_write_validate_canonical_h5ad.py"
mkdir -p "${STUDY_DIR}"

failed() {
  local line="$1"
  local status="$2"
  printf 'FAILED\nstudy_id=%s\njob_id=%s\narray_task_id=%s\nline=%s\nexit_status=%s\ntime=%s\n' \
    "${STUDY_ID}" "${SLURM_JOB_ID:-unknown}" "${SLURM_ARRAY_TASK_ID}" "${line}" "${status}" "$(date --iso-8601=seconds)" \
    > "${STUDY_DIR}/FAILED.txt"
}
trap 'failed "${LINENO}" "$?"' ERR

module purge
module load Bioinformatics
module load "${SEURAT5_MODULE:-r-seurat/5.1.0-R-4.4.1-c3m7yfq}"

{
  echo "started=$(date --iso-8601=seconds)"
  echo "host=$(hostname)"
  echo "study_id=${STUDY_ID}"
  echo "slurm_job_id=${SLURM_JOB_ID:-unknown}"
  echo "slurm_array_task_id=${SLURM_ARRAY_TASK_ID}"
  echo "loaded_modules=${LOADEDMODULES:-}"
  echo "R=$(Rscript --version 2>&1 | head -n 1)"
  echo "Python=$(${CONDA_ENV_BIN}/python --version 2>&1)"
} > "${STUDY_DIR}/environment.txt"

echo "Canonical input build started for ${STUDY_ID}: $(date --iso-8601=seconds)"
Rscript "${R_SCRIPT}" \
  --manifest "${PAPER2_CANONICAL_MANIFEST}" \
  --study-id "${STUDY_ID}" \
  --outdir "${STUDY_DIR}"

"${CONDA_ENV_BIN}/python" "${PY_SCRIPT}" \
  --study-id "${STUDY_ID}" \
  --study-dir "${STUDY_DIR}"

[[ -s "${STUDY_DIR}/${STUDY_ID}_minimal.rds" ]]
[[ -s "${STUDY_DIR}/${STUDY_ID}_minimal.h5ad" ]]
[[ -s "${STUDY_DIR}/SUCCESS.txt" ]]
[[ ! -e "${STUDY_DIR}/.bridge" ]]
rm -f "${STUDY_DIR}/FAILED.txt"
echo "Canonical input build completed for ${STUDY_ID}: $(date --iso-8601=seconds)"
