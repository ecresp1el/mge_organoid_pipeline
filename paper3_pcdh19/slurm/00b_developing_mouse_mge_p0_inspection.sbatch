#!/usr/bin/env bash
# Scheduler resources, array bounds, dependencies, and logs are assigned by the submission wrapper.
set -Eeuo pipefail

: "${PAPER3_CURATION_RUN_DIR:?submit through bin/submit_developing_mouse_mge_reference_curation.sh}"
: "${PAPER3_CURATION_SOURCE_ROOT:?source cache root is required}"
: "${PAPER3_CURATION_REGISTRY:?submitted source registry is required}"
: "${PAPER3_CURATION_PYTHON_BIN:?Python interpreter is required}"
: "${PAPER3_CURATION_R_MODULE:?R/Seurat module is required}"
: "${SLURM_ARRAY_TASK_ID:?this job must run as the submitted three-task array}"

RUN_DIR="${PAPER3_CURATION_RUN_DIR}"
SCRIPT="${RUN_DIR}/code/developing_mouse_mge_reference_curation.py"
R_HELPER="${RUN_DIR}/code/inspect_bandler_rds.R"
case "${SLURM_ARRAY_TASK_ID}" in
  0) STUDY=LaManno2021 ;;
  1) STUDY=Bandler2022 ;;
  2) STUDY=Mayer2018 ;;
  *) echo "Unexpected array task: ${SLURM_ARRAY_TASK_ID}" >&2; exit 2 ;;
esac

failed() {
  local line="$1" status="$2"
  printf 'FAILED\nstage=p0_inspection\nstudy=%s\njob_id=%s\narray_task=%s\nline=%s\nexit_status=%s\ntime=%s\n' \
    "${STUDY}" "${SLURM_JOB_ID:-unknown}" "${SLURM_ARRAY_TASK_ID}" "${line}" "${status}" "$(date --iso-8601=seconds)" \
    > "${RUN_DIR}/provenance/p0_${STUDY}_FAILED.txt"
  printf 'FAILED\nstage=p0_inspection\nstudy=%s\njob_id=%s\ntime=%s\n' \
    "${STUDY}" "${SLURM_JOB_ID:-unknown}" "$(date --iso-8601=seconds)" > "${RUN_DIR}/FAILED.txt"
}
trap 'failed "${LINENO}" "$?"' ERR

module purge
module load Bioinformatics
module load "${PAPER3_CURATION_R_MODULE}"

"${PAPER3_CURATION_PYTHON_BIN}" -c 'import anndata, h5py, matplotlib, pandas, scipy'
Rscript --version

{
  echo "started=$(date --iso-8601=seconds)"
  echo "host=$(hostname)"
  echo "slurm_job_id=${SLURM_JOB_ID:-unknown}"
  echo "array_task=${SLURM_ARRAY_TASK_ID}"
  echo "study=${STUDY}"
  echo "python=$(${PAPER3_CURATION_PYTHON_BIN} --version 2>&1)"
  echo "R=$(Rscript --version 2>&1 | head -n 1)"
  echo "loaded_modules=${LOADEDMODULES:-}"
} > "${RUN_DIR}/provenance/p0_${STUDY}_environment.txt"

"${PAPER3_CURATION_PYTHON_BIN}" "${SCRIPT}" inspect-study \
  --study "${STUDY}" \
  --registry "${PAPER3_CURATION_REGISTRY}" \
  --run-dir "${RUN_DIR}" \
  --source-root "${PAPER3_CURATION_SOURCE_ROOT}" \
  --r-helper "${R_HELPER}"

rm -f "${RUN_DIR}/provenance/p0_${STUDY}_FAILED.txt"
