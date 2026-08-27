#!/usr/bin/env bash

set -Eeuo pipefail

: "${PAPER2_LEGACY_RUN_DIR:?run directory is required}"
: "${CONDA_ENV_BIN:?Python environment bin directory is required}"

failed() {
  local line="$1"
  local status="$2"
  printf 'FAILED\njob_id=%s\nline=%s\nexit_status=%s\ntime=%s\n' \
    "${SLURM_JOB_ID:-unknown}" "${line}" "${status}" "$(date --iso-8601=seconds)" \
    > "${PAPER2_LEGACY_RUN_DIR}/FAILED.txt"
}
trap 'failed "${LINENO}" "$?"' ERR

{
  echo "started=$(date --iso-8601=seconds)"
  echo "host=$(hostname)"
  echo "slurm_job_id=${SLURM_JOB_ID:-unknown}"
  echo "Python=$(${CONDA_ENV_BIN}/python --version 2>&1)"
} > "${PAPER2_LEGACY_RUN_DIR}/provenance/environment.txt"

echo "Step 02b legacy gene-ID recovery audit started: $(date --iso-8601=seconds)"
"${CONDA_ENV_BIN}/python" \
  "${PAPER2_LEGACY_RUN_DIR}/code/02b_recover_legacy_gene_ids.py" \
  --run-dir "${PAPER2_LEGACY_RUN_DIR}"

[[ -s "${PAPER2_LEGACY_RUN_DIR}/SUCCESS.txt" ]]
[[ -s "${PAPER2_LEGACY_RUN_DIR}/tables/recovery_candidate_mapping_long.tsv.gz" ]]
[[ -s "${PAPER2_LEGACY_RUN_DIR}/tables/recovery_outcome_by_dataset.tsv" ]]
[[ -s "${PAPER2_LEGACY_RUN_DIR}/tables/proposed_six_way_intersection_summary.tsv" ]]
rm -f "${PAPER2_LEGACY_RUN_DIR}/FAILED.txt"
echo "Step 02b legacy gene-ID recovery audit completed: $(date --iso-8601=seconds)"
