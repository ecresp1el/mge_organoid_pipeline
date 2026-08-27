#!/usr/bin/env bash

set -Eeuo pipefail

: "${PAPER2_GENE_RUN_DIR:?run directory is required}"
: "${CONDA_ENV_BIN:?Python environment bin directory is required}"

failed() {
  local line="$1"
  local status="$2"
  printf 'FAILED\njob_id=%s\nline=%s\nexit_status=%s\ntime=%s\n' \
    "${SLURM_JOB_ID:-unknown}" "${line}" "${status}" "$(date --iso-8601=seconds)" \
    > "${PAPER2_GENE_RUN_DIR}/FAILED.txt"
}
trap 'failed "${LINENO}" "$?"' ERR

{
  echo "started=$(date --iso-8601=seconds)"
  echo "host=$(hostname)"
  echo "slurm_job_id=${SLURM_JOB_ID:-unknown}"
  echo "Python=$(${CONDA_ENV_BIN}/python --version 2>&1)"
} > "${PAPER2_GENE_RUN_DIR}/provenance/environment.txt"

echo "Step 02 gene identifier audit started: $(date --iso-8601=seconds)"
"${CONDA_ENV_BIN}/python" \
  "${PAPER2_GENE_RUN_DIR}/code/02_audit_harmonize_genes.py" \
  --run-dir "${PAPER2_GENE_RUN_DIR}"

[[ -s "${PAPER2_GENE_RUN_DIR}/SUCCESS.txt" ]]
[[ -s "${PAPER2_GENE_RUN_DIR}/tables/six_way_gene_intersection_strict.tsv" ]]
[[ -s "${PAPER2_GENE_RUN_DIR}/tables/feature_mapping_long.tsv.gz" ]]
rm -f "${PAPER2_GENE_RUN_DIR}/FAILED.txt"
echo "Step 02 gene identifier audit completed: $(date --iso-8601=seconds)"
