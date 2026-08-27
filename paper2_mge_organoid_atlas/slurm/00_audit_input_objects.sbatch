#!/usr/bin/env bash
# Resources and log paths are assigned by bin/submit.sh.

set -Eeuo pipefail

: "${PAPER2_RUN_DIR:?submit through paper2_mge_organoid_atlas/bin/submit.sh}"
: "${PAPER2_INPUT_MANIFEST:?input manifest is required}"

RUN_DIR="${PAPER2_RUN_DIR}"
SCRIPT="${RUN_DIR}/code/00_audit_input_objects.R"

failed() {
  local line="$1"
  local status="$2"
  printf 'FAILED\njob_id=%s\nline=%s\nexit_status=%s\ntime=%s\n' \
    "${SLURM_JOB_ID:-unknown}" "${line}" "${status}" "$(date --iso-8601=seconds)" \
    > "${RUN_DIR}/FAILED.txt"
}
trap 'failed "${LINENO}" "$?"' ERR

module purge
module load Bioinformatics
module load "${SEURAT5_MODULE:-r-seurat/5.1.0-R-4.4.1-c3m7yfq}"

{
  echo "started=$(date --iso-8601=seconds)"
  echo "host=$(hostname)"
  echo "slurm_job_id=${SLURM_JOB_ID:-unknown}"
  echo "loaded_modules=${LOADEDMODULES:-}"
  echo "R=$(Rscript --version 2>&1 | head -n 1)"
} > "${RUN_DIR}/provenance/environment.txt"

echo "Paper 2 six-study input audit started: $(date --iso-8601=seconds)"
echo "Run directory: ${RUN_DIR}"
echo "Input registry: ${PAPER2_INPUT_MANIFEST}"

Rscript "${SCRIPT}" \
  --manifest "${PAPER2_INPUT_MANIFEST}" \
  --outdir "${RUN_DIR}"

for required in \
  "${RUN_DIR}/tables/input_objects_audited.tsv" \
  "${RUN_DIR}/tables/object_summary.tsv" \
  "${RUN_DIR}/tables/assay_layer_inventory.tsv" \
  "${RUN_DIR}/tables/metadata_column_inventory.tsv" \
  "${RUN_DIR}/tables/reduction_inventory.tsv" \
  "${RUN_DIR}/tables/feature_overlap_pairwise.tsv" \
  "${RUN_DIR}/tables/feature_overlap_summary.tsv" \
  "${RUN_DIR}/provenance/sessionInfo.txt"; do
  [[ -s "${required}" ]] || { echo "Required output missing or empty: ${required}" >&2; exit 1; }
done

(
  cd "${RUN_DIR}"
  find README.md code config figures tables -type f -print0 \
    | sort -z \
    | xargs -0 sha256sum
) > "${RUN_DIR}/provenance/package_sha256_manifest.txt"
[[ -s "${RUN_DIR}/provenance/package_sha256_manifest.txt" ]]

printf 'PASS\njob_id=%s\ncompleted=%s\nstep=00_input_audit\n' \
  "${SLURM_JOB_ID:-unknown}" "$(date --iso-8601=seconds)" \
  > "${RUN_DIR}/SUCCESS.txt"
rm -f "${RUN_DIR}/FAILED.txt"

echo "Paper 2 six-study input audit completed: $(date --iso-8601=seconds)"
