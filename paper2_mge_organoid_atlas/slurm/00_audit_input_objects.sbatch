#!/usr/bin/env bash
# Resources and log paths are assigned by bin/submit.sh.

set -Eeuo pipefail

: "${PAPER2_RUN_DIR:?submit through paper2_mge_organoid_atlas/bin/submit.sh}"
: "${PAPER2_INPUT_MANIFEST:?input manifest is required}"
: "${PAPER2_CLUSTER_OVERRIDES:?cluster override table is required}"
: "${CONDA_ENV_BIN:?Python environment bin directory is required}"

RUN_DIR="${PAPER2_RUN_DIR}"
SCRIPT="${RUN_DIR}/code/00_audit_input_objects.R"
PLOT_SCRIPT="${RUN_DIR}/code/00_render_input_umap_inventory.py"

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
  --cluster-overrides "${PAPER2_CLUSTER_OVERRIDES}" \
  --outdir "${RUN_DIR}"

"${CONDA_ENV_BIN}/python" "${PLOT_SCRIPT}" \
  --umap-dir "${RUN_DIR}/tables/umap" \
  --outdir "${RUN_DIR}/figures" \
  --manifest-outdir "${RUN_DIR}/tables" \
  --png-dpi "${PNG_DPI:-300}" \
  --pdf-dpi "${PDF_DPI:-300}" \
  --make-svg "${MAKE_SVG:-false}" \
  --svg-dpi "${SVG_DPI:-300}"

for required in \
  "${RUN_DIR}/tables/input_objects_audited.tsv" \
  "${RUN_DIR}/tables/object_summary.tsv" \
  "${RUN_DIR}/tables/assay_layer_inventory.tsv" \
  "${RUN_DIR}/tables/metadata_column_inventory.tsv" \
  "${RUN_DIR}/tables/reduction_inventory.tsv" \
  "${RUN_DIR}/tables/umap_source_inventory.tsv" \
  "${RUN_DIR}/tables/umap_cluster_counts.tsv" \
  "${RUN_DIR}/tables/umap_figure_manifest.tsv" \
  "${RUN_DIR}/tables/umap_figure_panel_summary.tsv" \
  "${RUN_DIR}/tables/feature_overlap_pairwise.tsv" \
  "${RUN_DIR}/tables/feature_overlap_summary.tsv" \
  "${RUN_DIR}/provenance/sessionInfo.txt"; do
  [[ -s "${required}" ]] || { echo "Required output missing or empty: ${required}" >&2; exit 1; }
done

for figure_basename in \
  six_study_input_umap_cluster_inventory \
  varela_div30_div90_input_umap_cluster_inventory \
  varela_div30_input_umap_cluster_inventory \
  varela_div90_input_umap_cluster_inventory; do
  for extension in png pdf; do
    figure_path="${RUN_DIR}/figures/${extension}/${figure_basename}.${extension}"
    [[ -s "${figure_path}" ]] || { echo "Required UMAP figure missing or empty: ${figure_path}" >&2; exit 1; }
  done
  if [[ "${MAKE_SVG:-false}" == "true" ]]; then
    figure_path="${RUN_DIR}/figures/svg/${figure_basename}.svg"
    [[ -s "${figure_path}" ]] || { echo "Required UMAP figure missing or empty: ${figure_path}" >&2; exit 1; }
  fi
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
