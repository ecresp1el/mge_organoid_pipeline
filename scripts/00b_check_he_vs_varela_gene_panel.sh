#!/usr/bin/env bash

# Check completion of additive He-vs-Varela gene-panel pipeline stages.

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: 00b_check_he_vs_varela_gene_panel.sh --project-root /nfs/turbo/.../mgeo_neuron_scrnaseq_projectfolder [--run-label gene_panel_he_vs_varela_v1]
USAGE
}

project_root="${PROJECT_ROOT:-}"
run_label="gene_panel_he_vs_varela_v1"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -p|--project-root)
      project_root="$2"
      shift 2
      ;;
    -r|--run-label)
      run_label="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage
      exit 1
      ;;
  esac
done

if [[ -z "${project_root}" ]]; then
  echo "Error: --project-root or PROJECT_ROOT is required" >&2
  exit 1
fi

project_root="${project_root%/}"

yn() { [[ -f "$1" ]] && echo "Yes" || echo "No"; }

raw_h5ad="${project_root}/data/raw/he_et_al_zenodo/suppl/hnoca_allmeta.h5ad"
slice_counts="${project_root}/data/processed/he_et_al_scn8a_slice/scn8a_counts.mtx"
slice_genes="${project_root}/data/processed/he_et_al_scn8a_slice/genes.tsv"
slice_umap="${project_root}/data/processed/he_et_al_scn8a_slice/umap.tsv.gz"
he_seurat="${project_root}/results/he_et_al/he_et_al_scn8a_seurat.rds"
panel_png="${project_root}/results/${run_label}/plots/gene_panel_cross_study_log.png"
panel_status="${project_root}/results/${run_label}/plots/gene_panel_cross_study_status.tsv"

echo "| Stage | Artifact | Exists |"
echo "| --- | --- | --- |"
echo "| Download full V2 h5ad | ${raw_h5ad} | $(yn "${raw_h5ad}") |"
echo "| Extract slice counts | ${slice_counts} | $(yn "${slice_counts}") |"
echo "| Extract slice genes | ${slice_genes} | $(yn "${slice_genes}") |"
echo "| Extract slice UMAP | ${slice_umap} | $(yn "${slice_umap}") |"
echo "| Build He Seurat | ${he_seurat} | $(yn "${he_seurat}") |"
echo "| Gene panel image | ${panel_png} | $(yn "${panel_png}") |"
echo "| Gene panel status table | ${panel_status} | $(yn "${panel_status}") |"
