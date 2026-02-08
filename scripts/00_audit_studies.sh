#!/usr/bin/env bash
#
# Quick consistency/status audit across studies in PROJECT_ROOT.
# Prints a small Markdown table you can paste into WORKFLOW.md/README.md.

set -euo pipefail

usage() {
  cat <<'EOF'
Usage: 00_audit_studies.sh --project-root /nfs/turbo/.../mgeo_neuron_scrnaseq_projectfolder

Environment:
  PROJECT_ROOT   Alternative way to provide --project-root.
EOF
}

project_root="${PROJECT_ROOT:-}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    -p|--project-root)
      project_root="$2"
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
  usage
  exit 1
fi
project_root="${project_root%/}"

yn_file() { [[ -f "$1" ]] && echo "Yes" || echo "No"; }
yn_dir() { [[ -d "$1" ]] && echo "Yes" || echo "No"; }

note_if_missing() {
  local path="$1"
  [[ -e "$path" ]] || echo "missing: $path"
}

raw_dir="${project_root}/data/raw"
results_dir="${project_root}/results"

echo "| Study | Raw inputs available | Seurat object available | UMAP PNG available | Notes |"
echo "| --- | --- | --- | --- | --- |"

# Walsh (GSE250482)
walsh_raw="${raw_dir}/walsh_2025_geo_files/suppl/GSE250482_RAW.tar"
walsh_seu="${results_dir}/walsh_day75/walsh_day75_final.rds"
walsh_umap="${results_dir}/walsh_day75/plots/umap_by_cluster.png"
echo "| Walsh (GSE250482) | $(yn_file "$walsh_raw") | $(yn_file "$walsh_seu") | $(yn_file "$walsh_umap") | run label: walsh_day75 |"

# Bershteyn 2025 (GSE283775) provided Seurat object
bers25_raw="${raw_dir}/bershteyn_2025_geo_files/suppl/GSE283775_Seurat_scRNA_seq.rds.gz"
bers25_seu="${results_dir}/bershteyn_2025/bershteyn_2025_seurat.rds"
bers25_umap="${results_dir}/bershteyn_2025/plots/umap_by_cluster.png"
echo "| Bershteyn 2025 (GSE283775) | $(yn_file "$bers25_raw") | $(yn_file "$bers25_seu") | $(yn_file "$bers25_umap") | canonical path is results/bershteyn_2025 (legacy results/bershteyn kept for compatibility) |"

# Xiang 2018 (GSE98201 10x trio)
xiang_raw="${raw_dir}/xiang_2018_geo_files/suppl/GSE98201_matrix.mtx.gz"
xiang_seu="${results_dir}/xiang_2018/xiang_2018_seurat.rds"
xiang_umap="${results_dir}/xiang_2018/plots/umap_by_cluster.png"
echo "| Xiang 2018 (GSE98201) | $(yn_file "$xiang_raw") | $(yn_file "$xiang_seu") | $(yn_file "$xiang_umap") | built from 10x trio in suppl/ |"

# Bershteyn 2023 (GSE208672) provided Seurat object
bers23_raw="${raw_dir}/bershteyn_2023_geo_files/suppl/GSE208672_Seurat_allsamples.rds.gz"
bers23_seu="${results_dir}/bershteyn_2023/bershteyn_2023_seurat.rds"
bers23_umap="${results_dir}/bershteyn_2023/plots/umap_by_cluster.png"
echo "| Bershteyn 2023 (GSE208672) | $(yn_file "$bers23_raw") | $(yn_file "$bers23_seu") | $(yn_file "$bers23_umap") | canonical copy saved in results/bershteyn_2023 |"

# Samarasinghe 2021 (GSE165577) counts
samar_raw="${raw_dir}/samarasinghe_2021_geo_files/suppl/GSE165577_Filtered_counts_all_samples.csv.gz"
samar_seu="${results_dir}/samarasinghe_2021/samarasinghe_2021_seurat.rds"
samar_umap="${results_dir}/samarasinghe_2021/plots/umap_by_cluster.png"
echo "| Samarasinghe 2021 (GSE165577) | $(yn_file "$samar_raw") | $(yn_file "$samar_seu") | $(yn_file "$samar_umap") | seurat+umap not generated yet |"

# Siebert 2026 (NeMO placeholder)
siebert_meta="${raw_dir}/siebert_2026_nemo/metadata.json"
echo "| Siebert 2026 (NeMO nemo:dat-htzat9t) | $(yn_file "$siebert_meta") (metadata only) | No | No | waiting on NeMO manifests/files |"

echo
echo "Sanity checks:"
echo "- data/raw exists: $(yn_dir "$raw_dir")"
echo "- results exists: $(yn_dir "$results_dir")"
