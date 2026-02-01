#!/usr/bin/env bash

# Download GEO accessions GSE283775 and GSE250482 into the PROJECT_ROOT workspace.
# Idempotent: existing files are not re-downloaded, but integrity is re-checked and
# a fresh manifest with checksums is written on each run.

set -euo pipefail

usage() {
  cat <<'EOF'
Usage: 01_download_geo.sh --project-root /nfs/turbo/.../mgeo_neuron_scrnaseq_projectfolder

Required arguments:
  -p, --project-root   Path to runtime workspace (PROJECT_ROOT). Data will be
                       written under PROJECT_ROOT/data/raw/.

Environment:
  PROJECT_ROOT         Alternative way to provide --project-root.

Behavior:
  * Downloads GEO series files for GSE283775 (Bershteyn) and GSE250482 (Walsh).
  * Skips downloads when destination files already exist.
  * Verifies archives/compressed files (tar -tf, gunzip -t) even when skipped.
  * Writes per-dataset manifest TSVs with URL, destination, size, sha256.
EOF
}

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "Error: required command '$1' not found in PATH" >&2
    exit 1
  }
}

verify_file() {
  local path="$1"
  case "$path" in
    *.tar) tar -tf "$path" >/dev/null ;;
    *.tgz|*.gz) gunzip -t "$path" >/dev/null ;;
    *) ;; # nothing to do
  esac
}

download_and_record() {
  local url="$1"
  local dest="$2"
  local manifest="$3"
  local tmp

  mkdir -p "$(dirname "$dest")"

  if [[ -f "$dest" ]]; then
    echo "Found existing file, skipping download: $dest" >&2
  else
    echo "Downloading $url -> $dest" >&2
    tmp="${dest}.partial"
    rm -f "$tmp"
    curl -L --fail --retry 4 --retry-delay 5 -o "$tmp" "$url"
    mv "$tmp" "$dest"
  fi

  verify_file "$dest"

  local size_bytes sha256sum_value
  size_bytes=$(stat -c%s "$dest")
  sha256sum_value=$(sha256sum "$dest" | awk '{print $1}')
  printf "%s\t%s\t%s\t%s\n" "$url" "$dest" "$size_bytes" "$sha256sum_value" >> "$manifest"
}

main() {
  local project_root="${PROJECT_ROOT:-}"

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

  if [[ -z "$project_root" ]]; then
    echo "Error: --project-root or PROJECT_ROOT is required" >&2
    usage
    exit 1
  fi

  require_cmd curl
  require_cmd tar
  require_cmd gunzip
  require_cmd sha256sum

  local dest_root="${project_root%/}/data/raw"

  declare -A bases
  bases["GSE283775"]="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE283nnn/GSE283775"
  bases["GSE250482"]="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE250nnn/GSE250482"

  # Dataset definitions: dataset_dir_name|subdir|filename|url_suffix
  local gse283775_files=(
    "bershteyn_2025_geo_files|matrix|GSE283775_series_matrix.txt.gz|/matrix/GSE283775_series_matrix.txt.gz"
    "bershteyn_2025_geo_files|miniml|GSE283775_family.xml.tgz|/miniml/GSE283775_family.xml.tgz"
    "bershteyn_2025_geo_files|soft|GSE283775_family.soft.gz|/soft/GSE283775_family.soft.gz"
    "bershteyn_2025_geo_files|suppl|GSE283775_Seurat_scRNA_seq.rds.gz|/suppl/GSE283775_Seurat_scRNA_seq.rds.gz"
  )

  local gse250482_files=(
    "walsh_2025_geo_files|matrix|GSE250482_series_matrix.txt.gz|/matrix/GSE250482_series_matrix.txt.gz"
    "walsh_2025_geo_files|miniml|GSE250482_family.xml.tgz|/miniml/GSE250482_family.xml.tgz"
    "walsh_2025_geo_files|soft|GSE250482_family.soft.gz|/soft/GSE250482_family.soft.gz"
    "walsh_2025_geo_files|suppl|filelist.txt|/suppl/filelist.txt"
    "walsh_2025_geo_files|suppl|GSE250482_RAW.tar|/suppl/GSE250482_RAW.tar"
  )

  process_dataset() {
    local accession="$1"
    local dataset_dir_name="$2"
    local -n file_entries="$3" # nameref

    local base="${bases[$accession]}"
    local target_dir="${dest_root}/${dataset_dir_name}"
    local manifest="${target_dir}/download_manifest.tsv"

    mkdir -p "$target_dir"
    echo -e "url\tdestination\tsize_bytes\tsha256" > "$manifest"

    for entry in "${file_entries[@]}"; do
      IFS='|' read -r _ subdir filename url_suffix <<< "$entry"
      local dest="${target_dir}/${subdir}/${filename}"
      local url="${base}${url_suffix}"
      download_and_record "$url" "$dest" "$manifest"
    done
  }

  process_dataset "GSE283775" "bershteyn_2025_geo_files" gse283775_files
  process_dataset "GSE250482" "walsh_2025_geo_files" gse250482_files

  echo "Done."
}

main "$@"
