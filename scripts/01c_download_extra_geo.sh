#!/usr/bin/env bash

# Download GEO accessions:
#   * GSE97882  (Xiang_2018 organoid scRNA-seq)
#   * GSE165577 (Samarasinghe_2021 scRNA-seq)
# Follows the same pattern as 01_download_geo.sh, writing manifests and verifying archives.

set -euo pipefail

usage() {
  cat <<'EOF'
Usage: 01c_download_extra_geo.sh --project-root /nfs/turbo/.../mgeo_neuron_scrnaseq_projectfolder

Required arguments:
  -p, --project-root   Path to runtime workspace (PROJECT_ROOT). Data will be
                       written under PROJECT_ROOT/data/raw/.

Environment:
  PROJECT_ROOT         Alternative way to provide --project-root.

Behavior:
  * Downloads GEO series files for GSE97882 (Xiang_2018) and GSE165577 (Samarasinghe_2021).
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
  bases["GSE97882"]="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE97nnn/GSE97882"
  bases["GSE165577"]="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165577"

  # dataset_dir_name|subdir|filename|url_suffix
  local gse97882_files=(
    "xiang_2018_geo_files|matrix|GSE97882-GPL11154_series_matrix.txt.gz|/matrix/GSE97882-GPL11154_series_matrix.txt.gz"
    "xiang_2018_geo_files|matrix|GSE97882-GPL20301_series_matrix.txt.gz|/matrix/GSE97882-GPL20301_series_matrix.txt.gz"
    "xiang_2018_geo_files|miniml|GSE97882_family.xml.tgz|/miniml/GSE97882_family.xml.tgz"
    "xiang_2018_geo_files|soft|GSE97882_family.soft.gz|/soft/GSE97882_family.soft.gz"
    "xiang_2018_geo_files|suppl|filelist.txt|/suppl/filelist.txt"
    "xiang_2018_geo_files|suppl|GSE97882_RAW.tar|/suppl/GSE97882_RAW.tar"
  )

  local gse165577_files=(
    "samarasinghe_2021_geo_files|matrix|GSE165577_series_matrix.txt.gz|/matrix/GSE165577_series_matrix.txt.gz"
    "samarasinghe_2021_geo_files|miniml|GSE165577_family.xml.tgz|/miniml/GSE165577_family.xml.tgz"
    "samarasinghe_2021_geo_files|soft|GSE165577_family.soft.gz|/soft/GSE165577_family.soft.gz"
    "samarasinghe_2021_geo_files|suppl|GSE165577_Filtered_counts_all_samples.csv.gz|/suppl/GSE165577_Filtered_counts_all_samples.csv.gz"
    "samarasinghe_2021_geo_files|suppl|GSE165577_Normalized_counts_all_samples.csv.gz|/suppl/GSE165577_Normalized_counts_all_samples.csv.gz"
  )

  process_dataset() {
    local accession="$1"
    local dataset_dir_name="$2"
    local -n file_entries="$3"

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

  process_dataset "GSE97882" "xiang_2018_geo_files" gse97882_files
  process_dataset "GSE165577" "samarasinghe_2021_geo_files" gse165577_files

  echo "Done."
}

main "$@"
