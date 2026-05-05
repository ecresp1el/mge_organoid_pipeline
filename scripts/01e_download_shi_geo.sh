#!/usr/bin/env bash

# Download GEO accession GSE135827 (Shi et al, 2019) into PROJECT_ROOT.
# Idempotent: existing files are not re-downloaded, but integrity is re-checked and
# a fresh manifest with checksums is written on each run.

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: 01e_download_shi_geo.sh --project-root /nfs/turbo/.../mgeo_neuron_scrnaseq_projectfolder [options]

Required arguments:
  -p, --project-root     Path to runtime workspace (PROJECT_ROOT). Data will be
                         written under PROJECT_ROOT/data/raw/shi_2019_geo_files/.

Optional arguments:
  --skip-raw-tar         Skip downloading GSE135827_RAW.tar.
  -h, --help             Show this help and exit.

Environment:
  PROJECT_ROOT           Alternative way to provide --project-root.

Behavior:
  * Downloads GEO series files for GSE135827.
  * Skips downloads when destination files already exist.
  * Verifies archives/compressed files (tar -tf, gunzip -t) even when skipped.
  * Writes a manifest TSV with URL, destination, size, md5, sha256.
USAGE
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
    *) ;;
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

  local size_bytes md5_value sha256sum_value
  size_bytes=$(stat -c%s "$dest")
  md5_value=$(md5sum "$dest" | awk '{print $1}')
  sha256sum_value=$(sha256sum "$dest" | awk '{print $1}')

  printf "%s\t%s\t%s\t%s\t%s\n" \
    "$url" "$dest" "$size_bytes" "$md5_value" "$sha256sum_value" >> "$manifest"
}

main() {
  local project_root="${PROJECT_ROOT:-}"
  local skip_raw_tar="false"

  while [[ $# -gt 0 ]]; do
    case "$1" in
      -p|--project-root)
        project_root="$2"
        shift 2
        ;;
      --skip-raw-tar)
        skip_raw_tar="true"
        shift
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
  require_cmd md5sum
  require_cmd sha256sum

  local base="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135827"
  local dest_root="${project_root%/}/data/raw/shi_2019_geo_files"
  local manifest="${dest_root}/download_manifest.tsv"

  mkdir -p "$dest_root"
  echo -e "url\tdestination\tsize_bytes\tmd5\tsha256" > "$manifest"

  # matrix/miniml/soft descriptors
  local files=(
    "matrix|GSE135827_series_matrix.txt.gz|/matrix/GSE135827_series_matrix.txt.gz"
    "miniml|GSE135827_family.xml.tgz|/miniml/GSE135827_family.xml.tgz"
    "soft|GSE135827_family.soft.gz|/soft/GSE135827_family.soft.gz"
    "suppl|filelist.txt|/suppl/filelist.txt"
    "suppl|GSE135827_GE_mat_raw_count_with_week_info.txt.gz|/suppl/GSE135827_GE_mat_raw_count_with_week_info.txt.gz"
  )

  if [[ "$skip_raw_tar" != "true" ]]; then
    files+=("suppl|GSE135827_RAW.tar|/suppl/GSE135827_RAW.tar")
  fi

  local subdir filename suffix dest url
  for entry in "${files[@]}"; do
    IFS='|' read -r subdir filename suffix <<< "$entry"
    dest="${dest_root}/${subdir}/${filename}"
    url="${base}${suffix}"
    download_and_record "$url" "$dest" "$manifest"
  done

  echo "Done."
}

main "$@"
