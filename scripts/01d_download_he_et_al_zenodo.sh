#!/usr/bin/env bash

# Download He et al HNOCA full object (Version 2) from Zenodo into PROJECT_ROOT.
# Idempotent: existing files are not re-downloaded; checksum is always verified.

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: 01d_download_he_et_al_zenodo.sh --project-root /nfs/turbo/.../mgeo_neuron_scrnaseq_projectfolder [options]

Required arguments:
  -p, --project-root   Path to runtime workspace (PROJECT_ROOT). Data will be
                       written under PROJECT_ROOT/data/raw/he_et_al_zenodo/.

Optional arguments:
  --url URL            Direct download URL (default: Zenodo record 14160929, hnoca_allmeta.h5ad)
  --filename NAME      Output filename (default: hnoca_allmeta.h5ad)
  --expected-md5 HEX   Expected MD5 checksum (default: 2fb8b272b61646f314c0f933e4040d22)
  -h, --help           Show this help and exit

Environment:
  PROJECT_ROOT         Alternative way to provide --project-root.

Behavior:
  * Downloads the h5ad file under data/raw/he_et_al_zenodo/suppl/.
  * Skips download if file already exists.
  * Verifies md5 checksum and writes a download manifest TSV.
USAGE
}

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "Error: required command '$1' not found in PATH" >&2
    exit 1
  }
}

download_and_record() {
  local url="$1"
  local dest="$2"
  local expected_md5="$3"
  local manifest="$4"

  mkdir -p "$(dirname "$dest")"

  if [[ -f "$dest" ]]; then
    echo "Found existing file, skipping download: $dest" >&2
  else
    echo "Downloading $url -> $dest" >&2
    local tmp="${dest}.partial"
    rm -f "$tmp"
    curl -L --fail --retry 4 --retry-delay 5 -o "$tmp" "$url"
    mv "$tmp" "$dest"
  fi

  local file_md5
  file_md5=$(md5sum "$dest" | awk '{print $1}')
  if [[ "$file_md5" != "$expected_md5" ]]; then
    echo "Error: md5 mismatch for $dest" >&2
    echo "  expected: $expected_md5" >&2
    echo "  actual:   $file_md5" >&2
    exit 1
  fi

  local size_bytes sha256sum_value
  size_bytes=$(stat -c%s "$dest")
  sha256sum_value=$(sha256sum "$dest" | awk '{print $1}')

  printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$url" "$dest" "$size_bytes" "$file_md5" "$sha256sum_value" "$expected_md5" >> "$manifest"
}

main() {
  local project_root="${PROJECT_ROOT:-}"
  local url="https://zenodo.org/records/14160929/files/hnoca_allmeta.h5ad?download=1"
  local filename="hnoca_allmeta.h5ad"
  local expected_md5="2fb8b272b61646f314c0f933e4040d22"

  while [[ $# -gt 0 ]]; do
    case "$1" in
      -p|--project-root)
        project_root="$2"
        shift 2
        ;;
      --url)
        url="$2"
        shift 2
        ;;
      --filename)
        filename="$2"
        shift 2
        ;;
      --expected-md5)
        expected_md5="$2"
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
  require_cmd md5sum
  require_cmd sha256sum

  local dest_root="${project_root%/}/data/raw/he_et_al_zenodo"
  local dest="${dest_root}/suppl/${filename}"
  local manifest="${dest_root}/download_manifest.tsv"

  mkdir -p "$dest_root"
  echo -e "url\tdestination\tsize_bytes\tmd5\tsha256\texpected_md5" > "$manifest"

  download_and_record "$url" "$dest" "$expected_md5" "$manifest"

  echo "Done."
}

main "$@"
