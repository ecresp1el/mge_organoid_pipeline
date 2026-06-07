#!/usr/bin/env bash

# Download GEO files for GSE286235, the hnbMO/nucleus basalis organoid
# scRNA-seq study. The processed supplement is a TAR of Cell Ranger H5 files.

set -euo pipefail

usage() {
  cat <<'EOF'
Usage: 01f_download_gse286235_geo.sh --project-root /nfs/turbo/.../mgeo_neuron_scrnaseq_projectfolder

Required arguments:
  -p, --project-root   Runtime workspace. Data are written under
                       PROJECT_ROOT/data/raw/liu_2025_hnbmo_geo_files/.

Environment:
  PROJECT_ROOT         Alternative way to provide --project-root.

Behavior:
  * Downloads GSE286235 series metadata and processed Cell Ranger H5 archive.
  * Skips downloads when destination files already exist.
  * Verifies compressed/archive files.
  * Extracts H5 files from GSE286235_RAW.tar into the supplementary directory.
  * Writes download_manifest.tsv with URL, destination, size, and sha256.
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

  local base="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE286nnn/GSE286235"
  local target_dir="${project_root%/}/data/raw/liu_2025_hnbmo_geo_files"
  local manifest="${target_dir}/download_manifest.tsv"
  local suppl_dir="${target_dir}/suppl"
  local raw_tar="${suppl_dir}/GSE286235_RAW.tar"

  mkdir -p "$target_dir" "$suppl_dir"
  echo -e "url\tdestination\tsize_bytes\tsha256" > "$manifest"

  download_and_record "${base}/matrix/GSE286235_series_matrix.txt.gz" \
    "${target_dir}/matrix/GSE286235_series_matrix.txt.gz" "$manifest"
  download_and_record "${base}/miniml/GSE286235_family.xml.tgz" \
    "${target_dir}/miniml/GSE286235_family.xml.tgz" "$manifest"
  download_and_record "${base}/soft/GSE286235_family.soft.gz" \
    "${target_dir}/soft/GSE286235_family.soft.gz" "$manifest"
  download_and_record "${base}/suppl/filelist.txt" \
    "${suppl_dir}/filelist.txt" "$manifest"
  download_and_record "${base}/suppl/GSE286235_RAW.tar" \
    "$raw_tar" "$manifest"

  if ! compgen -G "${suppl_dir}/GSM872144*_raw_feature_bc_matrix.h5" >/dev/null; then
    echo "Extracting H5 files from: $raw_tar" >&2
    tar -xf "$raw_tar" -C "$suppl_dir"
  else
    echo "Found extracted H5 files, skipping extraction." >&2
  fi

  echo "Done. Files are under: $target_dir"
}

main "$@"
