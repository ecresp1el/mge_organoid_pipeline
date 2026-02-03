#!/usr/bin/env bash

# Download NeMO collection nemo:dat-htzat9t (Siebert 2026) landing-page assets
# following https://nemoarchive.org/resources/data-download#data-collection-landing-pages.
# Steps:
#   1) Fetch collection metadata JSON from assets.nemoarchive.org/api/collection/<id>.
#   2) Download any manifest_file_urls (typically BDBag .tgz archives).
#   3) Extract BDBags and filter fetch.txt entries for MGE organoid / Seurat content.
#   4) (Optional, default) download the filtered subset files to materialize raw/processed
#      organoid data and Seurat objects.
#
# Output layout under PROJECT_ROOT:
#   data/raw/siebert_2026_nemo/
#       metadata.json                 (API response)
#       manifest_download.tsv         (bag downloads + checksums)
#       bags/<bag>.tgz                (downloaded BDBags)
#       bags/<bag>/fetch.txt          (from bag)
#       subset_fetch.tsv              (filtered lines from fetch.txt)
#       subset_download_manifest.tsv  (if subset downloads run)
#       data/...                      (downloaded subset files follow fetch paths)

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: 01b_download_siebert_nemo.sh --project-root /nfs/.../mgeo_neuron_scrnaseq_projectfolder [options]

Required arguments:
  -p, --project-root   Path to runtime workspace (PROJECT_ROOT). Data lands under
                       PROJECT_ROOT/data/raw/siebert_2026_nemo/

Optional arguments:
  -i, --nemo-id ID     NeMO collection identifier (default: nemo:dat-htzat9t)
  -n, --dataset-name S Folder name under data/raw/ (default: siebert_2026_nemo)
  --no-download-subset Only build bag + filtered lists; skip downloading subset files
  -h, --help           Show this help and exit

Environment:
  PROJECT_ROOT         Alternative way to provide --project-root

Dependencies: curl, jq, tar, gunzip, sha256sum
USAGE
}

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || { echo "Error: required command '$1' not found" >&2; exit 1; }
}

verify_file() {
  local path="$1"
  case "$path" in
    *.tar|*.tgz) tar -tf "$path" >/dev/null ;;
    *.gz) gunzip -t "$path" >/dev/null ;;
    *) ;; # nothing
  esac
}

download_and_record() {
  local url="$1" dest="$2" manifest="$3"
  local tmp="${dest}.partial"
  mkdir -p "$(dirname "$dest")"

  if [[ -f "$dest" ]]; then
    echo "Found existing file, skipping download: $dest" >&2
  else
    echo "Downloading $url -> $dest" >&2
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
  local project_root="${PROJECT_ROOT:-}" nemo_id="nemo:dat-htzat9t" dataset_name="siebert_2026_nemo" download_subset=1

  while [[ $# -gt 0 ]]; do
    case "$1" in
      -p|--project-root) project_root="$2"; shift 2;;
      -i|--nemo-id) nemo_id="$2"; shift 2;;
      -n|--dataset-name) dataset_name="$2"; shift 2;;
      --no-download-subset) download_subset=0; shift 1;;
      -h|--help) usage; exit 0;;
      *) echo "Unknown option: $1" >&2; usage; exit 1;;
    esac
  done

  if [[ -z "$project_root" ]]; then
    echo "Error: --project-root or PROJECT_ROOT is required" >&2
    usage
    exit 1
  fi

  require_cmd curl
  require_cmd jq
  require_cmd tar
  require_cmd gunzip
  require_cmd sha256sum

  local dest_root="${project_root%/}/data/raw/${dataset_name}"
  mkdir -p "$dest_root"

  local meta_json="$dest_root/metadata.json"
  echo "Fetching metadata for $nemo_id" >&2
  curl -L --fail --silent "https://assets.nemoarchive.org/api/collection/${nemo_id}" | jq '.' > "$meta_json"

  local manifest_tsv="$dest_root/manifest_download.tsv"
  echo -e "url\tdestination\tsize_bytes\tsha256" > "$manifest_tsv"

  local bag_urls
  bag_urls=$(jq -r '.manifest_file_urls[]? | select(.protocol=="bdbag") | .url' "$meta_json")

  if [[ -z "$bag_urls" ]]; then
    echo "No manifest_file_urls (BDBags) published for $nemo_id yet. Nothing to download." >&2
    exit 2
  fi

  for url in $bag_urls; do
    local filename="${url##*/}"
    local dest="$dest_root/bags/$filename"
    download_and_record "$url" "$dest" "$manifest_tsv"

    local bag_dir="${dest%.tgz}"
    if [[ ! -d "$bag_dir" ]]; then
      echo "Extracting $dest" >&2
      mkdir -p "$bag_dir"
      tar -xzf "$dest" -C "$bag_dir/.."
    else
      echo "Bag already extracted: $bag_dir" >&2
    fi

    local fetch_file="$bag_dir/fetch.txt"
    if [[ ! -f "$fetch_file" ]]; then
      echo "Warning: fetch.txt missing in $bag_dir" >&2
      continue
    fi

    local subset_file="$dest_root/subset_fetch.tsv"
    : > "$subset_file"

    # Filter lines for organoid / MGE or Seurat RDS
    grep -Ei '\borganoid\b|\bmge\b|seurat|\.rds' "$fetch_file" >> "$subset_file" || true

    local subset_count
    subset_count=$(wc -l < "$subset_file" | tr -d ' ')
    echo "Filtered subset lines: $subset_count (written to $subset_file)" >&2

    if [[ $download_subset -eq 1 && $subset_count -gt 0 ]]; then
      local subset_manifest="$dest_root/subset_download_manifest.tsv"
      echo -e "url\tdestination\tsize_bytes\tsha256" > "$subset_manifest"

      while IFS=$'\t' read -r file_url file_size rel_path; do
        [[ -z "$file_url" || -z "$rel_path" ]] && continue
        local target="$bag_dir/$rel_path"
        download_and_record "$file_url" "$target" "$subset_manifest"
      done < "$subset_file"
    fi

    # Brief overview summary
    local summary="$dest_root/overview.txt"
    {
      echo "NeMO ID: $nemo_id"
      echo "Bag: $filename"
      echo "Filtered subset lines: $subset_count"
      echo "First 10 subset paths:"
      head -n 10 "$subset_file"
    } > "$summary"
  done

  echo "Done. Outputs in $dest_root" >&2
}

main "$@"
