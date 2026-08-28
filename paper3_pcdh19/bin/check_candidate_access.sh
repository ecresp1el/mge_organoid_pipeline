#!/usr/bin/env bash
# PURPOSE
#   Step 00, read-only source discovery. Report whether each path registered in
#   config/input_candidates.tsv exists and can be read/traversed.
#
# INPUTS
#   --manifest PATH       Candidate registry (default: input_candidates.tsv).
#   --max-depth INTEGER   Maximum directory depth for metadata aggregation
#                        (default: 3; zero is allowed).
#
# COMPUTATION
#   For a directory that can be traversed, find counts regular files and sums
#   their stat-reported byte sizes through the requested depth. For a regular
#   file, stat supplies its byte size. Permission tests use shell -r/-x.
#
# OUTPUTS AND SIDE EFFECTS
#   Writes one TSV to standard output with existence/access/type/count/size
#   fields. Diagnostics and usage go to standard error. It creates, copies,
#   hashes, opens for content inspection, and modifies no scientific file.
#
# SCIENTIFIC SCOPE
#   This is a storage/access inventory only. It cannot validate matrix content,
#   sample identity, tissue, condition, genotype, replicate structure, or QC.
#   The full contract is in ../PIPELINE_IO_AND_BIOLOGICAL_SCOPE.md.
set -Eeuo pipefail

BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
MANIFEST="${BUNDLE_DIR}/config/input_candidates.tsv"
MAX_DEPTH=3

usage() {
  # Print the command-line contract to stderr; this function changes no state.
  cat >&2 <<'EOF'
Usage: check_candidate_access.sh [--manifest PATH] [--max-depth INTEGER]

Print a metadata-only TSV access inventory to standard output. The command
does not open, checksum, copy, or modify scientific files.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest)
      [[ $# -ge 2 ]] || { usage; exit 2; }
      MANIFEST="$2"
      shift 2
      ;;
    --max-depth)
      [[ $# -ge 2 ]] || { usage; exit 2; }
      MAX_DEPTH="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 2
      ;;
  esac
done

[[ -f "${MANIFEST}" ]] || { echo "Missing manifest: ${MANIFEST}" >&2; exit 2; }
[[ "${MAX_DEPTH}" =~ ^[0-9]+$ ]] || {
  echo "--max-depth must be a non-negative integer; got: ${MAX_DEPTH}" >&2
  exit 2
}

printf 'candidate_id\texperiment_id\texists\treadable\ttraversable\tobject_type\tfiles_through_depth_%s\tbytes_through_depth_%s\tsource_root\n' \
  "${MAX_DEPTH}" "${MAX_DEPTH}"

{
  IFS= read -r _header
  while IFS=$'\t' read -r candidate_id experiment_id _display_name source_root _rest; do
    [[ -n "${candidate_id}" ]] || continue

    exists=false
    readable=false
    traversable=not_applicable
    object_type=missing
    file_count=not_scanned
    total_bytes=not_scanned

    if [[ -e "${source_root}" ]]; then
      exists=true
      [[ -r "${source_root}" ]] && readable=true

      if [[ -d "${source_root}" ]]; then
        object_type=directory
        traversable=false
        if [[ -x "${source_root}" ]]; then
          traversable=true
          file_count="$(find "${source_root}" -maxdepth "${MAX_DEPTH}" -type f -printf '.\n' 2>/dev/null | wc -l)"
          total_bytes="$(find "${source_root}" -maxdepth "${MAX_DEPTH}" -type f -printf '%s\n' 2>/dev/null | awk '{sum += $1} END {printf "%.0f", sum + 0}')"
        fi
      elif [[ -f "${source_root}" ]]; then
        object_type=file
        file_count=1
        total_bytes="$(stat -c '%s' "${source_root}")"
      else
        object_type=other
      fi
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${candidate_id}" "${experiment_id}" "${exists}" "${readable}" \
      "${traversable}" "${object_type}" "${file_count}" "${total_bytes}" \
      "${source_root}"
  done
} < "${MANIFEST}"
