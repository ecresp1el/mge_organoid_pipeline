#!/bin/bash

# Download healthy Wang/Liu 2025 hnbMO SRA runs and convert them to
# Cell Ranger-compatible FASTQ names. This reconstructs upstream FASTQs for
# GSE286235/PRJNA1206345 from SRA; GEO supplementary files only contain raw
# Cell Ranger matrices.

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  scripts/01g_gse286235_sra_to_10x_fastqs.sh --sample-id <sample_id> [options]

Required:
  --sample-id <sample_id>        BF_H9_D36, BF_H9_D63, or BFCO_IMR_D63

Options:
  --project-root <path>          Runtime project root
  --repo-root <path>             Repository root containing config manifest
  --manifest <path>              SRA run manifest TSV
  --threads <int>                fasterq-dump/pigz threads
  --force                        Recreate existing FASTQ outputs

Environment:
  PROJECT_ROOT                   Alternative project root
  REPO_ROOT                      Alternative repo root
  THREADS                        Alternative thread count
EOF
}

PROJECT_ROOT="${PROJECT_ROOT:-/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder}"
REPO_ROOT="${REPO_ROOT:-/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline}"
MANIFEST=""
SAMPLE_ID=""
THREADS="${THREADS:-8}"
FORCE="false"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample-id)
      SAMPLE_ID="$2"
      shift 2
      ;;
    --project-root|-p)
      PROJECT_ROOT="$2"
      shift 2
      ;;
    --repo-root)
      REPO_ROOT="$2"
      shift 2
      ;;
    --manifest)
      MANIFEST="$2"
      shift 2
      ;;
    --threads)
      THREADS="$2"
      shift 2
      ;;
    --force)
      FORCE="true"
      shift
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ -z "${SAMPLE_ID}" ]]; then
  echo "Missing --sample-id" >&2
  usage >&2
  exit 2
fi

if [[ -z "${MANIFEST}" ]]; then
  MANIFEST="${REPO_ROOT}/config/liu_2025_hnbmo_healthy_sra_runs.tsv"
fi

if [[ ! -f "${MANIFEST}" ]]; then
  echo "Manifest not found: ${MANIFEST}" >&2
  exit 1
fi

RAW_ROOT="${PROJECT_ROOT}/data/raw/liu_2025_hnbmo_sra"
SRA_DIR="${RAW_ROOT}/sra"
FASTQ_ROOT="${RAW_ROOT}/fastqs"
FASTQ_DIR="${FASTQ_ROOT}/${SAMPLE_ID}"
WORK_ROOT="${RAW_ROOT}/work/${SAMPLE_ID}"
TMP_ROOT="${WORK_ROOT}/tmp"

mkdir -p "${SRA_DIR}" "${FASTQ_DIR}" "${WORK_ROOT}" "${TMP_ROOT}"

mapfile -t RUN_ROWS < <(awk -v sample="${SAMPLE_ID}" 'BEGIN{FS="\t"} NR > 1 && $1 == sample {print $0}' "${MANIFEST}")

if [[ "${#RUN_ROWS[@]}" -eq 0 ]]; then
  echo "No runs found for sample ${SAMPLE_ID} in ${MANIFEST}" >&2
  exit 1
fi

echo "Sample: ${SAMPLE_ID}"
echo "Runs: ${#RUN_ROWS[@]}"
echo "SRA dir: ${SRA_DIR}"
echo "FASTQ dir: ${FASTQ_DIR}"
echo "Threads: ${THREADS}"

find_sra_file() {
  local run="$1"
  local direct="${SRA_DIR}/${run}/${run}.sra"
  if [[ -f "${direct}" ]]; then
    echo "${direct}"
    return 0
  fi
  local found
  found="$(find "${SRA_DIR}" -path "*/${run}*.sra" -type f 2>/dev/null | head -n 1 || true)"
  if [[ -n "${found}" ]]; then
    echo "${found}"
    return 0
  fi
  return 1
}

for row in "${RUN_ROWS[@]}"; do
  IFS=$'\t' read -r sample_id figure1c_sample sra_run lane sra_sample_name library_name size_mb <<<"${row}"

  out_r1="${FASTQ_DIR}/${SAMPLE_ID}_S1_${lane}_R1_001.fastq.gz"
  out_r2="${FASTQ_DIR}/${SAMPLE_ID}_S1_${lane}_R2_001.fastq.gz"
  if [[ "${FORCE}" != "true" && -s "${out_r1}" && -s "${out_r2}" ]]; then
    echo "FASTQs already exist for ${sra_run} (${lane}); skipping"
    continue
  fi

  if [[ "${FORCE}" == "true" ]]; then
    rm -f "${out_r1}" "${out_r2}"
  fi

  echo "Downloading ${sra_run} (${figure1c_sample}, ${lane}, expected SRA MB=${size_mb})"
  if ! find_sra_file "${sra_run}" >/dev/null 2>&1; then
    prefetch --output-directory "${SRA_DIR}" --max-size 250G "${sra_run}"
  else
    echo "SRA file already present for ${sra_run}"
  fi

  sra_file="$(find_sra_file "${sra_run}")"
  echo "Converting ${sra_run}: ${sra_file}"

  run_tmp="${WORK_ROOT}/${sra_run}_fastq_tmp"
  rm -rf "${run_tmp}"
  mkdir -p "${run_tmp}" "${TMP_ROOT}/${sra_run}"

  fasterq-dump \
    --split-files \
    --include-technical \
    --threads "${THREADS}" \
    --temp "${TMP_ROOT}/${sra_run}" \
    --outdir "${run_tmp}" \
    "${sra_file}"

  in_r1="${run_tmp}/${sra_run}_1.fastq"
  in_r2="${run_tmp}/${sra_run}_2.fastq"
  if [[ ! -s "${in_r1}" || ! -s "${in_r2}" ]]; then
    echo "Expected paired FASTQs were not produced for ${sra_run}; saw:" >&2
    find "${run_tmp}" -maxdepth 1 -type f -print >&2
    exit 1
  fi

  echo "Compressing/renaming ${sra_run} to ${SAMPLE_ID} ${lane}"
  pigz -p "${THREADS}" -c "${in_r1}" > "${out_r1}.tmp"
  pigz -p "${THREADS}" -c "${in_r2}" > "${out_r2}.tmp"
  mv "${out_r1}.tmp" "${out_r1}"
  mv "${out_r2}.tmp" "${out_r2}"
  rm -rf "${run_tmp}" "${TMP_ROOT:?}/${sra_run}"
done

ls -lh "${FASTQ_DIR}"
echo "FASTQ conversion complete for ${SAMPLE_ID}"
