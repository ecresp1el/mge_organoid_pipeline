#!/bin/bash

# Run Cell Ranger count on reconstructed GSE286235 healthy FASTQs.
# Uses the available Great Lakes Cell Ranger module rather than the paper's
# reported Cell Ranger v3.1, because v3.1 is not exposed as a module.

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  scripts/05m_liu_2025_hnbmo_cellranger_count.sh --sample-id <sample_id> [options]

Required:
  --sample-id <sample_id>        BF_H9_D36, BF_H9_D63, or BFCO_IMR_D63

Options:
  --project-root <path>          Runtime project root
  --fastq-root <path>            Parent directory containing per-sample FASTQ dirs
  --out-root <path>              Cell Ranger output root
  --transcriptome <path>         10x-compatible transcriptome reference
  --localcores <int>             Cell Ranger local cores
  --localmem <int>               Cell Ranger local memory in GB
  --force                        Delete existing sample output before rerun

Environment:
  PROJECT_ROOT                   Alternative project root
  CELLRANGER_REF                 Alternative transcriptome reference
  LOCALCORES                     Alternative localcores
  LOCALMEM                       Alternative localmem
EOF
}

PROJECT_ROOT="${PROJECT_ROOT:-/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder}"
SAMPLE_ID=""
FASTQ_ROOT=""
OUT_ROOT=""
TRANSCRIPTOME="${CELLRANGER_REF:-/nfs/turbo/umms-parent/Manny_human_ref/refdata-gex-GRCh38-2020-A}"
LOCALCORES="${LOCALCORES:-16}"
LOCALMEM="${LOCALMEM:-120}"
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
    --fastq-root)
      FASTQ_ROOT="$2"
      shift 2
      ;;
    --out-root)
      OUT_ROOT="$2"
      shift 2
      ;;
    --transcriptome)
      TRANSCRIPTOME="$2"
      shift 2
      ;;
    --localcores)
      LOCALCORES="$2"
      shift 2
      ;;
    --localmem)
      LOCALMEM="$2"
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

if [[ -z "${FASTQ_ROOT}" ]]; then
  FASTQ_ROOT="${PROJECT_ROOT}/data/raw/liu_2025_hnbmo_sra/fastqs"
fi
if [[ -z "${OUT_ROOT}" ]]; then
  OUT_ROOT="${PROJECT_ROOT}/results/liu_2025_hnbmo_cellranger_counts"
fi

FASTQ_DIR="${FASTQ_ROOT}/${SAMPLE_ID}"
OUT_DIR="${OUT_ROOT}/${SAMPLE_ID}"

if [[ ! -d "${FASTQ_DIR}" ]]; then
  echo "FASTQ directory not found: ${FASTQ_DIR}" >&2
  exit 1
fi
if [[ ! -d "${TRANSCRIPTOME}" ]]; then
  echo "Transcriptome reference not found: ${TRANSCRIPTOME}" >&2
  exit 1
fi

mkdir -p "${OUT_ROOT}"

if [[ "${FORCE}" == "true" && -d "${OUT_DIR}" ]]; then
  rm -rf "${OUT_DIR}"
fi

if [[ -s "${OUT_DIR}/outs/filtered_feature_bc_matrix.h5" ]]; then
  echo "Cell Ranger output already exists for ${SAMPLE_ID}: ${OUT_DIR}"
  exit 0
fi

echo "Sample: ${SAMPLE_ID}"
echo "FASTQs: ${FASTQ_DIR}"
echo "Output root: ${OUT_ROOT}"
echo "Transcriptome: ${TRANSCRIPTOME}"
echo "Cell Ranger: $(cellranger --version)"
echo "localcores=${LOCALCORES}, localmem=${LOCALMEM}"
ls -lh "${FASTQ_DIR}"

cd "${OUT_ROOT}"

cellranger count \
  --id "${SAMPLE_ID}" \
  --transcriptome "${TRANSCRIPTOME}" \
  --fastqs "${FASTQ_DIR}" \
  --sample "${SAMPLE_ID}" \
  --chemistry SC3Pv3 \
  --localcores "${LOCALCORES}" \
  --localmem "${LOCALMEM}" \
  --no-bam \
  --disable-ui

echo "Cell Ranger complete for ${SAMPLE_ID}: ${OUT_DIR}"
