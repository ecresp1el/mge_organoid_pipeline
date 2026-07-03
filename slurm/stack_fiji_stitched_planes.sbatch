#!/usr/bin/env bash
# Convert Fiji Grid/Collection plane-per-file output into one 16-bit BigTIFF
# stack per channel.

#SBATCH --job-name=fiji-stack-tiff
#SBATCH --account=parent0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/fiji-stack-tiff-%j.out
#SBATCH --error=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/fiji-stack-tiff-%j.err

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline}"
PROJECT_ROOT="${PROJECT_ROOT:-/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder}"
CONDA_BASE="${CONDA_BASE:-/home/elcrespo/miniconda3}"
CONDA_ENV="${CONDA_ENV:-imaris-careamics}"
INPUT_DIR="${INPUT_DIR:?set INPUT_DIR to the Fiji plane output directory}"
OUTPUT_DIR="${OUTPUT_DIR:-${INPUT_DIR}/tiff_stacks}"
PREFIX="${PREFIX:-cl32_bive4_pv_reporter_40x_realbive4_corrected_stage}"

mkdir -p "${PROJECT_ROOT}/logs" "${OUTPUT_DIR}"

echo "Stack conversion started: $(date)"
echo "Host: $(hostname)"
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-}"
echo "INPUT_DIR=${INPUT_DIR}"
echo "OUTPUT_DIR=${OUTPUT_DIR}"
echo "PREFIX=${PREFIX}"

source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

python "${REPO_ROOT}/imaging/imaris_careamics/fiji_planes_to_tiff_stacks.py" \
  --input-dir "${INPUT_DIR}" \
  --output-dir "${OUTPUT_DIR}" \
  --prefix "${PREFIX}" \
  --channels 1 2

echo "Stack conversion finished: $(date)"
find "${OUTPUT_DIR}" -maxdepth 1 -type f -name '*.tif' -printf '%p %s bytes\n' | sort
