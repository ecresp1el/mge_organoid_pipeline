#!/usr/bin/env bash
# Convert Fiji/ImageJ plane-per-file output into one BigTIFF OME-TIFF stack per
# channel.
#
# Inputs:
#   INPUT_DIR: plane directory containing files named img_t*_z*_c* or
#     img_t*_z*_c*.tif.
#   DTYPE: expected plane dtype, uint16 or uint8.
#   PREFIX: output filename prefix.
#
# Outputs:
#   OUTPUT_DIR/<PREFIX>_c1_stitched_stack.ome.tif
#   OUTPUT_DIR/<PREFIX>_c2_stitched_stack.ome.tif
#
# Submit examples:
#   INPUT_DIR=.../corrected_stage_full \
#   OUTPUT_DIR=.../corrected_stage_full/tiff_stacks \
#   PREFIX=cl32_bive4_pv_reporter_40x_realbive4_corrected_stage \
#   DTYPE=uint16 \
#   sbatch slurm/stack_fiji_stitched_planes.sbatch
#
#   INPUT_DIR=.../median_r1/planes_8bit \
#   OUTPUT_DIR=.../median_r1/tiff_stacks_8bit \
#   PREFIX=cl32_bive4_pv_reporter_40x_realbive4_corrected_stage_median_r1_8bit \
#   DTYPE=uint8 \
#   sbatch slurm/stack_fiji_stitched_planes.sbatch

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
DTYPE="${DTYPE:-uint16}"

mkdir -p "${PROJECT_ROOT}/logs" "${OUTPUT_DIR}"

echo "Stack conversion started: $(date)"
echo "Host: $(hostname)"
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-}"
echo "INPUT_DIR=${INPUT_DIR}"
echo "OUTPUT_DIR=${OUTPUT_DIR}"
echo "PREFIX=${PREFIX}"
echo "DTYPE=${DTYPE}"

source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

# The Python writer streams planes and writes OME metadata so Fiji/Bio-Formats
# opens one Z stack series instead of one series per plane.
python "${REPO_ROOT}/imaging/imaris_careamics/fiji_planes_to_tiff_stacks.py" \
  --input-dir "${INPUT_DIR}" \
  --output-dir "${OUTPUT_DIR}" \
  --prefix "${PREFIX}" \
  --channels 1 2 \
  --dtype "${DTYPE}"

echo "Stack conversion finished: $(date)"
find "${OUTPUT_DIR}" -maxdepth 1 -type f -name '*.tif' -printf '%p %s bytes\n' | sort
