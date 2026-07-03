#!/usr/bin/env bash
# Trim existing OME-TIFF stacks to a 1-based inclusive Z range without rerunning
# stitching, median filtering, or 8-bit conversion.
#
# Inputs:
#   INPUTS_FILE: text file containing one input OME-TIFF path per line.
#   Z_START/Z_END: 1-based inclusive slice range. For the current requested
#     BC43/realbive4 trim, use 25 and 430.
#
# Output:
#   OUTPUT_DIR: one trimmed OME-TIFF stack per input, preserving dtype.

#SBATCH --job-name=ome-ztrim
#SBATCH --account=parent0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#SBATCH --output=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/ome-ztrim-%j.out
#SBATCH --error=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/ome-ztrim-%j.err

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline}"
PROJECT_ROOT="${PROJECT_ROOT:-/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder}"
CONDA_BASE="${CONDA_BASE:-/home/elcrespo/miniconda3}"
CONDA_ENV="${CONDA_ENV:-imaris-careamics}"
INPUTS_FILE="${INPUTS_FILE:?set INPUTS_FILE to newline-delimited OME-TIFF paths}"
OUTPUT_DIR="${OUTPUT_DIR:?set OUTPUT_DIR for trimmed outputs}"
Z_START="${Z_START:-25}"
Z_END="${Z_END:-430}"
SUFFIX="${SUFFIX:-_z025_z430}"

mkdir -p "${PROJECT_ROOT}/logs" "${OUTPUT_DIR}"

echo "OME-TIFF Z trim started: $(date)"
echo "Host: $(hostname)"
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-}"
echo "INPUTS_FILE=${INPUTS_FILE}"
echo "OUTPUT_DIR=${OUTPUT_DIR}"
echo "Z_START=${Z_START}"
echo "Z_END=${Z_END}"
echo "SUFFIX=${SUFFIX}"

source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

mapfile -t INPUTS < "${INPUTS_FILE}"
python "${REPO_ROOT}/imaging/imaris_careamics/slice_ome_tiff_z_range.py" \
  --output-dir "${OUTPUT_DIR}" \
  --z-start "${Z_START}" \
  --z-end "${Z_END}" \
  --suffix "${SUFFIX}" \
  "${INPUTS[@]}"

echo "OME-TIFF Z trim finished: $(date)"
find "${OUTPUT_DIR}" -maxdepth 1 -type f -name '*.ome.tif' -printf '%p %s bytes\n' | sort
