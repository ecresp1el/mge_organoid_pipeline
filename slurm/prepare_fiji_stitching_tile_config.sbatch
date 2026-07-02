#!/usr/bin/env bash
# Generate Fiji Grid/Collection Stitching TileConfiguration.txt from a
# TeraStitcher-style XML file.
#
# Submit on Great Lakes:
#   INPUT_DIR=/nfs/turbo/umms-parent/andor_micropscope_data_dump/exp17_pv_reporter_with_biver3and4/20x \
#   XML_NAME=cl32_bive4_pv_reporter_40x_realbive4.xml \
#     sbatch slurm/prepare_fiji_stitching_tile_config.sbatch

#SBATCH --job-name=fiji-tilecfg
#SBATCH --account=parent0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/fiji-tilecfg-%j.out
#SBATCH --error=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/fiji-tilecfg-%j.err

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline}"
PROJECT_ROOT="${PROJECT_ROOT:-/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder}"
INPUT_DIR="${INPUT_DIR:-/nfs/turbo/umms-parent/andor_micropscope_data_dump/exp17_pv_reporter_with_biver3and4/20x}"
XML_NAME="${XML_NAME:-cl32_bive4_pv_reporter_40x_realbive4.xml}"
OUTPUT_NAME="${OUTPUT_NAME:-TileConfiguration.txt}"
PYTHON_BIN="${PYTHON_BIN:-python3}"

mkdir -p "${PROJECT_ROOT}/logs"

echo "Fiji tile configuration prep started: $(date)"
echo "Host: $(hostname)"
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-}"
echo "REPO_ROOT=${REPO_ROOT}"
echo "INPUT_DIR=${INPUT_DIR}"
echo "XML_NAME=${XML_NAME}"
echo "OUTPUT_NAME=${OUTPUT_NAME}"

"${PYTHON_BIN}" "${REPO_ROOT}/imaging/imaris_careamics/terastitcher_to_fiji_tile_config.py" \
  --xml "${INPUT_DIR}/${XML_NAME}" \
  --output "${INPUT_DIR}/${OUTPUT_NAME}"

echo
echo "Generated layout:"
sed -n '1,120p' "${INPUT_DIR}/${OUTPUT_NAME}"

echo
echo "Open Fiji from a Great Lakes desktop session with:"
echo "  module load fiji/1.5.4"
echo "  /sw/pkgs/med/fiji/1.5.4/ImageJ-linux64 &"
echo
echo "Use Plugins -> Stitching -> Grid/Collection Stitching:"
echo "  Type: Positions from file"
echo "  Layout file: ${INPUT_DIR}/${OUTPUT_NAME}"
echo "  Fusion method: Linear Blending"
echo "  Compute overlap: checked"
echo "  Regression threshold: 0.30"
echo "  Max/avg displacement threshold: 2.50"
echo "  Absolute displacement threshold: 3.50"
echo "  Image output: Fuse and display"
echo
echo "Do not use FusionStitcher output, do not split channels first, and trust XML tile order."
echo "Fiji tile configuration prep finished: $(date)"
