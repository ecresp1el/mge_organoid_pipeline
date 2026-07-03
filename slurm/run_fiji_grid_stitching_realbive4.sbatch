#!/usr/bin/env bash
# Stitch the BC43/realbive4 2x2 Imaris montage with Fiji Grid/Collection
# Stitching through the dialog-free Java runner.
#
# Inputs:
#   INPUT_DIR: directory containing the .ims tiles and the selected
#     TileConfiguration file.
#   LAYOUT_NAME: Fiji TileConfiguration file inside INPUT_DIR. For the corrected
#     BC43/realbive4 full stitch, use:
#       TileConfiguration.corrected_stage_mip_qc.txt
#   GENERATE_LAYOUT: true regenerates LAYOUT_NAME from XML; false uses an
#     existing reviewed layout. For corrected production output, use false.
#
# Output:
#   OUTPUT_DIR: Fiji plane-per-file fused output named img_t*_z*_c*.
#   This is an intermediate format. Run stack_fiji_stitched_planes.sbatch after
#   this if you need two single-series OME-TIFF stacks.
#
# Important:
#   This does not use FusionStitcher output and does not split channels before
#   stitching. Do not blindly use XML ABS_H/ABS_V for future montages without
#   checking MIP QC first.
#
# Submit:
#   GENERATE_LAYOUT=false \
#   LAYOUT_NAME=TileConfiguration.corrected_stage_mip_qc.txt \
#   OUTPUT_LABEL=cl32_bive4_pv_reporter_40x_realbive4_corrected_stage_full \
#   sbatch slurm/run_fiji_grid_stitching_realbive4.sbatch

#SBATCH --job-name=fiji-stitch-realbive4
#SBATCH --account=parent0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/fiji-stitch-realbive4-%j.out
#SBATCH --error=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/fiji-stitch-realbive4-%j.err

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline}"
PROJECT_ROOT="${PROJECT_ROOT:-/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder}"
INPUT_DIR="${INPUT_DIR:-/nfs/turbo/umms-parent/andor_micropscope_data_dump/exp17_pv_reporter_with_biver3and4/20x}"
XML_NAME="${XML_NAME:-cl32_bive4_pv_reporter_40x_realbive4.xml}"
LAYOUT_NAME="${LAYOUT_NAME:-TileConfiguration.txt}"
GENERATE_LAYOUT="${GENERATE_LAYOUT:-true}"
COMPUTE_OVERLAP="${COMPUTE_OVERLAP:-false}"
OUTPUT_LABEL="${OUTPUT_LABEL:-cl32_bive4_pv_reporter_40x_realbive4_xml_coords_trusted}"
OUTPUT_DIR="${OUTPUT_DIR:-${PROJECT_ROOT}/results/imaging/fiji_stitching/${OUTPUT_LABEL}}"
JAVA_BIN="${JAVA_BIN:-java}"
JAVA_HEAP="${JAVA_HEAP:-48g}"
CLASS_DIR="${CLASS_DIR:-${TMPDIR:-/tmp}/${USER:-elcrespo}/fiji_stitcher_classes_${SLURM_JOB_ID:-manual}}"

mkdir -p "${PROJECT_ROOT}/logs" "${OUTPUT_DIR}/logs"
exec > >(tee -a "${OUTPUT_DIR}/logs/run_fiji_grid_stitching_${SLURM_JOB_ID:-manual}.log") 2>&1

echo "Fiji Grid/Collection Stitching started: $(date)"
echo "Host: $(hostname)"
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-}"
echo "REPO_ROOT=${REPO_ROOT}"
echo "INPUT_DIR=${INPUT_DIR}"
echo "XML_NAME=${XML_NAME}"
echo "LAYOUT_NAME=${LAYOUT_NAME}"
echo "GENERATE_LAYOUT=${GENERATE_LAYOUT}"
echo "COMPUTE_OVERLAP=${COMPUTE_OVERLAP}"
echo "OUTPUT_DIR=${OUTPUT_DIR}"
echo "JAVA_BIN=${JAVA_BIN}"
echo "JAVA_HEAP=${JAVA_HEAP}"
echo "CLASS_DIR=${CLASS_DIR}"

if [[ "${GENERATE_LAYOUT}" == "true" ]]; then
  # Literal XML conversion. This is useful for quick prep, but final production
  # layouts should usually come from a reviewed MIP QC/audit step.
  python3 "${REPO_ROOT}/imaging/imaris_careamics/terastitcher_to_fiji_tile_config.py" \
    --xml "${INPUT_DIR}/${XML_NAME}" \
    --output "${INPUT_DIR}/${LAYOUT_NAME}"
else
  echo "Using existing tile configuration: ${INPUT_DIR}/${LAYOUT_NAME}"
fi

echo
echo "Tile configuration:"
sed -n '1,120p' "${INPUT_DIR}/${LAYOUT_NAME}"

# Compile the Java runner inside the job so the Slurm node uses the same Fiji
# jars that will execute the stitch.
mkdir -p "${CLASS_DIR}"
FIJI_CP="$(find /sw/pkgs/med/fiji/1.5.4/jars /sw/pkgs/med/fiji/1.5.4/plugins -name '*.jar' | paste -sd: -)"
javac -cp "${FIJI_CP}" \
  -d "${CLASS_DIR}" \
  "${REPO_ROOT}/imaging/imaris_careamics/java/plugin/FijiGridCollectionStitcher.java"

# Runs Fiji Grid/Collection Stitching directly. The last argument disables
# virtual input so Bio-Formats opens each IMS normally.
"${JAVA_BIN}" \
  -Djava.awt.headless=true \
  -Xmx"${JAVA_HEAP}" \
  -cp "${CLASS_DIR}:${FIJI_CP}" \
  plugin.FijiGridCollectionStitcher \
  "${INPUT_DIR}" \
  "${LAYOUT_NAME}" \
  "${OUTPUT_DIR}" \
  "${COMPUTE_OVERLAP}" \
  false

echo
echo "Fiji Grid/Collection Stitching finished: $(date)"
echo "Output inventory:"
find "${OUTPUT_DIR}" -maxdepth 3 -type f | sort
