#!/usr/bin/env bash
# Apply Fiji/ImageJ median radius 1 to corrected stitched planes, then package
# filtered 16-bit and 8-bit outputs as single-series OME-TIFF stacks.

#SBATCH --job-name=fiji-median-realbive4
#SBATCH --account=parent0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/fiji-median-realbive4-%j.out
#SBATCH --error=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/fiji-median-realbive4-%j.err

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline}"
PROJECT_ROOT="${PROJECT_ROOT:-/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder}"
CONDA_BASE="${CONDA_BASE:-/home/elcrespo/miniconda3}"
CONDA_ENV="${CONDA_ENV:-imaris-careamics}"
INPUT_DIR="${INPUT_DIR:-${PROJECT_ROOT}/results/imaging/fiji_stitching/cl32_bive4_pv_reporter_40x_realbive4_corrected_stage_full}"
OUTPUT_ROOT="${OUTPUT_ROOT:-${INPUT_DIR}/median_r1}"
RADIUS="${RADIUS:-1.0}"
JAVA_BIN="${JAVA_BIN:-java}"
JAVA_HEAP="${JAVA_HEAP:-24g}"
CLASS_DIR="${CLASS_DIR:-${TMPDIR:-/tmp}/${USER:-elcrespo}/fiji_median_classes_${SLURM_JOB_ID:-manual}}"

PLANES_16="${OUTPUT_ROOT}/planes_16bit"
PLANES_8="${OUTPUT_ROOT}/planes_8bit"
STACKS_16="${OUTPUT_ROOT}/tiff_stacks_16bit"
STACKS_8="${OUTPUT_ROOT}/tiff_stacks_8bit"

mkdir -p "${PROJECT_ROOT}/logs" "${PLANES_16}" "${PLANES_8}" "${STACKS_16}" "${STACKS_8}" "${OUTPUT_ROOT}/logs"
exec > >(tee -a "${OUTPUT_ROOT}/logs/run_fiji_median_filter_${SLURM_JOB_ID:-manual}.log") 2>&1

echo "Fiji median filtering started: $(date)"
echo "Host: $(hostname)"
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-}"
echo "INPUT_DIR=${INPUT_DIR}"
echo "OUTPUT_ROOT=${OUTPUT_ROOT}"
echo "RADIUS=${RADIUS}"

printf 'input channel 1 planes: '
find "${INPUT_DIR}" -maxdepth 1 -type f -name 'img_*_c1' | wc -l
printf 'input channel 2 planes: '
find "${INPUT_DIR}" -maxdepth 1 -type f -name 'img_*_c2' | wc -l

mkdir -p "${CLASS_DIR}"
FIJI_CP="$(find /sw/pkgs/med/fiji/1.5.4/jars /sw/pkgs/med/fiji/1.5.4/plugins -name '*.jar' | paste -sd: -)"
javac -cp "${FIJI_CP}" \
  -d "${CLASS_DIR}" \
  "${REPO_ROOT}/imaging/imaris_careamics/java/plugin/FijiMedianFilterPlanes.java"

"${JAVA_BIN}" \
  -Djava.awt.headless=true \
  -Xmx"${JAVA_HEAP}" \
  -cp "${CLASS_DIR}:${FIJI_CP}" \
  plugin.FijiMedianFilterPlanes \
  "${INPUT_DIR}" \
  "${PLANES_16}" \
  "${PLANES_8}" \
  "${RADIUS}"

source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

python "${REPO_ROOT}/imaging/imaris_careamics/fiji_planes_to_tiff_stacks.py" \
  --input-dir "${PLANES_16}" \
  --output-dir "${STACKS_16}" \
  --prefix cl32_bive4_pv_reporter_40x_realbive4_corrected_stage_median_r1_16bit \
  --channels 1 2 \
  --dtype uint16

python "${REPO_ROOT}/imaging/imaris_careamics/fiji_planes_to_tiff_stacks.py" \
  --input-dir "${PLANES_8}" \
  --output-dir "${STACKS_8}" \
  --prefix cl32_bive4_pv_reporter_40x_realbive4_corrected_stage_median_r1_8bit \
  --channels 1 2 \
  --dtype uint8

echo "Fiji median filtering finished: $(date)"
find "${OUTPUT_ROOT}" -maxdepth 3 -type f \( -name '*.ome.tif' -o -name '*.log' \) -printf '%p %s bytes\n' | sort
