#!/usr/bin/env bash
# Render the Fiji-saved display TIFFs into a synchronized 1x3 Z-scroll MP4.
#
# Inputs:
#   STACK_DIR: directory containing green_display.tif, magenta_display.tif,
#     and composite_display.tif.
#   FPS: movie frame rate.
#   SCALE: render scale after LUT/range mapping. 1.0 keeps native stitched size.
#   HEADER_HEIGHT/FONT_SIZE/FONT_STROKE_WIDTH: black label bar and thick white
#     panel text settings.
#   SCALE_BAR_UM/XY_UM_PER_PIXEL: optional bottom-right physical scale bar.
#
# Output:
#   OUTPUT: MP4 movie with panels PV-mNG | BiVe3-dTom | Merge.
#   RGB_TIFF: optional rendered RGB TIFF stack of the exact video frames.

#SBATCH --job-name=fiji-zmovie
#SBATCH --account=parent0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --output=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/fiji-zmovie-%j.out
#SBATCH --error=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/fiji-zmovie-%j.err

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline}"
PROJECT_ROOT="${PROJECT_ROOT:-/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder}"
CONDA_BASE="${CONDA_BASE:-/home/elcrespo/miniconda3}"
CONDA_ENV="${CONDA_ENV:-imaris-careamics}"
STACK_DIR="${STACK_DIR:-/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/fiji_stitching/cl32_bive4_pv_reporter_40x_realbive4_corrected_stage_full/z025_z430_stacks}"
OUTPUT="${OUTPUT:-${STACK_DIR}/panel_movie_10fps_bigtext.mp4}"
RGB_TIFF="${RGB_TIFF:-}"
FPS="${FPS:-10}"
BORDER="${BORDER:-10}"
SCALE="${SCALE:-1.0}"
HEADER_HEIGHT="${HEADER_HEIGHT:-180}"
FONT_SIZE="${FONT_SIZE:-120}"
FONT_STROKE_WIDTH="${FONT_STROKE_WIDTH:-5}"
SCALE_BAR_UM="${SCALE_BAR_UM:-0}"
XY_UM_PER_PIXEL="${XY_UM_PER_PIXEL:-}"
SCALE_BAR_LABEL="${SCALE_BAR_LABEL:-}"
SCALE_BAR_MARGIN="${SCALE_BAR_MARGIN:-100}"
SCALE_BAR_HEIGHT="${SCALE_BAR_HEIGHT:-16}"
SCALE_BAR_FONT_SIZE="${SCALE_BAR_FONT_SIZE:-72}"
SCALE_BAR_FONT_STROKE_WIDTH="${SCALE_BAR_FONT_STROKE_WIDTH:-3}"
SCALE_BAR_TEXT_GAP="${SCALE_BAR_TEXT_GAP:-16}"
QUALITY="${QUALITY:-6}"
PRESET="${PRESET:-ultrafast}"
FFMPEG_THREADS="${FFMPEG_THREADS:-4}"

mkdir -p "${PROJECT_ROOT}/logs" "$(dirname "${OUTPUT}")"

echo "Fiji Z-scroll panel movie started: $(date)"
echo "Host: $(hostname)"
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-}"
echo "STACK_DIR=${STACK_DIR}"
echo "OUTPUT=${OUTPUT}"
echo "RGB_TIFF=${RGB_TIFF}"
echo "FPS=${FPS}"
echo "SCALE=${SCALE}"
echo "HEADER_HEIGHT=${HEADER_HEIGHT}"
echo "FONT_SIZE=${FONT_SIZE}"
echo "FONT_STROKE_WIDTH=${FONT_STROKE_WIDTH}"
echo "SCALE_BAR_UM=${SCALE_BAR_UM}"
echo "XY_UM_PER_PIXEL=${XY_UM_PER_PIXEL}"
echo "SCALE_BAR_LABEL=${SCALE_BAR_LABEL}"
echo "SCALE_BAR_MARGIN=${SCALE_BAR_MARGIN}"
echo "SCALE_BAR_HEIGHT=${SCALE_BAR_HEIGHT}"
echo "SCALE_BAR_FONT_SIZE=${SCALE_BAR_FONT_SIZE}"
echo "SCALE_BAR_FONT_STROKE_WIDTH=${SCALE_BAR_FONT_STROKE_WIDTH}"
echo "SCALE_BAR_TEXT_GAP=${SCALE_BAR_TEXT_GAP}"
echo "QUALITY=${QUALITY}"
echo "PRESET=${PRESET}"
echo "FFMPEG_THREADS=${FFMPEG_THREADS}"

rm -f "${OUTPUT}"
if [[ -n "${RGB_TIFF}" ]]; then
  rm -f "${RGB_TIFF}"
fi

source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

CMD=(
  python "${REPO_ROOT}/imaging/imaris_careamics/fiji_z_scroll_panel_movie.py"
  --green "${STACK_DIR}/green_display.tif" \
  --magenta "${STACK_DIR}/magenta_display.tif" \
  --composite "${STACK_DIR}/composite_display.tif" \
  --fps "${FPS}" \
  --border "${BORDER}" \
  --scale "${SCALE}" \
  --header-height "${HEADER_HEIGHT}" \
  --font-size "${FONT_SIZE}" \
  --font-stroke-width "${FONT_STROKE_WIDTH}" \
  --scale-bar-um "${SCALE_BAR_UM}" \
  --scale-bar-margin "${SCALE_BAR_MARGIN}" \
  --scale-bar-height "${SCALE_BAR_HEIGHT}" \
  --scale-bar-font-size "${SCALE_BAR_FONT_SIZE}" \
  --scale-bar-font-stroke-width "${SCALE_BAR_FONT_STROKE_WIDTH}" \
  --scale-bar-text-gap "${SCALE_BAR_TEXT_GAP}" \
  --quality "${QUALITY}" \
  --preset "${PRESET}" \
  --ffmpeg-threads "${FFMPEG_THREADS}" \
  --output "${OUTPUT}"
)

if [[ -n "${RGB_TIFF}" ]]; then
  CMD+=(--rgb-tiff "${RGB_TIFF}")
fi
if [[ -n "${XY_UM_PER_PIXEL}" ]]; then
  CMD+=(--xy-um-per-pixel "${XY_UM_PER_PIXEL}")
fi
if [[ -n "${SCALE_BAR_LABEL}" ]]; then
  CMD+=(--scale-bar-label "${SCALE_BAR_LABEL}")
fi

"${CMD[@]}"

echo "Fiji Z-scroll panel movie finished: $(date)"
ls -lh "${OUTPUT}"
if [[ -n "${RGB_TIFF}" ]]; then
  ls -lh "${RGB_TIFF}"
fi
