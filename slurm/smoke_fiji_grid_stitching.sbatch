#!/usr/bin/env bash
# Smoke-test the Fiji Grid/Collection Stitching macro on four tiny synthetic
# TIFF tiles before attempting the large Imaris IMS montage.

#SBATCH --job-name=fiji-stitch-smoke
#SBATCH --account=parent0
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --output=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/fiji-stitch-smoke-%j.out
#SBATCH --error=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/fiji-stitch-smoke-%j.err

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline}"
PROJECT_ROOT="${PROJECT_ROOT:-/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder}"
SMOKE_ROOT="${SMOKE_ROOT:-${PROJECT_ROOT}/results/imaging/fiji_stitching/smoke_test}"
INPUT_DIR="${SMOKE_ROOT}/tiles"
OUTPUT_DIR="${SMOKE_ROOT}/stitched"
CONDA_BASE="${CONDA_BASE:-/home/elcrespo/miniconda3}"
CONDA_ENV="${CONDA_ENV:-imaris-careamics}"
JAVA_BIN="${JAVA_BIN:-java}"
CLASS_DIR="${CLASS_DIR:-${TMPDIR:-/tmp}/${USER:-elcrespo}/fiji_stitcher_smoke_classes_${SLURM_JOB_ID:-manual}}"

rm -rf "${SMOKE_ROOT}"
mkdir -p "${INPUT_DIR}" "${OUTPUT_DIR}" "${PROJECT_ROOT}/logs"

echo "Fiji stitching smoke test started: $(date)"
echo "Host: $(hostname)"
echo "SMOKE_ROOT=${SMOKE_ROOT}"

source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

python - <<'PY'
from pathlib import Path
import numpy as np
import tifffile

root = Path("/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/fiji_stitching/smoke_test/tiles")
root.mkdir(parents=True, exist_ok=True)

canvas = np.zeros((112, 112), dtype=np.uint16)
yy, xx = np.mgrid[:112, :112]
canvas += ((xx * 23 + yy * 31) % 2048).astype(np.uint16)
canvas[20:38, 16:72] += 12000
canvas[58:96, 38:54] += 18000
canvas[70:88, 70:104] += 9000

tiles = {
    "tile_F3.tif": (0, 0),
    "tile_F0.tif": (48, 0),
    "tile_F2.tif": (0, 48),
    "tile_F1.tif": (48, 48),
}
for name, (x, y) in tiles.items():
    tifffile.imwrite(root / name, canvas[y:y + 64, x:x + 64], photometric="minisblack")

(root / "TileConfiguration.txt").write_text(
    "# Define the number of dimensions we are working on\n"
    "dim = 2\n\n"
    "# Define the image coordinates\n"
    "tile_F3.tif; ; (0.0, 0.0)\n"
    "tile_F0.tif; ; (48.0, 0.0)\n"
    "tile_F2.tif; ; (0.0, 48.0)\n"
    "tile_F1.tif; ; (48.0, 48.0)\n",
    encoding="utf-8",
)
PY

mkdir -p "${CLASS_DIR}"
FIJI_CP="$(find /sw/pkgs/med/fiji/1.5.4/jars /sw/pkgs/med/fiji/1.5.4/plugins -name '*.jar' | paste -sd: -)"
javac -cp "${FIJI_CP}" \
  -d "${CLASS_DIR}" \
  "${REPO_ROOT}/imaging/imaris_careamics/java/plugin/FijiGridCollectionStitcher.java"

"${JAVA_BIN}" \
  -Djava.awt.headless=true \
  -Xmx4g \
  -cp "${CLASS_DIR}:${FIJI_CP}" \
  plugin.FijiGridCollectionStitcher \
  "${INPUT_DIR}" \
  "TileConfiguration.txt" \
  "${OUTPUT_DIR}" \
  true \
  false

echo
echo "Smoke output inventory:"
find "${SMOKE_ROOT}" -maxdepth 3 -type f | sort
echo "Fiji stitching smoke test finished: $(date)"
