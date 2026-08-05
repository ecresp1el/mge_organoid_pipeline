#!/usr/bin/env bash
# Resource requests and log paths are supplied by submit.sh so each run is isolated.

set -Eeuo pipefail

: "${PIPELINE_RUN_DIR:?submit this job through submit.sh}"
RUN_DIR="${PIPELINE_RUN_DIR}"
CONFIG="${RUN_DIR}/config/resolved.env"

mkdir -p "${RUN_DIR}/logs" "${RUN_DIR}/work/classes" "${RUN_DIR}/results/planes" \
  "${RUN_DIR}/results/stacks"
exec > >(tee -a "${RUN_DIR}/logs/pipeline.log") 2>&1

failed() {
  local line="$1"
  local status="$2"
  printf 'FAILED\njob_id=%s\nline=%s\nexit_status=%s\ntime=%s\n' \
    "${SLURM_JOB_ID:-unknown}" "${line}" "${status}" "$(date --iso-8601=seconds)" \
    > "${RUN_DIR}/FAILED.txt"
}
trap 'failed "${LINENO}" "$?"' ERR

source "${CONFIG}"
IFS=' ' read -r -a CHANNEL_LIST <<< "${CHANNELS}"
CODE_DIR="${RUN_DIR}/code"
INPUTS_DIR="${RUN_DIR}/inputs"
PLANES_DIR="${RUN_DIR}/results/planes"
STACKS_DIR="${RUN_DIR}/results/stacks"

echo "Portable Fiji pipeline started: $(date --iso-8601=seconds)"
echo "RUN_DIR=${RUN_DIR}"
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-}"
echo "MODE=${MODE}"
echo "Host=$(hostname)"

module purge
module load fiji/1.5.4
module load python3.11-anaconda/2024.02

{
  echo "timestamp=$(date --iso-8601=seconds)"
  echo "host=$(hostname)"
  echo "slurm_job_id=${SLURM_JOB_ID:-}"
  echo "loaded_modules=${LOADEDMODULES:-}"
  echo "java=$(java -version 2>&1 | head -1)"
  echo "javac=$(javac -version 2>&1)"
  echo "python=$(python --version 2>&1)"
  python - <<'PY'
import numpy, tifffile, PIL
print(f"numpy={numpy.__version__}")
print(f"tifffile={tifffile.__version__}")
print(f"Pillow={PIL.__version__}")
PY
} > "${RUN_DIR}/config/environment.txt"

if [[ "${MODE}" == "smoke" ]]; then
  python "${CODE_DIR}/src/pipeline_io.py" make-smoke --destination "${INPUTS_DIR}"
fi

if [[ ! -f "${INPUTS_DIR}/TileConfiguration.txt" ]]; then
  echo "Missing staged TileConfiguration.txt in ${INPUTS_DIR}" >&2
  exit 2
fi

echo "Tile configuration used:"
sed -n '1,120p' "${INPUTS_DIR}/TileConfiguration.txt"

FIJI_ROOT="$(dirname "$(command -v ImageJ-linux64)")"
FIJI_CP="$(find "${FIJI_ROOT}/jars" "${FIJI_ROOT}/plugins" -name '*.jar' -print | paste -sd: -)"
javac -cp "${FIJI_CP}" \
  -d "${RUN_DIR}/work/classes" \
  "${CODE_DIR}/src/FijiGridCollectionStitcher.java"

java -Djava.awt.headless=true -Xmx"${JAVA_HEAP}" \
  -cp "${RUN_DIR}/work/classes:${FIJI_CP}" \
  plugin.FijiGridCollectionStitcher \
  "${INPUTS_DIR}" TileConfiguration.txt "${PLANES_DIR}" false false

python "${CODE_DIR}/src/pipeline_io.py" stack \
  --input-dir "${PLANES_DIR}" \
  --output-dir "${STACKS_DIR}" \
  --prefix "${OUTPUT_PREFIX}" \
  --channels "${CHANNEL_LIST[@]}"

python "${CODE_DIR}/src/pipeline_io.py" validate \
  --mode "${MODE}" --inputs-dir "${INPUTS_DIR}" \
  --planes-dir "${PLANES_DIR}" --stacks-dir "${STACKS_DIR}" \
  --prefix "${OUTPUT_PREFIX}" --channels "${CHANNEL_LIST[@]}" \
  --output "${RUN_DIR}/validation.json"

printf 'PASS\njob_id=%s\ncompleted=%s\nvalidation=%s\n' \
  "${SLURM_JOB_ID:-unknown}" "$(date --iso-8601=seconds)" "${RUN_DIR}/validation.json" \
  > "${RUN_DIR}/SUCCESS.txt"

MANIFEST_ARGS=(
  python "${CODE_DIR}/src/pipeline_io.py" manifest
  --run-dir "${RUN_DIR}" --output "${RUN_DIR}/manifest.json"
)
if [[ "${CHECKSUM_OUTPUTS}" == "true" ]]; then
  MANIFEST_ARGS+=(--checksum-outputs)
fi
"${MANIFEST_ARGS[@]}"

rm -f "${RUN_DIR}/FAILED.txt"
echo "Portable Fiji pipeline finished successfully: $(date --iso-8601=seconds)"
echo "Validation: ${RUN_DIR}/validation.json"
