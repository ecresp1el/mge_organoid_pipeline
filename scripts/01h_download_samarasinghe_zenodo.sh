#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder}"
URL="${SAMARASINGHE_ZENODO_URL:-https://zenodo.org/records/5732813/files/Samarasinghe_2021_seurat_object.RData?download=1}"
EXPECTED_MD5="${SAMARASINGHE_ZENODO_MD5:-e124d3831855817e8649f05dee13f0ce}"
OUT_DIR="${PROJECT_ROOT}/data/raw/samarasinghe_2021_zenodo"
OUT_FILE="${OUT_DIR}/Samarasinghe_2021_seurat_object.RData"

mkdir -p "${OUT_DIR}"

echo "Downloading Samarasinghe 2021 Zenodo processed Seurat object"
echo "URL: ${URL}"
echo "Destination: ${OUT_FILE}"

curl -L --fail --continue-at - --progress-bar -o "${OUT_FILE}" "${URL}"

actual_md5="$(md5sum "${OUT_FILE}" | awk '{print $1}')"
if [[ "${actual_md5}" != "${EXPECTED_MD5}" ]]; then
  echo "ERROR: MD5 mismatch for ${OUT_FILE}" >&2
  echo "Expected: ${EXPECTED_MD5}" >&2
  echo "Observed: ${actual_md5}" >&2
  exit 1
fi

echo "MD5 OK: ${actual_md5}"
