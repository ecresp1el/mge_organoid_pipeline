#!/usr/bin/env bash
set -Eeuo pipefail

: "${PAPER3_MAPPING_RUN_DIR:?}"
: "${PAPER3_CELLRANGER_ROOT:?}"
: "${PAPER3_MAPPING_R_MODULE:?}"
: "${PAPER3_MAPPING_COMPILER_MODULE:?}"
: "${PAPER3_MAPPING_HDF5_MODULE:?}"
: "${PAPER3_MAPPING_R_ADDON_LIB:?}"
: "${PAPER3_MIND_THRESHOLD:?}"
: "${PAPER3_MAPPING_SEED:?}"

RUN_DIR="${PAPER3_MAPPING_RUN_DIR}"
trap 'printf "FAILED\nstage=bandler_mind\njob_id=%s\nline=%s\n" "${SLURM_JOB_ID:-unknown}" "${LINENO}" > "${RUN_DIR}/BANDLER_MIND_FAILED.txt"' ERR

module load "${PAPER3_MAPPING_COMPILER_MODULE}"
module load "${PAPER3_MAPPING_HDF5_MODULE}"
module load Bioinformatics
module load "${PAPER3_MAPPING_R_MODULE}"
export R_LIBS_USER="${PAPER3_MAPPING_R_ADDON_LIB}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
Rscript "${RUN_DIR}/code/transfer_bandler_mind_labels.R" \
  --cellranger-root "${PAPER3_CELLRANGER_ROOT}" \
  --sample-key "${RUN_DIR}/config/sample_key.csv" \
  --query-metadata "${RUN_DIR}/query/query_cell_metadata.tsv.gz" \
  --reference-matrix "${PAPER3_BANDLER_CA301_MATRIX}" \
  --reference-join "${RUN_DIR}/inputs/CA301_later_atlas_barcode_join.tsv" \
  --output-dir "${RUN_DIR}/Bandler_MIND" \
  --threshold "${PAPER3_MIND_THRESHOLD}" \
  --seed "${PAPER3_MAPPING_SEED}" \
  --cpus "${SLURM_CPUS_PER_TASK:-1}"

rm -f "${RUN_DIR}/BANDLER_MIND_FAILED.txt"
printf 'PASS\njob_id=%s\ncompleted=%s\n' "${SLURM_JOB_ID:-unknown}" "$(date --iso-8601=seconds)" > "${RUN_DIR}/BANDLER_MIND_SUCCESS.txt"
