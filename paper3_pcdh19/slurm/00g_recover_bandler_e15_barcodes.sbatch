#!/usr/bin/env bash
# Scheduler resources and logs are assigned by the standalone submission wrapper.
set -Eeuo pipefail

: "${PAPER3_CURATION_RUN_DIR:?completed curation run directory is required}"
: "${PAPER3_CURATION_SOURCE_ROOT:?Bandler source cache is required}"
: "${PAPER3_CURATION_PYTHON_BIN:?Python interpreter is required}"
: "${PAPER3_CURATION_R_MODULE:?R module is required}"

RUN_DIR="${PAPER3_CURATION_RUN_DIR}"
OUTPUT_DIR="${RUN_DIR}/Bandler2022/interactive_atlas/barcode_recovery"
PY_SCRIPT="${RUN_DIR}/code/recover_bandler_e15_barcodes.py"
R_SCRIPT="${RUN_DIR}/code/export_bandler_e15_fingerprints.R"
REPORT_SCRIPT="${RUN_DIR}/code/generate_reference_curation_report.py"

failed() {
  local line="$1" status="$2"
  printf 'FAILED\nstage=bandler_e15_barcode_recovery\njob_id=%s\nline=%s\nexit_status=%s\ntime=%s\n' \
    "${SLURM_JOB_ID:-unknown}" "${line}" "${status}" "$(date --iso-8601=seconds)" \
    > "${RUN_DIR}/BANDLER_E15_BARCODE_RECOVERY_FAILED.txt"
}
trap 'failed "${LINENO}" "$?"' ERR

module load Bioinformatics
module load "${PAPER3_CURATION_R_MODULE}"
Rscript "${R_SCRIPT}" \
  --source-root "${PAPER3_CURATION_SOURCE_ROOT}" \
  --genes "Nfib,Hes1,Fabp7,Top2a,Ccnd2,Nkx2-1,Lhx8,Npy,Nxph1,Nr2f2,Isl1,Foxp1,Gucy1a3,Ebf1,Sst,Maf,Lhx6,Rnd2,Zfp536,Pcsk1n,Nnat,Ackr3,Mef2c,Erbb4" \
  --output "${OUTPUT_DIR}/audit/deposited_e15_expression_fingerprints.tsv" \
  --inventory "${OUTPUT_DIR}/audit/deposited_e15_matrix_inventory.tsv"

"${PAPER3_CURATION_PYTHON_BIN}" "${PY_SCRIPT}" publish \
  --run-dir "${RUN_DIR}" \
  --source-root "${PAPER3_CURATION_SOURCE_ROOT}" \
  --fingerprints "${OUTPUT_DIR}/audit/deposited_e15_expression_fingerprints.tsv"
"${PAPER3_CURATION_PYTHON_BIN}" "${REPORT_SCRIPT}" --run-dir "${RUN_DIR}"

rm -f "${RUN_DIR}/BANDLER_E15_BARCODE_RECOVERY_FAILED.txt"
printf 'PASS\ncompleted=%s\njob_id=%s\n' "$(date --iso-8601=seconds)" "${SLURM_JOB_ID:-unknown}" \
  > "${RUN_DIR}/BANDLER_E15_BARCODE_RECOVERY_SUCCESS.txt"
