#!/usr/bin/env bash
# Freeze and submit a serial, verbose, scores-only Step 03 diagnostic.
set -Eeuo pipefail

BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
# shellcheck disable=SC1090
source "${BUNDLE_DIR}/config/greatlakes.env"
# shellcheck disable=SC1090
source "${BUNDLE_DIR}/config/primary_processing_step03.env"

: "${REPO_ROOT:?}" "${PAPER3_ROOT:?}" "${ACCOUNT:?}"

PRODUCTION_RUN_ID="03_scdblfinder_20260831_142124_d8a7242"
PRODUCTION_RUN_DIR="${PAPER3_ROOT}/results/primary_processing/03_scdblfinder/${PRODUCTION_RUN_ID}"
STEP02_RUN_DIR="${PAPER3_ROOT}/results/primary_processing/02_qc_filtering/${PRIMARY_PROCESSING_STEP03_STEP02_RUN_ID}"
INPUT_H5AD="${STEP02_RUN_DIR}/${PRIMARY_PROCESSING_STEP03_STEP02_CHECKPOINT}"
BRIDGE_H5="${PRODUCTION_RUN_DIR}/intermediate/step02_counts_for_scdblfinder.h5"
METADATA_TSV="${PRODUCTION_RUN_DIR}/intermediate/step03_cell_metadata.tsv"
PREPARE_CHECKS="${PRODUCTION_RUN_DIR}/intermediate/prepare_validation_checks.tsv"
R_SOURCE="${BUNDLE_DIR}/scripts/primary_processing/step03_scdblfinder_diagnostic.R"
SBATCH_SOURCE="${BUNDLE_DIR}/slurm/primary_processing_03_scdblfinder_diagnostic.sbatch"

for required in "${INPUT_H5AD}" "${BRIDGE_H5}" "${METADATA_TSV}" "${PREPARE_CHECKS}" "${R_SOURCE}" "${SBATCH_SOURCE}"; do
  [[ -f "${required}" ]] || { echo "Missing diagnostic asset: ${required}" >&2; exit 2; }
done

[[ -z "$(git -C "${REPO_ROOT}" status --porcelain)" ]] || { echo "Refusing diagnostic submission from a dirty repository" >&2; exit 2; }

module load "${PRIMARY_PROCESSING_STEP03_R_MODULE}"
export R_LIBS_USER="${PRIMARY_PROCESSING_STEP03_R_LIBRARY}"
Rscript - "${BUNDLE_DIR}/config/primary_processing_step03_r_packages.tsv" "${R_SOURCE}" <<'RS'
args <- commandArgs(trailingOnly=TRUE)
expected <- read.delim(args[[1]], stringsAsFactors=FALSE, check.names=FALSE)
observed <- vapply(expected$package, function(package) {
  if (package == "R") paste(R.version$major, R.version$minor, sep=".") else as.character(packageVersion(package))
}, character(1L))
if (!identical(unname(observed), expected$version)) stop("Diagnostic R environment mismatch")
invisible(parse(file=args[[2]]))
RS

DIAGNOSTIC_ROOT="${PAPER3_ROOT}/results/primary_processing/03_scdblfinder_diagnostics"
STAMP="$(date +%Y%m%d_%H%M%S)"
GIT_SHORT="$(git -C "${REPO_ROOT}" rev-parse --short HEAD)"
RUN_ID="03_scdblfinder_diagnostic_${STAMP}_${GIT_SHORT}"
RUN_DIR="${DIAGNOSTIC_ROOT}/${RUN_ID}"
[[ ! -e "${RUN_DIR}" ]] || { echo "Refusing existing diagnostic run: ${RUN_DIR}" >&2; exit 2; }
mkdir -p "${RUN_DIR}/code" "${RUN_DIR}/config" "${RUN_DIR}/logs"

cp -p "${R_SOURCE}" "${RUN_DIR}/code/"
cp -p "${SBATCH_SOURCE}" "${RUN_DIR}/code/"
cp -p "${BUNDLE_DIR}/config/primary_processing_step03.env" "${RUN_DIR}/config/"
cp -p "${BUNDLE_DIR}/config/primary_processing_step03_r_packages.tsv" "${RUN_DIR}/config/"
printf 'production_run_id=%s\nproduction_job_id=59389340\ninput_h5ad=%s\nbridge_h5=%s\nmetadata_tsv=%s\nprepare_checks=%s\nscientific_role=instrumentation_only\ncells_removed=0\npublished_as_step03=FALSE\n' \
  "${PRODUCTION_RUN_ID}" "${INPUT_H5AD}" "${BRIDGE_H5}" "${METADATA_TSV}" "${PREPARE_CHECKS}" > "${RUN_DIR}/config/resolved_inputs.txt"
printf 'git_commit=%s\ngit_status=clean\nsubmitted_at=%s\n' \
  "$(git -C "${REPO_ROOT}" rev-parse HEAD)" "$(date --iso-8601=seconds)" > "${RUN_DIR}/config/submission_provenance.txt"

SCHEDULER_STDOUT="${RUN_DIR}/logs/scheduler.stdout.log"
SCHEDULER_STDERR="${RUN_DIR}/logs/scheduler.stderr.log"
EXPORTS="ALL,PAPER3_DIAGNOSTIC_RUN_DIR=${RUN_DIR},PAPER3_DIAGNOSTIC_R_SCRIPT=${RUN_DIR}/code/$(basename "${R_SOURCE}"),PAPER3_DIAGNOSTIC_BRIDGE_H5=${BRIDGE_H5},PAPER3_DIAGNOSTIC_METADATA_TSV=${METADATA_TSV},PAPER3_DIAGNOSTIC_PREPARE_CHECKS=${PREPARE_CHECKS},PAPER3_DIAGNOSTIC_INPUT_H5AD=${INPUT_H5AD},PAPER3_DIAGNOSTIC_EXPECTED_INPUT_SHA256=${PRIMARY_PROCESSING_STEP03_INPUT_SHA256},PAPER3_DIAGNOSTIC_EXPECTED_CELLS=${PRIMARY_PROCESSING_STEP03_EXPECTED_CELLS},PAPER3_DIAGNOSTIC_EXPECTED_GENES=${PRIMARY_PROCESSING_STEP03_EXPECTED_GENES},PAPER3_DIAGNOSTIC_CAPTURE_ID=${PRIMARY_PROCESSING_STEP03_CAPTURE_ID},PAPER3_DIAGNOSTIC_PRIMARY_SEED=${PRIMARY_PROCESSING_STEP03_PRIMARY_SEED},PAPER3_DIAGNOSTIC_R_MODULE=${PRIMARY_PROCESSING_STEP03_R_MODULE},PAPER3_DIAGNOSTIC_R_LIBRARY=${PRIMARY_PROCESSING_STEP03_R_LIBRARY}"

JOB_ID="$(sbatch --parsable \
  --job-name=pcdh19-scdbl-diag \
  --account="${ACCOUNT}" \
  --partition="${PRIMARY_PROCESSING_STEP03_PARTITION}" \
  --nodes=1 --ntasks=1 --cpus-per-task=1 \
  --mem="${PRIMARY_PROCESSING_STEP03_MEMORY}" \
  --time="${PRIMARY_PROCESSING_STEP03_WALLTIME}" \
  --output="${SCHEDULER_STDOUT}" \
  --error="${SCHEDULER_STDERR}" \
  --export="${EXPORTS}" \
  "${RUN_DIR}/code/$(basename "${SBATCH_SOURCE}")")"

printf 'stage\tjob_id\nserial_observability\t%s\n' "${JOB_ID}" > "${RUN_DIR}/config/job_id.tsv"
printf 'DIAGNOSTIC_JOB_ID=%s\nDIAGNOSTIC_RUN_DIR=%s\nSCHEDULER_STDOUT=%s\nSCHEDULER_STDERR=%s\nR_STDOUT=%s\nR_STDERR=%s\nHEARTBEAT=%s\n' \
  "${JOB_ID}" "${RUN_DIR}" "${SCHEDULER_STDOUT}" "${SCHEDULER_STDERR}" \
  "${RUN_DIR}/logs/r.stdout.log" "${RUN_DIR}/logs/r.stderr.log" "${RUN_DIR}/logs/resource_heartbeat.tsv"
