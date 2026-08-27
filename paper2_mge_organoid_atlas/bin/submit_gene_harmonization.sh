#!/usr/bin/env bash

set -Eeuo pipefail

usage() {
  echo "Usage: $0 [--dry-run]" >&2
}

DRY_RUN=false
ORIGINAL_ARGS=("$@")
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run)
      DRY_RUN=true
      shift
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 2
      ;;
  esac
done

BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
GREATLAKES_CONFIG="${BUNDLE_DIR}/config/greatlakes.env"
GENE_CONFIG="${BUNDLE_DIR}/config/gene_harmonization.env"
PY_SCRIPT="${BUNDLE_DIR}/scripts/02_audit_harmonize_genes.py"
SBATCH_TEMPLATE="${BUNDLE_DIR}/slurm/02_harmonize_genes.sbatch"

# shellcheck disable=SC1090
source "${GREATLAKES_CONFIG}"
# shellcheck disable=SC1090
source "${GENE_CONFIG}"
: "${REPO_ROOT:?REPO_ROOT is required}"
: "${PAPER2_ROOT:?PAPER2_ROOT is required}"
: "${ACCOUNT:?ACCOUNT is required}"
: "${PARTITION:?PARTITION is required}"
: "${CONDA_ENV_BIN:?CONDA_ENV_BIN is required}"

CANONICAL_DIR="${PAPER2_ROOT}/inputs/canonical"
[[ -s "${CANONICAL_DIR}/FROZEN.txt" ]] || {
  echo "Frozen canonical inputs are missing: ${CANONICAL_DIR}" >&2
  exit 2
}
[[ -s "${CANONICAL_DIR}/provenance/canonical_dataset_manifest.tsv" ]] || exit 2
[[ -s "${CANONICAL_DIR}/provenance/canonical_file_checksums.tsv" ]] || exit 2

for required in "${GENE_CONFIG}" "${PY_SCRIPT}" "${SBATCH_TEMPLATE}"; do
  [[ -s "${required}" ]] || { echo "Missing required file: ${required}" >&2; exit 2; }
done

echo "Step: 02_harmonize_genes (report-only review stop)"
echo "Canonical input: ${CANONICAL_DIR}"
echo "Mapping authority: GENCODE ${GENCODE_RELEASE} (${GENCODE_ASSEMBLY}) plus frozen HGNC complete-set snapshot"
echo "No expression matrices will be loaded or changed."
if [[ "${DRY_RUN}" == "true" ]]; then
  echo "Dry run passed; nothing was downloaded, created, or submitted."
  exit 0
fi

"${BUNDLE_DIR}/bin/initialize_turbo.sh" >/dev/null
STAMP="$(date +%Y%m%d_%H%M%S)"
GIT_SHORT="$(git -C "${REPO_ROOT}" rev-parse --short HEAD)"
RUN_ID="02_harmonize_genes_${STAMP}_${GIT_SHORT}"
STEP_DIR="${PAPER2_ROOT}/results/02_harmonize_genes"
RUN_DIR="${STEP_DIR}/${RUN_ID}"
JOB_FILE="${PAPER2_ROOT}/jobs/${RUN_ID}.sbatch"
[[ ! -e "${RUN_DIR}" ]] || { echo "Run already exists: ${RUN_DIR}" >&2; exit 2; }

mkdir -p \
  "${RUN_DIR}/code" \
  "${RUN_DIR}/config" \
  "${RUN_DIR}/reference" \
  "${RUN_DIR}/tables" \
  "${RUN_DIR}/logs" \
  "${RUN_DIR}/provenance"

cp -p "${PY_SCRIPT}" "${SBATCH_TEMPLATE}" "${RUN_DIR}/code/"
cp -p "${GREATLAKES_CONFIG}" "${RUN_DIR}/config/submitted_greatlakes.env"
cp -p "${GENE_CONFIG}" "${RUN_DIR}/config/gene_harmonization.env"
cp -p "${CANONICAL_DIR}/provenance/canonical_dataset_manifest.tsv" \
  "${RUN_DIR}/config/canonical_dataset_manifest.tsv"
cp -p "${CANONICAL_DIR}/provenance/canonical_file_checksums.tsv" \
  "${RUN_DIR}/config/canonical_file_checksums.tsv"
cp -p "${CANONICAL_DIR}/FROZEN.txt" "${RUN_DIR}/config/canonical_FROZEN.txt"
cp -p "${SBATCH_TEMPLATE}" "${JOB_FILE}"

echo "Downloading frozen GENCODE ${GENCODE_RELEASE} reference..."
curl --fail --location --retry 3 --output \
  "${RUN_DIR}/reference/gencode.v${GENCODE_RELEASE}.annotation.gtf.gz" \
  "${GENCODE_GTF_URL}"
printf '%s  %s\n' "${GENCODE_GTF_MD5}" \
  "${RUN_DIR}/reference/gencode.v${GENCODE_RELEASE}.annotation.gtf.gz" | md5sum --check --status

echo "Downloading frozen HGNC complete-set snapshot..."
curl --fail --location --retry 3 --output \
  "${RUN_DIR}/reference/hgnc_complete_set.txt" \
  "${HGNC_COMPLETE_SET_URL}"

{
  echo -e "reference\turl\tretrieved_at\tsize_bytes\tsha256\tupstream_md5"
  for reference_file in \
    "${RUN_DIR}/reference/gencode.v${GENCODE_RELEASE}.annotation.gtf.gz" \
    "${RUN_DIR}/reference/hgnc_complete_set.txt"; do
    if [[ "$(basename "${reference_file}")" == gencode.* ]]; then
      reference_url="${GENCODE_GTF_URL}"
      upstream_md5="${GENCODE_GTF_MD5}"
    else
      reference_url="${HGNC_COMPLETE_SET_URL}"
      upstream_md5=""
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$(basename "${reference_file}")" "${reference_url}" "$(date --iso-8601=seconds)" \
      "$(stat -c %s "${reference_file}")" "$(sha256sum "${reference_file}" | cut -d ' ' -f 1)" \
      "${upstream_md5}"
  done
} > "${RUN_DIR}/provenance/reference_manifest.tsv"

{
  echo "REPO_ROOT=${REPO_ROOT}"
  echo "PAPER2_ROOT=${PAPER2_ROOT}"
  echo "CANONICAL_DIR=${CANONICAL_DIR}"
  echo "RUN_ID=${RUN_ID}"
  echo "RUN_DIR=${RUN_DIR}"
  echo "GENCODE_RELEASE=${GENCODE_RELEASE}"
  echo "GENCODE_ASSEMBLY=${GENCODE_ASSEMBLY}"
  echo "COMMON_GENE_NAMESPACE=${COMMON_GENE_NAMESPACE}"
  echo "CONDA_ENV_BIN=${CONDA_ENV_BIN}"
} > "${RUN_DIR}/config/resolved.env"

{
  echo "submitted_at=$(date --iso-8601=seconds)"
  echo "submit_host=$(hostname)"
  echo "submit_user=${USER:-unknown}"
  printf 'command='
  printf '%q ' "$0" "${ORIGINAL_ARGS[@]}"
  printf '\n'
} > "${RUN_DIR}/provenance/submission.txt"
git -C "${REPO_ROOT}" rev-parse HEAD > "${RUN_DIR}/provenance/git_commit.txt"
git -C "${REPO_ROOT}" status --short > "${RUN_DIR}/provenance/git_status.txt"
git -C "${REPO_ROOT}" diff --binary > "${RUN_DIR}/provenance/git_working_tree.patch"

JOB_ID="$(sbatch --parsable \
  --job-name="p2-gene-audit" \
  --account="${ACCOUNT}" \
  --partition="${PARTITION}" \
  --nodes=1 \
  --ntasks=1 \
  --cpus-per-task="${GENE_AUDIT_CPUS}" \
  --mem="${GENE_AUDIT_MEMORY}" \
  --time="${GENE_AUDIT_WALLTIME}" \
  --output="${RUN_DIR}/logs/slurm-%j.out" \
  --error="${RUN_DIR}/logs/slurm-%j.err" \
  --export="ALL,PAPER2_GENE_RUN_DIR=${RUN_DIR},CONDA_ENV_BIN=${CONDA_ENV_BIN}" \
  "${JOB_FILE}")"
printf '%s\n' "${JOB_ID}" > "${RUN_DIR}/provenance/slurm_job_id.txt"

echo "Submitted job ${JOB_ID}"
echo "Run directory: ${RUN_DIR}"
echo "Monitor: ${BUNDLE_DIR}/bin/status.sh ${RUN_DIR}"
