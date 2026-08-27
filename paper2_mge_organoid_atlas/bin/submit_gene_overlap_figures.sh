#!/usr/bin/env bash

set -Eeuo pipefail

usage() {
  echo "Usage: $0 RUN_DIR [--svg true|false]" >&2
}

[[ $# -ge 1 ]] || { usage; exit 2; }
RUN_DIR="${1%/}"
shift
MAKE_SVG=false
while [[ $# -gt 0 ]]; do
  case "$1" in
    --svg)
      [[ $# -ge 2 && ("$2" == "true" || "$2" == "false") ]] || { usage; exit 2; }
      MAKE_SVG="$2"
      shift 2
      ;;
    *)
      usage
      exit 2
      ;;
  esac
done

BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
# shellcheck disable=SC1090
source "${BUNDLE_DIR}/config/greatlakes.env"
: "${PAPER2_ROOT:?PAPER2_ROOT is required}"
: "${CONDA_ENV_BIN:?CONDA_ENV_BIN is required}"

EXPECTED_PARENT="${PAPER2_ROOT}/results/02_harmonize_genes"
[[ "$(dirname "${RUN_DIR}")" == "${EXPECTED_PARENT}" ]] || {
  echo "Run must be directly under ${EXPECTED_PARENT}" >&2
  exit 2
}
[[ -s "${RUN_DIR}/SUCCESS.txt" ]] || { echo "Step 02 run is not complete: ${RUN_DIR}" >&2; exit 2; }
[[ -s "${RUN_DIR}/provenance/package_checksums.tsv" ]] || exit 2

echo "Validating the completed report package before adding figures..."
(
  cd "${RUN_DIR}"
  awk -F '\t' 'NR>1 {print $3 "  " $1}' provenance/package_checksums.tsv | sha256sum --check --status
)

PY_SCRIPT="${BUNDLE_DIR}/scripts/02_render_gene_overlap_figures.py"
SBATCH_TEMPLATE="${BUNDLE_DIR}/slurm/02_render_gene_overlap_figures.sbatch"
cp -p "${PY_SCRIPT}" "${SBATCH_TEMPLATE}" "${RUN_DIR}/code/"
cp -p "${SBATCH_TEMPLATE}" "${PAPER2_ROOT}/jobs/$(basename "${RUN_DIR}")_plot_overlap.sbatch"

{
  echo "FIGURE_DPI=300"
  echo "MAKE_SVG=${MAKE_SVG}"
  echo "SOURCE_TABLES=cached_step02_mapping_tables_only"
  echo "CANONICAL_OBJECTS_READ=false"
  echo "EXPRESSION_MATRICES_READ=false"
} > "${RUN_DIR}/config/gene_overlap_figure_render.env"
{
  echo "submitted_at=$(date --iso-8601=seconds)"
  echo "submit_host=$(hostname)"
  echo "submit_user=${USER:-unknown}"
  printf 'command='
  printf '%q ' "$0" "${RUN_DIR}" --svg "${MAKE_SVG}"
  printf '\n'
} > "${RUN_DIR}/provenance/gene_overlap_plot_submission.txt"
git -C "${REPO_ROOT}" rev-parse HEAD > "${RUN_DIR}/provenance/gene_overlap_plot_git_commit.txt"
git -C "${REPO_ROOT}" status --short > "${RUN_DIR}/provenance/gene_overlap_plot_git_status.txt"
git -C "${REPO_ROOT}" diff --binary > "${RUN_DIR}/provenance/gene_overlap_plot_working_tree.patch"

JOB_ID="$(sbatch --parsable \
  --job-name="p2-gene-plots" \
  --account="${ACCOUNT}" \
  --partition="${PARTITION}" \
  --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=16G --time=02:00:00 \
  --output="${RUN_DIR}/logs/gene-overlap-figures-%j.out" \
  --error="${RUN_DIR}/logs/gene-overlap-figures-%j.err" \
  --export="ALL,PAPER2_GENE_RUN_DIR=${RUN_DIR},CONDA_ENV_BIN=${CONDA_ENV_BIN},MAKE_SVG=${MAKE_SVG},FIGURE_DPI=300" \
  "${PAPER2_ROOT}/jobs/$(basename "${RUN_DIR}")_plot_overlap.sbatch")"
printf '%s\n' "${JOB_ID}" > "${RUN_DIR}/provenance/gene_overlap_plot_job_id.txt"
echo "Submitted plot-only job ${JOB_ID}"
echo "Run directory: ${RUN_DIR}"
