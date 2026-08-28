#!/usr/bin/env bash
# Great Lakes wrapper for locked Step 02c. All count and validation logic is in
# the local runner/Python files documented in PIPELINE_IO_AND_BIOLOGICAL_SCOPE.md.
# No normalization, cell-type calling, or biological classification occurs.
#SBATCH --job-name=egfp_probe_audit
#SBATCH --account=parent0
#SBATCH --partition=standard
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/logs/%x-%j.out

set -Eeuo pipefail

REPO_ROOT=/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
"${REPO_ROOT}/paper3_pcdh19/bin/run_egfp_probe_audit_all.sh"
