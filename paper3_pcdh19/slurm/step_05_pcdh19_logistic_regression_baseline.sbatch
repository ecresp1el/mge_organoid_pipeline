#!/usr/bin/env bash
# Great Lakes wrapper for Step 05. All model, evaluation, plot, validation, and
# publication logic remains in the supported runner and Python implementation.
#SBATCH --job-name=pcdh19_logistic
#SBATCH --account=parent0
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/logs/%x-%j.out

set -Eeuo pipefail

REPO_ROOT=/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
"${REPO_ROOT}/paper3_pcdh19/bin/run_step_05_pcdh19_logistic_regression_baseline.sh"
