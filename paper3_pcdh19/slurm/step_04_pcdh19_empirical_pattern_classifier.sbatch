#!/usr/bin/env bash
# Great Lakes wrapper for Step 04. It allocates resources and invokes the
# supported runner; no scientific or plotting logic lives here.
#SBATCH --job-name=pcdh19_pattern_clf
#SBATCH --account=parent0
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/logs/%x-%j.out

set -Eeuo pipefail

REPO_ROOT=/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
"${REPO_ROOT}/paper3_pcdh19/bin/run_step_04_pcdh19_empirical_pattern_classifier.sh"
