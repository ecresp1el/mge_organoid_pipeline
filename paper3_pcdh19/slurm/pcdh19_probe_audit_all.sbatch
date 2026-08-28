#!/usr/bin/env bash
# PURPOSE
#   Great Lakes batch wrapper for Step 02a. It allocates resources and invokes
#   bin/run_pcdh19_probe_audit_all.sh; it contains no scientific computation.
#
# SCHEDULER CONTRACT
#   Account parent0, standard partition, one task using two CPUs, 8 GiB RAM,
#   and a four-hour wall-time limit. SLURM expands %x/%j in the scheduler log.
#
# INPUTS AND OUTPUTS
#   The fixed repository checkout supplies code/configuration. The runner reads
#   the checksum-locked references and source Cell Ranger assets configured in
#   greatlakes.env, writes the Step 02a result package, and creates its own
#   timestamped run log. SLURM stdout/stderr additionally goes to the declared
#   PAPER3_ROOT/logs/%x-%j.out path, whose parent must exist before submission.
#
# FAILURE/RESTART AND SCIENTIFIC SCOPE
#   Strict shell mode propagates failures to SLURM. Reruns reuse only validated
#   identical outputs. Scheduling does not change the runner's raw-count scope:
#   no normalization, biological labeling, cell typing, or inference occurs.
#   See ../PIPELINE_IO_AND_BIOLOGICAL_SCOPE.md for the full I/O contract.
#SBATCH --job-name=pcdh19_probe_audit
#SBATCH --account=parent0
#SBATCH --partition=standard
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/logs/%x-%j.out

set -Eeuo pipefail

REPO_ROOT=/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
"${REPO_ROOT}/paper3_pcdh19/bin/run_pcdh19_probe_audit_all.sh"
