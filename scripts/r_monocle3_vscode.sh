#!/usr/bin/env bash
set -euo pipefail

# R launcher for VS Code's R extension on Great Lakes compute nodes.
# VS Code does not inherit modules loaded in an arbitrary terminal, so this
# wrapper loads the Monocle3 module stack before starting R.

if ! type module >/dev/null 2>&1; then
  if [[ -f /etc/profile.d/modules.sh ]]; then
    # shellcheck disable=SC1091
    source /etc/profile.d/modules.sh
  elif [[ -f /usr/share/lmod/lmod/init/bash ]]; then
    # shellcheck disable=SC1091
    source /usr/share/lmod/lmod/init/bash
  fi
fi

module load Bioinformatics
module load Rmonocle3/1.3.7

export PROJECT_ROOT="${PROJECT_ROOT:-/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder}"

exec /sw/pkgs/arc/stacks/gcc/13.2.0/R/4.4.0/bin/R "$@"
