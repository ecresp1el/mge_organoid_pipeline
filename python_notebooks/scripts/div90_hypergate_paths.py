"""Keep code in the checkout and DIV90 runtime artifacts in PROJECT_ROOT."""
import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
PROJECT_ROOT = Path(os.environ.get(
    "PROJECT_ROOT", "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder"
)).expanduser().resolve()

if PROJECT_ROOT == REPO_ROOT or REPO_ROOT in PROJECT_ROOT.parents:
    raise ValueError("PROJECT_ROOT must be outside the Git checkout; use the Turbo runtime workspace.")
