"""External runtime root for the 2026-10-31 downstream exploratory work."""
import os
from pathlib import Path
RUNTIME_ROOT = Path(os.environ.get("PCDH19_RUNTIME_ROOT", "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19")).expanduser().resolve()
REPO_STUDY_ROOT = Path(__file__).resolve().parents[1]
if RUNTIME_ROOT == REPO_STUDY_ROOT or REPO_STUDY_ROOT in RUNTIME_ROOT.parents:
    raise ValueError("PCDH19_RUNTIME_ROOT must be outside the code checkout")
