"""Path helpers for repo code versus Great Lakes runtime data."""

import os
from pathlib import Path


DEFAULT_PROJECT_ROOT = Path("/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder")


def resolve_project_root(project_root=None):
    """Return the runtime PROJECT_ROOT path.

    Priority:
    1. explicit argument
    2. PROJECT_ROOT environment variable
    3. documented Great Lakes default
    """
    raw = project_root or os.environ.get("PROJECT_ROOT") or str(DEFAULT_PROJECT_ROOT)
    return Path(raw).expanduser().resolve()


def default_anndata_dir(project_root=None):
    """Return the canonical AnnData cache directory under PROJECT_ROOT."""
    root = resolve_project_root(project_root)
    return root / "results" / "python_anndata"


def ensure_under_path(path, parent):
    """Raise if path is not inside parent after resolution."""
    resolved = Path(path).expanduser().resolve()
    resolved_parent = Path(parent).expanduser().resolve()
    try:
        resolved.relative_to(resolved_parent)
    except ValueError:
        raise ValueError("Output path must remain under PROJECT_ROOT: {}".format(resolved))
    return resolved
