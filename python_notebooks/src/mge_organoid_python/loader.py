"""Load cached AnnData files produced by the Seurat conversion jobs."""

from pathlib import Path

import anndata as ad

from .paths import resolve_project_root
from .validation import validate_anndata


def cached_h5ad_path(study, project_root=None, output_dir=None):
    """Return the expected cached H5AD path for a study."""
    return study.h5ad_path(project_root=resolve_project_root(project_root), output_dir=output_dir)


def missing_cached_h5ads(studies, project_root=None, output_dir=None):
    """Return `(study_id, path)` entries for missing cached H5AD files."""
    missing = []
    for study in studies:
        path = cached_h5ad_path(study, project_root=project_root, output_dir=output_dir)
        if not path.exists():
            missing.append((study.study_id, str(path)))
    return missing


def load_cached_anndata(study, project_root=None, output_dir=None, backed="r"):
    """Load one cached H5AD file and return `(adata, report)`.

    By default this opens the file in backed read-only mode, which avoids loading
    the full expression matrix into memory. Use `backed=None` only when you
    intentionally want an in-memory AnnData object.
    """
    path = cached_h5ad_path(study, project_root=project_root, output_dir=output_dir)
    if not path.exists():
        raise FileNotFoundError(
            "Missing cached H5AD for {}: {}. Rebuild it with the Slurm batch "
            "template, then rerun this notebook.".format(study.study_id, path)
        )
    adata = ad.read_h5ad(path, backed=backed)
    report = validate_anndata(study, adata, Path(path))
    return adata, report


def load_cached_anndatas(studies, project_root=None, output_dir=None, backed="r"):
    """Load all cached H5AD files into a dict keyed by study_id."""
    adatas = {}
    reports = []
    for study in studies:
        adata, report = load_cached_anndata(
            study,
            project_root=project_root,
            output_dir=output_dir,
            backed=backed,
        )
        adatas[study.study_id] = adata
        reports.append(report)
    return adatas, reports
