"""Python entry points for notebook-based MGE organoid analysis."""

from .converter import SeuratToAnnDataConverter
from .loader import cached_h5ad_path, load_cached_anndata, load_cached_anndatas, missing_cached_h5ads
from .paths import DEFAULT_PROJECT_ROOT, resolve_project_root
from .studies import StudySpec, default_studies
from .validation import AnnDataReport, validate_anndata, validate_source_paths

__all__ = [
    "AnnDataReport",
    "DEFAULT_PROJECT_ROOT",
    "SeuratToAnnDataConverter",
    "StudySpec",
    "cached_h5ad_path",
    "default_studies",
    "load_cached_anndata",
    "load_cached_anndatas",
    "missing_cached_h5ads",
    "resolve_project_root",
    "validate_anndata",
    "validate_source_paths",
]
