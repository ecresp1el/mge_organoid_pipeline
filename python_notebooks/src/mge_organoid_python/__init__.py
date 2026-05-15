"""Python entry points for notebook-based MGE organoid analysis."""

from .converter import SeuratToAnnDataConverter
from .paths import DEFAULT_PROJECT_ROOT, resolve_project_root
from .studies import StudySpec, default_studies
from .validation import AnnDataReport, validate_anndata, validate_source_paths

__all__ = [
    "AnnDataReport",
    "DEFAULT_PROJECT_ROOT",
    "SeuratToAnnDataConverter",
    "StudySpec",
    "default_studies",
    "resolve_project_root",
    "validate_anndata",
    "validate_source_paths",
]
