"""Python entry points for notebook-based MGE organoid analysis."""

from .converter import SeuratToAnnDataConverter
from .data_sources import (
    DatasetLoadResult,
    Notebook00SourceConfig,
    cellbender_output_table,
    expected_cellbender_output_h5,
    load_dataset_result,
    sample_table,
    source_availability_table,
    summarize_source_availability,
)
from .loader import cached_h5ad_path, load_cached_anndata, load_cached_anndatas, missing_cached_h5ads
from .notebook00_plots import PlotConfig
from .notebook00_workflow import (
    ManualECCheckpointPaths,
    ManualECCheckpointSettings,
    ManualECFilterSettings,
    save_manual_ec_checkpoints,
    validate_manual_ec_checkpoint_pair,
)
from .paths import DEFAULT_PROJECT_ROOT, resolve_project_root
from .studies import StudySpec, default_studies
from .validation import AnnDataReport, validate_anndata, validate_source_paths

__all__ = [
    "AnnDataReport",
    "DEFAULT_PROJECT_ROOT",
    "DatasetLoadResult",
    "ManualECCheckpointPaths",
    "ManualECCheckpointSettings",
    "ManualECFilterSettings",
    "Notebook00SourceConfig",
    "PlotConfig",
    "SeuratToAnnDataConverter",
    "StudySpec",
    "cached_h5ad_path",
    "cellbender_output_table",
    "default_studies",
    "expected_cellbender_output_h5",
    "load_dataset_result",
    "load_cached_anndata",
    "load_cached_anndatas",
    "missing_cached_h5ads",
    "resolve_project_root",
    "sample_table",
    "save_manual_ec_checkpoints",
    "source_availability_table",
    "summarize_source_availability",
    "validate_anndata",
    "validate_manual_ec_checkpoint_pair",
    "validate_source_paths",
]
