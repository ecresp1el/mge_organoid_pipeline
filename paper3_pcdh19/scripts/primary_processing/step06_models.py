"""Immutable contracts and result containers for Step 06 diagnostics."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


@dataclass(frozen=True)
class Step06Settings:
    """Define every scientific and rendering parameter for Step 06."""

    expected_cells: int = 446_349
    expected_genes: int = 19_071
    sample_field: str = "technical_sample_id"
    genotype_field: str = "genotype"
    sex_field: str = "sex"
    design_field: str = "design_group"
    target_sum: float = 10_000.0
    hvg_flavor: str = "seurat_v3"
    n_top_genes: int = 3_000
    pca_components: int = 50
    centroid_components: int = 20
    n_neighbors: int = 30
    n_jobs: int = 1
    neighbor_metric: str = "cosine"
    umap_min_dist: float = 0.3
    umap_spread: float = 1.0
    leiden_resolution: float = 1.0
    random_seed: int = 20_260_903
    render_max_cells: int = 120_000
    plot_dpi: int = 180
    compression: str = "lzf"


@dataclass(frozen=True)
class Step06Paths:
    """Resolve approved inputs, bypass evidence, and run-scoped outputs."""

    input_h5ad: Path
    input_status: Path
    input_manifest: Path
    approval_ledger_snapshot: Path
    bypass_decision: Path
    run_dir: Path
    workflow_root: Path
    frozen_code_dir: Path
    expected_step02_run_id: str
    expected_input_bytes: int
    expected_input_sha256: str


@dataclass
class Step06Artifacts:
    """Carry full-cell analysis products into metrics and publication."""

    annotations: pd.DataFrame
    pca: np.ndarray
    umap: np.ndarray
    connectivities: sparse.csr_matrix
    distances: sparse.csr_matrix
    hvg_table: pd.DataFrame
    variance_ratio: np.ndarray
    cluster_program_scores: pd.DataFrame
    pseudobulk_logcpm: pd.DataFrame
    rendering_cells: pd.DataFrame


@dataclass
class Step06Results:
    """Bundle companion tables and the provisional diagnostic conclusion."""

    tables: dict[str, pd.DataFrame]
    outcome_code: int
    outcome_label: str
    outcome_rationale: str
