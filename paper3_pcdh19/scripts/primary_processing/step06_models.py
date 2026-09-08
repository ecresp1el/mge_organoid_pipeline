"""Define the data contracts passed between Step 06 components.

``Step06Settings`` contains environment/CLI-controlled settings, not every
literal algorithm choice. ``Step06Paths`` identifies the approved input and
run destination. ``Step06Artifacts`` carries arrays from analysis;
``Step06Results`` carries tables and the heuristic report headline.

Production values flow from config/primary_processing_step06.env through the
SLURM wrapper and CLI. Dataclass defaults apply to direct Python construction;
notably n_jobs defaults to 1 here but the production environment requests 32.
See STEP06_CODE_AND_TUNING_GUIDE.md for configurable and code-only choices.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


@dataclass(frozen=True)
class Step06Settings:
    """Store environment/CLI-controlled scientific and rendering settings.

    Notes
    -----
    Frozen values are passed unchanged between components. For each field, see the inline
    comment and environment parameter table in the tuning guide. Identity/dimension guards
    are not scientific tuning knobs. Additional hard-coded choices live in analysis and
    metrics.
    """

    # Identity guard: all 446,349 approved Step 02 cells must survive. Do not tune.
    expected_cells: int = 446_349
    # Identity guard: all 19,071 approved genes must survive. Do not tune.
    expected_genes: int = 19_071
    # Grouping key for HVGs, pseudobulk, mixing and sample summaries; changing it changes the question.
    sample_field: str = "technical_sample_id"
    # Registered genotype metadata used for descriptive grouping and colors, not PCA predictors.
    genotype_field: str = "genotype"
    # Registered sex metadata used for descriptive grouping and colors, not PCA predictors.
    sex_field: str = "sex"
    # Registered combined design group used for display; not an integration covariate.
    design_field: str = "design_group"
    # Per-cell count total before natural log1p; changes all downstream expression-based results.
    target_sum: float = 10_000.0
    # Called on RAW counts; seurat_v3 expects counts. Other flavors are not drop-in replacements.
    hvg_flavor: str = "seurat_v3"
    # Number of HVGs defining PCA and HVG-only pseudobulk similarity; all genes remain in output.
    n_top_genes: int = 3_000
    # Number of centered PCs fitted and supplied to neighbors; must be below n_top_genes.
    pca_components: int = 50
    # First PCs used for sample-center distances; affects review flags, not the neighbor graph.
    centroid_components: int = 20
    # Graph locality: larger neighborhoods average over broader structure; affects mixing and clustering.
    n_neighbors: int = 30
    # Scanpy worker setting. Production env supplies 32; does not guarantee every algorithm uses 32 workers.
    n_jobs: int = 1
    # Distance between cells in PCA space; changing metric changes graph, mixing, UMAP and Leiden.
    neighbor_metric: str = "cosine"
    # UMAP packing parameter; changes display geometry and UMAP displacement flags, not graph-based entropy.
    umap_min_dist: float = 0.3
    # UMAP scale parameter; affects layout and UMAP center flags, not the fitted PCA or graph.
    umap_spread: float = 1.0
    # Cluster granularity; changes marker-winner states and composition summaries, not graph construction.
    leiden_resolution: float = 1.0
    # Reproducibility control passed to PCA/graph/UMAP/Leiden/rendering, not a fit-quality knob.
    random_seed: int = 20_260_903
    # Display-only cap split across samples; scientific computations still use all cells.
    render_max_cells: int = 120_000
    # PNG rendering resolution; no effect on coordinates, scores, or decisions.
    plot_dpi: int = 180
    # H5AD storage codec; affects file size/time, not raw count values.
    compression: str = "lzf"


@dataclass(frozen=True)
class Step06Paths:
    """Resolve approved inputs, bypass evidence, and run-scoped outputs.

    Notes
    -----
    All paths are supplied explicitly; no latest-run discovery occurs.
    expected_input_bytes/SHA-256 and expected_step02_run_id bind the run to the approved
    checkpoint. workflow_root holds the shared approval ledger; frozen_code_dir is audited
    for provenance.
    """

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
    """Carry full-cell analysis products into metrics and publication.

    Notes
    -----
    Cell arrays and graphs retain input row order. annotations is cell-indexed; pca is
    cells x PCs, umap cells x 2, graphs cells x cells. hvg_table covers all input genes;
    pseudobulk_logcpm is samples x selected HVGs. rendering_cells is only a display index.
    """

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
    """Bundle companion tables and the provisional diagnostic conclusion.

    Notes
    -----
    tables maps output filenames to DataFrames. outcome_code/label/rationale are the
    automatic heuristic from Step06MetricBuilder._outcome, not approval or a validated
    biological classifier.
    """

    tables: dict[str, pd.DataFrame]
    outcome_code: int
    outcome_label: str
    outcome_rationale: str
