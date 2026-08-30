"""Immutable settings and paths for Step 01a MAD sensitivity analysis."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Step01aSettings:
    """Define all tunable, non-filtering Step 01a analysis choices."""

    expected_cells: int = 450_788
    expected_genes: int = 19_071
    expected_samples: int = 12
    sample_field: str = "technical_sample_id"
    stringencies: tuple[float, ...] = (3.0, 4.0, 5.0)
    mad_scale: float = 1.4826
    histogram_bins: int = 80
    plot_dpi: int = 180


@dataclass(frozen=True)
class Step01aPaths:
    """Resolve the exact Step 01 input and run-scoped Step 01a locations."""

    input_h5ad: Path
    input_step_status: Path
    input_output_manifest: Path
    approval_ledger_snapshot: Path
    run_dir: Path
    workflow_root: Path
    frozen_code_dir: Path
    expected_input_run_id: str
    expected_input_bytes: int
    expected_input_sha256: str

