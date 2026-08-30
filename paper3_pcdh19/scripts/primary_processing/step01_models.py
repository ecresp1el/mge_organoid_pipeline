"""Immutable settings and paths for primary-processing Step 01."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Step01Settings:
    """Define all tunable QC-metric and rendering choices for Step 01."""

    expected_cells: int = 450_788
    expected_genes: int = 19_071
    mitochondrial_prefix: str = "mt-"
    ribosomal_prefixes: tuple[str, ...] = ("Rpl", "Rps")
    percent_top: tuple[int, ...] | None = None
    histogram_bins: int = 80
    plot_dpi: int = 180
    compression: str = "lzf"


@dataclass(frozen=True)
class Step01Paths:
    """Resolve the approved input and run-scoped Step 01 locations."""

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
