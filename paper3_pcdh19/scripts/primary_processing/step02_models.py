"""Immutable settings and paths for approved Step 02 QC filtering."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Step02Settings:
    """Define the selected filtering policy and expected output dimensions."""

    expected_cells_before: int = 450_788
    expected_cells_after: int = 446_349
    expected_removed: int = 4_439
    expected_genes: int = 19_071
    selected_stringency: int = 5
    sample_field: str = "technical_sample_id"
    design_field: str = "design_group"
    plot_dpi: int = 180
    compression: str = "lzf"


@dataclass(frozen=True)
class Step02Paths:
    """Resolve approved inputs, registered evidence, and run-scoped outputs."""

    input_h5ad: Path
    input_step01_status: Path
    input_step01_manifest: Path
    input_flags: Path
    input_step01a_status: Path
    input_step01a_manifest: Path
    approval_ledger_snapshot: Path
    run_dir: Path
    workflow_root: Path
    frozen_code_dir: Path
    expected_step01_run_id: str
    expected_step01_bytes: int
    expected_step01_sha256: str
    expected_step01a_run_id: str
    expected_flags_bytes: int
    expected_flags_sha256: str
