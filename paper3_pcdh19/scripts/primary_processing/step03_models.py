"""Immutable contracts for Step 03 scDblFinder detection and review."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Step03Settings:
    """Define the approved capture, scDblFinder, and checkpoint contract."""

    expected_cells: int = 446_349
    expected_genes: int = 19_071
    capture_id: str = "GEX_1"
    sample_field: str = "technical_sample_id"
    design_field: str = "design_group"
    primary_seed: int = 20_260_830
    reproducibility_seed: int = 20_260_831
    plot_dpi: int = 180
    compression: str = "lzf"


@dataclass(frozen=True)
class Step03Paths:
    """Resolve the exact approved input and run-scoped Step 03 assets."""

    input_h5ad: Path
    input_status: Path
    input_manifest: Path
    approval_ledger_snapshot: Path
    run_dir: Path
    workflow_root: Path
    frozen_code_dir: Path
    expected_step02_run_id: str
    expected_input_bytes: int
    expected_input_sha256: str

    @property
    def intermediate_dir(self) -> Path:
        """Return the retained run-scoped cross-language intermediate directory."""

        return self.run_dir / "intermediate"

    @property
    def bridge_h5(self) -> Path:
        """Return the lossless 10x-schema bridge file consumed by native R."""

        return self.intermediate_dir / "step02_counts_for_scdblfinder.h5"

    @property
    def cell_metadata(self) -> Path:
        """Return the exact ordered metadata table consumed by native R."""

        return self.intermediate_dir / "step03_cell_metadata.tsv"

    @property
    def r_output_dir(self) -> Path:
        """Return the native-R output directory consumed by Python finalization."""

        return self.intermediate_dir / "r_output"
