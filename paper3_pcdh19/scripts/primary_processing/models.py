"""Immutable settings, paths, and validation records for primary processing."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import pandas as pd


@dataclass(frozen=True)
class Step00Settings:
    """Define the complete, tunable scientific contract for Step 00."""

    expected_samples: int = 12
    expected_cells: int = 450_788
    expected_genes: int = 19_071
    expected_genome: str = "GRCm39"
    feature_type: str = "Gene Expression"
    compression: str = "lzf"


@dataclass(frozen=True)
class Step00Paths:
    """Resolve all read-only inputs and run-scoped publication locations."""

    cellranger_root: Path
    sample_key: Path
    technical_manifest: Path
    run_dir: Path
    workflow_root: Path
    frozen_code_dir: Path

    def sample_directory(self, sample_id: str) -> Path:
        """Return the delivered per-sample Cell Ranger directory."""

        return self.cellranger_root / "per_sample_outs" / sample_id

    def filtered_matrix(self, sample_id: str) -> Path:
        """Return the canonical called-cell gene-expression H5 path."""

        return self.sample_directory(sample_id) / "sample_filtered_feature_bc_matrix.h5"

    def raw_matrix(self, sample_id: str) -> Path:
        """Return the unfiltered gene-expression H5 registered for Step 04."""

        return self.sample_directory(sample_id) / "sample_raw_feature_bc_matrix.h5"

    def filtered_barcodes(self, sample_id: str) -> Path:
        """Return the vendor filtered-barcode CSV used for exact validation."""

        return self.sample_directory(sample_id) / "sample_filtered_barcodes.csv"

    def metrics_summary(self, sample_id: str) -> Path:
        """Return the sample-level Cell Ranger metrics table."""

        return self.sample_directory(sample_id) / "metrics_summary.csv"


@dataclass(frozen=True)
class SampleStructure:
    """Describe one Cell Ranger H5 without loading its expression values."""

    sample_id: str
    filtered_path: Path
    raw_path: Path
    barcode_csv_path: Path
    metrics_path: Path
    n_cells: int
    n_features: int
    nnz: int
    count_dtype: str
    barcodes: tuple[str, ...]
    gene_ids: tuple[str, ...]
    gene_symbols: tuple[str, ...]
    feature_types: tuple[str, ...]
    genomes: tuple[str, ...]
    raw_n_barcodes: int
    raw_n_features: int
    raw_n_total_features: int
    raw_feature_types: tuple[str, ...]
    feature_signature: str
    raw_feature_signature: str
    raw_contains_filtered_features: bool
    raw_preserves_filtered_feature_order: bool
    vendor_metrics_cells: int


@dataclass(frozen=True)
class ValidationCheck:
    """Represent one explicit, machine-readable validation assertion."""

    check_id: str
    scope: str
    status: str
    observed: str
    expected: str
    message: str


@dataclass
class ValidationLedger:
    """Collect validation checks and raise once all diagnostics are recorded."""

    checks: list[ValidationCheck] = field(default_factory=list)

    def add(
        self,
        check_id: str,
        scope: str,
        passed: bool,
        observed: Any,
        expected: Any,
        message: str,
    ) -> None:
        """Append a PASS or FAIL record without stopping the current audit."""

        self.checks.append(
            ValidationCheck(
                check_id=check_id,
                scope=scope,
                status="PASS" if passed else "FAIL",
                observed=str(observed),
                expected=str(expected),
                message=message,
            )
        )

    def extend(self, checks: Iterable[ValidationCheck]) -> None:
        """Append an iterable of already constructed validation records."""

        self.checks.extend(checks)

    def to_frame(self) -> pd.DataFrame:
        """Return checks in stable column order for TSV publication."""

        return pd.DataFrame(
            [check.__dict__ for check in self.checks],
            columns=["check_id", "scope", "status", "observed", "expected", "message"],
        )

    def require_all_pass(self) -> None:
        """Raise a concise error after identifying every failed assertion."""

        failed = [check for check in self.checks if check.status != "PASS"]
        if failed:
            identifiers = ", ".join(f"{item.scope}:{item.check_id}" for item in failed)
            raise ValueError(f"Step 00 validation failed: {identifiers}")
