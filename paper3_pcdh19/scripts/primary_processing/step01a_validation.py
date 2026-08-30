"""Validate Step 01 identity, per-sample MAD logic, and the no-filter boundary."""

from __future__ import annotations

import hashlib

import anndata as ad
import numpy as np
import pandas as pd

from .models import ValidationLedger
from .step01a_models import Step01aPaths, Step01aSettings
from .step01a_outliers import Step01aResult


def _file_sha256(path) -> str:
    """Compute a streaming digest for the exact Step 01 H5AD input."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


class Step01aInputValidator:
    """Require the exact in-review Step 01 checkpoint authorized for amendment."""

    def __init__(self, paths: Step01aPaths, ledger: ValidationLedger):
        """Store frozen input evidence and the shared validation ledger."""

        self.paths = paths
        self.ledger = ledger

    def validate(self) -> None:
        """Validate run status, ledger identity, manifest, size, and SHA-256."""

        status = pd.read_csv(self.paths.input_step_status, sep="\t", dtype=str, keep_default_na=False)
        matching = status.loc[status["run_id"] == self.paths.expected_input_run_id]
        observed_status = matching.iloc[0]["status"] if len(matching) == 1 else "missing"
        self.ledger.add("input_step01_status", "input", observed_status == "IN_REVIEW", observed_status, "IN_REVIEW", "Step 01a is an authorized sensitivity amendment to the exact in-review Step 01 run.")
        approval = pd.read_csv(self.paths.approval_ledger_snapshot, sep="\t", dtype=str, keep_default_na=False)
        rows = approval.loc[(approval["run_id"] == self.paths.expected_input_run_id) & (approval["status"] == "IN_REVIEW")]
        self.ledger.add("approval_ledger_identity", "input", len(rows) == 1, len(rows), 1, "Frozen ledger must identify exactly one matching in-review Step 01 run.")
        observed_bytes = self.paths.input_h5ad.stat().st_size
        self.ledger.add("input_bytes", "input", observed_bytes == self.paths.expected_input_bytes, observed_bytes, self.paths.expected_input_bytes, "Step 01 input size must match the frozen contract.")
        manifest = pd.read_csv(self.paths.input_output_manifest, sep="\t", dtype=str, keep_default_na=False)
        row = manifest.loc[manifest["relative_path"] == "objects/pcdh19_step01_qc_metrics.h5ad"]
        manifest_sha = row.iloc[0]["sha256"] if len(row) == 1 else "missing"
        self.ledger.add("input_manifest_identity", "input", manifest_sha == self.paths.expected_input_sha256, manifest_sha, self.paths.expected_input_sha256, "The Step 01 output manifest must identify the exact sensitivity-analysis input.")
        observed_sha = _file_sha256(self.paths.input_h5ad)
        self.ledger.add("input_sha256", "input", observed_sha == self.paths.expected_input_sha256, observed_sha, self.paths.expected_input_sha256, "The Step 01 H5AD must be byte-identical to the frozen input contract.")


class Step01aValidator:
    """Enforce per-sample, one-sided, non-filtering candidate sensitivity logic."""

    REQUIRED_METRICS = ("total_counts", "n_genes_by_counts", "pct_counts_mt")

    def __init__(self, settings: Step01aSettings, ledger: ValidationLedger):
        """Store the exact analysis contract and shared validation ledger."""

        self.settings = settings
        self.ledger = ledger

    def validate_input_object(self, adata: ad.AnnData) -> pd.DataFrame:
        """Validate dimensions/metadata and return only permitted cell metadata."""

        self.ledger.add("input_shape", "object", adata.shape == (self.settings.expected_cells, self.settings.expected_genes), adata.shape, (self.settings.expected_cells, self.settings.expected_genes), "Step 01a cannot change the Step 01 cell/gene universe.")
        missing = [field for field in (self.settings.sample_field, *self.REQUIRED_METRICS) if field not in adata.obs]
        self.ledger.add("required_obs_fields", "object", not missing, missing or "none", "none", "All Step 01 metrics and the technical-sample ID must be present.")
        samples = adata.obs[self.settings.sample_field].astype(str)
        self.ledger.add("technical_sample_count", "object", samples.nunique() == self.settings.expected_samples, samples.nunique(), self.settings.expected_samples, "Candidate boundaries must be calculated across exactly 12 technical samples.")
        finite = all(np.isfinite(adata.obs[field].to_numpy(dtype=float)).all() for field in self.REQUIRED_METRICS)
        self.ledger.add("finite_metrics", "object", finite, finite, True, "All input QC metrics must be finite.")
        self.ledger.add("input_read_only_backed", "object", bool(adata.isbacked), adata.isbacked, True, "Step 01a reads the H5AD in backed read-only mode and writes no replacement object.")
        return adata.obs[[self.settings.sample_field, *self.REQUIRED_METRICS]].copy()

    def validate_results(self, result: Step01aResult) -> None:
        """Validate sample independence, one-sided logic, monotonicity, and joins."""

        expected_boundary_rows = self.settings.expected_samples * len(self.settings.stringencies) * 3
        self.ledger.add("boundary_row_count", "results", len(result.boundaries) == expected_boundary_rows, len(result.boundaries), expected_boundary_rows, "Every sample, stringency, and metric requires one boundary.")
        self.ledger.add("per_cell_row_count", "results", len(result.per_cell) == self.settings.expected_cells, len(result.per_cell), self.settings.expected_cells, "Every Step 01 cell must retain one candidate-sensitivity record.")
        self.ledger.add("unique_cell_ids", "results", result.per_cell["cell_id"].is_unique, result.per_cell["cell_id"].nunique(), self.settings.expected_cells, "Candidate records must remain one-to-one with cells.")
        self.ledger.add("no_design_group_definition", "results", "design_group" not in result.per_cell.columns and "design_group" not in result.boundaries.columns, "absent", "absent", "Genotype/design groups cannot influence Step 01a boundaries.")
        forbidden = [column for column in result.per_cell if column.startswith(("high_total_counts", "high_n_genes"))]
        self.ledger.add("no_upper_complexity_flags", "results", not forbidden, forbidden or "none", "none", "High counts and high detected genes are reserved for later Scrublet assessment.")
        for sample, frame in result.per_cell.groupby(self.settings.sample_field, sort=False, observed=True):
            for base in ("low_total_counts", "low_n_genes_by_counts", "high_pct_counts_mt", "any_candidate"):
                flags = [frame[f"{base}_{value:g}mad"].to_numpy(dtype=bool) for value in self.settings.stringencies]
                monotonic = all(np.logical_or(~later, earlier).all() for earlier, later in zip(flags, flags[1:]))
                self.ledger.add(f"monotonic_flags_{sample}_{base}", "results", monotonic, monotonic, True, "More stringent MAD boundaries must flag a subset of cells.")
        for stringency in self.settings.stringencies:
            suffix = f"{stringency:g}mad"
            individuals = result.per_cell[[f"low_total_counts_{suffix}", f"low_n_genes_by_counts_{suffix}", f"high_pct_counts_mt_{suffix}"]]
            union_matches = result.per_cell[f"any_candidate_{suffix}"].equals(individuals.any(axis=1))
            intersection_matches = result.per_cell[f"all_three_candidate_{suffix}"].equals(individuals.all(axis=1))
            self.ledger.add(f"joint_union_{suffix}", "results", union_matches, union_matches, True, "Any-candidate must equal the union of the three one-sided criteria.")
            self.ledger.add(f"joint_intersection_{suffix}", "results", intersection_matches, intersection_matches, True, "All-three must equal the intersection of the three criteria.")
        summary_totals = result.overlap_summary.groupby(["technical_sample_id", "stringency_mad"], observed=True)["n_cells"].sum()
        expected = result.per_cell.groupby(self.settings.sample_field, observed=True).size()
        overlap_complete = all(int(value) == int(expected.loc[sample]) for (sample, _), value in summary_totals.items())
        self.ledger.add("overlap_partition_complete", "results", overlap_complete, overlap_complete, True, "Exact overlap patterns must partition every sample at every stringency.")
