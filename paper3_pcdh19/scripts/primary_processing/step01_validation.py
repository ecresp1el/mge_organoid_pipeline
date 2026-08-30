"""Input approval, count invariance, and Step 01 scientific-boundary checks."""

from __future__ import annotations

import hashlib
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from .models import ValidationLedger
from .step01_metrics import QCFeatureSetResult
from .step01_models import Step01Paths, Step01Settings


def _file_sha256(path: Path) -> str:
    """Compute a streaming digest for the approved input checkpoint."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


class ApprovedStep00Validator:
    """Require the exact user-approved Step 00 checkpoint and manifest identity."""

    def __init__(self, paths: Step01Paths, ledger: ValidationLedger):
        """Store frozen approval evidence and the shared validation ledger."""

        self.paths = paths
        self.ledger = ledger

    def validate(self) -> None:
        """Validate status, run ID, size, manifest, and SHA-256 before loading."""

        status = pd.read_csv(self.paths.input_step_status, sep="\t", dtype=str, keep_default_na=False)
        matching = status.loc[status["run_id"] == self.paths.expected_input_run_id]
        approved = len(matching) == 1 and matching.iloc[0]["status"] == "APPROVED"
        self.ledger.add("approved_step00_status", "input", approved, matching.iloc[0]["status"] if len(matching) == 1 else "missing", "APPROVED", "Only the explicitly approved Step 00 run may feed Step 01.")
        approval_ledger = pd.read_csv(self.paths.approval_ledger_snapshot, sep="\t", dtype=str, keep_default_na=False)
        approved_rows = approval_ledger.loc[
            (approval_ledger["run_id"] == self.paths.expected_input_run_id)
            & (approval_ledger["status"] == "APPROVED")
            & (approval_ledger["approved_run_id"] == self.paths.expected_input_run_id)
        ]
        self.ledger.add("approval_ledger_identity", "input", len(approved_rows) == 1, len(approved_rows), 1, "Frozen approval-ledger snapshot must identify this exact approved Step 00 run.")
        observed_bytes = self.paths.input_h5ad.stat().st_size
        self.ledger.add("input_bytes", "input", observed_bytes == self.paths.expected_input_bytes, observed_bytes, self.paths.expected_input_bytes, "Approved checkpoint byte size must match the frozen contract.")
        manifest = pd.read_csv(self.paths.input_output_manifest, sep="\t", dtype=str, keep_default_na=False)
        manifest_row = manifest.loc[manifest["relative_path"] == "objects/pcdh19_step00_canonical_raw_counts.h5ad"]
        manifest_matches = len(manifest_row) == 1 and manifest_row.iloc[0]["sha256"] == self.paths.expected_input_sha256
        self.ledger.add("input_manifest_identity", "input", manifest_matches, manifest_row.iloc[0]["sha256"] if len(manifest_row) == 1 else "missing", self.paths.expected_input_sha256, "Frozen Step 00 manifest must identify the approved checkpoint.")
        observed_sha256 = _file_sha256(self.paths.input_h5ad)
        self.ledger.add("input_sha256", "input", observed_sha256 == self.paths.expected_input_sha256, observed_sha256, self.paths.expected_input_sha256, "Approved H5AD content must be byte-identical before Step 01.")


class H5ADMatrixFingerprint:
    """Hash logical CSR arrays to prove that Step 01 does not change counts."""

    DATASETS = ("data", "indices", "indptr")

    def calculate(self, path: Path) -> str:
        """Return a chunked digest of H5AD shape, dtype, data, indices, and indptr."""

        digest = hashlib.sha256()
        with h5py.File(path, "r") as handle:
            matrix = handle["X"]
            digest.update(str(tuple(matrix.attrs.get("shape", ()))).encode("utf-8"))
            for name in self.DATASETS:
                dataset = matrix[name]
                digest.update(name.encode("utf-8"))
                digest.update(str(dataset.dtype).encode("utf-8"))
                step = max(1, 8 * 1024 * 1024 // dataset.dtype.itemsize)
                for start in range(0, len(dataset), step):
                    digest.update(dataset[start : start + step].tobytes())
        return digest.hexdigest()


class Step01Validator:
    """Validate QC fields while enforcing the no-filter/no-transformation boundary."""

    REQUIRED_OBS_METRICS = (
        "total_counts",
        "log1p_total_counts",
        "n_genes_by_counts",
        "log1p_n_genes_by_counts",
        "total_counts_mt",
        "log1p_total_counts_mt",
        "pct_counts_mt",
    )

    def __init__(self, settings: Step01Settings, ledger: ValidationLedger):
        """Store dimensions and scientific boundaries for all Step 01 checks."""

        self.settings = settings
        self.ledger = ledger

    def validate_input_object(self, adata: ad.AnnData) -> None:
        """Require the approved raw-count state before Scanpy calculation."""

        self.ledger.add("input_shape", "object", adata.shape == (self.settings.expected_cells, self.settings.expected_genes), adata.shape, (self.settings.expected_cells, self.settings.expected_genes), "Step 01 must start with every approved cell and gene.")
        self.ledger.add("input_sparse", "object", sparse.issparse(adata.X), type(adata.X).__name__, "SciPy sparse", "Approved raw counts must be sparse.")
        self.ledger.add("input_integer", "object", np.issubdtype(adata.X.dtype, np.integer), adata.X.dtype, "integer", "Approved raw counts must remain integer-valued.")
        self.ledger.add("input_no_layers", "object", len(adata.layers) == 0, list(adata.layers.keys()), [], "Step 00 contains no derived expression layer.")
        self.ledger.add("input_no_reductions", "object", len(adata.obsm) == 0 and len(adata.obsp) == 0, sorted(set(adata.obsm.keys()) | set(adata.obsp.keys())), [], "Step 01 cannot inherit or create reductions/graphs.")

    def validate_metrics(self, adata: ad.AnnData, feature_sets: QCFeatureSetResult) -> None:
        """Require valid Scanpy fields and confirm that no filtering state exists."""

        missing = [field for field in self.REQUIRED_OBS_METRICS if field not in adata.obs]
        self.ledger.add("required_qc_metrics", "object", not missing, missing or "none", "none missing", "Scanpy must create all requested standard and mitochondrial fields.")
        for metric in self.REQUIRED_OBS_METRICS:
            if metric not in adata.obs:
                continue
            values = adata.obs[metric].to_numpy(dtype=float)
            self.ledger.add(f"finite_{metric}", "object", bool(np.isfinite(values).all()), int(np.isfinite(values).sum()), len(values), f"{metric} must be finite for every retained cell.")
            self.ledger.add(f"nonnegative_{metric}", "object", bool((values >= 0).all()), float(values.min()), ">=0", f"{metric} cannot be negative.")
        pct = adata.obs["pct_counts_mt"].to_numpy(dtype=float)
        self.ledger.add("mitochondrial_percentage_range", "object", bool((pct <= 100).all()), float(pct.max()), "<=100", "Mitochondrial percentage must lie between 0 and 100.")
        self.ledger.add("mitochondrial_feature_count", "features", len(feature_sets.mitochondrial_genes) == 13, len(feature_sets.mitochondrial_genes), 13, "The delivered Flex panel contains 13 mt- protein-coding genes.")
        self.ledger.add("ribosomal_panel_limitation", "features", len(feature_sets.ribosomal_genes) == 0, len(feature_sets.ribosomal_genes), 0, "No Rpl/Rps genes are present; ribosomal fraction is intentionally not fabricated.")
        forbidden = [column for column in adata.obs.columns if any(token in column.lower() for token in ("qc_pass", "qc_fail", "exclude", "filter_reason"))]
        self.ledger.add("no_filter_columns", "object", not forbidden, forbidden or "none", "none", "Step 01 cannot define pass/fail state or exclusion reasons.")
        self.ledger.add("output_shape", "object", adata.shape == (self.settings.expected_cells, self.settings.expected_genes), adata.shape, (self.settings.expected_cells, self.settings.expected_genes), "QC calculation cannot remove cells or genes.")
        self.ledger.add("output_integer_counts", "object", np.issubdtype(adata.X.dtype, np.integer), adata.X.dtype, "integer", "Expression counts must remain unchanged raw integers.")
        self.ledger.add("output_no_layers_raw", "object", len(adata.layers) == 0 and adata.raw is None, {"layers": list(adata.layers.keys()), "raw": adata.raw is not None}, {"layers": [], "raw": False}, "Step 01 adds metadata only, not expression representations.")
        self.ledger.add("output_no_reductions", "object", len(adata.obsm) == 0 and len(adata.obsp) == 0, sorted(set(adata.obsm.keys()) | set(adata.obsp.keys())), [], "Step 01 cannot create PCA, UMAP, neighbors, or graphs.")
