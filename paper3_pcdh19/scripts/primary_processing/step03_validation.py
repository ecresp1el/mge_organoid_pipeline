"""Validate Step 03 approvals, matrices, capture semantics, and R results."""

from __future__ import annotations

import hashlib
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from .models import ValidationLedger
from .step02_validation import Step02MatrixFingerprint
from .step03_models import Step03Paths, Step03Settings


def sha256_file(path: Path) -> str:
    """Calculate one file's streaming SHA-256 digest."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


class Step03InputValidator:
    """Require the exact approved Step 02 checkpoint and one-capture policy."""

    def __init__(self, paths: Step03Paths, settings: Step03Settings, ledger: ValidationLedger):
        """Store the frozen paths, settings, and shared validation ledger."""

        self.paths = paths
        self.settings = settings
        self.ledger = ledger

    def validate_approval(self) -> None:
        """Validate run-local and workflow-level Step 02 approval evidence."""

        status = pd.read_csv(self.paths.input_status, sep="\t", dtype=str, keep_default_na=False)
        rows = status.loc[
            (status["run_id"] == self.paths.expected_step02_run_id)
            & (status["status"] == "APPROVED")
            & (status["approved_run_id"] == self.paths.expected_step02_run_id)
        ]
        self.ledger.add("approved_step02_status", "input", len(rows) == 1, len(rows), 1, "Step 03 may consume only the exact approved Step 02 run.")
        approval = pd.read_csv(self.paths.approval_ledger_snapshot, sep="\t", dtype=str, keep_default_na=False)
        rows = approval.loc[
            (approval["run_id"] == self.paths.expected_step02_run_id)
            & (approval["status"] == "APPROVED")
            & (approval["approved_run_id"] == self.paths.expected_step02_run_id)
        ]
        self.ledger.add("approved_step02_ledger", "input", len(rows) == 1, len(rows), 1, "The frozen workflow ledger must independently record the same approval.")
        manifest = pd.read_csv(self.paths.input_manifest, sep="\t", dtype=str, keep_default_na=False)
        rows = manifest.loc[manifest["relative_path"] == "objects/pcdh19_step02_qc_filtered.h5ad"]
        manifest_sha = rows.iloc[0]["sha256"] if len(rows) == 1 else "missing"
        self.ledger.add("manifest_input_sha256", "input", manifest_sha == self.paths.expected_input_sha256, manifest_sha, self.paths.expected_input_sha256, "The upstream manifest must identify the approved checkpoint content.")
        observed_bytes = self.paths.input_h5ad.stat().st_size
        observed_sha = sha256_file(self.paths.input_h5ad)
        self.ledger.add("input_bytes", "input", observed_bytes == self.paths.expected_input_bytes, observed_bytes, self.paths.expected_input_bytes, "The approved input byte size must match.")
        self.ledger.add("input_sha256", "input", observed_sha == self.paths.expected_input_sha256, observed_sha, self.paths.expected_input_sha256, "The approved input must be byte-identical.")

    def validate_object(self, adata: ad.AnnData) -> None:
        """Require raw counts, expected dimensions, metadata, and no prior reduction."""

        required = {self.settings.sample_field, self.settings.design_field, "step02_qc_pass"}
        missing = sorted(required - set(adata.obs.columns))
        self.ledger.add("required_metadata", "object", not missing, missing or "none", "none", "Sample, design, and Step 02 disposition metadata must be present.")
        self.ledger.add("input_shape", "object", adata.shape == (self.settings.expected_cells, self.settings.expected_genes), adata.shape, (self.settings.expected_cells, self.settings.expected_genes), "Step 03 must retain every approved Step 02 cell and gene.")
        sparse_state = sparse.issparse(adata.X) or getattr(adata.X, "format", None) in {"csr", "csc"}
        raw_state = sparse_state and np.issubdtype(adata.X.dtype, np.integer)
        self.ledger.add("raw_sparse_integer_counts", "object", raw_state, {"sparse": sparse_state, "dtype": str(adata.X.dtype)}, {"sparse": True, "dtype": "integer"}, "scDblFinder must receive sparse raw integer gene counts.")
        empty_derived = len(adata.layers) == 0 and adata.raw is None and len(adata.obsm) == 0 and len(adata.obsp) == 0
        self.ledger.add("no_prior_analysis_representation", "object", empty_derived, {"layers": list(adata.layers), "raw": adata.raw is not None, "obsm": list(adata.obsm), "obsp": list(adata.obsp)}, "all empty", "No prior normalization, embedding, or graph may enter Step 03.")
        all_pass = bool(adata.obs["step02_qc_pass"].astype(bool).all())
        self.ledger.add("all_cells_pass_step02", "object", all_pass, int(adata.obs["step02_qc_pass"].astype(bool).sum()), adata.n_obs, "Every Step 03 cell must be retained by the approved Step 02 decision.")
        self.ledger.add("twelve_technical_samples", "object", adata.obs[self.settings.sample_field].nunique() == 12, adata.obs[self.settings.sample_field].nunique(), 12, "The 12 probe-barcode technical samples remain reporting metadata.")


class Step03BridgeValidator:
    """Validate lossless reuse of H5AD CSR arrays as transposed 10x CSC arrays."""

    def __init__(self, paths: Step03Paths, settings: Step03Settings, ledger: ValidationLedger):
        """Store bridge paths, expected dimensions, and the validation ledger."""

        self.paths = paths
        self.settings = settings
        self.ledger = ledger

    def validate(self) -> None:
        """Require exact array identity, dimensions, names, and one capture."""

        with h5py.File(self.paths.input_h5ad, "r") as source, h5py.File(self.paths.bridge_h5, "r") as target:
            source_x = source["X"]
            target_x = target["matrix"]
            expected_shape = (self.settings.expected_genes, self.settings.expected_cells)
            observed_shape = tuple(int(value) for value in target_x["shape"][:])
            self.ledger.add("bridge_shape", "bridge", observed_shape == expected_shape, observed_shape, expected_shape, "The R bridge must be genes by cells.")
            for name in ("data", "indices", "indptr"):
                source_dataset = source_x[name]
                target_dataset = target_x[name]
                same_shape = source_dataset.shape == target_dataset.shape
                same_dtype = source_dataset.dtype == target_dataset.dtype
                same_digest = self._dataset_digest(source_dataset) == self._dataset_digest(target_dataset)
                self.ledger.add(f"bridge_{name}_identity", "bridge", same_shape and same_dtype and same_digest, {"shape": target_dataset.shape, "dtype": str(target_dataset.dtype), "digest_equal": same_digest}, {"shape": source_dataset.shape, "dtype": str(source_dataset.dtype), "digest_equal": True}, "The bridge must reuse the exact sparse array without numerical conversion.")
            captures = pd.read_csv(self.paths.cell_metadata, sep="\t", usecols=["capture_id"])["capture_id"]
            self.ledger.add("one_capture", "bridge", captures.nunique() == 1 and captures.iloc[0] == self.settings.capture_id, sorted(captures.unique()), [self.settings.capture_id], "All cells must be passed as one independently processed GEX_1 capture.")

    @staticmethod
    def _dataset_digest(dataset: h5py.Dataset) -> str:
        """Hash one HDF5 numeric dataset without holding a duplicate in memory."""

        digest = hashlib.sha256(str(dataset.dtype).encode())
        step = max(1, 8 * 1024 * 1024 // dataset.dtype.itemsize)
        for start in range(0, len(dataset), step):
            digest.update(dataset[start : start + step].tobytes())
        return digest.hexdigest()


class Step03ResultValidator:
    """Require complete, aligned, non-filtering scDblFinder results."""

    def __init__(self, settings: Step03Settings, ledger: ValidationLedger):
        """Store the approved dimensions and shared validation ledger."""

        self.settings = settings
        self.ledger = ledger

    def validate(self, adata: ad.AnnData, results: pd.DataFrame, pca_ids: pd.Index, pca: np.ndarray, original_fingerprint: str) -> None:
        """Validate cell identity, scores, calls, diagnostic PCA, and raw counts."""

        required = {"cell_id", "primary_score", "primary_class", "primary_cluster", "replicate_score", "replicate_class"}
        missing = sorted(required - set(results.columns))
        self.ledger.add("r_result_fields", "result", not missing, missing or "none", "none", "Primary and reproducibility scores/calls plus generated clusters must be present.")
        result_ids = pd.Index(results["cell_id"].astype(str))
        object_ids = adata.obs_names.astype(str)
        self.ledger.add("r_result_rows", "result", len(results) == self.settings.expected_cells, len(results), self.settings.expected_cells, "Native R must return every input cell exactly once.")
        self.ledger.add("r_result_ids_unique", "result", result_ids.is_unique, result_ids.nunique(), self.settings.expected_cells, "Native R result cell IDs must be unique.")
        self.ledger.add("r_result_order", "result", result_ids.equals(object_ids), "exact" if result_ids.equals(object_ids) else "different", "exact", "Native R result order must match the approved checkpoint.")
        for field in ("primary_score", "replicate_score"):
            values = results[field].to_numpy(dtype=float)
            valid = bool(np.isfinite(values).all() and np.logical_and(values >= 0, values <= 1).all())
            self.ledger.add(f"{field}_range", "result", valid, {"min": float(np.nanmin(values)), "max": float(np.nanmax(values))}, "finite [0,1]", "scDblFinder scores must be finite probabilities/scores in the documented range.")
        for field in ("primary_class", "replicate_class"):
            observed = sorted(results[field].astype(str).unique())
            self.ledger.add(f"{field}_values", "result", set(observed) == {"singlet", "doublet"}, observed, ["doublet", "singlet"], "Both scDblFinder call classes must be represented and no undocumented call value is allowed.")
        self.ledger.add("pca_rows", "diagnostic", pca.shape[0] == self.settings.expected_cells, pca.shape, f"({self.settings.expected_cells}, dimensions)", "The preserved internal PCA must contain every real cell and no artificial doublets.")
        self.ledger.add("pca_ids", "diagnostic", pca_ids.equals(object_ids), "exact" if pca_ids.equals(object_ids) else "different", "exact", "Internal PCA rows must align exactly to the approved cells.")
        self.ledger.add("pca_finite", "diagnostic", bool(np.isfinite(pca).all()), int(np.isfinite(pca).sum()), pca.size, "Internal diagnostic PCA coordinates must be finite.")
        observed_fingerprint = Step02MatrixFingerprint.memory(adata.X)
        self.ledger.add("raw_matrix_unchanged", "object", observed_fingerprint == original_fingerprint, observed_fingerprint, original_fingerprint, "Attaching calls and diagnostic PCA must not alter a single raw count.")
        self.ledger.add("no_cell_or_gene_removal", "object", adata.shape == (self.settings.expected_cells, self.settings.expected_genes), adata.shape, (self.settings.expected_cells, self.settings.expected_genes), "Step 03 is detection-only and cannot remove cells or genes.")
