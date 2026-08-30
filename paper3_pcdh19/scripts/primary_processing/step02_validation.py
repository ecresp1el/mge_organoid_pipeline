"""Validate approved Step 02 inputs, filtering logic, and raw-count output."""

from __future__ import annotations

import hashlib
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from .models import ValidationLedger
from .step02_filtering import Step02FilterResult
from .step02_models import Step02Paths, Step02Settings


def _sha256(path: Path) -> str:
    """Calculate a streaming SHA-256 digest for one frozen input."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


class Step02InputValidator:
    """Require exact approved Step 01/01a identities."""

    def __init__(self, paths: Step02Paths, settings: Step02Settings, ledger: ValidationLedger):
        """Store frozen contracts and the shared validation ledger."""

        self.paths = paths
        self.settings = settings
        self.ledger = ledger

    def validate(self) -> None:
        """Validate approvals, manifests, byte sizes, and content hashes."""

        self._approved_status(self.paths.input_step01_status, self.paths.expected_step01_run_id, "step01")
        self._approved_status(self.paths.input_step01a_status, self.paths.expected_step01a_run_id, "step01a")
        approval = pd.read_csv(self.paths.approval_ledger_snapshot, sep="\t", dtype=str, keep_default_na=False)
        for label, run_id in (("step01", self.paths.expected_step01_run_id), ("step01a", self.paths.expected_step01a_run_id)):
            rows = approval.loc[(approval["run_id"] == run_id) & (approval["status"] == "APPROVED") & (approval["approved_run_id"] == run_id)]
            self.ledger.add(f"approval_ledger_{label}", "input", len(rows) == 1, len(rows), 1, f"The frozen ledger must contain exactly one approved {label} run.")
        self._manifest_identity(self.paths.input_step01_manifest, "objects/pcdh19_step01_qc_metrics.h5ad", self.paths.expected_step01_sha256, "step01")
        self._manifest_identity(self.paths.input_step01a_manifest, "tables/per_cell_mad_candidate_flags.tsv.gz", self.paths.expected_flags_sha256, "step01a_flags")
        self._file_identity(self.paths.input_h5ad, self.paths.expected_step01_bytes, self.paths.expected_step01_sha256, "step01_h5ad")
        self._file_identity(self.paths.input_flags, self.paths.expected_flags_bytes, self.paths.expected_flags_sha256, "step01a_flags")

    def _approved_status(self, path: Path, run_id: str, label: str) -> None:
        """Require one exact approved run-local status row."""

        status = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        rows = status.loc[(status["run_id"] == run_id) & (status["status"] == "APPROVED") & (status["approved_run_id"] == run_id)]
        self.ledger.add(f"approved_status_{label}", "input", len(rows) == 1, len(rows), 1, f"Step 02 may consume only the exact approved {label} run.")

    def _manifest_identity(self, path: Path, relative: str, expected_sha: str, label: str) -> None:
        """Require a manifested input digest to match the frozen contract."""

        manifest = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        rows = manifest.loc[manifest["relative_path"] == relative]
        observed = rows.iloc[0]["sha256"] if len(rows) == 1 else "missing"
        self.ledger.add(f"manifest_{label}", "input", observed == expected_sha, observed, expected_sha, f"The {label} input must match its upstream output manifest.")

    def _file_identity(self, path: Path, expected_bytes: int, expected_sha: str, label: str) -> None:
        """Require byte-size and content identity for one exact input file."""

        observed_bytes = path.stat().st_size
        observed_sha = _sha256(path)
        self.ledger.add(f"bytes_{label}", "input", observed_bytes == expected_bytes, observed_bytes, expected_bytes, f"The {label} byte size must match.")
        self.ledger.add(f"sha256_{label}", "input", observed_sha == expected_sha, observed_sha, expected_sha, f"The {label} content must be byte-identical.")


class Step02Validator:
    """Enforce selected-union counts, exclusion provenance, and raw-count state."""

    def __init__(self, settings: Step02Settings, ledger: ValidationLedger):
        """Store approved dimensions and the shared validation ledger."""

        self.settings = settings
        self.ledger = ledger

    def validate_input_object(self, adata: ad.AnnData, flags: pd.DataFrame) -> None:
        """Require the complete approved cell universe and exact flag alignment."""

        required = {
            "cell_id", "technical_sample_id", "total_counts", "n_genes_by_counts",
            "pct_counts_mt", "low_total_counts_5mad",
            "low_n_genes_by_counts_5mad", "high_pct_counts_mt_5mad",
        }
        missing = sorted(required - set(flags.columns))
        self.ledger.add("flag_columns", "flags", not missing, missing or "none", "none", "All approved 5-MAD flag and alignment fields must be present.")
        if missing:
            raise ValueError(f"Missing approved flag columns: {missing}")
        self.ledger.add("input_shape", "object", adata.shape == (self.settings.expected_cells_before, self.settings.expected_genes), adata.shape, (self.settings.expected_cells_before, self.settings.expected_genes), "Step 02 must start from the complete approved Step 01 checkpoint.")
        self.ledger.add("input_sparse_integer", "object", sparse.issparse(adata.X) and np.issubdtype(adata.X.dtype, np.integer), {"sparse": sparse.issparse(adata.X), "dtype": str(adata.X.dtype)}, {"sparse": True, "dtype": "integer"}, "Step 02 must start from sparse integer raw counts.")
        self.ledger.add("input_no_derived_expression", "object", len(adata.layers) == 0 and adata.raw is None and len(adata.obsm) == 0 and len(adata.obsp) == 0, {"layers": list(adata.layers), "raw": adata.raw is not None, "obsm": list(adata.obsm), "obsp": list(adata.obsp)}, "all empty", "No normalization, reductions, or graphs may enter filtering.")
        self.ledger.add("flag_rows", "flags", len(flags) == self.settings.expected_cells_before, len(flags), self.settings.expected_cells_before, "The candidate table must contain every original cell.")
        self.ledger.add("flag_cell_ids_unique", "flags", flags["cell_id"].is_unique, flags["cell_id"].nunique(), self.settings.expected_cells_before, "Candidate flags must remain one-to-one with cells.")
        flag_ids = pd.Index(flags["cell_id"].astype(str))
        object_ids = adata.obs_names.astype(str)
        self.ledger.add("flag_order_exact", "flags", flag_ids.equals(object_ids), "exact" if flag_ids.equals(object_ids) else "different", "exact", "The selected candidate table must preserve exact Step 01 cell order.")
        sample_match = np.array_equal(flags["technical_sample_id"].astype(str).to_numpy(), adata.obs[self.settings.sample_field].astype(str).to_numpy())
        self.ledger.add("flag_sample_alignment", "flags", sample_match, "exact" if sample_match else "different", "exact", "Candidate sample IDs must align to Step 01 metadata cell-by-cell.")
        counts_match = np.array_equal(flags["total_counts"].to_numpy(), adata.obs["total_counts"].to_numpy())
        genes_match = np.array_equal(flags["n_genes_by_counts"].to_numpy(), adata.obs["n_genes_by_counts"].to_numpy())
        mt_match = np.allclose(flags["pct_counts_mt"].to_numpy(dtype=float), adata.obs["pct_counts_mt"].to_numpy(dtype=float), rtol=0, atol=1e-12)
        self.ledger.add("flag_metric_alignment", "flags", counts_match and genes_match and mt_match, {"total_counts": counts_match, "n_genes": genes_match, "pct_mt": mt_match}, "all exact/tolerance matched", "Candidate metrics must align with the approved Step 01 object.")
        forbidden = [column for column in flags if column.startswith(("high_total_counts", "high_n_genes"))]
        self.ledger.add("no_upper_complexity_filter", "flags", not forbidden, forbidden or "none", "none", "Step 02 cannot use high counts or high genes as a doublet proxy.")

    def validate_result(self, result: Step02FilterResult, original_var_names: pd.Index) -> None:
        """Validate exact selected-union removal and complete reason preservation."""

        self.ledger.add("removed_count", "result", int(result.disposition["step02_excluded"].sum()) == self.settings.expected_removed, int(result.disposition["step02_excluded"].sum()), self.settings.expected_removed, "The approved 5-MAD union must remove exactly the reviewed candidate count.")
        self.ledger.add("filtered_shape", "result", result.filtered.shape == (self.settings.expected_cells_after, self.settings.expected_genes), result.filtered.shape, (self.settings.expected_cells_after, self.settings.expected_genes), "Step 02 must remove selected cells and no genes.")
        self.ledger.add("gene_order_preserved", "result", result.filtered.var_names.equals(original_var_names), result.filtered.n_vars, len(original_var_names), "All genes and exact order must be preserved.")
        excluded = result.disposition.loc[result.disposition["step02_excluded"]]
        retained = result.disposition.loc[~result.disposition["step02_excluded"]]
        self.ledger.add("excluded_reasons_complete", "result", bool(excluded["step02_exclusion_reasons"].str.len().gt(0).all()), int(excluded["step02_exclusion_reasons"].str.len().gt(0).sum()), len(excluded), "Every excluded cell must retain at least one exact reason.")
        self.ledger.add("retained_reasons_empty", "result", bool(retained["step02_exclusion_reasons"].eq("").all()), int(retained["step02_exclusion_reasons"].eq("").sum()), len(retained), "Retained cells cannot receive an exclusion reason.")
        retained_ids = pd.Index(retained["cell_id"].astype(str))
        self.ledger.add("retained_cell_order", "result", result.filtered.obs_names.astype(str).equals(retained_ids), "exact" if result.filtered.obs_names.astype(str).equals(retained_ids) else "different", "exact", "The filtered object must contain exactly the retained cells in original order.")
        self.ledger.add("no_low_count_at_5mad", "result", int(result.disposition["low_total_counts_5mad"].sum()) == 0, int(result.disposition["low_total_counts_5mad"].sum()), 0, "The reviewed 5-MAD sensitivity contained no low-count candidates.")
        self.ledger.add("low_gene_count", "result", int(result.disposition["low_n_genes_by_counts_5mad"].sum()) == 67, int(result.disposition["low_n_genes_by_counts_5mad"].sum()), 67, "The reviewed low-gene candidate count must be exact.")
        self.ledger.add("high_mt_count", "result", int(result.disposition["high_pct_counts_mt_5mad"].sum()) == 4_372, int(result.disposition["high_pct_counts_mt_5mad"].sum()), 4_372, "The reviewed high-mitochondrial candidate count must be exact.")
        self.ledger.add("no_overlapping_reasons", "result", int(result.by_sample["multiple_reasons_n"].sum()) == 0, int(result.by_sample["multiple_reasons_n"].sum()), 0, "The approved 5-MAD flag set contained no multi-criterion overlap.")
        expected_sample_removed = {
            "15662-JZ-1": 292, "15662-JZ-2": 577, "15662-JZ-3": 1_120,
            "15662-JZ-4": 112, "15662-JZ-5": 289, "15662-JZ-6": 342,
            "15662-JZ-7": 125, "15662-JZ-8": 122, "15662-JZ-9": 352,
            "15662-JZ-10": 452, "15662-JZ-11": 608, "15662-JZ-12": 48,
        }
        observed_sample_removed = result.by_sample.set_index("group")["cells_removed"].astype(int).to_dict()
        self.ledger.add("sample_removals_exact", "result", observed_sample_removed == expected_sample_removed, observed_sample_removed, expected_sample_removed, "Every technical-sample removal count must equal the reviewed Step 01a result.")
        expected_design_after = {"WT_M": 112_147, "WT_F": 118_674, "HET_F": 100_503, "KO_M": 115_025}
        observed_design_after = result.by_design.set_index("group")["cells_after"].astype(int).to_dict()
        self.ledger.add("design_after_exact", "result", observed_design_after == expected_design_after, observed_design_after, expected_design_after, "Design summaries are reporting-only but must total the reviewed retained cells exactly.")
        self.ledger.add("raw_output_state", "result", sparse.issparse(result.filtered.X) and np.issubdtype(result.filtered.X.dtype, np.integer) and len(result.filtered.layers) == 0 and result.filtered.raw is None and len(result.filtered.obsm) == 0 and len(result.filtered.obsp) == 0, {"sparse": sparse.issparse(result.filtered.X), "dtype": str(result.filtered.X.dtype), "layers": list(result.filtered.layers), "raw": result.filtered.raw is not None, "obsm": list(result.filtered.obsm), "obsp": list(result.filtered.obsp)}, "sparse integer; derived representations empty", "Filtering must preserve raw counts and create no analysis representation.")


class Step02MatrixFingerprint:
    """Hash in-memory and serialized CSR arrays for retained-count identity."""

    @staticmethod
    def memory(matrix) -> str:
        """Hash one in-memory CSR matrix's logical arrays."""

        csr = matrix.tocsr()
        digest = hashlib.sha256(str(tuple(int(value) for value in csr.shape)).encode())
        for name, array in (("data", csr.data), ("indices", csr.indices), ("indptr", csr.indptr)):
            digest.update(name.encode())
            digest.update(str(array.dtype).encode())
            view = memoryview(array).cast("B")
            step = 8 * 1024 * 1024
            for start in range(0, len(view), step):
                digest.update(view[start : start + step])
        return digest.hexdigest()

    @staticmethod
    def h5ad(path: Path) -> str:
        """Hash one serialized H5AD CSR matrix without loading it fully."""

        digest = hashlib.sha256()
        with h5py.File(path, "r") as handle:
            matrix = handle["X"]
            digest.update(str(tuple(int(value) for value in matrix.attrs["shape"])).encode())
            for name in ("data", "indices", "indptr"):
                dataset = matrix[name]
                digest.update(name.encode())
                digest.update(str(dataset.dtype).encode())
                step = max(1, 8 * 1024 * 1024 // dataset.dtype.itemsize)
                for start in range(0, len(dataset), step):
                    digest.update(dataset[start : start + step].tobytes())
        return digest.hexdigest()
