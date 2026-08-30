"""Scientific-boundary, matrix-integrity, and documentation validations."""

from __future__ import annotations

import ast
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import sparse

from .models import SampleStructure, Step00Settings, ValidationLedger


class DocumentationAuditor:
    """Require module, class, and function docstrings in frozen Python code."""

    def audit(self, source_dir: Path) -> list[dict[str, str]]:
        """Return one documentation-status record per Python definition."""

        records: list[dict[str, str]] = []
        for path in sorted(source_dir.glob("*.py")):
            tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
            records.append(self._record(path, "module", path.stem, ast.get_docstring(tree)))
            for node in ast.walk(tree):
                if isinstance(node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
                    kind = "class" if isinstance(node, ast.ClassDef) else "function"
                    records.append(self._record(path, kind, node.name, ast.get_docstring(node)))
        return records

    def _record(self, path: Path, kind: str, name: str, docstring: str | None) -> dict[str, str]:
        """Construct one stable row for the documentation audit table."""

        return {
            "file": path.name,
            "definition_type": kind,
            "definition_name": name,
            "status": "PASS" if docstring else "FAIL",
            "docstring_summary": docstring.splitlines()[0] if docstring else "",
        }


class Step00Validator:
    """Apply all structural assertions allowed before QC-metric calculation."""

    def __init__(self, settings: Step00Settings, ledger: ValidationLedger):
        """Store the immutable Step 00 contract and shared validation ledger."""

        self.settings = settings
        self.ledger = ledger

    def validate_structure(
        self,
        structure: SampleStructure,
        metadata_vendor_cells: int,
        reference_signature: str,
    ) -> None:
        """Validate one sample's files, dimensions, feature order, and cell totals."""

        scope = structure.sample_id
        unique_types = sorted(set(structure.feature_types))
        unique_genomes = sorted(set(structure.genomes))
        self.ledger.add("feature_count", scope, structure.n_features == self.settings.expected_genes, structure.n_features, self.settings.expected_genes, "Filtered matrix feature count must match the registered delivery.")
        self.ledger.add("feature_type", scope, unique_types == [self.settings.feature_type], unique_types, [self.settings.feature_type], "Only Gene Expression features may enter Step 00.")
        self.ledger.add("genome", scope, unique_genomes == [self.settings.expected_genome], unique_genomes, [self.settings.expected_genome], "All features must use the registered GRCm39 genome.")
        self.ledger.add("feature_order", scope, structure.feature_signature == reference_signature, structure.feature_signature, reference_signature, "Every sample must have the identical ordered feature definition.")
        self.ledger.add("gene_ids_unique", scope, len(set(structure.gene_ids)) == len(structure.gene_ids), len(set(structure.gene_ids)), len(structure.gene_ids), "Ensembl gene IDs must be unique before use as var_names.")
        self.ledger.add("barcodes_unique", scope, len(set(structure.barcodes)) == len(structure.barcodes), len(set(structure.barcodes)), len(structure.barcodes), "Cell Ranger barcodes must be unique within sample.")
        self.ledger.add("metrics_cell_count", scope, structure.vendor_metrics_cells == structure.n_cells, structure.vendor_metrics_cells, structure.n_cells, "Cell Ranger metrics and filtered H5 cell counts must agree.")
        self.ledger.add("manifest_cell_count", scope, metadata_vendor_cells == structure.n_cells, metadata_vendor_cells, structure.n_cells, "Registered technical manifest and filtered H5 cell counts must agree.")
        self.ledger.add("raw_matrix_registered", scope, structure.raw_n_barcodes >= structure.n_cells, structure.raw_n_barcodes, f">={structure.n_cells}", "Raw droplets are registered for Step 04 but do not enter this object.")
        self.ledger.add("raw_gene_expression_feature_superset", scope, structure.raw_n_features >= structure.n_features, structure.raw_n_features, f">={structure.n_features}", "The registered raw matrix may contain additional genes but cannot omit the filtered transcriptome features.")
        self.ledger.add("raw_contains_filtered_features", scope, structure.raw_contains_filtered_features, structure.raw_contains_filtered_features, True, "Every canonical filtered gene must be present in the registered raw matrix.")
        self.ledger.add("raw_preserves_filtered_feature_order", scope, structure.raw_preserves_filtered_feature_order, structure.raw_preserves_filtered_feature_order, True, "Canonical genes selected from the raw matrix must retain the filtered feature definition and order.")

    def validate_loaded_sample(self, adata: ad.AnnData, structure: SampleStructure) -> None:
        """Validate sparse raw integer counts after one sample is loaded."""

        scope = structure.sample_id
        self.ledger.add("loaded_shape", scope, adata.shape == (structure.n_cells, self.settings.expected_genes), adata.shape, (structure.n_cells, self.settings.expected_genes), "Loaded AnnData must retain every called cell and gene.")
        self.ledger.add("sparse_counts", scope, sparse.issparse(adata.X), type(adata.X).__name__, "SciPy sparse matrix", "Canonical counts must remain sparse.")
        self.ledger.add("integer_dtype", scope, np.issubdtype(adata.X.dtype, np.integer), adata.X.dtype, "integer", "Canonical Cell Ranger counts must remain integer-valued.")
        minimum = int(adata.X.data.min()) if adata.X.nnz else 0
        self.ledger.add("nonnegative_counts", scope, minimum >= 0, minimum, ">=0", "Raw UMI counts cannot be negative.")
        self.ledger.add("obs_names_unique", scope, adata.obs_names.is_unique, adata.n_obs, adata.n_obs, "Prefixed cell identifiers must be unique.")

    def validate_combined(self, adata: ad.AnnData) -> None:
        """Validate the complete canonical object before and after serialization."""

        self.ledger.add("combined_shape", "combined", adata.shape == (self.settings.expected_cells, self.settings.expected_genes), adata.shape, (self.settings.expected_cells, self.settings.expected_genes), "Combined object must contain every registered called cell and gene.")
        self.ledger.add("combined_obs_unique", "combined", adata.obs_names.is_unique, adata.n_obs, adata.n_obs, "Combined cell identifiers must be globally unique.")
        self.ledger.add("combined_var_unique", "combined", adata.var_names.is_unique, adata.n_vars, adata.n_vars, "Combined Ensembl gene identifiers must be unique.")
        self.ledger.add("combined_sparse", "combined", sparse.issparse(adata.X), type(adata.X).__name__, "SciPy sparse matrix", "Combined raw counts must remain sparse.")
        self.ledger.add("combined_integer", "combined", np.issubdtype(adata.X.dtype, np.integer), adata.X.dtype, "integer", "Combined raw counts must remain integer-valued.")
        forbidden = sorted(set(adata.obsm.keys()) | set(adata.obsp.keys()))
        self.ledger.add("no_reductions_or_graphs", "combined", not forbidden, forbidden or "none", "none", "Step 00 cannot contain embeddings, neighbors, or graph results.")
