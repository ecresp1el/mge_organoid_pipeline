"""Input-lineage, analysis-topology, and output-integrity validation."""

from __future__ import annotations

import ast
import hashlib
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from .models import ValidationLedger
from .step06_analysis import marker_programs
from .step06_models import Step06Artifacts, Step06Paths, Step06Settings


def sha256(path: Path) -> str:
    """Calculate a streaming SHA-256 checksum."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


class PythonRuntimeBoundaryValidator:
    """Prove that frozen Step 06 scientific modules cannot invoke an R runtime."""

    forbidden_imports = {"rpy2", "subprocess"}

    def validate(self, code_dir: Path, ledger: ValidationLedger) -> None:
        """Reject R files and imports capable of bridging to an external R process."""

        step_files = sorted(code_dir.glob("step06_*.py"))
        non_python = sorted(
            str(path.name)
            for path in code_dir.glob("step06_*")
            if path.is_file() and path.suffix != ".py"
        )
        imported: set[str] = set()
        for path in step_files:
            tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    imported.update(alias.name.split(".")[0] for alias in node.names)
                elif isinstance(node, ast.ImportFrom) and node.module:
                    imported.add(node.module.split(".")[0])
        forbidden_found = sorted(imported.intersection(self.forbidden_imports))
        passed = bool(step_files) and not non_python and not forbidden_found
        ledger.add(
            "python_only_scientific_runtime", "runtime", passed,
            {"step06_python_modules": len(step_files), "non_python": non_python, "forbidden_imports": forbidden_found},
            {"non_python": [], "forbidden_imports": []},
            "Step 06 computation must remain Python-only and must not invoke an R bridge or subprocess.",
        )


class Step06InputValidator:
    """Require the approved Step 02 input and explicit Steps 03-05 bypass."""

    def __init__(self, paths: Step06Paths, settings: Step06Settings, ledger: ValidationLedger):
        """Store immutable contracts and the shared validation ledger."""

        self.paths = paths
        self.settings = settings
        self.ledger = ledger

    def validate(self, adata: ad.AnnData) -> None:
        """Validate lineage, bytes, checksum, matrix state, and metadata."""

        parameter_contract = (
            1 <= self.settings.n_jobs
            and 2 <= self.settings.pca_components < self.settings.n_top_genes
            and 1 <= self.settings.centroid_components <= self.settings.pca_components
            and 2 <= self.settings.n_neighbors < self.settings.expected_cells
            and self.settings.render_max_cells > 0
        )
        self.ledger.add(
            "parameter_contract", "configuration", parameter_contract,
            {
                "n_jobs": self.settings.n_jobs,
                "pca_components": self.settings.pca_components,
                "centroid_components": self.settings.centroid_components,
                "n_top_genes": self.settings.n_top_genes,
                "n_neighbors": self.settings.n_neighbors,
                "render_max_cells": self.settings.render_max_cells,
            },
            "positive and internally consistent",
            "Frozen dimensional, graph, rendering, and worker settings must be valid.",
        )
        status = pd.read_csv(self.paths.input_status, sep="\t", dtype=str, keep_default_na=False)
        approved = len(status) == 1 and status.iloc[0]["run_id"] == self.paths.expected_step02_run_id and status.iloc[0]["status"] == "APPROVED"
        self.ledger.add("step02_status", "input", approved, status.to_dict("records"), "one approved expected run", "Step 06 can consume only the approved Step 02 checkpoint.")
        manifest = pd.read_csv(
            self.paths.input_manifest, sep="\t", dtype=str, keep_default_na=False
        )
        manifested = manifest.loc[
            manifest["relative_path"] == "objects/pcdh19_step02_qc_filtered.h5ad"
        ]
        manifest_matches = (
            len(manifested) == 1
            and int(manifested.iloc[0]["bytes"]) == self.paths.expected_input_bytes
            and manifested.iloc[0]["sha256"] == self.paths.expected_input_sha256
        )
        self.ledger.add(
            "step02_manifest_identity", "input", manifest_matches,
            manifested.to_dict("records"),
            {"bytes": self.paths.expected_input_bytes, "sha256": self.paths.expected_input_sha256},
            "The frozen upstream manifest must identify the exact Step 02 object.",
        )
        workflow = pd.read_csv(self.paths.approval_ledger_snapshot, sep="\t", dtype=str, keep_default_na=False)
        checks = {
            "step02_approved_in_ledger": ((workflow["run_id"] == self.paths.expected_step02_run_id) & (workflow["status"] == "APPROVED")).sum() == 1,
            "step03_rejected_in_ledger": ((workflow["step"] == "03_scdblfinder") & (workflow["status"] == "REJECTED")).sum() == 1,
            "step04_skipped_in_ledger": ((workflow["step"] == "04_ambient_rna_contamination_assessment") & (workflow["status"] == "SKIPPED")).sum() == 1,
            "step05_skipped_in_ledger": ((workflow["step"] == "05_broad_biological_contaminant_assessment") & (workflow["status"] == "SKIPPED")).sum() == 1,
        }
        for name, passed in checks.items():
            self.ledger.add(name, "lineage", passed, passed, True, "The frozen ledger must preserve the explicit Step 02-to-Step 06 lineage.")
        self.ledger.add("bypass_decision_present", "lineage", self.paths.bypass_decision.is_file(), self.paths.bypass_decision, "existing file", "The frozen Steps 04-05 bypass decision is required.")
        observed_bytes = self.paths.input_h5ad.stat().st_size
        self.ledger.add("input_bytes", "input", observed_bytes == self.paths.expected_input_bytes, observed_bytes, self.paths.expected_input_bytes, "Input byte size must match the approved checkpoint.")
        observed_sha = sha256(self.paths.input_h5ad)
        self.ledger.add("input_sha256", "input", observed_sha == self.paths.expected_input_sha256, observed_sha, self.paths.expected_input_sha256, "Input content must be byte-identical.")
        self.ledger.add("input_shape", "object", adata.shape == (self.settings.expected_cells, self.settings.expected_genes), adata.shape, (self.settings.expected_cells, self.settings.expected_genes), "Every approved Step 02 cell and gene must enter Step 06.")
        self.ledger.add("raw_sparse_integer", "object", sparse.issparse(adata.X) and np.issubdtype(adata.X.dtype, np.integer), {"sparse": sparse.issparse(adata.X), "dtype": str(adata.X.dtype)}, "sparse integer", "Step 06 must begin from raw sparse integer counts.")
        required = {self.settings.sample_field, "submitted_sample_name", self.settings.genotype_field, self.settings.sex_field, self.settings.design_field, "total_counts", "n_genes_by_counts", "pct_counts_mt"}
        self.ledger.add("required_metadata", "object", required.issubset(adata.obs.columns), sorted(required - set(adata.obs.columns)) or "none", "none", "All report grouping and QC fields must be present.")
        forbidden = {"capture_id", "scDblFinder_score", "scDblFinder_class"}
        self.ledger.add("no_step03_annotations", "lineage", forbidden.isdisjoint(adata.obs.columns), sorted(forbidden.intersection(adata.obs.columns)) or "none", "none", "Rejected Step 03 annotations must not enter Step 06.")
        symbols = set(adata.var["gene_symbol"].astype(str))
        missing = sorted({gene for genes in marker_programs().values() for gene in genes} - symbols)
        self.ledger.add("marker_program_genes", "biology", not missing, missing or "none", "none", "Every fixed broad-state marker must exist in the delivered panel.")


class Step06OutputValidator:
    """Validate full-cell diagnostics and raw-count preservation."""

    def __init__(self, settings: Step06Settings, ledger: ValidationLedger):
        """Store immutable settings and the shared ledger."""

        self.settings = settings
        self.ledger = ledger

    def validate_artifacts(self, artifacts: Step06Artifacts) -> None:
        """Validate embeddings, graphs, annotations, HVGs, and rendering scope."""

        expected = self.settings.expected_cells
        self.ledger.add("pca_shape", "analysis", artifacts.pca.shape == (expected, self.settings.pca_components), artifacts.pca.shape, (expected, self.settings.pca_components), "PCA must represent all cells.")
        self.ledger.add("umap_shape", "analysis", artifacts.umap.shape == (expected, 2), artifacts.umap.shape, (expected, 2), "The one global UMAP must represent all cells.")
        self.ledger.add("finite_coordinates", "analysis", np.isfinite(artifacts.pca).all() and np.isfinite(artifacts.umap).all(), "finite" if np.isfinite(artifacts.pca).all() and np.isfinite(artifacts.umap).all() else "nonfinite", "finite", "Diagnostic coordinates must be finite.")
        self.ledger.add("graph_shape", "analysis", artifacts.distances.shape == (expected, expected) and artifacts.connectivities.shape == (expected, expected), {"distances": artifacts.distances.shape, "connectivities": artifacts.connectivities.shape}, (expected, expected), "Neighbor graphs must represent all cells.")
        self.ledger.add("hvg_count", "analysis", int(artifacts.hvg_table["highly_variable"].sum()) == self.settings.n_top_genes, int(artifacts.hvg_table["highly_variable"].sum()), self.settings.n_top_genes, "HVG count must match the frozen setting.")
        self.ledger.add("annotation_rows", "analysis", len(artifacts.annotations) == expected, len(artifacts.annotations), expected, "Every cell must receive descriptive diagnostics.")
        self.ledger.add("rendering_is_subset", "rendering", len(artifacts.rendering_cells) <= self.settings.render_max_cells and artifacts.rendering_cells["cell_id"].is_unique, len(artifacts.rendering_cells), f"<= {self.settings.render_max_cells}", "Downsampling is permitted only for rendering.")

    def validate_serialized(self, output: ad.AnnData, original_fingerprint: str, output_fingerprint: str) -> None:
        """Validate the published review checkpoint after round trip."""

        self.ledger.add("output_shape", "output", output.shape == (self.settings.expected_cells, self.settings.expected_genes), output.shape, (self.settings.expected_cells, self.settings.expected_genes), "The output must retain all approved cells and genes.")
        self.ledger.add("raw_matrix_unchanged", "output", output_fingerprint == original_fingerprint, output_fingerprint, original_fingerprint, "No raw count may change.")
        self.ledger.add("output_raw_sparse_integer", "output", sparse.issparse(output.X) and np.issubdtype(output.X.dtype, np.integer), {"sparse": sparse.issparse(output.X), "dtype": str(output.X.dtype)}, "sparse integer", "The output X must remain raw counts.")
        self.ledger.add("output_no_scdblfinder", "output", not any(column.startswith("scDblFinder") or column == "capture_id" for column in output.obs.columns), "none", "none", "Rejected Step 03 fields cannot propagate.")
        self.ledger.add("output_embeddings", "output", {"X_pca", "X_umap"}.issubset(output.obsm.keys()), sorted(output.obsm.keys()), ["X_pca", "X_umap"], "PCA and the one shared UMAP must round-trip.")
        self.ledger.add("output_graphs", "output", {"distances", "connectivities"}.issubset(output.obsp.keys()), sorted(output.obsp.keys()), ["distances", "connectivities"], "The unintegrated graph must round-trip.")
