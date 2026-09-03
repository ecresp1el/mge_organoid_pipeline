"""Orchestrate, validate, and atomically publish Step 06 diagnostics."""

from __future__ import annotations

import gc
from pathlib import Path

import anndata as ad
import pandas as pd

from .models import ValidationLedger
from .step02_validation import Step02MatrixFingerprint
from .step06_analysis import Step06Analyzer
from .step06_metrics import Step06MetricBuilder
from .step06_models import Step06Paths, Step06Settings
from .step06_plots import Step06ReportPlotter
from .step06_publishing import (
    AtomicStep06Publisher,
    Step06ApprovalLedger,
    Step06ProvenancePublisher,
)
from .step06_validation import (
    PythonRuntimeBoundaryValidator,
    Step06InputValidator,
    Step06OutputValidator,
)
from .validation import DocumentationAuditor


class Step06Workflow:
    """Execute the additive unintegrated diagnostic workflow and stop."""

    def __init__(self, paths: Step06Paths, settings: Step06Settings):
        """Construct the immutable path/settings contract and validation ledger."""

        self.paths = paths
        self.settings = settings
        self.ledger = ValidationLedger()

    def run(self) -> Path:
        """Run all-cell diagnostics and publish only after every check passes."""

        publisher = AtomicStep06Publisher(self.paths.run_dir)
        try:
            original_fingerprint = Step02MatrixFingerprint.h5ad(self.paths.input_h5ad)
            adata = ad.read_h5ad(self.paths.input_h5ad)
            Step06InputValidator(self.paths, self.settings, self.ledger).validate(adata)
            self.ledger.require_all_pass()
            original_obs = adata.obs.copy()
            artifacts = Step06Analyzer(self.settings).run(adata)
            Step06OutputValidator(self.settings, self.ledger).validate_artifacts(artifacts)
            results = Step06MetricBuilder(self.settings).build(original_obs, artifacts)
            del adata
            gc.collect()

            output = ad.read_h5ad(self.paths.input_h5ad)
            for column in artifacts.annotations.columns:
                output.obs[column] = artifacts.annotations.loc[output.obs_names, column].to_numpy()
            output.obsm["X_pca"] = artifacts.pca
            output.obsm["X_umap"] = artifacts.umap
            output.obsp["connectivities"] = artifacts.connectivities
            output.obsp["distances"] = artifacts.distances
            hvg = artifacts.hvg_table.set_index("gene_id").reindex(output.var_names)
            output.var["highly_variable_step06"] = hvg["highly_variable"].to_numpy(bool)
            output.uns["neighbors"] = {
                "connectivities_key": "connectivities",
                "distances_key": "distances",
                "params": {
                    "n_neighbors": self.settings.n_neighbors,
                    "n_pcs": self.settings.pca_components,
                    "metric": self.settings.neighbor_metric,
                    "method": "umap",
                    "random_state": self.settings.random_seed,
                    "integrated": False,
                },
            }
            output.uns["step06_diagnostics"] = {
                "input_step02_run_id": self.paths.expected_step02_run_id,
                "input_sha256": self.paths.expected_input_sha256,
                "steps_03_05_bypassed": True,
                "normalization": f"normalize_total target_sum={self.settings.target_sum:g}; log1p",
                "hvg_flavor": self.settings.hvg_flavor,
                "n_top_genes": self.settings.n_top_genes,
                "pca_components": self.settings.pca_components,
                "n_neighbors": self.settings.n_neighbors,
                "neighbor_metric": self.settings.neighbor_metric,
                "umap_min_dist": self.settings.umap_min_dist,
                "umap_spread": self.settings.umap_spread,
                "leiden_resolution": self.settings.leiden_resolution,
                "random_seed": self.settings.random_seed,
                "batch_correction": "none",
                "integration": "none",
                "cells_removed": 0,
                "provisional_outcome_code": results.outcome_code,
                "provisional_outcome_label": results.outcome_label,
            }

            for filename, frame in results.tables.items():
                compression = "gzip" if filename.endswith(".gz") else None
                frame.to_csv(publisher.tables_dir / filename, sep="\t", index=False, compression=compression)
            artifacts.hvg_table.to_csv(
                publisher.tables_dir / "highly_variable_genes.tsv.gz",
                sep="\t",
                index=False,
                compression="gzip",
            )
            artifacts.cluster_program_scores.to_csv(
                publisher.tables_dir / "marker_program_scores_by_cluster.tsv",
                sep="\t",
                index=False,
            )
            artifacts.rendering_cells.to_csv(
                publisher.tables_dir / "rendering_cell_ids.tsv.gz",
                sep="\t",
                index=False,
                compression="gzip",
            )
            pd.DataFrame({
                "component": [f"PC{i + 1}" for i in range(len(artifacts.variance_ratio))],
                "variance_ratio": artifacts.variance_ratio,
            }).to_csv(publisher.tables_dir / "pca_variance_ratio.tsv", sep="\t", index=False)
            plot_manifest = Step06ReportPlotter(
                self.settings, publisher.figures_dir
            ).publish(original_obs, artifacts, results)
            plot_manifest.to_csv(publisher.tables_dir / "plot_manifest.tsv", sep="\t", index=False)

            documentation = pd.DataFrame(DocumentationAuditor().audit(self.paths.frozen_code_dir))
            documentation.to_csv(publisher.tables_dir / "documentation_audit.tsv", sep="\t", index=False)
            documented = bool((documentation["status"] == "PASS").all())
            self.ledger.add("documented_definitions", "frozen_code", documented, int((documentation["status"] == "PASS").sum()), len(documentation), "Every frozen Python module, class, and function must have a docstring.")
            PythonRuntimeBoundaryValidator().validate(self.paths.frozen_code_dir, self.ledger)
            Step06ProvenancePublisher().software_versions().to_csv(
                publisher.tables_dir / "software_versions.tsv", sep="\t", index=False
            )

            output_path = publisher.objects_dir / "pcdh19_step06_unintegrated_diagnostics.h5ad"
            output.write_h5ad(output_path, compression=self.settings.compression)
            del output
            gc.collect()
            output_fingerprint = Step02MatrixFingerprint.h5ad(output_path)
            reopened = ad.read_h5ad(output_path)
            Step06OutputValidator(self.settings, self.ledger).validate_serialized(
                reopened, original_fingerprint, output_fingerprint
            )
            del reopened
            gc.collect()

            checks = self.ledger.to_frame()
            checks.to_csv(publisher.tables_dir / "validation_checks.tsv", sep="\t", index=False)
            self.ledger.require_all_pass()
            run_id = self.paths.run_dir.name
            Step06ProvenancePublisher().status_frame(
                run_id, self.paths.expected_step02_run_id,
                self.settings.expected_cells, self.settings.expected_genes,
            ).to_csv(publisher.stage_dir / "STEP_STATUS.tsv", sep="\t", index=False)
            report = Step06ProvenancePublisher().report(
                run_id, self.paths.expected_step02_run_id, results, checks, self.settings
            )
            (publisher.stage_dir / "STEP06_DIAGNOSTIC_REPORT.md").write_text(report, encoding="utf-8")
            Step06ProvenancePublisher().output_manifest(publisher.stage_dir).to_csv(
                publisher.tables_dir / "output_manifest.tsv", sep="\t", index=False
            )
            publisher.publish()
            Step06ApprovalLedger().update(
                self.paths.workflow_root / "APPROVAL_LEDGER.tsv",
                run_id,
                self.paths.expected_step02_run_id,
                self.settings,
                results,
                self._code_version(),
            )
            return self.paths.run_dir / "objects" / "pcdh19_step06_unintegrated_diagnostics.h5ad"
        except Exception:
            publisher.discard()
            raise

    def _code_version(self) -> str:
        """Return the clean repository identity frozen at submission."""

        path = self.paths.run_dir / "provenance" / "repository_state.txt"
        return path.read_text(encoding="utf-8").strip().replace("\n", "; ") if path.exists() else "unavailable"
