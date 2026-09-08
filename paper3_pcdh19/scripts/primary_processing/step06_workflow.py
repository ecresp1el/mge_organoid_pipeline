"""Trace the complete Step 06 execution and its publication boundary.

Entry: ``step06_cli.main`` -> ``Step06Workflow.run``. Order: fingerprint and
load approved Step 02 -> validate lineage -> analyze a working object ->
validate artifacts -> build metrics -> reload raw Step 02 -> attach diagnostic
metadata/coordinates/graphs -> stage tables/report/H5AD -> reopen and validate
raw-count identity -> publish output groups -> record IN_REVIEW.

Progress records bracket major phases; their function labels can name logical
blocks (for example attach_diagnostics), not only standalone Python functions.
On an exception, the workflow records failure and removes remaining staging.
Publication uses individual atomic renames, not a transaction across all
files and the separate approval ledger. No downstream step is launched.
"""

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
from .step06_progress import Step06ProgressTracker
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
    """Execute the additive unintegrated diagnostic workflow and stop.

    Notes
    -----
    Own the ordered lifecycle from approved input to IN_REVIEW publication. Delegate
    analysis, metrics, rendering and validation to focused components; retain
    original QC metadata and reload raw counts before output assembly. No downstream
    workflow is started here.
    """

    def __init__(self, paths: Step06Paths, settings: Step06Settings):
        """Construct the immutable path/settings contract and validation ledger.

        Parameters
        ----------
        paths : Step06Paths
            Exact input identity, frozen lineage files, and new run destinations.
        settings : Step06Settings
            Resolved scientific/rendering controls; see the settings class and tuning
            guide.

        Notes
        -----
        Create the shared validation ledger and a new progress tracker. Tracker construction
        creates provenance_dir and refuses an existing event ledger.
        """

        self.paths = paths
        self.settings = settings
        self.ledger = ValidationLedger()
        self.progress = Step06ProgressTracker(paths.run_dir / "provenance")

    def run(self) -> Path:
        """Run all-cell diagnostics and publish only after every check passes.

        Returns
        -------
        pathlib.Path
            Published objects/pcdh19_step06_unintegrated_diagnostics.h5ad in this run.

        Notes
        -----
        This is the top-level trace of all analysis and file writes. The working
        AnnData is transformed by Step06Analyzer; publication instead starts by
        reloading the raw input, attaching annotations, PCA/UMAP, graphs, HVG flags
        and provenance. X therefore remains sparse integer raw counts; no normalized
        layer is published. Tables/report/checkpoint are staged and validated before
        output groups are renamed. A separate shared ledger update records IN_REVIEW.

        Validation failures raise through require_all_pass(). Any caught exception
        records FAILED, discards remaining staging, and propagates; already-published
        assets cannot be rolled back by discard(). Reading/hashing/reloading the
        large input and output is deliberate integrity work, not repeated analysis.
        No scientific thresholds live here: see analysis.py, metrics.py, and settings.
        """

        self.progress.note(
            "workflow",
            "Step06Workflow.run",
            {
                "input_h5ad": self.paths.input_h5ad,
                "expected_step02_run_id": self.paths.expected_step02_run_id,
                "expected_input_bytes": self.paths.expected_input_bytes,
                "expected_input_sha256": self.paths.expected_input_sha256,
                "run_dir": self.paths.run_dir,
                "settings": vars(self.settings),
            },
            {
                "status": "STARTED",
                "lineage": "approved Step 02 directly; Steps 03-05 bypassed",
                "integration": "none",
            },
        )
        publisher = AtomicStep06Publisher(self.paths.run_dir)
        try:
            with self.progress.track(
                "input.identity",
                "Step02MatrixFingerprint.h5ad",
                {
                    "path": self.paths.input_h5ad,
                    "bytes": self.paths.expected_input_bytes,
                    "purpose": "logical raw-count CSR identity before analysis",
                },
            ) as event:
                original_fingerprint = Step02MatrixFingerprint.h5ad(
                    self.paths.input_h5ad
                )
                event.outputs["raw_matrix_fingerprint"] = original_fingerprint
            with self.progress.track(
                "input.load",
                "anndata.read_h5ad",
                {
                    "path": self.paths.input_h5ad,
                    "backed": False,
                    "expected_shape": [
                        self.settings.expected_cells,
                        self.settings.expected_genes,
                    ],
                },
            ) as event:
                adata = ad.read_h5ad(self.paths.input_h5ad)
                event.outputs.update(
                    {
                        "shape": list(adata.shape),
                        "matrix_type": type(adata.X).__name__,
                        "matrix_dtype": str(adata.X.dtype),
                    }
                )
            with self.progress.track(
                "input.validation",
                "Step06InputValidator.validate",
                {
                    "object_shape": list(adata.shape),
                    "input_status": self.paths.input_status,
                    "input_manifest": self.paths.input_manifest,
                    "approval_ledger_snapshot": self.paths.approval_ledger_snapshot,
                    "bypass_decision": self.paths.bypass_decision,
                    "expected_sha256": self.paths.expected_input_sha256,
                },
            ) as event:
                Step06InputValidator(
                    self.paths, self.settings, self.ledger
                ).validate(adata)
                self.ledger.require_all_pass()
                event.outputs["validation_checks_passed"] = len(
                    self.ledger.to_frame()
                )
            # Keep original raw-count QC metadata; analyzer mutates its working X next.
            original_obs = adata.obs.copy()
            artifacts = Step06Analyzer(self.settings, self.progress).run(adata)
            with self.progress.track(
                "analysis.validation",
                "Step06OutputValidator.validate_artifacts",
                {
                    "cells": len(artifacts.annotations),
                    "pca_shape": list(artifacts.pca.shape),
                    "umap_shape": list(artifacts.umap.shape),
                    "distances_shape": list(artifacts.distances.shape),
                    "hvg_rows": len(artifacts.hvg_table),
                },
            ) as event:
                Step06OutputValidator(self.settings, self.ledger).validate_artifacts(
                    artifacts
                )
                self.ledger.require_all_pass()
                event.outputs["status"] = "PASS"
            with self.progress.track(
                "metrics",
                "Step06MetricBuilder.build",
                {
                    "cells": len(original_obs),
                    "samples": int(
                        original_obs[self.settings.sample_field].nunique()
                    ),
                    "sample_field": self.settings.sample_field,
                    "pca_components_for_centroids": self.settings.centroid_components,
                    "neighbor_graph_scope": "all cells",
                    "pseudobulk_scope": "all cells by sample across HVGs",
                    "inferential_tests": "none",
                },
            ) as event:
                results = Step06MetricBuilder(self.settings).build(
                    original_obs, artifacts
                )
                event.outputs.update(
                    {
                        "table_count": len(results.tables),
                        "table_names": sorted(results.tables),
                        "provisional_outcome_code": results.outcome_code,
                        "provisional_outcome_label": results.outcome_label,
                    }
                )
            # Release transformed counts; the published object is rebuilt from pristine raw input.
            del adata
            gc.collect()

            with self.progress.track(
                "output.reload_raw",
                "anndata.read_h5ad",
                {
                    "path": self.paths.input_h5ad,
                    "purpose": "reload pristine Step 02 raw counts for additive output",
                    "backed": False,
                },
            ) as event:
                # Crucial boundary: output.X comes from approved raw counts, not normalized/scaled X.
                output = ad.read_h5ad(self.paths.input_h5ad)
                event.outputs.update(
                    {
                        "shape": list(output.shape),
                        "matrix_type": type(output.X).__name__,
                        "matrix_dtype": str(output.X.dtype),
                    }
                )
            with self.progress.track(
                "output.assemble",
                "Step06Workflow.attach_diagnostics",
                {
                    "base_object": "fresh approved Step 02 AnnData",
                    "annotation_columns": list(artifacts.annotations.columns),
                    "obsm_keys": ["X_pca", "X_umap"],
                    "obsp_keys": ["connectivities", "distances"],
                    "raw_matrix_mutation": False,
                },
            ) as event:
                for column in artifacts.annotations.columns:
                    output.obs[column] = artifacts.annotations.loc[
                        output.obs_names, column
                    ].to_numpy()
                output.obsm["X_pca"] = artifacts.pca
                output.obsm["X_umap"] = artifacts.umap
                output.obsp["connectivities"] = artifacts.connectivities
                output.obsp["distances"] = artifacts.distances
                hvg = artifacts.hvg_table.set_index("gene_id").reindex(
                    output.var_names
                )
                output.var["highly_variable_step06"] = hvg[
                    "highly_variable"
                ].to_numpy(bool)
                event.outputs.update(
                    {
                        "output_shape": list(output.shape),
                        "obs_columns_added": len(artifacts.annotations.columns),
                        "cells_removed": 0,
                        "genes_removed": 0,
                    }
                )
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

            with self.progress.track(
                "output.tables",
                "pandas.DataFrame.to_csv",
                {
                    "destination": publisher.tables_dir,
                    "metric_tables": sorted(results.tables),
                    "additional_tables": [
                        "highly_variable_genes.tsv.gz",
                        "marker_program_scores_by_cluster.tsv",
                        "rendering_cell_ids.tsv.gz",
                        "pca_variance_ratio.tsv",
                    ],
                },
            ) as event:
                for filename, frame in results.tables.items():
                    compression = "gzip" if filename.endswith(".gz") else None
                    frame.to_csv(
                        publisher.tables_dir / filename,
                        sep="\t",
                        index=False,
                        compression=compression,
                    )
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
                pd.DataFrame(
                    {
                        "component": [
                            f"PC{i + 1}"
                            for i in range(len(artifacts.variance_ratio))
                        ],
                        "variance_ratio": artifacts.variance_ratio,
                    }
                ).to_csv(
                    publisher.tables_dir / "pca_variance_ratio.tsv",
                    sep="\t",
                    index=False,
                )
                event.outputs["files_written"] = len(results.tables) + 4
            with self.progress.track(
                "output.report",
                "Step06ReportPlotter.publish",
                {
                    "all_numerical_cells": len(original_obs),
                    "rendering_cells": len(artifacts.rendering_cells),
                    "rendering_seed": self.settings.random_seed,
                    "plot_dpi": self.settings.plot_dpi,
                    "panels": "A-L",
                    "shared_umap_coordinates": True,
                    "destination": publisher.figures_dir,
                },
            ) as event:
                plot_manifest = Step06ReportPlotter(
                    self.settings, publisher.figures_dir, self.progress
                ).publish(original_obs, artifacts, results)
                plot_manifest.to_csv(
                    publisher.tables_dir / "plot_manifest.tsv", sep="\t", index=False
                )
                event.outputs["assets"] = plot_manifest.to_dict("records")

            with self.progress.track(
                "output.provenance",
                "DocumentationAuditor.audit",
                {"frozen_code_directory": self.paths.frozen_code_dir},
            ) as event:
                documentation = pd.DataFrame(
                    DocumentationAuditor().audit(self.paths.frozen_code_dir)
                )
                documentation.to_csv(
                    publisher.tables_dir / "documentation_audit.tsv",
                    sep="\t",
                    index=False,
                )
                documented = bool((documentation["status"] == "PASS").all())
                self.ledger.add(
                    "documented_definitions",
                    "frozen_code",
                    documented,
                    int((documentation["status"] == "PASS").sum()),
                    len(documentation),
                    "Every frozen Python module, class, and function must have a docstring.",
                )
                PythonRuntimeBoundaryValidator().validate(
                    self.paths.frozen_code_dir, self.ledger
                )
                event.outputs.update(
                    {
                        "definitions_checked": len(documentation),
                        "definitions_documented": int(
                            (documentation["status"] == "PASS").sum()
                        ),
                        "python_only_runtime": True,
                    }
                )
            with self.progress.track(
                "output.provenance",
                "Step06ProvenancePublisher.software_versions",
                {
                    "requirements_file": self.paths.run_dir
                    / "config"
                    / "requirements.txt"
                },
            ) as event:
                versions = Step06ProvenancePublisher().software_versions()
                versions.to_csv(
                    publisher.tables_dir / "software_versions.tsv",
                    sep="\t",
                    index=False,
                )
                event.outputs["versions"] = versions.set_index("component")[
                    "version"
                ].to_dict()

            output_path = (
                publisher.objects_dir
                / "pcdh19_step06_unintegrated_diagnostics.h5ad"
            )
            with self.progress.track(
                "output.serialization",
                "anndata.AnnData.write_h5ad",
                {
                    "destination": output_path,
                    "shape": list(output.shape),
                    "matrix_type": type(output.X).__name__,
                    "matrix_dtype": str(output.X.dtype),
                    "compression": self.settings.compression,
                    "raw_matrix_expected_unchanged": True,
                },
            ) as event:
                output.write_h5ad(output_path, compression=self.settings.compression)
                event.outputs["bytes"] = output_path.stat().st_size
            del output
            gc.collect()
            with self.progress.track(
                "output.validation",
                "Step06OutputValidator.validate_serialized",
                {
                    "path": output_path,
                    "expected_shape": [
                        self.settings.expected_cells,
                        self.settings.expected_genes,
                    ],
                    "expected_raw_matrix_fingerprint": original_fingerprint,
                },
            ) as event:
                output_fingerprint = Step02MatrixFingerprint.h5ad(output_path)
                reopened = ad.read_h5ad(output_path)
                Step06OutputValidator(
                    self.settings, self.ledger
                ).validate_serialized(
                    reopened, original_fingerprint, output_fingerprint
                )
                event.outputs.update(
                    {
                        "observed_shape": list(reopened.shape),
                        "observed_raw_matrix_fingerprint": output_fingerprint,
                        "raw_matrix_identical": output_fingerprint
                        == original_fingerprint,
                    }
                )
            del reopened
            gc.collect()

            checks = self.ledger.to_frame()
            checks.to_csv(publisher.tables_dir / "validation_checks.tsv", sep="\t", index=False)
            # Numerical success is not approval: status remains IN_REVIEW after publication.
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
            with self.progress.track(
                "output.publication",
                "AtomicStep06Publisher.publish",
                {
                    "staging_directory": publisher.stage_dir,
                    "run_directory": self.paths.run_dir,
                    "validation_failures": int(
                        (checks["status"] != "PASS").sum()
                    ),
                    "target_status": "IN_REVIEW",
                },
            ) as event:
                publisher.publish()
                event.outputs["published"] = True
            with self.progress.track(
                "workflow.ledger",
                "Step06ApprovalLedger.update",
                {
                    "ledger": self.paths.workflow_root / "APPROVAL_LEDGER.tsv",
                    "run_id": run_id,
                    "input_run_id": self.paths.expected_step02_run_id,
                    "status": "IN_REVIEW",
                },
            ) as event:
                Step06ApprovalLedger().update(
                    self.paths.workflow_root / "APPROVAL_LEDGER.tsv",
                    run_id,
                    self.paths.expected_step02_run_id,
                    self.settings,
                    results,
                    self._code_version(),
                )
                event.outputs["ledger_updated"] = True
            self.progress.note(
                "workflow",
                "Step06Workflow.run",
                {},
                {
                    "status": "COMPLETE",
                    "output": self.paths.run_dir
                    / "objects"
                    / "pcdh19_step06_unintegrated_diagnostics.h5ad",
                    "review_status": "IN_REVIEW",
                },
            )
            return self.paths.run_dir / "objects" / "pcdh19_step06_unintegrated_diagnostics.h5ad"
        except Exception as error:
            self.progress.note(
                "workflow",
                "Step06Workflow.run",
                {},
                {
                    "status": "FAILED",
                    "error_type": type(error).__name__,
                    "error_message": str(error),
                },
            )
            publisher.discard()
            raise

    def _code_version(self) -> str:
        """Return the clean repository identity frozen at submission.

        Returns
        -------
        str
            Frozen repository_state.txt flattened with semicolons, or unavailable.

        Notes
        -----
        Reads the run snapshot, not the current Git working tree; no files are modified.
        """

        path = self.paths.run_dir / "provenance" / "repository_state.txt"
        return path.read_text(encoding="utf-8").strip().replace("\n", "; ") if path.exists() else "unavailable"
