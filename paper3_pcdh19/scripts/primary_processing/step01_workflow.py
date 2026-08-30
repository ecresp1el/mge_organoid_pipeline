"""Orchestration for Scanpy QC-metric calculation without filtering."""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import pandas as pd

from .models import ValidationLedger
from .step01_metrics import QCFeatureSetAnnotator, QCMetadataDictionary, QCSummaryBuilder, ScanpyQCMetricsCalculator
from .step01_models import Step01Paths, Step01Settings
from .step01_plots import QCPlotPublisher
from .step01_publishing import AtomicStep01Publisher, Step01ApprovalLedger, Step01ProvenancePublisher
from .step01_validation import ApprovedStep00Validator, H5ADMatrixFingerprint, Step01Validator
from .validation import DocumentationAuditor


class Step01Workflow:
    """Coordinate approved input loading, Scanpy metrics, plots, and publication."""

    def __init__(self, paths: Step01Paths, settings: Step01Settings):
        """Construct focused collaborators around immutable Step 01 settings."""

        self.paths = paths
        self.settings = settings
        self.ledger = ValidationLedger()
        self.input_validator = ApprovedStep00Validator(paths, self.ledger)
        self.validator = Step01Validator(settings, self.ledger)
        self.feature_annotator = QCFeatureSetAnnotator(settings)
        self.calculator = ScanpyQCMetricsCalculator(settings)
        self.summaries = QCSummaryBuilder()
        self.provenance = Step01ProvenancePublisher()
        self.fingerprint = H5ADMatrixFingerprint()

    def run(self) -> Path:
        """Calculate metrics, validate invariance, publish, and leave IN_REVIEW."""

        publisher = AtomicStep01Publisher(self.paths.run_dir)
        try:
            self.input_validator.validate()
            input_matrix_fingerprint = self.fingerprint.calculate(self.paths.input_h5ad)
            adata = ad.read_h5ad(self.paths.input_h5ad)
            self.validator.validate_input_object(adata)
            original_obs_names = adata.obs_names.copy()
            original_var_names = adata.var_names.copy()
            before_obs = tuple(adata.obs.columns)
            before_var = tuple(adata.var.columns)
            feature_sets = self.feature_annotator.annotate(adata)
            qc_vars = self.calculator.calculate(adata, feature_sets)
            self.validator.validate_metrics(adata, feature_sets)
            self.ledger.add("cell_order_preserved", "object", adata.obs_names.equals(original_obs_names), adata.n_obs, len(original_obs_names), "Step 01 must preserve exact cell identity and order.")
            self.ledger.add("gene_order_preserved", "object", adata.var_names.equals(original_var_names), adata.n_vars, len(original_var_names), "Step 01 must preserve exact gene identity and order.")
            primary = dict(adata.uns.get("primary_processing", {}))
            primary["latest_step"] = "01_qc_metrics"
            primary["step01_qc_metrics"] = {
                "calculator": "scanpy.pp.calculate_qc_metrics",
                "qc_vars": list(qc_vars),
                "percent_top": None if self.settings.percent_top is None else list(self.settings.percent_top),
                "mitochondrial_prefix": self.settings.mitochondrial_prefix,
                "ribosomal_prefixes": list(self.settings.ribosomal_prefixes),
                "ribosomal_metrics_created": bool(feature_sets.ribosomal_genes),
                "filtering": "none",
                "thresholds": "none",
            }
            adata.uns["primary_processing"] = primary
            pooled = self.summaries.pooled(adata.obs)
            by_sample = self.summaries.by_sample(adata.obs)
            by_design = self.summaries.by_design_group(adata.obs)
            metadata_dictionary = QCMetadataDictionary().build(before_obs, before_var, adata)
            documentation = pd.DataFrame(DocumentationAuditor().audit(self.paths.frozen_code_dir))
            self.ledger.add("documented_definitions", "frozen_code", bool((documentation["status"] == "PASS").all()), int((documentation["status"] == "PASS").sum()), len(documentation), "Every frozen Python module, class, and function must have a docstring.")
            plot_manifest = QCPlotPublisher(self.settings, publisher.figures_dir).publish_all(adata.obs)
            output_checkpoint = "objects/pcdh19_step01_qc_metrics.h5ad"
            output_h5ad = publisher.stage_dir / output_checkpoint
            adata.write_h5ad(output_h5ad, compression=self.settings.compression)
            output_matrix_fingerprint = self.fingerprint.calculate(output_h5ad)
            self.ledger.add("raw_matrix_fingerprint_preserved", "object", output_matrix_fingerprint == input_matrix_fingerprint, output_matrix_fingerprint, input_matrix_fingerprint, "Serialized Step 01 X/data, indices, and indptr must exactly match approved Step 00 counts.")
            reopened = ad.read_h5ad(output_h5ad, backed="r")
            self.ledger.add("h5ad_roundtrip_shape", "object", reopened.shape == adata.shape, reopened.shape, adata.shape, "Step 01 checkpoint must reopen with unchanged dimensions.")
            self.ledger.add("h5ad_roundtrip_obs", "object", reopened.obs_names.equals(adata.obs_names), reopened.n_obs, adata.n_obs, "Step 01 checkpoint must preserve ordered cells.")
            self.ledger.add("h5ad_roundtrip_var", "object", reopened.var_names.equals(adata.var_names), reopened.n_vars, adata.n_vars, "Step 01 checkpoint must preserve ordered genes.")
            reopened.file.close()
            feature_summary = pd.DataFrame([
                {"feature_set": "mitochondrial", "definition": f"gene_symbol starts with {self.settings.mitochondrial_prefix}", "n_genes": len(feature_sets.mitochondrial_genes), "metrics_created": True, "interpretation": "Fraction reflects the targeted 13 mitochondrial protein-coding genes."},
                {"feature_set": "cytosolic_ribosomal", "definition": f"gene_symbol starts with {self.settings.ribosomal_prefixes}", "n_genes": len(feature_sets.ribosomal_genes), "metrics_created": bool(feature_sets.ribosomal_genes), "interpretation": "No Rpl/Rps genes are present; ribosomal fraction is unavailable."},
            ])
            feature_genes = pd.DataFrame(
                [{"feature_set": "mitochondrial", "gene_symbol": gene} for gene in feature_sets.mitochondrial_genes]
                + [{"feature_set": "cytosolic_ribosomal", "gene_symbol": gene} for gene in feature_sets.ribosomal_genes]
            )
            pooled.to_csv(publisher.tables_dir / "qc_summary_pooled.tsv", sep="\t", index=False)
            by_sample.to_csv(publisher.tables_dir / "qc_summary_by_sample.tsv", sep="\t", index=False)
            by_design.to_csv(publisher.tables_dir / "qc_summary_by_design_group.tsv", sep="\t", index=False)
            metadata_dictionary.to_csv(publisher.tables_dir / "qc_metric_dictionary.tsv", sep="\t", index=False)
            feature_summary.to_csv(publisher.tables_dir / "qc_feature_set_summary.tsv", sep="\t", index=False)
            feature_genes.to_csv(publisher.tables_dir / "qc_feature_genes.tsv", sep="\t", index=False)
            plot_manifest.to_csv(publisher.tables_dir / "plot_manifest.tsv", sep="\t", index=False)
            documentation.to_csv(publisher.tables_dir / "documentation_audit.tsv", sep="\t", index=False)
            self.provenance.software_versions().to_csv(publisher.tables_dir / "software_versions.tsv", sep="\t", index=False)
            object_summary = pd.DataFrame([{
                "n_cells": adata.n_obs,
                "n_genes": adata.n_vars,
                "nnz": int(adata.X.nnz),
                "matrix_class": type(adata.X).__name__,
                "matrix_dtype": str(adata.X.dtype),
                "X_state": "unchanged sparse unnormalized integer Cell Ranger counts",
                "new_obs_fields": json.dumps([field for field in adata.obs.columns if field not in before_obs]),
                "new_var_fields": json.dumps([field for field in adata.var.columns if field not in before_var]),
                "layers": json.dumps(list(adata.layers.keys())),
                "raw": adata.raw is not None,
                "obsm": json.dumps(list(adata.obsm.keys())),
                "obsp": json.dumps(list(adata.obsp.keys())),
                "filtering": "none",
                "thresholds": "none",
            }])
            object_summary.to_csv(publisher.tables_dir / "step01_object_summary.tsv", sep="\t", index=False)
            checks = self.ledger.to_frame()
            checks.to_csv(publisher.tables_dir / "validation_checks.tsv", sep="\t", index=False)
            self.ledger.require_all_pass()
            run_id = self.paths.run_dir.name
            status = self.provenance.status_frame(run_id, self.paths.expected_input_run_id, output_checkpoint, adata.n_obs, adata.n_vars)
            status.to_csv(publisher.stage_dir / "STEP_STATUS.tsv", sep="\t", index=False)
            report = self.provenance.report(run_id, self.paths.expected_input_run_id, output_checkpoint, adata.n_vars, pooled, by_sample, checks, len(feature_sets.mitochondrial_genes), len(feature_sets.ribosomal_genes), len(plot_manifest))
            (publisher.stage_dir / "STEP01_QC_METRICS_REPORT.md").write_text(report, encoding="utf-8")
            manifest = self.provenance.output_manifest(publisher.stage_dir)
            manifest.to_csv(publisher.tables_dir / "output_manifest.tsv", sep="\t", index=False)
            publisher.publish()
            Step01ApprovalLedger().update_in_review(self.paths.workflow_root / "APPROVAL_LEDGER.tsv", run_id, self.paths.expected_input_run_id, output_checkpoint, adata.n_obs, adata.n_vars, self._code_version())
            return self.paths.run_dir / output_checkpoint
        except Exception:
            publisher.discard()
            raise

    def _code_version(self) -> str:
        """Return the clean repository identity frozen during submission."""

        path = self.paths.run_dir / "provenance" / "repository_state.txt"
        return path.read_text(encoding="utf-8").strip().replace("\n", "; ") if path.exists() else "unavailable"
