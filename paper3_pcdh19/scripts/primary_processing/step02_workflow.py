"""Orchestrate approved Step 02 filtering and stop before doublet detection."""

from __future__ import annotations

import anndata as ad
import pandas as pd

from .models import ValidationLedger
from .step02_filtering import ApprovedMADFilter
from .step02_models import Step02Paths, Step02Settings
from .step02_plots import Step02PlotPublisher
from .step02_publishing import AtomicStep02Publisher, Step02ApprovalLedger, Step02ProvenancePublisher
from .step02_validation import Step02InputValidator, Step02MatrixFingerprint, Step02Validator
from .validation import DocumentationAuditor


class Step02Workflow:
    """Coordinate exact input validation, filtering, reports, and publication."""

    def __init__(self, paths: Step02Paths, settings: Step02Settings):
        """Construct focused collaborators around the approved Step 02 policy."""

        self.paths = paths
        self.settings = settings
        self.ledger = ValidationLedger()
        self.input_validator = Step02InputValidator(paths, settings, self.ledger)
        self.validator = Step02Validator(settings, self.ledger)
        self.filter = ApprovedMADFilter(settings)
        self.provenance = Step02ProvenancePublisher()
        self.fingerprint = Step02MatrixFingerprint()

    def run(self):
        """Apply the approved union, publish a raw-count checkpoint, and stop."""

        publisher = AtomicStep02Publisher(self.paths.run_dir)
        try:
            self.input_validator.validate()
            adata = ad.read_h5ad(self.paths.input_h5ad)
            flags = pd.read_csv(self.paths.input_flags, sep="\t", compression="gzip", low_memory=False)
            self.validator.validate_input_object(adata, flags)
            original_var_names = adata.var_names.copy()
            result = self.filter.apply(adata, flags)
            self.validator.validate_result(result, original_var_names)
            documentation = pd.DataFrame(DocumentationAuditor().audit(self.paths.frozen_code_dir))
            self.ledger.add("documented_definitions", "frozen_code", bool((documentation["status"] == "PASS").all()), int((documentation["status"] == "PASS").sum()), len(documentation), "Every frozen Python module, class, and function must have a docstring.")
            plot_manifest = Step02PlotPublisher(self.settings, publisher.figures_dir).publish_all(result.by_sample, result.by_design)
            output_relative = "objects/pcdh19_step02_qc_filtered.h5ad"
            output_h5ad = publisher.stage_dir / output_relative
            memory_fingerprint = self.fingerprint.memory(result.filtered.X)
            result.filtered.write_h5ad(output_h5ad, compression=self.settings.compression)
            serialized_fingerprint = self.fingerprint.h5ad(output_h5ad)
            self.ledger.add("retained_raw_matrix_fingerprint", "object", memory_fingerprint == serialized_fingerprint, serialized_fingerprint, memory_fingerprint, "Serialized output raw counts must match the selected in-memory raw counts exactly.")
            reopened = ad.read_h5ad(output_h5ad, backed="r")
            self.ledger.add("h5ad_roundtrip_shape", "object", reopened.shape == result.filtered.shape, reopened.shape, result.filtered.shape, "The filtered checkpoint must reopen with the exact dimensions.")
            self.ledger.add("h5ad_roundtrip_cells", "object", reopened.obs_names.equals(result.filtered.obs_names), reopened.n_obs, result.filtered.n_obs, "The filtered checkpoint must preserve retained cell identity/order.")
            self.ledger.add("h5ad_roundtrip_genes", "object", reopened.var_names.equals(result.filtered.var_names), reopened.n_vars, result.filtered.n_vars, "The filtered checkpoint must preserve gene identity/order.")
            reopened.file.close()
            result.disposition.to_csv(publisher.tables_dir / "step02_per_cell_disposition.tsv.gz", sep="\t", index=False, compression="gzip")
            result.by_sample.to_csv(publisher.tables_dir / "step02_before_after_by_sample.tsv", sep="\t", index=False)
            result.by_design.to_csv(publisher.tables_dir / "step02_before_after_by_design_group.tsv", sep="\t", index=False)
            result.reason_summary.to_csv(publisher.tables_dir / "step02_exclusion_reason_summary.tsv", sep="\t", index=False)
            self._dictionary().to_csv(publisher.tables_dir / "step02_data_dictionary.tsv", sep="\t", index=False)
            plot_manifest.to_csv(publisher.tables_dir / "plot_manifest.tsv", sep="\t", index=False)
            documentation.to_csv(publisher.tables_dir / "documentation_audit.tsv", sep="\t", index=False)
            self.provenance.software_versions().to_csv(publisher.tables_dir / "software_versions.tsv", sep="\t", index=False)
            checks = self.ledger.to_frame()
            checks.to_csv(publisher.tables_dir / "validation_checks.tsv", sep="\t", index=False)
            self.ledger.require_all_pass()
            run_id = self.paths.run_dir.name
            self.provenance.status_frame(run_id, self.paths.expected_step01_run_id, adata.n_obs, result.filtered.n_obs, adata.n_vars).to_csv(publisher.stage_dir / "STEP_STATUS.tsv", sep="\t", index=False)
            filtering_report = self.provenance.filtering_report(run_id, self.paths.expected_step01_run_id, self.paths.expected_step01a_run_id, result.by_sample, result.by_design, checks)
            (publisher.stage_dir / "STEP02_QC_FILTERING_REPORT.md").write_text(filtering_report, encoding="utf-8")
            self.provenance.output_manifest(publisher.stage_dir).to_csv(publisher.tables_dir / "output_manifest.tsv", sep="\t", index=False)
            publisher.publish()
            Step02ApprovalLedger().update_in_review(self.paths.workflow_root / "APPROVAL_LEDGER.tsv", run_id, self.paths.expected_step01_run_id, output_relative, adata.n_obs, result.filtered.n_obs, adata.n_vars, self._code_version())
            return self.paths.run_dir / output_relative
        except Exception:
            publisher.discard()
            raise

    @staticmethod
    def _dictionary() -> pd.DataFrame:
        """Describe checkpoint and complete disposition semantics."""

        return pd.DataFrame([
            {"asset_or_field": "pcdh19_step02_qc_filtered.h5ad", "definition": "Raw integer counts for retained cells and all genes", "decision_status": "IN_REVIEW checkpoint"},
            {"asset_or_field": "step02_per_cell_disposition.tsv.gz", "definition": "All original cells with selected flags, retained/excluded state, and semicolon-delimited exact reasons", "decision_status": "complete audit"},
            {"asset_or_field": "step02_excluded", "definition": "5-MAD low counts OR low genes OR high mitochondrial percentage", "decision_status": "approved Step 02 filter"},
        ])

    def _code_version(self) -> str:
        """Return the clean repository identity frozen during submission."""

        path = self.paths.run_dir / "provenance" / "repository_state.txt"
        return path.read_text(encoding="utf-8").strip().replace("\n", "; ") if path.exists() else "unavailable"
