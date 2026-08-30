"""Orchestrate the read-only Step 01a per-sample MAD sensitivity analysis."""

from __future__ import annotations

import anndata as ad
import pandas as pd

from .models import ValidationLedger
from .step01a_models import Step01aPaths, Step01aSettings
from .step01a_outliers import RobustMADBoundaryCalculator
from .step01a_plots import Step01aPlotPublisher
from .step01a_publishing import AtomicStep01aPublisher, Step01aApprovalLedger, Step01aProvenancePublisher
from .step01a_validation import Step01aInputValidator, Step01aValidator
from .validation import DocumentationAuditor


class Step01aWorkflow:
    """Coordinate exact input validation, robust flags, plots, and publication."""

    def __init__(self, paths: Step01aPaths, settings: Step01aSettings):
        """Construct focused collaborators around immutable Step 01a settings."""

        self.paths = paths
        self.settings = settings
        self.ledger = ValidationLedger()
        self.input_validator = Step01aInputValidator(paths, self.ledger)
        self.validator = Step01aValidator(settings, self.ledger)
        self.calculator = RobustMADBoundaryCalculator(settings)
        self.provenance = Step01aProvenancePublisher()

    def run(self):
        """Calculate, validate, and publish candidate sensitivity without filtering."""

        publisher = AtomicStep01aPublisher(self.paths.run_dir)
        try:
            self.input_validator.validate()
            adata = ad.read_h5ad(self.paths.input_h5ad, backed="r")
            obs = self.validator.validate_input_object(adata)
            n_cells, n_genes = adata.shape
            adata.file.close()
            result = self.calculator.calculate(obs)
            self.validator.validate_results(result)
            documentation = pd.DataFrame(DocumentationAuditor().audit(self.paths.frozen_code_dir))
            self.ledger.add("documented_definitions", "frozen_code", bool((documentation["status"] == "PASS").all()), int((documentation["status"] == "PASS").sum()), len(documentation), "Every frozen Python module, class, and function must have a docstring.")
            plot_manifest = Step01aPlotPublisher(self.settings, publisher.figures_dir).publish_all(obs, result)
            result.boundaries.to_csv(publisher.tables_dir / "mad_candidate_boundaries_by_sample.tsv", sep="\t", index=False)
            result.flag_summary.to_csv(publisher.tables_dir / "mad_candidate_flag_summary_by_sample.tsv", sep="\t", index=False)
            result.overlap_summary.to_csv(publisher.tables_dir / "mad_candidate_overlap_by_sample.tsv", sep="\t", index=False)
            result.per_cell.to_csv(publisher.tables_dir / "per_cell_mad_candidate_flags.tsv.gz", sep="\t", index=False, compression="gzip")
            self._dictionary().to_csv(publisher.tables_dir / "step01a_data_dictionary.tsv", sep="\t", index=False)
            plot_manifest.to_csv(publisher.tables_dir / "plot_manifest.tsv", sep="\t", index=False)
            documentation.to_csv(publisher.tables_dir / "documentation_audit.tsv", sep="\t", index=False)
            self.provenance.software_versions().to_csv(publisher.tables_dir / "software_versions.tsv", sep="\t", index=False)
            checks = self.ledger.to_frame()
            checks.to_csv(publisher.tables_dir / "validation_checks.tsv", sep="\t", index=False)
            self.ledger.require_all_pass()
            run_id = self.paths.run_dir.name
            self.provenance.status_frame(run_id, self.paths.expected_input_run_id, n_cells, n_genes).to_csv(publisher.stage_dir / "STEP_STATUS.tsv", sep="\t", index=False)
            report = self.provenance.report(run_id, self.paths.expected_input_run_id, result.flag_summary, checks, len(plot_manifest))
            (publisher.stage_dir / "STEP01A_MAD_SENSITIVITY_REPORT.md").write_text(report, encoding="utf-8")
            self.provenance.output_manifest(publisher.stage_dir).to_csv(publisher.tables_dir / "output_manifest.tsv", sep="\t", index=False)
            publisher.publish()
            Step01aApprovalLedger().update_in_review(self.paths.workflow_root / "APPROVAL_LEDGER.tsv", run_id, self.paths.expected_input_run_id, n_cells, n_genes, self._code_version())
            return self.paths.run_dir / "tables" / "per_cell_mad_candidate_flags.tsv.gz"
        except Exception:
            publisher.discard()
            raise

    def _dictionary(self) -> pd.DataFrame:
        """Describe calculations and non-decision semantics for published fields."""

        rows = [
            {"field_or_concept": "scaled MAD", "definition": "1.4826 * median(abs(x - median(x)))", "decision_status": "sensitivity only"},
            {"field_or_concept": "low_total_counts", "definition": "total_counts below sample median(log1p counts) minus k scaled MAD, back-transformed", "decision_status": "candidate only"},
            {"field_or_concept": "low_n_genes_by_counts", "definition": "detected genes below sample median(log1p genes) minus k scaled MAD, back-transformed", "decision_status": "candidate only"},
            {"field_or_concept": "high_pct_counts_mt", "definition": "mitochondrial percentage above sample median plus k scaled MAD", "decision_status": "candidate only"},
            {"field_or_concept": "any_candidate", "definition": "union of the three candidate criteria", "decision_status": "candidate only"},
            {"field_or_concept": "all_three_candidate", "definition": "intersection of the three candidate criteria", "decision_status": "candidate only"},
            {"field_or_concept": "upper counts/genes", "definition": "not calculated; high-complexity evaluation reserved for Step 03 scDblFinder", "decision_status": "out of scope"},
            {"field_or_concept": "filtering", "definition": "none; source H5AD opened read-only and no output H5AD created", "decision_status": "prohibited"},
        ]
        return pd.DataFrame(rows)

    def _code_version(self) -> str:
        """Return the clean repository identity frozen during submission."""

        path = self.paths.run_dir / "provenance" / "repository_state.txt"
        return path.read_text(encoding="utf-8").strip().replace("\n", "; ") if path.exists() else "unavailable"
