"""Finalize native-R Step 03 outputs into a validated AnnData review package."""

from __future__ import annotations

import re
import shutil
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from .models import ValidationLedger
from .step02_validation import Step02MatrixFingerprint
from .step03_models import Step03Paths, Step03Settings
from .step03_plots import Step03PlotPublisher
from .step03_publishing import (
    AtomicStep03Publisher,
    Step03ApprovalLedger,
    Step03ProvenancePublisher,
    Step03ReportBuilder,
)
from .step03_validation import Step03ResultValidator
from .validation import DocumentationAuditor


class RDocumentationAuditor:
    """Require nearby roxygen documentation for each frozen R function definition."""

    PATTERN = re.compile(r"^\s*(?:[A-Za-z0-9_.]+\s*(?:<-|=)\s*)?function\s*\(|R6Class\s*\(")

    def audit(self, path: Path) -> pd.DataFrame:
        """Return one documentation record for every R function or R6 class line."""

        lines = path.read_text(encoding="utf-8").splitlines()
        records = []
        for index, line in enumerate(lines):
            if not self.PATTERN.search(line):
                continue
            context = lines[max(0, index - 12) : index]
            documented = any(candidate.lstrip().startswith("#'") for candidate in context)
            records.append({
                "file": path.name,
                "definition_type": "R function or R6 class",
                "definition_name": line.strip()[:100],
                "status": "PASS" if documented else "FAIL",
                "docstring_summary": next((candidate.lstrip()[2:].strip() for candidate in reversed(context) if candidate.lstrip().startswith("#'")), ""),
            })
        return pd.DataFrame(records)


class Step03FinalizationWorkflow:
    """Attach detector outputs, validate count identity, publish diagnostics, and stop."""

    def __init__(self, paths: Step03Paths, settings: Step03Settings):
        """Construct run-scoped settings, validators, and a shared ledger."""

        self.paths = paths
        self.settings = settings
        self.ledger = ValidationLedger()

    def run(self) -> Path:
        """Publish the non-filtering Step 03 checkpoint in IN_REVIEW state."""

        publisher = AtomicStep03Publisher(self.paths.run_dir)
        try:
            prepare_checks = pd.read_csv(self.paths.intermediate_dir / "prepare_validation_checks.tsv", sep="\t", dtype=str, keep_default_na=False)
            if (prepare_checks["status"] != "PASS").any():
                raise ValueError("Step 03 preparation contained failed validation checks")
            adata = ad.read_h5ad(self.paths.input_h5ad)
            original_fingerprint = Step02MatrixFingerprint.h5ad(self.paths.input_h5ad)
            results = pd.read_csv(self.paths.r_output_dir / "scdblfinder_per_cell_results.tsv.gz", sep="\t", compression="gzip", low_memory=False)
            pca_ids, pca = self._read_pca(self.paths.r_output_dir / "scdblfinder_internal_pca.tsv.gz")
            Step03ResultValidator(self.settings, self.ledger).validate(adata, results, pca_ids, pca, original_fingerprint)
            frame = self._attach_results(adata, results, pca)

            documentation = self._documentation_audit()
            documented = bool((documentation["status"] == "PASS").all())
            self.ledger.add("documented_definitions", "frozen_code", documented, int((documentation["status"] == "PASS").sum()), len(documentation), "Every frozen Python module/class/function and native-R function/R6 class must be documented.")

            plot_manifest = Step03PlotPublisher(self.settings, publisher.figures_dir).publish_all(frame, pca)
            report_builder = Step03ReportBuilder()
            summaries = report_builder.summaries(frame, self.settings.sample_field, self.settings.design_field)
            self._write_tables(publisher, frame, summaries, plot_manifest, documentation)
            self._copy_native_r_audit_objects(publisher)

            output_relative = "objects/pcdh19_step03_scdblfinder.h5ad"
            output_h5ad = publisher.stage_dir / output_relative
            adata.write_h5ad(output_h5ad, compression=self.settings.compression)
            output_fingerprint = Step02MatrixFingerprint.h5ad(output_h5ad)
            self.ledger.add("serialized_raw_matrix_unchanged", "object", output_fingerprint == original_fingerprint, output_fingerprint, original_fingerprint, "The published checkpoint must preserve every raw sparse count exactly.")
            reopened = ad.read_h5ad(output_h5ad, backed="r")
            self.ledger.add("h5ad_roundtrip_shape", "object", reopened.shape == adata.shape, reopened.shape, adata.shape, "The checkpoint must reopen with every approved cell and gene.")
            self.ledger.add("h5ad_roundtrip_pca", "object", "X_scdblfinder_pca" in reopened.obsm and reopened.obsm["X_scdblfinder_pca"].shape == pca.shape, list(reopened.obsm.keys()), f"X_scdblfinder_pca {pca.shape}", "The exact real-cell internal PCA must round-trip.")
            reopened.file.close()

            checks = pd.concat([prepare_checks, self.ledger.to_frame()], ignore_index=True)
            checks.to_csv(publisher.tables_dir / "validation_checks.tsv", sep="\t", index=False)
            if (checks["status"] != "PASS").any():
                failed = checks.loc[checks["status"] != "PASS", ["scope", "check_id"]].astype(str).agg(": ".join, axis=1).tolist()
                raise ValueError(f"Step 03 validation failed: {', '.join(failed)}")
            run_id = self.paths.run_dir.name
            Step03ProvenancePublisher().status_frame(run_id, self.paths.expected_step02_run_id, adata.n_obs, adata.n_vars).to_csv(publisher.stage_dir / "STEP_STATUS.tsv", sep="\t", index=False)
            report = report_builder.report(run_id, self.paths.expected_step02_run_id, summaries, checks, pca.shape[1])
            (publisher.stage_dir / "STEP03_SCDBLFINDER_REPORT.md").write_text(report, encoding="utf-8")
            Step03ProvenancePublisher().output_manifest(publisher.stage_dir).to_csv(publisher.tables_dir / "output_manifest.tsv", sep="\t", index=False)
            publisher.publish()
            overall = summaries["overall"].iloc[0]
            Step03ApprovalLedger().update_in_review(
                self.paths.workflow_root / "APPROVAL_LEDGER.tsv",
                run_id,
                self.paths.expected_step02_run_id,
                int(overall["primary_called_doublets"]),
                float(overall["primary_called_doublet_pct"]),
                self._code_version(),
            )
            return self.paths.run_dir / output_relative
        except Exception:
            publisher.discard()
            raise

    @staticmethod
    def _read_pca(path: Path) -> tuple[pd.Index, np.ndarray]:
        """Read the aligned real-cell PCA written by native R."""

        table = pd.read_csv(path, sep="\t", compression="gzip")
        identifiers = pd.Index(table.pop("cell_id").astype(str))
        return identifiers, table.to_numpy(dtype=float)

    def _attach_results(self, adata: ad.AnnData, results: pd.DataFrame, pca: np.ndarray) -> pd.DataFrame:
        """Attach review-only fields and describe the exact saved object state."""

        adata.obs["capture_id"] = pd.Categorical([self.settings.capture_id] * adata.n_obs)
        adata.obs["scDblFinder_score"] = results["primary_score"].to_numpy(float)
        adata.obs["scDblFinder_class"] = pd.Categorical(results["primary_class"], categories=["singlet", "doublet"])
        adata.obs["scDblFinder_cluster"] = pd.Categorical(results["primary_cluster"].astype(str))
        adata.obs["scDblFinder_replicate_score"] = results["replicate_score"].to_numpy(float)
        adata.obs["scDblFinder_replicate_class"] = pd.Categorical(results["replicate_class"], categories=["singlet", "doublet"])
        adata.obs["scDblFinder_call_reproduced"] = results["primary_class"].astype(str).to_numpy() == results["replicate_class"].astype(str).to_numpy()
        adata.obsm["X_scdblfinder_pca"] = np.asarray(pca, dtype=np.float32)
        adata.uns["step03_scdblfinder"] = {
            "capture_definition": "one independently processed capture",
            "capture_id": self.settings.capture_id,
            "samples_argument": "capture_id",
            "clusters": True,
            "dbr_sd": 1.0,
            "dbr": "not supplied",
            "other_model_parameters": "scDblFinder package defaults",
            "primary_seed": self.settings.primary_seed,
            "reproducibility_seed": self.settings.reproducibility_seed,
            "cells_removed": 0,
            "pca_state": "internal primary-run diagnostic PCA for real cells; not an integrated UMAP",
        }
        frame = pd.DataFrame({
            "cell_id": adata.obs_names.astype(str),
            self.settings.sample_field: adata.obs[self.settings.sample_field].astype(str).to_numpy(),
            self.settings.design_field: adata.obs[self.settings.design_field].astype(str).to_numpy(),
            "primary_score": results["primary_score"].to_numpy(float),
            "primary_class": results["primary_class"].astype(str).to_numpy(),
            "primary_cluster": results["primary_cluster"].astype(str).to_numpy(),
            "replicate_score": results["replicate_score"].to_numpy(float),
            "replicate_class": results["replicate_class"].astype(str).to_numpy(),
        })
        return frame

    def _documentation_audit(self) -> pd.DataFrame:
        """Audit every frozen Python and native-R definition for documentation."""

        python = pd.DataFrame(DocumentationAuditor().audit(self.paths.frozen_code_dir))
        r_path = self.paths.frozen_code_dir / "step03_scdblfinder.R"
        r_records = RDocumentationAuditor().audit(r_path)
        return pd.concat([python, r_records], ignore_index=True)

    def _write_tables(self, publisher: AtomicStep03Publisher, frame: pd.DataFrame, summaries: dict[str, pd.DataFrame], plot_manifest: pd.DataFrame, documentation: pd.DataFrame) -> None:
        """Write complete calls, diagnostics, dictionaries, and software tables."""

        shutil.copy2(self.paths.r_output_dir / "scdblfinder_per_cell_results.tsv.gz", publisher.tables_dir / "scdblfinder_per_cell_details.tsv.gz")
        for key, filename in (
            ("overall", "scdblfinder_overall_summary.tsv"),
            ("sample", "scdblfinder_by_technical_sample.tsv"),
            ("design", "scdblfinder_by_design_group.tsv"),
            ("cluster", "scdblfinder_by_generated_cluster.tsv"),
            ("score_quantiles", "scdblfinder_score_quantiles_by_sample.tsv"),
            ("confusion", "scdblfinder_reproducibility_confusion.tsv"),
            ("reproducibility_sample", "scdblfinder_reproducibility_by_technical_sample.tsv"),
            ("reproducibility_design", "scdblfinder_reproducibility_by_design_group.tsv"),
        ):
            summaries[key].to_csv(publisher.tables_dir / filename, sep="\t", index=False)
        frame[["cell_id", "primary_score", "primary_class", "primary_cluster", "replicate_score", "replicate_class"]].to_csv(publisher.tables_dir / "scdblfinder_canonical_calls.tsv.gz", sep="\t", index=False, compression="gzip")
        Step03ReportBuilder.data_dictionary().to_csv(publisher.tables_dir / "step03_data_dictionary.tsv", sep="\t", index=False)
        plot_manifest.to_csv(publisher.tables_dir / "plot_manifest.tsv", sep="\t", index=False)
        documentation.to_csv(publisher.tables_dir / "documentation_audit.tsv", sep="\t", index=False)
        Step03ProvenancePublisher().python_versions().to_csv(publisher.tables_dir / "python_software_versions.tsv", sep="\t", index=False)
        for filename in ("r_software_versions.tsv", "r_session_info.txt", "scdblfinder_method_contract.tsv"):
            shutil.copy2(self.paths.r_output_dir / filename, publisher.tables_dir / filename)

    def _copy_native_r_audit_objects(self, publisher: AtomicStep03Publisher) -> None:
        """Retain native-R model statistics and internal PCA evidence."""

        shutil.copy2(self.paths.r_output_dir / "scdblfinder_primary_stats.rds", publisher.objects_dir / "scdblfinder_primary_stats.rds")
        shutil.copy2(self.paths.r_output_dir / "scdblfinder_internal_pca.rds", publisher.objects_dir / "scdblfinder_internal_pca.rds")

    def _code_version(self) -> str:
        """Return the clean repository identity recorded during submission."""

        path = self.paths.run_dir / "provenance" / "repository_state.txt"
        return path.read_text(encoding="utf-8").strip().replace("\n", "; ") if path.exists() else "unavailable"
