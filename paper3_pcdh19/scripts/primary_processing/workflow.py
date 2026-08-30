"""Orchestration for Step 00 canonical AnnData construction and publication."""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import pandas as pd

from .loaders import CellRangerSampleReader, CellRangerStructureInspector, SampleMetadataRegistry
from .models import Step00Paths, Step00Settings, ValidationLedger
from .publishing import ApprovalLedger, AtomicRunPublisher, ProvenancePublisher
from .validation import DocumentationAuditor, Step00Validator


class Step00Workflow:
    """Coordinate registered inputs, validation, AnnData creation, and publication."""

    def __init__(self, paths: Step00Paths, settings: Step00Settings):
        """Construct independent collaborators around immutable paths/settings."""

        self.paths = paths
        self.settings = settings
        self.ledger = ValidationLedger()
        self.metadata_registry = SampleMetadataRegistry(paths, settings, self.ledger)
        self.inspector = CellRangerStructureInspector(paths)
        self.reader = CellRangerSampleReader(settings)
        self.validator = Step00Validator(settings, self.ledger)
        self.provenance = ProvenancePublisher()

    def run(self) -> Path:
        """Build, serialize, reopen, validate, and atomically publish Step 00."""

        publisher = AtomicRunPublisher(self.paths.run_dir)
        try:
            metadata = self.metadata_registry.load()
            documentation = pd.DataFrame(DocumentationAuditor().audit(self.paths.frozen_code_dir))
            documentation_passed = bool((documentation["status"] == "PASS").all())
            self.ledger.add("documented_definitions", "frozen_code", documentation_passed, int((documentation["status"] == "PASS").sum()), len(documentation), "Every frozen Python module, class, and function must have a docstring.")
            structures = []
            sample_objects = []
            reference_signature = ""
            sample_rows = []
            for row in metadata.itertuples(index=False):
                sample_id = row.technical_sample_id
                structure = self.inspector.inspect(sample_id)
                if not reference_signature:
                    reference_signature = structure.feature_signature
                self.validator.validate_structure(structure, int(row.vendor_called_cells), reference_signature)
                sample_metadata = metadata.loc[metadata["technical_sample_id"] == sample_id].iloc[0]
                sample_object = self.reader.read(structure, sample_metadata)
                self.validator.validate_loaded_sample(sample_object, structure)
                structures.append(structure)
                sample_objects.append(sample_object)
                sample_rows.append({
                    "technical_sample_id": sample_id,
                    "submitted_sample_name": row.submitted_sample_name,
                    "probe_barcode_id": row.probe_barcode_id,
                    "genotype": row.genotype,
                    "sex": row.sex,
                    "design_group": row.design_group,
                    "n_cells": structure.n_cells,
                    "n_genes": structure.n_features,
                    "nnz": structure.nnz,
                    "filtered_matrix": str(structure.filtered_path),
                    "raw_matrix_registered_for_step04": str(structure.raw_path),
                    "raw_matrix_barcodes": structure.raw_n_barcodes,
                    "raw_matrix_total_features": structure.raw_n_total_features,
                    "raw_matrix_gene_expression_features": structure.raw_n_features,
                    "raw_matrix_feature_types": ";".join(structure.raw_feature_types),
                    "feature_signature_sha256": structure.feature_signature,
                })
            combined = ad.concat(sample_objects, axis=0, join="inner", merge="same", uns_merge="same", index_unique=None)
            combined.var = sample_objects[0].var.copy()
            combined.uns["primary_processing"] = {
                "step": "00_input_validation_and_canonical_anndata",
                "matrix_state": "sparse_unnormalized_integer_cellranger_counts",
                "canonical_input": "per-sample sample_filtered_feature_bc_matrix.h5",
                "raw_matrix_policy": "registered_for_step04_not_loaded_into_canonical_object",
                "feature_policy": self.settings.feature_type,
                "genome": self.settings.expected_genome,
                "cell_id_policy": "<technical_sample_id>_<cellranger_barcode>",
                "var_names_policy": "Ensembl gene ID; gene symbol in var['gene_symbol']",
                "forbidden_operations": ["QC filtering", "normalization", "scaling", "PCA", "neighbors", "UMAP", "clustering", "annotation"],
            }
            self.validator.validate_combined(combined)
            per_sample_total = sum(item.n_cells for item in structures)
            self.ledger.add("sample_sum_matches_combined", "combined", per_sample_total == combined.n_obs, per_sample_total, combined.n_obs, "Combined cells must equal the sum of per-sample called cells.")
            checkpoint_name = "objects/pcdh19_step00_canonical_raw_counts.h5ad"
            checkpoint = publisher.stage_dir / checkpoint_name
            combined.write_h5ad(checkpoint, compression=self.settings.compression)
            reopened = ad.read_h5ad(checkpoint, backed="r")
            self.ledger.add("h5ad_roundtrip_shape", "combined", reopened.shape == combined.shape, reopened.shape, combined.shape, "Published H5AD must reopen with unchanged dimensions.")
            self.ledger.add("h5ad_roundtrip_obs", "combined", reopened.obs_names.equals(combined.obs_names), reopened.n_obs, combined.n_obs, "Published H5AD must preserve ordered cell identifiers.")
            self.ledger.add("h5ad_roundtrip_var", "combined", reopened.var_names.equals(combined.var_names), reopened.n_vars, combined.n_vars, "Published H5AD must preserve ordered Ensembl identifiers.")
            reopened.file.close()
            checks = self.ledger.to_frame()
            checks.to_csv(publisher.tables_dir / "validation_checks.tsv", sep="\t", index=False)
            documentation.to_csv(publisher.tables_dir / "documentation_audit.tsv", sep="\t", index=False)
            metadata.to_csv(publisher.tables_dir / "registered_sample_metadata.tsv", sep="\t", index=False)
            sample_summary = pd.DataFrame(sample_rows)
            sample_summary.to_csv(publisher.tables_dir / "sample_input_and_object_summary.tsv", sep="\t", index=False)
            self._write_metadata_dictionary(combined, publisher.tables_dir / "anndata_metadata_dictionary.tsv")
            self._write_object_summary(combined, structures, publisher.tables_dir / "canonical_object_summary.tsv")
            self.provenance.software_versions().to_csv(publisher.tables_dir / "software_versions.tsv", sep="\t", index=False)
            self.ledger.require_all_pass()
            run_id = self.paths.run_dir.name
            code_version = self._code_version()
            status = self.provenance.status_frame(run_id, checkpoint_name, combined.n_obs, combined.n_vars)
            status.to_csv(publisher.stage_dir / "STEP_STATUS.tsv", sep="\t", index=False)
            report = self.provenance.report(run_id, checkpoint_name, combined.n_obs, combined.n_vars, int(combined.X.nnz), sample_summary, checks)
            (publisher.stage_dir / "STEP00_VALIDATION_REPORT.md").write_text(report, encoding="utf-8")
            manifest = self.provenance.output_manifest(publisher.stage_dir)
            manifest.to_csv(publisher.tables_dir / "output_manifest.tsv", sep="\t", index=False)
            publisher.publish()
            ApprovalLedger().update_in_review(
                self.paths.workflow_root / "APPROVAL_LEDGER.tsv",
                run_id,
                checkpoint_name,
                combined.n_obs,
                combined.n_vars,
                code_version,
            )
            return self.paths.run_dir / checkpoint_name
        except Exception:
            publisher.discard()
            raise

    def _write_metadata_dictionary(self, adata: ad.AnnData, path: Path) -> None:
        """Describe every observation and variable field in the canonical object."""

        descriptions = {
            "cellranger_barcode": "Original Cell Ranger barcode, unique within technical sample.",
            "technical_sample_id": "Registered 15662-JZ technical sample identifier.",
            "submitted_sample_name": "Submitted sample name from the registered biological key.",
            "genotype": "Registered sample-level genotype; not inferred from expression.",
            "sex": "Registered sample-level sex.",
            "design_group": "Registered genotype/sex design group.",
            "gene_symbol": "Cell Ranger gene symbol corresponding to the Ensembl var index.",
            "feature_type": "Cell Ranger feature type; Step 00 retains Gene Expression only.",
            "genome": "Cell Ranger genome annotation.",
        }
        rows = []
        for axis, frame in (("obs", adata.obs), ("var", adata.var)):
            index_name = frame.index.name or f"{axis}_index"
            rows.append({"axis": axis, "field": index_name, "dtype": "index", "description": "Globally unique prefixed cell ID." if axis == "obs" else "Unique Ensembl gene identifier."})
            for column in frame.columns:
                rows.append({"axis": axis, "field": column, "dtype": str(frame[column].dtype), "description": descriptions.get(column, "Registered source metadata preserved without biological inference.")})
        pd.DataFrame(rows).to_csv(path, sep="\t", index=False)

    def _write_object_summary(self, adata: ad.AnnData, structures: list, path: Path) -> None:
        """Record the exact matrix state and prohibited Step 00 transformations."""

        summary = {
            "n_cells": adata.n_obs,
            "n_genes": adata.n_vars,
            "n_samples": len(structures),
            "nnz": int(adata.X.nnz),
            "matrix_class": type(adata.X).__name__,
            "matrix_dtype": str(adata.X.dtype),
            "X_state": "sparse unnormalized integer Cell Ranger counts",
            "counts_layer": "not duplicated; raw counts are in X",
            "raw_slot": "not populated",
            "obsm_keys": json.dumps(list(adata.obsm.keys())),
            "obsp_keys": json.dumps(list(adata.obsp.keys())),
            "normalization": "none",
            "filtering": "none",
            "annotation": "none",
        }
        pd.DataFrame([summary]).to_csv(path, sep="\t", index=False)

    def _code_version(self) -> str:
        """Read the frozen repository identity recorded at submission time."""

        provenance = self.paths.run_dir / "provenance" / "repository_state.txt"
        return provenance.read_text(encoding="utf-8").strip().replace("\n", "; ") if provenance.exists() else "unavailable"
