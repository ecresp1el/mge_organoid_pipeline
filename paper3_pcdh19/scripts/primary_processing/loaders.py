"""Readers for registered metadata and 10x Cell Ranger HDF5 matrices."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path
from typing import Iterable

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from .models import SampleStructure, Step00Paths, Step00Settings, ValidationLedger


def _decode(values: Iterable[bytes]) -> tuple[str, ...]:
    """Decode a one-dimensional HDF5 byte-string array as Unicode strings."""

    return tuple(value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values)


def _feature_signature(
    gene_ids: tuple[str, ...],
    gene_symbols: tuple[str, ...],
    feature_types: tuple[str, ...],
    genomes: tuple[str, ...],
) -> str:
    """Hash the ordered feature definition to compare matrices efficiently."""

    digest = hashlib.sha256()
    for record in zip(gene_ids, gene_symbols, feature_types, genomes):
        digest.update("\t".join(record).encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


class SampleMetadataRegistry:
    """Validate and combine the biological sample key and technical manifest."""

    REQUIRED_SAMPLE_KEY_COLUMNS = (
        "technical_sample_id",
        "submitted_sample_name",
        "target_cells",
        "organism",
        "tissue",
        "region",
        "genotype",
        "sex",
        "design_group",
        "metadata_source",
        "metadata_registered_date",
    )
    REQUIRED_TECHNICAL_COLUMNS = (
        "candidate_id",
        "experiment_id",
        "sample_id",
        "probe_barcode_id",
        "vendor_called_cells",
        "median_umi_counts_per_cell",
        "median_genes_per_cell",
        "n_features_in_delivered_matrix",
        "biological_identity_status",
    )

    def __init__(self, paths: Step00Paths, settings: Step00Settings, ledger: ValidationLedger):
        """Store immutable dependencies used by metadata validation."""

        self.paths = paths
        self.settings = settings
        self.ledger = ledger

    def load(self) -> pd.DataFrame:
        """Return one validated row per technical sample in numeric sample order."""

        sample_key = pd.read_csv(self.paths.sample_key, dtype=str, keep_default_na=False)
        technical = pd.read_csv(self.paths.technical_manifest, sep="\t", dtype=str, keep_default_na=False)
        self._validate_columns(sample_key, self.REQUIRED_SAMPLE_KEY_COLUMNS, "sample_key")
        self._validate_columns(technical, self.REQUIRED_TECHNICAL_COLUMNS, "technical_manifest")
        self._validate_unique(sample_key, "technical_sample_id", "sample_key")
        self._validate_unique(technical, "sample_id", "technical_manifest")
        merged = sample_key.merge(
            technical,
            left_on="technical_sample_id",
            right_on="sample_id",
            how="outer",
            validate="one_to_one",
            indicator=True,
        )
        complete = bool((merged["_merge"] == "both").all())
        self.ledger.add(
            "metadata_join_complete",
            "all_samples",
            complete,
            merged["_merge"].value_counts().to_dict(),
            {"both": self.settings.expected_samples},
            "Biological and technical registries must match one-to-one.",
        )
        merged = merged.loc[merged["_merge"] == "both"].drop(columns=["_merge", "sample_id"])
        merged["sample_number"] = merged["technical_sample_id"].str.rsplit("-", n=1).str[-1].astype(int)
        merged = merged.sort_values("sample_number").reset_index(drop=True)
        self.ledger.add(
            "expected_sample_count",
            "all_samples",
            len(merged) == self.settings.expected_samples,
            len(merged),
            self.settings.expected_samples,
            "Exactly the registered 12 samples must enter Step 00.",
        )
        return merged

    def _validate_columns(self, frame: pd.DataFrame, required: tuple[str, ...], scope: str) -> None:
        """Record whether a metadata table contains all mandatory columns."""

        missing = [column for column in required if column not in frame.columns]
        self.ledger.add(
            "required_columns",
            scope,
            not missing,
            "none" if not missing else ",".join(missing),
            "none missing",
            "Registered metadata schemas are explicit pipeline inputs.",
        )

    def _validate_unique(self, frame: pd.DataFrame, column: str, scope: str) -> None:
        """Record whether a registry key is populated and unique."""

        valid = column in frame.columns and frame[column].ne("").all() and frame[column].is_unique
        self.ledger.add(
            "unique_sample_key",
            scope,
            bool(valid),
            frame[column].nunique() if column in frame.columns else "column missing",
            len(frame),
            f"{column} must be non-empty and unique.",
        )


class CellRangerStructureInspector:
    """Inspect matrix structure, features, and barcodes before expression loading."""

    def __init__(self, paths: Step00Paths):
        """Store the source-path resolver used for all sample inspections."""

        self.paths = paths

    def inspect(self, sample_id: str) -> SampleStructure:
        """Read structural metadata for one filtered matrix and its raw companion."""

        filtered_path = self.paths.filtered_matrix(sample_id)
        raw_path = self.paths.raw_matrix(sample_id)
        barcode_path = self.paths.filtered_barcodes(sample_id)
        metrics_path = self.paths.metrics_summary(sample_id)
        for path in (filtered_path, raw_path, barcode_path, metrics_path):
            if not path.is_file():
                raise FileNotFoundError(f"Required Step 00 input is missing: {path}")
        filtered = self._inspect_h5(filtered_path, include_barcodes=True)
        raw = self._inspect_h5(raw_path, include_barcodes=False)
        csv_barcodes = self._read_barcode_csv(barcode_path)
        if csv_barcodes != filtered["barcodes"]:
            raise ValueError(f"Filtered barcode CSV does not match H5 order for {sample_id}")
        filtered_id_set = set(filtered["gene_ids"])
        raw_matching_indices = [
            index for index, gene_id in enumerate(raw["gene_ids"]) if gene_id in filtered_id_set
        ]
        raw_matched_signature = _feature_signature(
            tuple(raw["gene_ids"][index] for index in raw_matching_indices),
            tuple(raw["gene_symbols"][index] for index in raw_matching_indices),
            tuple(raw["feature_types"][index] for index in raw_matching_indices),
            tuple(raw["genomes"][index] for index in raw_matching_indices),
        )
        return SampleStructure(
            sample_id=sample_id,
            filtered_path=filtered_path,
            raw_path=raw_path,
            barcode_csv_path=barcode_path,
            metrics_path=metrics_path,
            n_cells=filtered["n_barcodes"],
            n_features=filtered["n_features"],
            nnz=filtered["nnz"],
            count_dtype=filtered["count_dtype"],
            barcodes=filtered["barcodes"],
            gene_ids=filtered["gene_ids"],
            gene_symbols=filtered["gene_symbols"],
            feature_types=filtered["feature_types"],
            genomes=filtered["genomes"],
            raw_n_barcodes=raw["n_barcodes"],
            raw_n_features=raw["n_gene_expression_features"],
            raw_n_total_features=raw["n_features"],
            raw_feature_types=tuple(sorted(set(raw["feature_types"]))),
            feature_signature=filtered["feature_signature"],
            raw_feature_signature=raw["feature_signature"],
            raw_contains_filtered_features=filtered_id_set.issubset(set(raw["gene_ids"])),
            raw_preserves_filtered_feature_order=raw_matched_signature == filtered["feature_signature"],
            vendor_metrics_cells=self._read_vendor_cell_count(metrics_path),
        )

    def _inspect_h5(self, path: Path, include_barcodes: bool) -> dict[str, object]:
        """Return shape and ordered feature definitions from one 10x HDF5 matrix."""

        with h5py.File(path, "r") as handle:
            matrix = handle["matrix"]
            features = matrix["features"]
            shape = tuple(int(value) for value in matrix["shape"][:])
            gene_ids = _decode(features["id"][:])
            gene_symbols = _decode(features["name"][:])
            feature_types = _decode(features["feature_type"][:])
            genomes = _decode(features["genome"][:])
            barcodes = _decode(matrix["barcodes"][:]) if include_barcodes else tuple()
            gene_expression_mask = [value == "Gene Expression" for value in feature_types]
            gene_expression_ids = tuple(value for value, keep in zip(gene_ids, gene_expression_mask) if keep)
            gene_expression_symbols = tuple(value for value, keep in zip(gene_symbols, gene_expression_mask) if keep)
            gene_expression_types = tuple(value for value, keep in zip(feature_types, gene_expression_mask) if keep)
            gene_expression_genomes = tuple(value for value, keep in zip(genomes, gene_expression_mask) if keep)
            return {
                "n_features": shape[0],
                "n_barcodes": shape[1],
                "nnz": int(matrix["data"].shape[0]),
                "count_dtype": str(matrix["data"].dtype),
                "barcodes": barcodes,
                "gene_ids": gene_ids,
                "gene_symbols": gene_symbols,
                "feature_types": feature_types,
                "genomes": genomes,
                "feature_signature": _feature_signature(gene_ids, gene_symbols, feature_types, genomes),
                "n_gene_expression_features": len(gene_expression_ids),
                "gene_expression_feature_signature": _feature_signature(
                    gene_expression_ids,
                    gene_expression_symbols,
                    gene_expression_types,
                    gene_expression_genomes,
                ),
            }

    def _read_barcode_csv(self, path: Path) -> tuple[str, ...]:
        """Read Cell Ranger's headerless genome/barcode validation file."""

        with path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.reader(handle))
        if any(len(row) != 2 for row in rows):
            raise ValueError(f"Unexpected filtered-barcode CSV schema: {path}")
        return tuple(row[1] for row in rows)

    def _read_vendor_cell_count(self, path: Path) -> int:
        """Extract the Gene Expression called-cell total from Cell Ranger metrics."""

        frame = pd.read_csv(path, dtype=str, keep_default_na=False)
        selected = frame.loc[
            (frame["Category"] == "Cells")
            & (frame["Library Type"] == "Gene Expression")
            & (frame["Metric Name"] == "Cells"),
            "Metric Value",
        ]
        if len(selected) != 1:
            raise ValueError(f"Could not identify one called-cell metric in {path}")
        return int(selected.iloc[0].replace(",", ""))


class CellRangerSampleReader:
    """Load one validated 10x H5 as a cells-by-genes AnnData object."""

    def __init__(self, settings: Step00Settings):
        """Store the feature-selection contract for expression loading."""

        self.settings = settings

    def read(self, structure: SampleStructure, metadata: pd.Series) -> ad.AnnData:
        """Load sparse counts and attach registered sample-level metadata."""

        with h5py.File(structure.filtered_path, "r") as handle:
            matrix = handle["matrix"]
            feature_mask = np.asarray(structure.feature_types) == self.settings.feature_type
            gene_by_cell = sparse.csc_matrix(
                (matrix["data"][:], matrix["indices"][:], matrix["indptr"][:]),
                shape=(structure.n_features, structure.n_cells),
            )
            if not bool(feature_mask.all()):
                gene_by_cell = gene_by_cell[feature_mask, :]
            cell_by_gene = gene_by_cell.transpose().tocsr()
        gene_ids = np.asarray(structure.gene_ids, dtype=object)[feature_mask]
        gene_symbols = np.asarray(structure.gene_symbols, dtype=object)[feature_mask]
        feature_types = np.asarray(structure.feature_types, dtype=object)[feature_mask]
        genomes = np.asarray(structure.genomes, dtype=object)[feature_mask]
        cell_ids = [f"{structure.sample_id}_{barcode}" for barcode in structure.barcodes]
        obs = pd.DataFrame(index=pd.Index(cell_ids, name="cell_id"))
        obs["cellranger_barcode"] = list(structure.barcodes)
        for column, value in metadata.items():
            obs[column] = str(value)
        var = pd.DataFrame(index=pd.Index(gene_ids, name="gene_id"))
        var["gene_symbol"] = gene_symbols
        var["feature_type"] = feature_types
        var["genome"] = genomes
        return ad.AnnData(X=cell_by_gene, obs=obs, var=var)
