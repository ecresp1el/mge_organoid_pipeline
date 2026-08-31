"""Losslessly bridge the approved H5AD counts to native R scDblFinder."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd

from .models import ValidationLedger
from .step03_models import Step03Paths, Step03Settings
from .step03_validation import Step03BridgeValidator, Step03InputValidator


class TenXBridgeWriter:
    """Write an exact 10x-schema view of H5AD CSR arrays for native R."""

    def __init__(self, source_h5ad: Path, destination_h5: Path):
        """Store the source checkpoint and run-scoped destination paths."""

        self.source_h5ad = source_h5ad
        self.destination_h5 = destination_h5

    def write(self, adata: ad.AnnData) -> None:
        """Copy sparse arrays chunkwise and write ordered cell/feature names."""

        if self.destination_h5.exists():
            raise FileExistsError(f"Refusing to replace bridge file: {self.destination_h5}")
        string_dtype = h5py.string_dtype(encoding="utf-8")
        with h5py.File(self.source_h5ad, "r") as source, h5py.File(self.destination_h5, "w") as target:
            matrix = target.create_group("matrix")
            for name in ("data", "indices", "indptr"):
                self._copy_dataset(source["X"][name], matrix, name)
            matrix.create_dataset("shape", data=np.asarray([adata.n_vars, adata.n_obs], dtype=np.int64))
            matrix.create_dataset("barcodes", data=adata.obs_names.astype(str).to_numpy(), dtype=string_dtype)
            features = matrix.create_group("features")
            features.create_dataset("id", data=adata.var_names.astype(str).to_numpy(), dtype=string_dtype)
            features.create_dataset("name", data=adata.var["gene_symbol"].astype(str).to_numpy(), dtype=string_dtype)
            features.create_dataset("feature_type", data=adata.var["feature_type"].astype(str).to_numpy(), dtype=string_dtype)
            features.create_dataset("genome", data=adata.var["genome"].astype(str).to_numpy(), dtype=string_dtype)
            features.create_dataset("_all_tag_keys", data=np.asarray(["genome"], dtype=object), dtype=string_dtype)
            matrix.attrs["encoding-type"] = "csc_matrix"
            matrix.attrs["encoding-version"] = "0.1.0"

    @staticmethod
    def _copy_dataset(source: h5py.Dataset, group: h5py.Group, name: str) -> None:
        """Copy one large sparse array in bounded chunks without conversion."""

        chunk = min(max(1, 8 * 1024 * 1024 // source.dtype.itemsize), len(source))
        destination = group.create_dataset(name, shape=source.shape, dtype=source.dtype, chunks=(chunk,), compression="lzf")
        for start in range(0, len(source), chunk):
            destination[start : start + chunk] = source[start : start + chunk]


class Step03PreparationWorkflow:
    """Validate approved input and prepare the frozen native-R bridge assets."""

    def __init__(self, paths: Step03Paths, settings: Step03Settings):
        """Construct focused validation and bridge collaborators."""

        self.paths = paths
        self.settings = settings
        self.ledger = ValidationLedger()

    def run(self) -> Path:
        """Create exact bridge inputs without changing cells, genes, or counts."""

        self.paths.intermediate_dir.mkdir(parents=True, exist_ok=False)
        validator = Step03InputValidator(self.paths, self.settings, self.ledger)
        validator.validate_approval()
        adata = ad.read_h5ad(self.paths.input_h5ad, backed="r")
        try:
            validator.validate_object(adata)
            metadata = adata.obs[[self.settings.sample_field, self.settings.design_field, "total_counts", "n_genes_by_counts", "pct_counts_mt"]].copy()
            metadata.insert(0, "cell_id", adata.obs_names.astype(str))
            metadata.insert(1, "capture_id", self.settings.capture_id)
            metadata.to_csv(self.paths.cell_metadata, sep="\t", index=False)
            TenXBridgeWriter(self.paths.input_h5ad, self.paths.bridge_h5).write(adata)
        finally:
            adata.file.close()
        Step03BridgeValidator(self.paths, self.settings, self.ledger).validate()
        checks = self.ledger.to_frame()
        checks.to_csv(self.paths.intermediate_dir / "prepare_validation_checks.tsv", sep="\t", index=False)
        self.ledger.require_all_pass()
        return self.paths.bridge_h5
