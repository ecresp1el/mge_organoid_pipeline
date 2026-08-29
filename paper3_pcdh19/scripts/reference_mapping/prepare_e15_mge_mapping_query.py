#!/usr/bin/env python3
"""Build a raw-count H5AD and exact metadata contract without recomputing biology."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse


class QueryContractError(RuntimeError):
    pass


def decode(values):
    return [v.decode("utf-8") if isinstance(v, bytes) else str(v) for v in values]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


@dataclass(frozen=True)
class SamplePaths:
    sample_id: str
    matrix: Path
    umap: Path
    clusters: Path


class CellRangerSample:
    def __init__(self, paths: SamplePaths, sample_row: pd.Series):
        self.paths = paths
        self.sample_row = sample_row

    def load(self, expected_ids=None):
        with h5py.File(self.paths.matrix, "r") as src:
            group = src["matrix"]
            gene_ids = decode(group["features/id"][:])
            gene_symbols = decode(group["features/name"][:])
            feature_types = decode(group["features/feature_type"][:])
            barcodes = decode(group["barcodes"][:])
            shape = tuple(int(x) for x in group["shape"][:])
            matrix = sparse.csc_matrix(
                (group["data"][:], group["indices"][:], group["indptr"][:]),
                shape=shape,
            ).transpose().tocsr()
        if len(set(gene_ids)) != len(gene_ids):
            raise QueryContractError(f"Duplicate Ensembl IDs in {self.paths.sample_id}")
        if set(feature_types) != {"Gene Expression"}:
            raise QueryContractError(f"Non-GEX features in {self.paths.sample_id}")
        if expected_ids is not None and gene_ids != expected_ids:
            raise QueryContractError(f"Feature identity/order mismatch in {self.paths.sample_id}")

        umap = pd.read_csv(self.paths.umap)
        clusters = pd.read_csv(self.paths.clusters)
        for table, label in ((umap, "UMAP"), (clusters, "cluster")):
            if table.iloc[:, 0].astype(str).tolist() != barcodes:
                raise QueryContractError(
                    f"{label} barcode identity/order mismatch in {self.paths.sample_id}"
                )
        stable_ids = [f"{self.paths.sample_id}_{barcode}" for barcode in barcodes]
        metadata = pd.DataFrame({
            "cell_id": stable_ids,
            "sample_id": self.paths.sample_id,
            "raw_barcode": barcodes,
            "submitted_sample_name": self.sample_row["submitted_sample_name"],
            "organism": self.sample_row["organism"],
            "tissue": self.sample_row["tissue"],
            "region": self.sample_row["region"],
            "genotype": self.sample_row["genotype"],
            "sex": self.sample_row["sex"],
            "design_group": self.sample_row["design_group"],
            "existing_cluster_numeric": clusters.iloc[:, 1].astype(str).values,
            "existing_cluster": [
                f"{self.paths.sample_id}:{value}" for value in clusters.iloc[:, 1].astype(str)
            ],
            "vendor_umap_1": pd.to_numeric(umap.iloc[:, 1]).values,
            "vendor_umap_2": pd.to_numeric(umap.iloc[:, 2]).values,
        })
        metadata.index = stable_ids
        return matrix, metadata, gene_ids, gene_symbols


class QueryPackageBuilder:
    def __init__(self, cellranger_root: Path, sample_key: Path, output_dir: Path):
        self.cellranger_root = cellranger_root
        self.sample_key = sample_key
        self.output_dir = output_dir

    def sample_paths(self, sample_id: str) -> SamplePaths:
        root = self.cellranger_root / "per_sample_outs" / sample_id
        return SamplePaths(
            sample_id=sample_id,
            matrix=root / "sample_filtered_feature_bc_matrix.h5",
            umap=root / "analysis/umap/gene_expression_2_components/projection.csv",
            clusters=root / "analysis/clustering/gene_expression_graphclust/clusters.csv",
        )

    def run(self):
        self.output_dir.mkdir(parents=True, exist_ok=True)
        sample_key = pd.read_csv(self.sample_key, dtype=str)
        matrices, metadata_frames, inventory = [], [], []
        expected_ids, expected_symbols = None, None
        for _, row in sample_key.iterrows():
            sample_id = row["technical_sample_id"]
            paths = self.sample_paths(sample_id)
            for path in (paths.matrix, paths.umap, paths.clusters):
                if not path.is_file():
                    raise QueryContractError(f"Missing required source: {path}")
            matrix, metadata, gene_ids, gene_symbols = CellRangerSample(paths, row).load(expected_ids)
            if expected_ids is None:
                expected_ids, expected_symbols = gene_ids, gene_symbols
            matrices.append(matrix)
            metadata_frames.append(metadata)
            inventory.append({
                "sample_id": sample_id,
                "cells": matrix.shape[0],
                "genes": matrix.shape[1],
                "nonzero_counts": matrix.nnz,
                "matrix_path": str(paths.matrix),
                "umap_path": str(paths.umap),
                "clusters_path": str(paths.clusters),
            })

        counts = sparse.vstack(matrices, format="csr")
        obs = pd.concat(metadata_frames, axis=0)
        if not obs.index.is_unique or counts.shape[0] != len(obs):
            raise QueryContractError("Combined query cell identity contract failed")
        var = pd.DataFrame({
            "gene_symbol": expected_symbols,
            "feature_type": "Gene Expression",
        }, index=pd.Index(expected_ids, name="ensembl_gene_id"))
        query = ad.AnnData(X=counts, obs=obs.copy(), var=var)
        query.obsm["X_cellranger_umap"] = obs[["vendor_umap_1", "vendor_umap_2"]].to_numpy()
        query.uns["coordinate_warning"] = (
            "Cell Ranger UMAPs were computed independently per sample; compare only within sample/facets."
        )
        query.uns["transformation"] = "none; X contains original integer UMI counts"

        metadata_path = self.output_dir / "query_cell_metadata.tsv.gz"
        h5ad_path = self.output_dir / "query_raw_counts_mapmycells.h5ad"
        inventory_path = self.output_dir / "query_sample_inventory.tsv"
        with gzip.open(metadata_path, "wt", encoding="utf-8", newline="") as handle:
            obs.to_csv(handle, sep="\t", index=False)
        pd.DataFrame(inventory).to_csv(inventory_path, sep="\t", index=False)
        query.write_h5ad(h5ad_path, compression="gzip")
        del query, counts, matrices

        manifest = {
            "cells": int(len(obs)),
            "genes": int(len(expected_ids)),
            "samples": int(len(sample_key)),
            "cell_ids_unique": bool(obs.index.is_unique),
            "raw_counts_h5ad": str(h5ad_path),
            "raw_counts_h5ad_bytes": h5ad_path.stat().st_size,
            "raw_counts_h5ad_sha256": sha256_file(h5ad_path),
            "metadata": str(metadata_path),
            "metadata_sha256": sha256_file(metadata_path),
            "existing_coordinates": "Cell Ranger per-sample UMAP; not recomputed",
            "existing_clusters": "Cell Ranger per-sample graph clusters; not recomputed",
        }
        (self.output_dir / "query_contract.json").write_text(
            json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cellranger-root", required=True, type=Path)
    parser.add_argument("--sample-key", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    QueryPackageBuilder(args.cellranger_root, args.sample_key, args.output_dir).run()


if __name__ == "__main__":
    main()
