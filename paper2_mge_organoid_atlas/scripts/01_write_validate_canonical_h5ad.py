#!/usr/bin/env python3
"""Write and validate one canonical H5AD from a reloaded minimal RDS bridge."""

from __future__ import annotations

import argparse
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import platform
import shutil

import anndata as ad
import numpy as np
import pandas as pd
from scipy import io as scipy_io
from scipy import sparse


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study-id", required=True)
    parser.add_argument("--study-dir", type=Path, required=True)
    return parser.parse_args()


def read_schema(path: Path, table: str) -> dict[str, str]:
    schema = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    rows = schema.loc[schema["table"] == table]
    return dict(zip(rows["column"], rows["r_class"], strict=True))


def read_metadata(path: Path, schema: dict[str, str], index_column: str) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if list(frame.columns) != list(schema):
        raise ValueError(f"Metadata columns do not match schema for {path}")
    for column, r_class in schema.items():
        if r_class in {"numeric", "integer"}:
            values = pd.to_numeric(frame[column].replace("", np.nan), errors="raise")
            frame[column] = values.astype("Int64") if r_class == "integer" else values.astype(float)
        elif r_class == "logical":
            mapping = {"TRUE": True, "FALSE": False, "": pd.NA}
            unknown = sorted(set(frame[column]).difference(mapping))
            if unknown:
                raise ValueError(f"Unexpected logical values in {column}: {unknown}")
            frame[column] = frame[column].map(mapping).astype("boolean")
        else:
            frame[column] = frame[column].astype(str)
    if frame[index_column].duplicated().any() or (frame[index_column] == "").any():
        raise ValueError(f"{index_column} must be nonempty and unique")
    return frame.set_index(index_column, drop=False)


def canonical_sparse(matrix: sparse.spmatrix) -> sparse.csr_matrix:
    out = matrix.tocsr(copy=True)
    out.sum_duplicates()
    out.sort_indices()
    out.eliminate_zeros()
    return out


def matrix_sha256(matrix: sparse.spmatrix) -> str:
    matrix = canonical_sparse(matrix)
    digest = hashlib.sha256()
    digest.update(np.asarray(matrix.shape, dtype="<i8").tobytes())
    digest.update(np.asarray(matrix.indptr, dtype="<i8").tobytes())
    digest.update(np.asarray(matrix.indices, dtype="<i8").tobytes())
    digest.update(np.asarray(matrix.data, dtype="<f8").tobytes())
    return digest.hexdigest()


def matrix_row(name: str, matrix: sparse.spmatrix) -> dict[str, object]:
    matrix = canonical_sparse(matrix)
    data = matrix.data.astype(np.float64, copy=False)
    return {
        "matrix": name,
        "n_cells": matrix.shape[0],
        "n_features": matrix.shape[1],
        "nnz": matrix.nnz,
        "value_sum": float(data.sum(dtype=np.float64)),
        "value_sum_squares": float(np.square(data).sum(dtype=np.float64)),
        "value_min_nonzero": float(data.min()) if data.size else np.nan,
        "value_max": float(data.max()) if data.size else 0.0,
        "matrix_sha256": matrix_sha256(matrix),
    }


def assert_sparse_equal(expected: sparse.spmatrix, observed: sparse.spmatrix, label: str) -> None:
    left = canonical_sparse(expected)
    right = canonical_sparse(observed)
    if left.shape != right.shape:
        raise ValueError(f"{label} shape mismatch: {left.shape} != {right.shape}")
    difference = left - right
    difference.eliminate_zeros()
    if difference.nnz:
        max_abs = float(np.max(np.abs(difference.data)))
        raise ValueError(f"{label} differs at {difference.nnz} entries; max abs diff={max_abs}")


def normalized_value(value: object) -> object:
    if value is None or value is pd.NA:
        return None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return float(value)
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    return str(value)


def metadata_sha256(frame: pd.DataFrame) -> str:
    digest = hashlib.sha256()
    digest.update("\t".join(map(str, frame.columns)).encode())
    digest.update(b"\n")
    for row in frame.itertuples(index=True, name=None):
        payload = json.dumps([normalized_value(x) for x in row], ensure_ascii=False, separators=(",", ":"))
        digest.update(payload.encode())
        digest.update(b"\n")
    return digest.hexdigest()


def main() -> None:
    args = parse_args()
    study_id = args.study_id
    study_dir = args.study_dir.resolve()
    bridge = study_dir / ".bridge"
    if bridge.parent != study_dir or bridge.name != ".bridge":
        raise ValueError("Unsafe bridge path")
    required = [
        bridge / "counts.mtx",
        bridge / "obs.tsv",
        bridge / "var.tsv",
        bridge / "metadata_schema.tsv",
        study_dir / f"{study_id}_minimal.rds",
        study_dir / "canonical_manifest.tsv",
        study_dir / "matrix_summary_rds.tsv",
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing canonical bridge inputs: {missing}")

    print(f"[canonical-h5ad] reading RDS-derived bridge for {study_id}", flush=True)
    counts_gene_cell = canonical_sparse(scipy_io.mmread(bridge / "counts.mtx"))
    counts = canonical_sparse(counts_gene_cell.transpose())
    lognorm_path = bridge / "lognorm.mtx"
    lognorm = canonical_sparse(scipy_io.mmread(lognorm_path).transpose()) if lognorm_path.is_file() else None

    obs_schema = read_schema(bridge / "metadata_schema.tsv", "obs")
    var_schema = read_schema(bridge / "metadata_schema.tsv", "var")
    obs = read_metadata(bridge / "obs.tsv", obs_schema, "canonical_cell_id")
    var = read_metadata(bridge / "var.tsv", var_schema, "canonical_feature_id")
    if counts.shape != (obs.shape[0], var.shape[0]):
        raise ValueError(f"Counts shape {counts.shape} does not match obs/var {(obs.shape[0], var.shape[0])}")
    if not obs.index.str.startswith(f"{study_id}::").all():
        raise ValueError("Canonical cell IDs do not carry the expected study prefix")

    adata = ad.AnnData(X=counts, obs=obs, var=var)
    if lognorm is not None:
        if lognorm.shape != counts.shape:
            raise ValueError("lognorm/counts shape mismatch")
        adata.layers["lognorm"] = lognorm
    software_versions = {
        "Python": platform.python_version(),
        "Scanpy": version("scanpy"),
        "AnnData": version("anndata"),
        "NumPy": version("numpy"),
        "SciPy": version("scipy"),
        "pandas": version("pandas"),
        "h5py": version("h5py"),
    }
    adata.uns["schema_name"] = "paper2_mge_organoid_canonical_anndata"
    adata.uns["schema_version"] = "1.0.0"
    adata.uns["study_id"] = study_id
    adata.uns["X_semantics"] = "raw_counts"
    adata.uns["lognorm_semantics"] = "source RNA data layer; absent when the source RNA assay has no data layer"
    adata.uns["analysis_specific_embeddings_removed"] = True
    adata.uns["software_versions"] = software_versions

    h5ad_path = study_dir / f"{study_id}_minimal.h5ad"
    print(f"[canonical-h5ad] writing {h5ad_path}", flush=True)
    adata.write_h5ad(h5ad_path, compression="gzip")

    print("[canonical-h5ad] reopening and validating exact equivalence", flush=True)
    observed = ad.read_h5ad(h5ad_path)
    if list(observed.obs_names) != list(obs.index):
        raise ValueError("H5AD cell IDs/order differ from the RDS bridge")
    if list(observed.var_names) != list(var.index):
        raise ValueError("H5AD feature IDs/order differ from the RDS bridge")
    pd.testing.assert_frame_equal(
        obs,
        observed.obs,
        check_dtype=False,
        check_categorical=False,
        check_exact=False,
        rtol=1e-14,
        atol=1e-14,
    )
    pd.testing.assert_frame_equal(
        var,
        observed.var,
        check_dtype=False,
        check_categorical=False,
        check_exact=False,
        rtol=1e-14,
        atol=1e-14,
    )
    if observed.obsm or observed.varm or observed.obsp or observed.varp:
        raise ValueError("Minimal H5AD unexpectedly contains embeddings or pairwise analysis artifacts")
    assert_sparse_equal(counts, observed.X, "counts/X")
    if lognorm is None:
        if "lognorm" in observed.layers:
            raise ValueError("Unexpected lognorm layer in H5AD")
    else:
        if "lognorm" not in observed.layers:
            raise ValueError("Expected lognorm layer is absent from H5AD")
        assert_sparse_equal(lognorm, observed.layers["lognorm"], "lognorm")

    rows = [matrix_row("counts", observed.X)]
    if lognorm is not None:
        rows.append(matrix_row("lognorm", observed.layers["lognorm"]))
    pd.DataFrame(rows).to_csv(study_dir / "matrix_summary_h5ad.tsv", sep="\t", index=False)
    validation = pd.DataFrame(
        [
            ("study_id", study_id),
            ("rds_exists", "PASS"),
            ("h5ad_reopened", "PASS"),
            ("cell_ids_and_order", "PASS"),
            ("feature_ids_and_order", "PASS"),
            ("cell_metadata_values", "PASS"),
            ("feature_metadata_values", "PASS"),
            ("counts_exact", "PASS"),
            ("lognorm_exact", "PASS" if lognorm is not None else "NOT_PRESENT_IN_SOURCE_RNA"),
            ("obs_rows", str(obs.shape[0])),
            ("obs_columns", str(obs.shape[1])),
            ("var_rows", str(var.shape[0])),
            ("var_columns", str(var.shape[1])),
            ("obs_sha256", metadata_sha256(observed.obs)),
            ("var_sha256", metadata_sha256(observed.var)),
            ("analysis_embeddings_absent", "PASS"),
            ("analysis_pairwise_arrays_absent", "PASS"),
        ],
        columns=["check", "value"],
    )
    validation.to_csv(study_dir / "rds_h5ad_equivalence.tsv", sep="\t", index=False)
    pd.DataFrame(
        [{"software": name, "version": package_version} for name, package_version in software_versions.items()]
    ).to_csv(study_dir / "python_environment.tsv", sep="\t", index=False)

    # The bridge is explicitly temporary, lives inside this one study package,
    # and is removed only after the final H5AD has reopened and validated.
    shutil.rmtree(bridge)
    (study_dir / "SUCCESS.txt").write_text(
        f"PASS\nstudy_id={study_id}\nrds_h5ad_equivalence=PASS\n",
        encoding="utf-8",
    )
    print(f"[canonical-h5ad] complete: {study_id}", flush=True)


if __name__ == "__main__":
    main()
