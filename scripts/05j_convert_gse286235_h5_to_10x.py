#!/usr/bin/env python3

"""Convert GSE286235 Cell Ranger H5 files to standard 10x directories.

Seurat::Read10X_h5 requires the R package hdf5r, which is not always installed
in cluster Seurat modules. This converter uses Python h5py/scipy and writes
matrix.mtx.gz, features.tsv.gz, and barcodes.tsv.gz so Seurat can use Read10X.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import h5py
import scipy.io
import scipy.sparse


SAMPLE_INFO = {
    "GSM8721440": {"sample_id": "BF_H9_D36", "disease_status": "healthy"},
    "GSM8721441": {"sample_id": "BF_H9_D63", "disease_status": "healthy"},
    "GSM8721442": {"sample_id": "BFCO_IMR_D63", "disease_status": "healthy"},
    "GSM8721443": {"sample_id": "BF_2DS3_D63", "disease_status": "DS"},
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5-dir", required=True, help="Directory containing GSE286235 H5 files.")
    parser.add_argument("--out-dir", required=True, help="Output root for 10x directories.")
    parser.add_argument(
        "--include-ds",
        action="store_true",
        help="Also convert the DS sample GSM8721443. Default converts healthy samples only.",
    )
    parser.add_argument("--force", action="store_true", help="Reconvert even if _SUCCESS exists.")
    return parser.parse_args()


def decode(values) -> list[str]:
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in values]


def accession_from_name(path: Path) -> str:
    return path.name.split("_", 1)[0]


def write_gzip_lines(path: Path, rows: list[str]) -> None:
    with gzip.open(path, "wt") as handle:
        for row in rows:
            handle.write(row)
            handle.write("\n")


def convert_one(h5_path: Path, sample_id: str, out_root: Path, force: bool) -> None:
    out_dir = out_root / sample_id
    success = out_dir / "_SUCCESS"
    if success.exists() and not force:
        print(f"Found existing conversion, skipping: {out_dir}")
        return

    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Converting {h5_path.name} -> {out_dir}")
    with h5py.File(h5_path, "r") as handle:
        group = handle["matrix"]
        shape = tuple(int(x) for x in group["shape"][:])
        data = group["data"][:]
        indices = group["indices"][:]
        indptr = group["indptr"][:]
        matrix = scipy.sparse.csc_matrix((data, indices, indptr), shape=shape)

        feature_group = group["features"]
        feature_ids = decode(feature_group["id"][:])
        feature_names = decode(feature_group["name"][:])
        if "feature_type" in feature_group:
            feature_types = decode(feature_group["feature_type"][:])
        else:
            feature_types = ["Gene Expression"] * len(feature_ids)
        barcodes = decode(group["barcodes"][:])

    with gzip.open(out_dir / "matrix.mtx.gz", "wb") as handle:
        scipy.io.mmwrite(handle, matrix, field="integer")
    write_gzip_lines(
        out_dir / "features.tsv.gz",
        [f"{fid}\t{name}\t{ftype}" for fid, name, ftype in zip(feature_ids, feature_names, feature_types)],
    )
    write_gzip_lines(out_dir / "barcodes.tsv.gz", barcodes)
    success.write_text("ok\n")


def main() -> None:
    args = parse_args()
    h5_dir = Path(args.h5_dir)
    out_dir = Path(args.out_dir)
    h5_files = sorted(h5_dir.glob("GSM872144*_raw_feature_bc_matrix.h5"))
    if not h5_files:
        raise SystemExit(f"No GSE286235 H5 files found in {h5_dir}")

    for h5_path in h5_files:
        accession = accession_from_name(h5_path)
        info = SAMPLE_INFO.get(accession)
        if info is None:
            print(f"Skipping unrecognized file: {h5_path.name}")
            continue
        if info["disease_status"] != "healthy" and not args.include_ds:
            print(f"Skipping DS sample by default: {h5_path.name}")
            continue
        convert_one(h5_path, info["sample_id"], out_dir, args.force)


if __name__ == "__main__":
    main()
