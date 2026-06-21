#!/usr/bin/env python3
"""Stage Siletti/CELLxGENE supercluster H5AD files for DIV90 label transfer."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path

import anndata as ad
import pandas as pd
import requests


PROJECT_ROOT_DEFAULT = "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder"
COLLECTION_ID = "283d65eb-dd53-496d-adb7-7570c7caa443"
TARGET_SUPERCLUSTERS = (
    "MGE interneuron",
    "CGE interneuron",
    "LAMP5-LHX6 and Chandelier",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", default=os.environ.get("PROJECT_ROOT", PROJECT_ROOT_DEFAULT))
    parser.add_argument(
        "--outdir",
        default=None,
        help="Default: PROJECT_ROOT/results/siletti_2023_whb_reference_label_transfer/source_cellxgene_superclusters",
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--chunk-size", type=int, default=16 * 1024 * 1024)
    parser.add_argument(
        "--target-superclusters",
        default=os.environ.get("SILETTI_TARGET_SUPERCLUSTERS", ",".join(TARGET_SUPERCLUSTERS)),
        help="Comma-separated CELLxGENE supercluster titles to stage.",
    )
    return parser.parse_args()


def safe_token(value: str) -> str:
    return "".join(ch if ch.isalnum() else "_" for ch in value).strip("_").lower()


def sha256_file(path: Path, chunk_size: int) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch_collection() -> dict:
    url = f"https://api.cellxgene.cziscience.com/curation/v1/collections/{COLLECTION_ID}"
    response = requests.get(url, timeout=60)
    response.raise_for_status()
    return response.json()


def download(url: str, path: Path, expected_size: int | None, chunk_size: int, overwrite: bool) -> dict[str, object]:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not overwrite:
        return {"download_status": "reused_existing", "bytes_written": path.stat().st_size}

    tmp = path.with_suffix(path.suffix + ".tmp")
    if tmp.exists():
        tmp.unlink()

    start = time.time()
    with requests.get(url, stream=True, timeout=120) as response:
        response.raise_for_status()
        with tmp.open("wb") as handle:
            for chunk in response.iter_content(chunk_size=chunk_size):
                if chunk:
                    handle.write(chunk)
    tmp.replace(path)
    elapsed = time.time() - start
    size = path.stat().st_size
    status = "downloaded"
    if expected_size is not None and size != expected_size:
        status = "downloaded_size_mismatch"
    return {"download_status": status, "bytes_written": size, "elapsed_seconds": round(elapsed, 3)}


def main() -> None:
    args = parse_args()
    target_superclusters = tuple(x.strip() for x in args.target_superclusters.split(",") if x.strip())
    if not target_superclusters:
        raise ValueError("No target superclusters requested.")
    project_root = Path(args.project_root)
    outdir = (
        Path(args.outdir)
        if args.outdir
        else project_root
        / "results"
        / "siletti_2023_whb_reference_label_transfer"
        / "source_cellxgene_superclusters"
    )
    source_dir = outdir / "h5ad"
    table_dir = outdir / "tables"
    report_dir = outdir / "reports"
    for directory in [source_dir, table_dir, report_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    collection = fetch_collection()
    (source_dir / "cellxgene_collection_metadata.json").write_text(json.dumps(collection, indent=2, sort_keys=True) + "\n")

    rows = []
    for dataset in collection.get("datasets", []):
        title = dataset.get("title", "")
        if not title.startswith("Supercluster:"):
            continue
        supercluster = title.removeprefix("Supercluster:").strip()
        if supercluster not in target_superclusters:
            continue
        h5ads = [asset for asset in dataset.get("assets", []) if asset.get("filetype") == "H5AD"]
        if len(h5ads) != 1:
            raise ValueError(f"Expected exactly one H5AD asset for {title}, found {len(h5ads)}")
        asset = h5ads[0]
        filename = f"siletti_whb_{safe_token(supercluster)}.h5ad"
        path = source_dir / filename
        row = {
            "collection_id": COLLECTION_ID,
            "collection_name": collection.get("name", ""),
            "collection_version_id": collection.get("collection_version_id", ""),
            "dataset_id": dataset.get("dataset_id", ""),
            "dataset_title": title,
            "supercluster": supercluster,
            "asset_filetype": asset.get("filetype", ""),
            "asset_url": asset.get("url", ""),
            "asset_filesize": asset.get("filesize", pd.NA),
            "local_path": str(path),
        }
        try:
            dl = download(asset["url"], path, int(asset["filesize"]), args.chunk_size, args.overwrite)
            row.update(dl)
            h5 = ad.read_h5ad(path, backed="r")
            try:
                row.update(
                    {
                        "h5ad_n_obs": int(h5.n_obs),
                        "h5ad_n_vars": int(h5.n_vars),
                        "obs_columns": ";".join(map(str, h5.obs.columns)),
                        "var_columns": ";".join(map(str, h5.var.columns)),
                    }
                )
            finally:
                h5.file.close()
            row["sha256"] = sha256_file(path, args.chunk_size)
        except Exception as exc:
            row["download_status"] = "failed"
            row["error"] = repr(exc)
            rows.append(row)
            pd.DataFrame(rows).to_csv(table_dir / "siletti_cellxgene_supercluster_download_manifest.tsv", sep="\t", index=False)
            raise
        rows.append(row)

    manifest = pd.DataFrame(rows)
    if manifest.empty:
        raise RuntimeError("No target Siletti supercluster H5AD assets were found in the CELLxGENE collection.")

    manifest.to_csv(table_dir / "siletti_cellxgene_supercluster_download_manifest.tsv", sep="\t", index=False)
    summary = (
        manifest.groupby("download_status", dropna=False)
        .agg(n_assets=("local_path", "count"), total_bytes=("bytes_written", "sum"))
        .reset_index()
    )
    summary.to_csv(table_dir / "siletti_cellxgene_supercluster_download_summary.tsv", sep="\t", index=False)

    report = [
        "# Siletti CELLxGENE Supercluster H5AD Staging",
        "",
        f"Collection: `{COLLECTION_ID}`",
        f"Output: `{outdir}`",
        "",
        "Target superclusters:",
        "",
    ]
    report.extend(f"- {x}" for x in target_superclusters)
    report.extend(
        [
            "",
            "Files:",
            "",
        ]
    )
    for row in manifest.itertuples():
        report.append(
            f"- `{row.supercluster}`: `{row.local_path}` "
            f"({int(row.bytes_written):,} bytes; {row.h5ad_n_obs:,} cells x {row.h5ad_n_vars:,} genes)"
        )
    (report_dir / "siletti_cellxgene_supercluster_staging_report.md").write_text("\n".join(report) + "\n")

    print(manifest.to_string(index=False))
    print(report_dir / "siletti_cellxgene_supercluster_staging_report.md")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"[siletti-stage] ERROR: {exc}", file=sys.stderr)
        raise
