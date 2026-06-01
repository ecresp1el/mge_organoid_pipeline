#!/usr/bin/env python
"""Prove expected CellBender output locations for Notebook 00 samples.

This proof is intentionally tolerant of missing outputs: it reports found and
missing files separately because CellBender may have been run for only a subset
of samples.
"""

from __future__ import annotations

import argparse
import os
import socket
import sys
from pathlib import Path

import h5py
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = REPO_ROOT / "python_notebooks" / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from mge_organoid_python.data_sources import Notebook00SourceConfig, cellbender_output_table


DEFAULT_SAMPLE_IDS = (
    "9853-MW-1",
    "9853-MW-2",
    "9853-MW-3",
    "9853-MW-4",
    "9853-MW-5",
    "9853-MW-6",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root",
        default=str(REPO_ROOT),
        help="Git checkout root containing metadata/. Default: current repository.",
    )
    parser.add_argument(
        "--data-root",
        default=os.environ.get("MGE_DATA_ROOT") or os.environ.get("PROJECT_ROOT"),
        help="Large runtime data root. Defaults to MGE_DATA_ROOT, then PROJECT_ROOT.",
    )
    parser.add_argument(
        "--sample-id",
        action="append",
        dest="sample_ids",
        help="run_sample_id to check. May be repeated. Default: Notebook 00 DIV30 sample IDs 1-6.",
    )
    parser.add_argument(
        "--target-div",
        action="append",
        dest="target_divs",
        help="DIV value to include. May be repeated. Default: DIV30.",
    )
    parser.add_argument(
        "--report-tsv",
        default=None,
        help="Optional TSV path for the proof summary.",
    )
    parser.add_argument(
        "--require-any-existing",
        action="store_true",
        help="Fail if none of the expected CellBender outputs exist.",
    )
    return parser.parse_args()


def inspect_h5(path: str) -> dict:
    """Return lightweight HDF5 metadata without loading the expression matrix."""
    output_path = Path(path)
    if not output_path.exists():
        return {
            "h5_opened": False,
            "h5_top_level_keys": "",
            "h5_matrix_shape": "",
            "h5_inspect_error": "",
        }

    try:
        with h5py.File(output_path, "r") as handle:
            top_level_keys = sorted(handle.keys())
            matrix_shape = ""
            if "matrix" in handle and "shape" in handle["matrix"]:
                matrix_shape = "x".join(str(int(value)) for value in handle["matrix"]["shape"][()])
            return {
                "h5_opened": True,
                "h5_top_level_keys": ",".join(top_level_keys),
                "h5_matrix_shape": matrix_shape,
                "h5_inspect_error": "",
            }
    except Exception as exc:
        return {
            "h5_opened": False,
            "h5_top_level_keys": "",
            "h5_matrix_shape": "",
            "h5_inspect_error": repr(exc),
        }


def main() -> int:
    args = parse_args()
    if not args.data_root:
        raise EnvironmentError("Set PROJECT_ROOT or MGE_DATA_ROOT, or pass --data-root.")

    config = Notebook00SourceConfig.from_defaults(
        data_source="cellbender_denoised",
        repo_root=args.repo_root,
        data_root=args.data_root,
        target_divs=tuple(args.target_divs or ("DIV30",)),
        target_run_sample_ids=tuple(args.sample_ids or DEFAULT_SAMPLE_IDS),
        strict_missing_matrix_dirs=False,
    )

    print("Notebook 00 CellBender location proof")
    print("host:", socket.gethostname())
    print("python:", sys.executable)
    print("repo_root:", config.repo_root)
    print("data_root:", config.data_root)
    print("clean_adata_dir:", config.data_root / "clean_adata")
    print("sample_ids:", list(config.target_run_sample_ids or []))
    print("target_divs:", list(config.target_divs))

    proof_df = cellbender_output_table(config)
    h5_inspection_df = pd.DataFrame([inspect_h5(path) for path in proof_df["cellbender_output_h5"]])
    proof_df = pd.concat([proof_df.reset_index(drop=True), h5_inspection_df], axis=1)

    found_n = int(proof_df["cellbender_output_exists"].sum())
    missing_n = int((~proof_df["cellbender_output_exists"]).sum())

    print("\nCellBender expected path summary:")
    print(proof_df.to_string(index=False))

    print("\nCounts:")
    print("cellbender_outputs_found:", found_n)
    print("cellbender_outputs_missing:", missing_n)

    if args.report_tsv:
        report_path = Path(args.report_tsv).expanduser().resolve()
        report_path.parent.mkdir(parents=True, exist_ok=True)
        proof_df.to_csv(report_path, sep="\t", index=False)
        print("\nWrote report:", report_path)

    failures = []
    unopened_existing = proof_df[
        proof_df["cellbender_output_exists"] & ~proof_df["h5_opened"]
    ][["run_sample_id", "cellbender_output_h5", "h5_inspect_error"]]
    if not unopened_existing.empty:
        failures.append("Some existing CellBender outputs could not be opened as HDF5.")
        print("\nExisting outputs that failed HDF5 inspection:")
        print(unopened_existing.to_string(index=False))

    if args.require_any_existing and found_n == 0:
        failures.append("No expected CellBender outputs exist.")

    if failures:
        print("\nFAILURES:")
        for failure in failures:
            print(" -", failure)
        return 1

    print("\nPASS: expected CellBender locations were derived and existing outputs were HDF5-readable.")
    if missing_n:
        print("NOTE: some expected outputs are missing; this is allowed by this proof.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
