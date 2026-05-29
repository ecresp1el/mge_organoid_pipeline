#!/usr/bin/env python
"""Prove Notebook 00 can switch between Cell Ranger raw and filtered sources.

This is intended to run on a Great Lakes compute node with the
`mge-organoid-python` conda environment active.
"""

from __future__ import annotations

import argparse
import json
import os
import socket
import sys
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = REPO_ROOT / "python_notebooks" / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from mge_organoid_python.data_sources import Notebook00SourceConfig, load_dataset, sample_table


DEFAULT_SAMPLE_IDS = ("9853-MW-1",)


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
        help="run_sample_id to prove. May be repeated. Default: 9853-MW-1.",
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
    return parser.parse_args()


def nnz_of(adata) -> int:
    x = adata.X
    return int(x.nnz) if hasattr(x, "nnz") else int((x != 0).sum())


def load_one_source(args: argparse.Namespace, data_source: str):
    config = Notebook00SourceConfig.from_defaults(
        data_source=data_source,
        repo_root=args.repo_root,
        data_root=args.data_root,
        target_divs=tuple(args.target_divs or ("DIV30",)),
        target_run_sample_ids=tuple(args.sample_ids or DEFAULT_SAMPLE_IDS),
    )
    selected_samples_df = sample_table(config)
    adata_names, adata_list, _ = load_dataset(config)
    records = []
    adatas_by_run_sample_id = {}

    for run_sample_id, adata in zip(adata_names, adata_list):
        adatas_by_run_sample_id[run_sample_id] = adata
        row = selected_samples_df.loc[selected_samples_df["run_sample_id"] == run_sample_id].iloc[0]
        records.append(
            {
                "data_source": data_source,
                "matrix_source": str(row["matrix_source"]),
                "run_sample_id": run_sample_id,
                "cell_line": str(row["cell_line"]),
                "matrix_dir": str(row["matrix_dir"]),
                "matrix_exists": bool(row["matrix_exists"]),
                "n_obs": int(adata.n_obs),
                "n_vars": int(adata.n_vars),
                "x_nnz": nnz_of(adata),
                "obs_names_unique": bool(adata.obs_names.is_unique),
                "has_run_sample_id_obs": "run_sample_id" in adata.obs.columns,
                "has_cell_line_obs": "cell_line" in adata.obs.columns,
            }
        )

    return pd.DataFrame(records), adatas_by_run_sample_id


def main() -> int:
    args = parse_args()
    if not args.data_root:
        raise EnvironmentError("Set PROJECT_ROOT or MGE_DATA_ROOT, or pass --data-root.")

    print("Notebook 00 raw-vs-filtered source proof")
    print("host:", socket.gethostname())
    print("python:", sys.executable)
    print("repo_root:", Path(args.repo_root).resolve())
    print("data_root:", Path(args.data_root).resolve())
    print("sample_ids:", args.sample_ids or list(DEFAULT_SAMPLE_IDS))
    print("target_divs:", args.target_divs or ["DIV30"])

    raw_df, raw_adatas_by_run_sample_id = load_one_source(args, "cellranger_raw")
    filtered_df, filtered_adatas_by_run_sample_id = load_one_source(args, "cellranger_filtered")
    proof_df = pd.concat([raw_df, filtered_df], ignore_index=True)

    print("\nLoaded source summary:")
    print(proof_df.to_string(index=False))

    comparison_records = []
    for run_sample_id in sorted(set(raw_df["run_sample_id"]).intersection(filtered_df["run_sample_id"])):
        raw_row = raw_df.loc[raw_df["run_sample_id"] == run_sample_id].iloc[0]
        filtered_row = filtered_df.loc[filtered_df["run_sample_id"] == run_sample_id].iloc[0]
        raw_var_names = set(raw_adatas_by_run_sample_id[run_sample_id].var_names)
        filtered_var_names = set(filtered_adatas_by_run_sample_id[run_sample_id].var_names)
        shared_var_names = raw_var_names.intersection(filtered_var_names)
        comparison_records.append(
            {
                "run_sample_id": run_sample_id,
                "raw_n_obs": int(raw_row["n_obs"]),
                "filtered_n_obs": int(filtered_row["n_obs"]),
                "raw_minus_filtered_n_obs": int(raw_row["n_obs"]) - int(filtered_row["n_obs"]),
                "raw_n_vars": int(raw_row["n_vars"]),
                "filtered_n_vars": int(filtered_row["n_vars"]),
                "same_n_vars": int(raw_row["n_vars"]) == int(filtered_row["n_vars"]),
                "shared_n_vars": len(shared_var_names),
                "filtered_vars_missing_from_raw": len(filtered_var_names.difference(raw_var_names)),
                "different_matrix_dirs": str(raw_row["matrix_dir"]) != str(filtered_row["matrix_dir"]),
            }
        )

    comparison_df = pd.DataFrame(comparison_records)
    print("\nRaw vs filtered comparison:")
    print(comparison_df.to_string(index=False))

    failures = []
    for _, row in comparison_df.iterrows():
        if row["raw_n_obs"] < row["filtered_n_obs"]:
            failures.append(f"{row['run_sample_id']}: raw cells < filtered cells")
        if int(row["filtered_vars_missing_from_raw"]) > 0:
            failures.append(
                f"{row['run_sample_id']}: filtered has {row['filtered_vars_missing_from_raw']} genes not found in raw"
            )
        if not bool(row["different_matrix_dirs"]):
            failures.append(f"{row['run_sample_id']}: raw and filtered resolved to the same matrix_dir")

    for _, row in proof_df.iterrows():
        if not bool(row["obs_names_unique"]):
            failures.append(f"{row['data_source']} {row['run_sample_id']}: obs_names are not unique")
        if not bool(row["has_run_sample_id_obs"]):
            failures.append(f"{row['data_source']} {row['run_sample_id']}: missing obs run_sample_id")
        if not bool(row["has_cell_line_obs"]):
            failures.append(f"{row['data_source']} {row['run_sample_id']}: missing obs cell_line")

    if args.report_tsv:
        report_path = Path(args.report_tsv).expanduser().resolve()
        report_path.parent.mkdir(parents=True, exist_ok=True)
        merged = proof_df.copy()
        merged["proof_json"] = json.dumps(comparison_records, sort_keys=True)
        merged.to_csv(report_path, sep="\t", index=False)
        print("\nWrote report:", report_path)

    if failures:
        print("\nFAILURES:")
        for failure in failures:
            print(" -", failure)
        return 1

    print("\nPASS: raw and filtered sources both loaded and resolved as distinct Cell Ranger inputs.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
