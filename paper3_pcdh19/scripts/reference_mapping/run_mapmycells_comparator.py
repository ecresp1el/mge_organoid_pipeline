#!/usr/bin/env python3
"""Run pinned local MapMyCells and publish every returned hierarchy metric."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd


class MapMyCellsContractError(RuntimeError):
    pass


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class MapMyCellsRunner:
    def __init__(self, args):
        self.args = args
        self.output_dir = args.output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.extended_json = self.output_dir / "mapmycells_extended.json"
        self.raw_csv = self.output_dir / "mapmycells_results.csv"
        self.log_path = self.output_dir / "mapmycells.log"

    def validate(self):
        for path in (self.args.query_h5ad, self.args.query_metadata,
                     self.args.marker_lookup, self.args.precomputed_stats):
            if not path.is_file():
                raise MapMyCellsContractError(f"Missing required input: {path}")
        if self.args.bootstrap_iteration < 2:
            raise MapMyCellsContractError("Bootstrapping must use at least two iterations")

    def command(self):
        return [
            sys.executable, "-m", "cell_type_mapper.cli.from_specified_markers",
            "--query_path", str(self.args.query_h5ad),
            "--extended_result_path", str(self.extended_json),
            "--csv_result_path", str(self.raw_csv),
            "--log_path", str(self.log_path),
            "--drop_level", self.args.drop_level,
            "--flatten", "False",
            "--max_gb", "100",
            "--cloud_safe", "False",
            "--verbose_csv", "True",
            "--query_markers.serialized_lookup", str(self.args.marker_lookup),
            "--query_markers.collapse_markers", "False",
            "--precomputed_stats.path", str(self.args.precomputed_stats),
            "--type_assignment.normalization", "raw",
            "--type_assignment.n_processors", str(self.args.cpus),
            "--type_assignment.bootstrap_iteration", str(self.args.bootstrap_iteration),
            "--type_assignment.bootstrap_factor", str(self.args.bootstrap_factor),
            "--type_assignment.rng_seed", str(self.args.seed),
            "--type_assignment.n_runners_up", "5",
            "--type_assignment.algorithm", "hierarchical",
        ]

    def execute(self):
        self.validate()
        command = self.command()
        (self.output_dir / "executed_command.txt").write_text(
            " ".join(subprocess.list2cmdline([item]) for item in command) + "\n",
            encoding="utf-8",
        )
        subprocess.run(command, check=True)
        if not self.extended_json.is_file() or not self.raw_csv.is_file():
            raise MapMyCellsContractError("MapMyCells did not produce both official outputs")

    @staticmethod
    def _scalar(value):
        if isinstance(value, (list, dict)):
            return json.dumps(value, separators=(",", ":"), sort_keys=True)
        return value

    def flatten(self):
        metadata = pd.read_csv(self.args.query_metadata, sep="\t", compression="gzip", dtype=str)
        if metadata["cell_id"].duplicated().any():
            raise MapMyCellsContractError("Query metadata cell IDs are not unique")
        csv_result = pd.read_csv(self.raw_csv, comment="#", dtype=str)
        if csv_result.columns[0] != "cell_id":
            csv_result = csv_result.rename(columns={csv_result.columns[0]: "cell_id"})
        csv_result = csv_result.rename(
            columns={name: "mmc_" + name for name in csv_result.columns if name != "cell_id"}
        )

        with self.extended_json.open("r", encoding="utf-8") as handle:
            extended = json.load(handle)
        results = extended.get("results")
        if not isinstance(results, list):
            raise MapMyCellsContractError("Extended JSON lacks a result list")
        flattened = []
        hierarchy_levels = set()
        for result in results:
            row = {"cell_id": str(result["cell_id"])}
            for level, payload in result.items():
                if level == "cell_id":
                    continue
                hierarchy_levels.add(level)
                if not isinstance(payload, dict):
                    row[f"mmc_extended__{level}"] = self._scalar(payload)
                    continue
                for field, value in payload.items():
                    row[f"mmc_extended__{level}__{field}"] = self._scalar(value)
            flattened.append(row)
        extended_frame = pd.DataFrame(flattened)

        expected = set(metadata["cell_id"])
        for frame, label in ((csv_result, "CSV"), (extended_frame, "extended JSON")):
            if frame["cell_id"].duplicated().any() or set(frame["cell_id"]) != expected:
                raise MapMyCellsContractError(f"MapMyCells {label} cell identity mismatch")
        combined = metadata.merge(csv_result, on="cell_id", validate="one_to_one")
        combined = combined.merge(extended_frame, on="cell_id", validate="one_to_one")
        output_path = self.output_dir / "mapmycells_per_cell.tsv.gz"
        combined.to_csv(output_path, sep="\t", index=False, compression="gzip", quoting=csv.QUOTE_MINIMAL)

        run_manifest = {
            "query_cells": len(metadata),
            "output_cells": len(combined),
            "hierarchy_machine_levels": sorted(hierarchy_levels),
            "drop_level": self.args.drop_level,
            "drop_level_reason": (
                "Official WMB level is not directly mapped during hierarchical traversal; "
                "MapMyCells retains it and backfills it from the lower assignment"
            ),
            "flatten": False,
            "normalization": "raw (MapMyCells internal log2(CPM+1))",
            "bootstrap_iteration": self.args.bootstrap_iteration,
            "bootstrap_factor": self.args.bootstrap_factor,
            "rng_seed": self.args.seed,
            "marker_lookup": str(self.args.marker_lookup),
            "marker_lookup_sha256": sha256_file(self.args.marker_lookup),
            "precomputed_stats": str(self.args.precomputed_stats),
            "precomputed_stats_sha256": sha256_file(self.args.precomputed_stats),
            "official_csv": str(self.raw_csv),
            "extended_json": str(self.extended_json) + ".gz",
            "flattened_tsv": str(output_path),
        }
        (self.output_dir / "mapmycells_contract.json").write_text(
            json.dumps(run_manifest, indent=2) + "\n", encoding="utf-8"
        )
        with self.extended_json.open("rb") as src, gzip.open(str(self.extended_json) + ".gz", "wb") as dst:
            shutil.copyfileobj(src, dst, length=8 * 1024 * 1024)
        self.extended_json.unlink()
        (self.output_dir / "SUCCESS.txt").write_text(
            "PASS\n"
            f"query_cells={len(metadata)}\n"
            f"output_cells={len(combined)}\n"
            "hierarchical=true\n"
            f"bootstrap_iteration={self.args.bootstrap_iteration}\n"
            "reclustered=false\nintegrated=false\nexisting_umap_recomputed=false\n",
            encoding="utf-8",
        )

    def run(self):
        self.execute()
        self.flatten()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--query-h5ad", required=True, type=Path)
    parser.add_argument("--query-metadata", required=True, type=Path)
    parser.add_argument("--marker-lookup", required=True, type=Path)
    parser.add_argument("--precomputed-stats", required=True, type=Path)
    parser.add_argument("--drop-level", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--cpus", required=True, type=int)
    parser.add_argument("--seed", required=True, type=int)
    parser.add_argument("--bootstrap-iteration", type=int, default=100)
    parser.add_argument("--bootstrap-factor", type=float, default=0.5)
    args = parser.parse_args()
    MapMyCellsRunner(args).run()


if __name__ == "__main__":
    main()
