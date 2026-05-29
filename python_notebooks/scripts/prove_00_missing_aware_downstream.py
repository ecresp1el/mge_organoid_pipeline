#!/usr/bin/env python
"""Prove downstream code can identify and plot missing source samples.

This proof intentionally does not load expression matrices. It exercises the
source-report layer that downstream QC, preprocessing, and plotting should use
before trying to operate on AnnData objects.
"""

from __future__ import annotations

import argparse
import os
import socket
import sys
from pathlib import Path

import matplotlib.pyplot as plt


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = REPO_ROOT / "python_notebooks" / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from mge_organoid_python.data_sources import Notebook00SourceConfig, load_dataset_result


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
        "--data-source",
        default="cellbender_denoised",
        choices=["cellranger_raw", "cellranger_filtered", "cellbender_denoised"],
        help="Source to check. Default: cellbender_denoised.",
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
    parser.add_argument("--report-tsv", required=True, help="TSV path for the source status report.")
    parser.add_argument("--summary-tsv", required=True, help="TSV path for the status count summary.")
    parser.add_argument("--plot-png", required=True, help="PNG path for the sample availability plot.")
    parser.add_argument(
        "--expect-missing-sample",
        action="append",
        dest="expected_missing_samples",
        help="Optional sample expected to be reported as missing. May be repeated.",
    )
    return parser.parse_args()


def write_availability_plot(source_table, output_png: Path) -> None:
    """Write a compact status plot for requested samples."""
    plot_df = source_table.copy()
    plot_df["is_available"] = plot_df["load_status"].eq("available")
    colors = plot_df["is_available"].map({True: "#2F855A", False: "#C53030"}).tolist()

    fig_width = max(7, 1.15 * len(plot_df))
    fig, ax = plt.subplots(figsize=(fig_width, 3.8))
    ax.bar(plot_df["run_sample_id"], [1] * len(plot_df), color=colors)
    ax.set_ylim(0, 1.25)
    ax.set_ylabel("source status")
    ax.set_xlabel("run_sample_id")
    ax.set_title(f"{plot_df['data_source'].iloc[0]} source availability")
    ax.set_yticks([])

    for idx, row in enumerate(plot_df.itertuples(index=False)):
        label = "available" if row.load_status == "available" else row.skip_reason
        ax.text(idx, 1.03, label, ha="center", va="bottom", rotation=35, fontsize=8)

    ax.spines[["left", "right", "top"]].set_visible(False)
    fig.tight_layout()
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=160)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    if not args.data_root:
        raise EnvironmentError("Set PROJECT_ROOT or MGE_DATA_ROOT, or pass --data-root.")

    config = Notebook00SourceConfig.from_defaults(
        data_source=args.data_source,
        repo_root=args.repo_root,
        data_root=args.data_root,
        target_divs=tuple(args.target_divs or ("DIV30",)),
        target_run_sample_ids=tuple(args.sample_ids or DEFAULT_SAMPLE_IDS),
        strict_missing_matrix_dirs=False,
    )

    result = load_dataset_result(config, load_matrices=False)
    source_table = result.source_table
    summary_df = result.availability_summary()

    print("Notebook 00 missing-aware downstream proof")
    print("host:", socket.gethostname())
    print("python:", sys.executable)
    print("repo_root:", config.repo_root)
    print("data_root:", config.data_root)
    print("data_source:", result.data_source)
    print("available_samples:", result.available_samples)
    print("skipped_samples:", result.skipped_samples)

    print("\nSource table:")
    print(
        source_table[
            ["run_sample_id", "data_source", "source_exists", "load_status", "skip_reason", "source_path"]
        ].to_string(index=False)
    )

    print("\nAvailability summary:")
    print(summary_df.to_string(index=False))

    report_tsv = Path(args.report_tsv).expanduser().resolve()
    summary_tsv = Path(args.summary_tsv).expanduser().resolve()
    plot_png = Path(args.plot_png).expanduser().resolve()
    report_tsv.parent.mkdir(parents=True, exist_ok=True)
    summary_tsv.parent.mkdir(parents=True, exist_ok=True)

    source_table.to_csv(report_tsv, sep="\t", index=False)
    summary_df.to_csv(summary_tsv, sep="\t", index=False)
    write_availability_plot(source_table, plot_png)

    print("\nWrote report:", report_tsv)
    print("Wrote summary:", summary_tsv)
    print("Wrote plot:", plot_png)

    failures = []
    expected_missing_samples = set(args.expected_missing_samples or [])
    actual_missing_samples = set(result.skipped_samples)
    missing_not_reported = sorted(expected_missing_samples.difference(actual_missing_samples))
    if missing_not_reported:
        failures.append(f"Expected missing samples were not reported as missing: {missing_not_reported}")

    if not plot_png.exists() or plot_png.stat().st_size == 0:
        failures.append(f"Availability plot was not written: {plot_png}")

    if failures:
        print("\nFAILURES:")
        for failure in failures:
            print(" -", failure)
        return 1

    print("\nPASS: downstream source reporting identifies available and skipped samples.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
