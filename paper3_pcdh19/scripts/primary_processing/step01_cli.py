"""Command-line entry point for frozen primary-processing Step 01."""

from __future__ import annotations

import argparse
from pathlib import Path

from .step01_models import Step01Paths, Step01Settings
from .step01_workflow import Step01Workflow


def _parse_percent_top(value: str) -> tuple[int, ...] | None:
    """Parse `none` or a comma-delimited positive integer sequence."""

    if value.lower() == "none":
        return None
    parsed = tuple(int(item) for item in value.split(",") if item)
    if not parsed or any(item <= 0 for item in parsed):
        raise argparse.ArgumentTypeError("percent-top must be none or positive comma-delimited integers")
    return parsed


def _parser() -> argparse.ArgumentParser:
    """Construct the explicit no-filter Step 01 command-line interface."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", type=Path, required=True)
    parser.add_argument("--input-step-status", type=Path, required=True)
    parser.add_argument("--input-output-manifest", type=Path, required=True)
    parser.add_argument("--approval-ledger-snapshot", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--workflow-root", type=Path, required=True)
    parser.add_argument("--frozen-code-dir", type=Path, required=True)
    parser.add_argument("--expected-input-run-id", required=True)
    parser.add_argument("--expected-input-bytes", type=int, required=True)
    parser.add_argument("--expected-input-sha256", required=True)
    parser.add_argument("--expected-cells", type=int, required=True)
    parser.add_argument("--expected-genes", type=int, required=True)
    parser.add_argument("--mitochondrial-prefix", required=True)
    parser.add_argument("--ribosomal-prefixes", required=True)
    parser.add_argument("--percent-top", type=_parse_percent_top, required=True)
    parser.add_argument("--histogram-bins", type=int, required=True)
    parser.add_argument("--plot-dpi", type=int, required=True)
    parser.add_argument("--compression", choices=("gzip", "lzf"), required=True)
    return parser


def main() -> int:
    """Run Step 01 only and print the IN_REVIEW checkpoint location."""

    args = _parser().parse_args()
    paths = Step01Paths(
        input_h5ad=args.input_h5ad.resolve(),
        input_step_status=args.input_step_status.resolve(),
        input_output_manifest=args.input_output_manifest.resolve(),
        approval_ledger_snapshot=args.approval_ledger_snapshot.resolve(),
        run_dir=args.run_dir.resolve(),
        workflow_root=args.workflow_root.resolve(),
        frozen_code_dir=args.frozen_code_dir.resolve(),
        expected_input_run_id=args.expected_input_run_id,
        expected_input_bytes=args.expected_input_bytes,
        expected_input_sha256=args.expected_input_sha256,
    )
    settings = Step01Settings(
        expected_cells=args.expected_cells,
        expected_genes=args.expected_genes,
        mitochondrial_prefix=args.mitochondrial_prefix,
        ribosomal_prefixes=tuple(item for item in args.ribosomal_prefixes.split(",") if item),
        percent_top=args.percent_top,
        histogram_bins=args.histogram_bins,
        plot_dpi=args.plot_dpi,
        compression=args.compression,
    )
    checkpoint = Step01Workflow(paths, settings).run()
    print(f"Published Step 01 checkpoint: {checkpoint}")
    print("Status: IN_REVIEW")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
