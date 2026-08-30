"""Command-line entry point for the frozen PCDH19 primary-processing package."""

from __future__ import annotations

import argparse
from pathlib import Path

from .models import Step00Paths, Step00Settings
from .workflow import Step00Workflow


def _parser() -> argparse.ArgumentParser:
    """Construct the explicit Step 00 command-line interface."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cellranger-root", type=Path, required=True)
    parser.add_argument("--sample-key", type=Path, required=True)
    parser.add_argument("--technical-manifest", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--workflow-root", type=Path, required=True)
    parser.add_argument("--frozen-code-dir", type=Path, required=True)
    parser.add_argument("--expected-samples", type=int, required=True)
    parser.add_argument("--expected-cells", type=int, required=True)
    parser.add_argument("--expected-genes", type=int, required=True)
    parser.add_argument("--expected-genome", required=True)
    parser.add_argument("--feature-type", required=True)
    parser.add_argument("--compression", choices=("gzip", "lzf"), required=True)
    return parser


def main() -> int:
    """Run only Step 00 and print the published checkpoint path."""

    args = _parser().parse_args()
    paths = Step00Paths(
        cellranger_root=args.cellranger_root.resolve(),
        sample_key=args.sample_key.resolve(),
        technical_manifest=args.technical_manifest.resolve(),
        run_dir=args.run_dir.resolve(),
        workflow_root=args.workflow_root.resolve(),
        frozen_code_dir=args.frozen_code_dir.resolve(),
    )
    settings = Step00Settings(
        expected_samples=args.expected_samples,
        expected_cells=args.expected_cells,
        expected_genes=args.expected_genes,
        expected_genome=args.expected_genome,
        feature_type=args.feature_type,
        compression=args.compression,
    )
    checkpoint = Step00Workflow(paths, settings).run()
    print(f"Published Step 00 checkpoint: {checkpoint}")
    print("Status: IN_REVIEW")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
