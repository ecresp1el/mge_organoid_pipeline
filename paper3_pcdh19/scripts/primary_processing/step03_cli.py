"""Command-line preparation and finalization for Step 03 scDblFinder."""

from __future__ import annotations

import argparse
from pathlib import Path

from .step03_io import Step03PreparationWorkflow
from .step03_models import Step03Paths, Step03Settings
from .step03_workflow import Step03FinalizationWorkflow


def parser() -> argparse.ArgumentParser:
    """Construct the explicit approved-input and one-capture command contract."""

    command = argparse.ArgumentParser(description=__doc__)
    command.add_argument("phase", choices=("prepare", "finalize"))
    command.add_argument("--input-h5ad", type=Path, required=True)
    command.add_argument("--input-status", type=Path, required=True)
    command.add_argument("--input-manifest", type=Path, required=True)
    command.add_argument("--approval-ledger-snapshot", type=Path, required=True)
    command.add_argument("--run-dir", type=Path, required=True)
    command.add_argument("--workflow-root", type=Path, required=True)
    command.add_argument("--frozen-code-dir", type=Path, required=True)
    command.add_argument("--expected-step02-run-id", required=True)
    command.add_argument("--expected-input-bytes", type=int, required=True)
    command.add_argument("--expected-input-sha256", required=True)
    command.add_argument("--expected-cells", type=int, default=446_349)
    command.add_argument("--expected-genes", type=int, default=19_071)
    command.add_argument("--capture-id", default="GEX_1")
    command.add_argument("--sample-field", default="technical_sample_id")
    command.add_argument("--design-field", default="design_group")
    command.add_argument("--primary-seed", type=int, default=20_260_830)
    command.add_argument("--reproducibility-seed", type=int, default=20_260_831)
    command.add_argument("--plot-dpi", type=int, default=180)
    command.add_argument("--compression", default="lzf")
    return command


def main() -> None:
    """Prepare the R bridge or finalize its non-filtering review checkpoint."""

    arguments = parser().parse_args()
    settings = Step03Settings(
        expected_cells=arguments.expected_cells,
        expected_genes=arguments.expected_genes,
        capture_id=arguments.capture_id,
        sample_field=arguments.sample_field,
        design_field=arguments.design_field,
        primary_seed=arguments.primary_seed,
        reproducibility_seed=arguments.reproducibility_seed,
        plot_dpi=arguments.plot_dpi,
        compression=arguments.compression,
    )
    paths = Step03Paths(
        input_h5ad=arguments.input_h5ad,
        input_status=arguments.input_status,
        input_manifest=arguments.input_manifest,
        approval_ledger_snapshot=arguments.approval_ledger_snapshot,
        run_dir=arguments.run_dir,
        workflow_root=arguments.workflow_root,
        frozen_code_dir=arguments.frozen_code_dir,
        expected_step02_run_id=arguments.expected_step02_run_id,
        expected_input_bytes=arguments.expected_input_bytes,
        expected_input_sha256=arguments.expected_input_sha256,
    )
    if arguments.phase == "prepare":
        output = Step03PreparationWorkflow(paths, settings).run()
        print(f"Prepared exact one-capture native-R bridge: {output}")
    else:
        output = Step03FinalizationWorkflow(paths, settings).run()
        print(f"Published Step 03 detection checkpoint: {output}")
        print("Status: IN_REVIEW; no predicted doublets removed")


if __name__ == "__main__":
    main()
