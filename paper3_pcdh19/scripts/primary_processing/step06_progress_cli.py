"""Display the live Step 06 function-level progress event stream."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path


def parser() -> argparse.ArgumentParser:
    """Construct the run-directory, follow, and detail display options."""

    command = argparse.ArgumentParser(description=__doc__)
    command.add_argument("--run-dir", type=Path, required=True)
    command.add_argument("--follow", action="store_true")
    command.add_argument("--details", action="store_true")
    command.add_argument("--interval", type=float, default=2.0)
    return command


class Step06ProgressViewer:
    """Read and print new records from one append-only Step 06 event ledger."""

    def __init__(self, run_dir: Path, details: bool):
        """Resolve the run-scoped ledger and display mode."""

        self.path = run_dir / "provenance" / "step06_progress_events.jsonl"
        self.details = details

    def show(self, follow: bool, interval: float) -> None:
        """Print all existing events and optionally follow new events."""

        if interval <= 0:
            raise ValueError("--interval must be positive")
        while not self.path.exists():
            if not follow:
                raise FileNotFoundError(f"Progress ledger not available: {self.path}")
            time.sleep(interval)
        with self.path.open("r", encoding="utf-8") as handle:
            while True:
                line = handle.readline()
                if line:
                    record = json.loads(line)
                    self._print(record)
                    if follow and self._terminal(record):
                        return
                    continue
                if not follow:
                    return
                time.sleep(interval)

    @staticmethod
    def _terminal(record: dict) -> bool:
        """Return whether a workflow-level event declares terminal state."""

        return (
            record.get("function") == "Step06Workflow.run"
            and record.get("outputs", {}).get("status") in {"COMPLETE", "FAILED"}
        )

    def _print(self, record: dict) -> None:
        """Print one compact event plus optional exact input/output JSON."""

        elapsed = record["workflow_elapsed_seconds"]
        duration = record["function_elapsed_seconds"]
        print(
            f"{record['sequence']:04d}  {record['timestamp_utc']}  "
            f"+{elapsed:>10.3f}s  {record['status']:<8}  "
            f"[{record['stage']}]  {record['function']}  "
            f"duration={duration:.3f}s  max_rss={record['max_rss_gb']:.3f}GB",
            flush=True,
        )
        if self.details:
            print(
                "      inputs=" + json.dumps(record["inputs"], sort_keys=True),
                flush=True,
            )
            print(
                "      outputs=" + json.dumps(record["outputs"], sort_keys=True),
                flush=True,
            )


def main() -> None:
    """Display one run's durable progress ledger."""

    arguments = parser().parse_args()
    Step06ProgressViewer(arguments.run_dir, arguments.details).show(
        arguments.follow, arguments.interval
    )


if __name__ == "__main__":
    main()
