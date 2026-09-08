"""Inspect a Step 06 JSONL progress stream without running analysis.

``Step06ProgressViewer`` prints old events and optionally follows new ones.
--details includes exact logged inputs/outputs; --interval sets polling seconds.
Follow mode stops on workflow COMPLETE/FAILED, not an individual function's
completion. It does not query SLURM or infer failure if a job is killed before
writing a terminal event. An older run may have no event stream at all.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path


def parser() -> argparse.ArgumentParser:
    """Construct the run-directory, follow, and detail display options.

    Returns
    -------
    argparse.ArgumentParser
        Parser for run directory, follow flag, details flag and polling interval.

    Notes
    -----
    Default interval is 2 seconds; display choices do not change analysis or scheduler
    state.
    """

    command = argparse.ArgumentParser(description=__doc__)
    command.add_argument("--run-dir", type=Path, required=True)
    command.add_argument("--follow", action="store_true")
    command.add_argument("--details", action="store_true")
    command.add_argument("--interval", type=float, default=2.0)
    return command


class Step06ProgressViewer:
    """Read and print new records from one append-only Step 06 event ledger.

    Notes
    -----
    Own read-only display of a run-scoped event ledger. show() can follow new events
    but does not query SLURM, cancel jobs or infer termination. The details switch
    controls whether exact recorded input/output metadata are printed.
    """

    def __init__(self, run_dir: Path, details: bool):
        """Resolve the run-scoped ledger and display mode.

        Parameters
        ----------
        run_dir : pathlib.Path
            Run-scoped directory; never the approved upstream checkpoint directory.
        details : bool
            Whether to print input/output JSON after each compact event line.

        Notes
        -----
        Resolve provenance/step06_progress_events.jsonl and store display mode. Does not
        create the file or run analysis.
        """

        self.path = run_dir / "provenance" / "step06_progress_events.jsonl"
        self.details = details

    def show(self, follow: bool, interval: float) -> None:
        """Print all existing events and optionally follow new events.

        Parameters
        ----------
        follow : bool
            Whether to wait for new events after reaching the current end of the file.
        interval : float
            Positive polling interval in seconds; affects display latency only.

        Returns
        -------
        None
            Effects are described below.

        Raises
        ------
        ValueError or FileNotFoundError
            Interval is nonpositive, or the ledger is missing without follow.

        Notes
        -----
        Print all existing records in order, then optionally poll for new records.
        Follow mode waits for a missing ledger and stops on a workflow terminal event.
        Without follow, missing files raise immediately and EOF ends reading. A killed
        job without a terminal event will keep follow mode waiting. Older job 59983804
        has no live ledger. Interrupt the viewer with Ctrl-C; it does not cancel jobs.
        """

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
        """Return whether a workflow-level event declares terminal state.

        Parameters
        ----------
        record : dict
            One decoded progress-event record.

        Returns
        -------
        bool
            True only for Step06Workflow.run with outputs.status COMPLETE or FAILED.

        Notes
        -----
        A single operation's COMPLETE event is not workflow completion.
        """

        return (
            record.get("function") == "Step06Workflow.run"
            and record.get("outputs", {}).get("status") in {"COMPLETE", "FAILED"}
        )

    def _print(self, record: dict) -> None:
        """Print one compact event plus optional exact input/output JSON.

        Parameters
        ----------
        record : dict
            One decoded progress-event record.

        Returns
        -------
        None
            Effects are described below.

        Notes
        -----
        Write sequence/time/status/function/duration/peak RSS to stdout, with optional sorted
        input/output JSON. No event files are changed.
        """

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
    """Display one run's durable progress ledger.

    Returns
    -------
    None
        Effects are described below.

    Notes
    -----
    Parse process arguments and display the chosen run's ledger. No scientific workflow is
    constructed.
    """

    arguments = parser().parse_args()
    Step06ProgressViewer(arguments.run_dir, arguments.details).show(
        arguments.follow, arguments.interval
    )


if __name__ == "__main__":
    main()
