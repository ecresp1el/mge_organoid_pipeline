"""Append-only live progress and exact-function provenance for Step 06."""

from __future__ import annotations

import json
import os
import resource
import time
from contextlib import contextmanager
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterator


def _utc_now() -> str:
    """Return a timezone-aware UTC timestamp for one progress event."""

    return datetime.now(timezone.utc).isoformat()


def _json_default(value):
    """Convert paths and scalar-like scientific values to JSON-safe values."""

    if isinstance(value, Path):
        return str(value)
    if hasattr(value, "item"):
        return value.item()
    return str(value)


@dataclass
class Step06ProgressContext:
    """Collect outputs that become part of a completed function event."""

    outputs: dict = field(default_factory=dict)


class NullStep06ProgressTracker:
    """Provide the progress API for isolated unit-level use without file output."""

    @contextmanager
    def track(
        self, stage: str, function: str, inputs: dict
    ) -> Iterator[Step06ProgressContext]:
        """Yield an in-memory context without publishing progress events."""

        del stage, function, inputs
        yield Step06ProgressContext()

    def note(self, stage: str, function: str, inputs: dict, outputs: dict) -> None:
        """Accept an instantaneous event without publishing it."""

        del stage, function, inputs, outputs


class Step06ProgressTracker:
    """Write durable JSONL events and an atomically refreshed latest-event file."""

    def __init__(self, provenance_dir: Path):
        """Initialize a new run-scoped append-only event stream."""

        self.provenance_dir = provenance_dir
        self.provenance_dir.mkdir(parents=True, exist_ok=True)
        self.events_path = provenance_dir / "step06_progress_events.jsonl"
        self.latest_path = provenance_dir / "step06_progress_latest.json"
        if self.events_path.exists():
            raise FileExistsError(
                f"Refusing to append a new execution to existing progress ledger: {self.events_path}"
            )
        self.sequence = 0
        self.workflow_started = time.monotonic()

    @contextmanager
    def track(
        self, stage: str, function: str, inputs: dict
    ) -> Iterator[Step06ProgressContext]:
        """Record START and COMPLETE/FAILED events around one major function."""

        started = time.monotonic()
        self._write(stage, "START", function, inputs, {}, 0.0)
        context = Step06ProgressContext()
        try:
            yield context
        except Exception as error:
            self._write(
                stage,
                "FAILED",
                function,
                inputs,
                {
                    **context.outputs,
                    "error_type": type(error).__name__,
                    "error_message": str(error),
                },
                time.monotonic() - started,
            )
            raise
        self._write(
            stage,
            "COMPLETE",
            function,
            inputs,
            context.outputs,
            time.monotonic() - started,
        )

    def note(self, stage: str, function: str, inputs: dict, outputs: dict) -> None:
        """Record one instantaneous informational event."""

        self._write(stage, "NOTE", function, inputs, outputs, 0.0)

    def _write(
        self,
        stage: str,
        status: str,
        function: str,
        inputs: dict,
        outputs: dict,
        function_elapsed_seconds: float,
    ) -> None:
        """Append, flush, and fsync one event before refreshing latest state."""

        self.sequence += 1
        record = {
            "sequence": self.sequence,
            "timestamp_utc": _utc_now(),
            "workflow_elapsed_seconds": round(time.monotonic() - self.workflow_started, 3),
            "function_elapsed_seconds": round(function_elapsed_seconds, 3),
            "pid": os.getpid(),
            "max_rss_gb": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**2, 3),
            "stage": stage,
            "status": status,
            "function": function,
            "inputs": inputs,
            "outputs": outputs,
        }
        encoded = json.dumps(record, sort_keys=True, default=_json_default)
        with self.events_path.open("a", encoding="utf-8") as handle:
            handle.write(encoded + "\n")
            handle.flush()
            os.fsync(handle.fileno())
        temporary = self.latest_path.with_suffix(".json.tmp")
        temporary.write_text(encoded + "\n", encoding="utf-8")
        os.replace(temporary, self.latest_path)
