"""Record durable function-level events without changing scientific results.

``Step06ProgressTracker.track`` writes START then COMPLETE/FAILED around a
block; callers add measured outputs to the yielded context. ``note`` writes
instantaneous events. JSONL is flushed/fsynced and the latest-event JSON is
replaced atomically. ``NullStep06ProgressTracker`` provides a no-I/O adapter.

Timing uses a monotonic clock. max_rss_gb is process-lifetime peak resident
memory (Linux ru_maxrss in KiB converted to GiB), not per-function allocation.
This is a single-writer run log. The September 3 job 59983804 predates these
live events; its frozen executable and original logs remain authoritative.
"""

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
    """Return a timezone-aware UTC timestamp for one progress event.

    Returns
    -------
    str
        Timezone-aware UTC ISO-8601 timestamp.

    Notes
    -----
    Wall-clock event timestamp; elapsed durations use a separate monotonic clock.
    """

    return datetime.now(timezone.utc).isoformat()


def _json_default(value):
    """Convert paths and scalar-like scientific values to JSON-safe values.

    Parameters
    ----------
    value : object
        Value requiring fallback JSON conversion.

    Returns
    -------
    object
        Path string, scalar .item() result, or string fallback.

    Notes
    -----
    Used only by JSON encoding. Does not serialize a full AnnData; callers should log
    compact metadata.
    """

    if isinstance(value, Path):
        return str(value)
    if hasattr(value, "item"):
        return value.item()
    return str(value)


@dataclass
class Step06ProgressContext:
    """Collect outputs that become part of a completed function event.

    Notes
    -----
    Each instance owns an independent outputs dictionary. A tracked block adds result
    metadata before its COMPLETE or FAILED event is encoded.
    """

    outputs: dict = field(default_factory=dict)


class NullStep06ProgressTracker:
    """Provide the progress API for isolated unit-level use without file output.

    Notes
    -----
    Provide the same track()/note() interface as the file-backed tracker for
    component-level use. Context outputs exist in memory only; no timing or event
    files are published. The production workflow supplies the real tracker.
    """

    @contextmanager
    def track(
        self, stage: str, function: str, inputs: dict
    ) -> Iterator[Step06ProgressContext]:
        """Yield an in-memory context without publishing progress events.

        Parameters
        ----------
        stage : str
            Readable workflow phase identifier.
        function : str
            Function or logical-block label to record; not resolved or invoked
            dynamically.
        inputs : dict
            Explicit provenance metadata for the operation, not the expression matrix
            itself.

        Yields
        ------
        Step06ProgressContext
            Mutable output context yielded by a context manager.

        Notes
        -----
        Accept the same arguments as the real tracker but write no events. Exceptions from the
        caller's block propagate unchanged.
        """

        del stage, function, inputs
        yield Step06ProgressContext()

    def note(self, stage: str, function: str, inputs: dict, outputs: dict) -> None:
        """Accept an instantaneous event without publishing it.

        Parameters
        ----------
        stage : str
            Readable workflow phase identifier.
        function : str
            Function or logical-block label to record; not resolved or invoked
            dynamically.
        inputs : dict
            Explicit provenance metadata for the operation, not the expression matrix
            itself.
        outputs : dict
            Result metadata for the operation; serialized into the event.

        Returns
        -------
        None
            Effects are described below.

        Notes
        -----
        Discard an informational event without file output; used for isolated component
        execution.
        """

        del stage, function, inputs, outputs


class Step06ProgressTracker:
    """Write durable JSONL events and an atomically refreshed latest-event file.

    Notes
    -----
    Own the single-writer JSONL stream and latest-event snapshot for one new
    execution. track() brackets operations and note() records instantaneous status.
    Existing logs are refused to keep runs distinct; errors propagate rather than
    silently losing provenance.
    """

    def __init__(self, provenance_dir: Path):
        """Initialize a new run-scoped append-only event stream.

        Parameters
        ----------
        provenance_dir : pathlib.Path
            Directory for the new JSONL event stream and latest-event JSON.

        Notes
        -----
        Create provenance_dir if needed; resolve JSONL/latest paths and start the elapsed
        timer. A preexisting events file raises FileExistsError to prevent mixing executions.
        """

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
        """Record START and COMPLETE/FAILED events around one major function.

        Parameters
        ----------
        stage : str
            Readable workflow phase identifier.
        function : str
            Function or logical-block label to record; not resolved or invoked
            dynamically.
        inputs : dict
            Explicit provenance metadata for the operation, not the expression matrix
            itself.

        Yields
        ------
        Step06ProgressContext
            Context-manager yield; add result metadata to context.outputs inside the with
            block.

        Notes
        -----
        Write START before yielding; write COMPLETE afterward, or FAILED with error
        type/message before re-raising a caught exception. Durations bracket the entire block.
        This is instrumentation, not an alternative implementation of the labeled function.
        """

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
        """Record one instantaneous informational event.

        Parameters
        ----------
        stage : str
            Readable workflow phase identifier.
        function : str
            Function or logical-block label to record; not resolved or invoked
            dynamically.
        inputs : dict
            Explicit provenance metadata for the operation, not the expression matrix
            itself.
        outputs : dict
            Result metadata for the operation; serialized into the event.

        Returns
        -------
        None
            Effects are described below.

        Notes
        -----
        Immediately write a NOTE event with duration zero; workflow termination is carried in
        outputs.status, not in the outer NOTE status.
        """

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
        """Append, flush, and fsync one event before refreshing latest state.

        Parameters
        ----------
        stage : str
            Readable workflow phase identifier.
        status : str
            Event state, normally START, COMPLETE, FAILED, or NOTE.
        function : str
            Function or logical-block label to record; not resolved or invoked
            dynamically.
        inputs : dict
            Explicit provenance metadata for the operation, not the expression matrix
            itself.
        outputs : dict
            Result metadata for the operation; serialized into the event.
        function_elapsed_seconds : float
            Measured wall-clock duration for this operation, in seconds.

        Returns
        -------
        None
            Effects are described below.

        Notes
        -----
        Increment sequence, append/flush/fsync JSONL, then replace latest JSON using
        a temporary sibling file. Each record includes PID, UTC time, elapsed time,
        and process-lifetime peak RSS. Linux ru_maxrss is divided by 1024**2, yielding
        GiB despite the legacy max_rss_gb field name. Single-writer use is assumed;
        there is no cross-process lock. Log I/O failures propagate.
        """

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
