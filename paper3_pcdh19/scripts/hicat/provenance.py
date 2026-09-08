"""Durable function events, content manifests, and JSON-safe metadata for HiCAT.

Each major operation prints START/COMPLETE/FAILED and writes the same event to
JSONL. This module manages provenance only; it contains no clustering rules.
"""
from contextlib import contextmanager
from collections.abc import Mapping
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import os
import resource
import time
import numpy as np
import pandas as pd


def json_value(value):
    """Convert Paths/NumPy scalars/arrays for provenance JSON, not expression data."""
    # AnnData 0.8 exposes uns as a Mapping wrapper rather than a native dict.
    # Preserve its nested structure so the inventory remains machine-readable.
    if isinstance(value, Mapping):
        return dict(value)
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return str(value)


def write_json(path, value):
    """Atomically write a JSON metadata file at path; replace only that file."""
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, default=json_value) + '\n')
    os.replace(temporary, path)


def sha256(path):
    """Return the streaming SHA-256 of a file without changing its contents."""
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def manifest(root):
    """Return relative paths, byte sizes and hashes of files below root.

    The manifest itself is excluded to avoid self-reference. The caller writes
    this table only after scientific outputs and their status are complete.
    """
    root = Path(root)
    return pd.DataFrame([
        {'relative_path': str(p.relative_to(root)), 'bytes': p.stat().st_size,
         'sha256': sha256(p)}
        for p in sorted(root.rglob('*')) if p.is_file() and p.name != 'output_manifest.tsv'
    ])


class Progress:
    """Record one execution's operations and their resolved inputs/results."""

    def __init__(self, directory):
        """Create a new single-writer stream; refuse to mix with an older stream."""
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.path = self.directory / 'hicat_progress_events.jsonl'
        if self.path.exists():
            raise FileExistsError(str(self.path))
        self.start = time.monotonic()
        self.sequence = 0

    def note(self, stage, status, inputs, outputs, elapsed=0):
        """Print and fsync one metadata event; RSS is lifetime peak GiB on Linux."""
        self.sequence += 1
        row = dict(sequence=self.sequence, utc=datetime.now(timezone.utc).isoformat(),
                   stage=stage, status=status, inputs=inputs, outputs=outputs,
                   operation_seconds=elapsed, total_seconds=time.monotonic()-self.start,
                   peak_rss_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2)
        with self.path.open('a') as handle:
            handle.write(json.dumps(row, default=json_value) + '\n')
            handle.flush()
            os.fsync(handle.fileno())
        write_json(self.directory / 'hicat_progress_latest.json', row)
        print(json.dumps(row, default=json_value), flush=True)

    @contextmanager
    def track(self, stage, **inputs):
        """Yield an output-metadata dict, then log completion or re-raise failure.

        Parameters
        ----------
        stage : str
            Function/phase name, including replicate and node when relevant.
        **inputs : dict
            Resolved settings, dimensions and input representation/path.

        Yields
        ------
        dict
            Populate with observed dimensions, counts and saved file names.
        """
        start = time.monotonic()
        self.note(stage, 'START', inputs, {})
        outputs = {}
        try:
            yield outputs
        except Exception as error:
            outputs.update(error_type=type(error).__name__, error=str(error))
            self.note(stage, 'FAILED', inputs, outputs, time.monotonic()-start)
            raise
        self.note(stage, 'COMPLETE', inputs, outputs, time.monotonic()-start)
