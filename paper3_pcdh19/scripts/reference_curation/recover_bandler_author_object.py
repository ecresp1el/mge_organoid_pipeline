#!/usr/bin/env python3
"""Recover an author-linked Bandler Seurat object with immutable provenance."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import platform
import sys
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping


def sha256_file(path: Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_tsv(path: Path, columns: Iterable[str], rows: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


@dataclass(frozen=True)
class AuthorArtifact:
    name: str
    url: str
    expected_bytes: int
    repository: str
    discovery_commit: str
    discovery_path: str
    current_status: str


class ImmutableHttpCache:
    """Materialize a remote artifact atomically, resuming a partial transfer."""

    def __init__(self, artifact: AuthorArtifact, destination: Path):
        self.artifact = artifact
        self.destination = destination.resolve()
        self.partial = self.destination.with_suffix(self.destination.suffix + ".partial")

    def _request(self, start: int = 0):
        headers = {"User-Agent": "mge-organoid-reference-curation/1.0"}
        if start:
            headers["Range"] = f"bytes={start}-"
        return urllib.request.urlopen(urllib.request.Request(self.artifact.url, headers=headers), timeout=120)

    def materialize(self) -> tuple[str, str]:
        self.destination.parent.mkdir(parents=True, exist_ok=True)
        if self.destination.exists():
            if self.destination.stat().st_size != self.artifact.expected_bytes:
                raise RuntimeError(f"Existing cache has unexpected size: {self.destination}")
            return "reused", sha256_file(self.destination)

        start = self.partial.stat().st_size if self.partial.exists() else 0
        response = self._request(start)
        status = getattr(response, "status", response.getcode())
        if start and status != 206:
            response.close()
            start = 0
            response = self._request(0)
        mode = "ab" if start and status == 206 else "wb"
        with response, self.partial.open(mode) as output:
            while True:
                block = response.read(8 * 1024 * 1024)
                if not block:
                    break
                output.write(block)
            output.flush()
            os.fsync(output.fileno())
        observed_bytes = self.partial.stat().st_size
        if observed_bytes != self.artifact.expected_bytes:
            raise RuntimeError(
                f"Incomplete author object: expected {self.artifact.expected_bytes}, observed {observed_bytes}"
            )
        digest = sha256_file(self.partial)
        os.replace(self.partial, self.destination)
        return "downloaded", digest


class BandlerAuthorObjectRecovery:
    """Publish discovery evidence and immutable-download provenance."""

    ARTIFACT = AuthorArtifact(
        name="STICR.seuratobject.RDS",
        url="https://keeper.mpdl.mpg.de/f/0d1a23ab36fe473481b0/?dl=1",
        expected_bytes=1_023_351_516,
        repository="https://github.com/mayer-lab/Bandler-et-al_lineage",
        discovery_commit="a8fa139a5ed6d8832b07b61f384982d630893c93",
        discovery_path="README.md",
        current_status="Historical author link remains live; current README replaced it with <seurat_obj_url>.",
    )

    def __init__(self, run_dir: Path, destination: Path):
        self.run_dir = run_dir.resolve()
        self.destination = destination.resolve()

    def run(self) -> None:
        cache = ImmutableHttpCache(self.ARTIFACT, self.destination)
        action, digest = cache.materialize()
        now = datetime.now(timezone.utc).isoformat(timespec="seconds")
        audit = self.run_dir / "Bandler2022" / "author_object_audit"
        atomic_tsv(audit / "author_object_discovery.tsv", (
            "artifact", "author_url", "repository", "discovery_commit", "discovery_path", "current_status",
        ), ({
            "artifact": self.ARTIFACT.name,
            "author_url": self.ARTIFACT.url,
            "repository": self.ARTIFACT.repository,
            "discovery_commit": self.ARTIFACT.discovery_commit,
            "discovery_path": self.ARTIFACT.discovery_path,
            "current_status": self.ARTIFACT.current_status,
        },))
        atomic_tsv(audit / "author_object_download.tsv", (
            "artifact", "resolved_path", "bytes", "sha256", "action", "completed_utc",
        ), ({
            "artifact": self.ARTIFACT.name,
            "resolved_path": self.destination,
            "bytes": self.destination.stat().st_size,
            "sha256": digest,
            "action": action,
            "completed_utc": now,
        },))
        atomic_tsv(audit / "recovery_environment.tsv", (
            "completed_utc", "python", "executed_script", "script_sha256", "command",
        ), ({
            "completed_utc": now,
            "python": platform.python_version(),
            "executed_script": str(Path(__file__).resolve()),
            "script_sha256": sha256_file(Path(__file__).resolve()),
            "command": " ".join(sys.argv),
        },))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--destination", required=True, type=Path)
    args = parser.parse_args()
    BandlerAuthorObjectRecovery(args.run_dir, args.destination).run()


if __name__ == "__main__":
    main()
