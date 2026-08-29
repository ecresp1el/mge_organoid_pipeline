#!/usr/bin/env python3
"""Recover public Bandler author objects and supplements with provenance."""

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
    relative_path: str
    artifact_role: str
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
    """Publish discovery evidence and cache each public artifact atomically."""

    ARTIFACTS = (
        AuthorArtifact(
            name="STICR.seuratobject.RDS",
            url="https://keeper.mpdl.mpg.de/f/0d1a23ab36fe473481b0/?dl=1",
            expected_bytes=1_023_351_516,
            relative_path="author_objects/STICR.seuratobject.RDS",
            artifact_role="Author Seurat object for the postnatal STICR reference.",
            repository="https://github.com/mayer-lab/Bandler-et-al_lineage",
            discovery_commit="a8fa139a5ed6d8832b07b61f384982d630893c93",
            discovery_path="README.md",
            current_status="Historical author link remains live; current README replaced it with <seurat_obj_url>.",
        ),
        AuthorArtifact(
            name="41586_2021_4237_MOESM2_ESM.xlsx",
            url="https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-021-04237-0/MediaObjects/41586_2021_4237_MOESM2_ESM.xlsx",
            expected_bytes=22_008,
            relative_path="metadata_sources/41586_2021_4237_MOESM2_ESM.xlsx",
            artifact_role="Supplementary Data 1: sample design and filtering ledger.",
            repository="https://www.nature.com/articles/s41586-021-04237-0",
            discovery_commit="NA",
            discovery_path="Supplementary Data 1",
            current_status="Public publisher supplement.",
        ),
        AuthorArtifact(
            name="41586_2021_4237_MOESM3_ESM.xlsx",
            url="https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-021-04237-0/MediaObjects/41586_2021_4237_MOESM3_ESM.xlsx",
            expected_bytes=21_586,
            relative_path="metadata_sources/41586_2021_4237_MOESM3_ESM.xlsx",
            artifact_role="Supplementary Data 2: postnatal broad-class markers.",
            repository="https://www.nature.com/articles/s41586-021-04237-0",
            discovery_commit="NA",
            discovery_path="Supplementary Data 2",
            current_status="Public publisher supplement.",
        ),
        AuthorArtifact(
            name="41586_2021_4237_MOESM4_ESM.xlsx",
            url="https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-021-04237-0/MediaObjects/41586_2021_4237_MOESM4_ESM.xlsx",
            expected_bytes=70_339,
            relative_path="metadata_sources/41586_2021_4237_MOESM4_ESM.xlsx",
            artifact_role="Supplementary Data 3: postnatal refined-cluster markers.",
            repository="https://www.nature.com/articles/s41586-021-04237-0",
            discovery_commit="NA",
            discovery_path="Supplementary Data 3",
            current_status="Public publisher supplement.",
        ),
        AuthorArtifact(
            name="41586_2021_4237_MOESM5_ESM.xlsx",
            url="https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-021-04237-0/MediaObjects/41586_2021_4237_MOESM5_ESM.xlsx",
            expected_bytes=46_056,
            relative_path="metadata_sources/41586_2021_4237_MOESM5_ESM.xlsx",
            artifact_role="Supplementary Data 4: embryonic-cluster markers.",
            repository="https://www.nature.com/articles/s41586-021-04237-0",
            discovery_commit="NA",
            discovery_path="Supplementary Data 4",
            current_status="Public publisher supplement.",
        ),
        AuthorArtifact(
            name="TrackerSeq.MUC28072.md.RDS",
            url="https://keeper.mpdl.mpg.de/f/9b9d33a7503b4658aad1/?dl=1",
            expected_bytes=179_390,
            relative_path="metadata_sources/TrackerSeq.MUC28072.md.RDS",
            artifact_role="Author-posted TrackerSeq cell-to-lineage-barcode metadata; not CA301 annotations.",
            repository="https://github.com/mayer-lab/Bandler-et-al_lineage/issues/1",
            discovery_commit="NA",
            discovery_path="Issue 1 author reply, 2022-06-01",
            current_status="Public author link remains live.",
        ),
    )

    def __init__(self, run_dir: Path, source_root: Path):
        self.run_dir = run_dir.resolve()
        self.source_root = source_root.resolve()

    def run(self) -> None:
        discovery_rows = []
        download_rows = []
        now = datetime.now(timezone.utc).isoformat(timespec="seconds")
        for artifact in self.ARTIFACTS:
            destination = self.source_root / "Bandler2022" / artifact.relative_path
            action, digest = ImmutableHttpCache(artifact, destination).materialize()
            discovery_rows.append({
                "artifact": artifact.name,
                "artifact_role": artifact.artifact_role,
                "author_url": artifact.url,
                "repository": artifact.repository,
                "discovery_commit": artifact.discovery_commit,
                "discovery_path": artifact.discovery_path,
                "current_status": artifact.current_status,
            })
            download_rows.append({
                "artifact": artifact.name,
                "artifact_role": artifact.artifact_role,
                "resolved_path": destination,
                "bytes": destination.stat().st_size,
                "sha256": digest,
                "action": action,
                "completed_utc": now,
            })
        audit = self.run_dir / "Bandler2022" / "author_object_audit"
        atomic_tsv(audit / "author_object_discovery.tsv", (
            "artifact", "artifact_role", "author_url", "repository", "discovery_commit", "discovery_path", "current_status",
        ), discovery_rows)
        atomic_tsv(audit / "author_object_download.tsv", (
            "artifact", "artifact_role", "resolved_path", "bytes", "sha256", "action", "completed_utc",
        ), download_rows)
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
    parser.add_argument("--source-root", required=True, type=Path)
    args = parser.parse_args()
    BandlerAuthorObjectRecovery(args.run_dir, args.source_root).run()


if __name__ == "__main__":
    main()
