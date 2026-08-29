#!/usr/bin/env python3
"""Recover deposited Bandler E15 barcodes in the public later MIND atlas.

The public app is queried only through its intended RNA-expression PDF download
handler.  Binary expression fingerprints and preserved plotting order are
matched to the deposited CA301/CA302/CA303 matrices.  Later-atlas labels are
kept explicitly distinct from the original Bandler 2022 annotation taxonomy.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import json
import math
import os
import random
import re
import shutil
import socket
import stat
import string
import subprocess
import tempfile
import time
import urllib.request
from collections import Counter, defaultdict
from contextlib import AbstractContextManager
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import requests
import websocket


ATLAS_PAGE = "https://mayerlab.net/mouse-inhibitory-neuron-development/"
SHINY_HTTP_ROOT = "http://141.5.108.55:3838/mind_shiny/"
SHINY_WS_ROOT = "ws://141.5.108.55:3838/mind_shiny/"
DATASET = "Embryonic Inhibitory and Excitatory"
GENES = (
    "Nfib", "Hes1", "Fabp7", "Top2a", "Ccnd2", "Nkx2-1", "Lhx8", "Npy",
    "Nxph1", "Nr2f2", "Isl1", "Foxp1", "Gucy1a3", "Ebf1", "Sst", "Maf",
    "Lhx6", "Rnd2", "Zfp536", "Pcsk1n", "Nnat", "Ackr3", "Mef2c", "Erbb4",
)
SAMPLES = (
    ("CA301", "GSM5684876", "MGE", "GSM5684876_CA301_filtered_RNA_counts.RDS.gz"),
    ("CA302", "GSM5684875", "CGE", "GSM5684875_CA302_filtered_RNA_counts.RDS.gz"),
    ("CA303", "GSM5684874", "LGE", "GSM5684874_CA303_filtered_RNA_counts.RDS.gz"),
)
GEO_ROOT = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5684nnn"

LABEL_COLORS = {
    "cluster": {
        "0.471 0.165 0.714": "Gas1_Ldha", "0.667 0.051 0.992": "Hes1_Fabp7",
        "0.333": "Fabp7_Mt3", "0.769 0.271 0.11": "Hist1h1b_Top2a",
        "0.969 0.878 0.624": "Ccnd2_Nudt4", "0.992 0 0.976": "Nkx2-1_Lhx8",
        "0.882": "Npy_Nxph1", "0.992 0.682 0.082": "Sst_Maf",
        "0.11 0.741 0.31": "Nr2f2_Nr2f1", "0.192 0.353 0.608": "Isl1_Zfp503",
        "0.973 0.627 0.624": "Foxp1_Gucy1a3", "0.871 0.624 0.988": "Ebf1_Foxp1",
        "0.522 0.4 0.051": "Neurog2_Rrm2", "0.082 1 0.192": "Neurog2_Eomes",
        "0.192 0.514 0.992": "Neurod2_Neurod6", "0.11 0.514 0.333": "Neurod6_Mef2c",
    },
    "class": {
        "0.667 0.051 0.992": "Mitotic", "0.522 0.4 0.051": "Inhibitory Neuron Precursor",
        "0.192 0.514 0.992": "Excitatory Neuron Precursor",
    },
    "stage": {
        "0.667 0.051 0.992": "E12", "0.192 0.514 0.992": "E13",
        "0.522 0.4 0.051": "E14", "0.471 0.165 0.714": "E15", "0.333": "E16",
    },
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_tsv(path: Path, columns: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent, text=True)
    temporary = Path(name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(columns), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in rows:
                writer.writerow({column: row.get(column, "") for column in columns})
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def atomic_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent, text=True)
    temporary = Path(name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


@dataclass(frozen=True)
class Panel:
    index: int
    x: float
    y: float
    width: float
    height: float


@dataclass(frozen=True)
class PlotPoint:
    panel_index: int
    x: float
    y: float
    color: str


class DepositedMatrixCache:
    """Atomically acquire immutable deposited E15 matrices."""

    def __init__(self, source_root: Path):
        self.source_root = source_root.resolve()

    def acquire(self) -> list[dict[str, object]]:
        rows = []
        self.source_root.mkdir(parents=True, exist_ok=True)
        for sample, gsm, region, filename in SAMPLES:
            path = self.source_root / filename
            url = f"{GEO_ROOT}/{gsm}/suppl/{filename}"
            reused = path.is_file()
            if not reused:
                descriptor, name = tempfile.mkstemp(prefix=".download.", dir=self.source_root)
                os.close(descriptor)
                temporary = Path(name)
                try:
                    request = urllib.request.Request(url, headers={"User-Agent": "paper3-pcdh19-bandler-barcode-recovery/1.0"})
                    with urllib.request.urlopen(request, timeout=240) as response, temporary.open("wb") as output:
                        shutil.copyfileobj(response, output, 8 * 1024 * 1024)
                        output.flush()
                        os.fsync(output.fileno())
                    if temporary.stat().st_size == 0:
                        raise RuntimeError(f"Empty GEO download: {url}")
                    os.replace(temporary, path)
                    path.chmod(stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH)
                finally:
                    temporary.unlink(missing_ok=True)
            rows.append({
                "sample": sample, "gsm": gsm, "age": "E15.5 (GEO shorthand E15)",
                "region": region, "url": url, "resolved_path": str(path),
                "bytes": path.stat().st_size, "sha256": sha256_file(path),
                "cache_action": "reused" if reused else "downloaded", "acquired_utc": utc_now(),
            })
        return rows


class FeatureShinySession(AbstractContextManager["FeatureShinySession"]):
    """Use the browser-facing SockJS protocol and intended feature download."""

    def __init__(self, gene: str, timeout: int = 240):
        self.gene = gene
        self.timeout = timeout
        self.socket: websocket.WebSocket | None = None
        self.session_id = ""
        self.download_href = ""

    @staticmethod
    def _token(length: int = 8) -> str:
        return "".join(random.SystemRandom().choice(string.ascii_lowercase + string.digits) for _ in range(length))

    def _send(self, message: str) -> None:
        assert self.socket is not None
        self.socket.send(json.dumps([message]))

    @staticmethod
    def _messages(payload: str) -> list[dict[str, object]]:
        if not payload.startswith("a"):
            return []
        return [json.loads(item[4:]) for item in json.loads(payload[1:]) if item.startswith("0|m|")]

    def __enter__(self) -> "FeatureShinySession":
        route = f"{SHINY_WS_ROOT}__sockjs__/000/{self._token()}/websocket"
        self.socket = websocket.create_connection(
            route, timeout=30, origin="http://141.5.108.55:3838", http_no_proxy=["141.5.108.55"],
        )
        self.socket.settimeout(self.timeout)
        if self.socket.recv() != "o":
            raise RuntimeError("Unexpected SockJS opening frame")
        self._send("0|o|mind_shiny")
        inputs = {
            "singletons": "", "allowDataUriScheme": True,
            ".clientdata_output_umap_width": 1200, ".clientdata_output_umap_height": 800,
            ".clientdata_output_umap_hidden": True, ".clientdata_output_umap_download_hidden": True,
            ".clientdata_output_feature_width": 1200, ".clientdata_output_feature_height": 800,
            ".clientdata_output_feature_hidden": False, ".clientdata_output_feature_download_hidden": False,
            ".clientdata_output_network_width": 1200, ".clientdata_output_network_height": 800,
            ".clientdata_output_network_hidden": True, ".clientdata_pixelratio": 1,
            ".clientdata_url_protocol": "http:", ".clientdata_url_hostname": "141.5.108.55",
            ".clientdata_url_pathname": "/mind_shiny/", ".clientdata_url_search": "",
            ".clientdata_url_hash_initial": "", ".clientdata_single_pixel_ratio": 1,
            "umap_dataset": DATASET, "umap_color_by": "cluster", "umap_split_by": "study", "umap_go": 0,
            "feature_dataset": DATASET, "gene": self.gene, "feature_split_by": "study", "feature_go": 1,
            "tf1": "Nfib", "tf2": "not chosen", "tf3": "not chosen", "only_TF": "Yes", "network_go": 0,
        }
        self._send("0|m|" + json.dumps({"method": "init", "data": inputs}, separators=(",", ":")))
        deadline = time.monotonic() + self.timeout
        while time.monotonic() < deadline and not (self.session_id and self.download_href):
            for message in self._messages(self.socket.recv()):
                config = message.get("config")
                if isinstance(config, dict):
                    self.session_id = str(config.get("sessionId", self.session_id))
                values = message.get("values")
                if isinstance(values, dict) and values.get("feature_download"):
                    self.download_href = str(values["feature_download"])
        if not self.session_id or not self.download_href:
            raise RuntimeError(f"Feature handler was not exposed for {self.gene}")
        return self

    def __exit__(self, exc_type: object, exc: object, traceback: object) -> None:
        if self.socket is not None:
            self.socket.close()
        self.socket = None

    def download(self, destination: Path) -> dict[str, object]:
        url = SHINY_HTTP_ROOT + self.download_href.lstrip("/")
        client = requests.Session()
        client.trust_env = False
        response = client.get(url, timeout=self.timeout, headers={"Referer": ATLAS_PAGE})
        response.raise_for_status()
        if not response.content.startswith(b"%PDF"):
            raise RuntimeError(f"Feature output for {self.gene} was not PDF")
        descriptor, name = tempfile.mkstemp(prefix=f".{destination.name}.", dir=destination.parent)
        temporary = Path(name)
        try:
            with os.fdopen(descriptor, "wb") as handle:
                handle.write(response.content)
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(temporary, destination)
        finally:
            temporary.unlink(missing_ok=True)
        return {
            "gene": self.gene, "captured_utc": utc_now(), "capture_host": socket.gethostname(),
            "atlas_page": ATLAS_PAGE, "public_download_url": url, "shiny_session_id": self.session_id,
            "bytes": destination.stat().st_size, "sha256": sha256_file(destination),
            "relative_pdf": destination.name,
        }


class AcquisitionWorkflow:
    def __init__(self, run_dir: Path, source_root: Path):
        self.run_dir = run_dir.resolve()
        self.output = self.run_dir / "Bandler2022" / "interactive_atlas" / "barcode_recovery"
        self.feature_dir = self.output / "public_features"
        self.source_root = source_root.resolve()

    def run(self) -> None:
        if not (self.run_dir / "SUCCESS.txt").is_file():
            raise RuntimeError(f"Completed curation run required: {self.run_dir}")
        self.feature_dir.mkdir(parents=True, exist_ok=True)
        source_rows = DepositedMatrixCache(self.source_root).acquire()
        captures = []
        for index, gene in enumerate(GENES, start=1):
            destination = self.feature_dir / f"{gene}.pdf"
            error: Exception | None = None
            for attempt in range(1, 4):
                try:
                    with FeatureShinySession(gene) as session:
                        row = session.download(destination)
                    row["attempt"] = attempt
                    captures.append(row)
                    print(f"Captured {index}/{len(GENES)}: {gene}", flush=True)
                    break
                except Exception as caught:
                    error = caught
                    if attempt < 3:
                        time.sleep(2 * attempt)
            else:
                raise RuntimeError(f"Could not capture public expression plot for {gene}: {error}")
        atomic_tsv(
            self.output / "audit" / "source_matrix_manifest.tsv",
            ("sample", "gsm", "age", "region", "url", "resolved_path", "bytes", "sha256", "cache_action", "acquired_utc"),
            source_rows,
        )
        atomic_tsv(
            self.output / "audit" / "public_feature_capture.tsv",
            ("gene", "captured_utc", "capture_host", "atlas_page", "public_download_url", "shiny_session_id", "bytes", "sha256", "relative_pdf", "attempt"),
            captures,
        )


class VectorPlotParser:
    RGB = re.compile(r"([0-9.]+) ([0-9.]+) ([0-9.]+) rg")
    GRAY = re.compile(r"([0-9.]+) g")
    MOVE = re.compile(r"(-?[0-9.]+) (-?[0-9.]+) m")
    CURVE = re.compile(r"(-?[0-9.]+) (-?[0-9.]+) (-?[0-9.]+) (-?[0-9.]+) (-?[0-9.]+) (-?[0-9.]+) c")
    CLIP = re.compile(r"(-?[0-9.]+) (-?[0-9.]+) ([0-9.]+) ([0-9.]+) re W n")

    def panels(self, postscript: Path) -> list[Panel]:
        unique = []
        for line in postscript.read_text(encoding="latin-1").splitlines():
            match = self.CLIP.fullmatch(line.strip())
            if match:
                rectangle = tuple(float(value) for value in match.groups())
                if rectangle[2] > 1000 and rectangle[3] > 4000 and rectangle not in unique:
                    unique.append(rectangle)
        if len(unique) != 3:
            raise RuntimeError(f"Expected three study panels in {postscript}; observed {len(unique)}")
        return [Panel(index, *values) for index, values in enumerate(sorted(unique))]

    def points(self, postscript: Path, accepted_colors: Mapping[str, str] | None = None) -> list[PlotPoint]:
        panels = self.panels(postscript)
        points = []
        color = ""
        start = None
        curves: list[tuple[float, ...]] = []
        for raw_line in postscript.read_text(encoding="latin-1").splitlines():
            line = raw_line.strip()
            rgb = self.RGB.fullmatch(line)
            gray = self.GRAY.fullmatch(line)
            if rgb:
                color = " ".join(rgb.groups())
                continue
            if gray:
                color = gray.group(1)
                continue
            move = self.MOVE.fullmatch(line)
            if move:
                start = (float(move.group(1)), float(move.group(2)))
                curves = []
                continue
            curve = self.CURVE.fullmatch(line)
            if curve:
                curves.append(tuple(float(value) for value in curve.groups()))
                continue
            if line == "f":
                if start is not None and len(curves) == 4 and (accepted_colors is None or color in accepted_colors):
                    for panel in panels:
                        if panel.x <= start[0] <= panel.x + panel.width and panel.y <= start[1] <= panel.y + panel.height:
                            points.append(PlotPoint(panel.index, curves[0][-2], (curves[0][-1] + curves[2][-1]) / 2, color))
                            break
                start, curves = None, []
            elif line in {"S", "Q", "q", "n"}:
                start, curves = None, []
        return points


def srgb_to_lab(rgb: Sequence[float]) -> np.ndarray:
    values = np.array(rgb, dtype=float)
    values = np.where(values <= 0.04045, values / 12.92, ((values + 0.055) / 1.055) ** 2.4)
    xyz = np.array([
        0.4124564 * values[0] + 0.3575761 * values[1] + 0.1804375 * values[2],
        0.2126729 * values[0] + 0.7151522 * values[1] + 0.0721750 * values[2],
        0.0193339 * values[0] + 0.1191920 * values[1] + 0.9503041 * values[2],
    ]) / np.array([0.95047, 1.0, 1.08883])
    delta = 6 / 29
    transformed = np.where(xyz > delta ** 3, np.cbrt(xyz), xyz / (3 * delta ** 2) + 4 / 29)
    return np.array([116 * transformed[1] - 16, 500 * (transformed[0] - transformed[1]), 200 * (transformed[1] - transformed[2])])


LOW_LAB = srgb_to_lab((1.0, 1.0, 191 / 255))
HIGH_LAB = srgb_to_lab((165 / 255, 0.0, 38 / 255))
LAB_DIRECTION = HIGH_LAB - LOW_LAB


def expression_fraction(color: str) -> float:
    values = [float(value) for value in color.split()]
    if len(values) == 1:
        values *= 3
    lab = srgb_to_lab(values)
    return float(np.clip(np.dot(lab - LOW_LAB, LAB_DIRECTION) / np.dot(LAB_DIRECTION, LAB_DIRECTION), 0, 1))


def affine_max_residual(first: Sequence[float], second: Sequence[float]) -> float:
    x = np.asarray(first, dtype=float)
    y = np.asarray(second, dtype=float)
    slope, intercept = np.polyfit(x, y, 1)
    return float(np.max(np.abs((slope * x + intercept) - y)))


class OrderedFingerprintMatcher:
    """Find an exact binary-expression subsequence, then resolve duplicates by continuous color."""

    def __init__(self, genes: Sequence[str], atlas: np.ndarray, original: np.ndarray):
        self.genes = tuple(genes)
        self.atlas = atlas
        self.original = original
        self.atlas_binary = atlas > 0.002
        self.original_binary = original > 0

    @staticmethod
    def signatures(binary: np.ndarray) -> list[bytes]:
        return [np.packbits(row.astype(np.uint8)).tobytes() for row in binary]

    def longest_prefix(self, atlas_start: int, original_indices: Sequence[int]) -> int:
        atlas_signatures = self.signatures(self.atlas_binary)
        original_signatures = self.signatures(self.original_binary)
        cursor = atlas_start
        for original_index in original_indices:
            if cursor < len(atlas_signatures) and original_signatures[original_index] == atlas_signatures[cursor]:
                cursor += 1
        return cursor - atlas_start

    def match_segment(self, atlas_indices: Sequence[int], original_indices: Sequence[int]) -> tuple[list[int], list[dict[str, object]]]:
        atlas_sig_all = self.signatures(self.atlas_binary)
        original_sig_all = self.signatures(self.original_binary)
        atlas_sig = [atlas_sig_all[index] for index in atlas_indices]
        positions_by_signature: dict[bytes, list[int]] = defaultdict(list)
        for index in original_indices:
            positions_by_signature[original_sig_all[index]].append(index)

        earliest, cursor = [], original_indices[0] - 1
        for signature in atlas_sig:
            positions = positions_by_signature[signature]
            slot = bisect.bisect_right(positions, cursor)
            if slot == len(positions):
                raise RuntimeError("Atlas signatures are not a forward subsequence of the deposited matrix")
            cursor = positions[slot]
            earliest.append(cursor)
        latest, cursor = [0] * len(atlas_sig), original_indices[-1] + 1
        for i in range(len(atlas_sig) - 1, -1, -1):
            positions = positions_by_signature[atlas_sig[i]]
            slot = bisect.bisect_left(positions, cursor) - 1
            if slot < 0:
                raise RuntimeError("Atlas signatures are not a reverse subsequence of the deposited matrix")
            cursor = positions[slot]
            latest[i] = cursor

        seed_pairs = [(atlas_indices[i], earliest[i]) for i in range(len(earliest)) if earliest[i] == latest[i]]
        maxima = np.ones(len(self.genes), dtype=float)
        for gene_index in range(len(self.genes)):
            pairs = [(self.atlas[a, gene_index], self.original[o, gene_index]) for a, o in seed_pairs if self.atlas[a, gene_index] > 0.01]
            if pairs:
                x = np.array([pair[0] for pair in pairs])
                y = np.array([pair[1] for pair in pairs])
                maxima[gene_index] = max(float(np.dot(x, y) / np.dot(x, x)), 1e-8)
        original_scaled = np.clip(self.original / maxima, 0, 1)

        candidate_layers: list[list[int]] = []
        local_costs: list[np.ndarray] = []
        for i, signature in enumerate(atlas_sig):
            positions = positions_by_signature[signature]
            left = bisect.bisect_left(positions, earliest[i])
            right = bisect.bisect_right(positions, latest[i])
            candidates = positions[left:right]
            candidate_layers.append(candidates)
            diff = original_scaled[candidates, :] - self.atlas[atlas_indices[i], :]
            local_costs.append(np.sqrt(np.sum(diff * diff, axis=1)))

        dp = local_costs[0].copy()
        backpointers: list[np.ndarray] = [np.full(len(candidate_layers[0]), -1, dtype=np.int32)]
        for layer in range(1, len(candidate_layers)):
            previous_positions = candidate_layers[layer - 1]
            current_positions = candidate_layers[layer]
            new_dp = np.full(len(current_positions), np.inf)
            back = np.full(len(current_positions), -1, dtype=np.int32)
            previous_cursor = 0
            best_value = np.inf
            best_index = -1
            for current_index, position in enumerate(current_positions):
                while previous_cursor < len(previous_positions) and previous_positions[previous_cursor] < position:
                    if dp[previous_cursor] < best_value:
                        best_value = dp[previous_cursor]
                        best_index = previous_cursor
                    previous_cursor += 1
                if best_index >= 0:
                    new_dp[current_index] = best_value + local_costs[layer][current_index]
                    back[current_index] = best_index
            if not np.isfinite(new_dp).any():
                raise RuntimeError(f"No monotonic continuous-fingerprint path at atlas segment row {layer + 1}")
            dp = new_dp
            backpointers.append(back)
        selected_slots = [0] * len(candidate_layers)
        selected_slots[-1] = int(np.nanargmin(dp))
        for layer in range(len(candidate_layers) - 1, 0, -1):
            selected_slots[layer - 1] = int(backpointers[layer][selected_slots[layer]])
        selected = [candidate_layers[i][slot] for i, slot in enumerate(selected_slots)]

        diagnostics = []
        for i, (candidates, costs, chosen) in enumerate(zip(candidate_layers, local_costs, selected)):
            ordered_costs = np.sort(costs)
            chosen_cost = float(costs[candidates.index(chosen)])
            gap = float(ordered_costs[1] - ordered_costs[0]) if len(ordered_costs) > 1 else math.inf
            unique = earliest[i] == latest[i]
            chosen_is_local_best = chosen == candidates[int(np.argmin(costs))]
            if unique:
                method, status = "exact_expression_and_order", "definitive"
            elif chosen_is_local_best and gap >= 0.03:
                method, status = "continuous_expression_with_order", "high_confidence"
            else:
                method, status = "order_constrained_fingerprint", "ambiguous_among_equivalent_candidates"
            diagnostics.append({
                "candidate_count": len(candidates), "earliest_original_index": earliest[i] + 1,
                "latest_original_index": latest[i] + 1, "continuous_distance": f"{chosen_cost:.8f}",
                "local_distance_gap": "Inf" if math.isinf(gap) else f"{gap:.8f}",
                "match_method": method, "assignment_status": status,
            })
        return selected, diagnostics


class BarcodeRecoveryPublisher:
    def __init__(self, run_dir: Path, source_root: Path, fingerprints: Path):
        self.run_dir = run_dir.resolve()
        self.source_root = source_root.resolve()
        self.fingerprints = fingerprints.resolve()
        self.atlas_dir = self.run_dir / "Bandler2022" / "interactive_atlas"
        self.output = self.atlas_dir / "barcode_recovery"
        self.parser = VectorPlotParser()

    def _to_postscript(self, pdf: Path, postscript: Path) -> None:
        subprocess.run([
            "gs", "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=ps2write",
            "-dLanguageLevel=3", f"-sOutputFile={postscript}", str(pdf),
        ], check=True)

    def _parse_pdf(self, pdf: Path, accepted: Mapping[str, str] | None = None) -> tuple[list[Panel], list[PlotPoint]]:
        descriptor, name = tempfile.mkstemp(suffix=".ps", dir=self.output / "audit")
        os.close(descriptor)
        postscript = Path(name)
        try:
            self._to_postscript(pdf, postscript)
            return self.parser.panels(postscript), self.parser.points(postscript, accepted)
        finally:
            postscript.unlink(missing_ok=True)

    def _atlas_metadata(self) -> tuple[list[dict[str, object]], dict[str, object]]:
        point_sets = {}
        panels_by_attribute = {}
        for attribute in ("cluster", "class", "stage"):
            pdf = self.atlas_dir / "figures" / f"public_umap_{attribute}_by_study.pdf"
            panels, points = self._parse_pdf(pdf, LABEL_COLORS[attribute])
            panels_by_attribute[attribute] = panels
            point_sets[attribute] = [point for point in points if point.panel_index == 1]
        count = len(point_sets["stage"])
        if count != 18424 or any(len(point_sets[name]) != count for name in point_sets):
            raise RuntimeError(f"Expected 18,424 Bandler points in each discrete plot: {[len(value) for value in point_sets.values()]}")
        residuals = {}
        base = point_sets["cluster"]
        for attribute in ("class", "stage"):
            other = point_sets[attribute]
            residuals[f"cluster_to_{attribute}_x_affine_residual"] = affine_max_residual([p.x for p in base], [p.x for p in other])
            residuals[f"cluster_to_{attribute}_y_affine_residual"] = affine_max_residual([p.y for p in base], [p.y for p in other])
        panel = panels_by_attribute["cluster"][1]
        rows = []
        for index in range(count):
            rows.append({
                "atlas_bandler_order": index + 1,
                "later_atlas_cluster": LABEL_COLORS["cluster"][point_sets["cluster"][index].color],
                "later_atlas_class": LABEL_COLORS["class"][point_sets["class"][index].color],
                "later_atlas_stage": LABEL_COLORS["stage"][point_sets["stage"][index].color],
                "public_vector_umap_x": (base[index].x - panel.x) / panel.width,
                "public_vector_umap_y": (base[index].y - panel.y) / panel.height,
            })
        return rows, residuals

    def _feature_matrix(self, atlas_rows: Sequence[Mapping[str, object]]) -> tuple[np.ndarray, dict[str, object]]:
        matrix = np.zeros((len(atlas_rows), len(GENES)), dtype=float)
        residuals = {}
        base_x = np.array([float(row["public_vector_umap_x"]) for row in atlas_rows])
        base_y = np.array([float(row["public_vector_umap_y"]) for row in atlas_rows])
        for gene_index, gene in enumerate(GENES):
            panels, points = self._parse_pdf(self.output / "public_features" / f"{gene}.pdf")
            bandler = [point for point in points if point.panel_index == 1]
            if len(bandler) != len(atlas_rows):
                raise RuntimeError(f"{gene} feature plot has {len(bandler)} Bandler points, expected {len(atlas_rows)}")
            panel = panels[1]
            x = np.array([(point.x - panel.x) / panel.width for point in bandler])
            y = np.array([(point.y - panel.y) / panel.height for point in bandler])
            residuals[f"{gene}_x_affine_residual"] = affine_max_residual(base_x, x)
            residuals[f"{gene}_y_affine_residual"] = affine_max_residual(base_y, y)
            matrix[:, gene_index] = [expression_fraction(point.color) for point in bandler]
        return matrix, residuals

    def _read_original(self) -> tuple[list[dict[str, str]], np.ndarray]:
        with self.fingerprints.open("r", encoding="utf-8", newline="") as handle:
            rows = [dict(row) for row in csv.DictReader(handle, delimiter="\t")]
        matrix = np.array([[float(row[gene]) for gene in GENES] for row in rows], dtype=float)
        return rows, matrix

    def _plot(self, rows: Sequence[Mapping[str, object]]) -> None:
        ca301 = [row for row in rows if row["sample"] == "CA301"]
        clusters = list(LABEL_COLORS["cluster"].values())
        palette = dict(zip(clusters, plt.get_cmap("tab20").colors[:len(clusters)]))
        fig, ax = plt.subplots(figsize=(9.4, 7.2))
        for cluster in clusters:
            selected = [row for row in ca301 if row["later_atlas_cluster"] == cluster]
            if selected:
                ax.scatter(
                    [float(row["public_vector_umap_x"]) for row in selected],
                    [float(row["public_vector_umap_y"]) for row in selected],
                    s=4, alpha=0.78, linewidths=0, color=palette[cluster], label=f"{cluster} ({len(selected):,})",
                )
        ax.set_title("Recovered CA301 WT E15.5 MGE cells in the later MIND atlas")
        ax.set_xlabel("Public vector-plot UMAP x (normalized)")
        ax.set_ylabel("Public vector-plot UMAP y (normalized)")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, fontsize=8)
        fig.tight_layout()
        figure_dir = self.output / "figures"
        figure_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(figure_dir / "ca301_later_atlas_clusters.png", dpi=220, bbox_inches="tight")
        fig.savefig(figure_dir / "ca301_later_atlas_clusters.pdf", bbox_inches="tight")
        plt.close(fig)

    def run(self) -> None:
        (self.output / "audit").mkdir(parents=True, exist_ok=True)
        (self.output / "metadata").mkdir(parents=True, exist_ok=True)
        atlas_rows, discrete_residuals = self._atlas_metadata()
        atlas_features, feature_residuals = self._feature_matrix(atlas_rows)
        e15_indices = [i for i, row in enumerate(atlas_rows) if row["later_atlas_stage"] == "E15"]
        if len(e15_indices) != 7420:
            raise RuntimeError(f"Expected 7,420 Bandler E15 cells; observed {len(e15_indices)}")
        original_rows, original_features = self._read_original()
        sample_indices = {
            sample: [i for i, row in enumerate(original_rows) if row["sample"] == sample]
            for sample, _, _, _ in SAMPLES
        }
        matcher = OrderedFingerprintMatcher(GENES, atlas_features[e15_indices, :], original_features)
        segment_rows = []
        atlas_cursor = 0
        all_selected: list[int] = []
        all_diagnostics: list[dict[str, object]] = []
        for sample, gsm, region, _ in SAMPLES:
            matched = matcher.longest_prefix(atlas_cursor, sample_indices[sample])
            if matched == 0:
                continue
            segment_atlas = list(range(atlas_cursor, atlas_cursor + matched))
            selected, diagnostics = matcher.match_segment(segment_atlas, sample_indices[sample])
            all_selected.extend(selected)
            all_diagnostics.extend(diagnostics)
            segment_rows.append({
                "sample": sample, "gsm": gsm, "region": region,
                "deposited_cells": len(sample_indices[sample]), "later_atlas_E15_cells_recovered": matched,
                "deposited_cells_not_in_later_atlas": len(sample_indices[sample]) - matched,
                "first_atlas_E15_order": atlas_cursor + 1, "last_atlas_E15_order": atlas_cursor + matched,
            })
            atlas_cursor += matched
        if atlas_cursor != len(e15_indices):
            raise RuntimeError(f"Only {atlas_cursor} of {len(e15_indices)} E15 atlas cells formed ordered sample subsequences")
        observed = {row["sample"]: row["later_atlas_E15_cells_recovered"] for row in segment_rows}
        expected = {"CA301": 4481, "CA302": 2937, "CA303": 2}
        if observed != expected:
            raise RuntimeError(f"Recovered sample composition changed from validated expectation: {observed}")

        join_rows = []
        for e15_order, (atlas_index, original_index, diagnostic) in enumerate(zip(e15_indices, all_selected, all_diagnostics), start=1):
            atlas = atlas_rows[atlas_index]
            original = original_rows[original_index]
            join_rows.append({
                "atlas_E15_order": e15_order, "atlas_bandler_order": atlas["atlas_bandler_order"],
                "sample": original["sample"], "gsm": original["gsm"], "region": original["region"],
                "age": "E15.5 (atlas/GEO shorthand E15)", "genotype": "WT",
                "original_sample_order": original["original_sample_order"], "cell_id": original["cell_id"],
                "barcode": original["cell_id"].split("_", 1)[-1],
                "later_atlas_stage": atlas["later_atlas_stage"], "later_atlas_class": atlas["later_atlas_class"],
                "later_atlas_cluster": atlas["later_atlas_cluster"],
                "public_vector_umap_x": f"{float(atlas['public_vector_umap_x']):.9f}",
                "public_vector_umap_y": f"{float(atlas['public_vector_umap_y']):.9f}",
                **diagnostic,
            })
        columns = (
            "atlas_E15_order", "atlas_bandler_order", "sample", "gsm", "region", "age", "genotype",
            "original_sample_order", "cell_id", "barcode", "later_atlas_stage", "later_atlas_class",
            "later_atlas_cluster", "public_vector_umap_x", "public_vector_umap_y", "candidate_count",
            "earliest_original_index", "latest_original_index", "continuous_distance", "local_distance_gap",
            "match_method", "assignment_status",
        )
        atomic_tsv(self.output / "metadata" / "bandler_e15_later_atlas_barcode_join.tsv", columns, join_rows)
        atomic_tsv(
            self.output / "metadata" / "CA301_later_atlas_barcode_join.tsv", columns,
            (row for row in join_rows if row["sample"] == "CA301"),
        )
        atomic_tsv(
            self.output / "metadata" / "sample_recovery_summary.tsv",
            ("sample", "gsm", "region", "deposited_cells", "later_atlas_E15_cells_recovered", "deposited_cells_not_in_later_atlas", "first_atlas_E15_order", "last_atlas_E15_order"),
            segment_rows,
        )
        cluster_counts = Counter((row["sample"], row["later_atlas_class"], row["later_atlas_cluster"]) for row in join_rows)
        atomic_tsv(
            self.output / "metadata" / "sample_class_cluster_composition.tsv",
            ("sample", "later_atlas_class", "later_atlas_cluster", "cells"),
            ({"sample": key[0], "later_atlas_class": key[1], "later_atlas_cluster": key[2], "cells": value} for key, value in sorted(cluster_counts.items())),
        )
        status_counts = Counter((row["sample"], row["assignment_status"]) for row in join_rows)
        atomic_tsv(
            self.output / "audit" / "match_status_summary.tsv", ("sample", "assignment_status", "cells"),
            ({"sample": key[0], "assignment_status": key[1], "cells": value} for key, value in sorted(status_counts.items())),
        )
        validation = {
            "fingerprint_genes": len(GENES), "fingerprint_gene_names": "|".join(GENES),
            "Bandler_later_atlas_cells": len(atlas_rows), "Bandler_E15_cells": len(e15_indices),
            "joined_E15_cells": len(join_rows), "CA301_MGE_cells_recovered": observed["CA301"],
            "CA302_CGE_cells_recovered": observed["CA302"], "CA303_LGE_cells_recovered": observed["CA303"],
            "binary_expression_order_subsequence": "PASS", "expected_sample_composition": "PASS",
            "maximum_discrete_affine_residual": f"{max(discrete_residuals.values()):.9f}",
            "maximum_feature_affine_residual": f"{max(feature_residuals.values()):.9f}",
        }
        atomic_tsv(self.output / "audit" / "barcode_recovery_validation.tsv", tuple(validation), (validation,))
        atomic_tsv(
            self.output / "audit" / "fingerprint_gene_manifest.tsv", ("gene", "public_feature_pdf", "matching_use"),
            ({"gene": gene, "public_feature_pdf": f"public_features/{gene}.pdf", "matching_use": "binary expression plus continuous gradient position"} for gene in GENES),
        )
        self._plot(join_rows)
        ca301_counts = Counter(row["later_atlas_cluster"] for row in join_rows if row["sample"] == "CA301")
        confidence = Counter(row["assignment_status"] for row in join_rows)
        table = [f"| `{label}` | {count:,} |" for label, count in ca301_counts.most_common()]
        report = "\n".join([
            "# Bandler E15 barcode recovery from the public later MIND atlas", "", "## Result", "",
            "The missing sample/barcode bridge has been recovered for every one of the **7,420 Bandler E15 cells** retained in the public later atlas. The preserved sequence resolves to **4,481 CA301 WT MGE**, **2,937 CA302 WT CGE**, and **2 CA303 WT LGE** cells.", "",
            "The deposited matrices contain 4,516 CA301, 2,948 CA302, and 2,763 CA303 cells. Thus the later inhibitory/excitatory atlas retained 4,481/4,516 CA301 cells, 2,937/2,948 CA302 cells, and only 2/2,763 CA303 cells. The former public E15 total must therefore not be interpreted as 7,420 MGE cells; its directly recovered MGE component is 4,481 cells.", "",
            f"Across all E15 rows, {confidence.get('definitive', 0):,} assignments are uniquely fixed by expression plus order, {confidence.get('high_confidence', 0):,} are resolved by continuous expression under the order constraint, and {confidence.get('ambiguous_among_equivalent_candidates', 0):,} retain an explicit ambiguity flag. No ambiguous row is silently presented as definitive.", "",
            "## CA301 later-atlas composition", "", "| Later-atlas cluster | CA301 cells |", "| --- | ---: |", *table, "",
            "![Recovered CA301 cells in the later atlas](figures/ca301_later_atlas_clusters.png)", "",
            "## What these labels mean", "",
            "`later_atlas_class` and `later_atlas_cluster` are the Mayer-lab MIND atlas reanalysis labels exposed by the current public app. They are not represented as the original Bandler 2022 21-cluster assignments. This step solves the missing deposited-barcode-to-later-atlas join; recovery of the original 2022 per-cell taxonomy would still require an original author barcode table/object.", "",
            "The coordinates are normalized coordinates decoded from the intended public vector UMAP export, not the unpublished native numeric UMAP columns. They preserve the public plot geometry and permit a reproducible CA301-only plot.", "",
            "## Evidence and reproducibility", "",
            "For each public expression plot, the Bandler panel contains 18,424 cell circles in the same order and affine-equivalent UMAP geometry as the independently captured stage/class/cluster plots. Twenty-four genes provide binary and continuous expression fingerprints. The E15 sequence is an exact ordered subsequence of the deposited CA301, then CA302, then CA303 matrices. Machine-readable joins, sample/cluster summaries, confidence fields, source hashes, public-download hashes, and validation metrics are stored under `metadata/` and `audit/`.", "",
            "The exact Python, R, shell, and SLURM files executed are frozen in the parent run's `code/` directory. Re-running this follow-up intentionally overwrites only `Bandler2022/interactive_atlas/barcode_recovery/` within this curation run and regenerates the parent report.", "",
        ])
        atomic_text(self.output / "BANDLER_E15_BARCODE_RECOVERY_REPORT.md", report)
        manifest = self.output / "audit" / "output_manifest.tsv"
        outputs = sorted(path for path in self.output.rglob("*") if path.is_file() and path != manifest)
        atomic_tsv(
            manifest, ("relative_path", "bytes", "sha256"),
            ({"relative_path": str(path.relative_to(self.run_dir)), "bytes": path.stat().st_size, "sha256": sha256_file(path)} for path in outputs),
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    acquire = subparsers.add_parser("acquire")
    acquire.add_argument("--run-dir", required=True, type=Path)
    acquire.add_argument("--source-root", required=True, type=Path)
    publish = subparsers.add_parser("publish")
    publish.add_argument("--run-dir", required=True, type=Path)
    publish.add_argument("--source-root", required=True, type=Path)
    publish.add_argument("--fingerprints", required=True, type=Path)
    args = parser.parse_args()
    if args.command == "acquire":
        AcquisitionWorkflow(args.run_dir, args.source_root).run()
    else:
        BarcodeRecoveryPublisher(args.run_dir, args.source_root, args.fingerprints).run()


if __name__ == "__main__":
    main()

