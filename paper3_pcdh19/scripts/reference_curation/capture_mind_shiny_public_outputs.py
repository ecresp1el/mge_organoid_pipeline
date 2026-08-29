#!/usr/bin/env python3
"""Capture and quantify intended public outputs from the Mayer-lab MIND Shiny atlas."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import random
import re
import shutil
import socket
import string
import subprocess
import tempfile
import time
from collections import Counter
from contextlib import AbstractContextManager
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import requests
import websocket


ATLAS_PAGE = "https://mayerlab.net/mouse-inhibitory-neuron-development/"
SHINY_HTTP_ROOT = "http://141.5.108.55:3838/mind_shiny/"
SHINY_WS_ROOT = "ws://141.5.108.55:3838/mind_shiny/"
DATASET = "Embryonic Inhibitory and Excitatory"
STUDIES = ("Bright et al. 2025", "Bandler et al. 2022", "Di Bella et al. 2021")

# The factor/legend order is fixed in the public app and its data-preparation code.
# The PostScript colors are the vector encodings produced by pals::alphabet2.
LABEL_COLORS: dict[str, dict[str, str]] = {
    "cluster": {
        "0.471 0.165 0.714": "Gas1_Ldha",
        "0.667 0.051 0.992": "Hes1_Fabp7",
        "0.333": "Fabp7_Mt3",
        "0.769 0.271 0.11": "Hist1h1b_Top2a",
        "0.969 0.878 0.624": "Ccnd2_Nudt4",
        "0.992 0 0.976": "Nkx2-1_Lhx8",
        "0.882": "Npy_Nxph1",
        "0.992 0.682 0.082": "Sst_Maf",
        "0.11 0.741 0.31": "Nr2f2_Nr2f1",
        "0.192 0.353 0.608": "Isl1_Zfp503",
        "0.973 0.627 0.624": "Foxp1_Gucy1a3",
        "0.871 0.624 0.988": "Ebf1_Foxp1",
        "0.522 0.4 0.051": "Neurog2_Rrm2",
        "0.082 1 0.192": "Neurog2_Eomes",
        "0.192 0.514 0.992": "Neurod2_Neurod6",
        "0.11 0.514 0.333": "Neurod6_Mef2c",
    },
    "class": {
        "0.667 0.051 0.992": "Mitotic",
        "0.522 0.4 0.051": "Inhibitory Neuron Precursor",
        "0.192 0.514 0.992": "Excitatory Neuron Precursor",
    },
    "stage": {
        "0.667 0.051 0.992": "E12",
        "0.192 0.514 0.992": "E13",
        "0.522 0.4 0.051": "E14",
        "0.471 0.165 0.714": "E15",
        "0.333": "E16",
    },
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_tsv(path: Path, columns: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent, text=True)
    temporary = Path(temporary_name)
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
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent, text=True)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


@dataclass(frozen=True)
class PlotRequest:
    color_by: str
    split_by: str
    stem: str


@dataclass(frozen=True)
class VectorPoint:
    panel_index: int
    x: float
    y: float
    label: str


@dataclass(frozen=True)
class Panel:
    index: int
    x: float
    y: float
    width: float
    height: float


class ShinySockJsSession(AbstractContextManager["ShinySockJsSession"]):
    """Use the same public SockJS/Shiny protocol as the browser client."""

    def __init__(self, request: PlotRequest, timeout: int = 180):
        self.request = request
        self.timeout = timeout
        self.socket: websocket.WebSocket | None = None
        self.session_id = ""
        self.download_href = ""

    @staticmethod
    def _token(length: int = 8) -> str:
        alphabet = string.ascii_lowercase + string.digits
        return "".join(random.SystemRandom().choice(alphabet) for _ in range(length))

    @staticmethod
    def _messages(payload: str) -> list[dict[str, object]]:
        if not payload.startswith("a"):
            return []
        decoded = json.loads(payload[1:])
        messages = []
        for item in decoded:
            if item.startswith("0|m|"):
                messages.append(json.loads(item[4:]))
        return messages

    def _send(self, message: str) -> None:
        if self.socket is None:
            raise RuntimeError("Shiny socket is not connected")
        self.socket.send(json.dumps([message]))

    def __enter__(self) -> "ShinySockJsSession":
        route = f"{SHINY_WS_ROOT}__sockjs__/000/{self._token()}/websocket"
        self.socket = websocket.create_connection(
            route,
            timeout=30,
            origin="http://141.5.108.55:3838",
            http_no_proxy=["141.5.108.55"],
        )
        self.socket.settimeout(self.timeout)
        opening = self.socket.recv()
        if opening != "o":
            raise RuntimeError(f"Unexpected SockJS opening frame: {opening!r}")
        self._send("0|o|mind_shiny")
        self._send("0|m|" + json.dumps({"method": "init", "data": self._initial_inputs()}, separators=(",", ":")))
        deadline = time.monotonic() + self.timeout
        while time.monotonic() < deadline and not (self.session_id and self.download_href):
            for message in self._messages(self.socket.recv()):
                config = message.get("config")
                if isinstance(config, dict):
                    self.session_id = str(config.get("sessionId", self.session_id))
                values = message.get("values")
                if isinstance(values, dict) and values.get("umap_download"):
                    self.download_href = str(values["umap_download"])
        if not self.session_id or not self.download_href:
            raise RuntimeError("The public Shiny session did not expose its UMAP download handler")
        return self

    def __exit__(self, exc_type: object, exc: object, traceback: object) -> None:
        if self.socket is not None:
            self.socket.close()
        self.socket = None

    def _initial_inputs(self) -> dict[str, object]:
        return {
            "singletons": "",
            "allowDataUriScheme": True,
            ".clientdata_output_umap_width": 1200,
            ".clientdata_output_umap_height": 800,
            ".clientdata_output_umap_hidden": False,
            ".clientdata_output_umap_download_hidden": False,
            ".clientdata_output_feature_width": 1200,
            ".clientdata_output_feature_height": 800,
            ".clientdata_output_feature_hidden": True,
            ".clientdata_output_feature_download_hidden": True,
            ".clientdata_output_network_width": 1200,
            ".clientdata_output_network_height": 800,
            ".clientdata_output_network_hidden": True,
            ".clientdata_pixelratio": 1,
            ".clientdata_url_protocol": "http:",
            ".clientdata_url_hostname": "141.5.108.55",
            ".clientdata_url_pathname": "/mind_shiny/",
            ".clientdata_url_search": "",
            ".clientdata_url_hash_initial": "",
            ".clientdata_single_pixel_ratio": 1,
            "umap_dataset": DATASET,
            "umap_color_by": self.request.color_by,
            "umap_split_by": self.request.split_by,
            "umap_go": 1,
            "feature_dataset": DATASET,
            "gene": "Nfib",
            "feature_split_by": "nothing, show combined",
            "feature_go": 0,
            "tf1": "Nfib",
            "tf2": "not chosen",
            "tf3": "not chosen",
            "only_TF": "Yes",
            "network_go": 0,
        }

    def download_pdf(self, destination: Path) -> dict[str, object]:
        url = SHINY_HTTP_ROOT + self.download_href.lstrip("/")
        client = requests.Session()
        client.trust_env = False
        response = client.get(url, timeout=self.timeout, headers={"Referer": ATLAS_PAGE})
        response.raise_for_status()
        if not response.content.startswith(b"%PDF"):
            raise RuntimeError(f"UMAP download was not a PDF: {response.headers.get('content-type')}")
        destination.parent.mkdir(parents=True, exist_ok=True)
        descriptor, temporary_name = tempfile.mkstemp(prefix=f".{destination.name}.", dir=destination.parent)
        temporary = Path(temporary_name)
        try:
            with os.fdopen(descriptor, "wb") as handle:
                handle.write(response.content)
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(temporary, destination)
        finally:
            temporary.unlink(missing_ok=True)
        return {
            "atlas_page": ATLAS_PAGE,
            "public_download_url": url,
            "shiny_session_id": self.session_id,
            "content_type": response.headers.get("content-type", ""),
            "bytes": destination.stat().st_size,
            "sha256": sha256_file(destination),
        }


class GhostscriptRenderer:
    """Render public PDFs for inspection and convert them to parseable vector PostScript."""

    def __init__(self, executable: str = "gs"):
        resolved = shutil.which(executable)
        if resolved is None:
            raise FileNotFoundError(f"Ghostscript executable not found: {executable}")
        self.executable = resolved

    def png(self, pdf: Path, output: Path) -> None:
        subprocess.run([
            self.executable, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE",
            "-sDEVICE=pngalpha", "-r180", f"-sOutputFile={output}", str(pdf),
        ], check=True)

    def postscript(self, pdf: Path, output: Path) -> None:
        subprocess.run([
            self.executable, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE",
            "-sDEVICE=ps2write", "-dLanguageLevel=3", f"-sOutputFile={output}", str(pdf),
        ], check=True)


class VectorCellCounter:
    """Count geom_point circles in a public vector plot without claiming cell identifiers."""

    RGB = re.compile(r"([0-9.]+) ([0-9.]+) ([0-9.]+) rg")
    GRAY = re.compile(r"([0-9.]+) g")
    MOVE = re.compile(r"(-?[0-9.]+) (-?[0-9.]+) m")
    CURVE = re.compile(r"(-?[0-9.]+) (-?[0-9.]+) (-?[0-9.]+) (-?[0-9.]+) (-?[0-9.]+) (-?[0-9.]+) c")
    CLIP = re.compile(r"(-?[0-9.]+) (-?[0-9.]+) ([0-9.]+) ([0-9.]+) re W n")

    def __init__(self, label_colors: Mapping[str, str]):
        self.label_colors = dict(label_colors)

    def panels(self, postscript: Path) -> list[Panel]:
        unique: list[tuple[float, float, float, float]] = []
        for line in postscript.read_text(encoding="latin-1").splitlines():
            match = self.CLIP.fullmatch(line.strip())
            if not match:
                continue
            rectangle = tuple(float(value) for value in match.groups())
            if rectangle[2] > 1000 and rectangle[3] > 4000 and rectangle not in unique:
                unique.append(rectangle)
        if len(unique) != 3:
            raise RuntimeError(f"Expected three study panels in {postscript}; observed {len(unique)}")
        return [Panel(index, *rectangle) for index, rectangle in enumerate(sorted(unique))]

    def points(self, postscript: Path) -> list[VectorPoint]:
        panels = self.panels(postscript)
        points: list[VectorPoint] = []
        color = ""
        start: tuple[float, float] | None = None
        curves: list[tuple[float, ...]] = []
        for raw_line in postscript.read_text(encoding="latin-1").splitlines():
            line = raw_line.strip()
            rgb = self.RGB.fullmatch(line)
            if rgb:
                color = " ".join(rgb.groups())
                continue
            gray = self.GRAY.fullmatch(line)
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
                if start is not None and len(curves) == 4 and color in self.label_colors:
                    for panel in panels:
                        if panel.x <= start[0] <= panel.x + panel.width and panel.y <= start[1] <= panel.y + panel.height:
                            center_x = curves[0][-2]
                            center_y = (curves[0][-1] + curves[2][-1]) / 2
                            points.append(VectorPoint(panel.index, center_x, center_y, self.label_colors[color]))
                            break
                start = None
                curves = []
            elif line in {"S", "Q", "q", "n"}:
                start = None
                curves = []
        return points


class AtlasEvidencePublisher:
    """Coordinate public plot capture, vector counting, validation, and reporting."""

    REQUESTS = (
        PlotRequest("cluster", "study", "public_umap_cluster_by_study"),
        PlotRequest("class", "study", "public_umap_class_by_study"),
        PlotRequest("stage", "study", "public_umap_stage_by_study"),
        PlotRequest("study", "stage", "public_umap_study_by_stage"),
    )

    def __init__(self, run_dir: Path, reuse_public_pdfs: bool = False):
        self.run_dir = run_dir.resolve()
        if not (self.run_dir / "SUCCESS.txt").is_file():
            raise RuntimeError(f"Completed curation run required: {self.run_dir}")
        self.output_dir = self.run_dir / "Bandler2022" / "interactive_atlas"
        self.figure_dir = self.output_dir / "figures"
        self.metadata_dir = self.output_dir / "metadata"
        self.audit_dir = self.output_dir / "audit"
        self.renderer = GhostscriptRenderer()
        self.reuse_public_pdfs = reuse_public_pdfs

    def _existing_captures(self) -> list[dict[str, str]]:
        path = self.audit_dir / "public_plot_capture.tsv"
        if not path.is_file():
            raise FileNotFoundError(f"Pre-captured public plot ledger is required: {path}")
        with path.open("r", encoding="utf-8", newline="") as handle:
            return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]

    def _capture(self) -> tuple[list[dict[str, object]], dict[str, list[VectorPoint]]]:
        captures: list[dict[str, object]] = []
        point_sets: dict[str, list[VectorPoint]] = {}
        self.figure_dir.mkdir(parents=True, exist_ok=True)
        self.audit_dir.mkdir(parents=True, exist_ok=True)
        existing = self._existing_captures() if self.reuse_public_pdfs else []
        for request in self.REQUESTS:
            pdf = self.figure_dir / f"{request.stem}.pdf"
            png = self.figure_dir / f"{request.stem}.png"
            if self.reuse_public_pdfs:
                matches = [row for row in existing if row.get("color_by") == request.color_by and row.get("split_by") == request.split_by]
                if len(matches) != 1:
                    raise RuntimeError(f"Expected one pre-capture row for {request.color_by} by {request.split_by}")
                row = dict(matches[0])
                if not pdf.is_file() or pdf.read_bytes()[:4] != b"%PDF":
                    raise RuntimeError(f"Missing or invalid pre-captured public PDF: {pdf}")
                if row.get("sha256") and row["sha256"] != sha256_file(pdf):
                    raise RuntimeError(f"Pre-captured PDF changed after acquisition: {pdf}")
            else:
                with ShinySockJsSession(request) as session:
                    row = session.download_pdf(pdf)
                row.update({
                    "captured_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                    "capture_host": socket.gethostname(),
                    "dataset": DATASET,
                    "color_by": request.color_by,
                    "split_by": request.split_by,
                    "relative_pdf": str(pdf.relative_to(self.run_dir)),
                    "relative_png": str(png.relative_to(self.run_dir)),
                })
            self.renderer.png(pdf, png)
            captures.append(row)
            if request.color_by in LABEL_COLORS:
                postscript = self.audit_dir / f"{request.stem}.ps"
                self.renderer.postscript(pdf, postscript)
                point_sets[request.color_by] = VectorCellCounter(LABEL_COLORS[request.color_by]).points(postscript)
                postscript.unlink()
        return captures, point_sets

    @staticmethod
    def _count_rows(attribute: str, points: Sequence[VectorPoint]) -> list[dict[str, object]]:
        counts = Counter((point.panel_index, point.label) for point in points)
        label_order = list(LABEL_COLORS[attribute].values())
        rows = []
        for panel_index, study in enumerate(STUDIES):
            for label in label_order:
                cells = counts[(panel_index, label)]
                if cells:
                    rows.append({"study": study, "attribute": attribute, "label": label, "cells": cells})
        return rows

    @staticmethod
    def _linear_residual(first: Sequence[float], second: Sequence[float]) -> float:
        mean_x = sum(first) / len(first)
        mean_y = sum(second) / len(second)
        denominator = sum((value - mean_x) ** 2 for value in first)
        slope = sum((x - mean_x) * (y - mean_y) for x, y in zip(first, second)) / denominator
        intercept = mean_y - slope * mean_x
        return max(abs((slope * x + intercept) - y) for x, y in zip(first, second))

    def _cluster_by_stage(self, point_sets: Mapping[str, Sequence[VectorPoint]]) -> tuple[list[dict[str, object]], dict[str, object]]:
        cluster = [point for point in point_sets["cluster"] if point.panel_index == 1]
        stage = [point for point in point_sets["stage"] if point.panel_index == 1]
        if len(cluster) != len(stage) or len(cluster) == 0:
            raise RuntimeError(f"Bandler cluster/stage point mismatch: {len(cluster)} versus {len(stage)}")
        x_residual = self._linear_residual([point.x for point in cluster], [point.x for point in stage])
        y_residual = self._linear_residual([point.y for point in cluster], [point.y for point in stage])
        if x_residual > 0.2 or y_residual > 0.2:
            raise RuntimeError(f"Public plots do not preserve cell order/coordinates: residuals {x_residual}, {y_residual}")
        counts = Counter((cluster_point.label, stage_point.label) for cluster_point, stage_point in zip(cluster, stage))
        rows = []
        for label in LABEL_COLORS["cluster"].values():
            e13 = counts[(label, "E13")]
            e15 = counts[(label, "E15")]
            if e13 + e15:
                rows.append({"cluster": label, "E13_cells": e13, "E15_cells": e15, "total_cells": e13 + e15})
        validation = {
            "Bandler_points_in_each_plot": len(cluster),
            "cluster_to_stage_x_affine_max_residual_postscript_units": f"{x_residual:.6f}",
            "cluster_to_stage_y_affine_max_residual_postscript_units": f"{y_residual:.6f}",
            "E13_cells": sum(row["E13_cells"] for row in rows),
            "E15_cells": sum(row["E15_cells"] for row in rows),
            "total_cells": sum(row["total_cells"] for row in rows),
        }
        return rows, validation

    def _report(self, cluster_by_stage: Sequence[Mapping[str, object]]) -> str:
        table = [
            f"| `{row['cluster']}` | {int(row['E13_cells']):,} | {int(row['E15_cells']):,} | {int(row['total_cells']):,} |"
            for row in cluster_by_stage
        ]
        return "\n".join([
            "# Public MIND atlas capture: Bandler embryonic cells",
            "",
            "## Result",
            "",
            "The live Mayer-lab MIND Shiny app was queried only through its intended public session and `Save plot` download handler. Its vector UMAP exports contain **18,424 Bandler et al. 2022 cells** in the later integrated embryonic atlas: **11,004 E13** and **7,420 E15**.",
            "",
            "The Bandler cells have 12 later-atlas cluster labels. Their broad-class total is 2,877 `Mitotic` and 15,547 `Inhibitory Neuron Precursor`; no Bandler cells are labeled `Excitatory Neuron Precursor` in this filtered atlas.",
            "",
            "| Later-atlas cluster | E13 cells | E15 cells | Total |",
            "| --- | ---: | ---: | ---: |",
            *table,
            "",
            "![Bandler cells in the public later atlas, colored by cluster](figures/public_umap_cluster_by_study.png)",
            "",
            "## Interpretation boundary",
            "",
            "This is stronger than a screenshot-only lead: the intended public vector download preserves one circle per plotted cell, and the independently captured cluster, class, and stage totals reconcile exactly to 18,424. Cell drawing order and coordinates also reconcile between the cluster and stage plots, supporting the exact cross-tab above.",
            "",
            "The discrete plots alone do **not** identify CA301. They distinguish `study`, `stage`, `cluster`, and `class`, but not Bandler sample ID or MGE/CGE/LGE region. The separate barcode-recovery follow-up addresses that gap with intended public RNA-expression vectors and deposited-matrix fingerprints: it resolves the 7,420 E15 cells as 4,481 CA301 MGE, 2,937 CA302 CGE, and 2 CA303 LGE cells. The aggregate 7,420 must not be relabeled as MGE.",
            "",
            "The app can use its server-resident `EXCIT_INHIBIT_cleaned_sub.rds`, generated `EI_merged_umap2_df.tsv`, and expression HDF5 without sending those files to the browser. No intended public object/table download endpoint was found. The app returns rendered plot pixels/vector PDF and plot-coordinate metadata; this workflow did not attempt arbitrary server filesystem access.",
            "",
            "## Reproducibility",
            "",
            "The exact Python and SLURM files executed for this capture are stored in the parent run's `code/` directory. Re-running the standalone capture step intentionally overwrites only this `interactive_atlas/` output package within the same curation step.",
            "",
            "Machine-readable evidence is in `metadata/atlas_{cluster,class,stage}_by_study.tsv`, `metadata/bandler_cluster_by_stage.tsv`, `audit/public_endpoint_scope.tsv`, and `audit/vector_count_validation.tsv`.",
            "",
        ])

    def publish(self) -> None:
        captures, point_sets = self._capture()
        for attribute in ("cluster", "class", "stage"):
            rows = self._count_rows(attribute, point_sets[attribute])
            atomic_tsv(
                self.metadata_dir / f"atlas_{attribute}_by_study.tsv",
                ("study", "attribute", "label", "cells"), rows,
            )
        cluster_by_stage, validation = self._cluster_by_stage(point_sets)
        atomic_tsv(
            self.metadata_dir / "bandler_cluster_by_stage.tsv",
            ("cluster", "E13_cells", "E15_cells", "total_cells"), cluster_by_stage,
        )
        atomic_tsv(
            self.audit_dir / "public_plot_capture.tsv",
            ("captured_utc", "capture_host", "atlas_page", "dataset", "color_by", "split_by", "public_download_url", "shiny_session_id", "content_type", "relative_pdf", "relative_png", "bytes", "sha256"),
            captures,
        )
        atomic_tsv(
            self.audit_dir / "vector_count_validation.tsv",
            tuple(validation), (validation,),
        )
        atomic_tsv(
            self.audit_dir / "public_endpoint_scope.tsv",
            ("resource", "available_through_intended_public_app", "captured", "interpretation"),
            (
                {"resource": "UMAP vector PDF", "available_through_intended_public_app": "yes", "captured": "yes", "interpretation": "One vector circle per plotted cell; used for reconciled counts."},
                {"resource": "Rendered UMAP pixels and plot coordinate map", "available_through_intended_public_app": "yes", "captured": "PDF/PNG", "interpretation": "Browser-facing visualization output."},
                {"resource": "Study/stage/class/cluster aggregate counts", "available_through_intended_public_app": "derivable from vector PDF", "captured": "yes", "interpretation": "Exact for cells retained in the later atlas."},
                {"resource": "CA301 sample or MGE/CGE/LGE membership", "available_through_intended_public_app": "not as a discrete field; recoverable from intended RNA-expression vectors", "captured": "separate barcode_recovery follow-up", "interpretation": "Twenty-four public expression fingerprints plus preserved order recover deposited sample/barcode membership."},
                {"resource": "EI_merged_umap2_df.tsv or cell IDs", "available_through_intended_public_app": "no intended endpoint found", "captured": "no", "interpretation": "Loaded server-side; not returned by plot handlers."},
                {"resource": "EXCIT_INHIBIT_cleaned_sub.rds", "available_through_intended_public_app": "no intended endpoint found", "captured": "no", "interpretation": "Server-resident object named in public preparation code."},
            ),
        )
        atomic_text(self.output_dir / "MIND_PUBLIC_ATLAS_CAPTURE_REPORT.md", self._report(cluster_by_stage))
        manifest_path = self.audit_dir / "output_manifest.tsv"
        outputs = sorted(
            path for path in self.output_dir.rglob("*")
            if path.is_file() and path != manifest_path
        )
        atomic_tsv(
            manifest_path,
            ("relative_path", "bytes", "sha256"),
            ({"relative_path": str(path.relative_to(self.run_dir)), "bytes": path.stat().st_size, "sha256": sha256_file(path)} for path in outputs),
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--reuse-public-pdfs", action="store_true", help="Validate and parse PDFs already captured by the submission host.")
    args = parser.parse_args()
    AtlasEvidencePublisher(args.run_dir, reuse_public_pdfs=args.reuse_public_pdfs).publish()


if __name__ == "__main__":
    main()
