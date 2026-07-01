#!/usr/bin/env python3
"""Recover user-created Loupe Browser cell tracks from a .cloupe file.

This is an unsupported but read-only extractor for Loupe v8 block-store files.
It follows the header chain, reads the base barcode table, and exports cell
track assignments from the latest saved Loupe Browser layer.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import struct
import zlib
from pathlib import Path

import pandas as pd


PROJECT_ROOT_DEFAULT = "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder"
CLOUPE_DEFAULT = (
    "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/"
    "data/raw/loupe_exports/div90_recluster_annotations/D90 combined SO multi reclustered.cloupe"
)
SELECTED_TRACKS = [
    "Subcluster Progen & Astro",
    "Subpallial only Recluster",
    "Cortical only recluster",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cloupe", default=CLOUPE_DEFAULT)
    parser.add_argument(
        "--outdir",
        default=None,
        help="Default: cloupe parent / recovered_from_cloupe",
    )
    parser.add_argument("--selected-tracks", default=",".join(SELECTED_TRACKS))
    return parser.parse_args()


def read_header(handle, offset: int) -> dict:
    handle.seek(offset)
    raw = handle.read(4096)
    text = raw.split(b"\0", 1)[0].decode("utf-8")
    return json.loads(text)


def read_index(handle, header: dict) -> dict:
    block = header["indexBlock"]
    handle.seek(block["Start"])
    raw = handle.read(block["End"] - block["Start"])
    return json.loads(gzip.decompress(raw))


def read_header_chain(path: Path) -> list[tuple[int, dict, dict]]:
    out = []
    with path.open("rb") as handle:
        offset = 0
        seen = set()
        while offset < path.stat().st_size and offset not in seen:
            seen.add(offset)
            header = read_header(handle, offset)
            out.append((offset, header, read_index(handle, header)))
            next_offset = int(header.get("nextHeaderOffset", 0) or 0)
            if next_offset <= offset or next_offset >= path.stat().st_size:
                break
            offset = next_offset
    return out


def read_block(path: Path, block: dict) -> bytes:
    with path.open("rb") as handle:
        handle.seek(block["Start"])
        raw = handle.read(block["End"] - block["Start"])
    comp = int(block.get("CompressionType", 0) or 0)
    if comp == 0:
        return raw
    if comp == 1:
        return gzip.decompress(raw)
    if comp == 2:
        for func in (gzip.decompress, zlib.decompress):
            try:
                return func(raw)
            except Exception:
                pass
        return raw
    return raw


def read_fixed_strings(path: Path, block: dict) -> list[str]:
    data = read_block(path, block)
    n = int(block["ArraySize"])
    width = int(block["ArrayWidth"])
    if width <= 0:
        raise ValueError("Fixed-width string block has non-positive ArrayWidth")
    return [
        data[i * width : (i + 1) * width].rstrip(b"\0").decode("utf-8", errors="replace")
        for i in range(n)
    ]


def read_uint16_assignments(path: Path, block: dict) -> list[int | None]:
    data = read_block(path, block)
    n = int(block["ArraySize"])
    if len(data) < n * 2:
        raise ValueError(f"Assignment block too short: {len(data)} bytes for {n} uint16 values")
    vals = struct.unpack("<" + "H" * n, data[: n * 2])
    return [None if v == 65535 else int(v) for v in vals]


def rgb_to_hex(rgb: list[int] | tuple[int, int, int]) -> str:
    return "#" + "".join(f"{int(x):02X}" for x in rgb[:3])


def main() -> None:
    opt = parse_args()
    cloupe = Path(opt.cloupe)
    outdir = Path(opt.outdir) if opt.outdir else cloupe.parent / "recovered_from_cloupe"
    outdir.mkdir(parents=True, exist_ok=True)
    selected = {x.strip() for x in opt.selected_tracks.split(",") if x.strip()}

    chain = read_header_chain(cloupe)
    if len(chain) < 1:
        raise RuntimeError(f"No readable header found in {cloupe}")
    base_offset, base_header, base_index = chain[0]
    latest_offset, latest_header, latest_index = chain[-1]

    barcodes = read_fixed_strings(cloupe, base_index["Matrices"][0]["Barcodes"])
    rows = []
    track_rows = []

    for track_i, track in enumerate(latest_index.get("CellTracks", []) or []):
        name = track["Name"]
        metadata = json.loads(read_block(cloupe, track["Metadata"]).decode("utf-8"))
        groups = list(metadata.get("groups", []))
        colors = [rgb_to_hex(x) for x in metadata.get("colors", [])]
        assignments = read_uint16_assignments(cloupe, track["Values"])
        if len(assignments) != len(barcodes):
            raise ValueError(f"{name}: {len(assignments)} assignments for {len(barcodes)} barcodes")
        non_null = sum(v is not None for v in assignments)
        track_rows.append(
            {
                "track_index": track_i,
                "track_name": name,
                "parent_uuid": track.get("ParentUuid", ""),
                "override_uuid": track.get("OverrideUuid", ""),
                "n_groups": len(groups),
                "n_assigned_cells": non_null,
                "groups": "; ".join(groups),
                "colors": "; ".join(colors),
                "user_created": bool(metadata.get("userCreated", False)),
            }
        )
        for barcode, idx in zip(barcodes, assignments):
            if idx is None:
                label = None
                color = None
            elif 0 <= idx < len(groups):
                label = groups[idx]
                color = colors[idx] if idx < len(colors) else None
            else:
                label = f"UNKNOWN_INDEX_{idx}"
                color = None
            rows.append(
                {
                    "cell_id": barcode,
                    "track_name": name,
                    "label_index": idx,
                    "label": label,
                    "color_hex": color,
                }
            )

    long = pd.DataFrame(rows)
    inventory = pd.DataFrame(track_rows)
    inventory.to_csv(outdir / "loupe_recovered_celltrack_inventory.tsv", sep="\t", index=False)
    long.to_csv(outdir / "loupe_recovered_celltracks_long.tsv.gz", sep="\t", index=False)

    selected_long = long.loc[long["track_name"].isin(selected) & long["label"].notna()].copy()
    selected_long.to_csv(outdir / "loupe_recovered_selected_celltracks.tsv.gz", sep="\t", index=False)

    selected_summary = (
        selected_long.groupby(["track_name", "label_index", "label", "color_hex"], dropna=False)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["track_name", "label_index"])
    )
    selected_summary.to_csv(outdir / "loupe_recovered_selected_celltrack_summary.tsv", sep="\t", index=False)

    provenance = {
        "cloupe": str(cloupe),
        "headers": [
            {
                "offset": offset,
                "version": header.get("version"),
                "index_start": header.get("indexBlock", {}).get("Start"),
                "index_end": header.get("indexBlock", {}).get("End"),
                "next_header_offset": header.get("nextHeaderOffset"),
            }
            for offset, header, _ in chain
        ],
        "base_header_offset": base_offset,
        "latest_header_offset": latest_offset,
        "selected_tracks": sorted(selected),
        "note": "Read-only unsupported extraction from Loupe v8 block-store cell tracks.",
    }
    (outdir / "loupe_recovery_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"Recovered {len(inventory)} cell tracks from {cloupe}")
    print(f"Selected assigned rows: {selected_long.shape[0]}")
    print(f"Output: {outdir}")


if __name__ == "__main__":
    main()
