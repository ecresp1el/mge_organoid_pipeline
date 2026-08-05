#!/usr/bin/env python3
"""Input staging, OME-TIFF packaging, and validation for Fiji stitching."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
from typing import Iterable

import numpy as np
import tifffile


PLANE_RE = re.compile(r"img_t(?P<t>\d+)_z(?P<z>\d+)_c(?P<c>\d+)(?:\.tif)?$")


def layout_tile_names(path: Path) -> list[str]:
    names: list[str] = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or line.startswith("dim"):
            continue
        name = line.split(";", 1)[0].strip()
        if not name or Path(name).name != name:
            raise ValueError(f"Unsafe or invalid tile name in {path}: {name!r}")
        names.append(name)
    if len(names) < 2:
        raise ValueError(f"Expected at least two tiles in {path}, found {len(names)}")
    if len(names) != len(set(names)):
        raise ValueError(f"Duplicate tile names in {path}")
    return names


def stage_inputs(args: argparse.Namespace) -> None:
    args.destination.mkdir(parents=True, exist_ok=True)
    shutil.copy2(args.layout, args.destination / "TileConfiguration.txt")
    staged = []
    for name in layout_tile_names(args.layout):
        source = (args.source_dir / name).resolve()
        if not source.is_file():
            raise FileNotFoundError(source)
        target = args.destination / name
        target.symlink_to(source)
        staged.append({"name": name, "source": str(source), "bytes": source.stat().st_size})
    (args.destination / "input_manifest.json").write_text(
        json.dumps({"layout": str(args.layout.resolve()), "tiles": staged}, indent=2) + "\n"
    )
    print(f"Staged {len(staged)} read-only tile symlinks in {args.destination}")


def make_smoke(args: argparse.Namespace) -> None:
    """Create an asymmetric, exactly known, two-channel 3D montage."""
    args.destination.mkdir(parents=True, exist_ok=True)
    z_count, channels, height, width = 8, 2, 112, 112
    zz, yy, xx = np.indices((z_count, height, width), dtype=np.uint16)
    global_data = np.empty((z_count, channels, height, width), dtype=np.uint16)
    global_data[:, 0] = 200 + zz * 90 + yy * 11 + xx * 7
    global_data[:, 1] = 400 + zz * 120 + yy * 5 + xx * 13
    global_data[:, 0, 12:34, 18:42] += 1200
    global_data[:, 1, 65:98, 73:105] += 1800
    specs = [
        ("synthetic_F1.tif", 0, 0),
        ("synthetic_F2.tif", 48, 0),
        ("synthetic_F0.tif", 0, 48),
        ("synthetic_F3.tif", 48, 48),
    ]
    for name, x, y in specs:
        tile = global_data[:, :, y : y + 64, x : x + 64]
        tifffile.imwrite(args.destination / name, tile, imagej=True, metadata={"axes": "ZCYX"})
    (args.destination / "TileConfiguration.txt").write_text(
        """# Synthetic 2x2 layout with 16-pixel overlaps
dim = 3

synthetic_F1.tif; ; (0.0, 0.0, 0.0)
synthetic_F2.tif; ; (48.0, 0.0, 0.0)
synthetic_F0.tif; ; (0.0, 48.0, 0.0)
synthetic_F3.tif; ; (48.0, 48.0, 0.0)
"""
    )
    np.save(args.destination / "smoke_ground_truth.npy", global_data)
    print(f"Wrote known volume and four overlapping tiles to {args.destination}")


def channel_planes(input_dir: Path, channel: int) -> list[tuple[int, Path]]:
    planes: list[tuple[int, Path]] = []
    for path in input_dir.iterdir():
        if not path.is_file():
            continue
        match = PLANE_RE.fullmatch(path.name)
        if match and int(match.group("c")) == channel:
            planes.append((int(match.group("z")), path))
    return sorted(planes)


def assert_contiguous(planes: list[tuple[int, Path]], channel: int) -> None:
    if not planes:
        raise ValueError(f"No Fiji planes found for channel {channel}")
    observed = [z for z, _ in planes]
    expected = list(range(observed[0], observed[-1] + 1))
    if observed != expected:
        missing = sorted(set(expected) - set(observed))
        raise ValueError(f"Channel {channel} has missing Z planes: {missing[:20]}")


def stack_planes(args: argparse.Namespace) -> None:
    """Stream Fiji plane files into one single-series OME-TIFF per channel."""
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary = []
    for channel in args.channels:
        planes = channel_planes(args.input_dir, channel)
        assert_contiguous(planes, channel)
        first = tifffile.imread(planes[0][1])
        if first.dtype != np.uint16:
            raise ValueError(f"{planes[0][1]} has dtype {first.dtype}; expected uint16")
        shape = (len(planes),) + first.shape

        def images() -> Iterable[np.ndarray]:
            for _, path in planes:
                image = tifffile.imread(path)
                if image.dtype != first.dtype or image.shape != first.shape:
                    raise ValueError(
                        f"Inconsistent plane {path}: {image.shape}/{image.dtype}; "
                        f"expected {first.shape}/{first.dtype}"
                    )
                yield image

        output = args.output_dir / f"{args.prefix}_c{channel}_stitched_stack.ome.tif"
        print(f"Writing {len(planes)} channel-{channel} planes to {output}", flush=True)
        tifffile.imwrite(
            output,
            images(),
            shape=shape,
            dtype=first.dtype,
            bigtiff=True,
            ome=True,
            metadata={"axes": "ZYX"},
            photometric="minisblack",
        )
        summary.append({"channel": channel, "path": str(output), "shape": shape})
    (args.output_dir / "stack_summary.json").write_text(json.dumps(summary, indent=2) + "\n")


def validate(args: argparse.Namespace) -> None:
    report: dict[str, object] = {"status": "PASS", "mode": args.mode, "stacks": []}
    plane_counts: dict[int, int] = {}
    stack_shapes: dict[int, tuple[int, ...]] = {}
    for channel in args.channels:
        planes = channel_planes(args.planes_dir, channel)
        assert_contiguous(planes, channel)
        plane_counts[channel] = len(planes)
        stack = args.stacks_dir / f"{args.prefix}_c{channel}_stitched_stack.ome.tif"
        with tifffile.TiffFile(stack) as tif:
            if len(tif.series) != 1 or tif.series[0].axes != "ZYX" or not tif.ome_metadata:
                raise ValueError(f"Invalid single-series OME-TIFF: {stack}")
            shape = tuple(tif.series[0].shape)
            dtype = str(tif.series[0].dtype)
            if shape[0] != len(planes) or dtype != "uint16":
                raise ValueError(f"Stack/plane mismatch for channel {channel}: {shape}, {len(planes)}, {dtype}")
            stack_shapes[channel] = shape
            report["stacks"].append(
                {
                    "channel": channel,
                    "path": str(stack),
                    "shape": shape,
                    "axes": "ZYX",
                    "dtype": dtype,
                    "ome": True,
                    "fiji_plane_count": len(planes),
                }
            )
    if len(set(plane_counts.values())) != 1 or len(set(stack_shapes.values())) != 1:
        raise ValueError(f"Channels disagree: planes={plane_counts}, shapes={stack_shapes}")

    if args.mode == "smoke":
        truth = np.load(args.inputs_dir / "smoke_ground_truth.npy")
        comparisons = []
        for channel in args.channels:
            stack = args.stacks_dir / f"{args.prefix}_c{channel}_stitched_stack.ome.tif"
            observed = tifffile.imread(stack)
            expected = truth[:, channel - 1]
            if observed.shape != expected.shape:
                raise ValueError(f"Synthetic stitched shape is {observed.shape}; expected {expected.shape}")
            difference = np.abs(observed.astype(np.int32) - expected.astype(np.int32))
            comparison = {
                "channel": channel,
                "max_abs_difference": int(difference.max()),
                "exact_fraction": float(np.mean(difference == 0)),
            }
            comparisons.append(comparison)
            if difference.max() != 0:
                raise ValueError(f"Synthetic channel {channel} was not reconstructed exactly")
        report["synthetic_ground_truth_comparison"] = comparisons

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def manifest(args: argparse.Namespace) -> None:
    rows = []
    for path in sorted(args.run_dir.rglob("*")):
        if not path.is_file() or path == args.output:
            continue
        relative = path.relative_to(args.run_dir)
        should_hash = str(relative).startswith(("code/", "config/", "results/stacks/"))
        if not args.checksum_outputs and str(relative).startswith("results/stacks/"):
            should_hash = False
        rows.append(
            {"path": str(relative), "bytes": path.stat().st_size, "sha256": sha256(path) if should_hash else None}
        )
    args.output.write_text(json.dumps({"files": rows}, indent=2) + "\n")
    print(f"Wrote manifest for {len(rows)} files to {args.output}")


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser()
    commands = root.add_subparsers(dest="command", required=True)
    stage = commands.add_parser("stage-inputs")
    stage.add_argument("--source-dir", type=Path, required=True)
    stage.add_argument("--layout", type=Path, required=True)
    stage.add_argument("--destination", type=Path, required=True)
    stage.set_defaults(function=stage_inputs)
    smoke = commands.add_parser("make-smoke")
    smoke.add_argument("--destination", type=Path, required=True)
    smoke.set_defaults(function=make_smoke)
    stack = commands.add_parser("stack")
    stack.add_argument("--input-dir", type=Path, required=True)
    stack.add_argument("--output-dir", type=Path, required=True)
    stack.add_argument("--prefix", required=True)
    stack.add_argument("--channels", type=int, nargs="+", required=True)
    stack.set_defaults(function=stack_planes)
    check = commands.add_parser("validate")
    check.add_argument("--mode", choices=["smoke", "production"], required=True)
    check.add_argument("--inputs-dir", type=Path, required=True)
    check.add_argument("--planes-dir", type=Path, required=True)
    check.add_argument("--stacks-dir", type=Path, required=True)
    check.add_argument("--prefix", required=True)
    check.add_argument("--channels", type=int, nargs="+", required=True)
    check.add_argument("--output", type=Path, required=True)
    check.set_defaults(function=validate)
    listing = commands.add_parser("manifest")
    listing.add_argument("--run-dir", type=Path, required=True)
    listing.add_argument("--output", type=Path, required=True)
    listing.add_argument("--checksum-outputs", action="store_true")
    listing.set_defaults(function=manifest)
    return root


def main() -> int:
    args = parser().parse_args()
    args.function(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
