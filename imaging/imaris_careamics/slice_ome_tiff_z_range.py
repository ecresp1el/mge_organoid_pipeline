#!/usr/bin/env python3
"""Copy a 1-based inclusive Z slice range from OME-TIFF stacks.

Inputs:
  One or more single-series OME-TIFF stacks with axes ZYX.

Output:
  One OME-TIFF stack per input, preserving dtype and XY dimensions while writing
  only the requested Z pages.

Use case:
  Trim the corrected BC43/realbive4 stitched deliverables to Fiji slices
  25..430 without rerunning stitching or median filtering.
"""

import argparse
from pathlib import Path

import tifffile


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Write Z-trimmed OME-TIFF stack copies.")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--z-start", required=True, type=int, help="1-based first Z slice to keep.")
    parser.add_argument("--z-end", required=True, type=int, help="1-based last Z slice to keep, inclusive.")
    parser.add_argument("--suffix", default="_z025_z430")
    parser.add_argument("inputs", nargs="+", type=Path)
    return parser.parse_args()


def output_name(input_path: Path, suffix: str) -> str:
    name = input_path.name
    if name.endswith(".ome.tif"):
        return name[:-8] + suffix + ".ome.tif"
    if name.endswith(".tif"):
        return name[:-4] + suffix + ".tif"
    return name + suffix + ".ome.tif"


def write_trimmed(input_path: Path, output_path: Path, z_start: int, z_end: int) -> None:
    z0 = z_start - 1
    z1 = z_end
    if z_start < 1 or z_end < z_start:
        raise ValueError(f"Invalid Z range {z_start}..{z_end}")

    with tifffile.TiffFile(input_path) as tif:
        if len(tif.series) != 1:
            raise ValueError(f"{input_path} has {len(tif.series)} series; expected 1")
        series = tif.series[0]
        if series.axes != "ZYX":
            raise ValueError(f"{input_path} axes are {series.axes}; expected ZYX")
        if len(series.shape) != 3:
            raise ValueError(f"{input_path} shape is {series.shape}; expected ZYX")
        if z_end > series.shape[0]:
            raise ValueError(f"{input_path} has only {series.shape[0]} Z slices")

        shape = (z_end - z_start + 1, series.shape[1], series.shape[2])
        dtype = series.dtype

        def planes():
            for page_index in range(z0, z1):
                yield tif.pages[page_index].asarray()

        output_path.parent.mkdir(parents=True, exist_ok=True)
        tifffile.imwrite(
            output_path,
            planes(),
            shape=shape,
            dtype=dtype,
            bigtiff=True,
            ome=True,
            metadata={"axes": "ZYX"},
            photometric="minisblack",
        )
        print(f"Wrote {output_path} shape={shape} dtype={dtype}", flush=True)


def main() -> int:
    args = parse_args()
    for input_path in args.inputs:
        output_path = args.output_dir / output_name(input_path, args.suffix)
        write_trimmed(input_path, output_path, args.z_start, args.z_end)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
