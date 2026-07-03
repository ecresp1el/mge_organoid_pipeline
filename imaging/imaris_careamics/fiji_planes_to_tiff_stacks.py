#!/usr/bin/env python3
"""Convert Fiji plane-per-file output into per-channel OME-TIFF stacks.

Inputs:
  --input-dir: directory containing Fiji plane files named like
    img_t1_z001_c1 or img_t1_z001_c1.tif.
  --channels: channel numbers to package, usually 1 2 for Fiji output.
  --dtype: expected plane dtype, either uint16 or uint8.

Outputs:
  One BigTIFF OME-TIFF stack per channel:
    <prefix>_c<channel>_stitched_stack.ome.tif

Important:
  The writer uses ome=True, axes=ZYX, and a known shape so Bio-Formats sees one
  stack series with 736 Z pages. Do not change this back to per-page writer
  calls; that created 736 Bio-Formats series and made Fiji ask for "series list"
  during open.
"""

import argparse
from pathlib import Path
import re

import tifffile


PLANE_RE = re.compile(r"img_t(?P<t>\d+)_z(?P<z>\d+)_c(?P<c>\d+)(?:\.tif)?$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Write one BigTIFF stack per channel from Fiji img_t*_z*_c* plane files."
    )
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--channels", type=int, nargs="+", default=[1, 2])
    parser.add_argument("--compression", default=None, choices=[None, "deflate"])
    parser.add_argument("--dtype", default="uint16", choices=["uint16", "uint8"])
    return parser.parse_args()


def channel_planes(input_dir: Path, channel: int) -> list[tuple[int, Path]]:
    """Return all planes for one channel sorted by numeric Z index."""
    planes = []
    for path in input_dir.iterdir():
        if not path.is_file():
            continue
        match = PLANE_RE.fullmatch(path.name)
        if not match:
            continue
        if int(match.group("c")) == channel:
            planes.append((int(match.group("z")), path))
    return sorted(planes)


def write_stack(
    planes: list[tuple[int, Path]],
    output_path: Path,
    compression: str | None,
    expected_dtype: str,
) -> None:
    """Write a single-series OME-TIFF Z stack from plane files.

    This streams planes from disk rather than loading the whole channel into
    memory. It also checks shape and dtype on every plane so partial or mixed
    outputs fail loudly.
    """
    if not planes:
        raise ValueError(f"No planes found for {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    first = tifffile.imread(planes[0][1])
    if first.dtype.name != expected_dtype:
        raise ValueError(f"{planes[0][1]} is {first.dtype}, expected {expected_dtype}")
    shape = (len(planes),) + first.shape

    def image_iter():
        for _, path in planes:
            image = tifffile.imread(path)
            if image.dtype != first.dtype:
                raise ValueError(f"{path} is {image.dtype}, expected {first.dtype}")
            if image.shape != first.shape:
                raise ValueError(f"{path} has shape {image.shape}, expected {first.shape}")
            yield image

    tifffile.imwrite(
        output_path,
        image_iter(),
        shape=shape,
        dtype=first.dtype,
        bigtiff=True,
        ome=True,
        metadata={"axes": "ZYX"},
        photometric="minisblack",
        compression=compression,
    )


def main() -> int:
    args = parse_args()
    for channel in args.channels:
        planes = channel_planes(args.input_dir, channel)
        if not planes:
            raise ValueError(f"No channel {channel} planes found in {args.input_dir}")
        expected = list(range(planes[0][0], planes[-1][0] + 1))
        observed = [z for z, _ in planes]
        if observed != expected:
            missing = sorted(set(expected) - set(observed))
            raise ValueError(f"Channel {channel} has missing z planes: {missing[:20]}")
        output_path = args.output_dir / f"{args.prefix}_c{channel}_stitched_stack.ome.tif"
        print(f"Writing channel {channel}: {len(planes)} planes -> {output_path}", flush=True)
        write_stack(planes, output_path, args.compression, args.dtype)
        print(f"Wrote {output_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
