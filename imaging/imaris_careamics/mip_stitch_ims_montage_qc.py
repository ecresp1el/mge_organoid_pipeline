#!/usr/bin/env python3
"""Make max-projection stitched QC images for a small IMS montage.

Purpose:
  Run this before full-volume fusion when XML, stage, or filename ordering might
  disagree. It reads each IMS tile directly, generates per-channel max
  projections, and writes candidate stitched MIP images for visual inspection.

Inputs:
  --xml: TeraStitcher-style XML with Stack records.
  --input-dir: directory containing the IMS tile files referenced by the XML.
  --output-dir: QC directory for per-tile MIPs, stitched MIPs, RGB displays,
    and position TSVs.

Outputs:
  stitched_mip_xml_*.tif: direct XML ABS_H/ABS_V placement.
  stitched_mip_imaris_extent_*.tif: placement derived from IMS physical extents.
  stitched_mip_corrected_stage_*.tif: candidate that preserves stage spacing
    but flips the row orientation used for BC43/realbive4.
  positions.tsv / positions_corrected_stage.tsv: numeric audit trail.

This script is QC/audit code, not the full Fiji stitcher. The selected layout
must still be written to a Fiji TileConfiguration file and used by the Slurm
stitching job.
"""

import argparse
from dataclasses import dataclass
import os
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
import tifffile


DATASET_TEMPLATE = "DataSet/ResolutionLevel {level}/TimePoint 0/Channel {channel}/Data"


def sanitize_hdf5_plugin_path() -> None:
    """Point HDF5 at a valid plugin path before importing h5py.

    Great Lakes may expose a missing default plugin directory. Imaris IMS files
    can require hdf5plugin filters, so each read function imports hdf5plugin
    before h5py after this environment cleanup.
    """
    plugin_path = os.environ.get("HDF5_PLUGIN_PATH")
    existing = []
    if plugin_path:
        existing = [path for path in plugin_path.split(os.pathsep) if path and Path(path).exists()]
    if existing:
        os.environ["HDF5_PLUGIN_PATH"] = os.pathsep.join(existing)
        return

    fallback = Path(os.environ.get("TMPDIR", "/tmp")) / "imaris_mip_qc_empty_hdf5_plugins"
    fallback.mkdir(parents=True, exist_ok=True)
    os.environ["HDF5_PLUGIN_PATH"] = str(fallback)


@dataclass
class Tile:
    image: str
    xml_h: float
    xml_v: float
    xml_d: float
    row: int
    col: int


@dataclass
class ImsMetadata:
    size_x: int
    size_y: int
    ext_min0: float
    ext_min1: float
    ext_max0: float
    ext_max1: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Max-project IMS tiles and stitch MIP QC images from XML and IMS metadata."
    )
    parser.add_argument("--xml", required=True, type=Path)
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--level", type=int, default=0)
    parser.add_argument("--channels", type=int, nargs="+", default=[0, 1])
    parser.add_argument("--z-step", type=int, default=32)
    return parser.parse_args()


def attr_text(attrs, name: str) -> str:
    value = attrs[name]
    if isinstance(value, bytes):
        return value.decode("utf-8")
    arr = np.asarray(value)
    if arr.dtype.kind == "S":
        return b"".join(arr.tolist()).decode("utf-8")
    return str(value)


def parse_tiles(xml_path: Path):
    root = ET.parse(xml_path).getroot()
    tiles = []
    for stack in root.findall(".//Stack"):
        image = stack.get("IMG_REGEX")
        if not image:
            raise ValueError("Stack missing IMG_REGEX")
        tiles.append(
            Tile(
                image=image,
                xml_h=float(stack.get("ABS_H", "0")),
                xml_v=float(stack.get("ABS_V", "0")),
                xml_d=float(stack.get("ABS_D", "0")),
                row=int(stack.get("ROW", "0")),
                col=int(stack.get("COL", "0")),
            )
        )
    if not tiles:
        raise ValueError(f"No Stack entries found in {xml_path}")
    return tiles


def read_metadata(path: Path) -> ImsMetadata:
    try:
        import hdf5plugin  # noqa: F401
    except ImportError:
        pass
    import h5py

    with h5py.File(path, "r") as handle:
        attrs = handle["DataSetInfo/Image"].attrs
        return ImsMetadata(
            size_x=int(float(attr_text(attrs, "X"))),
            size_y=int(float(attr_text(attrs, "Y"))),
            ext_min0=float(attr_text(attrs, "ExtMin0")),
            ext_min1=float(attr_text(attrs, "ExtMin1")),
            ext_max0=float(attr_text(attrs, "ExtMax0")),
            ext_max1=float(attr_text(attrs, "ExtMax1")),
        )


def max_project(path: Path, level: int, channel: int, size_y: int, size_x: int, z_step: int):
    """Read one IMS channel and return a Z max projection.

    The raw HDF5 arrays can be padded to 1024 x 1024, while IMS metadata and
    Bio-Formats report the valid image size as 1020 x 996 for this sample.
    size_y/size_x crop the projection to that valid extent.
    """
    try:
        import hdf5plugin  # noqa: F401
    except ImportError:
        pass
    import h5py

    dataset_path = DATASET_TEMPLATE.format(level=level, channel=channel)
    with h5py.File(path, "r") as handle:
        data = handle[dataset_path]
        mip = np.zeros((min(size_y, data.shape[1]), min(size_x, data.shape[2])), dtype=data.dtype)
        for z0 in range(0, data.shape[0], z_step):
            z1 = min(z0 + z_step, data.shape[0])
            block = data[z0:z1, : mip.shape[0], : mip.shape[1]]
            mip = np.maximum(mip, block.max(axis=0))
    return mip


def normalize_positions(positions):
    min_x = min(x for x, _ in positions.values())
    min_y = min(y for _, y in positions.values())
    return {name: (int(round(x - min_x)), int(round(y - min_y))) for name, (x, y) in positions.items()}


def xml_positions(tiles):
    """Candidate layout using the XML ABS_H/ABS_V values literally."""
    return normalize_positions({tile.image: (tile.xml_h, tile.xml_v) for tile in tiles})


def extent_positions(tiles, metadata):
    # Imaris ExtMin0 increases left-to-right, while ExtMin1 is larger for the
    # upper row in this dataset. Convert physical units back to pixel offsets.
    px = np.median(
        [
            (meta.ext_max0 - meta.ext_min0) / meta.size_x
            for meta in metadata.values()
        ]
    )
    py = np.median(
        [
            (meta.ext_max1 - meta.ext_min1) / meta.size_y
            for meta in metadata.values()
        ]
    )
    min_x = min(meta.ext_min0 for meta in metadata.values())
    max_y = max(meta.ext_min1 for meta in metadata.values())
    positions = {}
    for tile in tiles:
        meta = metadata[tile.image]
        positions[tile.image] = ((meta.ext_min0 - min_x) / px, (max_y - meta.ext_min1) / py)
    return normalize_positions(positions)


def corrected_stage_positions(tiles, metadata):
    """Candidate layout used for the corrected BC43/realbive4 stitch.

    This keeps the left/right correction implied by IMS ExtMin0 and uses ExtMin1
    increasing downward in the display image. For this sample, that visual row
    orientation best preserved the expected organoid morphology in MIP QC.
    """
    px = np.median(
        [
            (meta.ext_max0 - meta.ext_min0) / meta.size_x
            for meta in metadata.values()
        ]
    )
    py = np.median(
        [
            (meta.ext_max1 - meta.ext_min1) / meta.size_y
            for meta in metadata.values()
        ]
    )
    min_x = min(meta.ext_min0 for meta in metadata.values())
    min_y = min(meta.ext_min1 for meta in metadata.values())
    positions = {}
    for tile in tiles:
        meta = metadata[tile.image]
        positions[tile.image] = ((meta.ext_min0 - min_x) / px, (meta.ext_min1 - min_y) / py)
    return normalize_positions(positions)


def stitch(mips, positions, tile_shape):
    tile_h, tile_w = tile_shape
    width = max(x + tile_w for x, _ in positions.values())
    height = max(y + tile_h for _, y in positions.values())
    canvas = np.zeros((height, width), dtype=np.uint16)
    for name, image in mips.items():
        x, y = positions[name]
        canvas[y : y + tile_h, x : x + tile_w] = np.maximum(
            canvas[y : y + tile_h, x : x + tile_w], image
        )
    return canvas


def scale_to_u8(image):
    """Display-only scaling for MIP contact sheets.

    This percentile scaling is only for QC PNG/RGB previews from this script.
    It is not used for the final requested 8-bit filtered stacks.
    """
    lo, hi = np.percentile(image[image > 0], [0.5, 99.8]) if np.any(image > 0) else (0, 1)
    if hi <= lo:
        hi = lo + 1
    return np.clip((image.astype(np.float32) - lo) * 255.0 / (hi - lo), 0, 255).astype(np.uint8)


def write_rgb(path: Path, green, red):
    green8 = scale_to_u8(green)
    red8 = scale_to_u8(red)
    rgb = np.zeros((green.shape[0], green.shape[1], 3), dtype=np.uint8)
    rgb[..., 0] = red8
    rgb[..., 1] = green8
    tifffile.imwrite(path, rgb, photometric="rgb")


def main() -> int:
    args = parse_args()
    sanitize_hdf5_plugin_path()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    tiles = parse_tiles(args.xml)
    metadata = {tile.image: read_metadata(args.input_dir / tile.image) for tile in tiles}

    print("Tile metadata and XML positions:")
    for tile in tiles:
        meta = metadata[tile.image]
        print(
            f"{tile.image}: XML=({tile.xml_h},{tile.xml_v},{tile.xml_d}) "
            f"ROW={tile.row} COL={tile.col} "
            f"size={meta.size_x}x{meta.size_y} "
            f"ExtMin0={meta.ext_min0} ExtMin1={meta.ext_min1}"
        )

    xml_pos = xml_positions(tiles)
    extent_pos = extent_positions(tiles, metadata)
    corrected_pos = corrected_stage_positions(tiles, metadata)
    print("\nPixel positions used for MIP stitch:")
    for tile in tiles:
        print(
            f"{tile.image}: xml={xml_pos[tile.image]} "
            f"imaris_extent={extent_pos[tile.image]} "
            f"corrected_stage={corrected_pos[tile.image]}"
        )

    first_meta = metadata[tiles[0].image]
    tile_shape = (first_meta.size_y, first_meta.size_x)

    stitched_by_mode = {"xml": {}, "imaris_extent": {}, "corrected_stage": {}}
    for channel in args.channels:
        mips = {}
        for tile in tiles:
            path = args.input_dir / tile.image
            mip = max_project(path, args.level, channel, first_meta.size_y, first_meta.size_x, args.z_step)
            mips[tile.image] = mip
            stem = Path(tile.image).stem
            tifffile.imwrite(args.output_dir / f"{stem}_ch{channel}_mip.tif", mip)

        for mode, positions in [
            ("xml", xml_pos),
            ("imaris_extent", extent_pos),
            ("corrected_stage", corrected_pos),
        ]:
            canvas = stitch(mips, positions, tile_shape)
            stitched_by_mode[mode][channel] = canvas
            tifffile.imwrite(args.output_dir / f"stitched_mip_{mode}_ch{channel}.tif", canvas)
            tifffile.imwrite(args.output_dir / f"stitched_mip_{mode}_ch{channel}_display8.tif", scale_to_u8(canvas))

    if 0 in stitched_by_mode["xml"] and 1 in stitched_by_mode["xml"]:
        write_rgb(
            args.output_dir / "stitched_mip_xml_rgb_display.tif",
            stitched_by_mode["xml"][0],
            stitched_by_mode["xml"][1],
        )
        write_rgb(
            args.output_dir / "stitched_mip_imaris_extent_rgb_display.tif",
            stitched_by_mode["imaris_extent"][0],
            stitched_by_mode["imaris_extent"][1],
        )
        write_rgb(
            args.output_dir / "stitched_mip_corrected_stage_rgb_display.tif",
            stitched_by_mode["corrected_stage"][0],
            stitched_by_mode["corrected_stage"][1],
        )

    (args.output_dir / "positions.tsv").write_text(
        "image\txml_x\txml_y\timaris_extent_x\timaris_extent_y\t"
        "corrected_stage_x\tcorrected_stage_y\tExtMin0\tExtMin1\tExtMax0\tExtMax1\n"
        + "\n".join(
            [
                "\t".join(
                    [
                        tile.image,
                        str(xml_pos[tile.image][0]),
                        str(xml_pos[tile.image][1]),
                        str(extent_pos[tile.image][0]),
                        str(extent_pos[tile.image][1]),
                        str(corrected_pos[tile.image][0]),
                        str(corrected_pos[tile.image][1]),
                        str(metadata[tile.image].ext_min0),
                        str(metadata[tile.image].ext_min1),
                        str(metadata[tile.image].ext_max0),
                        str(metadata[tile.image].ext_max1),
                    ]
                )
                for tile in tiles
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(f"\nWrote MIP QC outputs to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
