#!/usr/bin/env python3
"""Convert a TeraStitcher XML Stack layout to Fiji TileConfiguration.txt."""

import argparse
from pathlib import Path
from typing import List, Tuple
import xml.etree.ElementTree as ET


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Read TeraStitcher-style <Stack> entries and write a Fiji "
            "Grid/Collection Stitching TileConfiguration.txt file."
        )
    )
    parser.add_argument(
        "--xml",
        required=True,
        type=Path,
        help="Input TeraStitcher XML file.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output TileConfiguration file. Defaults to <xml-dir>/TileConfiguration.txt.",
    )
    parser.add_argument(
        "--no-verify-images",
        action="store_true",
        help="Do not require each referenced tile image to exist next to the XML.",
    )
    return parser.parse_args()


def parse_float_attr(stack: ET.Element, name: str) -> float:
    value = stack.get(name)
    if value is None:
        raise ValueError(f"Stack entry is missing required attribute {name}")
    try:
        return float(value)
    except ValueError as exc:
        raise ValueError(f"Stack attribute {name} is not numeric: {value!r}") from exc


def stack_image_name(stack: ET.Element) -> str:
    image = stack.get("IMG_REGEX")
    if not image:
        raise ValueError("Stack entry is missing required attribute IMG_REGEX")

    dirname = stack.get("DIR_NAME", ".")
    if dirname in ("", "."):
        return image
    return str(Path(dirname) / image)


TileRecord = Tuple[str, float, float, float]


def read_stacks(xml_path: Path) -> List[TileRecord]:
    tree = ET.parse(xml_path)
    root = tree.getroot()
    stacks = root.findall(".//Stack")
    if not stacks:
        raise ValueError(f"No <Stack> entries found in {xml_path}")

    records = []  # type: List[TileRecord]
    for stack in stacks:
        image = stack_image_name(stack)
        abs_h = parse_float_attr(stack, "ABS_H")
        abs_v = parse_float_attr(stack, "ABS_V")
        abs_d = parse_float_attr(stack, "ABS_D")
        records.append((image, abs_h, abs_v, abs_d))
    return records


def format_coord(value: float) -> str:
    return f"{value:.1f}"


def render_tile_configuration(records: List[TileRecord]) -> str:
    lines = [
        "# Define the number of dimensions we are working on",
        "dim = 3",
        "",
        "# Define the image coordinates",
    ]
    for image, abs_h, abs_v, abs_d in records:
        lines.append(
            f"{image}; ; "
            f"({format_coord(abs_h)}, {format_coord(abs_v)}, {format_coord(abs_d)})"
        )
    return "\n".join(lines) + "\n"


def verify_images(xml_path: Path, records: List[TileRecord]) -> None:
    missing = []
    for image, _, _, _ in records:
        image_path = xml_path.parent / image
        if not image_path.exists():
            missing.append(str(image_path))

    if missing:
        joined = "\n  ".join(missing)
        raise FileNotFoundError(f"Referenced tile image(s) not found:\n  {joined}")


def main() -> int:
    args = parse_args()
    xml_path = args.xml.expanduser().resolve()
    if not xml_path.exists():
        raise FileNotFoundError(f"XML file does not exist: {xml_path}")

    records = read_stacks(xml_path)
    if not args.no_verify_images:
        verify_images(xml_path, records)

    output = args.output.expanduser() if args.output else xml_path.parent / "TileConfiguration.txt"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(render_tile_configuration(records), encoding="utf-8")

    print(f"Wrote {output}")
    print(f"Tiles: {len(records)}")
    for image, abs_h, abs_v, abs_d in records:
        print(
            f"  {image}: "
            f"ABS_H={format_coord(abs_h)} ABS_V={format_coord(abs_v)} ABS_D={format_coord(abs_d)}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
