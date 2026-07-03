#!/usr/bin/env python3
"""Create a 1x3 Fiji-style Z-scroll panel movie from saved TIFF stacks.

Inputs:
  1. green_display.tif: single-channel Z stack with Fiji display settings.
  2. magenta_display.tif: single-channel Z stack with Fiji display settings.
  3. composite_display.tif: two-channel composite hyperstack with per-channel
     Fiji display settings.

Output:
  An MP4 movie with synchronized Z slices:
    [ PV-mNG ] [ BiVe3-dTom ] [ Merge ]

Optional output:
  A lossless RGB OME-TIFF stack of the same panel frames.

Important:
  This is a Z-scroll movie, not a projection and not a time movie. The script
  does not auto-contrast or percentile-scale. It uses ImageJ/Fiji display
  ranges and LUTs parsed by tifffile when available. If display ranges cannot
  be parsed, provide --green-min/--green-max and --magenta-min/--magenta-max.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import re
from typing import Iterable

import imageio.v2 as imageio
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import tifffile


LABELS = ("PV-mNG", "BiVe3-dTom", "Merge")


@dataclass
class FijiStack:
    path: Path
    data: np.ndarray
    axes: str
    metadata: dict
    description: str
    ranges: list[tuple[float, float]]
    luts: list[np.ndarray]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render Fiji display stacks as a labeled green|magenta|merge Z-scroll movie."
    )
    parser.add_argument("--green", required=True, type=Path, help="Single-channel green TIFF stack.")
    parser.add_argument("--magenta", required=True, type=Path, help="Single-channel magenta TIFF stack.")
    parser.add_argument("--composite", required=True, type=Path, help="Two-channel composite TIFF stack.")
    parser.add_argument("--output", type=Path, default=Path("panel_movie.mp4"))
    parser.add_argument(
        "--rgb-tiff",
        type=Path,
        default=None,
        help="Optional lossless RGB OME-TIFF stack output for the rendered panel frames.",
    )
    parser.add_argument("--fps", type=float, default=10.0)
    parser.add_argument("--border", type=int, default=10, help="Black separator width between panels.")
    parser.add_argument("--scale", type=float, default=1.0, help="Resize rendered RGB panels after LUT mapping.")
    parser.add_argument("--header-height", type=int, default=72)
    parser.add_argument("--font-size", type=int, default=42)
    parser.add_argument(
        "--font-stroke-width",
        type=int,
        default=1,
        help="White text stroke width for thicker panel labels.",
    )
    parser.add_argument("--codec", default="libx264", help="FFmpeg codec for MP4 output.")
    parser.add_argument("--pixelformat", default="yuv420p", help="FFmpeg output pixel format.")
    parser.add_argument("--quality", type=int, default=6, help="ImageIO ffmpeg quality, 0-10.")
    parser.add_argument("--preset", default="veryfast", help="x264 preset passed to FFmpeg.")
    parser.add_argument(
        "--ffmpeg-threads",
        type=int,
        default=4,
        help="Limit FFmpeg encoder threads; full-resolution panel frames are very large.",
    )
    parser.add_argument("--green-min", type=float, default=None)
    parser.add_argument("--green-max", type=float, default=None)
    parser.add_argument("--magenta-min", type=float, default=None)
    parser.add_argument("--magenta-max", type=float, default=None)
    return parser.parse_args()


def load_font(size: int) -> ImageFont.ImageFont:
    candidates = [
        "/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf",
        "/usr/share/fonts/liberation-sans/LiberationSans-Bold.ttf",
        "/usr/share/fonts/gnu-free/FreeSansBold.ttf",
    ]
    for candidate in candidates:
        if Path(candidate).exists():
            return ImageFont.truetype(candidate, size=size)
    return ImageFont.load_default()


def parse_description_ranges(description: str) -> list[tuple[float, float]]:
    values = {}
    for line in description.splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        key = key.strip().lower()
        try:
            values[key] = float(value.strip())
        except ValueError:
            continue
    if "min" in values and "max" in values:
        return [(values["min"], values["max"])]
    return []


def normalize_ranges(raw_ranges, channel_count: int) -> list[tuple[float, float]]:
    if raw_ranges is None:
        return []
    values = [float(value) for value in np.asarray(raw_ranges).ravel()]
    if len(values) >= 2 * channel_count:
        return [(values[2 * i], values[2 * i + 1]) for i in range(channel_count)]
    if channel_count == 1 and len(values) >= 2:
        return [(values[0], values[1])]
    return []


def imagej_metadata_ranges(metadata: dict, description: str, channel_count: int) -> list[tuple[float, float]]:
    for key in ("Ranges", "ranges", "display_ranges"):
        ranges = normalize_ranges(metadata.get(key), channel_count)
        if ranges:
            return ranges
    single = parse_description_ranges(description)
    if channel_count == 1 and single:
        return single
    return []


def normalize_lut(raw_lut) -> np.ndarray | None:
    arr = np.asarray(raw_lut, dtype=np.uint8)
    if arr.shape == (256, 3):
        return arr
    if arr.shape == (3, 256):
        return arr.T.copy()
    if arr.size == 768:
        flat = arr.ravel()
        return np.stack([flat[:256], flat[256:512], flat[512:768]], axis=1).astype(np.uint8)
    return None


def imagej_metadata_luts(metadata: dict, channel_count: int) -> list[np.ndarray]:
    raw_luts = metadata.get("LUTs")
    if raw_luts is None:
        raw_luts = metadata.get("luts")
    if raw_luts is None:
        raw_luts = []
    luts = []
    for raw_lut in raw_luts:
        lut = normalize_lut(raw_lut)
        if lut is not None:
            luts.append(lut)
    return luts[:channel_count]


def fallback_lut(color: str) -> np.ndarray:
    x = np.arange(256, dtype=np.uint8)
    lut = np.zeros((256, 3), dtype=np.uint8)
    if color == "green":
        lut[:, 1] = x
    elif color == "magenta":
        lut[:, 0] = x
        lut[:, 2] = x
    else:
        lut[:] = x[:, None]
    return lut


def open_tiff_stack(path: Path) -> FijiStack:
    with tifffile.TiffFile(path) as tif:
        series = tif.series[0]
        axes = series.axes
        metadata = dict(tif.imagej_metadata or {})
        description = tif.pages[0].description or ""
    try:
        data = tifffile.memmap(path)
    except Exception:
        with tifffile.TiffFile(path) as tif:
            data = tif.asarray()
    data = np.asarray(data)
    data_zcyx, axes_zcyx = to_zcyx(data, axes, metadata)
    channel_count = data_zcyx.shape[1]
    ranges = imagej_metadata_ranges(metadata, description, channel_count)
    luts = imagej_metadata_luts(metadata, channel_count)
    return FijiStack(path, data_zcyx, axes_zcyx, metadata, description, ranges, luts)


def int_metadata(metadata: dict, key: str, default: int) -> int:
    value = metadata.get(key, default)
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def to_zcyx(data: np.ndarray, axes: str, metadata: dict) -> tuple[np.ndarray, str]:
    """Normalize tifffile output to Z,C,Y,X.

    Handles normal stack axes such as ZYX, ZCYX, CZYX, TZCYX, and ImageJ
    hyperstacks that tifffile exposes as a flattened page axis.
    """
    if axes and len(axes) == data.ndim:
        arr = data
        ax = list(axes)
        if "T" in ax:
            t_index = ax.index("T")
            arr = np.take(arr, 0, axis=t_index)
            ax.pop(t_index)
        if "S" in ax:
            s_index = ax.index("S")
            if arr.shape[s_index] == 1:
                arr = np.take(arr, 0, axis=s_index)
                ax.pop(s_index)
        if "C" not in ax:
            arr = np.expand_dims(arr, axis=0)
            ax.insert(0, "C")
        if "Z" not in ax:
            arr = np.expand_dims(arr, axis=0)
            ax.insert(0, "Z")
        wanted = ["Z", "C", "Y", "X"]
        if all(name in ax for name in wanted):
            order = [ax.index(name) for name in wanted]
            return np.transpose(arr, order), "ZCYX"

    if data.ndim == 2:
        return data[np.newaxis, np.newaxis, :, :], "ZCYX"

    if data.ndim == 3:
        channels = int_metadata(metadata, "channels", 1)
        slices = int_metadata(metadata, "slices", data.shape[0] // max(channels, 1))
        frames = int_metadata(metadata, "frames", 1)
        if channels * slices * frames == data.shape[0] and channels > 1:
            return data.reshape(frames, slices, channels, data.shape[-2], data.shape[-1])[0], "ZCYX"
        return data[:, np.newaxis, :, :], "ZCYX"

    raise ValueError(f"Cannot normalize stack with shape={data.shape} axes={axes!r}")


def require_range(
    stack: FijiStack,
    channel: int,
    manual_min: float | None,
    manual_max: float | None,
    label: str,
) -> tuple[float, float]:
    if channel < len(stack.ranges):
        lo, hi = stack.ranges[channel]
        if hi > lo:
            return lo, hi
    if manual_min is not None and manual_max is not None and manual_max > manual_min:
        return manual_min, manual_max
    raise ValueError(
        f"Could not parse Fiji display range for {label} from {stack.path}. "
        "Provide manual --green-min/--green-max and --magenta-min/--magenta-max."
    )


def render_channel_rgb(plane: np.ndarray, display_range: tuple[float, float], lut: np.ndarray) -> np.ndarray:
    lo, hi = display_range
    scaled = (plane.astype(np.float32) - lo) / (hi - lo)
    index = np.clip(np.round(scaled * 255.0), 0, 255).astype(np.uint8)
    return lut[index]


def render_composite_rgb(
    stack: FijiStack,
    z_index: int,
    green_range: tuple[float, float],
    magenta_range: tuple[float, float],
    green_lut: np.ndarray,
    magenta_lut: np.ndarray,
) -> np.ndarray:
    if stack.data.shape[1] < 2:
        raise ValueError(f"Composite stack {stack.path} has {stack.data.shape[1]} channel(s), expected 2")
    green_rgb = render_channel_rgb(stack.data[z_index, 0], green_range, green_lut).astype(np.uint16)
    magenta_rgb = render_channel_rgb(stack.data[z_index, 1], magenta_range, magenta_lut).astype(np.uint16)
    return np.clip(green_rgb + magenta_rgb, 0, 255).astype(np.uint8)


def resize_if_requested(rgb: np.ndarray, scale: float) -> np.ndarray:
    if scale == 1.0:
        return rgb
    if scale <= 0:
        raise ValueError("--scale must be > 0")
    image = Image.fromarray(rgb)
    width = max(1, int(round(image.width * scale)))
    height = max(1, int(round(image.height * scale)))
    return np.asarray(image.resize((width, height), Image.Resampling.LANCZOS))


def add_header(
    rgb: np.ndarray,
    label: str,
    header_height: int,
    font: ImageFont.ImageFont,
    font_stroke_width: int,
) -> np.ndarray:
    image = Image.fromarray(rgb)
    frame = Image.new("RGB", (image.width, image.height + header_height), "black")
    frame.paste(image, (0, header_height))
    draw = ImageDraw.Draw(frame)
    bbox = draw.textbbox((0, 0), label, font=font, stroke_width=font_stroke_width)
    text_width = bbox[2] - bbox[0]
    text_height = bbox[3] - bbox[1]
    x = (image.width - text_width) // 2
    y = (header_height - text_height) // 2
    draw.text(
        (x, y),
        label,
        fill="white",
        font=font,
        stroke_width=font_stroke_width,
        stroke_fill="white",
    )
    return np.asarray(frame)


def concatenate_panels(panels: Iterable[np.ndarray], border: int) -> np.ndarray:
    panel_list = list(panels)
    heights = {panel.shape[0] for panel in panel_list}
    if len(heights) != 1:
        raise ValueError(f"Panel heights differ: {sorted(heights)}")
    if border < 0:
        raise ValueError("--border must be >= 0")
    separator = np.zeros((panel_list[0].shape[0], border, 3), dtype=np.uint8)
    pieces = []
    for i, panel in enumerate(panel_list):
        if i:
            pieces.append(separator)
        pieces.append(panel)
    return np.concatenate(pieces, axis=1)


def validate_z_counts(stacks: list[FijiStack]) -> int:
    z_counts = {stack.data.shape[0] for stack in stacks}
    if len(z_counts) != 1:
        raise ValueError(f"Input stacks have different Z counts: {sorted(z_counts)}")
    return z_counts.pop()


def render_panel_frame(
    z_index: int,
    green: FijiStack,
    magenta: FijiStack,
    composite: FijiStack,
    green_range: tuple[float, float],
    magenta_range: tuple[float, float],
    composite_green_range: tuple[float, float],
    composite_magenta_range: tuple[float, float],
    green_lut: np.ndarray,
    magenta_lut: np.ndarray,
    composite_green_lut: np.ndarray,
    composite_magenta_lut: np.ndarray,
    scale: float,
    border: int,
    header_height: int,
    font: ImageFont.ImageFont,
    font_stroke_width: int,
) -> np.ndarray:
    """Render one Z slice as the final labeled 1x3 RGB panel frame."""
    green_rgb = render_channel_rgb(green.data[z_index, 0], green_range, green_lut)
    magenta_rgb = render_channel_rgb(magenta.data[z_index, 0], magenta_range, magenta_lut)
    merge_rgb = render_composite_rgb(
        composite,
        z_index,
        composite_green_range,
        composite_magenta_range,
        composite_green_lut,
        composite_magenta_lut,
    )

    panels = [
        add_header(resize_if_requested(green_rgb, scale), LABELS[0], header_height, font, font_stroke_width),
        add_header(resize_if_requested(magenta_rgb, scale), LABELS[1], header_height, font, font_stroke_width),
        add_header(resize_if_requested(merge_rgb, scale), LABELS[2], header_height, font, font_stroke_width),
    ]
    return concatenate_panels(panels, border)


def main() -> int:
    args = parse_args()
    green = open_tiff_stack(args.green)
    magenta = open_tiff_stack(args.magenta)
    composite = open_tiff_stack(args.composite)
    z_count = validate_z_counts([green, magenta, composite])

    green_range = require_range(green, 0, args.green_min, args.green_max, "green_display")
    magenta_range = require_range(magenta, 0, args.magenta_min, args.magenta_max, "magenta_display")
    composite_green_range = require_range(composite, 0, args.green_min, args.green_max, "composite channel 1")
    composite_magenta_range = require_range(composite, 1, args.magenta_min, args.magenta_max, "composite channel 2")

    green_lut = green.luts[0] if green.luts else fallback_lut("green")
    magenta_lut = magenta.luts[0] if magenta.luts else fallback_lut("magenta")
    composite_green_lut = composite.luts[0] if len(composite.luts) >= 1 else fallback_lut("green")
    composite_magenta_lut = composite.luts[1] if len(composite.luts) >= 2 else fallback_lut("magenta")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    font = load_font(args.font_size)

    def frame_at(z_index: int) -> np.ndarray:
        return render_panel_frame(
            z_index,
            green,
            magenta,
            composite,
            green_range,
            magenta_range,
            composite_green_range,
            composite_magenta_range,
            green_lut,
            magenta_lut,
            composite_green_lut,
            composite_magenta_lut,
            args.scale,
            args.border,
            args.header_height,
            font,
            args.font_stroke_width,
        )

    writer = imageio.get_writer(
        args.output,
        fps=args.fps,
        codec=args.codec,
        pixelformat=args.pixelformat,
        quality=args.quality,
        macro_block_size=2,
        output_params=[
            "-preset",
            args.preset,
            "-threads",
            str(args.ffmpeg_threads),
            "-movflags",
            "+faststart",
        ],
        ffmpeg_log_level="info",
    )
    try:
        for z_index in range(z_count):
            frame = frame_at(z_index)
            writer.append_data(frame)
            if (z_index + 1) % 25 == 0 or z_index + 1 == z_count:
                print(f"Rendered {z_index + 1} / {z_count} Z slices", flush=True)
    finally:
        writer.close()

    print(f"Wrote MP4: {args.output}", flush=True)
    if args.rgb_tiff is not None:
        args.rgb_tiff.parent.mkdir(parents=True, exist_ok=True)
        first = frame_at(0)

        def frames():
            yield first
            for z_index in range(1, z_count):
                yield frame_at(z_index)

        tifffile.imwrite(
            args.rgb_tiff,
            frames(),
            shape=(z_count,) + first.shape,
            dtype=np.uint8,
            bigtiff=True,
            ome=True,
            metadata={"axes": "ZYXS"},
            photometric="rgb",
        )
        print(f"Wrote RGB TIFF stack: {args.rgb_tiff}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
