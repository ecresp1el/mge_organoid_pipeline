#!/usr/bin/env python3
"""Create a native-resolution 1x3 denoised Z-stack panel movie."""

from __future__ import annotations

import argparse
from pathlib import Path

import imageio.v2 as imageio
import numpy as np
import tifffile
from PIL import Image, ImageDraw, ImageFont


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create a 1x3 movie from denoised green/red OME-TIFF stacks: "
            "magenta, green, merged."
        )
    )
    parser.add_argument("--green", required=True, type=Path, help="Green denoised OME-TIFF.")
    parser.add_argument("--red", required=True, type=Path, help="Red denoised OME-TIFF.")
    parser.add_argument("--output", required=True, type=Path, help="Output MP4 path.")
    parser.add_argument("--fps", type=float, default=12.0, help="Movie frames per second.")
    parser.add_argument("--lower-percentile", type=float, default=0.5)
    parser.add_argument("--upper-percentile", type=float, default=99.8)
    parser.add_argument(
        "--percentile-z-step",
        type=int,
        default=1,
        help="Use every Nth Z plane for whole-stack percentile calculation.",
    )
    parser.add_argument(
        "--z-step",
        type=int,
        default=1,
        help="Write every Nth Z plane to the movie.",
    )
    parser.add_argument(
        "--start-z",
        type=int,
        default=0,
        help="First Z plane to include, 0-indexed.",
    )
    parser.add_argument(
        "--end-z",
        type=int,
        default=None,
        help="Stop before this Z plane, 0-indexed. Defaults to stack depth.",
    )
    parser.add_argument("--label-red", default="BiVe3-dTom")
    parser.add_argument("--label-green", default="PV-mNG")
    parser.add_argument("--label-merged", default="Merged")
    parser.add_argument(
        "--quality",
        type=int,
        default=9,
        help="imageio/ffmpeg quality, 0-10. Higher is better/larger.",
    )
    return parser.parse_args()


def open_zyx(path: Path) -> np.ndarray:
    arr = tifffile.memmap(path)
    arr = np.asarray(arr)
    arr = np.squeeze(arr)
    if arr.ndim != 3:
        raise ValueError(f"Expected a 3D ZYX stack after squeezing {path}, got {arr.shape}")
    return arr


def percentile_limits(
    arr: np.ndarray,
    lower: float,
    upper: float,
    z_step: int,
) -> tuple[float, float]:
    sample = arr[:: max(1, z_step)]
    lo, hi = np.percentile(sample, (lower, upper))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo = float(np.nanmin(sample))
        hi = float(np.nanmax(sample))
    if hi <= lo:
        hi = lo + 1.0
    return float(lo), float(hi)


def scale_uint8(slice_2d: np.ndarray, lo: float, hi: float) -> np.ndarray:
    scaled = (slice_2d.astype(np.float32, copy=False) - lo) / (hi - lo)
    return (np.clip(scaled, 0.0, 1.0) * 255.0).astype(np.uint8)


def find_font(size: int) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    for candidate in (
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
        "/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf",
        "/Library/Fonts/Arial Bold.ttf",
    ):
        if Path(candidate).exists():
            return ImageFont.truetype(candidate, size=size)
    return ImageFont.load_default()


def draw_centered_label(
    draw: ImageDraw.ImageDraw,
    text: str,
    x0: int,
    x1: int,
    y0: int,
    label_h: int,
    font: ImageFont.FreeTypeFont | ImageFont.ImageFont,
) -> None:
    bbox = draw.textbbox((0, 0), text, font=font)
    text_w = bbox[2] - bbox[0]
    text_h = bbox[3] - bbox[1]
    x = x0 + ((x1 - x0) - text_w) // 2
    y = y0 + (label_h - text_h) // 2 - 2
    draw.text((x, y), text, fill=(245, 245, 245), font=font)


def make_frame(
    green_slice: np.ndarray,
    red_slice: np.ndarray,
    green_limits: tuple[float, float],
    red_limits: tuple[float, float],
    labels: tuple[str, str, str],
    font: ImageFont.FreeTypeFont | ImageFont.ImageFont,
    label_h: int,
    gutter: int,
) -> np.ndarray:
    green_u8 = scale_uint8(green_slice, *green_limits)
    red_u8 = scale_uint8(red_slice, *red_limits)

    h, w = green_u8.shape
    panel_h = h + label_h
    panel_w = (3 * w) + (2 * gutter)
    frame = np.zeros((panel_h, panel_w, 3), dtype=np.uint8)

    x_red = 0
    x_green = w + gutter
    x_merged = 2 * (w + gutter)

    frame[:h, x_red : x_red + w, 0] = red_u8
    frame[:h, x_red : x_red + w, 2] = red_u8

    frame[:h, x_green : x_green + w, 1] = green_u8

    frame[:h, x_merged : x_merged + w, 0] = red_u8
    frame[:h, x_merged : x_merged + w, 1] = green_u8
    frame[:h, x_merged : x_merged + w, 2] = red_u8

    img = Image.fromarray(frame)
    draw = ImageDraw.Draw(img)
    y_label = h
    draw.rectangle((0, y_label, panel_w, panel_h), fill=(0, 0, 0))
    draw_centered_label(draw, labels[0], x_red, x_red + w, y_label, label_h, font)
    draw_centered_label(draw, labels[1], x_green, x_green + w, y_label, label_h, font)
    draw_centered_label(draw, labels[2], x_merged, x_merged + w, y_label, label_h, font)
    return np.asarray(img)


def main() -> None:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    green = open_zyx(args.green)
    red = open_zyx(args.red)
    if green.shape != red.shape:
        raise ValueError(f"Green and red shapes differ: {green.shape} vs {red.shape}")

    z, y, x = green.shape
    start_z = max(0, args.start_z)
    end_z = z if args.end_z is None else min(z, args.end_z)
    if start_z >= end_z:
        raise ValueError(f"Invalid Z range: start={start_z}, end={end_z}, stack depth={z}")

    print(f"Green: {args.green} shape={green.shape} dtype={green.dtype}")
    print(f"Red:   {args.red} shape={red.shape} dtype={red.dtype}")
    print("Calculating whole-stack percentile limits...")
    green_limits = percentile_limits(
        green, args.lower_percentile, args.upper_percentile, args.percentile_z_step
    )
    red_limits = percentile_limits(
        red, args.lower_percentile, args.upper_percentile, args.percentile_z_step
    )
    print(f"Green limits: {green_limits[0]:.6g}, {green_limits[1]:.6g}")
    print(f"Red limits:   {red_limits[0]:.6g}, {red_limits[1]:.6g}")
    print(f"Writing native panel movie: {args.output}")
    print(f"Frame size: {(3 * x) + 32} x {y + 88}; Z frames: {len(range(start_z, end_z, args.z_step))}")

    font = find_font(size=38)
    labels = (args.label_red, args.label_green, args.label_merged)
    with imageio.get_writer(
        args.output,
        fps=args.fps,
        codec="libx264",
        quality=args.quality,
        macro_block_size=1,
        output_params=["-pix_fmt", "yuv420p", "-crf", "12", "-preset", "slow"],
    ) as writer:
        for idx, zi in enumerate(range(start_z, end_z, args.z_step), start=1):
            frame = make_frame(
                green[zi],
                red[zi],
                green_limits,
                red_limits,
                labels,
                font,
                label_h=88,
                gutter=16,
            )
            writer.append_data(frame)
            if idx == 1 or idx % 25 == 0:
                print(f"Wrote frame {idx} at Z={zi}")

    print(f"Done: {args.output}")


if __name__ == "__main__":
    main()
