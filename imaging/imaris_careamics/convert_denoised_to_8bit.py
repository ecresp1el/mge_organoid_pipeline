#!/usr/bin/env python3
"""Convert CAREamics denoised OME-TIFF outputs to fixed-scale 8-bit display files.

This is a post-processing/QC export step. It does not modify the raw `.ims`
input or the quantitative denoised OME-TIFFs.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path
from typing import Any, Iterator

import numpy as np
import tifffile
from PIL import Image


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert green_denoised.ome.tif and red_denoised.ome.tif to uint8 "
            "display copies using fixed whole-stack percentile scaling."
        )
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        type=Path,
        help="full_run directory containing green_denoised.ome.tif and red_denoised.ome.tif.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory. Defaults to <input-dir>/display_8bit/.",
    )
    parser.add_argument("--green-name", default="green_denoised.ome.tif")
    parser.add_argument("--red-name", default="red_denoised.ome.tif")
    parser.add_argument("--lower-percentile", type=float, default=0.5)
    parser.add_argument("--upper-percentile", type=float, default=99.8)
    parser.add_argument(
        "--percentile-z-step",
        type=int,
        default=2,
        help=(
            "Use every Nth Z plane to estimate display percentiles. "
            "Default 2 is deterministic and safer for large stacks; use 1 for exact all-Z limits."
        ),
    )
    parser.add_argument(
        "--percentile-yx-step",
        type=int,
        default=2,
        help=(
            "Use every Nth Y/X pixel to estimate display percentiles. "
            "Default 2 is deterministic and safer for large stacks; use 1 for exact all-pixel limits."
        ),
    )
    parser.add_argument(
        "--no-merged-rgb",
        action="store_true",
        help="Skip the RGB merged green/magenta display OME-TIFF.",
    )
    parser.add_argument(
        "--compression",
        default=None,
        choices=("none", "zlib", "lzw", "deflate"),
        help="Optional TIFF compression. Default writes uncompressed files for speed.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing 8-bit display outputs.",
    )
    return parser.parse_args()


def open_memmap(path: Path) -> tuple[np.ndarray, str, str]:
    with tifffile.TiffFile(path) as tif:
        series = tif.series[0]
        axes = series.axes
        dtype = str(series.dtype)
    arr = tifffile.memmap(path)
    arr = np.asarray(arr)
    if arr.ndim not in (3, 4):
        raise ValueError(f"Expected ZYX or TZYX data in {path}, got shape {arr.shape} axes={axes}")
    if axes not in ("ZYX", "TZYX"):
        if arr.ndim == 3:
            axes = "ZYX"
        elif arr.ndim == 4:
            axes = "TZYX"
    return arr, axes, dtype


def sampled_view(arr: np.ndarray, axes: str, z_step: int, yx_step: int) -> np.ndarray:
    z_step = max(1, int(z_step))
    yx_step = max(1, int(yx_step))
    if axes == "ZYX":
        return arr[::z_step, ::yx_step, ::yx_step]
    if axes == "TZYX":
        return arr[:, ::z_step, ::yx_step, ::yx_step]
    raise ValueError(f"Unsupported axes for sampling: {axes}")


def percentile_limits(
    arr: np.ndarray,
    axes: str,
    lower: float,
    upper: float,
    z_step: int,
    yx_step: int,
) -> tuple[float, float, dict[str, Any]]:
    sample = np.asarray(sampled_view(arr, axes, z_step, yx_step), dtype=np.float32)
    finite = np.isfinite(sample)
    n_sampled = int(sample.size)
    n_finite = int(np.count_nonzero(finite))
    if n_finite == 0:
        raise ValueError("Cannot calculate display limits: sampled pixels are all non-finite.")
    lo, hi = np.nanpercentile(sample, [lower, upper])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        finite_values = sample[finite]
        lo = float(np.min(finite_values))
        hi = float(np.max(finite_values))
    if hi <= lo:
        hi = lo + 1.0
    metadata = {
        "lower_percentile": float(lower),
        "upper_percentile": float(upper),
        "percentile_z_step": int(max(1, z_step)),
        "percentile_yx_step": int(max(1, yx_step)),
        "sample_shape": list(map(int, sample.shape)),
        "sampled_pixels": n_sampled,
        "sampled_finite_pixels": n_finite,
        "sampled_nonfinite_pixels": int(n_sampled - n_finite),
        "limits": [float(lo), float(hi)],
    }
    return float(lo), float(hi), metadata


def scale_plane_uint8(plane: np.ndarray, lo: float, hi: float) -> np.ndarray:
    plane_float = np.asarray(plane, dtype=np.float32)
    scaled = (np.nan_to_num(plane_float, nan=lo, posinf=hi, neginf=lo) - lo) / (hi - lo)
    return np.round(np.clip(scaled, 0.0, 1.0) * 255.0).astype(np.uint8)


def iter_indices(shape: tuple[int, ...], axes: str) -> Iterator[tuple[int, ...]]:
    if axes == "ZYX":
        for z in range(shape[0]):
            yield (z,)
    elif axes == "TZYX":
        for t in range(shape[0]):
            for z in range(shape[1]):
                yield (t, z)
    else:
        raise ValueError(f"Unsupported axes: {axes}")


def blank_max_projection(shape: tuple[int, ...], axes: str) -> np.ndarray:
    if axes == "ZYX":
        return np.zeros(shape[-2:], dtype=np.uint8)
    if axes == "TZYX":
        return np.zeros((shape[0], *shape[-2:]), dtype=np.uint8)
    raise ValueError(f"Unsupported axes: {axes}")


def update_max_projection(max_proj: np.ndarray, plane: np.ndarray, index: tuple[int, ...], axes: str) -> None:
    if axes == "ZYX":
        np.maximum(max_proj, plane, out=max_proj)
    else:
        np.maximum(max_proj[index[0]], plane, out=max_proj[index[0]])


def write_uint8_stack(
    path: Path,
    arr: np.ndarray,
    axes: str,
    limits: tuple[float, float],
    compression: str | None,
    overwrite: bool,
) -> dict[str, Any]:
    if path.exists() and not overwrite:
        raise FileExistsError(f"Output exists; pass --overwrite to replace it: {path}")

    lo, hi = limits
    stats = {
        "output": str(path),
        "axes": axes,
        "shape": list(map(int, arr.shape)),
        "dtype": "uint8",
        "input_nonfinite_pixels": 0,
        "output_zero_pixels": 0,
        "output_255_pixels": 0,
        "planes_written": 0,
    }
    max_proj = blank_max_projection(tuple(map(int, arr.shape)), axes)

    def planes() -> Iterator[np.ndarray]:
        for index in iter_indices(tuple(map(int, arr.shape)), axes):
            plane = arr[index]
            stats["input_nonfinite_pixels"] += int(plane.size - np.count_nonzero(np.isfinite(plane)))
            out = scale_plane_uint8(plane, lo, hi)
            stats["output_zero_pixels"] += int(np.count_nonzero(out == 0))
            stats["output_255_pixels"] += int(np.count_nonzero(out == 255))
            stats["planes_written"] += 1
            update_max_projection(max_proj, out, index, axes)
            yield out

    tifffile.imwrite(
        path,
        data=planes(),
        shape=arr.shape,
        dtype=np.uint8,
        bigtiff=True,
        ome=True,
        metadata={"axes": axes},
        photometric="minisblack",
        compression=None if compression in (None, "none") else compression,
    )
    stats["total_pixels"] = int(np.prod(arr.shape))
    stats["output_zero_percent"] = 100.0 * stats["output_zero_pixels"] / max(stats["total_pixels"], 1)
    stats["output_255_percent"] = 100.0 * stats["output_255_pixels"] / max(stats["total_pixels"], 1)
    stats["max_projection_shape"] = list(map(int, max_proj.shape))
    stats["_max_projection_array"] = max_proj
    return stats


def write_rgb_stack(
    path: Path,
    green: np.ndarray,
    red: np.ndarray,
    axes: str,
    green_limits: tuple[float, float],
    red_limits: tuple[float, float],
    compression: str | None,
    overwrite: bool,
) -> dict[str, Any]:
    if path.exists() and not overwrite:
        raise FileExistsError(f"Output exists; pass --overwrite to replace it: {path}")
    if green.shape != red.shape:
        raise ValueError(f"Green/red shapes differ: {green.shape} vs {red.shape}")

    rgb_shape = tuple(map(int, green.shape)) + (3,)
    rgb_axes = axes + "S"
    stats = {
        "output": str(path),
        "axes": rgb_axes,
        "shape": list(rgb_shape),
        "dtype": "uint8",
        "planes_written": 0,
    }
    max_green = blank_max_projection(tuple(map(int, green.shape)), axes)
    max_red = blank_max_projection(tuple(map(int, red.shape)), axes)

    def planes() -> Iterator[np.ndarray]:
        for index in iter_indices(tuple(map(int, green.shape)), axes):
            green_u8 = scale_plane_uint8(green[index], *green_limits)
            red_u8 = scale_plane_uint8(red[index], *red_limits)
            rgb = np.zeros((*green_u8.shape, 3), dtype=np.uint8)
            rgb[..., 0] = red_u8
            rgb[..., 1] = green_u8
            rgb[..., 2] = red_u8
            update_max_projection(max_green, green_u8, index, axes)
            update_max_projection(max_red, red_u8, index, axes)
            stats["planes_written"] += 1
            yield rgb

    tifffile.imwrite(
        path,
        data=planes(),
        shape=rgb_shape,
        dtype=np.uint8,
        bigtiff=True,
        ome=True,
        metadata={"axes": rgb_axes},
        photometric="rgb",
        compression=None if compression in (None, "none") else compression,
    )
    stats["max_green_projection_shape"] = list(map(int, max_green.shape))
    stats["max_red_projection_shape"] = list(map(int, max_red.shape))
    stats["_max_green_projection_array"] = max_green
    stats["_max_red_projection_array"] = max_red
    return stats


def save_projection_pngs(
    outdir: Path,
    axes: str,
    green_max: np.ndarray,
    red_max: np.ndarray,
    prefix: str = "denoised_8bit",
) -> list[str]:
    written: list[str] = []
    if axes == "TZYX":
        green_max_2d = np.max(green_max, axis=0)
        red_max_2d = np.max(red_max, axis=0)
    else:
        green_max_2d = green_max
        red_max_2d = red_max

    green_path = outdir / f"{prefix}_green_max_projection.png"
    red_path = outdir / f"{prefix}_red_as_magenta_max_projection.png"
    merged_path = outdir / f"{prefix}_merged_green_magenta_max_projection.png"

    Image.fromarray(green_max_2d).save(green_path)
    written.append(str(green_path))

    red_rgb = np.zeros((*red_max_2d.shape, 3), dtype=np.uint8)
    red_rgb[..., 0] = red_max_2d
    red_rgb[..., 2] = red_max_2d
    Image.fromarray(red_rgb).save(red_path)
    written.append(str(red_path))

    merged_rgb = np.zeros((*green_max_2d.shape, 3), dtype=np.uint8)
    merged_rgb[..., 0] = red_max_2d
    merged_rgb[..., 1] = green_max_2d
    merged_rgb[..., 2] = red_max_2d
    Image.fromarray(merged_rgb).save(merged_path)
    written.append(str(merged_path))
    return written


def public_stats(stats: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in stats.items() if not key.startswith("_")}


def main() -> int:
    args = parse_args()
    input_dir = args.input_dir.expanduser().resolve()
    output_dir = (args.output_dir or input_dir / "display_8bit").expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    green_path = input_dir / args.green_name
    red_path = input_dir / args.red_name
    if not green_path.exists():
        raise SystemExit(f"Missing green denoised input: {green_path}")
    if not red_path.exists():
        raise SystemExit(f"Missing red denoised input: {red_path}")

    print(f"Input directory: {input_dir}", flush=True)
    print(f"Output directory: {output_dir}", flush=True)
    print("Opening denoised OME-TIFFs read-only...", flush=True)
    green, green_axes, green_dtype = open_memmap(green_path)
    red, red_axes, red_dtype = open_memmap(red_path)
    if green.shape != red.shape:
        raise SystemExit(f"Green/red shapes differ: {green.shape} vs {red.shape}")
    if green_axes != red_axes:
        raise SystemExit(f"Green/red axes differ: {green_axes} vs {red_axes}")
    axes = green_axes
    print(f"Green: shape={green.shape} axes={green_axes} dtype={green_dtype}", flush=True)
    print(f"Red:   shape={red.shape} axes={red_axes} dtype={red_dtype}", flush=True)

    print("Calculating fixed whole-stack/channel display limits...", flush=True)
    green_lo, green_hi, green_limit_metadata = percentile_limits(
        green,
        axes,
        args.lower_percentile,
        args.upper_percentile,
        args.percentile_z_step,
        args.percentile_yx_step,
    )
    red_lo, red_hi, red_limit_metadata = percentile_limits(
        red,
        axes,
        args.lower_percentile,
        args.upper_percentile,
        args.percentile_z_step,
        args.percentile_yx_step,
    )
    print(f"Green limits: {green_lo:.6g}, {green_hi:.6g}", flush=True)
    print(f"Red limits:   {red_lo:.6g}, {red_hi:.6g}", flush=True)

    green_out = output_dir / "green_denoised_8bit_display.ome.tif"
    red_out = output_dir / "red_denoised_8bit_display.ome.tif"
    rgb_out = output_dir / "merged_green_magenta_8bit_display.ome.tif"

    print(f"Writing {green_out.name}", flush=True)
    green_stats = write_uint8_stack(
        green_out,
        green,
        axes,
        (green_lo, green_hi),
        args.compression,
        args.overwrite,
    )
    print(f"Writing {red_out.name}", flush=True)
    red_stats = write_uint8_stack(
        red_out,
        red,
        axes,
        (red_lo, red_hi),
        args.compression,
        args.overwrite,
    )

    rgb_stats: dict[str, Any] | None = None
    if not args.no_merged_rgb:
        print(f"Writing {rgb_out.name}", flush=True)
        rgb_stats = write_rgb_stack(
            rgb_out,
            green,
            red,
            axes,
            (green_lo, green_hi),
            (red_lo, red_hi),
            args.compression,
            args.overwrite,
        )

    png_paths = save_projection_pngs(
        output_dir,
        axes,
        green_stats["_max_projection_array"],
        red_stats["_max_projection_array"],
    )

    metadata = {
        "created_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "purpose": (
            "8-bit display/QC exports from CAREamics denoised outputs. "
            "Quantitative denoised OME-TIFFs remain unchanged."
        ),
        "input_dir": str(input_dir),
        "output_dir": str(output_dir),
        "input_files": {
            "green": str(green_path),
            "red": str(red_path),
        },
        "input": {
            "shape": list(map(int, green.shape)),
            "axes": axes,
            "green_dtype": green_dtype,
            "red_dtype": red_dtype,
        },
        "display_scaling": {
            "method": "fixed percentile scaling per channel from a deterministic whole-stack/movie sample",
            "green": green_limit_metadata,
            "red": red_limit_metadata,
            "uint8_mapping": "values <= low -> 0; values >= high -> 255",
        },
        "outputs": {
            "green_8bit": public_stats(green_stats),
            "red_8bit": public_stats(red_stats),
            "merged_rgb_8bit": public_stats(rgb_stats) if rgb_stats is not None else None,
            "max_projection_pngs": png_paths,
        },
    }
    metadata_path = output_dir / "display_scaling_8bit.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Wrote metadata: {metadata_path}", flush=True)
    print("Done. Original denoised OME-TIFFs were not modified.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
