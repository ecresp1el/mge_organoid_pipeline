#!/usr/bin/env python3
"""Denoise two-channel Imaris .ims confocal files with CAREamics N2V.

The raw .ims file is read only. Outputs are written as OME-TIFF previews,
model directories, and QC figures under smoke_test/ or full_run/.
"""

from __future__ import annotations

import argparse
import json
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np


CHANNELS = ((0, "green"), (1, "red"))


@dataclass
class ImsReadResult:
    data_tczyx: np.ndarray
    reader: str
    metadata: dict[str, Any]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Denoise green/red channels from an Imaris .ims file with CAREamics N2V."
    )
    parser.add_argument("--input", required=True, help="Input Imaris .ims file.")
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Base output directory. The run writes smoke_test/ or full_run/ below this directory.",
    )
    parser.add_argument(
        "--reader",
        choices=("auto", "h5py", "aicsimageio"),
        default="auto",
        help="Reader to use. h5py supports efficient .ims subset reads.",
    )
    parser.add_argument("--subset", action="store_true", help="Run a small smoke test.")
    parser.add_argument("--max-timepoints", type=int, default=2, help="Subset maximum timepoints.")
    parser.add_argument("--max-z", type=int, default=16, help="Subset maximum Z planes.")
    parser.add_argument("--crop-yx", type=int, default=512, help="Subset center crop size in Y and X.")
    parser.add_argument("--require-gpu", action="store_true", help="Exit if CUDA is unavailable.")
    parser.add_argument("--epochs", type=int, default=50, help="Training epochs per channel.")
    parser.add_argument("--batch-size", type=int, default=2, help="CAREamics training batch size.")
    parser.add_argument("--patch-size-z", type=int, default=16, help="N2V patch size in Z.")
    parser.add_argument("--patch-size-yx", type=int, default=64, help="N2V patch size in Y and X.")
    parser.add_argument(
        "--num-steps",
        type=int,
        default=None,
        help="Optional training steps per epoch. Defaults to 10 for smoke tests and CAREamics default for full runs.",
    )
    parser.add_argument("--seed", type=int, default=7, help="Random seed.")
    return parser.parse_args()


def print_header(title: str) -> None:
    print(f"\n=== {title} ===", flush=True)


def import_or_exit(module_name: str, install_hint: str | None = None) -> Any:
    try:
        return __import__(module_name)
    except ImportError as exc:
        hint = f" Install with: {install_hint}" if install_hint else ""
        raise SystemExit(f"Missing required Python module '{module_name}'.{hint}") from exc


def sanitize_hdf5_plugin_path() -> None:
    """Drop missing HDF5 plugin paths that can make h5py fail before reading data."""
    plugin_path = os.environ.get("HDF5_PLUGIN_PATH")
    fallback_plugin_dir = Path(os.environ.get("TMPDIR", "/tmp")) / "imaris_careamics_empty_hdf5_plugins"
    fallback_plugin_dir.mkdir(parents=True, exist_ok=True)
    if not plugin_path:
        os.environ["HDF5_PLUGIN_PATH"] = str(fallback_plugin_dir)
        return
    existing_paths = [part for part in plugin_path.split(os.pathsep) if part and Path(part).exists()]
    missing_paths = [part for part in plugin_path.split(os.pathsep) if part and not Path(part).exists()]
    if missing_paths:
        print(
            "WARNING: removing missing HDF5_PLUGIN_PATH entries: "
            + ", ".join(missing_paths),
            flush=True,
        )
    if existing_paths:
        os.environ["HDF5_PLUGIN_PATH"] = os.pathsep.join(existing_paths)
    else:
        os.environ["HDF5_PLUGIN_PATH"] = str(fallback_plugin_dir)


def sorted_numbered_keys(group: Any, prefix: str) -> list[str]:
    numbered: list[tuple[int, str]] = []
    for key in group.keys():
        match = re.fullmatch(rf"{re.escape(prefix)}\s*(\d+)", key)
        if match:
            numbered.append((int(match.group(1)), key))
    return [key for _, key in sorted(numbered)]


def numbered_key_map(group: Any, prefix: str) -> dict[int, str]:
    mapping: dict[int, str] = {}
    for key in group.keys():
        match = re.fullmatch(rf"{re.escape(prefix)}\s*(\d+)", key)
        if match:
            mapping[int(match.group(1))] = key
    return mapping


def center_crop_slice(size: int, requested: int | None) -> slice:
    if requested is None or requested <= 0 or requested >= size:
        return slice(0, size)
    start = max((size - requested) // 2, 0)
    return slice(start, start + requested)


def read_ims_h5py(
    path: Path,
    subset: bool,
    max_timepoints: int,
    max_z: int,
    crop_yx: int,
) -> ImsReadResult:
    sanitize_hdf5_plugin_path()
    try:
        import hdf5plugin  # noqa: F401
    except ImportError:
        print(
            "WARNING: hdf5plugin is not installed. Imaris .ims files using LZ4 or "
            "other external HDF5 filters may fail to read.",
            flush=True,
        )
    h5py = import_or_exit("h5py", "conda install h5py")

    with h5py.File(path, "r") as h5:
        if "DataSet" not in h5:
            raise ValueError("No /DataSet group found; this does not look like a standard .ims file.")
        dataset_group = h5["DataSet"]
        resolution_keys = sorted_numbered_keys(dataset_group, "ResolutionLevel")
        if not resolution_keys:
            raise ValueError("No /DataSet/ResolutionLevel N group found in .ims file.")
        resolution_key = resolution_keys[0]
        resolution_group = dataset_group[resolution_key]
        timepoint_keys = sorted_numbered_keys(resolution_group, "TimePoint")
        if not timepoint_keys:
            timepoint_keys = ["TimePoint 0"] if "TimePoint 0" in resolution_group else []
        if not timepoint_keys:
            raise ValueError("No TimePoint groups found in .ims file.")

        if subset:
            timepoint_keys = timepoint_keys[: max(1, max_timepoints)]

        first_tp = resolution_group[timepoint_keys[0]]
        channel_keys = sorted_numbered_keys(first_tp, "Channel")
        if len(channel_keys) < 2:
            raise ValueError(f"Expected at least 2 channels, found {len(channel_keys)}.")

        first_data = first_tp[channel_keys[0]]["Data"]
        if first_data.ndim != 3:
            raise ValueError(f"Expected ZYX channel data, found shape {first_data.shape}.")
        z_size, y_size, x_size = first_data.shape
        z_slice = slice(0, min(z_size, max_z)) if subset else slice(0, z_size)
        y_slice = center_crop_slice(y_size, crop_yx if subset else None)
        x_slice = center_crop_slice(x_size, crop_yx if subset else None)

        stacks_tczyx = []
        dataset_paths = []
        for timepoint_key in timepoint_keys:
            tp_group = resolution_group[timepoint_key]
            channel_key_by_index = numbered_key_map(tp_group, "Channel")
            stacks_czyx = []
            for channel_index, _channel_name in CHANNELS:
                channel_key = channel_key_by_index.get(channel_index)
                if channel_key is None:
                    raise ValueError(f"Missing Channel {channel_index} in {timepoint_key}.")
                data_path = f"/DataSet/{resolution_key}/{timepoint_key}/{channel_key}/Data"
                data = np.asarray(tp_group[channel_key]["Data"][z_slice, y_slice, x_slice])
                stacks_czyx.append(data)
                dataset_paths.append(data_path)
            stacks_tczyx.append(np.stack(stacks_czyx, axis=0))

        arr = np.stack(stacks_tczyx, axis=0)
        metadata = {
            "reader": "h5py",
            "resolution_level": resolution_key,
            "timepoints_read": len(timepoint_keys),
            "channels_read": [0, 1],
            "source_shape_zyx": [int(z_size), int(y_size), int(x_size)],
            "subset": bool(subset),
            "z_slice": [z_slice.start, z_slice.stop],
            "y_slice": [y_slice.start, y_slice.stop],
            "x_slice": [x_slice.start, x_slice.stop],
            "dataset_paths": dataset_paths,
        }
        return ImsReadResult(data_tczyx=arr, reader="h5py", metadata=metadata)


def read_ims_aicsimageio(path: Path, subset: bool, max_timepoints: int, max_z: int, crop_yx: int) -> ImsReadResult:
    try:
        from aicsimageio import AICSImage
    except ImportError as exc:
        raise SystemExit("Missing aicsimageio. Install with: conda install -c conda-forge aicsimageio") from exc

    img = AICSImage(path)
    arr = np.asarray(img.get_image_data("TCZYX", C=[0, 1]))
    if arr.ndim != 5:
        raise ValueError(f"AICSImage returned unexpected shape {arr.shape}; expected TCZYX.")
    if subset:
        arr = arr[: max(1, max_timepoints), :, : max(1, max_z), :, :]
        y_slice = center_crop_slice(arr.shape[-2], crop_yx)
        x_slice = center_crop_slice(arr.shape[-1], crop_yx)
        arr = arr[..., y_slice, x_slice]
    metadata = {
        "reader": "aicsimageio",
        "dims": str(img.dims),
        "subset": bool(subset),
        "note": "AICSImage may load more source data than h5py before applying subsets.",
    }
    return ImsReadResult(data_tczyx=arr, reader="aicsimageio", metadata=metadata)


def read_ims(path: Path, args: argparse.Namespace) -> ImsReadResult:
    if args.reader in ("auto", "h5py"):
        try:
            return read_ims_h5py(path, args.subset, args.max_timepoints, args.max_z, args.crop_yx)
        except Exception as exc:
            if args.reader == "h5py":
                raise
            print(f"WARNING: h5py .ims read failed; falling back to aicsimageio. Reason: {exc}", flush=True)
    return read_ims_aicsimageio(path, args.subset, args.max_timepoints, args.max_z, args.crop_yx)


def check_gpu(require_gpu: bool) -> bool:
    print_header("GPU check")
    torch = import_or_exit("torch", "install PyTorch with CUDA support")
    cuda_available = bool(torch.cuda.is_available())
    print(f"CUDA available: {cuda_available}", flush=True)
    if cuda_available:
        print(f"GPU name: {torch.cuda.get_device_name(0)}", flush=True)
        x = torch.ones((4, 4), device="cuda")
        y = (x @ x).sum().item()
        print(f"GPU tensor smoke test: {y:.1f}", flush=True)
    else:
        print("WARNING: CUDA is not available; CAREamics will run on CPU and may be very slow.", flush=True)
        if require_gpu:
            raise SystemExit("--require-gpu was passed, but PyTorch reports CUDA unavailable.")
    return cuda_available


def set_seeds(seed: int) -> None:
    np.random.seed(seed)
    try:
        import torch

        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
    except Exception:
        pass


def adjusted_patch_size(shape_zyx: tuple[int, int, int], requested_zyx: tuple[int, int, int]) -> tuple[int, int, int]:
    patch = []
    for dim, requested in zip(shape_zyx, requested_zyx):
        size = max(1, min(int(dim), int(requested)))
        if size >= 16:
            size = max(16, (size // 16) * 16)
        elif size >= 8:
            size = 8
        patch.append(size)
    return tuple(patch)  # type: ignore[return-value]


def create_n2v_config_for_channel(
    channel_name: str,
    patch_size: tuple[int, int, int],
    batch_size: int,
    epochs: int,
    num_steps: int | None,
    seed: int,
) -> Any:
    from careamics.config import create_advanced_n2v_config

    kwargs: dict[str, Any] = {
        "experiment_name": f"{channel_name}_n2v",
        "data_type": "array",
        "axes": "ZYX",
        "patch_size": patch_size,
        "batch_size": batch_size,
        "num_epochs": epochs,
        "seed": seed,
    }
    if num_steps is not None:
        kwargs["num_steps"] = num_steps
        kwargs["n_val_patches"] = 2
    return create_advanced_n2v_config(**kwargs)


def train_and_predict_channel(
    channel_name: str,
    raw_tzyx: np.ndarray,
    model_dir: Path,
    patch_size: tuple[int, int, int],
    batch_size: int,
    epochs: int,
    num_steps: int | None,
    seed: int,
) -> np.ndarray:
    from careamics import CAREamist

    model_dir.mkdir(parents=True, exist_ok=True)
    train_data = [np.asarray(raw_tzyx[t], dtype=np.float32) for t in range(raw_tzyx.shape[0])]
    config = create_n2v_config_for_channel(channel_name, patch_size, batch_size, epochs, num_steps, seed)

    print_header(f"Training {channel_name} N2V")
    print(f"Model directory: {model_dir}", flush=True)
    print(f"Training samples: {len(train_data)} timepoint(s)", flush=True)
    print(f"Patch size ZYX: {patch_size}", flush=True)
    print(f"Epochs: {epochs}", flush=True)
    print(f"Batch size: {batch_size}", flush=True)
    if num_steps is not None:
        print(f"Steps per epoch: {num_steps}", flush=True)

    careamist = CAREamist(config=config, work_dir=model_dir)
    careamist.train(train_data=train_data)

    print_header(f"Predicting {channel_name}")
    predictions = []
    tile_overlap = tuple(max(1, min(8, size // 4)) for size in patch_size)
    for t in range(raw_tzyx.shape[0]):
        pred = careamist.predict(
            pred_data=[np.asarray(raw_tzyx[t], dtype=np.float32)],
            tile_size=patch_size,
            tile_overlap=tile_overlap,
            batch_size=1,
        )
        if isinstance(pred, tuple):
            pred = pred[0]
        if isinstance(pred, list):
            pred_arr = np.asarray(pred[0])
        else:
            pred_arr = np.asarray(pred)
        pred_arr = np.squeeze(pred_arr).astype(np.float32, copy=False)
        if pred_arr.shape != raw_tzyx[t].shape:
            raise RuntimeError(
                f"{channel_name} prediction for T={t} has shape {pred_arr.shape}, "
                f"expected {raw_tzyx[t].shape}."
            )
        predictions.append(pred_arr)
    return np.stack(predictions, axis=0)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def save_ome_tiff(path: Path, arr: np.ndarray, axes: str) -> None:
    import tifffile

    path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(path, arr, ome=True, metadata={"axes": axes})


def percentile_scale_uint16(arr: np.ndarray, low: float = 0.5, high: float = 99.8) -> np.ndarray:
    arr_float = np.asarray(arr, dtype=np.float32)
    p_low, p_high = np.percentile(arr_float, [low, high])
    return scale_to_uint16(arr_float, p_low, p_high)


def scale_to_uint16(arr: np.ndarray, p_low: float, p_high: float) -> np.ndarray:
    arr_float = np.asarray(arr, dtype=np.float32)
    if not np.isfinite(p_low) or not np.isfinite(p_high) or p_high <= p_low:
        return np.zeros(arr.shape, dtype=np.uint16)
    scaled = np.clip((arr_float - p_low) / (p_high - p_low), 0, 1)
    return np.round(scaled * 65535).astype(np.uint16)


def percentile_scale_uint8(arr: np.ndarray, low: float = 0.5, high: float = 99.8) -> np.ndarray:
    scaled16 = percentile_scale_uint16(arr, low=low, high=high)
    return (scaled16 / 257).astype(np.uint8)


def scale_to_uint8(arr: np.ndarray, p_low: float, p_high: float) -> np.ndarray:
    return (scale_to_uint16(arr, p_low, p_high) / 257).astype(np.uint8)


def save_merged_preview(outdir: Path, green_tzyx: np.ndarray, red_tzyx: np.ndarray) -> None:
    green_scaled = percentile_scale_uint16(green_tzyx)
    red_scaled = percentile_scale_uint16(red_tzyx)
    merged_tczyx = np.stack([green_scaled, red_scaled], axis=1)
    save_ome_tiff(outdir / "merged_green_red_preview.ome.tif", merged_tczyx, axes="TCZYX")


def save_max_projection_mp4(outdir: Path, green_tzyx: np.ndarray, red_tzyx: np.ndarray) -> None:
    imageio = import_or_exit("imageio", "pip install imageio imageio-ffmpeg")

    green_u8 = percentile_scale_uint8(green_tzyx)
    red_u8 = percentile_scale_uint8(red_tzyx)
    frames = []
    for t in range(green_tzyx.shape[0]):
        rgb = np.zeros((*green_tzyx.shape[-2:], 3), dtype=np.uint8)
        rgb[..., 1] = np.max(green_u8[t], axis=0)
        rgb[..., 0] = np.max(red_u8[t], axis=0)
        frames.append(rgb)
    imageio.mimsave(outdir / "max_projection_preview.mp4", frames, fps=4, quality=8)


def write_saturation_report(outdir: Path, arrays: dict[str, np.ndarray]) -> None:
    lines = ["Clipping/saturation report", ""]
    warning = False
    for name, arr in arrays.items():
        finite = np.isfinite(arr)
        n = int(arr.size)
        n_nonfinite = int(n - np.count_nonzero(finite))
        lines.append(f"{name}:")
        lines.append(f"  dtype: {arr.dtype}")
        lines.append(f"  shape: {list(arr.shape)}")
        lines.append(f"  min: {float(np.nanmin(arr))}")
        lines.append(f"  max: {float(np.nanmax(arr))}")
        lines.append(f"  nonfinite_pixels: {n_nonfinite}")
        if np.issubdtype(arr.dtype, np.integer):
            dtype_info = np.iinfo(arr.dtype)
            n_min = int(np.count_nonzero(arr == dtype_info.min))
            n_max = int(np.count_nonzero(arr == dtype_info.max))
            pct_min = 100.0 * n_min / max(n, 1)
            pct_max = 100.0 * n_max / max(n, 1)
            lines.append(f"  pixels_at_dtype_min: {n_min} ({pct_min:.5f}%)")
            lines.append(f"  pixels_at_dtype_max: {n_max} ({pct_max:.5f}%)")
            if pct_max >= 0.1:
                lines.append("  WARNING: >=0.1% pixels are at dtype maximum.")
                warning = True
        if n_nonfinite > 0:
            lines.append("  WARNING: non-finite values found.")
            warning = True
        lines.append("")
    if not warning:
        lines.append("No obvious clipping/saturation warnings triggered.")
    (outdir / "qc" / "clipping_saturation_warnings.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def save_qc_figures(
    outdir: Path,
    raw_channels: dict[str, np.ndarray],
    denoised_channels: dict[str, np.ndarray],
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    qc_dir = outdir / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)

    names = ["green", "red"]
    t0 = 0
    z_mid = raw_channels["green"].shape[1] // 2

    fig, axes = plt.subplots(2, 2, figsize=(9, 8), constrained_layout=True)
    for row, name in enumerate(names):
        limits = np.percentile(
            np.concatenate(
                [
                    np.asarray(raw_channels[name], dtype=np.float32).ravel(),
                    np.asarray(denoised_channels[name], dtype=np.float32).ravel(),
                ]
            ),
            [0.5, 99.8],
        )
        raw = scale_to_uint8(raw_channels[name][t0, z_mid], float(limits[0]), float(limits[1]))
        den = scale_to_uint8(denoised_channels[name][t0, z_mid], float(limits[0]), float(limits[1]))
        axes[row, 0].imshow(raw, cmap="gray")
        axes[row, 0].set_title(f"{name} raw T0 Z{z_mid}")
        axes[row, 1].imshow(den, cmap="gray")
        axes[row, 1].set_title(f"{name} denoised T0 Z{z_mid}")
        for col in range(2):
            axes[row, col].axis("off")
    fig.savefig(qc_dir / "raw_vs_denoised_example_slices.png", dpi=200)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(9, 8), constrained_layout=True)
    for row, name in enumerate(names):
        limits = np.percentile(
            np.concatenate(
                [
                    np.asarray(raw_channels[name], dtype=np.float32).ravel(),
                    np.asarray(denoised_channels[name], dtype=np.float32).ravel(),
                ]
            ),
            [0.5, 99.8],
        )
        raw_max = np.max(raw_channels[name][t0], axis=0)
        den_max = np.max(denoised_channels[name][t0], axis=0)
        axes[row, 0].imshow(scale_to_uint8(raw_max, float(limits[0]), float(limits[1])), cmap="gray")
        axes[row, 0].set_title(f"{name} raw max projection")
        axes[row, 1].imshow(scale_to_uint8(den_max, float(limits[0]), float(limits[1])), cmap="gray")
        axes[row, 1].set_title(f"{name} denoised max projection")
        for col in range(2):
            axes[row, col].axis("off")
    fig.savefig(qc_dir / "raw_vs_denoised_max_projections.png", dpi=200)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    for col, name in enumerate(names):
        raw = raw_channels[name].ravel()
        den = denoised_channels[name].ravel()
        axes[col].hist(raw[np.isfinite(raw)], bins=128, alpha=0.55, label="raw", density=True)
        axes[col].hist(den[np.isfinite(den)], bins=128, alpha=0.55, label="denoised", density=True)
        axes[col].set_title(f"{name} intensity histogram")
        axes[col].set_xlabel("intensity")
        axes[col].set_ylabel("density")
        axes[col].legend()
    fig.savefig(qc_dir / "intensity_histograms.png", dpi=200)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    input_path = Path(args.input).expanduser().resolve()
    if not input_path.exists():
        raise SystemExit(f"Input file does not exist: {input_path}")

    run_mode = "smoke_test" if args.subset else "full_run"
    outdir = Path(args.output_dir).expanduser().resolve() / run_mode
    model_root = outdir / "models"
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "qc").mkdir(exist_ok=True)
    model_root.mkdir(exist_ok=True)

    epochs = args.epochs
    num_steps = args.num_steps
    if args.subset:
        epochs = min(epochs, 5)
        if num_steps is None:
            num_steps = 10

    print_header("Run setup")
    print(f"Input .ims: {input_path}", flush=True)
    print(f"Output directory: {outdir}", flush=True)
    print(f"Run mode: {run_mode}", flush=True)
    print("Raw .ims policy: read only; file is not modified.", flush=True)

    set_seeds(args.seed)
    cuda_available = check_gpu(args.require_gpu)

    print_header("Reading .ims")
    ims = read_ims(input_path, args)
    data = ims.data_tczyx
    if data.ndim != 5 or data.shape[1] < 2:
        raise RuntimeError(f"Internal data must be TCZYX with two channels; got {data.shape}.")
    print(f"Reader: {ims.reader}", flush=True)
    print(f"Internal shape TCZYX: {data.shape}", flush=True)
    print(f"Input dtype: {data.dtype}", flush=True)

    write_json(
        outdir / "run_metadata.json",
        {
            "input_path": str(input_path),
            "output_dir": str(outdir),
            "run_mode": run_mode,
            "cuda_available": cuda_available,
            "args": vars(args),
            "ims_metadata": ims.metadata,
            "internal_shape_tczyx": list(map(int, data.shape)),
            "dtype": str(data.dtype),
        },
    )

    raw_channels = {
        "green": np.asarray(data[:, 0, :, :, :]),
        "red": np.asarray(data[:, 1, :, :, :]),
    }

    print_header("Saving raw channel OME-TIFFs")
    save_ome_tiff(outdir / "green_raw.ome.tif", raw_channels["green"], axes="TZYX")
    save_ome_tiff(outdir / "red_raw.ome.tif", raw_channels["red"], axes="TZYX")

    requested_patch = (args.patch_size_z, args.patch_size_yx, args.patch_size_yx)
    patch_size = adjusted_patch_size(tuple(map(int, raw_channels["green"].shape[1:])), requested_patch)
    if patch_size != requested_patch:
        print(f"WARNING: adjusted patch size from {requested_patch} to {patch_size} for data shape.", flush=True)

    denoised_channels: dict[str, np.ndarray] = {}
    for _channel_index, name in CHANNELS:
        denoised_channels[name] = train_and_predict_channel(
            channel_name=name,
            raw_tzyx=raw_channels[name],
            model_dir=model_root / f"{name}_n2v",
            patch_size=patch_size,
            batch_size=args.batch_size,
            epochs=epochs,
            num_steps=num_steps,
            seed=args.seed,
        )

    print_header("Saving denoised OME-TIFFs and previews")
    save_ome_tiff(outdir / "green_denoised.ome.tif", denoised_channels["green"], axes="TZYX")
    save_ome_tiff(outdir / "red_denoised.ome.tif", denoised_channels["red"], axes="TZYX")
    save_merged_preview(outdir, denoised_channels["green"], denoised_channels["red"])
    save_max_projection_mp4(outdir, denoised_channels["green"], denoised_channels["red"])

    print_header("Saving QC")
    save_qc_figures(outdir, raw_channels, denoised_channels)
    write_saturation_report(
        outdir,
        {
            "green_raw": raw_channels["green"],
            "red_raw": raw_channels["red"],
            "green_denoised": denoised_channels["green"],
            "red_denoised": denoised_channels["red"],
        },
    )

    print_header("Done")
    print(f"Outputs written to: {outdir}", flush=True)
    print(f"Green model: {model_root / 'green_n2v'}", flush=True)
    print(f"Red model: {model_root / 'red_n2v'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
