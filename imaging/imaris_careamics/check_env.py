#!/usr/bin/env python3
"""Print environment details for the Imaris CAREamics pipeline."""

from __future__ import annotations

import importlib.metadata
import platform


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not installed"


def main() -> int:
    print(f"Python version: {platform.python_version()}")
    print(f"CAREamics version: {package_version('careamics')}")
    print(f"PyTorch version: {package_version('torch')}")

    try:
        import torch
    except ImportError:
        print("CUDA available: false")
        print("GPU tensor smoke test: skipped because torch is not installed")
        return 1

    cuda_available = bool(torch.cuda.is_available())
    print(f"CUDA available: {str(cuda_available).lower()}")
    if cuda_available:
        print(f"GPU name: {torch.cuda.get_device_name(0)}")
        x = torch.ones((4, 4), device="cuda")
        y = (x @ x).sum().item()
        print(f"GPU tensor smoke test: {y:.1f}")
    else:
        print("GPU name: none")
        print("GPU tensor smoke test: skipped because CUDA is unavailable")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
