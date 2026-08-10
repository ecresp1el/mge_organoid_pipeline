#!/usr/bin/env python3
"""Assemble marker-context figures above the receptor/MACS framework."""

from __future__ import annotations

import argparse
from datetime import datetime
import hashlib
from pathlib import Path
import shutil

import pandas as pd
from PIL import Image, ImageChops, ImageDraw, ImageFont


PROJECT_ROOT_DEFAULT = Path("/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder")
CONTEXT_ATLAS_RELATIVE = Path(
    "final_figures/fig_cross_study_marker_expression_pv_precursors_on_off_target_v12_candidate/"
    "figures/png/cross_study_marker_expression_pv_precursors_on_off_target.png"
)
CONTEXT_MATRIX_RELATIVE = Path(
    "final_figures/fig_cross_study_integrated_mge_marker_expression_v7_overlay_zero_baseline_candidate/"
    "figures/png/cross_study_canonical_mge_marker_expression_violin_positive_matrix.png"
)
FINAL_FOLDER_DEFAULT = "fig_cross_study_receptor_macs_enrichment_v1_candidate"
MACS_COMPONENT_NAME = "cross_study_receptor_macs_enrichment_composite_component.png"
OUTPUT_STEM = "cross_study_marker_context_and_receptor_macs_enrichment"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT_DEFAULT)
    parser.add_argument("--context-atlas", type=Path, default=None)
    parser.add_argument("--context-matrix", type=Path, default=None)
    parser.add_argument("--macs-component", type=Path, default=None)
    parser.add_argument("--final-dir", type=Path, default=None)
    parser.add_argument("--target-content-width", type=int, default=5500)
    parser.add_argument("--dpi", type=int, default=300)
    return parser.parse_args()


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    name = "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"
    candidates = [
        Path("/usr/share/fonts/dejavu-sans-fonts") / name,
        Path("/usr/share/fonts/dejavu") / name,
        Path("/usr/share/fonts/truetype/dejavu") / name,
    ]
    for path in candidates:
        if path.exists():
            return ImageFont.truetype(str(path), size=size)
    return ImageFont.load_default()


def trim_white(image: Image.Image, padding: int = 30) -> Image.Image:
    rgb = image.convert("RGB")
    background = Image.new("RGB", rgb.size, "white")
    difference = ImageChops.difference(rgb, background).convert("L")
    bbox = difference.point(lambda value: 255 if value > 8 else 0).getbbox()
    if bbox is None:
        return rgb
    left, top, right, bottom = bbox
    left = max(0, left - padding)
    top = max(0, top - padding)
    right = min(rgb.width, right + padding)
    bottom = min(rgb.height, bottom + padding)
    return rgb.crop((left, top, right, bottom))


def resize_width(image: Image.Image, width: int) -> Image.Image:
    height = int(round(image.height * width / image.width))
    return image.resize((width, height), Image.Resampling.LANCZOS)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    args = parse_args()
    project_root = args.project_root.resolve()
    final_dir = (args.final_dir or project_root / "final_figures" / FINAL_FOLDER_DEFAULT).resolve()
    atlas_path = (args.context_atlas or project_root / CONTEXT_ATLAS_RELATIVE).resolve()
    matrix_path = (args.context_matrix or project_root / CONTEXT_MATRIX_RELATIVE).resolve()
    macs_path = (
        args.macs_component
        or final_dir / "figures" / "png" / MACS_COMPONENT_NAME
    ).resolve()
    for path in [atlas_path, matrix_path, macs_path]:
        if not path.exists():
            raise FileNotFoundError(path)

    images = [resize_width(trim_white(Image.open(path)), args.target_content_width) for path in [atlas_path, matrix_path, macs_path]]
    left_margin = 155
    right_margin = 55
    header_height = 260
    gap = 105
    footer_height = 80
    canvas_width = left_margin + args.target_content_width + right_margin
    canvas_height = header_height + sum(image.height for image in images) + gap * 2 + footer_height
    canvas = Image.new("RGB", (canvas_width, canvas_height), "white")
    draw = ImageDraw.Draw(canvas)
    draw.text(
        (left_margin, 38),
        "Cross-study marker context and receptor-guided enrichment",
        fill="#17191B",
        font=font(88, bold=True),
    )
    draw.text(
        (left_margin, 148),
        "From on/off-target identity and sample-level expression to an RNA-proxy MACS composition model",
        fill="#555B61",
        font=font(38),
    )

    y = header_height
    source_rows = []
    for label, path, image in zip(["a", "b", ""], [atlas_path, matrix_path, macs_path], images, strict=True):
        if label:
            draw.text((42, y + 24), label, fill="#111416", font=font(78, bold=True))
        canvas.paste(image, (left_margin, y))
        source_rows.append(
            {
                "section": label or "c-e",
                "source_path": str(path),
                "source_sha256": sha256(path),
                "source_width_px": Image.open(path).width,
                "source_height_px": Image.open(path).height,
                "composite_width_px": image.width,
                "composite_height_px": image.height,
            }
        )
        y += image.height
        if path != macs_path:
            draw.line((left_margin, y + gap // 2, canvas_width - right_margin, y + gap // 2), fill="#D9DDE1", width=3)
            y += gap

    png_path = final_dir / "figures" / "png" / f"{OUTPUT_STEM}.png"
    pdf_path = final_dir / "figures" / "pdf" / f"{OUTPUT_STEM}.pdf"
    png_path.parent.mkdir(parents=True, exist_ok=True)
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(png_path, dpi=(args.dpi, args.dpi), optimize=True)
    canvas.save(pdf_path, "PDF", resolution=float(args.dpi))

    source_table = final_dir / "provenance" / "composite_source_manifest.tsv"
    pd.DataFrame(source_rows).to_csv(source_table, sep="\t", index=False)
    render_table = final_dir / "provenance" / "composite_render_manifest.tsv"
    pd.DataFrame(
        [
            ("rendered_at", datetime.now().astimezone().isoformat()),
            ("output_png", str(png_path)),
            ("output_pdf", str(pdf_path)),
            ("canvas_width_px", str(canvas.width)),
            ("canvas_height_px", str(canvas.height)),
            ("dpi", str(args.dpi)),
            ("layout", "full-width vertical context atlas; expression matrix; MACS panels c-e"),
        ],
        columns=["key", "value"],
    ).to_csv(render_table, sep="\t", index=False)
    shutil.copy2(Path(__file__).resolve(), final_dir / "code" / Path(__file__).name)

    checksum_path = final_dir / "provenance" / "sha256_manifest.txt"
    checksum_files = [
        *sorted((final_dir / "figures").glob("*/*")),
        *sorted((final_dir / "tables").glob("*")),
        *sorted((final_dir / "code").glob("*")),
        *sorted((final_dir / "provenance").glob("*")),
        final_dir / "README.md",
    ]
    checksum_path.write_text(
        "".join(
            f"{sha256(path)}  {path.relative_to(final_dir)}\n"
            for path in checksum_files
            if path.is_file() and path != checksum_path
        )
    )
    print(f"output_png={png_path}")
    print(f"output_pdf={pdf_path}")
    print(f"canvas={canvas.width}x{canvas.height}")


if __name__ == "__main__":
    main()
