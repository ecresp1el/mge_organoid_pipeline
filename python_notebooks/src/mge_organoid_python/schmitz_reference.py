"""Fetch Schmitz 2022 Cell Browser metadata and UMAP coordinates."""

from __future__ import annotations

import csv
import gzip
import json
import shutil
import tempfile
import urllib.request
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

from .paths import resolve_project_root


BASE_URL = "https://dev-inhibitory-neurons.cells.ucsc.edu/dev-inhibitory-neurons"
RUN_LABEL_DEFAULT = "schmitz_2022_reference_metadata_umap_v1"
RESULTS_DIRNAME = "schmitz_2022_reference_metadata_umap"


@dataclass(frozen=True)
class SchmitzDataset:
    dataset_id: str
    ucsc_path: str
    label: str
    species: str


DATASETS = (
    SchmitzDataset("schmitz_macaque_dev", "macaque", "Macaque developing inhibitory neurons", "Rhesus macaque"),
    SchmitzDataset("schmitz_mouse_dev", "mouse-dev", "Mouse developing inhibitory neurons", "Mouse"),
    SchmitzDataset("schmitz_mouse_adult", "mouse-adult", "Mouse adult inhibitory neurons", "Mouse"),
)


def run_dir(project_root: str | Path | None = None, run_label: str = RUN_LABEL_DEFAULT) -> Path:
    return resolve_project_root(project_root) / "results" / RESULTS_DIRNAME / run_label


def output_paths(project_root: str | Path | None = None, run_label: str = RUN_LABEL_DEFAULT) -> dict[str, Path]:
    root = run_dir(project_root, run_label)
    return {
        "run": root,
        "config": root / "config",
        "metadata": root / "metadata",
        "coords": root / "coords",
        "tables": root / "tables",
    }


def log(message: str) -> None:
    print(f"[Schmitz2022] {message}", flush=True)


def dataset_url(dataset: SchmitzDataset, filename: str) -> str:
    return f"{BASE_URL}/{dataset.ucsc_path}/{filename}"


def download(url: str, dest: Path, skip_existing: bool = True) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if skip_existing and dest.exists() and dest.stat().st_size > 0:
        return
    tmp_handle = tempfile.NamedTemporaryFile(delete=False, dir=str(dest.parent))
    tmp_path = Path(tmp_handle.name)
    try:
        with tmp_handle:
            with urllib.request.urlopen(url) as response:
                shutil.copyfileobj(response, tmp_handle)
        tmp_path.replace(dest)
    except Exception:
        if tmp_path.exists():
            tmp_path.unlink()
        raise


def stream_url_to_gzip(url: str, dest: Path, skip_existing: bool = True) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if skip_existing and dest.exists() and dest.stat().st_size > 0:
        return
    tmp_handle = tempfile.NamedTemporaryFile(delete=False, dir=str(dest.parent))
    tmp_path = Path(tmp_handle.name)
    try:
        with tmp_handle:
            with urllib.request.urlopen(url) as response:
                with gzip.GzipFile(fileobj=tmp_handle, mode="wb") as gz:
                    shutil.copyfileobj(response, gz)
        tmp_path.replace(dest)
    except Exception:
        if tmp_path.exists():
            tmp_path.unlink()
        raise


def read_json(path: Path) -> dict:
    with path.open() as handle:
        return json.load(handle)


def write_rows(path: Path, rows: Iterable[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def meta_path(paths: dict[str, Path], dataset: SchmitzDataset) -> Path:
    return paths["metadata"] / f"{dataset.dataset_id}.meta.tsv.gz"


def dataset_json_path(paths: dict[str, Path], dataset: SchmitzDataset) -> Path:
    return paths["config"] / f"{dataset.dataset_id}.dataset.json"


def desc_json_path(paths: dict[str, Path], dataset: SchmitzDataset) -> Path:
    return paths["config"] / f"{dataset.dataset_id}.desc.json"


def coord_bin_path(paths: dict[str, Path], dataset: SchmitzDataset, coord_name: str) -> Path:
    return paths["coords"] / dataset.dataset_id / f"{coord_name}.coords.bin"


def coord_tsv_path(paths: dict[str, Path], dataset: SchmitzDataset, coord_name: str) -> Path:
    return paths["coords"] / dataset.dataset_id / f"{coord_name}.coords.tsv.gz"


def metadata_umap_path(paths: dict[str, Path], dataset: SchmitzDataset) -> Path:
    return paths["metadata"] / f"{dataset.dataset_id}.metadata_with_umap.tsv.gz"


def iter_meta_rows(path: Path):
    with gzip.open(str(path), "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            yield row


def read_cell_ids(path: Path) -> list[str]:
    return [row["cellId"] for row in iter_meta_rows(path)]


def decode_coord_bin(bin_path: Path, cell_ids: list[str], coord_info: dict, out_path: Path) -> int:
    raw = np.fromfile(str(bin_path), dtype="<u2")
    expected = len(cell_ids) * 2
    if raw.size != expected:
        raise ValueError(f"{bin_path} has {raw.size} uint16 values; expected {expected}")
    coords = raw.reshape((len(cell_ids), 2))
    max_x = float(coord_info.get("maxX", 65535) or 65535)
    max_y = float(coord_info.get("maxY", 65535) or 65535)
    min_x = float(coord_info.get("minX", 0) or 0)
    min_y = float(coord_info.get("minY", 0) or 0)
    denom_x = max(max_x - min_x, 1.0)
    denom_y = max(max_y - min_y, 1.0)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(str(out_path), "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["cell_id", "umap_1_u16", "umap_2_u16", "umap_1_norm", "umap_2_norm"])
        for cell_id, (x, y) in zip(cell_ids, coords):
            x_norm = (float(x) - min_x) / denom_x
            y_norm = (float(y) - min_y) / denom_y
            writer.writerow([cell_id, int(x), int(y), f"{x_norm:.8f}", f"{y_norm:.8f}"])
    return len(cell_ids)


def write_metadata_with_umap(paths: dict[str, Path], dataset: SchmitzDataset, dataset_json: dict) -> int:
    source_meta = meta_path(paths, dataset)
    out_path = metadata_umap_path(paths, dataset)
    coord_names = [coord["name"] for coord in dataset_json.get("coords", [])]
    coord_arrays = {}
    for coord in dataset_json.get("coords", []):
        raw = np.fromfile(str(coord_bin_path(paths, dataset, coord["name"])), dtype="<u2").reshape((-1, 2))
        max_x = float(coord.get("maxX", 65535) or 65535)
        max_y = float(coord.get("maxY", 65535) or 65535)
        min_x = float(coord.get("minX", 0) or 0)
        min_y = float(coord.get("minY", 0) or 0)
        coord_arrays[coord["name"]] = (
            raw,
            min_x,
            max(max_x - min_x, 1.0),
            min_y,
            max(max_y - min_y, 1.0),
        )

    n_rows = 0
    with gzip.open(str(source_meta), "rt") as in_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        for coord_name in coord_names:
            fieldnames.extend(
                [
                    f"{coord_name}_umap_1_u16",
                    f"{coord_name}_umap_2_u16",
                    f"{coord_name}_umap_1_norm",
                    f"{coord_name}_umap_2_norm",
                ]
            )
        with gzip.open(str(out_path), "wt", newline="") as out_handle:
            writer = csv.DictWriter(out_handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for idx, row in enumerate(reader):
                for coord_name in coord_names:
                    raw, min_x, denom_x, min_y, denom_y = coord_arrays[coord_name]
                    x, y = raw[idx]
                    row[f"{coord_name}_umap_1_u16"] = int(x)
                    row[f"{coord_name}_umap_2_u16"] = int(y)
                    row[f"{coord_name}_umap_1_norm"] = f"{(float(x) - min_x) / denom_x:.8f}"
                    row[f"{coord_name}_umap_2_norm"] = f"{(float(y) - min_y) / denom_y:.8f}"
                writer.writerow(row)
                n_rows += 1
    return n_rows


def write_dataset_table(paths: dict[str, Path], datasets_by_id: dict[str, dict]) -> None:
    rows = []
    for dataset in DATASETS:
        conf = datasets_by_id[dataset.dataset_id]
        rows.append(
            {
                "dataset_id": dataset.dataset_id,
                "study_label": dataset.label,
                "species": dataset.species,
                "ucsc_path": dataset.ucsc_path,
                "sample_count": conf.get("sampleCount", ""),
                "label_field": conf.get("labelField", ""),
                "default_color_field": conf.get("defColorField", ""),
                "coords": ",".join(coord["name"] for coord in conf.get("coords", [])),
                "metadata_path": str(meta_path(paths, dataset)),
                "metadata_with_umap_path": str(metadata_umap_path(paths, dataset)),
                "dataset_json_path": str(dataset_json_path(paths, dataset)),
                "desc_json_path": str(desc_json_path(paths, dataset)),
            }
        )
    write_rows(
        paths["tables"] / "schmitz_2022_datasets.tsv",
        rows,
        [
            "dataset_id",
            "study_label",
            "species",
            "ucsc_path",
            "sample_count",
            "label_field",
            "default_color_field",
            "coords",
            "metadata_path",
            "metadata_with_umap_path",
            "dataset_json_path",
            "desc_json_path",
        ],
    )


def write_metadata_fields(paths: dict[str, Path], datasets_by_id: dict[str, dict]) -> None:
    rows = []
    for dataset in DATASETS:
        conf = datasets_by_id[dataset.dataset_id]
        for field in conf.get("metaFields", []):
            rows.append(
                {
                    "dataset_id": dataset.dataset_id,
                    "field_name": field.get("name", ""),
                    "field_label": field.get("label", ""),
                    "field_type": field.get("type", ""),
                    "array_type": field.get("arrType", ""),
                    "value_count_entries": len(field.get("valCounts") or []),
                }
            )
    write_rows(
        paths["tables"] / "schmitz_2022_metadata_fields.tsv",
        rows,
        ["dataset_id", "field_name", "field_label", "field_type", "array_type", "value_count_entries"],
    )


def write_metadata_value_counts(paths: dict[str, Path]) -> None:
    rows = []
    for dataset in DATASETS:
        counters: dict[str, Counter] = defaultdict(Counter)
        total = 0
        for row in iter_meta_rows(meta_path(paths, dataset)):
            total += 1
            for field, value in row.items():
                if field == "cellId":
                    continue
                counters[field][value] += 1
        for field, counter in counters.items():
            for value, count in counter.most_common():
                rows.append(
                    {
                        "dataset_id": dataset.dataset_id,
                        "field_name": field,
                        "value": value,
                        "count": count,
                        "pct_of_dataset": f"{100.0 * count / total:.4f}",
                    }
                )
    write_rows(
        paths["tables"] / "schmitz_2022_metadata_value_counts.tsv",
        rows,
        ["dataset_id", "field_name", "value", "count", "pct_of_dataset"],
    )


def write_class_composition(paths: dict[str, Path]) -> None:
    rows = []
    group_fields = ["region", "timepoint", "batch_name", "dataset_name", "phase", "leiden", "hires_leiden"]
    for dataset in DATASETS:
        totals = Counter()
        counts: dict[tuple[str, str], Counter] = defaultdict(Counter)
        for row in iter_meta_rows(meta_path(paths, dataset)):
            label = row.get("class", "")
            totals[label] += 1
            for field in group_fields:
                if field in row:
                    counts[(label, field)][row[field]] += 1
        for label, total in totals.most_common():
            for field in group_fields:
                if (label, field) not in counts:
                    continue
                for value, count in counts[(label, field)].most_common():
                    rows.append(
                        {
                            "dataset_id": dataset.dataset_id,
                            "class": label,
                            "group_field": field,
                            "group_value": value,
                            "count": count,
                            "pct_of_class": f"{100.0 * count / total:.4f}",
                        }
                    )
    write_rows(
        paths["tables"] / "schmitz_2022_class_composition.tsv",
        rows,
        ["dataset_id", "class", "group_field", "group_value", "count", "pct_of_class"],
    )


def write_reference_label_summary(paths: dict[str, Path]) -> None:
    rows = []
    for dataset in DATASETS:
        class_counter = Counter()
        region_counter = Counter()
        mge_class_counter = Counter()
        total = 0
        mge_like = 0
        for row in iter_meta_rows(meta_path(paths, dataset)):
            total += 1
            label = row.get("class", "")
            region = row.get("region", "")
            class_counter[label] += 1
            region_counter[region] += 1
            if label.startswith("MGE_") or label.startswith("Ctx_LHX6") or label.startswith("Ctx_SST"):
                mge_like += 1
                mge_class_counter[label] += 1
        rows.append(
            {
                "dataset_id": dataset.dataset_id,
                "summary_type": "TOTAL",
                "label": "all_cells",
                "count": total,
                "pct_of_dataset": "100.0000",
            }
        )
        rows.append(
            {
                "dataset_id": dataset.dataset_id,
                "summary_type": "MGE_LIKE_BY_CLASS_PREFIX",
                "label": "MGE_/Ctx_LHX6/Ctx_SST",
                "count": mge_like,
                "pct_of_dataset": f"{100.0 * mge_like / total:.4f}",
            }
        )
        for label, count in class_counter.most_common():
            rows.append(
                {
                    "dataset_id": dataset.dataset_id,
                    "summary_type": "CLASS",
                    "label": label,
                    "count": count,
                    "pct_of_dataset": f"{100.0 * count / total:.4f}",
                }
            )
        for label, count in mge_class_counter.most_common():
            rows.append(
                {
                    "dataset_id": dataset.dataset_id,
                    "summary_type": "MGE_LIKE_CLASS",
                    "label": label,
                    "count": count,
                    "pct_of_dataset": f"{100.0 * count / total:.4f}",
                }
            )
        for label, count in region_counter.most_common():
            rows.append(
                {
                    "dataset_id": dataset.dataset_id,
                    "summary_type": "REGION",
                    "label": label,
                    "count": count,
                    "pct_of_dataset": f"{100.0 * count / total:.4f}",
                }
            )
    write_rows(
        paths["tables"] / "schmitz_2022_reference_label_summary.tsv",
        rows,
        ["dataset_id", "summary_type", "label", "count", "pct_of_dataset"],
    )


def write_coord_summary(paths: dict[str, Path], datasets_by_id: dict[str, dict], coord_rows_by_key: dict[tuple[str, str], int]) -> None:
    rows = []
    for dataset in DATASETS:
        for coord in datasets_by_id[dataset.dataset_id].get("coords", []):
            coord_name = coord["name"]
            rows.append(
                {
                    "dataset_id": dataset.dataset_id,
                    "coord_name": coord_name,
                    "coord_label": coord.get("shortLabel", ""),
                    "coord_bin_path": str(coord_bin_path(paths, dataset, coord_name)),
                    "coord_tsv_path": str(coord_tsv_path(paths, dataset, coord_name)),
                    "n_coordinate_rows": coord_rows_by_key[(dataset.dataset_id, coord_name)],
                    "coord_storage_type": coord.get("type", ""),
                    "source_text_filename": coord.get("textFname", ""),
                    "notes": "Decoded from UCSC Cell Browser UInt16 coords.bin; normalized coordinates are browser-scale 0..1, not original Scanpy UMAP units.",
                }
            )
    write_rows(
        paths["tables"] / "schmitz_2022_umap_files.tsv",
        rows,
        [
            "dataset_id",
            "coord_name",
            "coord_label",
            "coord_bin_path",
            "coord_tsv_path",
            "n_coordinate_rows",
            "coord_storage_type",
            "source_text_filename",
            "notes",
        ],
    )


def write_previews(paths: dict[str, Path], n_rows: int = 25) -> None:
    preview_dir = paths["tables"] / "previews"
    preview_dir.mkdir(parents=True, exist_ok=True)
    for dataset in DATASETS:
        out_path = preview_dir / f"{dataset.dataset_id}.metadata_with_umap.first_{n_rows}_rows.tsv"
        with gzip.open(str(metadata_umap_path(paths, dataset)), "rt") as in_handle, out_path.open("w") as out_handle:
            for idx, line in enumerate(in_handle):
                if idx > n_rows:
                    break
                out_handle.write(line)


def write_completion(paths: dict[str, Path], datasets_by_id: dict[str, dict], metadata_rows: dict[str, int]) -> None:
    rows = [
        {"metric": "run_dir", "value": str(paths["run"])},
        {"metric": "base_url", "value": BASE_URL},
        {"metric": "n_datasets", "value": str(len(DATASETS))},
        {"metric": "download_scope", "value": "metadata/config/UMAP only; expression matrices intentionally not downloaded"},
    ]
    for dataset in DATASETS:
        rows.append({"metric": f"dataset_cells:{dataset.dataset_id}", "value": str(datasets_by_id[dataset.dataset_id].get("sampleCount", ""))})
        rows.append({"metric": f"metadata_with_umap_rows:{dataset.dataset_id}", "value": str(metadata_rows[dataset.dataset_id])})
    write_rows(paths["tables"] / "schmitz_2022_reference_metadata_umap_complete.tsv", rows, ["metric", "value"])


def write_readme(paths: dict[str, Path]) -> None:
    content = """# Schmitz 2022 Reference Metadata And UMAP

Source: UCSC Cell Browser collection `dev-inhibitory-neurons`, associated with
Schmitz et al. 2022, Nature, `10.1038/s41586-022-04510-w`.

This run stages Cell Browser configs, cell-level metadata, decoded UMAP
coordinates, and metadata-with-UMAP tables for the three available Schmitz
datasets:

- `schmitz_macaque_dev`
- `schmitz_mouse_dev`
- `schmitz_mouse_adult`

It intentionally does not download the expression matrices. The Cell Browser
matrix files are Float32 processed expression matrices and include negative
centered values, so raw-count or Seurat-reference use needs a separate decision.

Key tables:

- `tables/schmitz_2022_datasets.tsv`
- `tables/schmitz_2022_metadata_fields.tsv`
- `tables/schmitz_2022_metadata_value_counts.tsv`
- `tables/schmitz_2022_class_composition.tsv`
- `tables/schmitz_2022_reference_label_summary.tsv`
- `tables/schmitz_2022_umap_files.tsv`
- `metadata/*metadata_with_umap.tsv.gz`
"""
    (paths["run"] / "README_schmitz_2022_reference_metadata_umap.md").write_text(content)


def fetch_inventory(
    project_root: str | Path | None = None,
    run_label: str = RUN_LABEL_DEFAULT,
    skip_existing: bool = True,
) -> Path:
    paths = output_paths(project_root, run_label)
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)

    datasets_by_id: dict[str, dict] = {}
    coord_rows_by_key: dict[tuple[str, str], int] = {}
    metadata_rows: dict[str, int] = {}
    for dataset in DATASETS:
        log(f"Fetching config, metadata, and UMAP for {dataset.dataset_id}")
        download(dataset_url(dataset, "dataset.json"), dataset_json_path(paths, dataset), skip_existing=skip_existing)
        download(dataset_url(dataset, "desc.json"), desc_json_path(paths, dataset), skip_existing=skip_existing)
        stream_url_to_gzip(dataset_url(dataset, "meta.tsv"), meta_path(paths, dataset), skip_existing=skip_existing)
        conf = read_json(dataset_json_path(paths, dataset))
        datasets_by_id[dataset.dataset_id] = conf
        cell_ids = read_cell_ids(meta_path(paths, dataset))
        for coord in conf.get("coords", []):
            coord_name = coord["name"]
            download(
                dataset_url(dataset, f"coords/{coord_name}/coords.bin"),
                coord_bin_path(paths, dataset, coord_name),
                skip_existing=skip_existing,
            )
            n_rows = decode_coord_bin(coord_bin_path(paths, dataset, coord_name), cell_ids, coord, coord_tsv_path(paths, dataset, coord_name))
            coord_rows_by_key[(dataset.dataset_id, coord_name)] = n_rows
        metadata_rows[dataset.dataset_id] = write_metadata_with_umap(paths, dataset, conf)

    log("Writing Schmitz inventory tables")
    write_dataset_table(paths, datasets_by_id)
    write_metadata_fields(paths, datasets_by_id)
    write_metadata_value_counts(paths)
    write_class_composition(paths)
    write_reference_label_summary(paths)
    write_coord_summary(paths, datasets_by_id, coord_rows_by_key)
    write_previews(paths)
    write_completion(paths, datasets_by_id, metadata_rows)
    write_readme(paths)
    log(f"Complete: {paths['run']}")
    return paths["run"]


def main(argv: list[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", default=None)
    parser.add_argument("--run-label", default=RUN_LABEL_DEFAULT)
    parser.add_argument("--overwrite", action="store_true", help="Re-download existing files")
    args = parser.parse_args(argv)
    fetch_inventory(args.project_root, args.run_label, skip_existing=not args.overwrite)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
