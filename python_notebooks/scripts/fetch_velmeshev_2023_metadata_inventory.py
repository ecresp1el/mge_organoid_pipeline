#!/usr/bin/env python3
"""Fetch Velmeshev 2023 UCSC Cell Browser metadata/UMAP inventory.

This stages lightweight metadata needed to judge whether Velmeshev et al. 2023
is usable as a reference or comparison source. It deliberately does not
download expression matrices or the multi-GB H5AD.
"""

import csv
import gzip
import json
import os
import re
import shutil
import sys
import tempfile
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path


PROJECT_ROOT_DEFAULT = "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder"
RESULTS_DIRNAME = "velmeshev_2023_metadata_inventory"
RUN_LABEL_DEFAULT = "velmeshev_2023_metadata_inventory_v1"
BASE_URL = "https://pre-postnatal-cortex.cells.ucsc.edu/pre-postnatal-cortex"

VIEWS = [
    ("all_rna", "all/rna", "All cells snRNA-seq"),
    ("ex_neu_rna", "ex-neu/rna", "Excitatory neuron lineage snRNA-seq"),
    ("in_rna", "in/rna", "Interneuron lineage snRNA-seq"),
    ("glia_rna", "glia/rna", "Glial lineage snRNA-seq"),
    ("mg", "mg", "Microglial cells"),
    ("per", "per", "Pericytes"),
    ("end", "end", "Vascular cells"),
]

CLUSTER_LABEL_FIELDS = {
    "all_rna": ("Seurat_clusters", "Lineage"),
    "ex_neu_rna": ("Seurat_Clusters", "Cell_Type"),
    "in_rna": ("Seurat_Clusters", "Cell_Type"),
    "glia_rna": ("Seurat_Clusters", "Cell_Type"),
}

COMPOSITION_FIELDS = [
    "Lineage",
    "Cell_Type",
    "Age_Range",
    "Region_Broad",
    "Dataset",
    "Individual",
    "Sample",
    "Sex",
]


def log(message):
    print("[Velmeshev2023] {}".format(message), flush=True)


def project_root():
    return Path(os.environ.get("PROJECT_ROOT", PROJECT_ROOT_DEFAULT)).expanduser().resolve()


def run_dir():
    label = os.environ.get("VELMESHEV_2023_RUN_LABEL", RUN_LABEL_DEFAULT)
    return project_root() / "results" / RESULTS_DIRNAME / label


def url_for(view_path, filename):
    return "{}/{}/{}".format(BASE_URL, view_path, filename)


def download(url, dest, skip_existing=True):
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


def stream_url_to_gzip(url, dest, skip_existing=True):
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


def read_json(path):
    with path.open() as handle:
        return json.load(handle)


def open_meta(path):
    return gzip.open(str(path), "rt")


def natural_cluster_key(value):
    parts = []
    for part in re.split(r"([0-9]+)", str(value)):
        if part.isdigit():
            parts.append((0, int(part)))
        else:
            parts.append((1, part.lower()))
    return parts


def get_indices(header, names):
    return {name: header.index(name) for name in names if name in header}


def view_records(paths, datasets):
    rows = []
    for view_id, view_path, view_label in VIEWS:
        dataset = datasets[view_id]
        rows.append(
            {
                "view_id": view_id,
                "view_label": view_label,
                "ucsc_path": view_path,
                "sample_count": dataset.get("sampleCount", ""),
                "metadata_path": str(paths["metadata"] / "{}.meta.tsv.gz".format(view_id)),
                "umap_path": str(paths["umap"] / "{}.UMAP.coords.tsv.gz".format(view_id)),
                "dataset_json_path": str(paths["config"] / "{}.dataset.json".format(view_id)),
                "desc_json_path": str(paths["config"] / "{}.desc.json".format(view_id)),
            }
        )
    return rows


def write_rows(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_metadata_fields(paths, datasets):
    rows = []
    for view_id, _, _ in VIEWS:
        dataset = datasets[view_id]
        for field in dataset.get("metaFields", []):
            val_counts = field.get("valCounts") or []
            top_values = "; ".join("{}:{}".format(value, count) for value, count in val_counts[:15])
            rows.append(
                {
                    "view_id": view_id,
                    "sample_count": dataset.get("sampleCount", ""),
                    "field_name": field.get("name", ""),
                    "field_label": field.get("label", ""),
                    "field_type": field.get("type", ""),
                    "diff_val_count": field.get("diffValCount", ""),
                    "value_count_entries": len(val_counts),
                    "top_values": top_values,
                }
            )
    write_rows(
        paths["tables"] / "velmeshev_2023_metadata_fields.tsv",
        rows,
        [
            "view_id",
            "sample_count",
            "field_name",
            "field_label",
            "field_type",
            "diff_val_count",
            "value_count_entries",
            "top_values",
        ],
    )


def write_metadata_value_counts(paths, datasets):
    rows = []
    for view_id, _, _ in VIEWS:
        dataset = datasets[view_id]
        for field in dataset.get("metaFields", []):
            for value, count in field.get("valCounts") or []:
                rows.append(
                    {
                        "view_id": view_id,
                        "field_name": field.get("name", ""),
                        "value": value,
                        "count": count,
                    }
                )
    write_rows(
        paths["tables"] / "velmeshev_2023_metadata_value_counts.tsv",
        rows,
        ["view_id", "field_name", "value", "count"],
    )


def write_previews(paths, n_rows=25):
    preview_dir = paths["tables"] / "previews"
    preview_dir.mkdir(parents=True, exist_ok=True)
    for view_id, _, _ in VIEWS:
        meta_path = paths["metadata"] / "{}.meta.tsv.gz".format(view_id)
        out_path = preview_dir / "{}.first_{}_rows.tsv".format(view_id, n_rows)
        with open_meta(meta_path) as in_handle, out_path.open("w") as out_handle:
            for idx, line in enumerate(in_handle):
                if idx > n_rows:
                    break
                out_handle.write(line)


def write_cluster_label_crosswalk(paths):
    rows = []
    for view_id, fields in CLUSTER_LABEL_FIELDS.items():
        cluster_field, label_field = fields
        meta_path = paths["metadata"] / "{}.meta.tsv.gz".format(view_id)
        counts = defaultdict(Counter)
        with open_meta(meta_path) as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader)
            idx = get_indices(header, [cluster_field, label_field])
            for row in reader:
                counts[row[idx[cluster_field]]][row[idx[label_field]]] += 1
        for cluster_id in sorted(counts, key=natural_cluster_key):
            label_counts = counts[cluster_id]
            dominant_label, dominant_count = label_counts.most_common(1)[0]
            cluster_count = sum(label_counts.values())
            rows.append(
                {
                    "view_id": view_id,
                    "cluster_field": cluster_field,
                    "cluster_id": cluster_id,
                    "label_field": label_field,
                    "dominant_label": dominant_label,
                    "dominant_count": dominant_count,
                    "cluster_count": cluster_count,
                    "dominant_pct": "{:.2f}".format(100.0 * dominant_count / cluster_count),
                }
            )
    write_rows(
        paths["tables"] / "velmeshev_2023_cluster_label_crosswalk.tsv",
        rows,
        [
            "view_id",
            "cluster_field",
            "cluster_id",
            "label_field",
            "dominant_label",
            "dominant_count",
            "cluster_count",
            "dominant_pct",
        ],
    )


def write_cluster_composition(paths):
    rows = []
    for view_id, _, _ in VIEWS:
        meta_path = paths["metadata"] / "{}.meta.tsv.gz".format(view_id)
        cluster_field = "Seurat_clusters" if view_id == "all_rna" else "Seurat_Clusters"
        totals = Counter()
        counts = defaultdict(Counter)
        with open_meta(meta_path) as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader)
            fields = [cluster_field] + [field for field in COMPOSITION_FIELDS if field in header]
            idx = get_indices(header, fields)
            if cluster_field not in idx:
                continue
            for row in reader:
                cluster_id = row[idx[cluster_field]]
                totals[cluster_id] += 1
                for group_field in fields[1:]:
                    counts[(cluster_id, group_field)][row[idx[group_field]]] += 1
        for cluster_id in sorted(totals, key=natural_cluster_key):
            for group_field in [field for field in COMPOSITION_FIELDS if (cluster_id, field) in counts]:
                for group_value, count in counts[(cluster_id, group_field)].most_common():
                    rows.append(
                        {
                            "view_id": view_id,
                            "cluster_field": cluster_field,
                            "cluster_id": cluster_id,
                            "group_field": group_field,
                            "group_value": group_value,
                            "count": count,
                            "pct_of_cluster": "{:.2f}".format(100.0 * count / totals[cluster_id]),
                        }
                    )
    write_rows(
        paths["tables"] / "velmeshev_2023_cluster_composition.tsv",
        rows,
        ["view_id", "cluster_field", "cluster_id", "group_field", "group_value", "count", "pct_of_cluster"],
    )


def write_mge_summary(paths):
    rows = []
    for view_id, _, _ in VIEWS:
        meta_path = paths["metadata"] / "{}.meta.tsv.gz".format(view_id)
        total = 0
        counters = defaultdict(Counter)
        with open_meta(meta_path) as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader)
            wanted = [
                "Region_Broad",
                "Age_Range",
                "Lineage",
                "Cell_Type",
                "Seurat_clusters",
                "Seurat_Clusters",
                "Dataset",
                "Individual",
                "Sample",
                "Sex",
            ]
            idx = get_indices(header, wanted)
            if "Region_Broad" not in idx:
                continue
            for row in reader:
                if row[idx["Region_Broad"]] != "MGE":
                    continue
                total += 1
                for field_name, col_idx in idx.items():
                    if field_name == "Region_Broad":
                        continue
                    counters[field_name][row[col_idx]] += 1
        rows.append({"view_id": view_id, "group_field": "TOTAL", "group_value": "MGE", "count": total})
        for field_name in sorted(counters):
            for value, count in counters[field_name].most_common():
                rows.append({"view_id": view_id, "group_field": field_name, "group_value": value, "count": count})
    write_rows(
        paths["tables"] / "velmeshev_2023_mge_subset_summary.tsv",
        rows,
        ["view_id", "group_field", "group_value", "count"],
    )


def write_umap_summary(paths):
    rows = []
    for view_id, _, _ in VIEWS:
        umap_path = paths["umap"] / "{}.UMAP.coords.tsv.gz".format(view_id)
        n_rows = 0
        first_row = ""
        with gzip.open(str(umap_path), "rt") as handle:
            for line in handle:
                if n_rows == 0:
                    first_row = line.rstrip("\n")
                n_rows += 1
        rows.append(
            {
                "view_id": view_id,
                "umap_path": str(umap_path),
                "n_coordinate_rows": n_rows,
                "has_header": "no",
                "columns": "cell_id,umap_1,umap_2",
                "first_row": first_row,
            }
        )
    write_rows(
        paths["tables"] / "velmeshev_2023_umap_files.tsv",
        rows,
        ["view_id", "umap_path", "n_coordinate_rows", "has_header", "columns", "first_row"],
    )


def write_completion(paths, datasets):
    rows = [
        {"metric": "run_dir", "value": str(paths["run"])},
        {"metric": "base_url", "value": BASE_URL},
        {"metric": "n_views", "value": str(len(VIEWS))},
        {"metric": "download_scope", "value": "metadata/config/UMAP only; no expression matrix or H5AD"},
    ]
    for view_id, _, _ in VIEWS:
        rows.append({"metric": "view_cells:{}".format(view_id), "value": str(datasets[view_id].get("sampleCount", ""))})
    write_rows(paths["tables"] / "velmeshev_2023_metadata_inventory_complete.tsv", rows, ["metric", "value"])


def write_readme(paths):
    content = """# Velmeshev 2023 Metadata Inventory

Source: UCSC Cell Browser mirror of Velmeshev et al. 2023, Science, DOI `10.1126/science.adf0834`.

This run stages metadata, UMAP coordinates, and Cell Browser config files only.
It does not download the expression matrix or H5AD.

Key tables:

- `tables/velmeshev_2023_views.tsv`
- `tables/velmeshev_2023_metadata_fields.tsv`
- `tables/velmeshev_2023_metadata_value_counts.tsv`
- `tables/velmeshev_2023_cluster_label_crosswalk.tsv`
- `tables/velmeshev_2023_cluster_composition.tsv`
- `tables/velmeshev_2023_mge_subset_summary.tsv`
- `tables/velmeshev_2023_umap_files.tsv`
- `tables/previews/*.tsv`

UMAP files are headerless TSV.GZ files with columns:

```text
cell_id    umap_1    umap_2
```
"""
    (paths["run"] / "README_velmeshev_2023_metadata_inventory.md").write_text(content)


def main():
    paths = {
        "run": run_dir(),
        "config": run_dir() / "config",
        "metadata": run_dir() / "metadata",
        "umap": run_dir() / "umap",
        "tables": run_dir() / "tables",
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)

    datasets = {}
    for view_id, view_path, _ in VIEWS:
        log("Fetching config and metadata for {} ({})".format(view_id, view_path))
        dataset_json = paths["config"] / "{}.dataset.json".format(view_id)
        desc_json = paths["config"] / "{}.desc.json".format(view_id)
        meta_gz = paths["metadata"] / "{}.meta.tsv.gz".format(view_id)
        umap_gz = paths["umap"] / "{}.UMAP.coords.tsv.gz".format(view_id)
        download(url_for(view_path, "dataset.json"), dataset_json)
        download(url_for(view_path, "desc.json"), desc_json)
        stream_url_to_gzip(url_for(view_path, "meta.tsv"), meta_gz)
        download(url_for(view_path, "UMAP.coords.tsv.gz"), umap_gz)
        datasets[view_id] = read_json(dataset_json)

    write_rows(
        paths["tables"] / "velmeshev_2023_views.tsv",
        view_records(paths, datasets),
        [
            "view_id",
            "view_label",
            "ucsc_path",
            "sample_count",
            "metadata_path",
            "umap_path",
            "dataset_json_path",
            "desc_json_path",
        ],
    )
    log("Writing inventory tables")
    write_metadata_fields(paths, datasets)
    write_metadata_value_counts(paths, datasets)
    write_previews(paths)
    write_cluster_label_crosswalk(paths)
    write_cluster_composition(paths)
    write_mge_summary(paths)
    write_umap_summary(paths)
    write_completion(paths, datasets)
    write_readme(paths)
    log("Complete: {}".format(paths["run"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
