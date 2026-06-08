"""Cross-study marker-expression preparation and plotting helpers.

This module is the Python-side home for the first major cross-study marker
expression figure.  It deliberately works from a compact per-cell table
contract so notebooks can iterate on plot design without reloading large
Seurat objects.
"""

from __future__ import annotations

import argparse
import gzip
from dataclasses import dataclass
from pathlib import Path
import re
from typing import Iterable, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
import numpy as np
import pandas as pd
from scipy import sparse

from .paths import resolve_project_root


MARKER_EXPRESSION_SCHEMA_VERSION = "cross_study_marker_expression_v8"

GREY_TO_RED_CMAP = LinearSegmentedColormap.from_list("grey_to_red_expression", ["#d0d0d0", "#d7301f"])

ON_TARGET_GENES = [
    "DCX",
    "GAD2",
    "DLX5",
    "LHX6",
    "MAF",
    "SST",
    "PVALB",
    "ERBB4",
    "MEF2C",
    "MAFB",
    "LHX8",
    "NKX2-1",
]

OFF_TARGET_GENES = ["SP8", "PAX6", "NEUROD2", "ISL1", "ACHE", "NKX6-2", "MKI67"]
GENE_ORDER = ON_TARGET_GENES + OFF_TARGET_GENES

GENE_ALIASES = {
    "NKX2-1": ["NKX2-1", "NKX2.1"],
    "NKX6-2": ["NKX6-2", "NKX6.2"],
}

BASE_COLUMNS = [
    "cell_id",
    "study_id",
    "study_label",
    "sample",
    "cluster",
    "umap_1",
    "umap_2",
]


@dataclass(frozen=True)
class CrossStudyMarkerSpec:
    """Description of one study for cross-study marker-expression plotting."""

    study_id: str
    study_label: str
    seurat_path: str
    h5ad_path: str = ""
    reduction: str = "umap"
    assay: str = "RNA"
    expression_layer: str = "data"
    expression_transform: str = "none"
    sample_col: str = ""
    cluster_col: str = "seurat_clusters"
    feature_map_path: str = ""
    include_in_first_plot: bool = True
    note: str = ""

    def resolved_seurat_path(self, project_root: str | Path | None = None) -> Path:
        return resolve_project_path(self.seurat_path, project_root)

    def resolved_h5ad_path(self, project_root: str | Path | None = None) -> Path | None:
        if not self.h5ad_path:
            return None
        return resolve_project_path(self.h5ad_path, project_root)

    def resolved_feature_map_path(self, project_root: str | Path | None = None) -> Path | None:
        if not self.feature_map_path:
            return None
        return resolve_project_path(self.feature_map_path, project_root)


def resolve_project_path(path: str | Path, project_root: str | Path | None = None) -> Path:
    """Resolve an absolute path or a PROJECT_ROOT-relative path."""
    raw = Path(path).expanduser()
    if raw.is_absolute():
        return raw
    return resolve_project_root(project_root) / raw


def safe_token(value: object) -> str:
    """Return a stable filename/column token."""
    token = re.sub(r"[^A-Za-z0-9]+", "_", str(value).strip())
    token = re.sub(r"_+", "_", token).strip("_")
    return token or "value"


def natural_sort_key(value: object) -> list[object]:
    """Sort strings with embedded numbers in human order."""
    parts = re.split(r"(\d+)", str(value))
    return [int(part) if part.isdigit() else part for part in parts]


def default_cross_study_marker_specs(
    project_root: str | Path | None = None,
    *,
    include_xiang: bool = False,
) -> list[CrossStudyMarkerSpec]:
    """Return the canonical study set for the first Python marker plot."""
    root = resolve_project_root(project_root)
    return [
        CrossStudyMarkerSpec(
            study_id="varela_div30",
            study_label="This Study, DIV 30",
            seurat_path="results/varela_this_paper/varela_this_paper_seurat.rds",
            h5ad_path="results/python_anndata/varela_div30.h5ad",
            sample_col="orig.ident",
            cluster_col="seurat_clusters",
        ),
        CrossStudyMarkerSpec(
            study_id="varela_div90",
            study_label="This Study, DIV 90",
            seurat_path="/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/umap_props_output/clustered_day90_with_cluster_names_2.rds",
            h5ad_path="results/python_anndata/varela_div90.h5ad",
            sample_col="orig.ident",
            cluster_col="seurat_clusters",
        ),
        CrossStudyMarkerSpec(
            study_id="siebert_2026",
            study_label="Siebert et al. 2026",
            seurat_path="results/siebert_2026/siebert_2026_seurat.rds",
            expression_layer="counts",
            expression_transform="log1p_cp10k",
            sample_col="orig.ident",
            cluster_col="seurat_clusters",
            note="RNA counts are transformed to log1p(CP10K) during export.",
        ),
        CrossStudyMarkerSpec(
            study_id="walsh",
            study_label="Walsh et al. 2025",
            seurat_path="results/walsh_day75/walsh_day75_final_annotated.rds",
            reduction="umap_sel",
            sample_col="sample_id",
            cluster_col="seurat_clusters",
        ),
        CrossStudyMarkerSpec(
            study_id="bershteyn_2025",
            study_label="Bershteyn et al. 2025",
            seurat_path="results/bershteyn_2025/bershteyn_2025_seurat.rds",
            sample_col="sample",
            cluster_col="seurat_clusters",
        ),
        CrossStudyMarkerSpec(
            study_id="bershteyn_2023",
            study_label="Bershteyn et al. 2023",
            seurat_path="results/bershteyn_2023/bershteyn_2023_seurat.rds",
            sample_col="orig.ident",
            cluster_col="seurat_clusters",
        ),
        CrossStudyMarkerSpec(
            study_id="samarasinghe_2021",
            study_label="Samarasinghe et al. 2021",
            seurat_path="results/samarasinghe_2021_zenodo_processed_object/samarasinghe_2021_zenodo_seurat.rds",
            sample_col="orig.ident",
            cluster_col="seurat_clusters",
            note="Use the official Zenodo processed object.",
        ),
        CrossStudyMarkerSpec(
            study_id="xiang_2018",
            study_label="Xiang",
            seurat_path="results/xiang_2018/xiang_2018_seurat.rds",
            sample_col="orig.ident",
            cluster_col="seurat_clusters",
            feature_map_path="data/raw/xiang_2018_geo_files/suppl/GSE98201_genes.tsv.gz",
            include_in_first_plot=include_xiang,
            note="Excluded by default until sample biology/metadata are validated.",
        ),
    ]


def included_specs(specs: Sequence[CrossStudyMarkerSpec]) -> list[CrossStudyMarkerSpec]:
    """Return studies included in the first major plot."""
    return [spec for spec in specs if spec.include_in_first_plot]


def exclude_specs(
    specs: Sequence[CrossStudyMarkerSpec],
    exclude_study_ids: Sequence[str] | None = None,
) -> list[CrossStudyMarkerSpec]:
    """Return specs after dropping selected studies from plotting."""
    excluded = {str(study_id).strip() for study_id in (exclude_study_ids or []) if str(study_id).strip()}
    if not excluded:
        return list(specs)
    return [spec for spec in specs if spec.study_id not in excluded]


def run_dir(project_root: str | Path | None = None, run_label: str = "cross_study_marker_expression_v8") -> Path:
    """Return the run directory for this workflow."""
    return resolve_project_root(project_root) / "results" / "cross_study_marker_expression" / run_label


def table_dir(project_root: str | Path | None = None, run_label: str = "cross_study_marker_expression_v8") -> Path:
    """Return the table directory for this workflow."""
    return run_dir(project_root, run_label) / "tables"


def plot_dir(project_root: str | Path | None = None, run_label: str = "cross_study_marker_expression_v8") -> Path:
    """Return the plot directory for this workflow."""
    return run_dir(project_root, run_label) / "plots"


def per_study_table_path(
    study_id: str,
    project_root: str | Path | None = None,
    run_label: str = "cross_study_marker_expression_v8",
) -> Path:
    """Return the standardized per-study marker table path."""
    return table_dir(project_root, run_label) / "per_study" / f"{study_id}_marker_expression.tsv.gz"


def ensure_output_dirs(project_root: str | Path | None = None, run_label: str = "cross_study_marker_expression_v8") -> dict[str, Path]:
    """Create and return output directories."""
    paths = {
        "run_dir": run_dir(project_root, run_label),
        "table_dir": table_dir(project_root, run_label),
        "per_study_dir": table_dir(project_root, run_label) / "per_study",
        "plot_dir": plot_dir(project_root, run_label),
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    return paths


def gene_group(gene: str) -> str:
    """Return ON-target, OFF-target, or ungrouped for a marker."""
    if gene in ON_TARGET_GENES:
        return "ON-target"
    if gene in OFF_TARGET_GENES:
        return "OFF-target"
    return "ungrouped"


def marker_gene_table(genes: Sequence[str] = GENE_ORDER) -> pd.DataFrame:
    """Return the canonical marker gene table."""
    return pd.DataFrame(
        {
            "gene": list(genes),
            "gene_order": np.arange(1, len(genes) + 1),
            "gene_group": [gene_group(gene) for gene in genes],
            "schema_version": MARKER_EXPRESSION_SCHEMA_VERSION,
        }
    )


def study_table(specs: Sequence[CrossStudyMarkerSpec], project_root: str | Path | None = None) -> pd.DataFrame:
    """Return the canonical study table for notebook display/audit."""
    rows = []
    for idx, spec in enumerate(specs, start=1):
        h5ad_path = spec.resolved_h5ad_path(project_root)
        feature_map_path = spec.resolved_feature_map_path(project_root)
        rows.append(
            {
                "study_id": spec.study_id,
                "study_label": spec.study_label,
                "study_order": idx,
                "include_in_first_plot": spec.include_in_first_plot,
                "seurat_path": str(spec.resolved_seurat_path(project_root)),
                "h5ad_path": "" if h5ad_path is None else str(h5ad_path),
                "reduction": spec.reduction,
                "assay": spec.assay,
                "expression_layer": spec.expression_layer,
                "expression_transform": spec.expression_transform,
                "sample_col": spec.sample_col,
                "cluster_col": spec.cluster_col,
                "feature_map_path": "" if feature_map_path is None else str(feature_map_path),
                "note": spec.note,
            }
        )
    return pd.DataFrame(rows)


def _case_insensitive_lookup(values: Iterable[str]) -> tuple[dict[str, str], set[str]]:
    lookup: dict[str, str] = {}
    ambiguous: set[str] = set()
    for value in values:
        text = str(value)
        key = text.upper()
        if key in lookup and lookup[key] != text:
            ambiguous.add(key)
        else:
            lookup[key] = text
    return lookup, ambiguous


def match_genes_to_var_names(genes: Sequence[str], var_names: Sequence[str]) -> tuple[dict[str, str], pd.DataFrame]:
    """Match requested genes against AnnData var names with alias fallback."""
    var_values = [str(value) for value in var_names]
    var_set = set(var_values)
    upper_lookup, ambiguous = _case_insensitive_lookup(var_values)
    rows = []
    matches = {}
    for gene in genes:
        candidates = GENE_ALIASES.get(gene, [gene])
        matched = ""
        match_type = "missing"
        for candidate in candidates:
            key = candidate.upper()
            if candidate in var_set:
                matched = candidate
                match_type = "exact" if candidate == gene else "alias_exact"
                break
            if key in ambiguous:
                match_type = "ambiguous_case_insensitive"
                continue
            if key in upper_lookup:
                matched = upper_lookup[key]
                match_type = "case_insensitive" if candidate == gene else "alias_case_insensitive"
                break
        if matched:
            matches[gene] = matched
        rows.append(
            {
                "gene": gene,
                "matched_feature": matched,
                "matched": bool(matched),
                "match_type": match_type,
                "aliases_considered": ",".join(candidates),
            }
        )
    return matches, pd.DataFrame(rows)


def _matrix_to_dense_array(matrix) -> np.ndarray:
    if sparse.issparse(matrix):
        return np.asarray(matrix.toarray())
    return np.asarray(matrix)


def _h5_attr_text(value: object) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def _decode_h5_values(values: np.ndarray) -> np.ndarray:
    if values.dtype.kind in {"S", "O", "U"}:
        decoded = []
        for value in values:
            if isinstance(value, bytes):
                decoded.append(value.decode("utf-8"))
            else:
                decoded.append(str(value))
        return np.asarray(decoded, dtype=object)
    return values


def _read_h5ad_string_or_array(node) -> np.ndarray:
    """Read a simple H5AD dataset as strings/numbers."""
    try:
        values = node.asstr()[:]
    except (AttributeError, TypeError):
        values = node[:]
    return _decode_h5_values(np.asarray(values))


def _read_h5ad_dataframe_value(frame_group, column: str) -> np.ndarray | None:
    """Read one AnnData dataframe column from an H5AD group."""
    import h5py

    if column not in frame_group:
        return None
    node = frame_group[column]
    if isinstance(node, h5py.Dataset):
        return _read_h5ad_string_or_array(node)
    encoding = _h5_attr_text(node.attrs.get("encoding-type", ""))
    if encoding == "categorical" or {"categories", "codes"}.issubset(node.keys()):
        categories = _read_h5ad_string_or_array(node["categories"])
        codes = np.asarray(node["codes"][:], dtype=np.int64)
        values = np.full(codes.shape[0], "", dtype=object)
        valid = codes >= 0
        values[valid] = categories[codes[valid]]
        return values
    return None


def _read_h5ad_dataframe_index(frame_group) -> np.ndarray:
    index_name = _h5_attr_text(frame_group.attrs.get("_index", ""))
    if index_name and index_name in frame_group:
        values = _read_h5ad_dataframe_value(frame_group, index_name)
        if values is not None:
            return values
    for key in frame_group.keys():
        values = _read_h5ad_dataframe_value(frame_group, key)
        if values is not None:
            return values
    raise ValueError("Could not read H5AD dataframe index.")


def _read_h5ad_obsm(handle, obsm_keys: Sequence[str]) -> tuple[np.ndarray, str]:
    for key in obsm_keys:
        path = f"obsm/{key}"
        if path in handle:
            values = np.asarray(handle[path][:], dtype=float)
            if values.ndim == 2 and values.shape[1] >= 2:
                return values[:, :2], key
    raise ValueError(f"H5AD lacks one of these UMAP keys: {', '.join(obsm_keys)}")


def _read_h5ad_matrix_columns(
    h5ad_path: str | Path,
    matrix_key: str,
    target_indices: Sequence[int],
    shape: tuple[int, int],
    chunk_rows: int = 2048,
) -> np.ndarray:
    """Read selected matrix columns from an H5AD without materializing all genes."""
    import h5py

    target_indices = list(map(int, target_indices))
    n_obs, n_vars = shape
    result = np.zeros((n_obs, len(target_indices)), dtype=np.float32)
    with h5py.File(h5ad_path, "r") as handle:
        if matrix_key not in handle:
            raise ValueError(f"H5AD lacks matrix key: {matrix_key}")
        node = handle[matrix_key]
        if isinstance(node, h5py.Group):
            encoding = _h5_attr_text(node.attrs.get("encoding-type", ""))
            if encoding in {"csr_matrix", "csc_matrix"} or {"data", "indices", "indptr"}.issubset(node.keys()):
                data_ds = node["data"]
                indices_ds = node["indices"]
                indptr = np.asarray(node["indptr"][:], dtype=np.int64)
                if encoding == "csc_matrix":
                    for out_idx, col_idx in enumerate(target_indices):
                        start = int(indptr[col_idx])
                        stop = int(indptr[col_idx + 1])
                        rows = np.asarray(indices_ds[start:stop], dtype=np.int64)
                        values = np.asarray(data_ds[start:stop], dtype=np.float32)
                        result[rows, out_idx] = values
                    return result

                # AnnData normally writes X as CSR: rows are cells, columns are genes.
                for row_start in range(0, n_obs, chunk_rows):
                    row_stop = min(row_start + chunk_rows, n_obs)
                    data_start = int(indptr[row_start])
                    data_stop = int(indptr[row_stop])
                    chunk_indptr = indptr[row_start : row_stop + 1] - data_start
                    chunk_indices = np.asarray(indices_ds[data_start:data_stop], dtype=np.int64)
                    chunk_data = np.asarray(data_ds[data_start:data_stop], dtype=np.float32)
                    chunk = sparse.csr_matrix(
                        (chunk_data, chunk_indices, chunk_indptr),
                        shape=(row_stop - row_start, n_vars),
                    )
                    result[row_start:row_stop, :] = chunk[:, target_indices].toarray()
                return result
            raise ValueError(f"Unsupported H5AD matrix group encoding for {matrix_key}: {encoding}")

        # Dense matrix path. h5py requires fancy-indexed columns to be sorted.
        sorted_pairs = sorted(enumerate(target_indices), key=lambda pair: pair[1])
        sorted_indices = [pair[1] for pair in sorted_pairs]
        dense = np.asarray(node[:, sorted_indices], dtype=np.float32)
        for sorted_pos, (out_pos, _) in enumerate(sorted_pairs):
            result[:, out_pos] = dense[:, sorted_pos]
        return result


def _extract_obs_column(obs: pd.DataFrame, preferred: str, fallbacks: Sequence[str]) -> pd.Series:
    candidates = [preferred] if preferred else []
    candidates.extend([value for value in fallbacks if value and value not in candidates])
    for col in candidates:
        if col in obs.columns:
            return obs[col].astype(str)
    return pd.Series([""] * obs.shape[0], index=obs.index, dtype="object")


def extract_marker_expression_from_h5ad(
    spec: CrossStudyMarkerSpec,
    output_path: str | Path,
    project_root: str | Path | None = None,
    genes: Sequence[str] = GENE_ORDER,
    layer: str | None = None,
    obsm_keys: Sequence[str] = ("X_umap_seurat", "X_umap"),
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Write one standardized marker-expression table from an AnnData cache.

    Returns `(table, gene_match_table)`.
    """
    import h5py

    h5ad_path = spec.resolved_h5ad_path(project_root)
    if h5ad_path is None:
        raise ValueError(f"{spec.study_id} has no h5ad_path configured.")
    if not h5ad_path.exists():
        raise FileNotFoundError(f"Missing H5AD for {spec.study_id}: {h5ad_path}")

    with h5py.File(h5ad_path, "r") as handle:
        if "var" not in handle or "obs" not in handle:
            raise ValueError(f"{spec.study_id} H5AD lacks obs/var groups.")
        var_names = _read_h5ad_dataframe_index(handle["var"]).astype(str)
        cell_ids = _read_h5ad_dataframe_index(handle["obs"]).astype(str)
        sample_values = _read_h5ad_dataframe_value(handle["obs"], spec.sample_col) if spec.sample_col else None
        if sample_values is None:
            for fallback in ["sample_id", "sample", "orig.ident"]:
                sample_values = _read_h5ad_dataframe_value(handle["obs"], fallback)
                if sample_values is not None:
                    break
        cluster_values = _read_h5ad_dataframe_value(handle["obs"], spec.cluster_col) if spec.cluster_col else None
        if cluster_values is None:
            cluster_values = _read_h5ad_dataframe_value(handle["obs"], "seurat_clusters")
        umap, _ = _read_h5ad_obsm(handle, obsm_keys)
        x_node = handle["X"]
        shape_attr = x_node.attrs.get("shape", None) if hasattr(x_node, "attrs") else None
        if shape_attr is None:
            shape = tuple(x_node.shape)
        else:
            shape = tuple(map(int, shape_attr))

    matches, match_table = match_genes_to_var_names(genes, var_names)
    missing = [gene for gene in genes if gene not in matches]
    if missing:
        raise ValueError(f"{spec.study_id} is missing marker genes in H5AD: {', '.join(missing)}")

    n_obs, n_vars = int(shape[0]), int(shape[1])
    if cell_ids.shape[0] != n_obs:
        raise ValueError(f"{spec.study_id} cell_id length mismatch: {cell_ids.shape[0]} != {n_obs}")
    if umap.shape[0] != n_obs:
        raise ValueError(f"{spec.study_id} UMAP row mismatch: {umap.shape[0]} != {n_obs}")
    if sample_values is None:
        sample_values = np.full(n_obs, "", dtype=object)
    if cluster_values is None:
        cluster_values = np.full(n_obs, "", dtype=object)

    table = pd.DataFrame(
        {
            "cell_id": cell_ids,
            "study_id": spec.study_id,
            "study_label": spec.study_label,
            "sample": sample_values.astype(str),
            "cluster": cluster_values.astype(str),
            "umap_1": umap[:, 0],
            "umap_2": umap[:, 1],
        }
    )
    var_positions = {str(name): idx for idx, name in enumerate(var_names)}
    target_indices = [var_positions[matches[gene]] for gene in genes]
    matrix_key = "X" if not layer else f"layers/{layer}"
    expression = _read_h5ad_matrix_columns(
        h5ad_path=h5ad_path,
        matrix_key=matrix_key,
        target_indices=target_indices,
        shape=(n_obs, n_vars),
    )
    for idx, gene in enumerate(genes):
        table[gene] = expression[:, idx]

    out = Path(output_path).expanduser()
    out.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(out, sep="\t", index=False)
    return table, match_table


def count_tsv_rows(path: str | Path) -> int:
    """Count data rows in a TSV/TSV.GZ file without loading it into pandas."""
    p = Path(path)
    opener = gzip.open if p.suffix == ".gz" else open
    with opener(p, "rt") as handle:
        return max(sum(1 for _ in handle) - 1, 0)


def validate_marker_expression_table(
    path: str | Path,
    genes: Sequence[str] = GENE_ORDER,
) -> dict[str, object]:
    """Return validation status for a standardized marker-expression table."""
    p = Path(path).expanduser()
    status: dict[str, object] = {
        "marker_table_path": str(p),
        "marker_table_exists": p.exists(),
        "python_ready": False,
        "n_cells_marker_table": np.nan,
        "missing_required_columns": "",
        "missing_marker_columns": "",
    }
    if not p.exists():
        return status
    header = pd.read_csv(p, sep="\t", nrows=0).columns.tolist()
    missing_required = [col for col in BASE_COLUMNS if col not in header]
    missing_markers = [gene for gene in genes if gene not in header]
    status["missing_required_columns"] = ",".join(missing_required)
    status["missing_marker_columns"] = ",".join(missing_markers)
    status["python_ready"] = not missing_required and not missing_markers
    status["n_cells_marker_table"] = count_tsv_rows(p)
    return status


def readiness_table(
    specs: Sequence[CrossStudyMarkerSpec] | None = None,
    project_root: str | Path | None = None,
    run_label: str = "cross_study_marker_expression_v8",
    genes: Sequence[str] = GENE_ORDER,
) -> pd.DataFrame:
    """Return readiness for source paths, H5AD caches, and Python marker tables."""
    specs = list(specs or default_cross_study_marker_specs(project_root))
    rows = []
    for spec in specs:
        h5ad_path = spec.resolved_h5ad_path(project_root)
        marker_path = per_study_table_path(spec.study_id, project_root, run_label)
        validation = validate_marker_expression_table(marker_path, genes)
        needs = []
        seurat_path = spec.resolved_seurat_path(project_root)
        if not seurat_path.exists():
            needs.append("missing Seurat source")
        if not validation["python_ready"]:
            if h5ad_path is not None and h5ad_path.exists():
                needs.append("extract marker table from H5AD")
            else:
                needs.append("export marker table from Seurat")
        rows.append(
            {
                "study_id": spec.study_id,
                "study_label": spec.study_label,
                "include_in_first_plot": spec.include_in_first_plot,
                "seurat_path": str(seurat_path),
                "seurat_exists": seurat_path.exists(),
                "h5ad_path": "" if h5ad_path is None else str(h5ad_path),
                "h5ad_exists": False if h5ad_path is None else h5ad_path.exists(),
                "sample_col": spec.sample_col,
                "cluster_col": spec.cluster_col,
                "reduction": spec.reduction,
                "assay": spec.assay,
                "expression_layer": spec.expression_layer,
                "expression_transform": spec.expression_transform,
                "note": spec.note,
                **validation,
                "next_action": "; ".join(needs),
            }
        )
    return pd.DataFrame(rows)


def write_setup_tables(
    project_root: str | Path | None = None,
    run_label: str = "cross_study_marker_expression_v8",
    include_xiang: bool = False,
) -> dict[str, Path]:
    """Write canonical setup/readiness tables for the notebook."""
    paths = ensure_output_dirs(project_root, run_label)
    specs = default_cross_study_marker_specs(project_root, include_xiang=include_xiang)
    study_path = paths["table_dir"] / "cross_study_marker_expression_studies.tsv"
    gene_path = paths["table_dir"] / "cross_study_marker_expression_genes.tsv"
    readiness_path = paths["table_dir"] / "cross_study_marker_expression_readiness.tsv"
    study_table(specs, project_root).to_csv(study_path, sep="\t", index=False)
    marker_gene_table().to_csv(gene_path, sep="\t", index=False)
    readiness_table(specs, project_root, run_label).to_csv(readiness_path, sep="\t", index=False)
    return {
        "study_table": study_path,
        "gene_table": gene_path,
        "readiness_table": readiness_path,
    }


def extract_available_h5ad_marker_tables(
    project_root: str | Path | None = None,
    run_label: str = "cross_study_marker_expression_v8",
    study_ids: Sequence[str] | None = None,
    include_xiang: bool = False,
) -> pd.DataFrame:
    """Extract marker tables for studies that already have configured H5AD files."""
    ensure_output_dirs(project_root, run_label)
    specs = default_cross_study_marker_specs(project_root, include_xiang=include_xiang)
    if study_ids:
        wanted = set(study_ids)
        specs = [spec for spec in specs if spec.study_id in wanted]
    rows = []
    for spec in specs:
        h5ad_path = spec.resolved_h5ad_path(project_root)
        if h5ad_path is None or not h5ad_path.exists():
            rows.append(
                {
                    "study_id": spec.study_id,
                    "status": "skipped_no_h5ad",
                    "output_path": str(per_study_table_path(spec.study_id, project_root, run_label)),
                    "n_cells": np.nan,
                    "message": "No configured existing H5AD cache.",
                }
            )
            continue
        output_path = per_study_table_path(spec.study_id, project_root, run_label)
        table, match_table = extract_marker_expression_from_h5ad(spec, output_path, project_root)
        match_path = output_path.with_name(output_path.name.replace("_marker_expression.tsv.gz", "_gene_matches.tsv"))
        match_table.to_csv(match_path, sep="\t", index=False)
        rows.append(
            {
                "study_id": spec.study_id,
                "status": "ok",
                "output_path": str(output_path),
                "n_cells": table.shape[0],
                "message": "",
            }
        )
    summary = pd.DataFrame(rows)
    summary.to_csv(table_dir(project_root, run_label) / "cross_study_marker_expression_h5ad_extract_summary.tsv", sep="\t", index=False)
    readiness_table(default_cross_study_marker_specs(project_root, include_xiang=include_xiang), project_root, run_label).to_csv(
        table_dir(project_root, run_label) / "cross_study_marker_expression_readiness.tsv",
        sep="\t",
        index=False,
    )
    return summary


def load_marker_expression_tables(
    specs: Sequence[CrossStudyMarkerSpec] | None = None,
    project_root: str | Path | None = None,
    run_label: str = "cross_study_marker_expression_v8",
    genes: Sequence[str] = GENE_ORDER,
    require_all: bool = True,
) -> pd.DataFrame:
    """Load standardized per-study marker-expression tables."""
    specs = included_specs(list(specs or default_cross_study_marker_specs(project_root)))
    tables = []
    missing = []
    usecols = BASE_COLUMNS + list(genes)
    for spec in specs:
        path = per_study_table_path(spec.study_id, project_root, run_label)
        validation = validate_marker_expression_table(path, genes)
        if not validation["python_ready"]:
            missing.append(f"{spec.study_id}: {validation['marker_table_path']}")
            continue
        tables.append(pd.read_csv(path, sep="\t", usecols=usecols))
    if missing and require_all:
        raise FileNotFoundError("Missing or incomplete marker-expression tables:\n" + "\n".join(missing))
    if not tables:
        raise FileNotFoundError("No marker-expression tables were loaded.")
    data = pd.concat(tables, ignore_index=True)
    data["study_id"] = data["study_id"].astype(str)
    data["study_label"] = data["study_label"].astype(str)
    return data


def downsample_by_study(data: pd.DataFrame, max_cells_per_study: int | None, random_state: int = 0) -> pd.DataFrame:
    """Downsample each study for faster plotting while keeping deterministic output."""
    if max_cells_per_study is None or max_cells_per_study <= 0:
        return data
    parts = []
    rng = np.random.default_rng(random_state)
    for _, group in data.groupby("study_id", sort=False):
        if group.shape[0] <= max_cells_per_study:
            parts.append(group)
            continue
        idx = rng.choice(group.index.to_numpy(), size=max_cells_per_study, replace=False)
        parts.append(group.loc[np.sort(idx)])
    return pd.concat(parts, ignore_index=True)


def expression_limits(
    data: pd.DataFrame,
    genes: Sequence[str],
    quantile: float = 0.99,
    minimum: float = 1e-6,
) -> dict[str, float]:
    """Return per-gene upper limits for marker-expression color scales."""
    limits = {}
    for gene in genes:
        values = pd.to_numeric(data[gene], errors="coerce").to_numpy(dtype=float)
        values = values[np.isfinite(values)]
        if values.size == 0:
            limits[gene] = minimum
            continue
        positive = values[values > 0]
        if positive.size == 0:
            limits[gene] = minimum
            continue
        limits[gene] = max(float(np.quantile(positive, quantile)), minimum)
    return limits


def expression_scale_groups(genes: Sequence[str]) -> list[tuple[str, list[str]]]:
    """Return marker groups that should share one expression color scale."""
    requested = list(genes)
    groups: list[tuple[str, list[str]]] = []
    on_genes = [gene for gene in ON_TARGET_GENES if gene in requested]
    off_genes = [gene for gene in OFF_TARGET_GENES if gene in requested]
    if on_genes:
        groups.append(("ON-target", on_genes))
    if off_genes:
        groups.append(("OFF-target", off_genes))
    grouped = set(on_genes + off_genes)
    other_genes = [gene for gene in requested if gene not in grouped]
    if other_genes:
        groups.append(("Other", other_genes))
    return groups


def expression_group_limits(
    data: pd.DataFrame,
    scale_groups: Sequence[tuple[str, Sequence[str]]],
    quantile: float = 0.99,
    minimum: float = 1e-6,
) -> dict[str, float]:
    """Return one upper expression limit per marker group."""
    limits = {}
    for group_label, group_genes in scale_groups:
        columns = [gene for gene in group_genes if gene in data.columns]
        if not columns:
            limits[group_label] = minimum
            continue
        values = data[columns].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float).ravel()
        values = values[np.isfinite(values)]
        positive = values[values > 0]
        if positive.size == 0:
            limits[group_label] = minimum
            continue
        limits[group_label] = max(float(np.quantile(positive, quantile)), minimum)
    return limits


def expression_summary_table(
    data: pd.DataFrame,
    genes: Sequence[str] = GENE_ORDER,
) -> pd.DataFrame:
    """Summarize marker-expression values by study for validation."""
    rows = []
    for (study_id, study_label), group in data.groupby(["study_id", "study_label"], sort=False):
        values = group[list(genes)].stack()
        values = pd.to_numeric(values, errors="coerce")
        values = values[np.isfinite(values)]
        rows.append(
            {
                "study_id": study_id,
                "study_label": study_label,
                "n_cells": group.shape[0],
                "n_gene_values": values.shape[0],
                "min": float(values.min()) if values.shape[0] else np.nan,
                "q50": float(values.quantile(0.50)) if values.shape[0] else np.nan,
                "q95": float(values.quantile(0.95)) if values.shape[0] else np.nan,
                "q99": float(values.quantile(0.99)) if values.shape[0] else np.nan,
                "max": float(values.max()) if values.shape[0] else np.nan,
                "pct_nonzero": float((values > 0).mean() * 100.0) if values.shape[0] else np.nan,
            }
        )
    return pd.DataFrame(rows)


def _ordered_study_labels(data: pd.DataFrame, specs: Sequence[CrossStudyMarkerSpec] | None = None) -> list[str]:
    if specs:
        ordered = [spec.study_label for spec in specs if spec.include_in_first_plot]
        present = set(data["study_label"].astype(str))
        return [label for label in ordered if label in present]
    return list(pd.unique(data["study_label"].astype(str)))


def _apply_internal_umap_plot_filters(data: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Apply hard-coded biological filters that define what each UMAP row shows."""
    summaries = []
    filtered_parts = []
    for study_id, group in data.groupby("study_id", sort=False):
        rule = "none"
        kept = group
        if study_id == "samarasinghe_2021":
            rule = "samarasinghe_2021_controls_only"
            control = group["sample"].astype(str).str.contains("Ctrl", case=False, na=False)
            kept = group.loc[control].copy()
            if kept.empty:
                raise ValueError("Samarasinghe control-only UMAP filter removed all cells.")
            non_control_samples = sorted(kept.loc[~kept["sample"].astype(str).str.contains("Ctrl", case=False, na=False), "sample"].unique())
            if non_control_samples:
                raise ValueError("Samarasinghe control-only UMAP filter retained non-control samples: " + ", ".join(non_control_samples))
        summaries.append(
            {
                "study_id": study_id,
                "study_label": str(group["study_label"].iloc[0]),
                "plot_filter": rule,
                "plot_filter_label": "Controls only" if rule == "samarasinghe_2021_controls_only" else "",
                "n_cells_before_plot_filter": int(group.shape[0]),
                "n_cells_after_plot_filter": int(kept.shape[0]),
                "n_cells_removed_by_plot_filter": int(group.shape[0] - kept.shape[0]),
                "sample_values_after_plot_filter": ";".join(sorted(kept["sample"].astype(str).unique())),
            }
        )
        filtered_parts.append(kept)
    filtered = pd.concat(filtered_parts, ignore_index=True)
    summary = pd.DataFrame(summaries)
    filtered_total = int(summary["n_cells_after_plot_filter"].sum())
    if filtered_total != filtered.shape[0]:
        raise ValueError(f"Internal UMAP plot filter count mismatch: summary={filtered_total}, data={filtered.shape[0]}")
    return filtered, summary


def plot_marker_umap_grid(
    data: pd.DataFrame,
    output_path: str | Path,
    genes: Sequence[str],
    specs: Sequence[CrossStudyMarkerSpec] | None = None,
    title: str = "Cross-study marker expression",
    max_cells_per_study: int | None = None,
    random_state: int = 0,
    point_size: float = 0.36,
    background_point_size: float = 0.15,
    cmap=GREY_TO_RED_CMAP,
    vmax_quantile: float = 0.99,
) -> pd.DataFrame:
    """Plot a study-by-gene UMAP grid from standardized marker tables."""
    genes = list(genes)
    data, filter_summary = _apply_internal_umap_plot_filters(data)
    filter_by_study = filter_summary.set_index("study_label").to_dict(orient="index")
    plot_data = downsample_by_study(data, max_cells_per_study, random_state=random_state)
    study_labels = _ordered_study_labels(plot_data, specs)
    if not study_labels:
        raise ValueError("No study labels available for plotting.")
    missing_genes = [gene for gene in genes if gene not in plot_data.columns]
    if missing_genes:
        raise ValueError("Plot data missing marker columns: " + ", ".join(missing_genes))

    vmax_by_gene = expression_limits(plot_data, genes, quantile=vmax_quantile)
    n_rows = len(study_labels)
    n_cols = len(genes)
    fig_width = max(1.05 * n_cols + 2.2, 7.0)
    fig_height = max(1.65 * n_rows + 1.2, 5.0)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_width, fig_height), squeeze=False)
    manifest_rows = []
    study_subsets: dict[str, pd.DataFrame] = {}
    study_counts: dict[str, int] = {}
    gene_to_col = {gene: idx for idx, gene in enumerate(genes)}
    divider_after_col = None
    if any(gene in ON_TARGET_GENES for gene in genes) and any(gene in OFF_TARGET_GENES for gene in genes):
        on_cols = [gene_to_col[gene] for gene in genes if gene in ON_TARGET_GENES]
        divider_after_col = max(on_cols) if on_cols else None

    for study_label in study_labels:
        subset = plot_data.loc[plot_data["study_label"].astype(str) == study_label].copy()
        x = subset["umap_1"].to_numpy(dtype=float)
        y = subset["umap_2"].to_numpy(dtype=float)
        finite = np.isfinite(x) & np.isfinite(y)
        subset = subset.loc[finite].copy()
        study_subsets[study_label] = subset
        study_counts[study_label] = int(subset.shape[0])

    for row_idx, study_label in enumerate(study_labels):
        subset = study_subsets[study_label]
        x = subset["umap_1"].to_numpy(dtype=float)
        y = subset["umap_2"].to_numpy(dtype=float)
        for col_idx, gene in enumerate(genes):
            ax = axes[row_idx, col_idx]
            vmax = vmax_by_gene[gene]
            norm = Normalize(vmin=0.0, vmax=vmax)
            expr = pd.to_numeric(subset[gene], errors="coerce").to_numpy(dtype=float)
            x = subset["umap_1"].to_numpy(dtype=float)
            y = subset["umap_2"].to_numpy(dtype=float)
            ax.scatter(x, y, s=background_point_size, c="#d0d0d0", linewidths=0, rasterized=True)
            positive = np.isfinite(expr) & (expr > 0)
            if positive.any():
                order = np.argsort(expr[positive])
                ax.scatter(
                    x[positive][order],
                    y[positive][order],
                    s=point_size,
                    c=np.clip(expr[positive][order], 0.0, vmax),
                    cmap=cmap,
                    norm=norm,
                    linewidths=0,
                    rasterized=True,
                )
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_aspect("equal", adjustable="box")
            for spine in ax.spines.values():
                spine.set_visible(False)
            if row_idx == 0:
                ax.set_title(gene, fontsize=8)
            manifest_rows.append(
                {
                    "gene": gene,
                    "gene_group": gene_group(gene),
                    "scale_group": gene,
                    "scale_basis": "per_gene",
                    "study_label": study_label,
                    "plot_filter": filter_by_study[study_label]["plot_filter"],
                    "plot_filter_label": filter_by_study[study_label]["plot_filter_label"],
                    "n_cells_before_plot_filter": filter_by_study[study_label]["n_cells_before_plot_filter"],
                    "n_cells_after_plot_filter": filter_by_study[study_label]["n_cells_after_plot_filter"],
                    "n_cells_removed_by_plot_filter": filter_by_study[study_label]["n_cells_removed_by_plot_filter"],
                    "sample_values_after_plot_filter": filter_by_study[study_label]["sample_values_after_plot_filter"],
                    "n_cells_plotted": study_counts[study_label],
                    "n_positive_cells": int(positive.sum()),
                    "vmax": vmax,
                    "vmax_quantile": vmax_quantile,
                }
            )
    fig.suptitle(title, fontsize=12, y=0.995)
    fig.tight_layout(rect=(0.15, 0.13, 0.995, 0.965))
    fig.canvas.draw()

    for row_idx, study_label in enumerate(study_labels):
        position = axes[row_idx, 0].get_position()
        filter_label = filter_by_study[study_label]["plot_filter_label"]
        row_label = f"{study_label}\n(n = {study_counts[study_label]:,} cells)"
        if filter_label:
            row_label = f"{study_label}\n{filter_label}\n(n = {study_counts[study_label]:,} cells)"
        fig.text(
            position.x0 - 0.012,
            (position.y0 + position.y1) / 2,
            row_label,
            ha="right",
            va="center",
            fontsize=8,
        )

    if divider_after_col is not None and divider_after_col + 1 < n_cols:
        left_position = axes[0, divider_after_col].get_position()
        right_position = axes[0, divider_after_col + 1].get_position()
        x_line = (left_position.x1 + right_position.x0) / 2
        y_top = axes[0, 0].get_position().y1
        y_bottom = axes[-1, 0].get_position().y0
        fig.add_artist(
            plt.Line2D(
                [x_line, x_line],
                [y_bottom, y_top],
                transform=fig.transFigure,
                color="#1f1f1f",
                linewidth=1.2,
            )
        )

    bottom_position = axes[-1, 0].get_position()
    cbar_y = max(0.035, bottom_position.y0 - 0.055)
    for col_idx, gene in enumerate(genes):
        position = axes[-1, col_idx].get_position()
        cax = fig.add_axes([position.x0, cbar_y, position.x1 - position.x0, 0.010])
        vmax = vmax_by_gene[gene]
        sm = ScalarMappable(norm=Normalize(vmin=0.0, vmax=vmax), cmap=cmap)
        sm.set_array([])
        cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
        cbar.ax.tick_params(labelsize=5.5, length=1.5, pad=1)
        cbar.set_ticks([0.0, vmax])
        cbar.set_ticklabels(["0", f"{vmax:.2g}"])
        cbar.outline.set_linewidth(0.4)

    fig.text(
        0.5,
        max(0.012, cbar_y - 0.030),
        "log1p(CP10K); each gene column has its own scale",
        ha="center",
        va="center",
        fontsize=7,
    )

    out = Path(output_path).expanduser()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    if out.suffix.lower() != ".pdf":
        fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    return pd.DataFrame(manifest_rows)


def plot_default_marker_grids(
    project_root: str | Path | None = None,
    run_label: str = "cross_study_marker_expression_v8",
    max_cells_per_study: int | None = None,
    include_xiang: bool = False,
    exclude_study_ids: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Load marker tables and write ON-target, OFF-target, and combined grids."""
    ensure_output_dirs(project_root, run_label)
    specs = default_cross_study_marker_specs(project_root, include_xiang=include_xiang)
    data = load_marker_expression_tables(specs, project_root, run_label)
    expression_summary_table(data).to_csv(
        table_dir(project_root, run_label) / "cross_study_marker_expression_value_summary.tsv",
        sep="\t",
        index=False,
    )
    plot_specs = exclude_specs(specs, exclude_study_ids)
    excluded = {spec.study_id for spec in specs} - {spec.study_id for spec in plot_specs}
    if excluded:
        data = data.loc[~data["study_id"].astype(str).isin(excluded)].copy()
    if data.empty:
        raise ValueError("No cells remain after applying plot study exclusions.")
    _, plot_filter_summary = _apply_internal_umap_plot_filters(data)
    plot_filter_summary.to_csv(
        table_dir(project_root, run_label) / "cross_study_marker_expression_internal_plot_filter_summary.tsv",
        sep="\t",
        index=False,
    )
    outputs = [
        ("on_target", ON_TARGET_GENES, "Cross-study ON-target marker expression"),
        ("off_target", OFF_TARGET_GENES, "Cross-study OFF-target marker expression"),
        ("on_off_target", GENE_ORDER, "Cross-study ON/OFF marker expression"),
    ]
    manifests = []
    for token, genes, title in outputs:
        output_path = plot_dir(project_root, run_label) / f"cross_study_marker_expression_{token}.png"
        manifest = plot_marker_umap_grid(
            data=data,
            output_path=output_path,
            genes=genes,
            specs=plot_specs,
            title=title,
            max_cells_per_study=max_cells_per_study,
        )
        manifest.insert(0, "plot_token", token)
        manifest.insert(1, "plot_path", str(output_path))
        manifest["excluded_study_ids"] = ",".join(sorted(excluded))
        manifests.append(manifest)
    combined = pd.concat(manifests, ignore_index=True)
    combined.to_csv(table_dir(project_root, run_label) / "cross_study_marker_expression_plot_manifest.tsv", sep="\t", index=False)
    return combined


def _parse_study_ids(raw: Sequence[str] | None) -> list[str] | None:
    if not raw:
        return None
    values: list[str] = []
    for item in raw:
        values.extend([part for part in re.split(r"[,;\s]+", item) if part and part.lower() not in {"none", "null", "na"}])
    return values or None


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", default=None)
    parser.add_argument("--run-label", default="cross_study_marker_expression_v8")
    parser.add_argument("--include-xiang", action="store_true")
    subparsers = parser.add_subparsers(dest="command", required=True)

    subparsers.add_parser("setup", help="Write study/gene/readiness setup tables.")

    extract = subparsers.add_parser("extract-h5ad", help="Extract marker tables from configured H5AD caches.")
    extract.add_argument("--study-id", action="append", default=[], help="Study ID to extract. Repeatable or comma-separated.")

    plot = subparsers.add_parser("plot", help="Plot default marker-expression grids from prepared tables.")
    plot.add_argument("--max-cells-per-study", type=int, default=None)
    plot.add_argument("--exclude-study-id", action="append", default=[], help="Study ID to exclude from plots. Repeatable or comma-separated.")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    if args.command == "setup":
        paths = write_setup_tables(args.project_root, args.run_label, include_xiang=args.include_xiang)
        for key, path in paths.items():
            print(f"{key}\t{path}", flush=True)
    elif args.command == "extract-h5ad":
        summary = extract_available_h5ad_marker_tables(
            args.project_root,
            args.run_label,
            study_ids=_parse_study_ids(args.study_id),
            include_xiang=args.include_xiang,
        )
        print(summary.to_string(index=False), flush=True)
    elif args.command == "plot":
        manifest = plot_default_marker_grids(
            args.project_root,
            args.run_label,
            max_cells_per_study=args.max_cells_per_study,
            include_xiang=args.include_xiang,
            exclude_study_ids=_parse_study_ids(args.exclude_study_id),
        )
        print(manifest[["plot_token", "plot_path"]].drop_duplicates().to_string(index=False), flush=True)
    else:
        raise ValueError(f"Unhandled command: {args.command}")


if __name__ == "__main__":
    main()
