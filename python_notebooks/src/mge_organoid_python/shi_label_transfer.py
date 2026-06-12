"""Shi reference label-transfer helpers for standalone notebook workflows.

The kNN label-transfer functions in this module are legacy provenance for the
completed Scanpy-side DIV30 run. Current Shi prediction scores must come from
Seurat TransferData exports, not from neighbor-vote fractions.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import math
import re
from typing import Iterable

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler, normalize

from .shi_prediction_schema import sanitize_shi_label_token


DEFAULT_NON_NEURAL_LABELS = {"Microglia", "OPC", "Endothelial"}
DEFAULT_REGIONAL_LABELS = {"MGE", "LGE", "CGE"}


@dataclass(frozen=True)
class LabelTransferConfig:
    """Runtime settings for one Shi label-transfer run."""

    query_h5ad: Path
    reference_h5ad: Path
    shi_table_s2_xlsx: Path
    run_dir: Path
    n_top_variable_genes: int = 3000
    n_pcs: int = 50
    n_neighbors: int = 31
    random_state: int = 0
    n_jobs: int = 1
    min_reference_label_match_fraction: float = 0.8
    write_h5ad: bool = True
    save_plots: bool = True


def safe_token(value: object) -> str:
    """Return a stable token suitable for filenames and column names."""
    return sanitize_shi_label_token(value)


def natural_sort_key(value: object) -> list[object]:
    """Sort cluster-like strings with embedded numbers in human order."""
    parts = re.split(r"(\d+)", str(value))
    return [int(part) if part.isdigit() else part.lower() for part in parts]


def parse_csv_set(raw: str | None) -> set[str]:
    """Parse comma/semicolon/colon-delimited values into a string set."""
    if raw is None or str(raw).strip() == "":
        return set()
    return {part.strip() for part in re.split(r"[,;:]+", str(raw)) if part.strip()}


def normalize_shi_barcode(value: object) -> str:
    """Normalize Shi cell IDs so Table S2 and cached AnnData cells can join."""
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    if not text:
        return ""
    text = re.sub(r"-GW\d+(?:_[A-Za-z0-9]+)?$", "", text)
    text = re.sub(r"-1$", "", text)
    return text


def shi_barcode_base(value: object) -> str:
    """Return the 10x barcode base before any sample suffix."""
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    if not text:
        return ""
    match = re.match(r"^([ACGT]+)-", text)
    if match:
        return match.group(1)
    return normalize_shi_barcode(text)


def shi_reference_suffix(value: object) -> str:
    """Return the cached Shi suffix after the 10x barcode base."""
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    match = re.match(r"^[ACGT]+-(.+)$", text)
    return match.group(1) if match else ""


def shi_table_numeric_suffix(value: object) -> str:
    """Return the numeric Cell Ranger suffix from a Table S2 cell ID."""
    if value is None or pd.isna(value):
        return ""
    match = re.search(r"-(\d+)$", str(value).strip())
    return match.group(1) if match else ""


def load_shi_table_s2_labels(path: str | Path) -> pd.DataFrame:
    """Read Shi Table S2 cell-level major-type labels from the workbook."""
    import openpyxl

    table_path = Path(path).expanduser()
    if not table_path.exists():
        raise FileNotFoundError(f"Missing Shi Table S2 workbook: {table_path}")

    workbook = openpyxl.load_workbook(table_path, read_only=True, data_only=True)
    try:
        worksheet = workbook["UMAP coordinates"]
    except KeyError:
        worksheet = workbook.worksheets[0]

    rows = worksheet.iter_rows(values_only=True)
    try:
        next(rows)
        headers = [str(value).strip() if value is not None else "" for value in next(rows)]
    except StopIteration as exc:
        raise ValueError(f"Shi Table S2 appears empty: {table_path}") from exc

    if "Cells" not in headers or "Major types" not in headers:
        raise ValueError(
            "Shi Table S2 must contain 'Cells' and 'Major types' columns; found "
            + ", ".join(headers)
        )

    cell_idx = headers.index("Cells")
    label_idx = headers.index("Major types")
    records = []
    for row in rows:
        if row is None:
            continue
        cell = row[cell_idx] if cell_idx < len(row) else None
        label = row[label_idx] if label_idx < len(row) else None
        if cell is None or label is None:
            continue
        cell_text = str(cell).strip()
        label_text = str(label).strip()
        barcode = shi_barcode_base(cell_text)
        numeric_suffix = shi_table_numeric_suffix(cell_text)
        cell_key = f"{barcode}-{numeric_suffix}" if barcode and numeric_suffix else ""
        if barcode and numeric_suffix and label_text:
            records.append(
                {
                    "shi_table_s2_cell": cell_text,
                    "shi_table_s2_barcode": barcode,
                    "shi_table_s2_numeric_suffix": numeric_suffix,
                    "shi_table_s2_cell_key": cell_key,
                    "shi_label": label_text,
                }
            )
    workbook.close()

    labels = pd.DataFrame.from_records(records)
    if labels.empty:
        raise ValueError(f"No usable cell labels found in Shi Table S2: {table_path}")
    labels = labels.drop_duplicates("shi_table_s2_cell_key", keep="first").reset_index(drop=True)
    return labels


def _candidate_reference_keys(obs: pd.DataFrame) -> pd.DataFrame:
    candidates = pd.DataFrame(index=obs.index)
    candidates["obs_name"] = [normalize_shi_barcode(value) for value in obs.index]
    for column in ["barcode", "raw_cell_id", "cell_id"]:
        if column in obs.columns:
            candidates[column] = [normalize_shi_barcode(value) for value in obs[column]]
    return candidates


def _reference_cell_id_series(obs: pd.DataFrame) -> pd.Series:
    if "cell_id" in obs.columns:
        return obs["cell_id"].astype(str)
    if "raw_cell_id" in obs.columns:
        return obs["raw_cell_id"].astype(str)
    return pd.Series(obs.index.astype(str), index=obs.index)


def infer_reference_to_table_suffix_map(
    reference_obs: pd.DataFrame,
    table_s2_labels: pd.DataFrame,
) -> pd.DataFrame:
    """Infer cached Shi suffix -> Table S2 numeric suffix using barcode overlap."""
    ref_cell_ids = _reference_cell_id_series(reference_obs)
    ref = pd.DataFrame(
        {
            "reference_suffix": [shi_reference_suffix(value) for value in ref_cell_ids],
            "barcode_base": [shi_barcode_base(value) for value in ref_cell_ids],
        },
        index=reference_obs.index,
    )
    table = table_s2_labels[["shi_table_s2_barcode", "shi_table_s2_numeric_suffix"]].drop_duplicates()
    table_suffixes_by_barcode = (
        table.groupby("shi_table_s2_barcode")["shi_table_s2_numeric_suffix"]
        .apply(lambda values: sorted(set(values.dropna().astype(str))))
        .to_dict()
    )

    records = []
    for reference_suffix, sub in ref.groupby("reference_suffix", observed=True, dropna=False):
        counts: dict[str, int] = {}
        for barcode in sub["barcode_base"].astype(str):
            for suffix in table_suffixes_by_barcode.get(barcode, []):
                counts[suffix] = counts.get(suffix, 0) + 1
        ranked = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
        best_suffix = ranked[0][0] if ranked else ""
        best_count = ranked[0][1] if ranked else 0
        second_count = ranked[1][1] if len(ranked) > 1 else 0
        unambiguous = bool(best_suffix) and best_count > second_count
        records.append(
            {
                "reference_suffix": reference_suffix,
                "n_reference_cells": int(sub.shape[0]),
                "inferred_table_suffix": best_suffix if unambiguous else "",
                "best_table_suffix": best_suffix,
                "best_suffix_overlap_count": int(best_count),
                "second_suffix_overlap_count": int(second_count),
                "suffix_mapping_unambiguous": unambiguous,
            }
        )
    return pd.DataFrame.from_records(records)


def attach_reference_labels(
    reference: ad.AnnData,
    table_s2_labels: pd.DataFrame,
    min_match_fraction: float = 0.8,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Join Shi Table S2 major-type labels to reference AnnData obs."""
    suffix_map = infer_reference_to_table_suffix_map(reference.obs, table_s2_labels)
    suffix_to_table = suffix_map.set_index("reference_suffix")["inferred_table_suffix"].to_dict()
    label_map = table_s2_labels.set_index("shi_table_s2_cell_key")["shi_label"].to_dict()
    table_cell_map = table_s2_labels.set_index("shi_table_s2_cell_key")["shi_table_s2_cell"].to_dict()
    keys_by_barcode = (
        table_s2_labels.groupby("shi_table_s2_barcode")["shi_table_s2_cell_key"]
        .apply(lambda values: sorted(set(values.dropna().astype(str))))
        .to_dict()
    )
    barcode_label_map = (
        table_s2_labels[table_s2_labels["shi_table_s2_numeric_suffix"] == "1"]
        .set_index("shi_table_s2_barcode")["shi_label"]
        .to_dict()
    )
    barcode_cell_map = (
        table_s2_labels[table_s2_labels["shi_table_s2_numeric_suffix"] == "1"]
        .set_index("shi_table_s2_barcode")["shi_table_s2_cell"]
        .to_dict()
    )
    candidates = _candidate_reference_keys(reference.obs)
    ref_cell_ids = _reference_cell_id_series(reference.obs)

    labels = []
    join_keys = []
    join_sources = []
    table_cells = []
    reference_suffixes = []
    inferred_table_suffixes = []
    for idx, row in candidates.iterrows():
        chosen_label = ""
        chosen_key = ""
        chosen_source = ""
        reference_suffix = shi_reference_suffix(ref_cell_ids.loc[idx])
        barcode_base = shi_barcode_base(ref_cell_ids.loc[idx])
        inferred_table_suffix = suffix_to_table.get(reference_suffix, "")
        barcode_keys = keys_by_barcode.get(barcode_base, [])
        if len(barcode_keys) == 1:
            table_key = barcode_keys[0]
            if table_key in label_map:
                chosen_label = label_map[table_key]
                chosen_key = table_key
                chosen_source = "unique_table_s2_barcode_key"
        if not chosen_label and barcode_base and inferred_table_suffix:
            table_key = f"{barcode_base}-{inferred_table_suffix}"
            if table_key in label_map:
                chosen_label = label_map[table_key]
                chosen_key = table_key
                chosen_source = "inferred_table_s2_cell_key"
        for source, value in row.items():
            if chosen_label:
                break
            if value and value in label_map:
                chosen_label = label_map[value]
                chosen_key = value
                chosen_source = source
                break
            if value and value in barcode_label_map:
                chosen_label = barcode_label_map[value]
                chosen_key = value
                chosen_source = f"{source}_table_suffix_1_fallback"
                break
        labels.append(chosen_label)
        join_keys.append(chosen_key)
        join_sources.append(chosen_source)
        table_cells.append(table_cell_map.get(chosen_key, barcode_cell_map.get(chosen_key, "")))
        reference_suffixes.append(reference_suffix)
        inferred_table_suffixes.append(inferred_table_suffix)

    label_obs = pd.DataFrame(
        {
            "reference_obs_name": reference.obs_names.astype(str),
            "reference_suffix": reference_suffixes,
            "inferred_table_suffix": inferred_table_suffixes,
            "shi_table_s2_join_key": join_keys,
            "shi_table_s2_join_source": join_sources,
            "shi_table_s2_cell": table_cells,
            "shi_label": labels,
        },
        index=reference.obs_names,
    )
    label_obs["shi_label"] = label_obs["shi_label"].replace("", pd.NA)

    matched = int(label_obs["shi_label"].notna().sum())
    total = int(label_obs.shape[0])
    fraction = matched / total if total else 0.0
    summary = pd.DataFrame(
        [
            {"metric": "reference_cells", "value": total},
            {"metric": "table_s2_label_rows", "value": int(table_s2_labels.shape[0])},
            {"metric": "matched_reference_cells", "value": matched},
            {"metric": "unmatched_reference_cells", "value": total - matched},
            {"metric": "matched_reference_fraction", "value": fraction},
            {"metric": "min_required_match_fraction", "value": float(min_match_fraction)},
            {"metric": "unambiguous_suffix_mappings", "value": int(suffix_map["suffix_mapping_unambiguous"].sum())},
            {"metric": "total_reference_suffixes", "value": int(suffix_map.shape[0])},
        ]
    )
    if fraction < min_match_fraction:
        raise ValueError(
            "Only {:.1%} of Shi reference cells matched Table S2 labels; expected at least {:.1%}".format(
                fraction, min_match_fraction
            )
        )
    return label_obs, summary


def broad_region_from_label(label: object) -> str:
    """Map a Shi label to MGE/LGE/CGE/other."""
    text = str(label).strip().lower()
    if "mge" in text:
        return "MGE"
    if "lge" in text:
        return "LGE"
    if "cge" in text:
        return "CGE"
    return "other"


def developmental_class_from_label(label: object) -> str:
    """Map a Shi label to a broad developmental class."""
    text = str(label).strip().lower()
    if "ipc" in text:
        return "IPC"
    if "rgc" in text or "radial" in text or "progenitor" in text:
        return "RGC"
    if "neuron" in text or "mge" in text or "lge" in text or "cge" in text:
        return "neuron"
    return "other"


def comparison_label_sets(
    labels: Iterable[str],
    non_neural_labels: set[str] | None = None,
    regional_labels: set[str] | None = None,
) -> dict[str, set[str]]:
    """Return full-developmental and restricted regional reference label sets."""
    observed = {str(label).strip() for label in labels if str(label).strip()}
    excluded = non_neural_labels if non_neural_labels is not None else set()
    regional = regional_labels if regional_labels is not None else DEFAULT_REGIONAL_LABELS
    return {
        "full_relevant": {label for label in observed if label not in excluded},
        "mge_lge_cge_only": {label for label in observed if label in regional},
    }


def shi_week_label_from_cell_id(value: object) -> str:
    """Infer a Shi gestational-week label from cached cell IDs."""
    suffix = shi_reference_suffix(value)
    match = re.search(r"GW(\d+)", suffix)
    if not match:
        return ""
    week = int(match.group(1))
    return f"GW{week:02d}"


def shi_week_numeric_from_label(value: object) -> float:
    """Return numeric gestational week from a Shi week label or cell suffix."""
    if value is None or pd.isna(value):
        return np.nan
    match = re.search(r"GW\s*0*(\d+)", str(value))
    return float(match.group(1)) if match else np.nan


def reference_week_metadata(reference_obs: pd.DataFrame) -> pd.DataFrame:
    """Return filled Shi reference week labels/numeric values."""
    cell_ids = _reference_cell_id_series(reference_obs)
    inferred_label = pd.Series(
        [shi_week_label_from_cell_id(value) for value in cell_ids],
        index=reference_obs.index,
        dtype="string",
    )
    if "week_label" in reference_obs.columns:
        week_label = reference_obs["week_label"].astype("string").replace("", pd.NA)
        week_label = week_label.fillna(inferred_label.replace("", pd.NA))
    else:
        week_label = inferred_label.replace("", pd.NA)

    if "week_numeric" in reference_obs.columns:
        week_numeric = pd.to_numeric(reference_obs["week_numeric"], errors="coerce")
        inferred_numeric = week_label.map(shi_week_numeric_from_label)
        week_numeric = week_numeric.fillna(inferred_numeric)
    else:
        week_numeric = week_label.map(shi_week_numeric_from_label)

    return pd.DataFrame(
        {
            "reference_obs_name": reference_obs.index.astype(str),
            "shi_week_label": week_label.astype("string"),
            "shi_week_numeric": week_numeric.astype(float),
        },
        index=reference_obs.index,
    )


def match_common_genes(
    query_var_names: Iterable[object],
    reference_var_names: Iterable[object],
) -> tuple[list[int], list[int], pd.DataFrame]:
    """Match query/reference gene names with exact or unambiguous uppercase matching."""
    query_values = [str(value) for value in query_var_names]
    ref_values = [str(value) for value in reference_var_names]

    ref_exact = {value: idx for idx, value in enumerate(ref_values)}
    ref_upper = {}
    ambiguous = set()
    for idx, value in enumerate(ref_values):
        key = value.upper()
        if key in ref_upper and ref_upper[key][1] != value:
            ambiguous.add(key)
        else:
            ref_upper[key] = (idx, value)

    records = []
    query_indices = []
    reference_indices = []
    used_ref = set()
    for q_idx, q_gene in enumerate(query_values):
        r_idx = None
        r_gene = ""
        match_type = "missing"
        if q_gene in ref_exact:
            r_idx = ref_exact[q_gene]
            r_gene = q_gene
            match_type = "exact"
        else:
            key = q_gene.upper()
            if key in ambiguous:
                match_type = "ambiguous_case_insensitive"
            elif key in ref_upper:
                r_idx, r_gene = ref_upper[key]
                match_type = "case_insensitive"

        used = r_idx is not None and r_idx not in used_ref
        if used:
            query_indices.append(q_idx)
            reference_indices.append(r_idx)
            used_ref.add(r_idx)
        records.append(
            {
                "query_gene": q_gene,
                "reference_gene": r_gene,
                "match_type": match_type,
                "used_for_transfer": used,
            }
        )

    detail = pd.DataFrame.from_records(records)
    if not query_indices:
        raise ValueError("No common genes found between query and Shi reference")
    return query_indices, reference_indices, detail


def _as_csr_float32(matrix) -> sparse.csr_matrix:
    if sparse.issparse(matrix):
        return matrix.tocsr().astype(np.float32)
    return sparse.csr_matrix(np.asarray(matrix, dtype=np.float32))


def subset_expression_matrix(adata_obj: ad.AnnData, row_mask, gene_indices: list[int]) -> sparse.csr_matrix:
    """Return a CSR float32 expression matrix for selected cells and genes."""
    subset = adata_obj[row_mask, gene_indices].X
    return _as_csr_float32(subset)


def select_variable_gene_indices(matrix: sparse.csr_matrix, n_top: int) -> tuple[np.ndarray, pd.DataFrame]:
    """Pick top variable genes from a sparse log-expression matrix."""
    if n_top <= 0 or n_top >= matrix.shape[1]:
        idx = np.arange(matrix.shape[1], dtype=int)
    else:
        mean = np.asarray(matrix.mean(axis=0)).ravel()
        mean_sq = np.asarray(matrix.multiply(matrix).mean(axis=0)).ravel()
        variance = mean_sq - np.square(mean)
        variance[~np.isfinite(variance)] = -np.inf
        idx = np.argsort(variance)[::-1][:n_top]
    stats = pd.DataFrame(
        {
            "harmonized_gene_index": np.arange(matrix.shape[1], dtype=int),
            "selected_for_transfer": np.isin(np.arange(matrix.shape[1]), idx),
        }
    )
    return np.sort(idx), stats


def _fit_project_reference_query(
    reference_matrix: sparse.csr_matrix,
    query_matrix: sparse.csr_matrix,
    n_pcs: int,
    random_state: int,
) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    max_components = min(n_pcs, reference_matrix.shape[0] - 1, reference_matrix.shape[1] - 1)
    if max_components < 2:
        raise ValueError("Need at least two cells and genes for transfer PCA/SVD")

    scaler = StandardScaler(with_mean=False, copy=True)
    reference_scaled = scaler.fit_transform(reference_matrix)
    query_scaled = scaler.transform(query_matrix)

    svd = TruncatedSVD(n_components=max_components, random_state=random_state)
    reference_embedding = svd.fit_transform(reference_scaled)
    query_embedding = svd.transform(query_scaled)
    reference_embedding = normalize(reference_embedding)
    query_embedding = normalize(query_embedding)

    variance = pd.DataFrame(
        {
            "component": np.arange(1, max_components + 1),
            "explained_variance_ratio": svd.explained_variance_ratio_,
        }
    )
    return reference_embedding, query_embedding, variance


def _label_vote_predictions(
    neighbor_indices: np.ndarray,
    reference_labels: np.ndarray,
    output_prefix: str = "",
) -> pd.DataFrame:
    labels, inverse = np.unique(reference_labels.astype(str), return_inverse=True)
    neighbor_codes = inverse[neighbor_indices]
    n_labels = len(labels)
    norm_base = math.log(n_labels) if n_labels > 1 else 1.0

    pred_labels = []
    scores = []
    entropy_values = []
    for codes in neighbor_codes:
        counts = np.bincount(codes, minlength=n_labels).astype(float)
        probs = counts / counts.sum()
        top_idx = int(np.argmax(probs))
        nonzero = probs[probs > 0]
        entropy = -float(np.sum(nonzero * np.log(nonzero))) / norm_base if norm_base else 0.0
        pred_labels.append(labels[top_idx])
        scores.append(float(probs[top_idx]))
        entropy_values.append(entropy)

    predictions = pd.DataFrame(
        {
            f"{output_prefix}predicted_shi_label": pred_labels,
            f"{output_prefix}prediction_score": scores,
            f"{output_prefix}uncertainty_score": [1.0 - value for value in scores],
            f"{output_prefix}prediction_entropy": entropy_values,
        }
    )
    return predictions


def run_knn_label_transfer(
    query: ad.AnnData,
    reference: ad.AnnData,
    reference_label_obs: pd.DataFrame,
    query_gene_indices: list[int],
    reference_gene_indices: list[int],
    allowed_labels: set[str],
    comparison_name: str,
    n_top_variable_genes: int = 3000,
    n_pcs: int = 50,
    n_neighbors: int = 31,
    random_state: int = 0,
    n_jobs: int = 1,
    reference_week_obs: pd.DataFrame | None = None,
    allow_legacy_knn_scoring: bool = False,
) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    """Run one reference-to-query kNN label-transfer comparison."""
    if not allow_legacy_knn_scoring:
        raise RuntimeError(
            "run_knn_label_transfer() uses legacy Scanpy-side neighbor-vote "
            "prediction scores. Current Shi prediction scoring must use Seurat "
            "TransferData exports: prediction.score.max and prediction.score.<label>. "
            "Pass allow_legacy_knn_scoring=True only to reproduce the completed "
            "shi_reference_div30_label_transfer_v2 provenance run."
        )

    label_series = reference_label_obs["shi_label"].astype("string")
    ref_mask = label_series.isin(sorted(allowed_labels)).to_numpy()
    if int(ref_mask.sum()) < 10:
        raise ValueError(f"{comparison_name}: fewer than 10 reference cells after label filtering")

    reference_matrix = subset_expression_matrix(reference, ref_mask, reference_gene_indices)
    query_matrix = subset_expression_matrix(query, slice(None), query_gene_indices)

    selected_idx, variable_gene_stats = select_variable_gene_indices(reference_matrix, n_top_variable_genes)
    reference_matrix = reference_matrix[:, selected_idx]
    query_matrix = query_matrix[:, selected_idx]

    ref_embedding, query_embedding, svd_variance = _fit_project_reference_query(
        reference_matrix,
        query_matrix,
        n_pcs=n_pcs,
        random_state=random_state,
    )

    k = min(int(n_neighbors), ref_embedding.shape[0])
    nn = NearestNeighbors(n_neighbors=k, metric="cosine", n_jobs=n_jobs)
    nn.fit(ref_embedding)
    distances, indices = nn.kneighbors(query_embedding, return_distance=True)

    ref_labels = label_series[ref_mask].astype(str).to_numpy()
    predictions = _label_vote_predictions(indices, ref_labels)
    predictions["broad_region_class"] = predictions["predicted_shi_label"].map(broad_region_from_label)
    predictions["developmental_class"] = predictions["predicted_shi_label"].map(developmental_class_from_label)
    predictions.index = query.obs_names
    predictions.insert(0, "obs_name", query.obs_names.astype(str))
    predictions.insert(1, "comparison", comparison_name)
    predictions["n_neighbors_used"] = k
    predictions["mean_neighbor_distance"] = distances.mean(axis=1)

    if reference_week_obs is not None:
        week_labels = reference_week_obs.loc[reference_label_obs.index, "shi_week_label"].astype("string")
        week_labels = week_labels[ref_mask].fillna("unknown").astype(str).to_numpy()
        week_predictions = _label_vote_predictions(
            indices,
            week_labels,
            output_prefix="week_",
        )
        week_predictions = week_predictions.rename(
            columns={
                "week_predicted_shi_label": "predicted_shi_week_label",
                "week_prediction_score": "week_prediction_score",
                "week_uncertainty_score": "week_uncertainty_score",
                "week_prediction_entropy": "week_prediction_entropy",
            }
        )
        week_predictions.index = predictions.index
        week_numeric = pd.to_numeric(
            reference_week_obs.loc[reference_label_obs.index, "shi_week_numeric"],
            errors="coerce",
        )[ref_mask].to_numpy(dtype=float)
        neighbor_week = week_numeric[indices]
        predictions = pd.concat([predictions, week_predictions], axis=1)
        predictions["mean_neighbor_week_numeric"] = np.nanmean(neighbor_week, axis=1)
        predictions["median_neighbor_week_numeric"] = np.nanmedian(neighbor_week, axis=1)
        predictions["std_neighbor_week_numeric"] = np.nanstd(neighbor_week, axis=1)

    ref_label_counts = (
        pd.Series(ref_labels, name="shi_label")
        .value_counts()
        .rename_axis("shi_label")
        .reset_index(name="n_reference_cells")
    )
    ref_label_counts["comparison"] = comparison_name

    tables = {
        "reference_label_counts": ref_label_counts,
        "variable_gene_stats": variable_gene_stats,
        "svd_variance": svd_variance,
    }
    return predictions, tables


def normalized_entropy_from_counts(counts: pd.Series) -> float:
    total = float(counts.sum())
    if total <= 0:
        return np.nan
    probs = counts.astype(float) / total
    probs = probs[probs > 0]
    if probs.shape[0] <= 1:
        return 0.0
    return -float((probs * np.log(probs)).sum()) / math.log(probs.shape[0])


def summarize_predictions_by_cluster(
    obs: pd.DataFrame,
    label_col: str,
    score_col: str,
    uncertainty_col: str,
    cluster_col: str = "seurat_clusters",
    label_categories: Iterable[str] | None = None,
) -> dict[str, pd.DataFrame]:
    """Summarize transferred labels by existing query Seurat clusters."""
    if cluster_col not in obs.columns:
        raise ValueError(f"Missing post hoc cluster column: {cluster_col}")
    work = obs[[cluster_col, label_col, score_col, uncertainty_col]].copy()
    work[cluster_col] = work[cluster_col].astype(str)
    work[label_col] = work[label_col].astype(str)

    cluster_counts = work.groupby(cluster_col, observed=True).size().reset_index(name="n_cells")
    cluster_values = sorted(cluster_counts[cluster_col].astype(str), key=natural_sort_key)
    if label_categories is None:
        label_values = sorted(work[label_col].dropna().astype(str).unique())
    else:
        label_values = [str(value) for value in label_categories]

    label_counts = (
        work.groupby([cluster_col, label_col], observed=True)
        .size()
        .reset_index(name="n_cells")
    )
    full_index = pd.MultiIndex.from_product(
        [cluster_values, label_values],
        names=[cluster_col, label_col],
    )
    label_counts = (
        label_counts.set_index([cluster_col, label_col])
        .reindex(full_index, fill_value=0)
        .reset_index()
    )
    cluster_size_map = cluster_counts.set_index(cluster_col)["n_cells"].to_dict()
    label_counts["n_cells_cluster"] = label_counts[cluster_col].map(cluster_size_map).astype(int)
    label_counts["fraction"] = label_counts["n_cells"] / label_counts["n_cells_cluster"]

    dominant = (
        label_counts.sort_values([cluster_col, "n_cells", "fraction"], ascending=[True, False, False])
        .groupby(cluster_col, observed=True)
        .head(1)
        .rename(columns={label_col: "dominant_label", "fraction": "dominant_fraction"})
    )
    entropy = (
        label_counts.groupby(cluster_col, observed=True)["n_cells"]
        .apply(normalized_entropy_from_counts)
        .reset_index(name="label_entropy")
    )
    confidence = (
        work.groupby(cluster_col, observed=True)
        .agg(mean_prediction_score=(score_col, "mean"), mean_uncertainty=(uncertainty_col, "mean"))
        .reset_index()
    )
    summary = (
        cluster_counts.merge(dominant[[cluster_col, "dominant_label", "dominant_fraction"]], on=cluster_col)
        .merge(entropy, on=cluster_col)
        .merge(confidence, on=cluster_col)
    )
    return {
        "cluster_counts": cluster_counts,
        "label_counts": label_counts,
        "cluster_summary": summary,
    }


def get_umap_coordinates(adata_obj: ad.AnnData) -> tuple[np.ndarray, str]:
    """Return query UMAP coordinates and the obsm key used."""
    for key in ["X_umap_seurat", "X_umap"]:
        if key in adata_obj.obsm:
            return np.asarray(adata_obj.obsm[key]), key
    raise ValueError("Query AnnData lacks .obsm['X_umap_seurat'] or .obsm['X_umap']")


def plot_umap_categorical(
    coords: np.ndarray,
    values: pd.Series,
    title: str,
    path: str | Path,
    point_size: float = 1.0,
    category_order: Iterable[str] | None = None,
) -> None:
    if category_order is None:
        categories = pd.Series(values, dtype="category")
    else:
        categories = pd.Series(
            pd.Categorical(values.astype(str), categories=[str(value) for value in category_order])
        )
    codes = categories.cat.codes.to_numpy()
    labels = list(categories.cat.categories)
    cmap = plt.get_cmap("tab20", max(len(labels), 1))
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    ax.scatter(
        coords[:, 0],
        coords[:, 1],
        c=codes,
        s=point_size,
        cmap=cmap,
        vmin=0,
        vmax=max(len(labels) - 1, 0),
        linewidths=0,
        rasterized=True,
    )
    ax.set_title(title)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_xticks([])
    ax.set_yticks([])
    handles = [
        plt.Line2D([0], [0], marker="o", linestyle="", color=cmap(i), markersize=5, label=label)
        for i, label in enumerate(labels)
    ]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, title=values.name)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_overlaid_density_by_group(
    data: pd.DataFrame,
    value_col: str,
    group_col: str,
    title: str,
    path: str | Path,
    bins: int | Iterable[float] = 40,
    x_label: str | None = None,
) -> None:
    """Plot overlaid density histograms for one value column by group."""
    fig, ax = plt.subplots(figsize=(8.0, 5.5))
    groups = sorted(data[group_col].dropna().astype(str).unique(), key=natural_sort_key)
    cmap = plt.get_cmap("tab10", max(len(groups), 1))
    for idx, group in enumerate(groups):
        values = pd.to_numeric(
            data.loc[data[group_col].astype(str) == group, value_col],
            errors="coerce",
        ).dropna()
        if values.empty:
            continue
        ax.hist(
            values,
            bins=bins,
            density=True,
            histtype="step",
            linewidth=1.8,
            color=cmap(idx),
            label=f"{group} (n={values.shape[0]})",
        )
    ax.set_title(title)
    ax.set_xlabel(x_label or value_col)
    ax.set_ylabel("density")
    ax.legend(frameon=False, loc="center left", bbox_to_anchor=(1.02, 0.5), title=group_col)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def summarize_week_predictions_by_sample(
    obs: pd.DataFrame,
    week_label_col: str,
    week_score_col: str,
    week_numeric_col: str,
    sample_col: str = "orig.ident",
) -> pd.DataFrame:
    """Summarize transferred Shi gestational week predictions by query sample."""
    if sample_col not in obs.columns:
        raise ValueError(f"Missing sample column: {sample_col}")
    work = obs[[sample_col, week_label_col, week_score_col, week_numeric_col]].copy()
    work[sample_col] = work[sample_col].astype(str)
    work[week_label_col] = work[week_label_col].astype(str)
    counts = (
        work.groupby([sample_col, week_label_col], observed=True)
        .size()
        .reset_index(name="n_cells")
    )
    total = work.groupby(sample_col, observed=True).size().rename("n_cells_total")
    counts = counts.merge(total.reset_index(), on=sample_col)
    counts["fraction"] = counts["n_cells"] / counts["n_cells_total"]
    dominant = (
        counts.sort_values([sample_col, "n_cells", "fraction"], ascending=[True, False, False])
        .groupby(sample_col, observed=True)
        .head(1)
        .rename(columns={week_label_col: "dominant_week_label", "fraction": "dominant_week_fraction"})
    )
    numeric = (
        work.groupby(sample_col, observed=True)
        .agg(
            n_cells=(week_label_col, "size"),
            mean_week_prediction_score=(week_score_col, "mean"),
            mean_neighbor_week_numeric=(week_numeric_col, "mean"),
            median_neighbor_week_numeric=(week_numeric_col, "median"),
        )
        .reset_index()
    )
    return numeric.merge(
        dominant[[sample_col, "dominant_week_label", "dominant_week_fraction"]],
        on=sample_col,
    )


def plot_umap_continuous(
    coords: np.ndarray,
    values: pd.Series,
    title: str,
    path: str | Path,
    point_size: float = 1.0,
    cmap: str = "viridis",
) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    artist = ax.scatter(coords[:, 0], coords[:, 1], c=values.astype(float), s=point_size, cmap=cmap, linewidths=0, rasterized=True)
    ax.set_title(title)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_xticks([])
    ax.set_yticks([])
    fig.colorbar(artist, ax=ax, label=values.name)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_stacked_bar(
    label_counts: pd.DataFrame,
    cluster_col: str,
    label_col: str,
    path: str | Path,
    title: str,
) -> None:
    wide = label_counts.pivot(index=cluster_col, columns=label_col, values="fraction").fillna(0.0)
    wide = wide.loc[sorted(wide.index, key=natural_sort_key)]
    fig, ax = plt.subplots(figsize=(max(8, 0.45 * wide.shape[0]), 5.5))
    wide.plot(kind="bar", stacked=True, ax=ax, width=0.85, colormap="tab20")
    ax.set_title(title)
    ax.set_xlabel("seurat_clusters")
    ax.set_ylabel("fraction of cells")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, title=label_col)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_heatmap(
    label_counts: pd.DataFrame,
    cluster_col: str,
    label_col: str,
    path: str | Path,
    title: str,
) -> None:
    wide = label_counts.pivot(index=cluster_col, columns=label_col, values="fraction").fillna(0.0)
    wide = wide.loc[sorted(wide.index, key=natural_sort_key)]
    fig, ax = plt.subplots(figsize=(max(7, 0.45 * wide.shape[1]), max(5, 0.3 * wide.shape[0])))
    im = ax.imshow(wide.to_numpy(), aspect="auto", vmin=0.0, vmax=1.0, cmap="magma")
    ax.set_title(title)
    ax.set_xlabel("Shi predicted label")
    ax.set_ylabel("seurat_clusters")
    ax.set_xticks(np.arange(wide.shape[1]))
    ax.set_xticklabels(wide.columns, rotation=45, ha="right")
    ax.set_yticks(np.arange(wide.shape[0]))
    ax.set_yticklabels(wide.index)
    fig.colorbar(im, ax=ax, label="fraction")
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
