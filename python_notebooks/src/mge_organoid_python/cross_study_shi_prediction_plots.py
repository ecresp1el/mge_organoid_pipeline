"""Cross-study Shi Seurat label-transfer table preparation and plotting."""

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
from matplotlib.cm import ScalarMappable
from matplotlib.colors import BoundaryNorm, ListedColormap, Normalize
import numpy as np
import pandas as pd

from .paths import resolve_project_root


RUN_LABEL_DEFAULT = "cross_study_shi_seurat_label_transfer_v1"
RESULTS_DIRNAME = "cross_study_shi_seurat_label_transfer"

STUDY_ORDER = [
    "varela_div30",
    "varela_div90",
    "siebert_2026",
    "walsh",
    "bershteyn_2025",
    "bershteyn_2023",
    "samarasinghe_2021",
]

STUDY_LABELS = {
    "varela_div30": "This Study, DIV 30",
    "varela_div90": "This Study, DIV 90",
    "siebert_2026": "Siebert et al. 2026",
    "walsh": "Walsh et al. 2025",
    "bershteyn_2025": "Bershteyn et al. 2025",
    "bershteyn_2023": "Bershteyn et al. 2023",
    "samarasinghe_2021": "Samarasinghe et al. 2021",
}

VARELA_H5AD_PATHS = {
    "varela_div30": "results/python_anndata/varela_div30.h5ad",
    "varela_div90": "results/python_anndata/varela_div90.h5ad",
}

SHI_LABEL_ORDER = [
    "MGE",
    "LGE",
    "CGE",
    "progenitor",
    "Excitatory IPC",
    "Excitatory neuron",
    "Thalamic neurons",
    "Microglia",
    "OPC",
    "Endothelial",
]

LABEL_COLORS = {
    "MGE": "#16697a",
    "LGE": "#4895ef",
    "CGE": "#52b788",
    "progenitor": "#f4a261",
    "Excitatory IPC": "#e76f51",
    "Excitatory neuron": "#d62828",
    "Thalamic neurons": "#7b2cbf",
    "Microglia": "#6c757d",
    "OPC": "#a7c957",
    "Endothelial": "#2a9d8f",
}

LABEL_SCORE_MAP = {
    "MGE": "shi_seurat_full_prediction_score_MGE",
    "LGE": "shi_seurat_full_prediction_score_LGE",
    "CGE": "shi_seurat_full_prediction_score_CGE",
    "progenitor": "shi_seurat_full_prediction_score_progenitor",
    "Excitatory IPC": "shi_seurat_full_prediction_score_Excitatory_IPC",
    "Excitatory neuron": "shi_seurat_full_prediction_score_Excitatory_neuron",
    "Thalamic neurons": "shi_seurat_full_prediction_score_Thalamic_neurons",
    "Microglia": "shi_seurat_full_prediction_score_Microglia",
    "OPC": "shi_seurat_full_prediction_score_OPC",
    "Endothelial": "shi_seurat_full_prediction_score_Endothelial",
}

BASE_COLUMNS = [
    "cell_id",
    "study_id",
    "study_label",
    "sample",
    "sample_label",
    "cluster",
    "umap_1",
    "umap_2",
    "shi_seurat_full_predicted_shi_label",
    "shi_seurat_full_prediction_score",
    "shi_seurat_full_uncertainty_score",
    "shi_seurat_full_predicted_shi_week_label",
    "shi_seurat_full_week_prediction_score",
    "shi_seurat_full_week_uncertainty_score",
    "shi_seurat_full_expected_shi_gw_numeric",
    "shi_seurat_full_expected_shi_gw_even",
]


@dataclass(frozen=True)
class OutputPaths:
    run_dir: Path
    seurat_dir: Path
    table_dir: Path
    per_study_table_dir: Path
    plot_dir: Path
    umap_grid_dir: Path
    summary_plot_dir: Path
    h5ad_dir: Path
    diagnostics_dir: Path


def safe_token(value: object) -> str:
    token = re.sub(r"[^A-Za-z0-9]+", "_", str(value).strip())
    token = re.sub(r"_+", "_", token).strip("_")
    return token or "value"


def natural_sort_key(value: object) -> list[object]:
    parts = re.split(r"(\d+)", str(value))
    return [int(part) if part.isdigit() else part.lower() for part in parts]


def parse_study_ids(raw: Sequence[str] | None) -> list[str] | None:
    if not raw:
        return None
    values: list[str] = []
    for item in raw:
        values.extend(part for part in re.split(r"[,;\s]+", item) if part and part.lower() not in {"none", "null", "na"})
    return values or None


def run_dir(project_root: str | Path | None = None, run_label: str = RUN_LABEL_DEFAULT) -> Path:
    return resolve_project_root(project_root) / "results" / RESULTS_DIRNAME / run_label


def output_paths(project_root: str | Path | None = None, run_label: str = RUN_LABEL_DEFAULT) -> OutputPaths:
    root = run_dir(project_root, run_label)
    return OutputPaths(
        run_dir=root,
        seurat_dir=root / "seurat" / "per_study",
        table_dir=root / "tables",
        per_study_table_dir=root / "tables" / "per_study",
        plot_dir=root / "plots",
        umap_grid_dir=root / "plots" / "umap_grids",
        summary_plot_dir=root / "plots" / "summary",
        h5ad_dir=root / "h5ad",
        diagnostics_dir=root / "diagnostics",
    )


def ensure_output_dirs(project_root: str | Path | None = None, run_label: str = RUN_LABEL_DEFAULT) -> OutputPaths:
    paths = output_paths(project_root, run_label)
    for path in paths.__dict__.values():
        path.mkdir(parents=True, exist_ok=True)
    return paths


def count_tsv_rows(path: str | Path) -> int:
    p = Path(path)
    opener = gzip.open if p.suffix == ".gz" else open
    with opener(p, "rt") as handle:
        return max(sum(1 for _ in handle) - 1, 0)


def parse_gw_numeric(value: object) -> float:
    nums = re.findall(r"\d+", str(value))
    if not nums:
        return np.nan
    return float(np.mean([int(num) for num in nums]))


def study_table(project_root: str | Path | None = None, run_label: str = RUN_LABEL_DEFAULT) -> pd.DataFrame:
    paths = output_paths(project_root, run_label)
    existing = paths.table_dir / "cross_study_shi_seurat_label_transfer_studies.tsv"
    if existing.exists():
        return pd.read_csv(existing, sep="\t")
    return pd.DataFrame(
        {
            "study_id": STUDY_ORDER,
            "study_label": [STUDY_LABELS[x] for x in STUDY_ORDER],
            "study_order": np.arange(1, len(STUDY_ORDER) + 1),
        }
    )


def per_study_obs_path(study_id: str, project_root: str | Path | None = None, run_label: str = RUN_LABEL_DEFAULT) -> Path:
    return output_paths(project_root, run_label).per_study_table_dir / f"{study_id}_shi_seurat_label_transfer_obs.tsv.gz"


def setup_tables(project_root: str | Path | None = None, run_label: str = RUN_LABEL_DEFAULT) -> pd.DataFrame:
    paths = ensure_output_dirs(project_root, run_label)
    table = study_table(project_root, run_label)
    if "study_order" not in table.columns:
        table["study_order"] = table["study_id"].map({study_id: idx for idx, study_id in enumerate(STUDY_ORDER, start=1)})
    rows = []
    for _, row in table.iterrows():
        study_id = str(row["study_id"])
        path = per_study_obs_path(study_id, project_root, run_label)
        rows.append(
            {
                "study_id": study_id,
                "study_label": row.get("study_label", STUDY_LABELS.get(study_id, study_id)),
                "study_order": row.get("study_order", STUDY_ORDER.index(study_id) + 1 if study_id in STUDY_ORDER else np.nan),
                "per_study_obs_path": str(path),
                "per_study_obs_exists": path.exists(),
                "n_cells_per_study_obs": count_tsv_rows(path) if path.exists() else np.nan,
            }
        )
    ready = pd.DataFrame(rows).sort_values("study_order")
    ready.to_csv(paths.table_dir / "cross_study_shi_seurat_label_transfer_readiness.tsv", sep="\t", index=False)
    return ready


def label_score_columns(data: pd.DataFrame) -> list[str]:
    known = [col for col in LABEL_SCORE_MAP.values() if col in data.columns]
    other = [
        col
        for col in data.columns
        if col.startswith("shi_seurat_full_prediction_score_")
        and col not in known
        and col != "shi_seurat_full_prediction_score"
    ]
    return known + sorted(other, key=natural_sort_key)


def week_score_columns(data: pd.DataFrame) -> list[str]:
    cols = [
        col
        for col in data.columns
        if col.startswith("shi_seurat_full_week_prediction_score_")
        and col != "shi_seurat_full_week_prediction_score"
    ]
    return sorted(cols, key=lambda col: (parse_gw_numeric(col), natural_sort_key(col)))


def label_from_score_col(col: str) -> str:
    reverse = {value: key for key, value in LABEL_SCORE_MAP.items()}
    if col in reverse:
        return reverse[col]
    return col.removeprefix("shi_seurat_full_prediction_score_").replace("_", " ")


def week_label_from_score_col(col: str) -> str:
    return col.removeprefix("shi_seurat_full_week_prediction_score_")


def normalize_obs_table(data: pd.DataFrame) -> pd.DataFrame:
    data = data.copy()
    if "sample_label" not in data.columns:
        data["sample_label"] = data.get("sample", "")
    if "cluster" not in data.columns and "seurat_clusters" in data.columns:
        data["cluster"] = data["seurat_clusters"]
    if "sample" not in data.columns and "orig.ident" in data.columns:
        data["sample"] = data["orig.ident"]
    for col in BASE_COLUMNS:
        if col not in data.columns:
            data[col] = np.nan if col.startswith("shi_") or col.startswith("umap") else ""
    data["study_id"] = data["study_id"].astype(str)
    data["study_label"] = data["study_label"].astype(str)
    data["sample"] = data["sample"].astype(str)
    data["sample_label"] = data["sample_label"].astype(str)
    data["cluster"] = data["cluster"].astype(str)
    data["shi_seurat_full_uncertainty_score"] = pd.to_numeric(
        data["shi_seurat_full_uncertainty_score"], errors="coerce"
    )
    if data["shi_seurat_full_uncertainty_score"].isna().all():
        data["shi_seurat_full_uncertainty_score"] = 1.0 - pd.to_numeric(
            data["shi_seurat_full_prediction_score"], errors="coerce"
        )
    data["shi_seurat_full_week_uncertainty_score"] = pd.to_numeric(
        data["shi_seurat_full_week_uncertainty_score"], errors="coerce"
    )
    if data["shi_seurat_full_week_uncertainty_score"].isna().all():
        data["shi_seurat_full_week_uncertainty_score"] = 1.0 - pd.to_numeric(
            data["shi_seurat_full_week_prediction_score"], errors="coerce"
        )
    return data


def resolve_project_path(path: str | Path, project_root: str | Path | None = None) -> Path:
    raw = Path(path).expanduser()
    if raw.is_absolute():
        return raw
    return resolve_project_root(project_root) / raw


def augment_reused_varela_from_h5ad(
    data: pd.DataFrame,
    study_id: str,
    project_root: str | Path | None = None,
) -> pd.DataFrame:
    """Add UMAP/sample/cluster metadata to reused Varela Seurat prediction tables."""
    if study_id not in VARELA_H5AD_PATHS or {"umap_1", "umap_2"}.issubset(data.columns):
        return data
    h5ad_path = resolve_project_path(VARELA_H5AD_PATHS[study_id], project_root)
    if not h5ad_path.exists():
        raise FileNotFoundError(f"{study_id} reused predictions need UMAP coordinates; missing H5AD: {h5ad_path}")
    from .cross_study_marker_expression import (
        _read_h5ad_dataframe_index,
        _read_h5ad_dataframe_value,
        _read_h5ad_obsm,
    )
    import h5py

    with h5py.File(h5ad_path, "r") as handle:
        obs_names = pd.Index(_read_h5ad_dataframe_index(handle["obs"]).astype(str), name="_h5ad_cell_id")
        umap, _ = _read_h5ad_obsm(handle, ("X_umap_seurat", "X_umap"))
        sample_values = _read_h5ad_dataframe_value(handle["obs"], "orig.ident")
        if sample_values is None:
            sample_values = np.full(obs_names.shape[0], "", dtype=object)
        cluster_values = _read_h5ad_dataframe_value(handle["obs"], "seurat_clusters")
        if cluster_values is None:
            cluster_values = np.full(obs_names.shape[0], "", dtype=object)
    meta = pd.DataFrame(
        {
            "_h5ad_cell_id": obs_names.astype(str),
            "_h5ad_umap_1": umap[:, 0],
            "_h5ad_umap_2": umap[:, 1],
            "_h5ad_sample": sample_values.astype(str),
            "_h5ad_cluster": cluster_values.astype(str),
        }
    )
    out = data.copy()
    join_col = ""
    for candidate in ["cell_id", "cell_id_for_join"]:
        if candidate in out.columns and out[candidate].astype(str).isin(set(meta["_h5ad_cell_id"])).all():
            join_col = candidate
            break
    if not join_col:
        raise ValueError(f"{study_id} reused prediction cell IDs could not be matched to H5AD obs names: {h5ad_path}")
    out = out.merge(meta, left_on=join_col, right_on="_h5ad_cell_id", how="left", sort=False)
    if out["_h5ad_umap_1"].isna().any() or out["_h5ad_umap_2"].isna().any():
        raise ValueError(f"{study_id} H5AD UMAP join left missing coordinates.")
    out["umap_1"] = out["_h5ad_umap_1"]
    out["umap_2"] = out["_h5ad_umap_2"]
    if "sample" not in out.columns:
        out["sample"] = out["_h5ad_sample"]
    if "cluster" not in out.columns:
        out["cluster"] = out["_h5ad_cluster"]
    drop_cols = [col for col in out.columns if col.startswith("_h5ad_")]
    return out.drop(columns=drop_cols)


def compute_expected_gw(data: pd.DataFrame) -> pd.DataFrame:
    data = data.copy()
    cols = week_score_columns(data)
    if not cols:
        return data
    labels = [week_label_from_score_col(col) for col in cols]
    numeric = np.asarray([parse_gw_numeric(label) for label in labels], dtype=float)
    if np.isnan(numeric).any():
        bad = [label for label, value in zip(labels, numeric, strict=True) if np.isnan(value)]
        raise ValueError("Could not parse numeric GW labels from score columns: " + ", ".join(bad))
    even_values = np.asarray([sorted(set(numeric)).index(value) + 1 for value in numeric], dtype=float)
    scores = data[cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    denom = np.nansum(scores, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        expected_numeric = np.nansum(scores * numeric[None, :], axis=1) / denom
        expected_even = np.nansum(scores * even_values[None, :], axis=1) / denom
    data["shi_seurat_full_expected_shi_gw_numeric"] = np.where(np.isfinite(expected_numeric), expected_numeric, np.nan)
    data["shi_seurat_full_expected_shi_gw_even"] = np.where(np.isfinite(expected_even), expected_even, np.nan)
    return data


def validate_combined(data: pd.DataFrame, require_studies: Sequence[str] | None = None) -> None:
    missing_cols = [col for col in BASE_COLUMNS if col not in data.columns]
    if missing_cols:
        raise ValueError("Combined obs table missing required columns: " + ", ".join(missing_cols))
    score_cols = (
        label_score_columns(data)
        + week_score_columns(data)
        + ["shi_seurat_full_prediction_score", "shi_seurat_full_week_prediction_score"]
    )
    for col in score_cols:
        values = pd.to_numeric(data[col], errors="coerce")
        bad = values.notna() & ((values < -1e-8) | (values > 1.0 + 1e-8))
        if bad.any():
            raise ValueError(f"{col} contains scores outside [0,1].")
    if require_studies:
        missing = [study_id for study_id in require_studies if study_id not in set(data["study_id"].astype(str))]
        if missing:
            raise ValueError("Missing requested studies in combined obs table: " + ", ".join(missing))
    gw = pd.to_numeric(data["shi_seurat_full_expected_shi_gw_numeric"], errors="coerce")
    if gw.notna().any() and ((gw.dropna() < 0).any() or (gw.dropna() > 40).any()):
        raise ValueError("Expected GW numeric values look outside a plausible range.")


def load_per_study_tables(
    project_root: str | Path | None = None,
    run_label: str = RUN_LABEL_DEFAULT,
    study_ids: Sequence[str] | None = None,
) -> pd.DataFrame:
    wanted = list(study_ids or STUDY_ORDER)
    tables = []
    missing = []
    for study_id in wanted:
        path = per_study_obs_path(study_id, project_root, run_label)
        if not path.exists():
            missing.append(f"{study_id}: {path}")
            continue
        table = pd.read_csv(path, sep="\t")
        table = augment_reused_varela_from_h5ad(table, study_id, project_root)
        tables.append(normalize_obs_table(table))
    if missing:
        raise FileNotFoundError("Missing per-study Shi obs tables:\n" + "\n".join(missing))
    combined = pd.concat(tables, ignore_index=True)
    combined = compute_expected_gw(combined)
    combined["study_id"] = pd.Categorical(combined["study_id"], categories=STUDY_ORDER, ordered=True)
    combined = combined.sort_values(["study_id", "sample", "cell_id"]).reset_index(drop=True)
    combined["study_id"] = combined["study_id"].astype(str)
    validate_combined(combined, wanted)
    return combined


def write_long_score_tables(data: pd.DataFrame, paths: OutputPaths) -> tuple[pd.DataFrame, pd.DataFrame]:
    id_cols = ["cell_id", "study_id", "study_label", "sample", "sample_label", "cluster"]
    label_cols = label_score_columns(data)
    label_long = data[id_cols + label_cols].melt(
        id_vars=id_cols,
        value_vars=label_cols,
        var_name="score_column",
        value_name="prediction_score",
    )
    label_long["shi_label"] = label_long["score_column"].map(label_from_score_col)
    label_long["shi_label"] = pd.Categorical(label_long["shi_label"], categories=SHI_LABEL_ORDER, ordered=True)
    label_long = label_long.sort_values(["study_id", "cell_id", "shi_label"])
    label_long.to_csv(paths.table_dir / "cross_study_shi_seurat_label_scores_long.tsv.gz", sep="\t", index=False)

    week_cols = week_score_columns(data)
    week_long = data[id_cols + week_cols].melt(
        id_vars=id_cols,
        value_vars=week_cols,
        var_name="score_column",
        value_name="prediction_score",
    )
    week_long["shi_week_label"] = week_long["score_column"].map(week_label_from_score_col)
    week_long["shi_week_numeric"] = week_long["shi_week_label"].map(parse_gw_numeric)
    week_long = week_long.sort_values(["study_id", "cell_id", "shi_week_numeric"])
    week_long.to_csv(paths.table_dir / "cross_study_shi_seurat_week_scores_long.tsv.gz", sep="\t", index=False)
    return label_long, week_long


def write_summary_tables(data: pd.DataFrame, paths: OutputPaths) -> dict[str, Path]:
    outputs: dict[str, Path] = {}
    label = "shi_seurat_full_predicted_shi_label"
    label_cat = pd.Categorical(data[label], categories=SHI_LABEL_ORDER, ordered=True)
    tmp = data.assign(_label=label_cat)

    counts_study = pd.crosstab(tmp["study_id"], tmp["_label"]).reindex(index=STUDY_ORDER, columns=SHI_LABEL_ORDER, fill_value=0)
    counts_study.index.name = "study_id"
    counts_study = counts_study.reset_index()
    counts_study.insert(1, "study_label", counts_study["study_id"].map(STUDY_LABELS))
    outputs["label_counts_by_study"] = paths.table_dir / "cross_study_shi_label_counts_by_study.tsv"
    counts_study.to_csv(outputs["label_counts_by_study"], sep="\t", index=False)
    frac_study = counts_study.copy()
    denom = frac_study[SHI_LABEL_ORDER].sum(axis=1).replace(0, np.nan)
    frac_study[SHI_LABEL_ORDER] = frac_study[SHI_LABEL_ORDER].div(denom, axis=0).fillna(0)
    outputs["label_fractions_by_study"] = paths.table_dir / "cross_study_shi_label_fractions_by_study.tsv"
    frac_study.to_csv(outputs["label_fractions_by_study"], sep="\t", index=False)

    sample_index = ["study_id", "study_label", "sample", "sample_label"]
    counts_sample = tmp.groupby(sample_index + ["_label"], observed=False).size().unstack("_label", fill_value=0)
    counts_sample = counts_sample.reindex(columns=SHI_LABEL_ORDER, fill_value=0).reset_index()
    outputs["label_counts_by_sample"] = paths.table_dir / "cross_study_shi_label_counts_by_sample.tsv"
    counts_sample.to_csv(outputs["label_counts_by_sample"], sep="\t", index=False)
    frac_sample = counts_sample.copy()
    denom = frac_sample[SHI_LABEL_ORDER].sum(axis=1).replace(0, np.nan)
    frac_sample[SHI_LABEL_ORDER] = frac_sample[SHI_LABEL_ORDER].div(denom, axis=0).fillna(0)
    outputs["label_fractions_by_sample"] = paths.table_dir / "cross_study_shi_label_fractions_by_sample.tsv"
    frac_sample.to_csv(outputs["label_fractions_by_sample"], sep="\t", index=False)

    score_cols = label_score_columns(data)
    mean_study = data.groupby(["study_id", "study_label"], sort=False)[score_cols].mean(numeric_only=True).reset_index()
    outputs["mean_label_scores_by_study"] = paths.table_dir / "cross_study_shi_mean_label_scores_by_study.tsv"
    mean_study.to_csv(outputs["mean_label_scores_by_study"], sep="\t", index=False)
    mean_sample = data.groupby(sample_index, sort=False)[score_cols].mean(numeric_only=True).reset_index()
    outputs["mean_label_scores_by_sample"] = paths.table_dir / "cross_study_shi_mean_label_scores_by_sample.tsv"
    mean_sample.to_csv(outputs["mean_label_scores_by_sample"], sep="\t", index=False)

    def _gw_summary(group_cols: list[str], name: str) -> None:
        summary = (
            data.groupby(group_cols, sort=False)
            .agg(
                n_cells=("cell_id", "size"),
                expected_gw_numeric_mean=("shi_seurat_full_expected_shi_gw_numeric", "mean"),
                expected_gw_numeric_median=("shi_seurat_full_expected_shi_gw_numeric", "median"),
                expected_gw_numeric_q10=("shi_seurat_full_expected_shi_gw_numeric", lambda x: x.quantile(0.10)),
                expected_gw_numeric_q90=("shi_seurat_full_expected_shi_gw_numeric", lambda x: x.quantile(0.90)),
                expected_gw_even_mean=("shi_seurat_full_expected_shi_gw_even", "mean"),
                expected_gw_even_median=("shi_seurat_full_expected_shi_gw_even", "median"),
                max_gw_score_mean=("shi_seurat_full_week_prediction_score", "mean"),
            )
            .reset_index()
        )
        outputs[name] = paths.table_dir / f"cross_study_shi_gw_summary_by_{'sample' if 'sample' in name else 'study'}.tsv"
        summary.to_csv(outputs[name], sep="\t", index=False)

    _gw_summary(["study_id", "study_label"], "gw_summary_by_study")
    _gw_summary(sample_index, "gw_summary_by_sample")

    diag_src = paths.diagnostics_dir / "cross_study_shi_transfer_diagnostics_summary.tsv"
    diag_dst = paths.table_dir / "cross_study_shi_transfer_diagnostics_summary.tsv"
    if diag_src.exists():
        diag = pd.read_csv(diag_src, sep="\t")
        diag.to_csv(diag_dst, sep="\t", index=False)
        outputs["diagnostics_summary"] = diag_dst
    return outputs


def combine_tables(
    project_root: str | Path | None = None,
    run_label: str = RUN_LABEL_DEFAULT,
    study_ids: Sequence[str] | None = None,
) -> pd.DataFrame:
    paths = ensure_output_dirs(project_root, run_label)
    data = load_per_study_tables(project_root, run_label, study_ids)
    output = paths.table_dir / "cross_study_shi_seurat_label_transfer_obs.tsv.gz"
    data.to_csv(output, sep="\t", index=False)
    write_long_score_tables(data, paths)
    write_summary_tables(data, paths)
    return data


def load_combined_table(project_root: str | Path | None = None, run_label: str = RUN_LABEL_DEFAULT) -> pd.DataFrame:
    path = output_paths(project_root, run_label).table_dir / "cross_study_shi_seurat_label_transfer_obs.tsv.gz"
    if not path.exists():
        return combine_tables(project_root, run_label)
    data = pd.read_csv(path, sep="\t")
    data = compute_expected_gw(normalize_obs_table(data))
    validate_combined(data)
    return data


def ordered_study_ids(data: pd.DataFrame) -> list[str]:
    present = set(data["study_id"].astype(str))
    return [study_id for study_id in STUDY_ORDER if study_id in present]


def downsample_by_study(data: pd.DataFrame, max_cells_per_study: int | None, random_state: int = 0) -> tuple[pd.DataFrame, pd.DataFrame]:
    if max_cells_per_study is None or max_cells_per_study <= 0:
        return data.copy(), pd.DataFrame(columns=["study_id", "n_cells_input", "n_cells_plotted", "downsampled"])
    rng = np.random.default_rng(random_state)
    parts = []
    rows = []
    for study_id, group in data.groupby("study_id", sort=False):
        if group.shape[0] > max_cells_per_study:
            idx = rng.choice(group.index.to_numpy(), size=max_cells_per_study, replace=False)
            sub = group.loc[np.sort(idx)].copy()
        else:
            sub = group.copy()
        parts.append(sub)
        rows.append(
            {
                "study_id": study_id,
                "n_cells_input": group.shape[0],
                "n_cells_plotted": sub.shape[0],
                "downsampled": group.shape[0] != sub.shape[0],
            }
        )
    return pd.concat(parts, ignore_index=True), pd.DataFrame(rows)


def _prep_axes(ax: plt.Axes) -> None:
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("equal", adjustable="box")
    for spine in ax.spines.values():
        spine.set_visible(False)


def _study_subsets(data: pd.DataFrame) -> tuple[list[str], dict[str, pd.DataFrame], dict[str, int]]:
    ids = ordered_study_ids(data)
    subsets = {}
    counts = {}
    for study_id in ids:
        subset = data.loc[data["study_id"].astype(str) == study_id].copy()
        finite = np.isfinite(pd.to_numeric(subset["umap_1"], errors="coerce")) & np.isfinite(
            pd.to_numeric(subset["umap_2"], errors="coerce")
        )
        subset = subset.loc[finite].copy()
        subsets[study_id] = subset
        counts[study_id] = subset.shape[0]
    return ids, subsets, counts


def save_figure(fig: plt.Figure, output_path: str | Path, *, also_pdf: bool = True, also_svg: bool = False) -> list[Path]:
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    paths = [out]
    fig.savefig(out, dpi=300, bbox_inches="tight")
    if also_pdf:
        pdf = out.with_suffix(".pdf")
        fig.savefig(pdf, bbox_inches="tight")
        paths.append(pdf)
    if also_svg:
        svg = out.with_suffix(".svg")
        fig.savefig(svg, bbox_inches="tight")
        paths.append(svg)
    plt.close(fig)
    return paths


def plot_categorical_umap_grid(
    data: pd.DataFrame,
    output_path: str | Path,
    category_col: str,
    title: str,
    *,
    category_order: Sequence[str],
    colors: dict[str, str] | None = None,
    point_size: float = 0.12,
    background_point_size: float = 0.05,
) -> pd.DataFrame:
    study_ids, subsets, counts = _study_subsets(data)
    n_rows = len(study_ids)
    fig, axes = plt.subplots(n_rows, 1, figsize=(4.4, max(1.65 * n_rows + 1.0, 5.0)), squeeze=False)
    color_map = colors or {cat: plt.get_cmap("tab20")(idx % 20) for idx, cat in enumerate(category_order)}
    rows = []
    for row_idx, study_id in enumerate(study_ids):
        ax = axes[row_idx, 0]
        subset = subsets[study_id]
        x = pd.to_numeric(subset["umap_1"], errors="coerce").to_numpy()
        y = pd.to_numeric(subset["umap_2"], errors="coerce").to_numpy()
        ax.scatter(x, y, s=background_point_size, c="#d7d7d7", linewidths=0, rasterized=True)
        values = subset[category_col].astype(str)
        for cat in category_order:
            mask = values == str(cat)
            if mask.any():
                ax.scatter(x[mask], y[mask], s=point_size, c=[color_map.get(str(cat), "#333333")], linewidths=0, rasterized=True)
            rows.append({"study_id": study_id, "category": cat, "n_cells_plotted": int(mask.sum())})
        _prep_axes(ax)
        ax.set_ylabel(f"{STUDY_LABELS.get(study_id, study_id)}\n(n = {counts[study_id]:,} cells)", fontsize=8, rotation=0, ha="right", va="center")
    axes[0, 0].set_title("Predicted Shi label", fontsize=9)
    handles = [
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=color_map.get(cat, "#333333"), markersize=5, label=cat)
        for cat in category_order
    ]
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(0.83, 0.5), fontsize=7, frameon=False)
    fig.suptitle(title, fontsize=12)
    fig.tight_layout(rect=(0.20, 0.03, 0.82, 0.96))
    save_figure(fig, output_path, also_pdf=True, also_svg=True)
    return pd.DataFrame(rows)


def plot_continuous_umap_grid(
    data: pd.DataFrame,
    output_path: str | Path,
    features: Sequence[tuple[str, str, float, float, str]],
    title: str,
    *,
    point_size: float = 0.12,
    background_point_size: float = 0.05,
) -> pd.DataFrame:
    study_ids, subsets, counts = _study_subsets(data)
    n_rows = len(study_ids)
    n_cols = len(features)
    fig_width = max(1.35 * n_cols + 2.0, 6.0)
    fig_height = max(1.65 * n_rows + 1.2, 5.0)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_width, fig_height), squeeze=False)
    rows = []
    for row_idx, study_id in enumerate(study_ids):
        subset = subsets[study_id]
        x = pd.to_numeric(subset["umap_1"], errors="coerce").to_numpy()
        y = pd.to_numeric(subset["umap_2"], errors="coerce").to_numpy()
        for col_idx, (col, label, vmin, vmax, cmap) in enumerate(features):
            ax = axes[row_idx, col_idx]
            values = pd.to_numeric(subset[col], errors="coerce").to_numpy(dtype=float)
            finite = np.isfinite(values)
            ax.scatter(x, y, s=background_point_size, c="#d7d7d7", linewidths=0, rasterized=True)
            if finite.any():
                order = np.argsort(values[finite])
                ax.scatter(
                    x[finite][order],
                    y[finite][order],
                    s=point_size,
                    c=np.clip(values[finite][order], vmin, vmax),
                    cmap=cmap,
                    norm=Normalize(vmin=vmin, vmax=vmax),
                    linewidths=0,
                    rasterized=True,
                )
            _prep_axes(ax)
            if row_idx == 0:
                ax.set_title(label, fontsize=8)
            rows.append(
                {
                    "study_id": study_id,
                    "feature": col,
                    "n_cells_plotted": counts[study_id],
                    "n_finite_values": int(finite.sum()),
                    "vmin": vmin,
                    "vmax": vmax,
                }
            )
        axes[row_idx, 0].set_ylabel(
            f"{STUDY_LABELS.get(study_id, study_id)}\n(n = {counts[study_id]:,} cells)",
            fontsize=8,
            rotation=0,
            ha="right",
            va="center",
        )
    fig.suptitle(title, fontsize=12, y=0.995)
    fig.tight_layout(rect=(0.16, 0.12, 0.995, 0.965))
    fig.canvas.draw()
    cbar_y = max(0.04, axes[-1, 0].get_position().y0 - 0.055)
    for col_idx, (_, label, vmin, vmax, cmap) in enumerate(features):
        pos = axes[-1, col_idx].get_position()
        cax = fig.add_axes([pos.x0, cbar_y, pos.x1 - pos.x0, 0.010])
        sm = ScalarMappable(norm=Normalize(vmin=vmin, vmax=vmax), cmap=cmap)
        sm.set_array([])
        cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
        cbar.ax.tick_params(labelsize=5.5, length=1.5, pad=1)
        cbar.set_ticks([vmin, vmax])
        cbar.set_ticklabels([f"{vmin:g}", f"{vmax:g}"])
        cbar.outline.set_linewidth(0.4)
    save_figure(fig, output_path, also_pdf=True, also_svg=False)
    return pd.DataFrame(rows)


def plot_gw_umap_grid(data: pd.DataFrame, output_path: str | Path) -> pd.DataFrame:
    study_ids, subsets, counts = _study_subsets(data)
    week_labels = sorted(pd.unique(data["shi_seurat_full_predicted_shi_week_label"].dropna().astype(str)), key=lambda x: parse_gw_numeric(x))
    cmap_week = plt.get_cmap("viridis", max(len(week_labels), 2))
    week_colors = {label: cmap_week(idx / max(len(week_labels) - 1, 1)) for idx, label in enumerate(week_labels)}
    features = [
        ("shi_seurat_full_expected_shi_gw_numeric", "Expected Shi GW", "continuous_gw"),
        ("shi_seurat_full_predicted_shi_week_label", "Predicted Shi GW", "categorical_week"),
        ("shi_seurat_full_week_prediction_score", "Max GW score", "continuous_score"),
    ]
    n_rows = len(study_ids)
    n_cols = len(features)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6.3, max(1.65 * n_rows + 1.2, 5.0)), squeeze=False)
    rows = []
    for row_idx, study_id in enumerate(study_ids):
        subset = subsets[study_id]
        x = pd.to_numeric(subset["umap_1"], errors="coerce").to_numpy()
        y = pd.to_numeric(subset["umap_2"], errors="coerce").to_numpy()
        for col_idx, (col, label, kind) in enumerate(features):
            ax = axes[row_idx, col_idx]
            ax.scatter(x, y, s=0.05, c="#d7d7d7", linewidths=0, rasterized=True)
            if kind == "categorical_week":
                values = subset[col].astype(str)
                for week in week_labels:
                    mask = values == week
                    if mask.any():
                        ax.scatter(x[mask], y[mask], s=0.12, c=[week_colors[week]], linewidths=0, rasterized=True)
                finite = values.notna().to_numpy()
            else:
                values = pd.to_numeric(subset[col], errors="coerce").to_numpy(dtype=float)
                finite = np.isfinite(values)
                if finite.any():
                    if kind == "continuous_gw":
                        norm = Normalize(vmin=9, vmax=18)
                        cmap = "viridis"
                        clipped = np.clip(values[finite], 9, 18)
                    else:
                        norm = Normalize(vmin=0, vmax=1)
                        cmap = "magma"
                        clipped = np.clip(values[finite], 0, 1)
                    order = np.argsort(values[finite])
                    ax.scatter(x[finite][order], y[finite][order], s=0.12, c=clipped[order], cmap=cmap, norm=norm, linewidths=0, rasterized=True)
            _prep_axes(ax)
            if row_idx == 0:
                ax.set_title(label, fontsize=8)
            rows.append({"study_id": study_id, "feature": col, "n_cells_plotted": counts[study_id], "n_finite_values": int(np.sum(finite))})
        axes[row_idx, 0].set_ylabel(
            f"{STUDY_LABELS.get(study_id, study_id)}\n(n = {counts[study_id]:,} cells)",
            fontsize=8,
            rotation=0,
            ha="right",
            va="center",
        )
    fig.suptitle("Cross-study Shi gestational-week predictions", fontsize=12, y=0.995)
    fig.tight_layout(rect=(0.18, 0.15, 0.99, 0.965))
    fig.canvas.draw()
    y = max(0.04, axes[-1, 0].get_position().y0 - 0.055)
    for col_idx, (_, _, kind) in enumerate(features):
        pos = axes[-1, col_idx].get_position()
        cax = fig.add_axes([pos.x0, y, pos.x1 - pos.x0, 0.010])
        if kind == "continuous_gw":
            sm = ScalarMappable(norm=Normalize(vmin=9, vmax=18), cmap="viridis")
            cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
            cbar.set_ticks([9, 18])
            cbar.set_ticklabels(["GW9", "GW18"])
        elif kind == "continuous_score":
            sm = ScalarMappable(norm=Normalize(vmin=0, vmax=1), cmap="magma")
            cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
            cbar.set_ticks([0, 1])
            cbar.set_ticklabels(["0", "1"])
        else:
            cmap = ListedColormap([week_colors[label] for label in week_labels])
            norm = BoundaryNorm(np.arange(len(week_labels) + 1), cmap.N)
            sm = ScalarMappable(norm=norm, cmap=cmap)
            cbar = fig.colorbar(sm, cax=cax, orientation="horizontal", ticks=np.arange(len(week_labels)) + 0.5)
            cbar.set_ticklabels(week_labels)
        cbar.ax.tick_params(labelsize=5.5, length=1.5, pad=1)
        cbar.outline.set_linewidth(0.4)
    save_figure(fig, output_path, also_pdf=True, also_svg=False)
    return pd.DataFrame(rows)


def plot_stacked_bar(table: pd.DataFrame, output_path: str | Path, group_cols: Sequence[str], title: str, *, fraction: bool) -> None:
    labels = [label for label in SHI_LABEL_ORDER if label in table.columns]
    plot = table.copy()
    if fraction:
        values = plot[labels].to_numpy(dtype=float)
    else:
        denom = plot[labels].sum(axis=1).replace(0, np.nan)
        values = plot[labels].div(denom, axis=0).fillna(0).to_numpy(dtype=float)
    names = [" | ".join(map(str, row)) for row in plot[list(group_cols)].to_numpy()]
    fig, ax = plt.subplots(figsize=(max(7.0, 0.32 * len(names) + 2), 4.2))
    bottom = np.zeros(len(plot))
    for idx, label in enumerate(labels):
        vals = values[:, idx]
        ax.bar(np.arange(len(plot)), vals, bottom=bottom, color=LABEL_COLORS.get(label, "#999999"), width=0.85, label=label)
        bottom += vals
    ax.set_xticks(np.arange(len(plot)))
    ax.set_xticklabels(names, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("Fraction of cells")
    ax.set_title(title)
    ax.set_ylim(0, 1)
    ax.legend(frameon=False, fontsize=7, ncol=2, bbox_to_anchor=(1.01, 1.0), loc="upper left")
    fig.tight_layout()
    save_figure(fig, output_path, also_pdf=True)


def plot_heatmap(table: pd.DataFrame, output_path: str | Path, group_cols: Sequence[str], title: str) -> None:
    score_cols = [col for col in LABEL_SCORE_MAP.values() if col in table.columns]
    names = [" | ".join(map(str, row)) for row in table[list(group_cols)].to_numpy()]
    matrix = table[score_cols].to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(max(7, 0.55 * len(score_cols) + 3), max(3.5, 0.26 * len(names) + 1.5)))
    im = ax.imshow(matrix, aspect="auto", vmin=0, vmax=1, cmap="magma")
    ax.set_xticks(np.arange(len(score_cols)))
    ax.set_xticklabels([label_from_score_col(col) for col in score_cols], rotation=45, ha="right", fontsize=7)
    ax.set_yticks(np.arange(len(names)))
    ax.set_yticklabels(names, fontsize=7)
    ax.set_title(title)
    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("Mean Seurat prediction score")
    fig.tight_layout()
    save_figure(fig, output_path, also_pdf=True)


def _ridge_density(ax: plt.Axes, values: np.ndarray, y: float, color: str, xlim: tuple[float, float]) -> None:
    values = values[np.isfinite(values)]
    if values.size < 2:
        return
    bins = np.linspace(xlim[0], xlim[1], 80)
    hist, edges = np.histogram(values, bins=bins, density=True)
    centers = (edges[:-1] + edges[1:]) / 2
    if hist.max() > 0:
        hist = hist / hist.max() * 0.75
    ax.fill_between(centers, y, y + hist, color=color, alpha=0.75, linewidth=0)
    ax.plot(centers, y + hist, color=color, linewidth=0.8)


def plot_ridge(data: pd.DataFrame, output_path: str | Path, value_col: str, group_cols: Sequence[str], title: str, xlim: tuple[float, float]) -> None:
    groups = list(data.groupby(list(group_cols), sort=False))
    fig, ax = plt.subplots(figsize=(7.0, max(3.8, 0.34 * len(groups) + 1.4)))
    colors = plt.get_cmap("tab20")(np.linspace(0, 1, max(len(groups), 2)))
    labels = []
    for idx, (key, group) in enumerate(groups):
        label = " | ".join(map(str, key if isinstance(key, tuple) else (key,)))
        labels.append(label)
        values = pd.to_numeric(group[value_col], errors="coerce").to_numpy(dtype=float)
        _ridge_density(ax, values, float(idx), colors[idx % len(colors)], xlim)
    ax.set_yticks(np.arange(len(labels)) + 0.35)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlim(*xlim)
    ax.set_xlabel(value_col)
    ax.set_title(title)
    ax.spines[["top", "right", "left"]].set_visible(False)
    fig.tight_layout()
    save_figure(fig, output_path, also_pdf=True)


def plot_violin(data: pd.DataFrame, output_path: str | Path, value_col: str, title: str) -> None:
    groups = [(study_id, data.loc[data["study_id"].astype(str) == study_id, value_col].dropna().astype(float).to_numpy()) for study_id in ordered_study_ids(data)]
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    ax.violinplot([vals for _, vals in groups], positions=np.arange(len(groups)), showmeans=True, showextrema=False)
    ax.set_xticks(np.arange(len(groups)))
    ax.set_xticklabels([STUDY_LABELS.get(study_id, study_id) for study_id, _ in groups], rotation=35, ha="right", fontsize=8)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Seurat prediction score")
    ax.set_title(title)
    fig.tight_layout()
    save_figure(fig, output_path, also_pdf=True)


def make_umap_grids(
    data: pd.DataFrame,
    paths: OutputPaths,
    *,
    max_cells_per_study: int | None = None,
) -> pd.DataFrame:
    plot_data, downsample = downsample_by_study(data, max_cells_per_study)
    if not downsample.empty:
        downsample.to_csv(paths.diagnostics_dir / "cross_study_shi_plot_downsampling.tsv", sep="\t", index=False)
    manifests = []
    manifest = plot_categorical_umap_grid(
        plot_data,
        paths.umap_grid_dir / "cross_study_umap_shi_seurat_full_predicted_shi_label_grid.png",
        "shi_seurat_full_predicted_shi_label",
        "Cross-study predicted Shi labels",
        category_order=SHI_LABEL_ORDER,
        colors=LABEL_COLORS,
    )
    manifest.insert(0, "plot_name", "predicted_shi_label_grid")
    manifests.append(manifest)

    main_features = [
        (LABEL_SCORE_MAP["MGE"], "MGE score", 0, 1, "magma"),
        (LABEL_SCORE_MAP["LGE"], "LGE score", 0, 1, "magma"),
        (LABEL_SCORE_MAP["CGE"], "CGE score", 0, 1, "magma"),
        (LABEL_SCORE_MAP["progenitor"], "progenitor score", 0, 1, "magma"),
        ("shi_seurat_full_prediction_score", "Max label score", 0, 1, "magma"),
    ]
    main_features = [feat for feat in main_features if feat[0] in plot_data.columns]
    manifest = plot_continuous_umap_grid(
        plot_data,
        paths.umap_grid_dir / "cross_study_umap_shi_seurat_full_label_score_grid.png",
        main_features,
        "Cross-study Shi label prediction scores",
    )
    manifest.insert(0, "plot_name", "label_score_grid")
    manifests.append(manifest)

    all_features = [(col, label_from_score_col(col), 0, 1, "magma") for col in label_score_columns(plot_data)]
    manifest = plot_continuous_umap_grid(
        plot_data,
        paths.umap_grid_dir / "cross_study_umap_shi_seurat_full_all_label_scores_grid.png",
        all_features,
        "Cross-study all Shi label prediction scores",
    )
    manifest.insert(0, "plot_name", "all_label_scores_grid")
    manifests.append(manifest)

    manifest = plot_gw_umap_grid(
        plot_data,
        paths.umap_grid_dir / "cross_study_umap_shi_seurat_full_gw_prediction_grid.png",
    )
    manifest.insert(0, "plot_name", "gw_prediction_grid")
    manifests.append(manifest)

    gw_features = [(col, week_label_from_score_col(col), 0, 1, "magma") for col in week_score_columns(plot_data)]
    manifest = plot_continuous_umap_grid(
        plot_data,
        paths.umap_grid_dir / "cross_study_umap_shi_seurat_full_individual_gw_scores_grid.png",
        gw_features,
        "Cross-study individual Shi GW prediction scores",
    )
    manifest.insert(0, "plot_name", "individual_gw_scores_grid")
    manifests.append(manifest)

    combined = pd.concat(manifests, ignore_index=True)
    combined.to_csv(paths.table_dir / "cross_study_shi_umap_plot_manifest.tsv", sep="\t", index=False)
    return combined


def make_summary_plots(data: pd.DataFrame, paths: OutputPaths) -> list[Path]:
    write_summary_tables(data, paths)
    outputs: list[Path] = []
    counts_study = pd.read_csv(paths.table_dir / "cross_study_shi_label_counts_by_study.tsv", sep="\t")
    counts_sample = pd.read_csv(paths.table_dir / "cross_study_shi_label_counts_by_sample.tsv", sep="\t")
    mean_study = pd.read_csv(paths.table_dir / "cross_study_shi_mean_label_scores_by_study.tsv", sep="\t")
    mean_sample = pd.read_csv(paths.table_dir / "cross_study_shi_mean_label_scores_by_sample.tsv", sep="\t")
    specs = [
        (counts_study, "cross_study_shi_label_fractions_by_study_stacked_bar.png", ["study_label"], "Predicted Shi label fractions by study"),
        (counts_sample, "cross_study_shi_label_fractions_by_sample_stacked_bar.png", ["study_label", "sample_label"], "Predicted Shi label fractions by sample"),
    ]
    for table, filename, group_cols, title in specs:
        out = paths.summary_plot_dir / filename
        plot_stacked_bar(table, out, group_cols, title, fraction=False)
        outputs.append(out)
    heat_specs = [
        (mean_study, "cross_study_shi_mean_label_scores_by_study_heatmap.png", ["study_label"], "Mean Shi label score by study"),
        (mean_sample, "cross_study_shi_mean_label_scores_by_sample_heatmap.png", ["study_label", "sample_label"], "Mean Shi label score by sample"),
    ]
    for table, filename, group_cols, title in heat_specs:
        out = paths.summary_plot_dir / filename
        plot_heatmap(table, out, group_cols, title)
        outputs.append(out)
    ridge_specs = [
        ("shi_seurat_full_expected_shi_gw_numeric", ["study_label"], "cross_study_shi_expected_gw_numeric_by_study_ridge.png", "Expected Shi GW numeric by study", (9, 18)),
        ("shi_seurat_full_expected_shi_gw_numeric", ["study_label", "sample_label"], "cross_study_shi_expected_gw_numeric_by_sample_ridge.png", "Expected Shi GW numeric by sample", (9, 18)),
        ("shi_seurat_full_expected_shi_gw_even", ["study_label"], "cross_study_shi_expected_gw_even_by_study_ridge.png", "Expected Shi GW even-spacing by study", (1, max(2, data["shi_seurat_full_expected_shi_gw_even"].max()))),
        ("shi_seurat_full_expected_shi_gw_even", ["study_label", "sample_label"], "cross_study_shi_expected_gw_even_by_sample_ridge.png", "Expected Shi GW even-spacing by sample", (1, max(2, data["shi_seurat_full_expected_shi_gw_even"].max()))),
    ]
    for value_col, group_cols, filename, title, xlim in ridge_specs:
        out = paths.summary_plot_dir / filename
        plot_ridge(data, out, value_col, group_cols, title, xlim)
        outputs.append(out)
    violin_specs = [
        ("shi_seurat_full_prediction_score", "cross_study_shi_max_label_prediction_score_by_study_violin.png", "Max label prediction score by study"),
        ("shi_seurat_full_week_prediction_score", "cross_study_shi_max_gw_prediction_score_by_study_violin.png", "Max GW prediction score by study"),
    ]
    for value_col, filename, title in violin_specs:
        out = paths.summary_plot_dir / filename
        plot_violin(data, out, value_col, title)
        outputs.append(out)
    return outputs


def plot_all(
    project_root: str | Path | None = None,
    run_label: str = RUN_LABEL_DEFAULT,
    max_cells_per_study: int | None = None,
) -> dict[str, object]:
    paths = ensure_output_dirs(project_root, run_label)
    data = load_combined_table(project_root, run_label)
    umap_manifest = make_umap_grids(data, paths, max_cells_per_study=max_cells_per_study)
    summary_outputs = make_summary_plots(data, paths)
    complete = pd.DataFrame(
        {
            "output_root": [str(paths.run_dir)],
            "n_cells": [data.shape[0]],
            "n_studies": [data["study_id"].nunique()],
            "n_umap_grid_png": [len(list(paths.umap_grid_dir.glob("*.png")))],
            "n_summary_png": [len(list(paths.summary_plot_dir.glob("*.png")))],
        }
    )
    complete.to_csv(paths.table_dir / "cross_study_shi_seurat_label_transfer_complete.tsv", sep="\t", index=False)
    return {"data": data, "umap_manifest": umap_manifest, "summary_outputs": summary_outputs, "paths": paths}


def print_final_report(paths: OutputPaths) -> None:
    plot_list = sorted([*paths.umap_grid_dir.glob("*.png"), *paths.summary_plot_dir.glob("*.png")])
    table_list = sorted(paths.table_dir.glob("cross_study_shi*.tsv*"))
    print(f"output_root\t{paths.run_dir}", flush=True)
    print("completed_plot_list", flush=True)
    for path in plot_list:
        print(f"  {path}", flush=True)
    print("completed_table_list", flush=True)
    for path in table_list:
        print(f"  {path}", flush=True)
    diag = paths.table_dir / "cross_study_shi_transfer_diagnostics_summary.tsv"
    if diag.exists():
        print("diagnostics_summary", flush=True)
        print(pd.read_csv(diag, sep="\t").to_string(index=False), flush=True)
    print("rsync_commands", flush=True)
    print("mkdir -p /Users/ecrespo/Downloads/cross_study_shi_seurat_label_transfer_v1", flush=True)
    print(
        "rsync -avh --progress "
        "elcrespo@greatlakes.arc-ts.umich.edu:"
        f"{paths.run_dir}/plots/ "
        "/Users/ecrespo/Downloads/cross_study_shi_seurat_label_transfer_v1/plots/",
        flush=True,
    )
    print(
        "rsync -avh --progress "
        "elcrespo@greatlakes.arc-ts.umich.edu:"
        f"{paths.run_dir}/tables/ "
        "/Users/ecrespo/Downloads/cross_study_shi_seurat_label_transfer_v1/tables/",
        flush=True,
    )


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", default=None)
    parser.add_argument("--run-label", default=RUN_LABEL_DEFAULT)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("setup", help="Create output directories and readiness table.")
    combine = subparsers.add_parser("combine", help="Combine per-study obs tables and write long/summary tables.")
    combine.add_argument("--study-id", action="append", default=[])
    plot = subparsers.add_parser("plot", help="Render UMAP grids and summary plots.")
    plot.add_argument("--max-cells-per-study", type=int, default=None)
    subparsers.add_parser("report", help="Print output root, plot/table lists, diagnostics, and rsync commands.")
    all_cmd = subparsers.add_parser("all", help="Run setup, combine, plot, and report.")
    all_cmd.add_argument("--max-cells-per-study", type=int, default=None)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    if args.command == "setup":
        ready = setup_tables(args.project_root, args.run_label)
        print(ready.to_string(index=False), flush=True)
    elif args.command == "combine":
        data = combine_tables(args.project_root, args.run_label, parse_study_ids(args.study_id))
        print(f"combined_obs_rows\t{data.shape[0]}", flush=True)
    elif args.command == "plot":
        result = plot_all(args.project_root, args.run_label, max_cells_per_study=args.max_cells_per_study)
        print(result["umap_manifest"][["plot_name"]].drop_duplicates().to_string(index=False), flush=True)
    elif args.command == "report":
        print_final_report(output_paths(args.project_root, args.run_label))
    elif args.command == "all":
        setup_tables(args.project_root, args.run_label)
        plot_all(args.project_root, args.run_label, max_cells_per_study=args.max_cells_per_study)
        print_final_report(output_paths(args.project_root, args.run_label))
    else:
        raise ValueError(f"Unhandled command: {args.command}")


if __name__ == "__main__":
    main()
