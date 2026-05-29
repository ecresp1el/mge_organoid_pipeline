"""Data-source helpers for Notebook 00 raw/filtered and CellBender workflows."""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

import anndata as ad
import pandas as pd
import scanpy as sc


SAMPLE_MAP_NAME = "metadata/div30_div90_sample_id_to_biolabel_map.tsv"
DEFAULT_DATA_ROOT = Path("/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder")

CELLRANGER_DATA_SOURCES = {
    "cellranger_raw": "raw",
    "cellranger_filtered": "filtered",
}
CELLBENDER_DATA_SOURCE = "cellbender_denoised"
SUPPORTED_DATA_SOURCES = tuple([*CELLRANGER_DATA_SOURCES.keys(), CELLBENDER_DATA_SOURCE])


def find_repo_root(start: Path | str | None = None) -> Path:
    """Return the repo root containing the DIV30/DIV90 sample map."""
    start_path = Path(start or Path.cwd()).expanduser().resolve()
    for candidate in [start_path, *start_path.parents]:
        if (candidate / SAMPLE_MAP_NAME).exists():
            return candidate
    raise FileNotFoundError(f"Could not find repo root containing {SAMPLE_MAP_NAME}")


def resolve_data_root(data_root: Path | str | None = None) -> Path:
    """Return the large runtime data root, separate from the git checkout."""
    raw = data_root or os.environ.get("MGE_DATA_ROOT") or os.environ.get("PROJECT_ROOT") or DEFAULT_DATA_ROOT
    return Path(raw).expanduser().resolve()


def normalize_data_source(data_source: str) -> str:
    """Normalize and validate a supported Notebook 00 data source name."""
    normalized = str(data_source).strip().lower()
    if normalized not in SUPPORTED_DATA_SOURCES:
        valid = sorted(SUPPORTED_DATA_SOURCES)
        raise ValueError(f"data_source must be one of {valid}; got {data_source!r}")
    return normalized


def normalize_matrix_source(matrix_source: str) -> str:
    """Normalize and validate a Cell Ranger matrix source."""
    normalized = str(matrix_source).strip().lower()
    valid_sources = {"raw", "filtered", "auto"}
    if normalized not in valid_sources:
        raise ValueError(f"matrix_source must be one of {sorted(valid_sources)}; got {matrix_source!r}")
    return normalized


def resolve_matrix_dir(per_sample_metrics_csv: str | Path, matrix_source: str) -> Path:
    """Resolve a per-sample Cell Ranger raw or filtered 10x matrix directory."""
    count_dir = Path(per_sample_metrics_csv).parent / "count"
    resolved_source = normalize_matrix_source(matrix_source)

    raw_candidates = [
        count_dir / "sample_raw_feature_bc_matrix",
        count_dir / "raw_feature_bc_matrix",
    ]
    filtered_candidates = [
        count_dir / "sample_filtered_feature_bc_matrix",
        count_dir / "filtered_feature_bc_matrix",
    ]

    if resolved_source == "raw":
        candidates = raw_candidates + filtered_candidates
    elif resolved_source == "filtered":
        candidates = filtered_candidates + raw_candidates
    else:
        candidates = raw_candidates + filtered_candidates

    for candidate in candidates:
        if candidate.exists():
            return candidate

    return candidates[0]


def _cell_line_from_biological_label(label: object) -> Optional[str]:
    match = re.match(r"^(H9|79B|2E)", str(label))
    return match.group(1) if match else None


def safe_sample_filename_stem(run_sample_id: object) -> str:
    """Return the raw `.h5ad` basename stem used by Notebook 00 raw export."""
    return str(run_sample_id).replace(" ", "_")


def expected_raw_h5ad_path(data_root: Path | str, run_sample_id: object) -> Path:
    """Return the primary raw AnnData path used as CellBender input."""
    return Path(data_root).expanduser().resolve() / "raw_adata" / f"{safe_sample_filename_stem(run_sample_id)}.h5ad"


def expected_cellbender_output_h5(data_root: Path | str, run_sample_id: object) -> Path:
    """Return the primary CellBender denoised H5 path from `scripts/cellbender.sh`."""
    stem = safe_sample_filename_stem(run_sample_id)
    return Path(data_root).expanduser().resolve() / "clean_adata" / f"{stem}_cellbender_denoised.h5"


@dataclass(frozen=True)
class Notebook00SourceConfig:
    """Configuration for loading Notebook 00 Cell Ranger raw/filtered sources."""

    repo_root: Path
    data_root: Path
    data_source: str = "cellranger_filtered"
    target_divs: Sequence[str] = ("DIV30",)
    target_run_sample_ids: Optional[Sequence[str]] = None
    sample_map_name: str = SAMPLE_MAP_NAME
    strict_missing_matrix_dirs: bool = True

    @classmethod
    def from_defaults(
        cls,
        data_source: str = "cellranger_filtered",
        repo_root: Path | str | None = None,
        data_root: Path | str | None = None,
        target_divs: Sequence[str] = ("DIV30",),
        target_run_sample_ids: Optional[Sequence[str]] = None,
        strict_missing_matrix_dirs: bool = True,
    ) -> "Notebook00SourceConfig":
        """Build config from current repo plus HPC data-root conventions."""
        return cls(
            repo_root=find_repo_root(repo_root),
            data_root=resolve_data_root(data_root),
            data_source=normalize_data_source(data_source),
            target_divs=tuple(target_divs),
            target_run_sample_ids=tuple(target_run_sample_ids) if target_run_sample_ids is not None else None,
            strict_missing_matrix_dirs=strict_missing_matrix_dirs,
        )

    @property
    def sample_map_tsv(self) -> Path:
        return self.repo_root / self.sample_map_name

    @property
    def matrix_source(self) -> str:
        data_source = normalize_data_source(self.data_source)
        if data_source not in CELLRANGER_DATA_SOURCES:
            raise ValueError(f"{data_source!r} does not have a Cell Ranger matrix_source.")
        return CELLRANGER_DATA_SOURCES[data_source]


def sample_table(config: Notebook00SourceConfig) -> pd.DataFrame:
    """Return selected samples with resolved raw/filtered matrix directories."""
    if normalize_data_source(config.data_source) not in CELLRANGER_DATA_SOURCES:
        raise ValueError("sample_table() is only for Cell Ranger raw/filtered sources. Use cellbender_output_table().")

    sample_metadata_df = pd.read_csv(config.sample_map_tsv, sep="\t")
    required = {"DIV", "run_sample_id", "biological_label", "per_sample_metrics_csv"}
    missing = required.difference(sample_metadata_df.columns)
    if missing:
        raise ValueError(f"Sample map is missing required columns: {sorted(missing)}")

    sample_metadata_df = sample_metadata_df[sample_metadata_df["DIV"].isin(config.target_divs)].copy()

    if config.target_run_sample_ids is not None:
        sample_metadata_df = sample_metadata_df[
            sample_metadata_df["run_sample_id"].isin(config.target_run_sample_ids)
        ].copy()
        found = set(sample_metadata_df["run_sample_id"].astype(str))
        missing_ids = [sid for sid in config.target_run_sample_ids if sid not in found]
        if missing_ids:
            raise ValueError(f"Missing requested run_sample_id values in sample map: {missing_ids}")

    if sample_metadata_df.empty:
        raise ValueError("No samples matched target_divs/target_run_sample_ids.")

    sample_metadata_df["cell_line"] = sample_metadata_df["biological_label"].map(_cell_line_from_biological_label)
    invalid_rows = sample_metadata_df[sample_metadata_df["cell_line"].isna()][["run_sample_id", "biological_label"]]
    if not invalid_rows.empty:
        examples = invalid_rows.head(10).to_dict(orient="records")
        raise ValueError(f"Could not map biological_label values to cell_line prefixes. Examples: {examples}")

    sample_metadata_df["DIV"] = pd.Categorical(sample_metadata_df["DIV"], categories=config.target_divs, ordered=True)
    if config.target_run_sample_ids is not None:
        sample_metadata_df["run_sample_id"] = pd.Categorical(
            sample_metadata_df["run_sample_id"],
            categories=config.target_run_sample_ids,
            ordered=True,
        )

    sample_metadata_df = sample_metadata_df.sort_values(["DIV", "run_sample_id"]).reset_index(drop=True)
    sample_metadata_df["run_sample_id"] = sample_metadata_df["run_sample_id"].astype(str)
    sample_metadata_df["data_source"] = normalize_data_source(config.data_source)
    sample_metadata_df["matrix_source"] = config.matrix_source
    sample_metadata_df["matrix_dir"] = sample_metadata_df["per_sample_metrics_csv"].map(
        lambda p: str(resolve_matrix_dir(p, matrix_source=config.matrix_source))
    )
    sample_metadata_df["matrix_exists"] = sample_metadata_df["matrix_dir"].map(lambda p: Path(p).exists())

    return sample_metadata_df[
        [
            "DIV",
            "run_sample_id",
            "biological_label",
            "cell_line",
            "data_source",
            "matrix_source",
            "matrix_dir",
            "matrix_exists",
        ]
    ]


def selected_sample_metadata(config: Notebook00SourceConfig) -> pd.DataFrame:
    """Return selected sample metadata shared by raw/filtered and CellBender tables."""
    sample_metadata_df = pd.read_csv(config.sample_map_tsv, sep="\t")
    required = {"DIV", "run_sample_id", "biological_label", "per_sample_metrics_csv"}
    missing = required.difference(sample_metadata_df.columns)
    if missing:
        raise ValueError(f"Sample map is missing required columns: {sorted(missing)}")

    sample_metadata_df = sample_metadata_df[sample_metadata_df["DIV"].isin(config.target_divs)].copy()

    if config.target_run_sample_ids is not None:
        sample_metadata_df = sample_metadata_df[
            sample_metadata_df["run_sample_id"].isin(config.target_run_sample_ids)
        ].copy()
        found = set(sample_metadata_df["run_sample_id"].astype(str))
        missing_ids = [sid for sid in config.target_run_sample_ids if sid not in found]
        if missing_ids:
            raise ValueError(f"Missing requested run_sample_id values in sample map: {missing_ids}")

    if sample_metadata_df.empty:
        raise ValueError("No samples matched target_divs/target_run_sample_ids.")

    sample_metadata_df["cell_line"] = sample_metadata_df["biological_label"].map(_cell_line_from_biological_label)
    invalid_rows = sample_metadata_df[sample_metadata_df["cell_line"].isna()][["run_sample_id", "biological_label"]]
    if not invalid_rows.empty:
        examples = invalid_rows.head(10).to_dict(orient="records")
        raise ValueError(f"Could not map biological_label values to cell_line prefixes. Examples: {examples}")

    sample_metadata_df["DIV"] = pd.Categorical(sample_metadata_df["DIV"], categories=config.target_divs, ordered=True)
    if config.target_run_sample_ids is not None:
        sample_metadata_df["run_sample_id"] = pd.Categorical(
            sample_metadata_df["run_sample_id"],
            categories=config.target_run_sample_ids,
            ordered=True,
        )

    sample_metadata_df = sample_metadata_df.sort_values(["DIV", "run_sample_id"]).reset_index(drop=True)
    sample_metadata_df["run_sample_id"] = sample_metadata_df["run_sample_id"].astype(str)
    return sample_metadata_df[["DIV", "run_sample_id", "biological_label", "cell_line", "per_sample_metrics_csv"]]


def cellbender_output_table(config: Notebook00SourceConfig) -> pd.DataFrame:
    """Return expected CellBender input/output paths and existence status per sample."""
    sample_metadata_df = selected_sample_metadata(config)
    sample_metadata_df["data_source"] = CELLBENDER_DATA_SOURCE
    sample_metadata_df["raw_h5ad_path"] = sample_metadata_df["run_sample_id"].map(
        lambda run_sample_id: str(expected_raw_h5ad_path(config.data_root, run_sample_id))
    )
    sample_metadata_df["cellbender_output_h5"] = sample_metadata_df["run_sample_id"].map(
        lambda run_sample_id: str(expected_cellbender_output_h5(config.data_root, run_sample_id))
    )
    sample_metadata_df["raw_h5ad_exists"] = sample_metadata_df["raw_h5ad_path"].map(lambda path: Path(path).exists())
    sample_metadata_df["cellbender_output_exists"] = sample_metadata_df["cellbender_output_h5"].map(
        lambda path: Path(path).exists()
    )
    sample_metadata_df["cellbender_output_size_bytes"] = sample_metadata_df["cellbender_output_h5"].map(
        lambda path: Path(path).stat().st_size if Path(path).exists() else 0
    )
    return sample_metadata_df[
        [
            "DIV",
            "run_sample_id",
            "biological_label",
            "cell_line",
            "data_source",
            "raw_h5ad_path",
            "raw_h5ad_exists",
            "cellbender_output_h5",
            "cellbender_output_exists",
            "cellbender_output_size_bytes",
        ]
    ]


def read_cellranger_sample(row: pd.Series) -> ad.AnnData:
    """Read one resolved Cell Ranger sample row into AnnData with stable obs metadata."""
    matrix_dir = Path(row["matrix_dir"])
    run_sample_id = str(row["run_sample_id"])
    if not matrix_dir.exists():
        raise FileNotFoundError(f"Missing 10x matrix directory for {run_sample_id}: {matrix_dir}")

    one_sample_adata = sc.read_10x_mtx(matrix_dir, var_names="gene_symbols", make_unique=True)
    one_sample_adata.obs_names = [f"{run_sample_id}:{barcode}" for barcode in one_sample_adata.obs_names]
    one_sample_adata.obs["DIV"] = str(row["DIV"])
    one_sample_adata.obs["run_sample_id"] = run_sample_id
    one_sample_adata.obs["biological_label"] = str(row["biological_label"])
    one_sample_adata.obs["cell_line"] = str(row["cell_line"])
    one_sample_adata.obs["data_source"] = str(row["data_source"])
    one_sample_adata.obs["matrix_source"] = str(row["matrix_source"])
    one_sample_adata.obs["matrix_dir"] = str(matrix_dir)
    one_sample_adata.uns["data_source"] = str(row["data_source"])
    one_sample_adata.uns["matrix_source"] = str(row["matrix_source"])
    one_sample_adata.uns["matrix_dir"] = str(matrix_dir)
    return one_sample_adata


def load_dataset(config: Notebook00SourceConfig) -> tuple[list[str], list[ad.AnnData], pd.DataFrame]:
    """Load all selected samples for a Notebook 00 raw/filtered data source."""
    selected_samples_df = sample_table(config)
    missing = selected_samples_df.loc[~selected_samples_df["matrix_exists"], "matrix_dir"].tolist()
    if missing and config.strict_missing_matrix_dirs:
        message = "Missing per-sample 10x matrix directories:\n" + "\n".join(f" - {path}" for path in missing)
        raise FileNotFoundError(message)

    adata_names = []
    adata_list = []
    for _, row in selected_samples_df.iterrows():
        if not bool(row["matrix_exists"]):
            continue
        adata_names.append(str(row["run_sample_id"]))
        adata_list.append(read_cellranger_sample(row))

    if not adata_list:
        raise FileNotFoundError("No selected matrix directories were found.")

    return adata_names, adata_list, selected_samples_df
