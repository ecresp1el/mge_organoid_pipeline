#!/usr/bin/env python3
"""Export a DIV90 Matrix Market bundle for the Jia-lineage URD smoke run.

This exporter is the DIV90 counterpart to the DIV30 first-URD input exporter.
It keeps the major processing decisions outside the URD R step:

  - read the existing DIV90 AnnData-derived Matrix Market count export;
  - retain the Jia-lineage neuronal/root clusters for the first smoke test;
  - compute Jia RGC1/RGC2/IPC scores from the same selected Jia marker table
    used in the DIV30 workflow;
  - compute the DIV90 Jia RootScore within the retained smoke cells;
  - select top RootScore cells inside cluster 12 as the URD root;
  - write explicit tip groups approximating Jia lineage endpoints.

The output file names intentionally use the legacy `div30_first_urd_*` names
because scripts/14_div30_first_urd.R consumes that input bundle format. The
manifest records that this is a DIV90 Jia-lineage run.
"""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread, mmwrite
from scipy import sparse


DEFAULT_PROJECT_ROOT = Path("/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder")
DEFAULT_RUN_LABEL = "div90_urd_jia_lineage_smoke5k_knn100_v3_glia_cells"
DEFAULT_MAX_CELLS = 5_000
DEFAULT_ROOT_TOP_PERCENT = 2.0
DEFAULT_ROOT_MIN_CELLS = 8
DEFAULT_RETAIN_CLUSTERS = "0,1,2,3,4,5,8,9,10,11,12"


TIP_BY_CLUSTER = {
    0: "tip_lhx8_isl1",
    5: "tip_lhx8_isl1",
    8: "tip_lhx8_isl1",
    1: "tip_lhx6_nfia",
    2: "tip_crabp1_angpt2",
}

TIP_DISPLAY = {
    "tip_lhx8_isl1": "LHX8/ISL1-like lineage",
    "tip_lhx6_nfia": "LHX6/NFIA-like lineage",
    "tip_crabp1_angpt2": "CRABP1/ANGPT2-like lineage",
}

RETAINED_UNASSIGNED_CANDIDATE_CLUSTERS = {3, 11}
RETAINED_GLIA_NON_TIP_CLUSTERS = {4, 9, 10}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path(os.environ.get("PROJECT_ROOT", DEFAULT_PROJECT_ROOT)))
    parser.add_argument("--run-label", default=os.environ.get("RUN_LABEL", DEFAULT_RUN_LABEL))
    parser.add_argument("--anndata-dir", type=Path, default=None)
    parser.add_argument("--jia-marker-table", type=Path, default=None)
    parser.add_argument("--outdir", type=Path, default=None)
    parser.add_argument("--max-cells", type=int, default=int(os.environ.get("MAX_CELLS", DEFAULT_MAX_CELLS)))
    parser.add_argument("--seed", type=int, default=int(os.environ.get("SEED", "7")))
    parser.add_argument("--retain-clusters", default=os.environ.get("DIV90_RETAIN_CLUSTERS", DEFAULT_RETAIN_CLUSTERS))
    parser.add_argument("--root-cluster", type=int, default=int(os.environ.get("DIV90_ROOT_CLUSTER", "12")))
    parser.add_argument("--root-top-percent", type=float, default=float(os.environ.get("DIV90_ROOT_TOP_PERCENT", DEFAULT_ROOT_TOP_PERCENT)))
    parser.add_argument("--root-min-cells", type=int, default=int(os.environ.get("DIV90_ROOT_MIN_CELLS", DEFAULT_ROOT_MIN_CELLS)))
    return parser.parse_args()


def default_paths(project_root: Path, run_label: str) -> dict[str, Path]:
    return {
        "anndata_dir": project_root / "results" / "python_anndata" / "varela_div90_srznf5gg",
        "jia_marker_table": project_root
        / "results"
        / "jia_program_div30_scoring"
        / "jia_program_div30_scoring_v1"
        / "tables"
        / "jia_program_markers_selected.tsv",
        "outdir": project_root / "results" / "div90_jia_lineage_urd" / run_label / "inputs",
    }


def parse_cluster_list(text: str) -> list[int]:
    out = []
    for part in text.split(","):
        part = part.strip()
        if part:
            out.append(int(part))
    return sorted(set(out))


def read_single_column_tsv(path: Path, column_name: str) -> pd.Series:
    df = pd.read_csv(path, sep="\t")
    if column_name not in df.columns:
        raise ValueError(f"{path} does not contain column {column_name!r}")
    return df[column_name].astype(str)


def load_program_markers(path: Path, available_genes: set[str]) -> dict[str, list[str]]:
    markers = pd.read_csv(path, sep="\t")
    needed = {"gene", "cluster"}
    missing = needed.difference(markers.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
    programs: dict[str, list[str]] = {}
    upper_available = {g.upper(): g for g in available_genes}
    for program in ["RGC1", "RGC2", "IPC"]:
        genes = []
        for gene in markers.loc[markers["cluster"].astype(str) == program, "gene"].astype(str):
            hit = upper_available.get(gene.upper())
            if hit is not None and hit not in genes:
                genes.append(hit)
        if not genes:
            raise ValueError(f"No {program} genes from {path} were present in DIV90 features")
        programs[program] = genes
    return programs


def allocate_stratified(metadata: pd.DataFrame, max_cells: int, seed: int, root_cluster: int) -> pd.DataFrame:
    """Keep all root-cluster cells, then stratify the remaining smoke cells."""

    if max_cells <= 0 or len(metadata) <= max_cells:
        out = metadata.copy()
        out["urd_sampling_status"] = "all_retained_cells"
        return out

    root_pool = metadata[metadata["cluster_id_numeric"] == root_cluster].copy()
    nonroot = metadata[metadata["cluster_id_numeric"] != root_cluster].copy()
    remaining = max(max_cells - len(root_pool), 0)
    if remaining <= 0:
        out = root_pool.sort_values("matrix_col_index").copy()
        out["urd_sampling_status"] = f"all_root_pool_only_max_{max_cells}"
        return out

    rng = np.random.default_rng(seed)
    group_sizes = nonroot["cluster_id_numeric"].value_counts().sort_index()
    allocation = (group_sizes / group_sizes.sum() * remaining).round().astype(int)
    allocation = allocation.clip(lower=1, upper=group_sizes)
    while allocation.sum() > remaining:
        candidates = allocation[allocation > 1].sort_values(ascending=False)
        allocation.loc[candidates.index[0]] -= 1
    while allocation.sum() < remaining:
        candidates = (group_sizes - allocation).sort_values(ascending=False)
        candidates = candidates[candidates > 0]
        if candidates.empty:
            break
        allocation.loc[candidates.index[0]] += 1

    sampled = [root_pool]
    for cluster_id, n_cells in allocation.items():
        group = nonroot[nonroot["cluster_id_numeric"] == cluster_id]
        chosen = rng.choice(group.index.to_numpy(), size=int(n_cells), replace=False)
        sampled.append(group.loc[np.sort(chosen)])
    out = pd.concat(sampled, axis=0).sort_values("matrix_col_index").copy()
    out["urd_sampling_status"] = f"root_pool_plus_stratified_max_{max_cells}"
    return out


def zscore(values: np.ndarray) -> np.ndarray:
    values = values.astype(float)
    mu = np.nanmean(values)
    sigma = np.nanstd(values)
    if not np.isfinite(sigma) or sigma == 0:
        return np.zeros_like(values, dtype=float)
    return (values - mu) / sigma


def score_programs_and_root(
    counts: sparse.spmatrix,
    features: pd.Series,
    selected: pd.DataFrame,
    marker_programs: dict[str, list[str]],
    root_cluster: int,
    root_top_percent: float,
    root_min_cells: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    counts_sel = counts[:, selected["matrix_col_index"].to_numpy()].tocsc()
    lib_size = np.asarray(counts_sel.sum(axis=0)).ravel().astype(float)
    lib_size[lib_size <= 0] = np.nan
    scale = 1e4 / lib_size
    log_norm = counts_sel.multiply(scale)
    log_norm.data = np.log1p(log_norm.data)
    log_norm = log_norm.tocsr()

    gene_index = {gene.upper(): i for i, gene in enumerate(features.astype(str))}
    for program, genes in marker_programs.items():
        idx = [gene_index[g.upper()] for g in genes if g.upper() in gene_index]
        selected[f"jia_score_{program}"] = np.asarray(log_norm[idx, :].mean(axis=0)).ravel()
        selected[f"jia_score_{program}_n_genes"] = len(idx)

    marker_genes = ["HES1", "VIM", "NES", "DLX1", "DLX2", "ASCL1", "CACNA1E", "DCX"]
    for gene in marker_genes:
        idx = gene_index.get(gene.upper())
        if idx is None:
            selected[f"lognorm_{gene}"] = np.nan
        else:
            selected[f"lognorm_{gene}"] = np.asarray(log_norm[idx, :].todense()).ravel()

    components = {
        "z_RGC1": selected["jia_score_RGC1"].to_numpy(),
        "z_RGC2": selected["jia_score_RGC2"].to_numpy(),
        "z_HES1": selected["lognorm_HES1"].to_numpy(),
        "z_VIM": selected["lognorm_VIM"].to_numpy(),
        "z_NES": selected["lognorm_NES"].to_numpy(),
        "z_IPC": selected["jia_score_IPC"].to_numpy(),
        "z_DLX1": selected["lognorm_DLX1"].to_numpy(),
        "z_DLX2": selected["lognorm_DLX2"].to_numpy(),
        "z_ASCL1": selected["lognorm_ASCL1"].to_numpy(),
    }
    for name, values in components.items():
        selected[f"div90_rootscore_{name}"] = zscore(values)

    selected["div90_jia_rootscore"] = (
        selected["div90_rootscore_z_RGC1"]
        + selected["div90_rootscore_z_RGC2"]
        + selected["div90_rootscore_z_HES1"]
        + selected["div90_rootscore_z_VIM"]
        + selected["div90_rootscore_z_NES"]
        - selected["div90_rootscore_z_IPC"]
        - selected["div90_rootscore_z_DLX1"]
        - selected["div90_rootscore_z_DLX2"]
        - selected["div90_rootscore_z_ASCL1"]
    )

    root_pool = selected[selected["cluster_id_numeric"] == root_cluster].copy()
    n_root = max(root_min_cells, int(math.ceil(len(root_pool) * root_top_percent / 100.0)))
    n_root = min(n_root, len(root_pool))
    root_cells = set(root_pool.sort_values("div90_jia_rootscore", ascending=False).head(n_root)["cell_id"])
    selected["urd_root_candidate"] = selected["cell_id"].isin(root_cells)
    selected["div90_root_pool_cluster12"] = selected["cluster_id_numeric"] == root_cluster
    selected["div90_root_top_percent"] = np.where(selected["urd_root_candidate"], root_top_percent, np.nan)

    root_summary = pd.DataFrame(
        [
            {"metric": "root_cluster", "value": str(root_cluster)},
            {"metric": "root_top_percent", "value": str(root_top_percent)},
            {"metric": "root_min_cells", "value": str(root_min_cells)},
            {"metric": "n_root_pool_cells_selected", "value": str(len(root_pool))},
            {"metric": "n_root_cells", "value": str(n_root)},
            {"metric": "rootscore_min_selected_root", "value": str(selected.loc[selected["urd_root_candidate"], "div90_jia_rootscore"].min())},
            {"metric": "rootscore_median_selected_root", "value": str(selected.loc[selected["urd_root_candidate"], "div90_jia_rootscore"].median())},
            {"metric": "rootscore_max_selected_root", "value": str(selected.loc[selected["urd_root_candidate"], "div90_jia_rootscore"].max())},
        ]
    )
    return selected, root_summary


def add_lineage_roles(metadata: pd.DataFrame) -> pd.DataFrame:
    metadata["div90_jia_tip_group"] = metadata["cluster_id_numeric"].map(TIP_BY_CLUSTER)
    metadata["div90_jia_tip_display"] = metadata["div90_jia_tip_group"].map(TIP_DISPLAY)
    metadata["div90_jia_urd_role"] = "excluded_or_unassigned"
    metadata.loc[metadata["cluster_id_numeric"] == 12, "div90_jia_urd_role"] = "root_pool_cluster12_dividing_cells"
    metadata.loc[metadata["div90_jia_tip_group"].notna(), "div90_jia_urd_role"] = metadata.loc[
        metadata["div90_jia_tip_group"].notna(), "div90_jia_tip_group"
    ]
    metadata.loc[
        metadata["cluster_id_numeric"].isin(RETAINED_UNASSIGNED_CANDIDATE_CLUSTERS),
        "div90_jia_urd_role",
    ] = "retained_unassigned_candidate"
    metadata.loc[
        metadata["cluster_id_numeric"].isin(RETAINED_GLIA_NON_TIP_CLUSTERS),
        "div90_jia_urd_role",
    ] = "retained_glia_non_tip"
    metadata["paper_cluster_annotation"] = metadata["cluster_number_name"].astype(str)
    return metadata


def main() -> None:
    args = parse_args()
    project_root = args.project_root.expanduser().resolve()
    paths = default_paths(project_root, args.run_label)
    anndata_dir = (args.anndata_dir or paths["anndata_dir"]).expanduser().resolve()
    marker_path = (args.jia_marker_table or paths["jia_marker_table"]).expanduser().resolve()
    outdir = (args.outdir or paths["outdir"]).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    counts_path = anndata_dir / "matrix_counts.mtx"
    barcodes_path = anndata_dir / "barcodes.tsv"
    features_path = anndata_dir / "features.tsv"
    obs_path = anndata_dir / "obs.tsv"
    umap_path = anndata_dir / "umap.tsv"
    for path in [counts_path, barcodes_path, features_path, obs_path, umap_path, marker_path]:
        if not path.exists():
            raise FileNotFoundError(path)

    retain_clusters = parse_cluster_list(args.retain_clusters)
    print(f"[Div90JiaURDExport] counts={counts_path}", flush=True)
    print(f"[Div90JiaURDExport] obs={obs_path}", flush=True)
    print(f"[Div90JiaURDExport] markers={marker_path}", flush=True)
    print(f"[Div90JiaURDExport] retain_clusters={retain_clusters}", flush=True)
    print(f"[Div90JiaURDExport] outdir={outdir}", flush=True)

    barcodes = read_single_column_tsv(barcodes_path, "cell_id")
    features = read_single_column_tsv(features_path, "feature_id")
    obs = pd.read_csv(obs_path, sep="\t")
    umap = pd.read_csv(umap_path, sep="\t")
    metadata = obs.merge(umap, on="cell_id", how="inner", validate="one_to_one")
    metadata["cluster_id_numeric"] = pd.to_numeric(metadata["cluster_id"], errors="raise").astype(int)
    metadata = add_lineage_roles(metadata)

    barcode_index = pd.Series(np.arange(len(barcodes)), index=barcodes.values)
    metadata["matrix_col_index"] = metadata["cell_id"].map(barcode_index)
    if metadata["matrix_col_index"].isna().any():
        raise ValueError("Some DIV90 metadata cells were absent from barcodes.tsv")
    metadata["matrix_col_index"] = metadata["matrix_col_index"].astype(int)
    retained = metadata[metadata["cluster_id_numeric"].isin(retain_clusters)].copy()
    selected = allocate_stratified(retained, args.max_cells, args.seed, args.root_cluster)

    print("[Div90JiaURDExport] selected cluster counts:", flush=True)
    print(selected["cluster_number_name"].value_counts().sort_index().to_string(), flush=True)

    print("[Div90JiaURDExport] reading Matrix Market counts", flush=True)
    counts = mmread(counts_path).tocsc()
    if counts.shape != (len(features), len(barcodes)):
        raise ValueError(f"Counts shape {counts.shape} does not match features/barcodes ({len(features)}, {len(barcodes)})")

    marker_programs = load_program_markers(marker_path, set(features.astype(str)))
    selected, root_summary = score_programs_and_root(
        counts,
        features,
        selected,
        marker_programs,
        args.root_cluster,
        args.root_top_percent,
        args.root_min_cells,
    )
    if not bool(selected["urd_root_candidate"].any()):
        raise ValueError("RootScore selection produced zero root cells")

    selected_counts = counts[:, selected["matrix_col_index"].to_numpy()].tocoo()
    selected_counts_path = outdir / "div30_first_urd_counts.mtx"
    selected_features_path = outdir / "div30_first_urd_features.tsv"
    selected_barcodes_path = outdir / "div30_first_urd_barcodes.tsv"
    selected_metadata_path = outdir / "div30_first_urd_cell_metadata.tsv"
    manifest_path = outdir / "div30_first_urd_input_manifest.tsv"

    print(f"[Div90JiaURDExport] writing {selected_counts_path}", flush=True)
    mmwrite(selected_counts_path, selected_counts)
    pd.DataFrame({"feature_id": features}).to_csv(selected_features_path, sep="\t", index=False)
    pd.DataFrame({"cell_id": selected["cell_id"].to_numpy()}).to_csv(selected_barcodes_path, sep="\t", index=False)

    metadata_cols = [
        "cell_id",
        "orig.ident",
        "nCount_RNA",
        "nFeature_RNA",
        "percent.mt",
        "S.Score",
        "G2M.Score",
        "Phase",
        "old.ident",
        "CC.Difference",
        "RNA_snn_res.0.5",
        "seurat_clusters",
        "cluster_id",
        "cluster_id_numeric",
        "cluster_number_name",
        "paper_cluster_annotation",
        "div90_jia_tip_group",
        "div90_jia_tip_display",
        "div90_jia_urd_role",
        "jia_score_RGC1",
        "jia_score_RGC2",
        "jia_score_IPC",
        "jia_score_RGC1_n_genes",
        "jia_score_RGC2_n_genes",
        "jia_score_IPC_n_genes",
        "lognorm_HES1",
        "lognorm_VIM",
        "lognorm_NES",
        "lognorm_DLX1",
        "lognorm_DLX2",
        "lognorm_ASCL1",
        "lognorm_CACNA1E",
        "lognorm_DCX",
        "div90_jia_rootscore",
        "div90_rootscore_z_RGC1",
        "div90_rootscore_z_RGC2",
        "div90_rootscore_z_HES1",
        "div90_rootscore_z_VIM",
        "div90_rootscore_z_NES",
        "div90_rootscore_z_IPC",
        "div90_rootscore_z_DLX1",
        "div90_rootscore_z_DLX2",
        "div90_rootscore_z_ASCL1",
        "div90_root_pool_cluster12",
        "div90_root_top_percent",
        "UMAP_1",
        "UMAP_2",
        "urd_root_candidate",
        "urd_sampling_status",
        "matrix_col_index",
    ]
    export_cols = [col for col in metadata_cols if col in selected.columns]
    export_metadata = selected[export_cols].copy()
    for bool_col in ["urd_root_candidate", "div90_root_pool_cluster12"]:
        if bool_col in export_metadata.columns:
            export_metadata[bool_col] = export_metadata[bool_col].map({True: "TRUE", False: "FALSE"})
    export_metadata.to_csv(selected_metadata_path, sep="\t", index=False)

    root_summary.to_csv(outdir / "div90_jia_rootscore_root_summary.tsv", sep="\t", index=False)
    selected.sort_values("div90_jia_rootscore", ascending=False).to_csv(
        outdir / "div90_jia_rootscore_scored_metadata.tsv",
        sep="\t",
        index=False,
    )
    marker_rows = []
    for program, genes in marker_programs.items():
        for rank, gene in enumerate(genes, start=1):
            marker_rows.append({"program": program, "rank": rank, "gene": gene})
    pd.DataFrame(marker_rows).to_csv(outdir / "div90_jia_program_markers_used.tsv", sep="\t", index=False)

    manifest = pd.DataFrame(
        [
            {"key": "dataset", "value": "varela_div90"},
            {"key": "run_label", "value": args.run_label},
            {"key": "source_counts", "value": str(counts_path)},
            {"key": "source_barcodes", "value": str(barcodes_path)},
            {"key": "source_features", "value": str(features_path)},
            {"key": "source_obs", "value": str(obs_path)},
            {"key": "source_umap", "value": str(umap_path)},
            {"key": "source_jia_marker_table", "value": str(marker_path)},
            {"key": "matrix_orientation", "value": "features_by_cells"},
            {"key": "metadata_join_key", "value": "cell_id"},
            {"key": "retained_clusters", "value": ",".join(map(str, retain_clusters))},
            {"key": "excluded_clusters_first_smoke", "value": "6,7"},
            {"key": "cell_selection_logic", "value": "keep all root-cluster cells then stratify remaining retained clusters to max_cells"},
            {"key": "root_score_formula", "value": "z(RGC1)+z(RGC2)+z(HES1)+z(VIM)+z(NES)-z(IPC)-z(DLX1)-z(DLX2)-z(ASCL1)"},
            {"key": "root_cluster", "value": str(args.root_cluster)},
            {"key": "root_top_percent", "value": str(args.root_top_percent)},
            {"key": "root_min_cells", "value": str(args.root_min_cells)},
            {"key": "tip_lhx8_isl1_clusters", "value": "0,5,8"},
            {"key": "tip_lhx6_nfia_clusters", "value": "1"},
            {"key": "tip_crabp1_angpt2_clusters", "value": "2"},
            {"key": "retained_unassigned_candidate_clusters", "value": "3,11"},
            {"key": "retained_glia_non_tip_clusters", "value": "4,9,10"},
            {"key": "post_tree_candidate_marker_genes", "value": "MEF2C,EPHA5,LHX6,CRABP1,LHX8,NR2F1,NR2F2"},
            {"key": "max_cells", "value": str(args.max_cells)},
            {"key": "seed", "value": str(args.seed)},
            {"key": "n_selected_cells", "value": str(len(selected))},
            {"key": "n_root_cells", "value": str(int(selected["urd_root_candidate"].sum()))},
            {"key": "compatibility_note", "value": "Output bundle uses div30_first_urd_* file names for scripts/14_div30_first_urd.R compatibility."},
            {"key": "counts_mtx", "value": str(selected_counts_path)},
            {"key": "features_tsv", "value": str(selected_features_path)},
            {"key": "barcodes_tsv", "value": str(selected_barcodes_path)},
            {"key": "cell_metadata_tsv", "value": str(selected_metadata_path)},
        ]
    )
    manifest.to_csv(manifest_path, sep="\t", index=False)
    print(f"[Div90JiaURDExport] root_cells={int(selected['urd_root_candidate'].sum())}", flush=True)
    print(f"[Div90JiaURDExport] complete: {outdir}", flush=True)


if __name__ == "__main__":
    main()
