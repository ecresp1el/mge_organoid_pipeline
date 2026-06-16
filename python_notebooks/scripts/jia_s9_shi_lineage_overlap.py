#!/usr/bin/env python3
"""Compare Jia Science Data S9 modules with Shi et al. marker workbooks.

This is a workbook-only overlap workflow. It reads Jia full-module and TF-only
module sheets from science.adw1803_data_s9.xlsx, reads one or more Shi marker
tables, and writes overlap tables, crosswalk checks, plots, and a short report.

The biological crosswalk is intentionally kept separate from the computed
overlaps. If a Shi workbook does not contain the expected M2/M3/M5/M6/M7 labels,
the report marks those expected labels as missing instead of forcing the match.
"""

from __future__ import annotations

import argparse
import math
import os
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import hypergeom


PROJECT_ROOT_DEFAULT = "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder"
DEFAULT_JIA_XLSX = "reference/science.adw1803_data_s9.xlsx"
DEFAULT_SHI_XLSX = "reference/shi_2021_tables_s2_to_s9/science.abj6641_table_s6.xlsx"
DEFAULT_RUN_LABEL = "jia_s9_shi_table_s6_overlap_v1"
CURATED_M_LABELS = ["M2", "M3", "M4", "M5", "M6", "M7"]


JIA_MODULE_METADATA = {
    "Module 1": {
        "jia_module": "module_1",
        "jia_exact_name": "LHX8/ISL1",
        "collapsed_lineage": "Early subpallial cholinergic",
        "expected_shi_names": "M5/M6",
        "jia_logic": "VZ RGC-derived, early, adult subpallial cholinergic; ZIC4/HMGA1-associated.",
        "anchor_genes": ["LHX8", "ISL1", "GBX2", "ZIC1"],
    },
    "Module 2": {
        "jia_module": "module_2",
        "jia_exact_name": "NR2F1/NR2F2",
        "collapsed_lineage": "Early subpallial GABAergic",
        "expected_shi_names": "M4/M7",
        "jia_logic": "VZ RGC-derived, early, subpallial GABAergic.",
        "anchor_genes": ["NR2F1", "NR2F2", "ZFHX3"],
    },
    "Module 3": {
        "jia_module": "module_3",
        "jia_exact_name": "EPHA5/MEF2C",
        "collapsed_lineage": "Cortical MGE interneuron output",
        "expected_shi_names": "M2",
        "jia_logic": "SVZ RGC-derived, cortex-bound; broadly maps to adult cortical GABAergic types.",
        "anchor_genes": ["EPHA5", "MEF2C"],
    },
    "Module 4": {
        "jia_module": "module_4",
        "jia_exact_name": "LHX6/NFIA",
        "collapsed_lineage": "Cortical MGE interneuron output, chandelier-biased",
        "expected_shi_names": "M2",
        "jia_logic": "SVZ RGC-derived, cortex-bound; especially adult chandelier-associated in Jia.",
        "anchor_genes": ["LHX6", "NFIA"],
    },
    "Module 5": {
        "jia_module": "module_5",
        "jia_exact_name": "CRABP1/ANGPT2",
        "collapsed_lineage": "CRABP1 bridge / contested lineage",
        "expected_shi_names": "M3",
        "jia_logic": "SVZ RGC-derived subpallial GABAergic, EPHA5-low/non-DEN; ETV1-associated.",
        "anchor_genes": ["CRABP1", "ANGPT2", "ETV1", "NFIA"],
    },
}


@dataclass
class Paths:
    outdir: Path
    table_dir: Path
    plot_dir: Path
    report_dir: Path
    source_dir: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", default=os.environ.get("PROJECT_ROOT", PROJECT_ROOT_DEFAULT))
    parser.add_argument("--jia-xlsx", default=DEFAULT_JIA_XLSX)
    parser.add_argument(
        "--shi-xlsx",
        action="append",
        default=None,
        help="Shi workbook to compare. Repeatable. Defaults to Shi Table S6.",
    )
    parser.add_argument("--run-label", default=DEFAULT_RUN_LABEL)
    parser.add_argument("--outdir", default=None)
    parser.add_argument("--min-shi-avg-logfc", type=float, default=None)
    parser.add_argument("--max-shi-p-adj", type=float, default=None)
    parser.add_argument("--top-n-shi-per-cluster", type=int, default=None)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--copy-inputs", action="store_true")
    return parser.parse_args()


def resolve_path(path: str | Path, project_root: Path) -> Path:
    candidate = Path(path).expanduser()
    if not candidate.is_absolute():
        candidate = project_root / candidate
    return candidate.resolve()


def prepare_paths(args: argparse.Namespace, project_root: Path) -> Paths:
    outdir = (
        Path(args.outdir).expanduser().resolve()
        if args.outdir
        else project_root / "results" / "jia_s9_shi_lineage_overlap" / args.run_label
    )
    if outdir.exists() and any(outdir.iterdir()):
        if not args.overwrite:
            raise FileExistsError(f"Output directory is not empty; use --overwrite: {outdir}")
        for path in outdir.iterdir():
            if path.is_file() or path.is_symlink():
                path.unlink()
            else:
                shutil.rmtree(path)
    table_dir = outdir / "tables"
    plot_dir = outdir / "plots"
    report_dir = outdir / "reports"
    source_dir = outdir / "source"
    for directory in [table_dir, plot_dir, report_dir, source_dir]:
        directory.mkdir(parents=True, exist_ok=True)
    return Paths(outdir, table_dir, plot_dir, report_dir, source_dir)


def safe_token(value: object) -> str:
    token = re.sub(r"[^A-Za-z0-9]+", "_", str(value).strip())
    return re.sub(r"_+", "_", token).strip("_").lower() or "value"


def normalize_gene(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip().upper()


def is_gene_symbol(value: str) -> bool:
    """Return True for symbol-like gene values and False for numeric spillover."""
    if not value or value in {"NAN", "NONE"}:
        return False
    return bool(re.search(r"[A-Z]", value)) and bool(re.fullmatch(r"[A-Z0-9_.-]+", value))


def gene_from_row(row: pd.Series) -> tuple[str, str, int]:
    """Find the first symbol-like value in a Jia row.

    Module 5 contains a shifted block where symbols move from the `gene` column
    into the next column. This recovers those rows while recording provenance.
    """
    for idx, (column, value) in enumerate(row.items()):
        gene = normalize_gene(value)
        if is_gene_symbol(gene):
            return gene, str(column), idx
    return "", "", -1


def numeric_after(row: pd.Series, source_idx: int, offset: int, fallback_col: str) -> object:
    if source_idx == 0:
        return row.get(fallback_col, math.nan)
    target_idx = source_idx + offset
    if target_idx < len(row.index):
        return row.iloc[target_idx]
    return math.nan


def bh_adjust(p_values: Iterable[float]) -> list[float]:
    values = list(p_values)
    n = len(values)
    indexed = sorted(enumerate(values), key=lambda item: (math.inf if pd.isna(item[1]) else item[1]))
    adjusted = [math.nan] * n
    prev = 1.0
    for rank, (idx, p_val) in enumerate(reversed(indexed), start=1):
        if pd.isna(p_val):
            adjusted[idx] = math.nan
            continue
        original_rank = n - rank + 1
        q_val = min(prev, float(p_val) * n / original_rank)
        adjusted[idx] = min(q_val, 1.0)
        prev = q_val
    return adjusted


def read_jia_modules(path: Path) -> pd.DataFrame:
    rows = []
    xl = pd.ExcelFile(path)
    for sheet in xl.sheet_names:
        df = pd.read_excel(path, sheet_name=sheet)
        if "gene" not in df.columns:
            continue
        full_sheet = sheet.replace("（", "(").replace("）", ")")
        is_tf = "TF" in full_sheet.upper()
        base_sheet = re.sub(r"\s*\(TFS?\)\s*$", "", full_sheet, flags=re.IGNORECASE)
        meta = JIA_MODULE_METADATA.get(base_sheet)
        if meta is None:
            continue
        for idx, row in df.iterrows():
            gene, gene_source_column, gene_source_index = gene_from_row(row)
            if not gene:
                continue
            rows.append(
                {
                    "source": "Jia_Science_Data_S9",
                    "workbook": str(path),
                    "sheet": sheet,
                    "gene": gene,
                    "module_set": "tf_only" if is_tf else "full",
                    "jia_module": meta["jia_module"],
                    "jia_exact_name": meta["jia_exact_name"],
                    "collapsed_lineage": meta["collapsed_lineage"],
                    "expected_shi_names": meta["expected_shi_names"],
                    "gene_rank_in_sheet": idx + 1,
                    "gene_source_column": gene_source_column,
                    "pval": numeric_after(row, gene_source_index, 1, "pval"),
                    "qval": numeric_after(row, gene_source_index, 2, "qval"),
                    "vst_mean": numeric_after(row, gene_source_index, 3, "vst.mean"),
                    "vst_variance_standardized": row.get("vst.variance.standardized", math.nan),
                }
            )
    if not rows:
        raise ValueError(f"No Jia module genes found in {path}")
    return pd.DataFrame(rows).drop_duplicates(["sheet", "gene", "module_set"])


def detect_header_row(raw: pd.DataFrame) -> int | None:
    for idx, row in raw.iterrows():
        values = {str(value).strip().lower() for value in row.tolist() if not pd.isna(value)}
        if {"gene", "cluster"}.issubset(values):
            return int(idx)
    return None


def read_shi_workbook(path: Path, args: argparse.Namespace) -> pd.DataFrame:
    rows = []
    xl = pd.ExcelFile(path)
    for sheet in xl.sheet_names:
        raw = pd.read_excel(path, sheet_name=sheet, header=None)
        header_row = detect_header_row(raw)
        if header_row is None:
            continue
        df = pd.read_excel(path, sheet_name=sheet, header=header_row)
        df.columns = [str(col).strip() for col in df.columns]
        if "gene" not in df.columns or "cluster" not in df.columns:
            continue
        work = df.copy()
        work["gene"] = work["gene"].map(normalize_gene)
        work["cluster"] = work["cluster"].astype(str).str.strip()
        work = work.loc[work["gene"].ne("") & work["cluster"].ne("")]
        work = work.loc[work["gene"].map(is_gene_symbol)]
        if args.min_shi_avg_logfc is not None and "avg_logFC" in work.columns:
            work = work.loc[pd.to_numeric(work["avg_logFC"], errors="coerce") >= args.min_shi_avg_logfc]
        if args.max_shi_p_adj is not None and "p_val_adj" in work.columns:
            work = work.loc[pd.to_numeric(work["p_val_adj"], errors="coerce") <= args.max_shi_p_adj]
        if args.top_n_shi_per_cluster is not None:
            sort_cols = [col for col in ["cluster", "p_val_adj", "p_val", "avg_logFC"] if col in work.columns]
            ascending = [True] + [True if col != "avg_logFC" else False for col in sort_cols[1:]]
            work = work.sort_values(sort_cols, ascending=ascending).groupby("cluster", sort=False).head(args.top_n_shi_per_cluster)
        for idx, row in work.iterrows():
            rows.append(
                {
                    "source": "Shi_2021",
                    "workbook": str(path),
                    "workbook_name": path.name,
                    "sheet": sheet,
                    "cluster": str(row["cluster"]).strip(),
                    "shi_set_id": f"{path.stem}:{sheet}:{row['cluster']}",
                    "shi_table": path.stem.replace("science.abj6641_table_", "").upper(),
                    "gene": row["gene"],
                    "gene_rank_in_sheet": idx + 1,
                    "p_val": row.get("p_val", math.nan),
                    "avg_logFC": row.get("avg_logFC", math.nan),
                    "pct_1": row.get("pct.1", math.nan),
                    "pct_2": row.get("pct.2", math.nan),
                    "p_val_adj": row.get("p_val_adj", math.nan),
                }
            )
    if not rows:
        raise ValueError(f"No Shi marker genes found in {path}")
    return pd.DataFrame(rows).drop_duplicates(["workbook_name", "sheet", "cluster", "gene"])


def shi_context(row: pd.Series) -> str:
    workbook = str(row.get("workbook_name", ""))
    sheet = str(row.get("sheet", ""))
    cluster = str(row.get("cluster", ""))
    if "table_s3" in workbook:
        return "major_cell_type"
    if "table_s4" in workbook:
        return "progenitor_vs_postmitotic"
    if "table_s5" in workbook:
        return "rgc_vs_ipc"
    if "table_s6" in workbook:
        if cluster.startswith("pM"):
            return "ge_progenitor_subcluster_MGE"
        if cluster.startswith("pL"):
            return "ge_progenitor_subcluster_LGE"
        if cluster.startswith("pC"):
            return "ge_progenitor_subcluster_CGE"
        return "ge_progenitor_subcluster"
    if "table_s7" in workbook:
        return "branchpoint_1" if "branch point 1" in sheet else "branchpoint_2"
    if "table_s8" in workbook:
        return "lge_ob_vs_striatal"
    if "table_s9" in workbook:
        return "integrated_human_mouse_subcluster"
    return "shi_marker_set"


def add_shi_context(shi: pd.DataFrame) -> pd.DataFrame:
    out = shi.copy()
    out["shi_context"] = out.apply(shi_context, axis=1)
    out["shi_label"] = (
        out["workbook_name"].str.replace(".xlsx", "", regex=False)
        + " | "
        + out["sheet"].astype(str)
        + " | "
        + out["cluster"].astype(str)
    )
    return out


def set_summary(df: pd.DataFrame, group_cols: list[str], label: str) -> pd.DataFrame:
    rows = []
    for keys, group in df.groupby(group_cols, dropna=False, sort=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        row = dict(zip(group_cols, keys))
        genes = sorted(set(group["gene"]))
        row.update(
            {
                "set_label": label,
                "n_genes": len(genes),
                "genes": ",".join(genes),
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


def pairwise_overlap(jia: pd.DataFrame, shi: pd.DataFrame) -> pd.DataFrame:
    background = sorted(set(jia["gene"]) | set(shi["gene"]))
    background_size = len(background)
    rows = []
    jia_groups = [
        ("full", group)
        for _, group in jia.loc[jia["module_set"].eq("full")].groupby("jia_module", sort=False)
    ] + [
        ("tf_only", group)
        for _, group in jia.loc[jia["module_set"].eq("tf_only")].groupby("jia_module", sort=False)
    ]
    shi_groups = list(shi.groupby(["workbook_name", "sheet", "cluster", "shi_set_id", "shi_context"], sort=False))
    for _, jia_group in jia_groups:
        jia_genes = set(jia_group["gene"])
        jia_meta = jia_group.iloc[0].to_dict()
        for (workbook_name, sheet, cluster, shi_set_id, context), shi_group in shi_groups:
            shi_genes = set(shi_group["gene"])
            overlap = sorted(jia_genes & shi_genes)
            union = jia_genes | shi_genes
            p_val = hypergeom.sf(len(overlap) - 1, background_size, len(jia_genes), len(shi_genes))
            rows.append(
                {
                    "jia_module": jia_meta["jia_module"],
                    "jia_exact_name": jia_meta["jia_exact_name"],
                    "collapsed_lineage": jia_meta["collapsed_lineage"],
                    "jia_module_set": jia_meta["module_set"],
                    "expected_shi_names": jia_meta["expected_shi_names"],
                    "shi_workbook": workbook_name,
                    "shi_sheet": sheet,
                    "shi_cluster": cluster,
                    "shi_set_id": shi_set_id,
                    "shi_context": context,
                    "n_jia_genes": len(jia_genes),
                    "n_shi_genes": len(shi_genes),
                    "n_overlap": len(overlap),
                    "jaccard": len(overlap) / len(union) if union else 0.0,
                    "overlap_coefficient": len(overlap) / min(len(jia_genes), len(shi_genes)) if jia_genes and shi_genes else 0.0,
                    "fraction_jia_in_shi": len(overlap) / len(jia_genes) if jia_genes else 0.0,
                    "fraction_shi_in_jia": len(overlap) / len(shi_genes) if shi_genes else 0.0,
                    "hypergeom_p": p_val,
                    "overlap_genes": ",".join(overlap),
                }
            )
    out = pd.DataFrame(rows)
    out["hypergeom_q"] = bh_adjust(out["hypergeom_p"])
    return out.sort_values(["jia_module_set", "jia_module", "hypergeom_q", "n_overlap"], ascending=[True, True, True, False])


def module_internal_overlap(jia: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for module_set, subset in jia.groupby("module_set", sort=False):
        groups = list(subset.groupby(["jia_module", "jia_exact_name", "collapsed_lineage"], sort=False))
        for i, ((left_module, left_name, left_lineage), left_df) in enumerate(groups):
            left_genes = set(left_df["gene"])
            for right_module, right_name, right_lineage in [g[0] for g in groups[i + 1 :]]:
                right_df = subset.loc[subset["jia_module"].eq(right_module)]
                right_genes = set(right_df["gene"])
                overlap = sorted(left_genes & right_genes)
                union = left_genes | right_genes
                rows.append(
                    {
                        "module_set": module_set,
                        "left_module": left_module,
                        "left_name": left_name,
                        "left_lineage": left_lineage,
                        "right_module": right_module,
                        "right_name": right_name,
                        "right_lineage": right_lineage,
                        "n_left_genes": len(left_genes),
                        "n_right_genes": len(right_genes),
                        "n_overlap": len(overlap),
                        "jaccard": len(overlap) / len(union) if union else 0.0,
                        "overlap_coefficient": len(overlap) / min(len(left_genes), len(right_genes)) if left_genes and right_genes else 0.0,
                        "overlap_genes": ",".join(overlap),
                    }
                )
    return pd.DataFrame(rows)


def jia_gene_membership(jia: pd.DataFrame, shi: pd.DataFrame) -> pd.DataFrame:
    jia_gene_to_rows = {
        gene: group
        for gene, group in jia.groupby("gene", sort=True)
    }
    shi_gene_to_rows = {
        gene: group
        for gene, group in shi.groupby("gene", sort=True)
    }
    rows = []
    for gene in sorted(set(jia_gene_to_rows) | set(shi_gene_to_rows)):
        jia_rows = jia_gene_to_rows.get(gene, pd.DataFrame())
        shi_rows = shi_gene_to_rows.get(gene, pd.DataFrame())
        in_jia = not jia_rows.empty
        in_shi = not shi_rows.empty
        rows.append(
            {
                "gene": gene,
                "membership_class": "shared" if in_jia and in_shi else ("jia_only" if in_jia else "shi_only"),
                "in_jia": in_jia,
                "in_shi": in_shi,
                "jia_modules": ",".join(sorted(set(jia_rows.get("jia_module", [])))),
                "jia_names": ",".join(sorted(set(jia_rows.get("jia_exact_name", [])))),
                "jia_module_sets": ",".join(sorted(set(jia_rows.get("module_set", [])))),
                "shi_workbooks": ",".join(sorted(set(shi_rows.get("workbook_name", [])))),
                "shi_contexts": ",".join(sorted(set(shi_rows.get("shi_context", [])))),
                "shi_clusters": ",".join(sorted(set(shi_rows.get("cluster", [])))),
                "shi_set_ids": ";".join(sorted(set(shi_rows.get("shi_set_id", [])))),
            }
        )
    return pd.DataFrame(rows)


def module_presence_summary(jia: pd.DataFrame, shi: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    shi_genes = set(shi["gene"])
    rows = []
    shared_rows = []
    private_rows = []
    for keys, group in jia.groupby(["jia_module", "jia_exact_name", "module_set", "collapsed_lineage", "expected_shi_names"], sort=False):
        module, name, module_set, lineage, expected = keys
        genes = sorted(set(group["gene"]))
        shared = [gene for gene in genes if gene in shi_genes]
        private = [gene for gene in genes if gene not in shi_genes]
        rows.append(
            {
                "jia_module": module,
                "jia_exact_name": name,
                "module_set": module_set,
                "collapsed_lineage": lineage,
                "expected_shi_names": expected,
                "n_jia_genes": len(genes),
                "n_shared_with_any_shi_set": len(shared),
                "n_not_seen_in_any_shi_set": len(private),
                "fraction_shared_with_any_shi_set": len(shared) / len(genes) if genes else 0.0,
                "shared_genes": ",".join(shared),
                "not_seen_in_shi_genes": ",".join(private),
            }
        )
        shared_rows.extend([{"jia_module": module, "jia_exact_name": name, "module_set": module_set, "gene": gene} for gene in shared])
        private_rows.extend([{"jia_module": module, "jia_exact_name": name, "module_set": module_set, "gene": gene} for gene in private])
    return pd.DataFrame(rows), pd.DataFrame(shared_rows), pd.DataFrame(private_rows)


def shi_presence_summary(shi: pd.DataFrame, jia: pd.DataFrame) -> pd.DataFrame:
    jia_full = set(jia.loc[jia["module_set"].eq("full"), "gene"])
    jia_tf = set(jia.loc[jia["module_set"].eq("tf_only"), "gene"])
    rows = []
    for keys, group in shi.groupby(["workbook_name", "sheet", "cluster", "shi_set_id", "shi_context"], sort=False):
        workbook, sheet, cluster, set_id, context = keys
        genes = sorted(set(group["gene"]))
        full = [gene for gene in genes if gene in jia_full]
        tf = [gene for gene in genes if gene in jia_tf]
        rows.append(
            {
                "workbook_name": workbook,
                "sheet": sheet,
                "cluster": cluster,
                "shi_set_id": set_id,
                "shi_context": context,
                "n_shi_genes": len(genes),
                "n_in_jia_full_modules": len(full),
                "n_in_jia_tf_modules": len(tf),
                "fraction_in_jia_full_modules": len(full) / len(genes) if genes else 0.0,
                "fraction_in_jia_tf_modules": len(tf) / len(genes) if genes else 0.0,
                "genes_in_jia_full_modules": ",".join(full),
                "genes_in_jia_tf_modules": ",".join(tf),
            }
        )
    return pd.DataFrame(rows)


def marker_audit_table(jia: pd.DataFrame, shi: pd.DataFrame) -> pd.DataFrame:
    marker_genes = sorted({gene.upper() for meta in JIA_MODULE_METADATA.values() for gene in meta["anchor_genes"]})
    rows = []
    for gene in marker_genes:
        jia_hits = jia.loc[jia["gene"].eq(gene)]
        shi_hits = shi.loc[shi["gene"].eq(gene)]
        rows.append(
            {
                "marker_gene": gene,
                "in_jia": not jia_hits.empty,
                "jia_modules": ",".join(sorted(set(jia_hits["jia_module"]))),
                "jia_module_sets": ",".join(sorted(set(jia_hits["module_set"]))),
                "in_shi": not shi_hits.empty,
                "shi_contexts": ",".join(sorted(set(shi_hits["shi_context"]))),
                "shi_clusters": ",".join(sorted(set(shi_hits["cluster"]))),
                "shi_set_ids": ";".join(sorted(set(shi_hits["shi_set_id"]))),
            }
        )
    return pd.DataFrame(rows)


def convention_audit_table(shi: pd.DataFrame) -> pd.DataFrame:
    rows = []
    clusters = set(shi["cluster"].astype(str))
    all_genes = set(shi["gene"].astype(str))
    for label in CURATED_M_LABELS:
        rows.append(
            {
                "term": label,
                "term_role": "curated_shi_lineage_label_from_notes",
                "present_as_cluster_label": label in clusters,
                "present_as_gene_symbol": label in all_genes,
                "note": "Expected from narrative crosswalk, but not necessarily encoded as workbook cluster labels.",
            }
        )
    for cluster in sorted(clusters):
        if re.fullmatch(r"p[CLM][0-9]+", cluster):
            rows.append(
                {
                    "term": cluster,
                    "term_role": "observed_table_s6_progenitor_subcluster_label",
                    "present_as_cluster_label": True,
                    "present_as_gene_symbol": cluster in all_genes,
                    "note": "Computed overlap uses this actual workbook label when Table S6 is included.",
                }
            )
    return pd.DataFrame(rows)


def split_expected_shi_names(value: str) -> list[str]:
    tokens = re.split(r"[/,;|]\s*", str(value))
    return [token.strip() for token in tokens if token.strip()]


def crosswalk_table(jia: pd.DataFrame, shi: pd.DataFrame, overlap: pd.DataFrame) -> pd.DataFrame:
    rows = []
    shi_clusters = set(shi["cluster"].astype(str))
    shi_gene_sets = shi.groupby("cluster")["gene"].apply(lambda s: set(s)).to_dict()
    for base_sheet, meta in JIA_MODULE_METADATA.items():
        module = meta["jia_module"]
        full_genes = set(jia.loc[jia["jia_module"].eq(module) & jia["module_set"].eq("full"), "gene"])
        tf_genes = set(jia.loc[jia["jia_module"].eq(module) & jia["module_set"].eq("tf_only"), "gene"])
        expected = split_expected_shi_names(meta["expected_shi_names"])
        expected_present = [name for name in expected if name in shi_clusters]
        expected_missing = [name for name in expected if name not in shi_clusters]
        expected_shi_genes = set()
        for cluster in expected_present:
            expected_shi_genes.update(shi_gene_sets.get(cluster, set()))
        anchor_genes = [gene.upper() for gene in meta["anchor_genes"]]
        best = overlap.loc[overlap["jia_module"].eq(module) & overlap["jia_module_set"].eq("full")].head(1)
        best_by_context = []
        for context, ctx_df in overlap.loc[overlap["jia_module"].eq(module) & overlap["jia_module_set"].eq("full")].groupby("shi_context", sort=False):
            ctx_best = ctx_df.sort_values(["hypergeom_q", "n_overlap"], ascending=[True, False]).head(1)
            if not ctx_best.empty:
                best_by_context.append(f"{context}:{ctx_best.iloc[0]['shi_cluster']}({int(ctx_best.iloc[0]['n_overlap'])})")
        rows.append(
            {
                "collapsed_lineage": meta["collapsed_lineage"],
                "jia_exact_name": meta["jia_exact_name"],
                "jia_module": module,
                "jia_sheet": base_sheet,
                "expected_shi_names": meta["expected_shi_names"],
                "expected_shi_names_present_in_input": ",".join(expected_present),
                "expected_shi_names_missing_from_input": ",".join(expected_missing),
                "anchor_genes": ",".join(anchor_genes),
                "anchor_genes_in_jia_full": ",".join([gene for gene in anchor_genes if gene in full_genes]),
                "anchor_genes_in_jia_tf": ",".join([gene for gene in anchor_genes if gene in tf_genes]),
                "anchor_genes_in_expected_shi_sets": ",".join([gene for gene in anchor_genes if gene in expected_shi_genes]),
                "best_observed_shi_cluster": "" if best.empty else str(best.iloc[0]["shi_cluster"]),
                "best_observed_shi_context": "" if best.empty else str(best.iloc[0]["shi_context"]),
                "best_observed_shi_sheet": "" if best.empty else str(best.iloc[0]["shi_sheet"]),
                "best_observed_n_overlap": "" if best.empty else int(best.iloc[0]["n_overlap"]),
                "best_observed_overlap_coefficient": "" if best.empty else float(best.iloc[0]["overlap_coefficient"]),
                "best_observed_by_context": "; ".join(best_by_context),
                "jia_logic": meta["jia_logic"],
            }
        )
    return pd.DataFrame(rows)


def collapse_guidance_table(overlap: pd.DataFrame, internal: pd.DataFrame, presence: pd.DataFrame, crosswalk: pd.DataFrame) -> pd.DataFrame:
    rows = []
    full_internal = internal.loc[internal["module_set"].eq("full")].copy()
    module3_4 = full_internal.loc[
        ((full_internal["left_module"].eq("module_3") & full_internal["right_module"].eq("module_4"))
         | (full_internal["left_module"].eq("module_4") & full_internal["right_module"].eq("module_3")))
    ]
    module1_2 = full_internal.loc[
        ((full_internal["left_module"].eq("module_1") & full_internal["right_module"].eq("module_2"))
         | (full_internal["left_module"].eq("module_2") & full_internal["right_module"].eq("module_1")))
    ]
    def best(module: str) -> str:
        hit = overlap.loc[overlap["jia_module"].eq(module) & overlap["jia_module_set"].eq("full")].sort_values(["hypergeom_q", "n_overlap"], ascending=[True, False]).head(1)
        if hit.empty:
            return "none"
        row = hit.iloc[0]
        return f"{row['shi_context']} / {row['shi_cluster']} / {int(row['n_overlap'])} genes"

    def shared_fraction(module: str) -> str:
        row = presence.loc[presence["jia_module"].eq(module) & presence["module_set"].eq("full")]
        if row.empty:
            return ""
        row = row.iloc[0]
        return f"{int(row['n_shared_with_any_shi_set'])}/{int(row['n_jia_genes'])} ({float(row['fraction_shared_with_any_shi_set']):.1%})"

    rows.append(
        {
            "question": "Collapse EPHA5/MEF2C and LHX6/NFIA?",
            "modules": "module_3 + module_4",
            "curated_reason": "Both map to Shi M2/cortical MGE interneuron output in the notes.",
            "within_jia_overlap": "" if module3_4.empty else f"{int(module3_4.iloc[0]['n_overlap'])} genes; Jaccard {float(module3_4.iloc[0]['jaccard']):.3g}",
            "best_observed_shi_hits": f"module_3: {best('module_3')} | module_4: {best('module_4')}",
            "shi_coverage": f"module_3: {shared_fraction('module_3')} | module_4: {shared_fraction('module_4')}",
            "recommendation": "Collapse for a reader-facing cortical-output tier, but keep separate for Jia-resolution/chandelier-biased interpretation.",
        }
    )
    rows.append(
        {
            "question": "Collapse LHX8/ISL1 and NR2F1/NR2F2?",
            "modules": "module_1 + module_2",
            "curated_reason": "Both are early VZ-derived subpallial outputs, but one is cholinergic-biased and one GABAergic-biased.",
            "within_jia_overlap": "" if module1_2.empty else f"{int(module1_2.iloc[0]['n_overlap'])} genes; Jaccard {float(module1_2.iloc[0]['jaccard']):.3g}",
            "best_observed_shi_hits": f"module_1: {best('module_1')} | module_2: {best('module_2')}",
            "shi_coverage": f"module_1: {shared_fraction('module_1')} | module_2: {shared_fraction('module_2')}",
            "recommendation": "Do not collapse biologically; group under early subpallial outputs only as a parent tier.",
        }
    )
    rows.append(
        {
            "question": "Collapse CRABP1/ANGPT2 into subpallial-only?",
            "modules": "module_5",
            "curated_reason": "Jia calls it subpallial GABAergic; Shi notes support striatal plus cortical/PV-associated CRABP1 biology.",
            "within_jia_overlap": "see module_5 pairwise rows in jia_internal_module_overlap.tsv",
            "best_observed_shi_hits": f"module_5: {best('module_5')}",
            "shi_coverage": f"module_5: {shared_fraction('module_5')}",
            "recommendation": "Do not collapse into subpallial-only; present as a CRABP1 bridge/contested lineage.",
        }
    )
    return pd.DataFrame(rows)


def best_overlap_by_module(overlap: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "jia_module",
        "jia_exact_name",
        "collapsed_lineage",
        "jia_module_set",
        "expected_shi_names",
        "shi_workbook",
        "shi_sheet",
        "shi_cluster",
        "shi_context",
        "n_jia_genes",
        "n_shi_genes",
        "n_overlap",
        "jaccard",
        "overlap_coefficient",
        "fraction_jia_in_shi",
        "fraction_shi_in_jia",
        "hypergeom_q",
        "overlap_genes",
    ]
    ranked = overlap.sort_values(
        ["jia_module_set", "jia_module", "hypergeom_q", "n_overlap", "overlap_coefficient"],
        ascending=[True, True, True, False, False],
    )
    return ranked.groupby(["jia_module_set", "jia_module"], sort=False).head(1)[cols].reset_index(drop=True)


def plot_heatmap(overlap: pd.DataFrame, path: Path, module_set: str, value_col: str, title: str) -> None:
    data = overlap.loc[overlap["jia_module_set"].eq(module_set)].copy()
    if data.empty:
        return
    data["jia_label"] = data["jia_module"] + " " + data["jia_exact_name"]
    data["shi_label"] = data["shi_workbook"].str.replace(".xlsx", "", regex=False) + "\n" + data["shi_cluster"].astype(str)
    matrix = data.pivot_table(index="jia_label", columns="shi_label", values=value_col, aggfunc="max", fill_value=0)
    fig, ax = plt.subplots(figsize=(max(8, 0.45 * matrix.shape[1] + 4), max(4, 0.55 * matrix.shape[0] + 2)))
    sns.heatmap(matrix, cmap="viridis", linewidths=0.4, linecolor="white", ax=ax)
    ax.set_title(title)
    ax.set_xlabel("Shi marker set")
    ax.set_ylabel("Jia module")
    fig.tight_layout()
    fig.savefig(path, dpi=240)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_dotplot(overlap: pd.DataFrame, path: Path, module_set: str) -> None:
    data = overlap.loc[overlap["jia_module_set"].eq(module_set)].copy()
    if data.empty:
        return
    data["jia_label"] = data["jia_module"] + " " + data["jia_exact_name"]
    data["shi_label"] = data["shi_workbook"].str.replace(".xlsx", "", regex=False) + " | " + data["shi_cluster"].astype(str)
    data["minus_log10_q"] = data["hypergeom_q"].map(lambda q: -math.log10(max(float(q), 1e-300)) if not pd.isna(q) else 0.0)
    fig_width = max(8, 0.42 * data["shi_label"].nunique() + 4)
    fig, ax = plt.subplots(figsize=(fig_width, 5.5))
    sns.scatterplot(
        data=data,
        x="shi_label",
        y="jia_label",
        size="n_overlap",
        hue="minus_log10_q",
        sizes=(20, 260),
        palette="mako",
        edgecolor="black",
        linewidth=0.3,
        ax=ax,
    )
    ax.set_title(f"Jia {module_set} module overlap with Shi marker sets")
    ax.set_xlabel("Shi marker set")
    ax.set_ylabel("Jia module")
    ax.tick_params(axis="x", rotation=70)
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)
    fig.tight_layout()
    fig.savefig(path, dpi=240)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def simple_markdown_table(df: pd.DataFrame, max_rows: int = 12) -> str:
    if df.empty:
        return "_No rows._"
    display = df.head(max_rows).copy()
    for col in display.columns:
        if pd.api.types.is_float_dtype(display[col]):
            display[col] = display[col].map(lambda value: "" if pd.isna(value) else f"{value:.3g}")
        else:
            display[col] = display[col].map(
                lambda value: "" if pd.isna(value) else str(value).replace("|", "\\|").replace("\n", " ")
            )
    lines = [
        "| " + " | ".join(display.columns) + " |",
        "| " + " | ".join(["---"] * len(display.columns)) + " |",
    ]
    for row in display.to_numpy().tolist():
        lines.append("| " + " | ".join(row) + " |")
    return "\n".join(lines)


def write_report(paths: Paths, args: argparse.Namespace, jia: pd.DataFrame, shi: pd.DataFrame, overlap: pd.DataFrame, crosswalk: pd.DataFrame) -> None:
    presence, _, _ = module_presence_summary(jia, shi)
    internal = module_internal_overlap(jia)
    marker_audit = marker_audit_table(jia, shi)
    guidance = collapse_guidance_table(overlap, internal, presence, crosswalk)
    membership = jia_gene_membership(jia, shi)
    convention = convention_audit_table(shi)
    best_by_module = best_overlap_by_module(overlap)
    top_cols = [
        "jia_module",
        "jia_exact_name",
        "jia_module_set",
        "shi_workbook",
        "shi_cluster",
        "n_overlap",
        "overlap_coefficient",
        "hypergeom_q",
        "overlap_genes",
    ]
    top = overlap.sort_values(["jia_module_set", "jia_module", "hypergeom_q", "n_overlap"], ascending=[True, True, True, False])
    missing_expected = crosswalk.loc[crosswalk["expected_shi_names_missing_from_input"].astype(str).ne("")]
    shared_count = int((membership["membership_class"] == "shared").sum())
    jia_only_count = int((membership["membership_class"] == "jia_only").sum())
    shi_only_count = int((membership["membership_class"] == "shi_only").sum())
    curated_labels_present = int(
        convention.loc[
            convention["term_role"].eq("curated_shi_lineage_label_from_notes"),
            "present_as_cluster_label",
        ].sum()
    )
    lines = [
        "# Jia S9 vs Shi Marker Workbook Overlap",
        "",
        "## Purpose",
        "",
        "Compare Jia Science Data S9 full modules and TF-only modules with Shi et al. marker workbooks.",
        "The biological crosswalk is reported separately from the computed workbook overlap.",
        "",
        "## Executive Summary",
        "",
        f"- Jia module genes represented: `{jia['gene'].nunique()}` unique genes across full and TF-only modules.",
        f"- Shi marker genes represented: `{shi['gene'].nunique()}` unique genes across `{shi['shi_set_id'].nunique()}` marker sets.",
        f"- Shared Jia/Shi vocabulary: `{shared_count}` genes.",
        f"- Jia-only vocabulary in this comparison: `{jia_only_count}` genes.",
        f"- Shi-only vocabulary in this comparison: `{shi_only_count}` genes.",
        f"- Literal curated Shi M labels present as workbook cluster labels: `{curated_labels_present}` of `{len(CURATED_M_LABELS)}`.",
        "- The clean reader-facing collapse is Jia 5 -> Shi 4-ish, with EPHA5/MEF2C and LHX6/NFIA collapsed only at the cortical-output tier.",
        "- CRABP1/ANGPT2 should not be collapsed into subpallial-only; it remains the bridge/contested class.",
        "",
        "## Inputs",
        "",
        f"- Jia workbook: `{resolve_path(args.jia_xlsx, Path(args.project_root))}`",
        "- Shi workbook(s):",
        *[f"  - `{path}`" for path in args.shi_xlsx],
        "",
        "## Important Caveat",
        "",
        "When Shi Table S6 is included, it contains GE progenitor subcluster markers such as `pM1-pM4`.",
        "The user's biological crosswalk references Shi MGE lineage labels such as `M2`, `M3`, `M5/M6`, and `M4/M7`.",
        "If those labels are absent from the input workbook, this workflow still computes gene overlap but marks expected crosswalk labels as missing.",
        "In the S3-S9 exhaustive run, the curated M labels are not literal workbook cluster labels; the comparison is therefore a gene-vocabulary/conceptual crosswalk, not a direct same-label merge.",
        "",
        "## Crosswalk Check",
        "",
        simple_markdown_table(
            crosswalk[
                [
                    "collapsed_lineage",
                    "jia_exact_name",
                    "jia_module",
                    "expected_shi_names",
                    "expected_shi_names_present_in_input",
                    "expected_shi_names_missing_from_input",
                    "best_observed_shi_cluster",
                    "best_observed_n_overlap",
                ]
            ],
            max_rows=10,
        ),
        "",
        "## Best Observed Hit Per Jia Module",
        "",
        simple_markdown_table(
            best_by_module.loc[best_by_module["jia_module_set"].eq("full")][
                [
                    "jia_module",
                    "jia_exact_name",
                    "collapsed_lineage",
                    "shi_workbook",
                    "shi_context",
                    "shi_cluster",
                    "n_overlap",
                    "overlap_coefficient",
                    "fraction_jia_in_shi",
                    "hypergeom_q",
                ]
            ],
            max_rows=10,
        ),
        "",
        "## Shared vs Private Gene Counts",
        "",
        simple_markdown_table(
            presence[
                [
                    "jia_module",
                    "jia_exact_name",
                    "module_set",
                    "n_jia_genes",
                    "n_shared_with_any_shi_set",
                    "n_not_seen_in_any_shi_set",
                    "fraction_shared_with_any_shi_set",
                ]
            ],
            max_rows=12,
        ),
        "",
        "## Collapse Guidance",
        "",
        simple_markdown_table(guidance, max_rows=10),
        "",
        "## Anchor Marker Audit",
        "",
        simple_markdown_table(marker_audit, max_rows=20),
        "",
        "## Jia Internal Similarity",
        "",
        simple_markdown_table(
            internal.loc[internal["module_set"].eq("full")][
                ["left_module", "left_name", "right_module", "right_name", "n_overlap", "jaccard", "overlap_coefficient"]
            ],
            max_rows=20,
        ),
        "",
        "## Top Computed Overlaps",
        "",
        simple_markdown_table(top[top_cols], max_rows=20),
        "",
        "## Missing Expected Crosswalk Labels",
        "",
        simple_markdown_table(
            missing_expected[["jia_exact_name", "expected_shi_names", "expected_shi_names_missing_from_input"]],
            max_rows=10,
        ),
        "",
        "## Main Output Tables",
        "",
        "- `tables/jia_s9_module_genes.tsv`",
        "- `tables/shi_marker_genes_long.tsv`",
        "- `tables/jia_s9_shi_pairwise_overlap.tsv`",
        "- `tables/jia_shi_crosswalk_marker_presence.tsv`",
        "- `tables/jia_s9_gene_set_summary.tsv`",
        "- `tables/shi_gene_set_summary.tsv`",
        "- `tables/jia_gene_membership_all_sources.tsv`",
        "- `tables/jia_module_shared_private_summary.tsv`",
        "- `tables/jia_shared_genes_by_module.tsv`",
        "- `tables/jia_private_genes_by_module.tsv`",
        "- `tables/shi_marker_set_recovery_by_jia.tsv`",
        "- `tables/jia_internal_module_overlap.tsv`",
        "- `tables/jia_shi_anchor_marker_audit.tsv`",
        "- `tables/jia_shi_naming_convention_audit.tsv`",
        "- `tables/jia_shi_collapse_guidance.tsv`",
        "- `tables/jia_shi_best_observed_overlap_by_module.tsv`",
        "",
        "## Main Plots",
        "",
        "- `plots/jia_full_vs_shi_overlap_coefficient_heatmap.png`",
        "- `plots/jia_tf_vs_shi_overlap_coefficient_heatmap.png`",
        "- `plots/jia_full_vs_shi_overlap_dotplot.png`",
        "- `plots/jia_tf_vs_shi_overlap_dotplot.png`",
    ]
    (paths.report_dir / "jia_s9_shi_lineage_overlap_report.md").write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    project_root = Path(args.project_root).expanduser().resolve()
    args.jia_xlsx = str(resolve_path(args.jia_xlsx, project_root))
    args.shi_xlsx = [str(resolve_path(path, project_root)) for path in (args.shi_xlsx or [DEFAULT_SHI_XLSX])]
    paths = prepare_paths(args, project_root)

    jia_path = Path(args.jia_xlsx)
    shi_paths = [Path(path) for path in args.shi_xlsx]
    for path in [jia_path, *shi_paths]:
        if not path.exists():
            raise FileNotFoundError(path)

    if args.copy_inputs:
        shutil.copy2(jia_path, paths.source_dir / jia_path.name)
        for path in shi_paths:
            shutil.copy2(path, paths.source_dir / path.name)

    print(f"Reading Jia workbook: {jia_path}", flush=True)
    jia = read_jia_modules(jia_path)
    print(f"Reading {len(shi_paths)} Shi workbook(s)", flush=True)
    shi = pd.concat([read_shi_workbook(path, args) for path in shi_paths], ignore_index=True)
    shi = add_shi_context(shi)
    overlap = pairwise_overlap(jia, shi)
    crosswalk = crosswalk_table(jia, shi, overlap)
    gene_membership = jia_gene_membership(jia, shi)
    module_presence, shared_genes, private_genes = module_presence_summary(jia, shi)
    shi_recovery = shi_presence_summary(shi, jia)
    internal_overlap = module_internal_overlap(jia)
    marker_audit = marker_audit_table(jia, shi)
    convention_audit = convention_audit_table(shi)
    collapse_guidance = collapse_guidance_table(overlap, internal_overlap, module_presence, crosswalk)

    jia_summary = set_summary(jia, ["jia_module", "jia_exact_name", "module_set", "collapsed_lineage", "expected_shi_names"], "jia_module")
    shi_summary = set_summary(shi, ["workbook_name", "sheet", "cluster", "shi_set_id", "shi_context"], "shi_marker_set")

    jia.to_csv(paths.table_dir / "jia_s9_module_genes.tsv", sep="\t", index=False)
    shi.to_csv(paths.table_dir / "shi_marker_genes_long.tsv", sep="\t", index=False)
    overlap.to_csv(paths.table_dir / "jia_s9_shi_pairwise_overlap.tsv", sep="\t", index=False)
    crosswalk.to_csv(paths.table_dir / "jia_shi_crosswalk_marker_presence.tsv", sep="\t", index=False)
    jia_summary.to_csv(paths.table_dir / "jia_s9_gene_set_summary.tsv", sep="\t", index=False)
    shi_summary.to_csv(paths.table_dir / "shi_gene_set_summary.tsv", sep="\t", index=False)
    gene_membership.to_csv(paths.table_dir / "jia_gene_membership_all_sources.tsv", sep="\t", index=False)
    module_presence.to_csv(paths.table_dir / "jia_module_shared_private_summary.tsv", sep="\t", index=False)
    shared_genes.to_csv(paths.table_dir / "jia_shared_genes_by_module.tsv", sep="\t", index=False)
    private_genes.to_csv(paths.table_dir / "jia_private_genes_by_module.tsv", sep="\t", index=False)
    shi_recovery.to_csv(paths.table_dir / "shi_marker_set_recovery_by_jia.tsv", sep="\t", index=False)
    internal_overlap.to_csv(paths.table_dir / "jia_internal_module_overlap.tsv", sep="\t", index=False)
    marker_audit.to_csv(paths.table_dir / "jia_shi_anchor_marker_audit.tsv", sep="\t", index=False)
    convention_audit.to_csv(paths.table_dir / "jia_shi_naming_convention_audit.tsv", sep="\t", index=False)
    collapse_guidance.to_csv(paths.table_dir / "jia_shi_collapse_guidance.tsv", sep="\t", index=False)
    best_overlap_by_module(overlap).to_csv(paths.table_dir / "jia_shi_best_observed_overlap_by_module.tsv", sep="\t", index=False)

    plot_heatmap(
        overlap,
        paths.plot_dir / "jia_full_vs_shi_overlap_coefficient_heatmap.png",
        "full",
        "overlap_coefficient",
        "Jia full-module overlap coefficient vs Shi marker sets",
    )
    plot_heatmap(
        overlap,
        paths.plot_dir / "jia_tf_vs_shi_overlap_coefficient_heatmap.png",
        "tf_only",
        "overlap_coefficient",
        "Jia TF-module overlap coefficient vs Shi marker sets",
    )
    plot_dotplot(overlap, paths.plot_dir / "jia_full_vs_shi_overlap_dotplot.png", "full")
    plot_dotplot(overlap, paths.plot_dir / "jia_tf_vs_shi_overlap_dotplot.png", "tf_only")
    write_report(paths, args, jia, shi, overlap, crosswalk)

    manifest = pd.DataFrame(
        [
            {"key": "jia_xlsx", "value": str(jia_path)},
            {"key": "shi_xlsx", "value": ";".join(str(path) for path in shi_paths)},
            {"key": "outdir", "value": str(paths.outdir)},
            {"key": "min_shi_avg_logfc", "value": "" if args.min_shi_avg_logfc is None else str(args.min_shi_avg_logfc)},
            {"key": "max_shi_p_adj", "value": "" if args.max_shi_p_adj is None else str(args.max_shi_p_adj)},
            {"key": "top_n_shi_per_cluster", "value": "" if args.top_n_shi_per_cluster is None else str(args.top_n_shi_per_cluster)},
            {"key": "n_jia_gene_rows", "value": str(len(jia))},
            {"key": "n_shi_gene_rows", "value": str(len(shi))},
            {"key": "n_overlap_rows", "value": str(len(overlap))},
        ]
    )
    manifest.to_csv(paths.table_dir / "jia_s9_shi_overlap_run_manifest.tsv", sep="\t", index=False)
    print(f"Wrote Jia/Shi overlap workflow outputs to {paths.outdir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
