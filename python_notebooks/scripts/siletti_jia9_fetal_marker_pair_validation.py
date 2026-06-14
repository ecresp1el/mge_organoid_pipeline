#!/usr/bin/env python3
"""Validate adult Siletti/Jia-style groups against Jia fetal marker pairs.

This is a lightweight workbook-only workflow. It does not download WHB
expression matrices, submit Slurm jobs, or require a compute node.
"""

import argparse
import os
import sys
import pandas as pd
import numpy as np
import re
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import shutil
import textwrap
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, fcluster

# ============================================================
# Siletti/Linnarsson Excel → Jia-style 9 adult groups
# + fetal marker-pair validation layer
# ============================================================
#
# User-requested addition:
#   Validate the adult-group reconstruction by checking the marker pairs
#   used to name the 5 fetal clusters in the Jia figure:
#
#   EPHA5 / MEF2C
#   LHX6 / NFIA
#   CRABP1 / ANGPT2
#   NR2F1 / NR2F2
#   LHX8 / ISL1
#
# Important:
#   This is NOT label transfer.
#   This is a workbook-only marker-overlap audit using Siletti's
#   `Top enriched genes` field.
#
# Default input:
#   PROJECT_ROOT/results/siletti_2023_whb_reference_metadata/
#     siletti_2023_whb_metadata_inventory_v1/source/subcluster_annotation.xlsx
#
# Default output:
#   PROJECT_ROOT/results/siletti_2023_whb_reference_metadata/
#     siletti_jia9_fetal_marker_pair_validation_v1
# ============================================================

PROJECT_ROOT_DEFAULT = "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder"
DEFAULT_WORKBOOK_REL = (
    "results/siletti_2023_whb_reference_metadata/"
    "siletti_2023_whb_metadata_inventory_v1/source/subcluster_annotation.xlsx"
)
DEFAULT_OUTDIR_REL = (
    "results/siletti_2023_whb_reference_metadata/"
    "siletti_jia9_fetal_marker_pair_validation_v1"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", default=os.environ.get("PROJECT_ROOT", PROJECT_ROOT_DEFAULT))
    parser.add_argument("--workbook", default=None, help="Path to Linnarsson/Siletti subcluster_annotation.xlsx.")
    parser.add_argument("--outdir", default=None, help="Output run directory. Defaults under PROJECT_ROOT/results.")
    parser.add_argument(
        "--previous-folder",
        default=None,
        help="Optional prior Jia 9-group output folder to copy into source/ for continuity.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Delete an existing output directory before writing.",
    )
    parser.add_argument(
        "--make-zip",
        action="store_true",
        help="Also write a zip archive beside the output directory.",
    )
    return parser.parse_args()


args = parse_args()
project_root = Path(args.project_root)
xlsx_path = Path(args.workbook) if args.workbook else project_root / DEFAULT_WORKBOOK_REL
previous_folder = Path(args.previous_folder) if args.previous_folder else None
outdir = Path(args.outdir) if args.outdir else project_root / DEFAULT_OUTDIR_REL

if outdir.exists() and any(outdir.iterdir()):
    if not args.overwrite:
        raise FileExistsError(f"Output directory is not empty; use --overwrite to replace it: {outdir}")
    for p in outdir.iterdir():
        if p.is_file() or p.is_symlink():
            p.unlink()
        elif p.is_dir():
            shutil.rmtree(p)

outdir.mkdir(parents=True, exist_ok=True)
tables_dir = outdir / "tables"
plots_dir = outdir / "plots"
reports_dir = outdir / "reports"
source_dir = outdir / "source"
for directory in [tables_dir, plots_dir, reports_dir, source_dir]:
    directory.mkdir(parents=True, exist_ok=True)

# Copy prior v2 outputs for continuity, if explicitly supplied.
prev_out = source_dir / "previous_jia9_gene_signature_v2_outputs"
if previous_folder and previous_folder.exists():
    shutil.copytree(previous_folder, prev_out)

# -----------------------------
# Load workbook
# -----------------------------
df = pd.read_excel(xlsx_path)

# Expected workbook columns
super_col = "Supercluster"
cluster_col = "Cluster"
subcluster_col = "Subcluster"
mtg_col = "Transferred MTG Label (Transferred from cluster level)"
class_col = "Class"
nt_col = "Neurotransmitter"
roi_group_col = "Top ROIGroupFine"
roi_col = "Top ROI"
genes_col = "Top enriched genes"
cells_col = "Number of cells"

required = [
    subcluster_col, cluster_col, super_col, mtg_col, class_col,
    nt_col, roi_group_col, roi_col, genes_col, cells_col
]
missing = [c for c in required if c not in df.columns]
if missing:
    raise ValueError(f"Missing expected columns: {missing}")

# -----------------------------
# Helpers
# -----------------------------
def parse_genes(x):
    """Parse gene symbols from workbook Top enriched genes field."""
    if pd.isna(x):
        return []
    s = str(x)
    genes = re.findall(r"'([^']+)'", s)
    if len(genes) == 0:
        genes = re.split(r"[\s,;\[\]\n]+", s)
    genes = [g.strip().strip('"').strip("'") for g in genes]
    genes = [g for g in genes if len(g) > 0 and g.lower() != "nan"]
    return genes

def primary_before_colon(x):
    """Take first ranked ROI before ':'."""
    if pd.isna(x):
        return "unknown"
    s = str(x).strip()
    if len(s) == 0:
        return "unknown"
    return s.split(":")[0].strip()

def anatomy_bin_from_roi(primary_roi_group, primary_roi):
    """
    Explicit assumption:
      Hippocampus is grouped with pallium/cortex.

    CerebralCortex or Hippocampus -> Pallial/cortical
    Amygdala / striatum / pallidum / basal forebrain terms -> Subpallial
    Other -> Other/uncertain
    """
    if primary_roi_group in ["CerebralCortex", "Hippocampus"]:
        return "Pallial/cortical"

    roi_text = f"{primary_roi_group} {primary_roi}".lower()
    subpallial_terms = [
        "striat", "caudate", "putamen", "pallid", "globus", "accumbens",
        "basal", "sept", "amygdala", "claustr", "endopiriform",
        "substantia innominata", "bed nucleus"
    ]
    if any(term in roi_text for term in subpallial_terms):
        return "Subpallial"

    return "Other/uncertain"

def clean_gene_for_clustering(g):
    """Remove transcript-like IDs for clustering/visualization only."""
    g = str(g)
    if g.startswith(("AC", "AL", "AP")):
        return False
    if g.startswith("LINC"):
        return False
    if "-AS" in g:
        return False
    return True

def has_any(genes, markers):
    s = set(genes)
    return any(m in s for m in markers)

def has_gene(genes, marker):
    return marker in set(genes)

# -----------------------------
# Prepare workbook fields
# -----------------------------
work = df.copy()
work["gene_list"] = work[genes_col].apply(parse_genes)
work["clean_gene_list"] = work["gene_list"].apply(lambda L: [g for g in L if clean_gene_for_clustering(g)])
work["primary_roi_group"] = work[roi_group_col].apply(primary_before_colon)
work["primary_roi"] = work[roi_col].apply(primary_before_colon)
work["anatomy_bin"] = work.apply(lambda r: anatomy_bin_from_roi(r["primary_roi_group"], r["primary_roi"]), axis=1)
work["nt_string"] = work[nt_col].fillna("").astype(str)

# ============================================================
# Rebuild Jia-style 9 adult group assignments from v2
# ============================================================

# Split Sst parent clusters into 3 branches by Top enriched genes.
sst = work[work[mtg_col].eq("Sst")].copy()

rows = []
for cl, sub in sst.groupby(cluster_col):
    lists = sub["clean_gene_list"].tolist()
    n = len(lists)
    counts = {}
    for L in lists:
        for g in set(L):
            counts[g] = counts.get(g, 0) + 1
    for gene, c in counts.items():
        rows.append({
            "Cluster": int(cl),
            "gene": gene,
            "fraction_subclusters": c / n,
        })

sst_gene_freq = pd.DataFrame(rows)
if len(sst_gene_freq) > 0:
    gene_scores = (
        sst_gene_freq.groupby("gene")
        .agg(n_clusters=("Cluster", "nunique"), total_fraction=("fraction_subclusters", "sum"))
        .query("n_clusters >= 2")
        .sort_values("total_fraction", ascending=False)
    )
    top_genes = gene_scores.head(40).index.tolist()
    sst_mat = (
        sst_gene_freq[sst_gene_freq["gene"].isin(top_genes)]
        .pivot_table(index="Cluster", columns="gene", values="fraction_subclusters", fill_value=0, aggfunc="max")
    )
else:
    sst_mat = pd.DataFrame()

if sst_mat.shape[0] >= 3:
    dist = pdist(sst_mat.values, metric="jaccard")
    Z = linkage(dist, method="average")
    branch_ids = fcluster(Z, 3, criterion="maxclust")
    branch_map = pd.DataFrame({"Cluster": sst_mat.index.astype(int), "sst_branch_id": branch_ids})
else:
    branch_map = pd.DataFrame({"Cluster": sst[cluster_col].dropna().astype(int).unique(), "sst_branch_id": 1})

sst_branch_annot = sst.merge(branch_map, on="Cluster", how="left")
branch_summary = (
    sst_branch_annot.groupby("sst_branch_id")
    .agg(
        n_parent_clusters=("Cluster", "nunique"),
        n_subclusters=("Subcluster", "nunique"),
        cells=(cells_col, "sum"),
        parent_clusters=("Cluster", lambda x: ",".join(map(str, sorted(set(x.astype(int)))))),
        top_genes=("clean_gene_list", lambda lists: ", ".join(pd.Series([g for L in lists for g in L]).value_counts().head(18).index.astype(str))),
        top_roi_groups=("primary_roi_group", lambda x: "; ".join(x.value_counts().head(3).index.astype(str))),
    )
    .reset_index()
    .sort_values("cells", ascending=False)
)

ordered_branches = list(branch_summary["sst_branch_id"])
label_by_branch = {}
if len(ordered_branches) >= 1:
    label_by_branch[ordered_branches[0]] = "Cortical SST+ Mt neurons"
if len(ordered_branches) >= 2:
    label_by_branch[ordered_branches[1]] = "Cortical SST+ nMt neurons"
if len(ordered_branches) >= 3:
    label_by_branch[ordered_branches[2]] = "Cortical SST+ LRP neurons"
branch_summary["candidate_jia_label"] = branch_summary["sst_branch_id"].map(label_by_branch)
branch_summary.to_csv(tables_dir / "sst_branch_model_from_workbook.csv", index=False)

def is_sst_lrp_like(genes):
    s = set(genes)
    return ("NOS1" in s) or (("NPY" in s) and ("CORT" in s))

def assign_candidate_jia_group(row):
    anatomy = row["anatomy_bin"]
    mtg = row[mtg_col]
    sup = row[super_col]
    genes = row["gene_list"]
    nt = row["nt_string"]

    # Subpallial cholinergic first
    if anatomy == "Subpallial" and ("NT-CHOL" in nt) and has_any(genes, ["CHAT", "SLC5A7"]):
        return "Subpallial Cholinergic neurons"

    # Pallial/cortical side
    if anatomy == "Pallial/cortical":
        if mtg == "Chandelier":
            return "Cortical PV+ Chandelier neurons"
        if mtg == "Pvalb" and sup == "MGE interneuron":
            return "Cortical PV+ basket neurons"
        if mtg == "Sst" and sup == "MGE interneuron":
            cl = int(row[cluster_col])
            match = branch_map.loc[branch_map["Cluster"].eq(cl), "sst_branch_id"]
            if len(match) > 0:
                return label_by_branch.get(int(match.iloc[0]), "Cortical SST+ unassigned neurons")
            return "Cortical SST+ unassigned neurons"

    # Subpallial side
    if anatomy == "Subpallial":
        if (mtg == "Pvalb") or has_gene(genes, "PVALB"):
            return "Subpallial PV+ neurons"
        if (mtg == "Sst") or has_gene(genes, "SST"):
            if is_sst_lrp_like(genes):
                return "Subpallial SST+ LRP neurons"
            return "Subpallial SST+ neurons"

    return "Excluded / not assigned to Jia-style 9 groups"

work["candidate_jia_group"] = work.apply(assign_candidate_jia_group, axis=1)
work["is_candidate_jia_group"] = ~work["candidate_jia_group"].eq("Excluded / not assigned to Jia-style 9 groups")

def candidate_side(g):
    if g.startswith("Cortical"):
        return "Cortical GABAergic inhibitory neurons"
    if g.startswith("Subpallial"):
        return "Subpallial inhibitory neurons"
    return "Excluded"

work["jia_side"] = work["candidate_jia_group"].apply(candidate_side)
selected = work[work["is_candidate_jia_group"]].copy()

# ============================================================
# Fetal marker pair layer
# ============================================================

fetal_signatures = {
    "EPHA5/MEF2C fetal neurons": ["EPHA5", "MEF2C"],
    "LHX6/NFIA fetal neurons": ["LHX6", "NFIA"],
    "CRABP1/ANGPT2 fetal neurons": ["CRABP1", "ANGPT2"],
    "NR2F1/NR2F2 fetal neurons": ["NR2F1", "NR2F2"],
    "LHX8/ISL1 fetal neurons": ["LHX8", "ISL1"],
}

fetal_signature_def = pd.DataFrame([
    {
        "fetal_label_from_jia_figure": k,
        "gene_1": v[0],
        "gene_2": v[1],
        "n_genes": len(v),
    }
    for k, v in fetal_signatures.items()
])
fetal_signature_def.to_csv(tables_dir / "jia_fetal_marker_pair_definitions.csv", index=False)

# Score row-level pair overlap.
for fetal_label, genes in fetal_signatures.items():
    safe = re.sub(r"[^A-Za-z0-9]+", "_", fetal_label).strip("_")
    hit_col = f"fetal_pair_hits__{safe}"
    n_hit_col = f"fetal_pair_n_hits__{safe}"
    score_col = f"fetal_pair_score__{safe}"
    both_col = f"fetal_pair_both_genes__{safe}"
    either_col = f"fetal_pair_either_gene__{safe}"

    def hits(L, genes=genes):
        s = set(L)
        return [g for g in genes if g in s]

    work[hit_col] = work["gene_list"].apply(lambda L: ",".join(hits(L)))
    work[n_hit_col] = work["gene_list"].apply(lambda L: len(hits(L)))
    work[score_col] = work[n_hit_col] / len(genes)
    work[both_col] = work[n_hit_col].eq(len(genes))
    work[either_col] = work[n_hit_col].gt(0)

    selected[hit_col] = selected["gene_list"].apply(lambda L: ",".join(hits(L)))
    selected[n_hit_col] = selected["gene_list"].apply(lambda L: len(hits(L)))
    selected[score_col] = selected[n_hit_col] / len(genes)
    selected[both_col] = selected[n_hit_col].eq(len(genes))
    selected[either_col] = selected[n_hit_col].gt(0)

score_cols = [f"fetal_pair_score__{re.sub(r'[^A-Za-z0-9]+', '_', k).strip('_')}" for k in fetal_signatures]
n_hit_cols = [f"fetal_pair_n_hits__{re.sub(r'[^A-Za-z0-9]+', '_', k).strip('_')}" for k in fetal_signatures]
either_cols = [f"fetal_pair_either_gene__{re.sub(r'[^A-Za-z0-9]+', '_', k).strip('_')}" for k in fetal_signatures]
both_cols = [f"fetal_pair_both_genes__{re.sub(r'[^A-Za-z0-9]+', '_', k).strip('_')}" for k in fetal_signatures]
score_to_fetal = dict(zip(score_cols, fetal_signatures.keys()))

# Best fetal pair per row.
if len(score_cols) > 0:
    score_matrix = work[score_cols].copy()
    work["best_fetal_pair_score"] = score_matrix.max(axis=1)
    work["best_fetal_pair"] = score_matrix.idxmax(axis=1).map(score_to_fetal)
    # If all scores are zero, make it explicit.
    work.loc[work["best_fetal_pair_score"].eq(0), "best_fetal_pair"] = "No fetal marker-pair hit"

    score_matrix_sel = selected[score_cols].copy()
    selected["best_fetal_pair_score"] = score_matrix_sel.max(axis=1)
    selected["best_fetal_pair"] = score_matrix_sel.idxmax(axis=1).map(score_to_fetal)
    selected.loc[selected["best_fetal_pair_score"].eq(0), "best_fetal_pair"] = "No fetal marker-pair hit"

# ============================================================
# QC / summary tables
# ============================================================

stage_qc = pd.DataFrame([
    {
        "stage": "0_all_workbook",
        "rows": len(work),
        "unique_subclusters": work[subcluster_col].nunique(),
        "unique_parent_clusters": work[cluster_col].nunique(),
        "cells_sum": int(work[cells_col].sum()),
    },
    {
        "stage": "1_candidate_jia_9_groups",
        "rows": len(selected),
        "unique_subclusters": selected[subcluster_col].nunique(),
        "unique_parent_clusters": selected[cluster_col].nunique(),
        "cells_sum": int(selected[cells_col].sum()),
    },
    {
        "stage": "2_excluded",
        "rows": int((~work["is_candidate_jia_group"]).sum()),
        "unique_subclusters": work.loc[~work["is_candidate_jia_group"], subcluster_col].nunique(),
        "unique_parent_clusters": work.loc[~work["is_candidate_jia_group"], cluster_col].nunique(),
        "cells_sum": int(work.loc[~work["is_candidate_jia_group"], cells_col].sum()),
    },
])
stage_qc.to_csv(tables_dir / "jia9_fetal_marker_stage_qc_totals.csv", index=False)

# Fetal marker occurrence across entire workbook.
gene_occurrence_rows = []
for fetal_label, genes in fetal_signatures.items():
    for gene in genes:
        mask = work["gene_list"].apply(lambda L, g=gene: g in set(L))
        gene_occurrence_rows.append({
            "fetal_label": fetal_label,
            "gene": gene,
            "all_rows_with_gene": int(mask.sum()),
            "all_cells_with_gene": int(work.loc[mask, cells_col].sum()),
            "selected_rows_with_gene": int(mask[work["is_candidate_jia_group"]].sum()),
            "selected_cells_with_gene": int(work.loc[mask & work["is_candidate_jia_group"], cells_col].sum()),
            "top_superclusters_all": "; ".join(work.loc[mask, super_col].value_counts().head(8).index.astype(str)),
            "top_candidate_groups_selected": "; ".join(work.loc[mask & work["is_candidate_jia_group"], "candidate_jia_group"].value_counts().head(8).index.astype(str)),
        })
gene_occurrence = pd.DataFrame(gene_occurrence_rows)
gene_occurrence.to_csv(tables_dir / "jia_fetal_pair_gene_occurrence_in_workbook.csv", index=False)

# Adult candidate group × fetal pair matrix.
group_fetal_rows = []
for group, sub in selected.groupby("candidate_jia_group"):
    row = {
        "candidate_jia_group": group,
        "rows": len(sub),
        "unique_subclusters": sub[subcluster_col].nunique(),
        "unique_parent_clusters": sub[cluster_col].nunique(),
        "cells_sum": int(sub[cells_col].sum()),
    }
    for fetal_label in fetal_signatures:
        safe = re.sub(r"[^A-Za-z0-9]+", "_", fetal_label).strip("_")
        score_col = f"fetal_pair_score__{safe}"
        n_hit_col = f"fetal_pair_n_hits__{safe}"
        either_col = f"fetal_pair_either_gene__{safe}"
        both_col = f"fetal_pair_both_genes__{safe}"

        row[f"{fetal_label}__mean_pair_score"] = sub[score_col].mean()
        row[f"{fetal_label}__fraction_rows_either_gene"] = sub[either_col].mean()
        row[f"{fetal_label}__fraction_rows_both_genes"] = sub[both_col].mean()
        row[f"{fetal_label}__cells_weighted_pair_score"] = float((sub[score_col] * sub[cells_col]).sum() / sub[cells_col].sum()) if sub[cells_col].sum() > 0 else 0
        row[f"{fetal_label}__evidence_cell_mass"] = float((sub[score_col] * sub[cells_col]).sum())
    group_fetal_rows.append(row)

group_fetal_summary = pd.DataFrame(group_fetal_rows).sort_values("cells_sum", ascending=False)
group_fetal_summary.to_csv(tables_dir / "adult_candidate_group_by_fetal_marker_pair_summary.csv", index=False)

# Long form group x fetal
long_rows = []
for _, row in group_fetal_summary.iterrows():
    group = row["candidate_jia_group"]
    for fetal_label in fetal_signatures:
        long_rows.append({
            "candidate_jia_group": group,
            "fetal_label": fetal_label,
            "rows": int(row["rows"]),
            "cells_sum": int(row["cells_sum"]),
            "mean_pair_score": row[f"{fetal_label}__mean_pair_score"],
            "fraction_rows_either_gene": row[f"{fetal_label}__fraction_rows_either_gene"],
            "fraction_rows_both_genes": row[f"{fetal_label}__fraction_rows_both_genes"],
            "cells_weighted_pair_score": row[f"{fetal_label}__cells_weighted_pair_score"],
            "evidence_cell_mass": row[f"{fetal_label}__evidence_cell_mass"],
        })
group_fetal_long = pd.DataFrame(long_rows)
group_fetal_long.to_csv(tables_dir / "adult_candidate_group_by_fetal_marker_pair_long.csv", index=False)

# Parent cluster × fetal pair.
selected["source_cluster"] = selected.apply(lambda r: f"{r[super_col]} / {int(r[cluster_col])}", axis=1)
cluster_rows = []
for (group, source_cluster), sub in selected.groupby(["candidate_jia_group", "source_cluster"]):
    row = {
        "candidate_jia_group": group,
        "source_cluster": source_cluster,
        "rows": len(sub),
        "cells_sum": int(sub[cells_col].sum()),
    }
    for fetal_label in fetal_signatures:
        safe = re.sub(r"[^A-Za-z0-9]+", "_", fetal_label).strip("_")
        score_col = f"fetal_pair_score__{safe}"
        row[f"{fetal_label}__mean_pair_score"] = sub[score_col].mean()
        row[f"{fetal_label}__cells_weighted_pair_score"] = float((sub[score_col] * sub[cells_col]).sum() / sub[cells_col].sum()) if sub[cells_col].sum() > 0 else 0
        row[f"{fetal_label}__evidence_cell_mass"] = float((sub[score_col] * sub[cells_col]).sum())
    cluster_rows.append(row)

cluster_fetal_summary = pd.DataFrame(cluster_rows).sort_values(["candidate_jia_group", "source_cluster"])
cluster_fetal_summary.to_csv(tables_dir / "parent_cluster_by_fetal_marker_pair_summary.csv", index=False)

# Best fetal pair by candidate group
best_by_group = []
for group, sub in group_fetal_long.groupby("candidate_jia_group"):
    tmp = sub.sort_values("cells_weighted_pair_score", ascending=False)
    best = tmp.iloc[0]
    best_by_group.append({
        "candidate_jia_group": group,
        "best_fetal_marker_pair": best["fetal_label"],
        "best_cells_weighted_pair_score": best["cells_weighted_pair_score"],
        "best_evidence_cell_mass": best["evidence_cell_mass"],
        "best_fraction_rows_either_gene": best["fraction_rows_either_gene"],
        "best_fraction_rows_both_genes": best["fraction_rows_both_genes"],
        "cells_sum": best["cells_sum"],
    })
best_by_group = pd.DataFrame(best_by_group).sort_values("best_cells_weighted_pair_score", ascending=False)
best_by_group.to_csv(tables_dir / "best_fetal_marker_pair_by_adult_candidate_group.csv", index=False)

# Save row-level tables
base_cols = [
    subcluster_col, cluster_col, super_col, mtg_col, nt_col,
    roi_group_col, roi_col, cells_col,
    "primary_roi_group", "primary_roi", "anatomy_bin",
    "jia_side", "candidate_jia_group", "is_candidate_jia_group",
    "best_fetal_pair", "best_fetal_pair_score",
    genes_col
]
work[base_cols + score_cols + n_hit_cols + either_cols + both_cols].to_csv(tables_dir / "row_level_fetal_marker_pair_scores_all_rows.csv", index=False)
selected[base_cols + score_cols + n_hit_cols + either_cols + both_cols].to_csv(tables_dir / "row_level_fetal_marker_pair_scores_selected_only.csv", index=False)

# ============================================================
# Plots
# ============================================================

def draw_river(df_edges, levels, value_col, title, outpath, figsize=(20, 12), min_flow_to_label=0):
    level_values = []
    level_y = []
    for lev in levels:
        vals = df_edges.groupby(lev, dropna=False)[value_col].sum().sort_values(ascending=False)
        labels = list(vals.index.astype(str))
        y_positions = {}
        if len(labels) == 1:
            y_positions[labels[0]] = 0.5
        else:
            for lab, yy in zip(labels, np.linspace(0.94, 0.06, len(labels))):
                y_positions[lab] = yy
        level_values.append(vals)
        level_y.append(y_positions)

    xs = np.linspace(0.04, 0.88, len(levels))
    fig_h = max(figsize[1], 0.35 * max(len(v) for v in level_values))
    fig, ax = plt.subplots(figsize=(figsize[0], fig_h))
    ax.axis("off")

    max_flow = max(df_edges[value_col].max(), 1)
    def lw(v):
        return 0.4 + 13 * (v / max_flow)

    for i in range(len(levels) - 1):
        grouped = (
            df_edges.groupby([levels[i], levels[i+1]], dropna=False)[value_col]
            .sum()
            .reset_index()
        )
        for _, r in grouped.iterrows():
            a = str(r[levels[i]])
            b = str(r[levels[i+1]])
            v = r[value_col]
            if v <= 0:
                continue
            ax.plot(
                [xs[i] + 0.13, xs[i+1] - 0.035],
                [level_y[i][a], level_y[i+1][b]],
                linewidth=lw(v),
                alpha=0.30,
                solid_capstyle="round",
            )

    for i, lev in enumerate(levels):
        vals = level_values[i]
        total = vals.sum()
        ax.text(xs[i], 0.995, lev, ha="left", va="bottom", fontsize=12, fontweight="bold")
        for lab, v in vals.items():
            if v < min_flow_to_label:
                continue
            lab_str = str(lab)
            pct = 100 * v / total if total else 0
            fs = 8 if len(vals) > 14 else 10
            ax.text(
                xs[i],
                level_y[i][lab_str],
                f"{lab_str}\n{v:,.1f} evidence ({pct:.1f}%)",
                ha="left",
                va="center",
                fontsize=fs,
            )

    ax.set_title(title, fontsize=14, pad=18)
    fig.savefig(outpath, dpi=220, bbox_inches="tight")
    plt.close(fig)

# Plot 1: filter accounting
fig, ax = plt.subplots(figsize=(12, 5))
ax.axis("off")
x_positions = [0.08, 0.48, 0.86]
for x, (_, row) in zip(x_positions, stage_qc.iterrows()):
    text = (
        f"{row['stage']}\n\n"
        f"rows: {int(row['rows']):,}\n"
        f"subclusters: {int(row['unique_subclusters']):,}\n"
        f"parent clusters: {int(row['unique_parent_clusters']):,}\n"
        f"cell-count sum: {int(row['cells_sum']):,}"
    )
    ax.text(
        x, 0.55, text, ha="center", va="center", fontsize=11,
        bbox=dict(boxstyle="round,pad=0.5", fc="white", ec="black", lw=1)
    )
ax.annotate("", xy=(0.36, 0.55), xytext=(0.20, 0.55), arrowprops=dict(arrowstyle="->", lw=2))
ax.annotate("", xy=(0.74, 0.55), xytext=(0.60, 0.55), arrowprops=dict(arrowstyle="->", lw=2))
ax.set_title(
    "Workbook-to-Jia-style 9-group filter accounting\n"
    "Fetal marker-pair validation added after adult group assignment",
    fontsize=14,
)
qc_plot = plots_dir / "fetal_marker_pair_filter_accounting.png"
fig.savefig(qc_plot, dpi=220, bbox_inches="tight")
plt.close(fig)

# Plot 2: adult group x fetal pair heatmap, cell-weighted score
heat = group_fetal_long.pivot_table(
    index="candidate_jia_group",
    columns="fetal_label",
    values="cells_weighted_pair_score",
    fill_value=0
)
# order adult groups like Jia-ish
adult_order = [
    "Cortical SST+ LRP neurons",
    "Cortical SST+ nMt neurons",
    "Cortical SST+ Mt neurons",
    "Cortical PV+ basket neurons",
    "Cortical PV+ Chandelier neurons",
    "Subpallial SST+ LRP neurons",
    "Subpallial SST+ neurons",
    "Subpallial PV+ neurons",
    "Subpallial Cholinergic neurons",
]
heat = heat.reindex([g for g in adult_order if g in heat.index])
heat = heat[list(fetal_signatures.keys())]

fig, ax = plt.subplots(figsize=(10, 7))
im = ax.imshow(heat.values, aspect="auto", interpolation="nearest")
ax.set_xticks(np.arange(heat.shape[1]))
ax.set_xticklabels(heat.columns, rotation=45, ha="right", fontsize=8)
ax.set_yticks(np.arange(heat.shape[0]))
ax.set_yticklabels(heat.index, fontsize=8)
ax.set_title(
    "Adult candidate group × Jia fetal marker-pair evidence\n"
    "Value = cell-weighted pair score from Siletti Top enriched genes"
)
cbar = fig.colorbar(im, ax=ax, shrink=0.8)
cbar.set_label("Cell-weighted pair score")
for i in range(heat.shape[0]):
    for j in range(heat.shape[1]):
        ax.text(j, i, f"{heat.values[i,j]:.2f}", ha="center", va="center", fontsize=7)
group_fetal_heatmap = plots_dir / "adult_group_by_jia_fetal_marker_pair_heatmap.png"
fig.savefig(group_fetal_heatmap, dpi=220, bbox_inches="tight")
plt.close(fig)

# Plot 3: fraction rows with either gene
frac_heat = group_fetal_long.pivot_table(
    index="candidate_jia_group",
    columns="fetal_label",
    values="fraction_rows_either_gene",
    fill_value=0
).reindex(heat.index)[list(fetal_signatures.keys())]

fig, ax = plt.subplots(figsize=(10, 7))
im = ax.imshow(frac_heat.values, aspect="auto", interpolation="nearest")
ax.set_xticks(np.arange(frac_heat.shape[1]))
ax.set_xticklabels(frac_heat.columns, rotation=45, ha="right", fontsize=8)
ax.set_yticks(np.arange(frac_heat.shape[0]))
ax.set_yticklabels(frac_heat.index, fontsize=8)
ax.set_title(
    "Adult candidate group × Jia fetal marker-pair occurrence\n"
    "Value = fraction of selected rows with either marker from the pair"
)
cbar = fig.colorbar(im, ax=ax, shrink=0.8)
cbar.set_label("Fraction rows with either gene")
for i in range(frac_heat.shape[0]):
    for j in range(frac_heat.shape[1]):
        ax.text(j, i, f"{frac_heat.values[i,j]:.2f}", ha="center", va="center", fontsize=7)
group_fetal_fraction_heatmap = plots_dir / "adult_group_by_jia_fetal_marker_pair_fraction_rows_heatmap.png"
fig.savefig(group_fetal_fraction_heatmap, dpi=220, bbox_inches="tight")
plt.close(fig)

# Plot 4: parent cluster × fetal pair heatmap
cluster_heat_rows = []
for _, r in cluster_fetal_summary.iterrows():
    row = {
        "candidate_jia_group": r["candidate_jia_group"],
        "source_cluster": r["source_cluster"],
        "row_label": f"{r['candidate_jia_group']} | {r['source_cluster']}",
    }
    for fetal_label in fetal_signatures:
        row[fetal_label] = r[f"{fetal_label}__cells_weighted_pair_score"]
    cluster_heat_rows.append(row)
cluster_heat = pd.DataFrame(cluster_heat_rows)
cluster_heat = cluster_heat.sort_values(["candidate_jia_group", "source_cluster"])
cluster_mat = cluster_heat.set_index("row_label")[list(fetal_signatures.keys())]

fig, ax = plt.subplots(figsize=(10, max(8, 0.25 * len(cluster_mat))))
im = ax.imshow(cluster_mat.values, aspect="auto", interpolation="nearest")
ax.set_xticks(np.arange(cluster_mat.shape[1]))
ax.set_xticklabels(cluster_mat.columns, rotation=45, ha="right", fontsize=8)
ax.set_yticks(np.arange(cluster_mat.shape[0]))
ax.set_yticklabels(cluster_mat.index, fontsize=6)
ax.set_title("Parent cluster × Jia fetal marker-pair evidence\nValue = cell-weighted pair score")
cbar = fig.colorbar(im, ax=ax, shrink=0.8)
cbar.set_label("Cell-weighted pair score")
cluster_fetal_heatmap = plots_dir / "parent_cluster_by_jia_fetal_marker_pair_heatmap.png"
fig.savefig(cluster_fetal_heatmap, dpi=220, bbox_inches="tight")
plt.close(fig)

# Plot 5: Fetal-pair evidence river
# left fetal pair -> right adult group, weighted by evidence cell mass.
river_edges = group_fetal_long[group_fetal_long["evidence_cell_mass"] > 0].copy()
fetal_river = river_edges[["fetal_label", "candidate_jia_group", "evidence_cell_mass"]].copy()
fetal_river_png = plots_dir / "jia_fetal_marker_pair_to_adult_group_evidence_riverplot.png"
draw_river(
    fetal_river,
    levels=["fetal_label", "candidate_jia_group"],
    value_col="evidence_cell_mass",
    title=(
        "Jia fetal marker-pair → adult candidate group evidence river\n"
        "Edge width = sum(Number of cells × fetal pair score); marker overlap only, not label transfer"
    ),
    outpath=fetal_river_png,
    figsize=(18, 10),
    min_flow_to_label=0,
)

# Plot 6: individual fetal genes occurrence bar
gene_occ_plot = gene_occurrence.sort_values("selected_cells_with_gene", ascending=True)
fig, ax = plt.subplots(figsize=(10, 6))
y = np.arange(len(gene_occ_plot))
ax.barh(y, gene_occ_plot["selected_cells_with_gene"].values)
ax.set_yticks(y)
ax.set_yticklabels(gene_occ_plot["fetal_label"] + " | " + gene_occ_plot["gene"], fontsize=8)
ax.set_xlabel("Selected candidate Jia-style cells where gene appears in Top enriched genes")
ax.set_title("Individual Jia fetal-label genes found in selected Siletti workbook rows")
for yi, val, rows in zip(y, gene_occ_plot["selected_cells_with_gene"].values, gene_occ_plot["selected_rows_with_gene"].values):
    ax.text(val, yi, f" {int(val):,} cells; {int(rows)} rows", va="center", fontsize=7)
fetal_gene_bar = plots_dir / "individual_jia_fetal_label_gene_occurrence_barplot.png"
fig.savefig(fetal_gene_bar, dpi=220, bbox_inches="tight")
plt.close(fig)

# Plot 7: best fetal pair by adult group
best_plot = best_by_group.sort_values("best_cells_weighted_pair_score", ascending=True)
fig, ax = plt.subplots(figsize=(10, 6))
y = np.arange(len(best_plot))
ax.barh(y, best_plot["best_cells_weighted_pair_score"].values)
ax.set_yticks(y)
ax.set_yticklabels(best_plot["candidate_jia_group"], fontsize=8)
ax.set_xlabel("Best cell-weighted fetal pair score")
ax.set_title("Best Jia fetal marker-pair match for each adult candidate group")
for yi, score, label in zip(y, best_plot["best_cells_weighted_pair_score"].values, best_plot["best_fetal_marker_pair"].values):
    ax.text(score, yi, f" {label}; score={score:.2f}", va="center", fontsize=7)
best_fetal_bar = plots_dir / "best_fetal_marker_pair_by_adult_group_barplot.png"
fig.savefig(best_fetal_bar, dpi=220, bbox_inches="tight")
plt.close(fig)

# ============================================================
# README / interpretation
# ============================================================

selected_cells = int(selected[cells_col].sum())
excluded_cells = int(work.loc[~work["is_candidate_jia_group"], cells_col].sum())
total_cells = int(work[cells_col].sum())
assert selected_cells + excluded_cells == total_cells

readme = reports_dir / "README_fetal_marker_pair_validation.md"
readme.write_text(textwrap.dedent(f"""
# Jia fetal marker-pair validation of Siletti adult 9-group model

This folder adds a fetal-marker validation layer to the previous Jia-style
adult 9-group workbook model.

## Workflow provenance

Workflow script:

`{Path(__file__).resolve() if "__file__" in globals() else "unknown"}`

Python executable:

`{sys.executable}`

Output directory:

`{outdir}`

Previous Jia v2 outputs copied for provenance:

`{previous_folder if previous_folder else "not supplied"}`

## Input

`{xlsx_path}`

Rows loaded:

- all workbook rows: {len(work):,}
- unique subclusters: {work[subcluster_col].nunique():,}
- unique parent clusters: {work[cluster_col].nunique():,}
- cell-count sum: {total_cells:,}

Candidate Jia-style adult groups selected:

- rows: {len(selected):,}
- unique subclusters: {selected[subcluster_col].nunique():,}
- unique parent clusters: {selected[cluster_col].nunique():,}
- cell-count sum: {selected_cells:,}

Excluded:

- rows: {(~work["is_candidate_jia_group"]).sum():,}
- cell-count sum: {excluded_cells:,}

Total check:

{selected_cells:,} + {excluded_cells:,} = {total_cells:,}

## Anatomical assumption

Hippocampus is grouped with pallium/cortex.

- CerebralCortex + Hippocampus -> Pallial/cortical
- Amygdala / basal ganglia / basal forebrain / striatal / pallidal terms -> Subpallial
- everything else -> Other/uncertain

## Derived adult candidate labels

The Jia-style adult labels are operational group assignments made from workbook
fields. In particular, `Cortical PV+ basket neurons` is assigned when a row is
`Pallial/cortical`, has transferred MTG label `Pvalb`, and belongs to the
`MGE interneuron` supercluster. The word `basket` is not a literal label in the
Siletti/Linnarsson workbook.

## Jia fetal marker pairs tested

The five fetal labels in the Jia figure are represented as two-gene signatures:

1. EPHA5/MEF2C fetal neurons
2. LHX6/NFIA fetal neurons
3. CRABP1/ANGPT2 fetal neurons
4. NR2F1/NR2F2 fetal neurons
5. LHX8/ISL1 fetal neurons

## Score definition

For each row and each fetal pair:

`fetal pair score = number of pair genes found in Top enriched genes / 2`

So:

- 0.0 = neither gene found
- 0.5 = one gene found
- 1.0 = both genes found

For group summaries:

`cell-weighted pair score = sum(Number of cells × fetal pair score) / sum(Number of cells)`

For the river plot:

`evidence cell mass = sum(Number of cells × fetal pair score)`

## Important interpretation warning

This is marker-pair overlap, not label transfer.

The fetal marker pairs were used by Jia to name fetal populations. Seeing those
genes among Siletti adult `Top enriched genes` can support a thematic connection,
but absence does not rule out a lineage relationship because adult endpoints may
not retain fetal marker genes.
""").strip() + "\n")

interpretation = reports_dir / "INTERPRETATION_fetal_marker_pair_validation.md"
# Identify best pair text
best_lines = []
for _, r in best_by_group.iterrows():
    best_lines.append(
        f"- {r['candidate_jia_group']}: best = {r['best_fetal_marker_pair']} "
        f"(cell-weighted score {r['best_cells_weighted_pair_score']:.3f})"
    )

gene_lines = []
for _, r in gene_occurrence.iterrows():
    gene_lines.append(
        f"- {r['fetal_label']} / {r['gene']}: "
        f"selected rows {int(r['selected_rows_with_gene'])}, "
        f"selected cells {int(r['selected_cells_with_gene'])}"
    )

interpretation.write_text(textwrap.dedent("""
# Interpretation: fetal marker-pair validation

This analysis asks whether the Siletti adult candidate groups retain any of the
two-gene marker themes used to name the five Jia fetal populations.

## Stronger evidence means

A candidate adult group has repeated occurrence of one or both fetal marker-pair
genes in its `Top enriched genes`.

## Weak or absent evidence means

The adult workbook top markers do not retain that fetal naming gene pair.
This does not disprove a developmental relationship, because fetal markers can
be downregulated in adult terminal states.

## Operational label caveat

`Cortical PV+ basket neurons` is not a direct Siletti/Linnarsson workbook label.
It is the operational adult-group name used here for rows that are
Pallial/cortical, transferred MTG label `Pvalb`, and `MGE interneuron`
supercluster. A case-insensitive search of the workbook metadata found no
literal `basket` label.

## Best fetal marker-pair match by adult candidate group

""").strip() + "\n\n" + "\n".join(best_lines) + "\n\n## Individual fetal-label gene occurrence in selected adult groups\n\n" + "\n".join(gene_lines) + "\n")

# Save exact executed code
(source_dir / "siletti_jia9_fetal_marker_pair_validation.py").write_text(Path(__file__).read_text() if "__file__" in globals() and Path(__file__).exists() else "# Code was generated from ChatGPT execution cell. See zipped folder.\n")

zip_path = shutil.make_archive(str(outdir), "zip", root_dir=outdir) if args.make_zip else None

print("Created fetal marker-pair validation outputs.")
print("Output folder:", outdir)
print("ZIP:", zip_path if zip_path else "not requested")
print("\nStage QC:")
print(stage_qc.to_string(index=False))
print("\nFetal marker gene occurrence:")
print(gene_occurrence.to_string(index=False))
print("\nBest fetal marker-pair by adult candidate group:")
print(best_by_group.to_string(index=False))
print("\nKey files:")
for p in [
    qc_plot,
    group_fetal_heatmap,
    group_fetal_fraction_heatmap,
    cluster_fetal_heatmap,
    fetal_river_png,
    fetal_gene_bar,
    best_fetal_bar,
    tables_dir / "jia_fetal_marker_pair_definitions.csv",
    tables_dir / "jia_fetal_pair_gene_occurrence_in_workbook.csv",
    tables_dir / "adult_candidate_group_by_fetal_marker_pair_summary.csv",
    tables_dir / "adult_candidate_group_by_fetal_marker_pair_long.csv",
    tables_dir / "parent_cluster_by_fetal_marker_pair_summary.csv",
    tables_dir / "best_fetal_marker_pair_by_adult_candidate_group.csv",
    tables_dir / "row_level_fetal_marker_pair_scores_all_rows.csv",
    tables_dir / "row_level_fetal_marker_pair_scores_selected_only.csv",
    readme,
    interpretation,
]:
    print(p)
if zip_path:
    print(Path(zip_path))
