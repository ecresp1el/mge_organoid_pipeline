"""Gene-program scoring helpers for notebook-driven analyses."""

from importlib.metadata import version as package_version
from pathlib import Path
import re

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd
import scanpy as sc
from scanpy.tools._score_genes import _check_score_genes_args, _nan_means, _score_genes_bins


def parse_csv_list(raw):
    """Return a comma/colon/semicolon/space separated string as a clean list."""
    if raw is None:
        return []
    if isinstance(raw, (list, tuple)):
        values = raw
    else:
        values = re.split(r"[,;:\s]+", str(raw).strip())
    return [str(value).strip() for value in values if str(value).strip()]


def parse_optional_int(raw):
    """Parse an optional integer environment/config value."""
    if raw is None or str(raw).strip() == "":
        return None
    return int(raw)


def parse_optional_float(raw):
    """Parse an optional float environment/config value."""
    if raw is None or str(raw).strip() == "":
        return None
    return float(raw)


def safe_token(value):
    """Return a stable filename/column token."""
    token = re.sub(r"[^A-Za-z0-9]+", "_", str(value).strip())
    token = re.sub(r"_+", "_", token).strip("_")
    return token or "value"


def natural_sort_key(value):
    """Sort strings with embedded numbers in human order."""
    parts = re.split(r"(\d+)", str(value))
    return [int(part) if part.isdigit() else part for part in parts]


def read_marker_program_csv(path, gene_col="gene", program_col="cluster"):
    """Read a marker table and normalize the gene/program columns."""
    marker_path = Path(path).expanduser()
    if not marker_path.exists():
        raise FileNotFoundError("Missing marker program CSV: {}".format(marker_path))
    markers = pd.read_csv(marker_path, encoding="utf-8-sig")
    missing = [col for col in [gene_col, program_col] if col not in markers.columns]
    if missing:
        raise ValueError("Missing required marker columns: {}".format(", ".join(missing)))
    markers = markers.copy()
    markers[gene_col] = markers[gene_col].astype(str).str.strip()
    markers[program_col] = markers[program_col].astype(str).str.strip()
    markers = markers[(markers[gene_col] != "") & (markers[program_col] != "")].copy()
    markers.insert(0, "source_marker_row", np.arange(1, len(markers) + 1))
    markers["source_program_gene_rank"] = markers.groupby(program_col, sort=False).cumcount() + 1
    return markers


def select_program_markers(
    markers,
    gene_col="gene",
    program_col="cluster",
    top_n=None,
    min_avg_log2fc=None,
    max_p_val_adj=None,
    sort_by="avg_log2FC",
):
    """Filter and optionally take the top N markers per program."""
    selected = markers.copy()
    if min_avg_log2fc is not None and "avg_log2FC" in selected.columns:
        selected = selected[selected["avg_log2FC"] >= float(min_avg_log2fc)].copy()
    if max_p_val_adj is not None and "p_val_adj" in selected.columns:
        selected = selected[selected["p_val_adj"] <= float(max_p_val_adj)].copy()
    if sort_by and sort_by in selected.columns:
        ascending = sort_by.lower() in {"p_val", "p_val_adj", "padj", "qval", "q_value"}
        selected = selected.sort_values([program_col, sort_by], ascending=[True, ascending])
    if top_n is not None:
        selected = selected.groupby(program_col, sort=False).head(int(top_n)).copy()
    selected["selected_program_gene_rank"] = selected.groupby(program_col, sort=False).cumcount() + 1
    return selected


def programs_from_marker_table(
    markers,
    gene_col="gene",
    program_col="cluster",
    program_order=None,
):
    """Return an ordered dict-like mapping from program name to unique genes."""
    programs_present = list(pd.unique(markers[program_col]))
    ordered = []
    for program in program_order or []:
        if program in programs_present and program not in ordered:
            ordered.append(program)
    for program in programs_present:
        if program not in ordered:
            ordered.append(program)

    programs = {}
    for program in ordered:
        values = markers.loc[markers[program_col] == program, gene_col].astype(str)
        genes = []
        seen = set()
        for gene in values:
            clean = gene.strip()
            if clean and clean not in seen:
                genes.append(clean)
                seen.add(clean)
        programs[program] = genes
    return programs


def match_program_genes(programs, var_names):
    """Match requested program genes to AnnData var_names with exact/case-insensitive fallback."""
    var_values = [str(value) for value in var_names]
    var_set = set(var_values)
    upper_to_var = {}
    ambiguous_upper = set()
    for value in var_values:
        key = value.upper()
        if key in upper_to_var and upper_to_var[key] != value:
            ambiguous_upper.add(key)
        else:
            upper_to_var[key] = value

    detail_records = []
    matched_programs = {}
    summary_records = []
    for program, genes in programs.items():
        used = []
        used_set = set()
        unique_input = []
        seen_input = set()
        for gene in genes:
            gene = str(gene).strip()
            if not gene or gene in seen_input:
                continue
            unique_input.append(gene)
            seen_input.add(gene)

        for gene in unique_input:
            matched_gene = ""
            match_type = "missing"
            key = gene.upper()
            if gene in var_set:
                matched_gene = gene
                match_type = "exact"
            elif key in ambiguous_upper:
                match_type = "ambiguous_case_insensitive"
            elif key in upper_to_var:
                matched_gene = upper_to_var[key]
                match_type = "case_insensitive"

            used_for_scoring = bool(matched_gene) and matched_gene not in used_set
            if used_for_scoring:
                used.append(matched_gene)
                used_set.add(matched_gene)
            detail_records.append(
                {
                    "program": program,
                    "input_gene": gene,
                    "matched_gene": matched_gene,
                    "matched": bool(matched_gene),
                    "match_type": match_type,
                    "used_for_scoring": used_for_scoring,
                }
            )

        matched_programs[program] = used
        n_input = len(unique_input)
        n_matched = sum(
            1
            for record in detail_records
            if record["program"] == program and record["matched"]
        )
        summary_records.append(
            {
                "program": program,
                "n_input_genes": n_input,
                "n_matched_genes": n_matched,
                "n_used_genes": len(used),
                "n_missing_genes": n_input - n_matched,
                "pct_input_genes_used": 100.0 * len(used) / n_input if n_input else 0.0,
            }
        )

    return matched_programs, pd.DataFrame(summary_records), pd.DataFrame(detail_records)


def score_programs_scanpy(
    adata,
    matched_programs,
    score_prefix="jia_score_",
    ctrl_size=50,
    random_state=0,
):
    """Score each gene program with scanpy.tl.score_genes on adata.X."""
    score_columns = {}
    for program, genes in matched_programs.items():
        if not genes:
            continue
        score_col = "{}{}".format(score_prefix, safe_token(program))
        print(
            "[GeneProgramScoring] scoring {} with {} genes -> {}".format(
                program,
                len(genes),
                score_col,
            ),
            flush=True,
        )
        sc.tl.score_genes(
            adata,
            gene_list=list(genes),
            score_name=score_col,
            ctrl_size=int(ctrl_size),
            random_state=int(random_state),
            use_raw=False,
            copy=False,
        )
        score_columns[program] = score_col
    if not score_columns:
        raise ValueError("No programs had matched genes for scoring.")
    return score_columns


def scanpy_score_genes_control_audit(
    adata,
    matched_programs,
    score_columns=None,
    ctrl_size=50,
    random_state=0,
    n_bins=25,
    ctrl_as_ref=True,
    gene_pool=None,
    use_raw=False,
    layer=None,
):
    """Return the exact Scanpy score_genes control genes for each program.

    This intentionally calls Scanpy's private helper functions to mirror the
    installed Scanpy version used by `scanpy.tl.score_genes`.
    """
    scanpy_version = package_version("scanpy")
    score_columns = score_columns or {}
    control_records = []
    program_bin_records = []
    summary_records = []

    for program, genes in matched_programs.items():
        if not genes:
            continue
        if random_state is not None:
            np.random.seed(random_state)
        gene_list, resolved_gene_pool, get_subset = _check_score_genes_args(
            adata,
            genes,
            gene_pool,
            use_raw=bool(use_raw),
            layer=layer,
        )

        obs_avg = pd.Series(
            _nan_means(get_subset(resolved_gene_pool), axis=0),
            index=resolved_gene_pool,
        )
        obs_avg = obs_avg[np.isfinite(obs_avg)]
        n_items = int(np.round(len(obs_avg) / (int(n_bins) - 1)))
        obs_cut = obs_avg.rank(method="min") // n_items

        for gene in gene_list:
            expression_bin = obs_cut.get(gene, np.nan)
            mean_expression = obs_avg.get(gene, np.nan)
            program_bin_records.append(
                {
                    "program": program,
                    "score_col": score_columns.get(program, ""),
                    "program_gene": gene,
                    "mean_expression": float(mean_expression) if pd.notna(mean_expression) else np.nan,
                    "expression_bin": int(expression_bin) if pd.notna(expression_bin) else np.nan,
                    "scanpy_version": scanpy_version,
                    "ctrl_size": int(ctrl_size),
                    "n_bins": int(n_bins),
                    "ctrl_as_ref": bool(ctrl_as_ref),
                    "random_state": random_state,
                }
            )

        control_genes = pd.Index([], dtype="string")
        for r_genes in _score_genes_bins(
            gene_list,
            resolved_gene_pool,
            ctrl_as_ref=bool(ctrl_as_ref),
            ctrl_size=int(ctrl_size),
            n_bins=int(n_bins),
            get_subset=get_subset,
        ):
            for gene in r_genes:
                expression_bin = obs_cut.get(gene, np.nan)
                mean_expression = obs_avg.get(gene, np.nan)
                control_records.append(
                    {
                        "program": program,
                        "score_col": score_columns.get(program, ""),
                        "control_gene": gene,
                        "mean_expression": (
                            float(mean_expression) if pd.notna(mean_expression) else np.nan
                        ),
                        "expression_bin": (
                            int(expression_bin) if pd.notna(expression_bin) else np.nan
                        ),
                        "scanpy_version": scanpy_version,
                        "ctrl_size": int(ctrl_size),
                        "n_bins": int(n_bins),
                        "ctrl_as_ref": bool(ctrl_as_ref),
                        "random_state": random_state,
                    }
                )
            control_genes = control_genes.union(r_genes)

        program_bins = sorted(
            [int(value) for value in pd.unique(obs_cut.loc[gene_list]) if pd.notna(value)]
        )
        summary_records.append(
            {
                "program": program,
                "score_col": score_columns.get(program, ""),
                "scanpy_version": scanpy_version,
                "n_program_genes": int(len(gene_list)),
                "n_gene_pool_genes": int(len(resolved_gene_pool)),
                "n_expression_bins": int(n_bins),
                "n_program_expression_bins": int(len(program_bins)),
                "program_expression_bins": ",".join(str(value) for value in program_bins),
                "ctrl_size": int(ctrl_size),
                "ctrl_as_ref": bool(ctrl_as_ref),
                "random_state": random_state,
                "n_control_genes": int(len(control_genes)),
            }
        )

    return (
        pd.DataFrame(summary_records),
        pd.DataFrame(control_records),
        pd.DataFrame(program_bin_records),
    )


def choose_first_existing(columns, candidates):
    """Return the first candidate column present in a pandas Index/list."""
    column_set = set(columns)
    for candidate in candidates:
        if candidate in column_set:
            return candidate
    return None


def choose_umap_key(adata, candidates=("X_umap_seurat", "X_umap")):
    """Return the first available UMAP key in adata.obsm."""
    keys = set(adata.obsm.keys())
    for candidate in candidates:
        if candidate in keys:
            return candidate
    raise KeyError("No UMAP key found. Checked: {}".format(", ".join(candidates)))


def attach_resolution_assignments(adata, assignments_path, id_col="cell_id"):
    """Attach precomputed resolution-sweep cluster columns to adata.obs."""
    path = Path(assignments_path).expanduser()
    if not path.exists():
        print("[GeneProgramScoring] resolution assignment table not found: {}".format(path), flush=True)
        return []

    assignments = pd.read_csv(path, sep="\t", dtype=str)
    if id_col not in assignments.columns:
        raise ValueError("Resolution assignment table is missing {}: {}".format(id_col, path))
    resolution_cols = [col for col in assignments.columns if col != id_col]
    assignments = assignments.set_index(id_col)

    if id_col in adata.obs.columns:
        obs_key = adata.obs[id_col].astype(str)
        matched_by_obs_col = int(obs_key.isin(assignments.index).sum())
    else:
        obs_key = pd.Series(adata.obs_names.astype(str), index=adata.obs_names)
        matched_by_obs_col = 0

    index_key = pd.Series(adata.obs_names.astype(str), index=adata.obs_names)
    matched_by_index = int(index_key.isin(assignments.index).sum())
    align_key = obs_key if matched_by_obs_col >= matched_by_index else index_key

    aligned = assignments.reindex(align_key.values)
    aligned.index = adata.obs_names
    for col in resolution_cols:
        adata.obs[col] = aligned[col].astype("category")

    print(
        "[GeneProgramScoring] attached {} resolution columns from {}".format(
            len(resolution_cols),
            path,
        ),
        flush=True,
    )
    return resolution_cols


def score_threshold_table(obs, score_columns, quantile=0.9):
    """Return global high-score thresholds for each program score."""
    records = []
    for program, score_col in score_columns.items():
        values = pd.to_numeric(obs[score_col], errors="coerce")
        threshold = float(values.quantile(float(quantile)))
        records.append(
            {
                "program": program,
                "score_col": score_col,
                "high_score_quantile": float(quantile),
                "high_score_threshold": threshold,
            }
        )
    return pd.DataFrame(records)


def gene_set_overlap(left_genes, right_genes, label, left_label="left", right_label="right"):
    """Return summary/detail tables for overlap between two gene sets."""
    left = pd.Index([str(gene) for gene in left_genes]).drop_duplicates()
    right = pd.Index([str(gene) for gene in right_genes]).drop_duplicates()
    overlap = left.intersection(right)
    union = left.union(right)
    summary = pd.DataFrame(
        [
            {
                "comparison": label,
                "left_label": left_label,
                "right_label": right_label,
                "n_left_genes": int(len(left)),
                "n_right_genes": int(len(right)),
                "n_overlap_genes": int(len(overlap)),
                "n_union_genes": int(len(union)),
                "pct_left_overlapping": 100.0 * len(overlap) / len(left) if len(left) else np.nan,
                "pct_right_overlapping": 100.0 * len(overlap) / len(right) if len(right) else np.nan,
                "jaccard_index": len(overlap) / len(union) if len(union) else np.nan,
            }
        ]
    )
    detail = pd.DataFrame(
        [
            {
                "comparison": label,
                "overlap_gene": gene,
                "left_label": left_label,
                "right_label": right_label,
            }
            for gene in sorted(overlap, key=natural_sort_key)
        ]
    )
    return summary, detail


def score_pair_correlations(obs, pairs):
    """Compute Pearson and Spearman correlations for named score-column pairs."""
    records = []
    for pair in pairs:
        left_col = pair["left_col"]
        right_col = pair["right_col"]
        frame = obs[[left_col, right_col]].apply(pd.to_numeric, errors="coerce").dropna()
        records.append(
            {
                "comparison": pair.get("comparison", ""),
                "comparison_display": pair.get("comparison_display", pair.get("comparison", "")),
                "top_n": pair.get("top_n", ""),
                "left_label": pair.get("left_label", left_col),
                "right_label": pair.get("right_label", right_col),
                "left_col": left_col,
                "right_col": right_col,
                "n_cells": int(len(frame)),
                "pearson_r": (
                    float(frame[left_col].corr(frame[right_col], method="pearson"))
                    if len(frame) > 1
                    else np.nan
                ),
                "spearman_r": (
                    float(frame[left_col].corr(frame[right_col], method="spearman"))
                    if len(frame) > 1
                    else np.nan
                ),
            }
        )
    return pd.DataFrame(records)


def add_contrast_score(adata, left_col, right_col, contrast_col):
    """Add left minus right contrast score to adata.obs."""
    adata.obs[contrast_col] = (
        pd.to_numeric(adata.obs[left_col], errors="coerce")
        - pd.to_numeric(adata.obs[right_col], errors="coerce")
    )
    return contrast_col


def summarize_columns_by_group(obs, score_columns, group_col):
    """Summarize arbitrary score columns by group."""
    if group_col not in obs.columns:
        raise KeyError("Missing group column: {}".format(group_col))
    frame = obs[[group_col] + list(score_columns.values())].copy()
    frame[group_col] = frame[group_col].astype(str)
    frame = frame[frame[group_col].notna() & (frame[group_col] != "nan")]
    records = []
    groups = sorted(pd.unique(frame[group_col]), key=natural_sort_key)
    for group in groups:
        sub = frame[frame[group_col] == group]
        for label, score_col in score_columns.items():
            values = pd.to_numeric(sub[score_col], errors="coerce").dropna()
            records.append(
                {
                    "groupby": group_col,
                    "group": group,
                    "score_label": label,
                    "score_col": score_col,
                    "n_cells": int(len(values)),
                    "mean_score": float(values.mean()) if len(values) else np.nan,
                    "median_score": float(values.median()) if len(values) else np.nan,
                    "std_score": float(values.std(ddof=0)) if len(values) else np.nan,
                    "q25_score": float(values.quantile(0.25)) if len(values) else np.nan,
                    "q75_score": float(values.quantile(0.75)) if len(values) else np.nan,
                    "fraction_positive": (
                        float((values > 0).mean()) if len(values) else np.nan
                    ),
                }
            )
    return pd.DataFrame(records)


def summarize_scores_by_group(obs, score_columns, group_col, thresholds=None):
    """Summarize program scores by one categorical grouping column."""
    if group_col not in obs.columns:
        raise KeyError("Missing group column: {}".format(group_col))
    threshold_map = {}
    if thresholds is not None and len(thresholds) > 0:
        threshold_map = dict(zip(thresholds["score_col"], thresholds["high_score_threshold"]))

    frame = obs[[group_col] + list(score_columns.values())].copy()
    frame[group_col] = frame[group_col].astype(str)
    frame = frame[frame[group_col].notna() & (frame[group_col] != "nan")]

    records = []
    groups = sorted(pd.unique(frame[group_col]), key=natural_sort_key)
    for group in groups:
        sub = frame[frame[group_col] == group]
        for program, score_col in score_columns.items():
            values = pd.to_numeric(sub[score_col], errors="coerce").dropna()
            threshold = threshold_map.get(score_col, np.nan)
            records.append(
                {
                    "groupby": group_col,
                    "group": group,
                    "program": program,
                    "score_col": score_col,
                    "n_cells": int(len(values)),
                    "mean_score": float(values.mean()) if len(values) else np.nan,
                    "median_score": float(values.median()) if len(values) else np.nan,
                    "std_score": float(values.std(ddof=0)) if len(values) else np.nan,
                    "q25_score": float(values.quantile(0.25)) if len(values) else np.nan,
                    "q75_score": float(values.quantile(0.75)) if len(values) else np.nan,
                    "high_score_threshold": threshold,
                    "fraction_high_score": (
                        float((values >= threshold).mean())
                        if len(values) and not pd.isna(threshold)
                        else np.nan
                    ),
                }
            )
    return pd.DataFrame(records)


def best_matches(summary, top_n=3):
    """Return top scoring groups per program and grouping column."""
    if summary.empty:
        return summary.copy()
    ranked = summary.copy()
    ranked["mean_score_rank"] = ranked.groupby(["groupby", "program"])["mean_score"].rank(
        method="first",
        ascending=False,
    )
    ranked["fraction_high_score_rank"] = ranked.groupby(["groupby", "program"])[
        "fraction_high_score"
    ].rank(method="first", ascending=False)
    ranked = ranked[
        (ranked["mean_score_rank"] <= int(top_n))
        | (ranked["fraction_high_score_rank"] <= int(top_n))
    ].copy()
    return ranked.sort_values(["groupby", "program", "mean_score_rank", "fraction_high_score_rank"])


def score_output_obs_table(adata, score_columns, extra_columns=None, id_col="cell_id"):
    """Build a compact per-cell score table from adata.obs."""
    columns = []
    if id_col in adata.obs.columns:
        columns.append(id_col)
    for col in extra_columns or []:
        if col in adata.obs.columns and col not in columns:
            columns.append(col)
    for col in score_columns.values():
        if col not in columns:
            columns.append(col)
    table = adata.obs[columns].copy()
    table.insert(0, "obs_name", adata.obs_names.astype(str))
    return table


def dense_vector(values):
    """Return a one-dimensional dense numpy vector from dense/sparse matrix slices."""
    if hasattr(values, "toarray"):
        values = values.toarray()
    array = np.asarray(values)
    return array.reshape(-1)


def summarize_marker_expression_by_group(
    adata,
    matched_marker_groups,
    group_col,
    expression_cutoff=0.0,
):
    """Summarize target marker expression by one categorical grouping column."""
    if group_col not in adata.obs.columns:
        raise KeyError("Missing group column: {}".format(group_col))

    groups = adata.obs[group_col].astype(str)
    valid_groups = sorted(pd.unique(groups[groups != "nan"]), key=natural_sort_key)
    records = []
    for marker_group, genes in matched_marker_groups.items():
        for gene in genes:
            if gene not in adata.var_names:
                continue
            values = dense_vector(adata[:, [gene]].X)
            for group in valid_groups:
                mask = groups.values == group
                group_values = values[mask]
                detected = group_values > float(expression_cutoff)
                records.append(
                    {
                        "groupby": group_col,
                        "group": group,
                        "marker_group": marker_group,
                        "marker_gene": gene,
                        "n_cells": int(mask.sum()),
                        "mean_expression": (
                            float(np.mean(group_values)) if len(group_values) else np.nan
                        ),
                        "median_expression": (
                            float(np.median(group_values)) if len(group_values) else np.nan
                        ),
                        "fraction_expressing": (
                            float(np.mean(detected)) if len(group_values) else np.nan
                        ),
                        "pct_expressing": (
                            100.0 * float(np.mean(detected)) if len(group_values) else np.nan
                        ),
                        "expression_cutoff": float(expression_cutoff),
                    }
                )
    return pd.DataFrame(records)


def marker_group_support_table(marker_summary):
    """Collapse per-gene marker summaries into marker-group support by cluster."""
    if marker_summary.empty:
        return marker_summary.copy()
    grouped = (
        marker_summary.groupby(["groupby", "group", "marker_group"], sort=False)
        .agg(
            n_marker_genes=("marker_gene", "nunique"),
            mean_marker_expression=("mean_expression", "mean"),
            mean_fraction_expressing=("fraction_expressing", "mean"),
            max_marker_expression=("mean_expression", "max"),
            max_fraction_expressing=("fraction_expressing", "max"),
        )
        .reset_index()
    )
    grouped["_group_sort"] = grouped["group"].map(lambda value: tuple(natural_sort_key(value)))
    grouped = grouped.sort_values(["groupby", "marker_group", "_group_sort"]).drop(
        columns=["_group_sort"]
    )
    return grouped


def cluster_marker_interpretation_table(base_summary, marker_support):
    """Create a compact score-plus-marker support table by cluster."""
    if base_summary.empty:
        return pd.DataFrame()

    score = base_summary.copy()
    score["mean_score_rank"] = score.groupby("program")["mean_score"].rank(
        method="first",
        ascending=False,
    )
    score["fraction_high_score_rank"] = score.groupby("program")["fraction_high_score"].rank(
        method="first",
        ascending=False,
    )
    score_records = []
    for group, sub in score.groupby("group", sort=False):
        top = sub.sort_values("mean_score", ascending=False).iloc[0]
        record = {
            "groupby": top["groupby"],
            "group": str(group),
            "n_cells": int(top["n_cells"]),
            "top_jia_program_by_mean_score": top["program"],
            "top_jia_program_mean_score": float(top["mean_score"]),
            "top_jia_program_fraction_high_score": float(top["fraction_high_score"]),
        }
        for _, row in sub.iterrows():
            token = safe_token(row["program"]).lower()
            record["{}_mean_score".format(token)] = float(row["mean_score"])
            record["{}_mean_score_rank".format(token)] = float(row["mean_score_rank"])
            record["{}_fraction_high_score".format(token)] = float(row["fraction_high_score"])
            record["{}_fraction_high_score_rank".format(token)] = float(row["fraction_high_score_rank"])
        score_records.append(record)
    result = pd.DataFrame(score_records)

    if marker_support.empty:
        result["_group_sort"] = result["group"].map(lambda value: tuple(natural_sort_key(value)))
        return result.sort_values("_group_sort").drop(columns=["_group_sort"])

    marker = marker_support.copy()
    marker_records = []
    for group, sub in marker.groupby("group", sort=False):
        record = {"group": str(group)}
        for _, row in sub.iterrows():
            token = safe_token(row["marker_group"]).lower()
            record["{}_n_marker_genes".format(token)] = int(row["n_marker_genes"])
            record["{}_mean_marker_expression".format(token)] = float(row["mean_marker_expression"])
            record["{}_mean_fraction_expressing".format(token)] = float(
                row["mean_fraction_expressing"]
            )
        marker_records.append(record)
    marker_table = pd.DataFrame(marker_records)
    result = result.merge(marker_table, on="group", how="left")
    result["_group_sort"] = result["group"].map(lambda value: tuple(natural_sort_key(value)))
    return result.sort_values("_group_sort").drop(columns=["_group_sort"])


def _plot_base_categorical(ax, x, y, values, title, point_size):
    values = pd.Series(values).astype(str)
    levels = sorted(pd.unique(values), key=natural_sort_key)
    if len(levels) <= 20:
        cmap = plt.get_cmap("tab20")
        colors = [cmap(i % 20) for i in range(len(levels))]
    else:
        cmap = plt.get_cmap("hsv")
        colors = [cmap(i / max(len(levels), 1)) for i in range(len(levels))]
    for level, color in zip(levels, colors):
        mask = values.values == level
        ax.scatter(
            x[mask],
            y[mask],
            s=point_size,
            c=[color],
            linewidths=0,
            rasterized=True,
            label=str(level),
        )
    ax.set_title(title, fontsize=10)
    if len(levels) <= 15:
        ax.legend(
            title="cluster",
            loc="center left",
            bbox_to_anchor=(1.01, 0.5),
            frameon=False,
            fontsize=6,
            title_fontsize=7,
            markerscale=4,
        )


def _plot_score_overlay(
    ax,
    x,
    y,
    values,
    title,
    point_size,
    background_color,
    cmap="viridis",
    center_zero=False,
):
    scores = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy()
    ax.scatter(
        x,
        y,
        s=point_size,
        c=background_color,
        linewidths=0,
        rasterized=True,
    )
    mask = np.isfinite(scores)
    if not np.any(mask):
        ax.set_title(title, fontsize=10)
        return None
    order = np.argsort(scores[mask])
    x_plot = x[mask][order]
    y_plot = y[mask][order]
    score_plot = scores[mask][order]
    norm = None
    if center_zero:
        vmax_abs = float(np.nanquantile(np.abs(score_plot), 0.99))
        if np.isclose(vmax_abs, 0):
            vmax_abs = float(np.nanmax(np.abs(score_plot)))
        if not np.isclose(vmax_abs, 0):
            norm = TwoSlopeNorm(vmin=-vmax_abs, vcenter=0.0, vmax=vmax_abs)
        vmin = vmax = None
    else:
        vmin, vmax = np.nanquantile(score_plot, [0.01, 0.99])
        if np.isclose(vmin, vmax):
            vmin = float(np.nanmin(score_plot))
            vmax = float(np.nanmax(score_plot))
    scatter = ax.scatter(
        x_plot,
        y_plot,
        s=point_size,
        c=score_plot,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        norm=norm,
        linewidths=0,
        alpha=0.9,
        rasterized=True,
    )
    ax.set_title(title, fontsize=10)
    return scatter


def plot_marker_umap_feature_grid(
    adata,
    matched_marker_groups,
    output_path,
    umap_key=None,
    point_size=None,
    background_color="#d4d4d4",
):
    """Save a marker feature UMAP grid grouped by marker class."""
    selected_umap_key = umap_key or choose_umap_key(adata)
    coords = np.asarray(adata.obsm[selected_umap_key])
    x = coords[:, 0]
    y = coords[:, 1]
    if point_size is None:
        point_size = 1.0 if adata.n_obs >= 50000 else 3.0

    marker_groups = [(group, genes) for group, genes in matched_marker_groups.items() if genes]
    if not marker_groups:
        raise ValueError("No matched marker genes to plot.")
    n_rows = len(marker_groups)
    n_cols = max(len(genes) for _, genes in marker_groups)
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(max(3.2 * n_cols, 6), max(3.0 * n_rows, 3.5)),
        squeeze=False,
        constrained_layout=True,
    )

    for row_idx, (marker_group, genes) in enumerate(marker_groups):
        for col_idx in range(n_cols):
            ax = axes[row_idx, col_idx]
            if col_idx >= len(genes):
                ax.axis("off")
                continue
            gene = genes[col_idx]
            values = dense_vector(adata[:, [gene]].X)
            scatter = _plot_score_overlay(
                ax,
                x,
                y,
                values,
                gene,
                point_size,
                background_color,
            )
            if col_idx == 0:
                ax.set_ylabel(marker_group, fontsize=9)
            if scatter is not None:
                cbar = fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.02)
                cbar.ax.tick_params(labelsize=6)
    output = Path(output_path).expanduser()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    return fig


def plot_umap_score_overlay_panel(
    adata,
    score_columns,
    output_path,
    base_color_col,
    umap_key=None,
    title_prefix="DIV30",
    point_size=None,
    background_color="#d4d4d4",
    score_cmap="viridis",
    center_zero_scores=False,
):
    """Save a 1 x N panel: base cluster UMAP, then score overlays."""
    if base_color_col not in adata.obs.columns:
        raise KeyError("Missing base UMAP color column: {}".format(base_color_col))
    selected_umap_key = umap_key or choose_umap_key(adata)
    coords = np.asarray(adata.obsm[selected_umap_key])
    if coords.shape[1] < 2:
        raise ValueError("UMAP coordinates must have at least two columns: {}".format(selected_umap_key))
    x = coords[:, 0]
    y = coords[:, 1]
    if point_size is None:
        point_size = 1.0 if adata.n_obs >= 50000 else 3.0

    n_panels = 1 + len(score_columns)
    fig, axes = plt.subplots(
        1,
        n_panels,
        figsize=(max(4.2 * n_panels, 8), 4.2),
        constrained_layout=True,
    )
    if n_panels == 1:
        axes = [axes]

    _plot_base_categorical(
        axes[0],
        x,
        y,
        adata.obs[base_color_col],
        "{} {}".format(title_prefix, base_color_col),
        point_size,
    )
    for ax, (program, score_col) in zip(axes[1:], score_columns.items()):
        scatter = _plot_score_overlay(
            ax,
            x,
            y,
            adata.obs[score_col],
            "{} {}".format(title_prefix, program),
            point_size,
            background_color,
            cmap=score_cmap,
            center_zero=center_zero_scores,
        )
        if scatter is not None:
            cbar = fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.02)
            cbar.ax.tick_params(labelsize=7)
            cbar.set_label(score_col, fontsize=7)

    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_aspect("equal", adjustable="box")
        for spine in ax.spines.values():
            spine.set_visible(False)

    output = Path(output_path).expanduser()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    return fig


def plot_score_heatmap(summary, output_path, value_col="mean_score", title=None):
    """Save a compact group x program score heatmap."""
    if summary.empty:
        raise ValueError("Cannot plot heatmap from an empty summary table.")
    matrix = summary.pivot(index="group", columns="program", values=value_col)
    matrix = matrix.reindex(sorted(matrix.index, key=natural_sort_key))
    fig_height = max(3.5, min(14, 0.28 * len(matrix.index) + 1.5))
    fig_width = max(5.0, 1.4 * len(matrix.columns) + 2.0)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), constrained_layout=True)
    image = ax.imshow(matrix.to_numpy(dtype=float), aspect="auto", cmap="viridis")
    ax.set_xticks(np.arange(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(np.arange(matrix.shape[0]))
    ax.set_yticklabels(matrix.index, fontsize=7)
    ax.set_xlabel("program")
    ax.set_ylabel(summary["groupby"].iloc[0])
    ax.set_title(title or value_col, fontsize=10)
    cbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.02)
    cbar.set_label(value_col, fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    output = Path(output_path).expanduser()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=200, bbox_inches="tight")
    return fig


def plot_program_proportion_dotplot(
    summary,
    output_path,
    proportion_col="fraction_high_score",
    color_col="mean_score",
    title=None,
    min_dot_size=20,
    max_dot_size=360,
):
    """Save a cluster x program dot plot.

    Dot area shows the proportion/fraction column. Dot color shows the
    continuous score summary column.
    """
    if summary.empty:
        raise ValueError("Cannot plot dotplot from an empty summary table.")
    required = {"group", "program", proportion_col, color_col}
    missing = sorted(required.difference(summary.columns))
    if missing:
        raise ValueError("Summary table is missing columns: {}".format(", ".join(missing)))

    plot_df = summary.copy()
    plot_df["group"] = plot_df["group"].astype(str)
    groups = sorted(pd.unique(plot_df["group"]), key=natural_sort_key)
    programs = list(pd.unique(plot_df["program"]))

    group_to_y = {group: i for i, group in enumerate(groups)}
    program_to_x = {program: i for i, program in enumerate(programs)}
    plot_df["x"] = plot_df["program"].map(program_to_x)
    plot_df["y"] = plot_df["group"].map(group_to_y)
    proportions = pd.to_numeric(plot_df[proportion_col], errors="coerce").fillna(0.0).clip(0, 1)
    colors = pd.to_numeric(plot_df[color_col], errors="coerce")
    sizes = min_dot_size + proportions * (max_dot_size - min_dot_size)

    fig_height = max(3.8, 0.38 * len(groups) + 1.6)
    fig_width = max(5.2, 1.25 * len(programs) + 2.8)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), constrained_layout=True)
    scatter = ax.scatter(
        plot_df["x"],
        plot_df["y"],
        s=sizes,
        c=colors,
        cmap="viridis",
        edgecolors="#333333",
        linewidths=0.25,
    )
    ax.set_xticks(range(len(programs)))
    ax.set_xticklabels(programs, fontsize=9)
    ax.set_yticks(range(len(groups)))
    ax.set_yticklabels(groups, fontsize=8)
    ax.set_xlabel("Jia program")
    ax.set_ylabel(summary["groupby"].iloc[0])
    ax.set_title(title or "{} by cluster".format(proportion_col), fontsize=10)
    ax.grid(axis="both", color="#e5e5e5", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.invert_yaxis()

    cbar = fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label(color_col, fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    legend_values = [0.05, 0.10, 0.25, 0.50]
    handles = []
    labels = []
    for value in legend_values:
        handles.append(
            ax.scatter(
                [],
                [],
                s=min_dot_size + value * (max_dot_size - min_dot_size),
                c="#bdbdbd",
                edgecolors="#333333",
                linewidths=0.25,
            )
        )
        labels.append("{:.0f}%".format(value * 100))
    ax.legend(
        handles,
        labels,
        title=proportion_col,
        loc="center left",
        bbox_to_anchor=(1.14, 0.5),
        frameon=False,
        fontsize=7,
        title_fontsize=8,
    )

    output = Path(output_path).expanduser()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    return fig


def plot_marker_expression_dotplot(
    marker_summary,
    output_path,
    marker_order=None,
    group_order=None,
    title=None,
    min_dot_size=20,
    max_dot_size=360,
):
    """Save a cluster x marker dot plot.

    Dot area shows percent expressing. Dot color shows mean expression.
    """
    if marker_summary.empty:
        raise ValueError("Cannot plot marker dotplot from an empty summary table.")
    plot_df = marker_summary.copy()
    plot_df["group"] = plot_df["group"].astype(str)
    if marker_order is None:
        marker_order = list(pd.unique(plot_df["marker_gene"]))
    if group_order is None:
        group_order = sorted(pd.unique(plot_df["group"]), key=natural_sort_key)

    group_to_y = {group: i for i, group in enumerate(group_order)}
    marker_to_x = {marker: i for i, marker in enumerate(marker_order)}
    plot_df = plot_df[
        plot_df["group"].isin(group_to_y) & plot_df["marker_gene"].isin(marker_to_x)
    ].copy()
    plot_df["x"] = plot_df["marker_gene"].map(marker_to_x)
    plot_df["y"] = plot_df["group"].map(group_to_y)
    fractions = pd.to_numeric(plot_df["fraction_expressing"], errors="coerce").fillna(0).clip(0, 1)
    colors = pd.to_numeric(plot_df["mean_expression"], errors="coerce")
    sizes = min_dot_size + fractions * (max_dot_size - min_dot_size)

    fig_height = max(3.8, 0.38 * len(group_order) + 1.6)
    fig_width = max(8.0, 0.75 * len(marker_order) + 3.5)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), constrained_layout=True)
    scatter = ax.scatter(
        plot_df["x"],
        plot_df["y"],
        s=sizes,
        c=colors,
        cmap="viridis",
        edgecolors="#333333",
        linewidths=0.25,
    )
    ax.set_xticks(range(len(marker_order)))
    ax.set_xticklabels(marker_order, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(group_order)))
    ax.set_yticklabels(group_order, fontsize=8)
    ax.set_xlabel("marker gene")
    ax.set_ylabel(marker_summary["groupby"].iloc[0])
    ax.set_title(title or "Marker expression by cluster", fontsize=10)
    ax.grid(axis="both", color="#e5e5e5", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.invert_yaxis()

    cbar = fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label("mean_expression", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    legend_values = [0.05, 0.10, 0.25, 0.50]
    handles = []
    labels = []
    for value in legend_values:
        handles.append(
            ax.scatter(
                [],
                [],
                s=min_dot_size + value * (max_dot_size - min_dot_size),
                c="#bdbdbd",
                edgecolors="#333333",
                linewidths=0.25,
            )
        )
        labels.append("{:.0f}%".format(value * 100))
    ax.legend(
        handles,
        labels,
        title="pct expressing",
        loc="center left",
        bbox_to_anchor=(1.14, 0.5),
        frameon=False,
        fontsize=7,
        title_fontsize=8,
    )

    output = Path(output_path).expanduser()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    return fig
