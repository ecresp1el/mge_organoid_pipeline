#!/usr/bin/env python3
"""Merge per-parent DIV90 subclustering audit outputs."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import pandas as pd


PROJECT_ROOT_DEFAULT = "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run-dir",
        default=f"{PROJECT_ROOT_DEFAULT}/results/div90_parent_cluster_subclustering_audit/div90_parent_cluster_subclustering_audit_v1",
    )
    parser.add_argument("--job-id", default="", help="Optional Slurm array job ID to record in the merged report.")
    return parser.parse_args()


def read_tables(run_dir: Path, name: str, gz: bool = False) -> pd.DataFrame:
    paths = sorted((run_dir / "per_parent").glob(f"cluster_*/tables/{name}"))
    frames = []
    for path in paths:
        frames.append(pd.read_csv(path, sep="\t"))
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def write_tsv_gz(df: pd.DataFrame, path: Path) -> None:
    with gzip.open(path, "wt") as handle:
        df.to_csv(handle, sep="\t", index=False)


def main() -> int:
    args = parse_args()
    run_dir = Path(args.run_dir).expanduser().resolve()
    table_dir = run_dir / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)

    rec = read_tables(run_dir, "div90_parent_cluster_subclustering_recommendations.tsv")
    summary = read_tables(run_dir, "div90_parent_cluster_subclustering_resolution_summary.tsv")
    marker_support = read_tables(run_dir, "div90_parent_cluster_subclustering_marker_support.tsv")
    top_marker_paths = sorted((run_dir / "per_parent").glob("cluster_*/tables/div90_parent_cluster_subclustering_top_markers.tsv.gz"))
    top_markers = pd.concat([pd.read_csv(path, sep="\t") for path in top_marker_paths], ignore_index=True) if top_marker_paths else pd.DataFrame()

    rec.to_csv(table_dir / "div90_parent_cluster_subclustering_recommendations.tsv", sep="\t", index=False)
    summary.to_csv(table_dir / "div90_parent_cluster_subclustering_resolution_summary.tsv", sep="\t", index=False)
    marker_support.to_csv(table_dir / "div90_parent_cluster_subclustering_marker_support.tsv", sep="\t", index=False)
    write_tsv_gz(top_markers, table_dir / "div90_parent_cluster_subclustering_top_markers.tsv.gz")

    lines = [
        "# DIV90 Parent Cluster Subclustering Audit",
        "",
        "Merged from per-parent Slurm array outputs. Parent cluster `12` dividing cells was intentionally excluded.",
        "",
        "## Slurm Array",
        "",
        f"- Array job: `{args.job_id or 'not recorded'}`",
        "- Parent clusters audited: `0-11`",
        "",
        "## Recommendation Counts",
        "",
        rec["recommendation"].value_counts().to_string() if not rec.empty else "No recommendations found.",
        "",
        "## Candidate Subcluster Calls",
        "",
    ]
    candidates = rec.loc[rec["recommendation"] == "candidate_subcluster"].copy() if not rec.empty else pd.DataFrame()
    if candidates.empty:
        lines.append("No parent cluster passed the current marker-support and confounding filters.")
    else:
        for row in candidates.sort_values("parent_cluster_id").itertuples():
            lines.append(
                f"- Parent `{row.parent_cluster_id}`: resolution `{row.recommended_resolution}`, "
                f"{row.recommended_n_subclusters} candidate subclusters; marker-supported fraction "
                f"{row.best_marker_supported_fraction:.2f}; max technical effect {row.best_max_technical_effect:.2f}."
            )
    lines.extend(["", "## Keep Parent Only", ""])
    keep = rec.loc[rec["recommendation"] == "keep_parent_only"].copy() if not rec.empty else pd.DataFrame()
    for row in keep.sort_values("parent_cluster_id").itertuples():
        lines.append(f"- Parent `{row.parent_cluster_id}`: {row.reason}")

    (run_dir / "div90_parent_cluster_subclustering_audit_merged_report.md").write_text("\n".join(lines) + "\n")
    print(f"Wrote merged outputs under {table_dir}")
    print(run_dir / "div90_parent_cluster_subclustering_audit_merged_report.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
