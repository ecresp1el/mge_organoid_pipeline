#!/usr/bin/env python3
"""Write a reproducible output manifest for DIV30/DIV90 URD runs.

The URD scripts intentionally reuse some historical filenames, including
``div30_first_urd_*`` for the shared Matrix Market bundle consumed by the R
pipeline. This audit layer makes the run contents explicit without renaming
live pipeline artifacts.
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class ExpectedArtifact:
    rel_path: str
    stage: str
    artifact_type: str
    required: bool
    description: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument(
        "--mode",
        required=True,
        choices=("div30_first", "div30_root10", "div90_jia"),
        help="Expected artifact set to audit.",
    )
    parser.add_argument("--outdir", type=Path, default=None)
    parser.add_argument("--selected-pct", default="10")
    parser.add_argument(
        "--fail-on-missing-required",
        action="store_true",
        help="Exit 1 if any required expected artifact is absent.",
    )
    return parser.parse_args()


def common_initial_urd() -> list[ExpectedArtifact]:
    return [
        ExpectedArtifact("div30_first_urd_object.rds", "initial_urd", "rds", True, "Saved URD object after PCA, diffusion map, and flood pseudotime."),
        ExpectedArtifact("tables/div30_first_urd_parameters.tsv", "initial_urd", "table", True, "Explicit imposed URD parameters."),
        ExpectedArtifact("tables/div30_first_urd_pseudotime.tsv", "initial_urd", "table", True, "Per-cell pseudotime and metadata export."),
        ExpectedArtifact("tables/div30_first_urd_summary.tsv", "initial_urd", "table", True, "Initial URD summary with key output paths."),
        ExpectedArtifact("tables/div30_first_urd_marker_correlations.tsv", "initial_urd", "table", True, "Marker/program correlations with URD pseudotime."),
        ExpectedArtifact("tables/div30_first_urd_stage_timings.tsv", "initial_urd", "table", True, "Stage timing and memory breadcrumbs."),
        ExpectedArtifact("tables/slurm_resource_heartbeat.tsv", "slurm", "table", False, "Slurm resource heartbeat if submitted through a monitored template."),
        ExpectedArtifact("plots/div30_first_urd_umap_pseudotime.png", "initial_urd", "plot", True, "UMAP colored by URD pseudotime."),
        ExpectedArtifact("plots/div30_first_urd_pseudotime_by_paper_cluster.png", "initial_urd", "plot", False, "Pseudotime by paper cluster when paper annotations are available."),
    ]


def decision_report(prefix: str = "lineage_decision_report") -> list[ExpectedArtifact]:
    base = prefix.rstrip("/")
    tables = [
        "root_annotation_composition.tsv",
        "pseudotime_ordering_by_annotation.tsv",
        "top_negative_pseudotime_genes.tsv",
        "top_positive_pseudotime_genes.tsv",
        "branch_structure_status.tsv",
        "decision_genes_between_branches.tsv",
        "flood_stability_summary.tsv",
        "gene_cascade_heatmap_matrix.tsv",
    ]
    plots = [
        "umap_pseudotime.png",
        "diffusion_map_pseudotime.png",
        "diffusion_map_annotation.png",
        "diffusion_map_annotation_grid.png",
        "flood_stability.png",
        "tree_visualization.png",
        "gene_cascade_heatmap.png",
        "lineage_decision_tree.png",
    ]
    out = [
        ExpectedArtifact(f"{base}/lineage_decision_tree_report.md", "decision_report", "report", True, "Post-URD lineage decision report."),
    ]
    out.extend(
        ExpectedArtifact(f"{base}/tables/{name}", "decision_report", "table", True, "Decision-report table.")
        for name in tables
    )
    out.extend(
        ExpectedArtifact(f"{base}/plots/{name}", "decision_report", "plot", True, "Decision-report plot.")
        for name in plots
    )
    return out


def lineage_tree(prefix: str) -> list[ExpectedArtifact]:
    base = prefix.rstrip("/")
    return [
        ExpectedArtifact(f"{base}/div30_urd_lineage_tree_object.rds", "lineage_tree", "rds", True, "URD object after tree construction."),
        ExpectedArtifact(f"{base}/urd_lineage_tree_report.md", "lineage_tree", "report", True, "Tree report."),
        ExpectedArtifact(f"{base}/tables/lineage_tree_parameters.tsv", "lineage_tree", "table", True, "Tree-construction parameters."),
        ExpectedArtifact(f"{base}/tables/lineage_tree_stage_timings.tsv", "lineage_tree", "table", True, "Tree stage timings."),
        ExpectedArtifact(f"{base}/tables/tree_tip_mapping.tsv", "lineage_tree", "table", True, "Requested tip ID to label mapping."),
        ExpectedArtifact(f"{base}/tables/tree_status.tsv", "lineage_tree", "table", True, "Tree status summary."),
        ExpectedArtifact(f"{base}/tables/tree_segment_joins.tsv", "lineage_tree", "table", True, "Tree segment joins."),
        ExpectedArtifact(f"{base}/tables/tree_tip_composition.tsv", "lineage_tree", "table", True, "Tip composition table."),
        ExpectedArtifact(f"{base}/tables/branch_specific_genes.tsv", "lineage_tree", "table", True, "Tip-vs-other branch marker table."),
        ExpectedArtifact(f"{base}/tables/pseudotime_logistic_parameters.tsv", "lineage_tree", "table", True, "URD pseudotime logistic parameters."),
        ExpectedArtifact(f"{base}/plots/pseudotime_logistic.png", "lineage_tree", "plot", True, "Pseudotime transition logistic plot."),
        ExpectedArtifact(f"{base}/plots/urd_tree_annotation.png", "lineage_tree", "plot", True, "URD tree colored by annotation."),
        ExpectedArtifact(f"{base}/plots/urd_tree_annotation_grid.png", "lineage_tree", "plot", True, "Faceted tree annotation overlays, one group highlighted per panel."),
        ExpectedArtifact(f"{base}/plots/urd_tree_pseudotime.png", "lineage_tree", "plot", True, "URD tree colored by pseudotime."),
    ]


def finalized_tree_report(prefix: str, stage: str) -> list[ExpectedArtifact]:
    base = prefix.rstrip("/")
    return [
        ExpectedArtifact(f"{base}/urd_lineage_tree_report.md", stage, "report", True, "Finalized tree report."),
        ExpectedArtifact(f"{base}/tables/tree_status.tsv", stage, "table", True, "Tree status summary."),
        ExpectedArtifact(f"{base}/tables/tree_segment_joins.tsv", stage, "table", True, "Tree segment joins."),
        ExpectedArtifact(f"{base}/tables/tree_tip_composition.tsv", stage, "table", True, "Tip composition table."),
        ExpectedArtifact(f"{base}/tables/branch_specific_genes.tsv", stage, "table", True, "Tip-vs-other branch marker table."),
        ExpectedArtifact(f"{base}/plots/urd_tree_annotation.png", stage, "plot", True, "URD tree colored by annotation."),
        ExpectedArtifact(f"{base}/plots/urd_tree_annotation_grid.png", stage, "plot", True, "Faceted tree annotation overlays, one group highlighted per panel."),
        ExpectedArtifact(f"{base}/plots/urd_tree_pseudotime.png", stage, "plot", True, "URD tree colored by pseudotime."),
    ]


def div90_inputs() -> list[ExpectedArtifact]:
    return [
        ExpectedArtifact("inputs/div30_first_urd_counts.mtx", "input_export", "matrix", True, "Matrix Market counts consumed by shared URD R script."),
        ExpectedArtifact("inputs/div30_first_urd_features.tsv", "input_export", "table", True, "Feature table consumed by shared URD R script."),
        ExpectedArtifact("inputs/div30_first_urd_barcodes.tsv", "input_export", "table", True, "Selected cell barcode table."),
        ExpectedArtifact("inputs/div30_first_urd_cell_metadata.tsv", "input_export", "table", True, "Per-cell metadata including root candidates."),
        ExpectedArtifact("inputs/div30_first_urd_input_manifest.tsv", "input_export", "table", True, "Input export manifest."),
        ExpectedArtifact("inputs/div90_jia_rootscore_root_summary.tsv", "input_export", "table", True, "DIV90 root selection summary."),
        ExpectedArtifact("inputs/div90_jia_rootscore_scored_metadata.tsv", "input_export", "table", True, "DIV90 scored metadata with root flags."),
        ExpectedArtifact("inputs/div90_jia_program_markers_used.tsv", "input_export", "table", True, "Jia program markers used for scoring."),
        ExpectedArtifact("inputs/root_score_program_marker_summary.tsv", "input_export", "table", True, "Selected roots versus pool/all-cell marker/program summary."),
    ]


def jia_fig_s11(prefix: str = "jia_fig_s11_style_marker_validation_v1") -> list[ExpectedArtifact]:
    base = prefix.rstrip("/")
    genes = ["HES1", "CACNA1E", "DLX2", "DCX", "LHX8", "NR2F1", "EPHA5", "MEF2C", "CRABP1"]
    artifacts = [
        ExpectedArtifact(f"{base}/jia_fig_s11_style_marker_validation_report.md", "marker_validation", "report", True, "Jia Fig. S11-style marker validation report."),
        ExpectedArtifact(f"{base}/tables/jia_fig_s11_marker_order.tsv", "marker_validation", "table", True, "Jia marker order/specification."),
        ExpectedArtifact(f"{base}/tables/jia_fig_s11_marker_expression_summary.tsv", "marker_validation", "table", True, "Jia marker expression summary."),
        ExpectedArtifact(f"{base}/tables/jia_fig_s11_marker_expression_by_cell.tsv.gz", "marker_validation", "table", True, "Jia marker expression by tree cell."),
        ExpectedArtifact(f"{base}/plots/jia_fig_s11_style_urd_marker_validation.png", "marker_validation", "plot", True, "Combined Jia marker validation PNG."),
        ExpectedArtifact(f"{base}/plots/jia_fig_s11_style_urd_marker_validation.pdf", "marker_validation", "plot", True, "Combined Jia marker validation PDF."),
        ExpectedArtifact(f"{base}/plots/jia_fig_s11_panel_a_developmental_markers.png", "marker_validation", "plot", True, "Jia panel A developmental marker PNG."),
        ExpectedArtifact(f"{base}/plots/jia_fig_s11_panel_b_lineage_markers.png", "marker_validation", "plot", True, "Jia panel B lineage marker PNG."),
    ]
    for gene in genes:
        artifacts.append(ExpectedArtifact(f"{base}/plots/jia_fig_s11_marker_tree_overlay_{gene}.png", "marker_validation", "plot", True, f"Individual Jia marker tree overlay for {gene}, PNG."))
        artifacts.append(ExpectedArtifact(f"{base}/plots/jia_fig_s11_marker_tree_overlay_{gene}.pdf", "marker_validation", "plot", True, f"Individual Jia marker tree overlay for {gene}, PDF."))
    return artifacts


def div90_candidate_markers(prefix: str = "candidate_pv_marker_projection_v1") -> list[ExpectedArtifact]:
    base = prefix.rstrip("/")
    genes = ["MEF2C", "EPHA5", "LHX6", "CRABP1", "LHX8", "NR2F1", "NR2F2"]
    artifacts = [
        ExpectedArtifact(f"{base}/div90_candidate_marker_projection_report.md", "candidate_marker_projection", "report", True, "DIV90 candidate marker projection report."),
        ExpectedArtifact(f"{base}/tables/div90_candidate_marker_expression_by_cell.tsv.gz", "candidate_marker_projection", "table", True, "Candidate marker expression by tree cell."),
        ExpectedArtifact(f"{base}/tables/div90_candidate_marker_expression_by_cell_with_metadata.tsv.gz", "candidate_marker_projection", "table", True, "Candidate marker expression plus metadata."),
        ExpectedArtifact(f"{base}/tables/div90_marker_expression_summary_by_cluster.tsv", "candidate_marker_projection", "table", True, "Marker expression summary by cluster."),
        ExpectedArtifact(f"{base}/tables/div90_marker_expression_summary_by_tree_segment.tsv", "candidate_marker_projection", "table", True, "Marker expression summary by tree segment."),
        ExpectedArtifact(f"{base}/tables/div90_candidate_marker_expression_by_cluster_segment.tsv", "candidate_marker_projection", "table", True, "Candidate-cluster expression by segment."),
        ExpectedArtifact(f"{base}/tables/div90_candidate_lineage_proxy_scores_by_cell.tsv.gz", "candidate_marker_projection", "table", True, "Lineage proxy scores by tree cell."),
        ExpectedArtifact(f"{base}/tables/div90_candidate_lineage_proxy_scores_by_cluster_segment.tsv", "candidate_marker_projection", "table", True, "Lineage proxy scores by cluster and segment."),
        ExpectedArtifact(f"{base}/tables/div90_candidate_and_tip_marker_profiles.tsv", "candidate_marker_projection", "table", True, "Candidate and tip marker profiles."),
        ExpectedArtifact(f"{base}/tables/div90_candidate_cluster_to_tip_profile_correlations.tsv", "candidate_marker_projection", "table", True, "Candidate-to-tip profile correlations."),
        ExpectedArtifact(f"{base}/tables/div90_candidate_cluster_lineage_assignment.tsv", "candidate_marker_projection", "table", True, "Guardrailed candidate assignment summary."),
        ExpectedArtifact(f"{base}/plots/div90_candidate_marker_tree_overlays.png", "candidate_marker_projection", "plot", True, "Candidate marker overlays on tree PNG."),
        ExpectedArtifact(f"{base}/plots/div90_candidate_marker_tree_overlays.pdf", "candidate_marker_projection", "plot", True, "Candidate marker overlays on tree PDF."),
        ExpectedArtifact(f"{base}/plots/div90_candidate_and_tip_marker_profile_heatmap.png", "candidate_marker_projection", "plot", True, "Candidate/tip profile heatmap PNG."),
        ExpectedArtifact(f"{base}/plots/div90_candidate_and_tip_marker_profile_heatmap.pdf", "candidate_marker_projection", "plot", True, "Candidate/tip profile heatmap PDF."),
    ]
    for gene in genes:
        artifacts.append(ExpectedArtifact(f"{base}/plots/div90_candidate_marker_tree_overlay_{gene}.png", "candidate_marker_projection", "plot", True, f"Individual DIV90 candidate marker tree overlay for {gene}, PNG."))
        artifacts.append(ExpectedArtifact(f"{base}/plots/div90_candidate_marker_tree_overlay_{gene}.pdf", "candidate_marker_projection", "plot", True, f"Individual DIV90 candidate marker tree overlay for {gene}, PDF."))
    return artifacts


def root_score_outputs(selected_pct: str) -> list[ExpectedArtifact]:
    pct = str(selected_pct)
    reflood = f"reflood_top{pct}pct_v1"
    out = [
        ExpectedArtifact("div30_urd_jia_rootscore_object.rds", "root_score", "rds", True, "URD object annotated with Jia RootScore metadata."),
        ExpectedArtifact("jia_rootscore_candidate_report.md", "root_score", "report", True, "RootScore candidate report."),
        ExpectedArtifact("tables/root_score_all_cells.tsv", "root_score", "table", True, "RootScore values for all scored cells."),
        ExpectedArtifact(f"tables/root_score_candidate_cells_top{pct}pct.tsv", "root_score", "table", True, "Selected top-percent root candidate cells."),
        ExpectedArtifact("tables/root_score_candidate_counts.tsv", "root_score", "table", True, "Candidate counts by top-percent threshold."),
        ExpectedArtifact("tables/root_score_distribution.tsv", "root_score", "table", True, "RootScore distribution summary."),
        ExpectedArtifact("tables/root_score_marker_expression_by_candidate_set.tsv", "root_score", "table", True, "Marker expression by candidate set."),
        ExpectedArtifact("tables/root_score_marker_genes_used.tsv", "root_score", "table", True, "Marker genes used for RootScore."),
        ExpectedArtifact("tables/root_score_program_marker_summary.tsv", "root_score", "table", True, "Selected roots versus pool/all-cell marker/program summary."),
        ExpectedArtifact("tables/root_score_proliferation_genes_used.tsv", "root_score", "table", True, "Proliferation genes used in RootScore context."),
        ExpectedArtifact("tables/root_score_pseudotime_by_candidate_set.tsv", "root_score", "table", True, "Pseudotime summary by candidate set."),
        ExpectedArtifact("tables/root_score_root_pool_cells.tsv", "root_score", "table", True, "Cells in the designated root pool."),
        ExpectedArtifact("plots/root_score_by_pseudotime.png", "root_score", "plot", True, "RootScore by pseudotime."),
        ExpectedArtifact("plots/root_score_candidate_umap.png", "root_score", "plot", True, "Selected root candidates on UMAP."),
        ExpectedArtifact("plots/root_score_distribution.png", "root_score", "plot", True, "RootScore distribution."),
        ExpectedArtifact("plots/root_score_umap.png", "root_score", "plot", True, "RootScore on UMAP."),
        ExpectedArtifact(f"{reflood}/div30_urd_reflood_object.rds", "reflood", "rds", True, "URD object reflooded from selected RootScore roots."),
        ExpectedArtifact(f"{reflood}/tables/reflood_root_summary.tsv", "reflood", "table", True, "Reflood root and parameter summary."),
        ExpectedArtifact(f"{reflood}/tables/reflood_root_cells.tsv", "reflood", "table", True, "Selected root cells used for reflood."),
    ]
    out.extend(decision_report(f"{reflood}/lineage_decision_report"))
    return out


def expected_for_mode(mode: str, selected_pct: str) -> list[ExpectedArtifact]:
    if mode == "div30_first":
        return common_initial_urd() + decision_report() + lineage_tree("lineage_tree_radial_glia_tips_v1")
    if mode == "div30_root10":
        return root_score_outputs(selected_pct)
    if mode == "div90_jia":
        return (
            div90_inputs()
            + common_initial_urd()
            + decision_report()
            + lineage_tree("lineage_tree_jia_endpoint_tips_v1")
            + finalized_tree_report("lineage_tree_cluster_number_name_v1", "cluster_name_tree_report")
            + jia_fig_s11()
            + div90_candidate_markers()
        )
    raise ValueError(mode)


def file_record(run_root: Path, artifact: ExpectedArtifact) -> dict[str, object]:
    path = run_root / artifact.rel_path
    exists = path.exists()
    stat = path.stat() if exists else None
    return {
        "modeled_stage": artifact.stage,
        "artifact_type": artifact.artifact_type,
        "required": str(artifact.required).upper(),
        "status": "present" if exists else "missing",
        "relative_path": artifact.rel_path,
        "size_bytes": stat.st_size if stat else "",
        "modified_time": datetime.fromtimestamp(stat.st_mtime).isoformat(timespec="seconds") if stat else "",
        "description": artifact.description,
    }


def discover_files(run_root: Path) -> Iterable[dict[str, object]]:
    if not run_root.exists():
        return []
    rows = []
    for path in sorted(p for p in run_root.rglob("*") if p.is_file()):
        stat = path.stat()
        rel = path.relative_to(run_root).as_posix()
        rows.append(
            {
                "relative_path": rel,
                "size_bytes": stat.st_size,
                "modified_time": datetime.fromtimestamp(stat.st_mtime).isoformat(timespec="seconds"),
                "suffix": "".join(path.suffixes),
            }
        )
    return rows


def write_tsv(rows: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if rows:
        fieldnames = list(rows[0].keys())
    else:
        fieldnames = ["status"]
        rows = [{"status": "no_rows"}]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def summarize(rows: list[dict[str, object]], discovered_rows: list[dict[str, object]], mode: str, run_root: Path) -> list[dict[str, object]]:
    required_rows = [row for row in rows if row["required"] == "TRUE"]
    missing_required = [row for row in required_rows if row["status"] != "present"]
    by_stage = {}
    for row in rows:
      stage = row["modeled_stage"]
      if stage not in by_stage:
          by_stage[stage] = {"expected": 0, "present": 0, "missing_required": 0}
      by_stage[stage]["expected"] += 1
      if row["status"] == "present":
          by_stage[stage]["present"] += 1
      if row["required"] == "TRUE" and row["status"] != "present":
          by_stage[stage]["missing_required"] += 1
    summary = [
        {"metric": "mode", "value": mode},
        {"metric": "run_root", "value": str(run_root)},
        {"metric": "expected_artifacts", "value": len(rows)},
        {"metric": "present_expected_artifacts", "value": sum(1 for row in rows if row["status"] == "present")},
        {"metric": "missing_required_artifacts", "value": len(missing_required)},
        {"metric": "discovered_files", "value": len(discovered_rows)},
    ]
    for stage, counts in sorted(by_stage.items()):
        summary.append({"metric": f"stage_{stage}_expected", "value": counts["expected"]})
        summary.append({"metric": f"stage_{stage}_present", "value": counts["present"]})
        summary.append({"metric": f"stage_{stage}_missing_required", "value": counts["missing_required"]})
    return summary


def main() -> int:
    args = parse_args()
    run_root = args.run_root.expanduser().resolve()
    outdir = (args.outdir.expanduser().resolve() if args.outdir else run_root / "tables")
    expected = expected_for_mode(args.mode, args.selected_pct)
    rows = [file_record(run_root, artifact) for artifact in expected]
    discovered_rows = list(discover_files(run_root))
    summary_rows = summarize(rows, discovered_rows, args.mode, run_root)

    write_tsv(rows, outdir / "urd_output_manifest.tsv")
    write_tsv(summary_rows, outdir / "urd_output_manifest_summary.tsv")
    write_tsv(discovered_rows, outdir / "urd_output_discovered_files.tsv")

    missing_required = [row for row in rows if row["required"] == "TRUE" and row["status"] != "present"]
    print(f"[URDOutputAudit] mode={args.mode} run_root={run_root}")
    print(f"[URDOutputAudit] expected={len(rows)} present={sum(1 for row in rows if row['status'] == 'present')} missing_required={len(missing_required)}")
    print(f"[URDOutputAudit] manifest={outdir / 'urd_output_manifest.tsv'}")
    if missing_required:
        print("[URDOutputAudit] missing required artifacts:", file=sys.stderr)
        for row in missing_required[:50]:
            print(f"  - {row['relative_path']}", file=sys.stderr)
        if len(missing_required) > 50:
            print(f"  - ... {len(missing_required) - 50} more", file=sys.stderr)
    return 1 if args.fail_on_missing_required and missing_required else 0


if __name__ == "__main__":
    raise SystemExit(main())
