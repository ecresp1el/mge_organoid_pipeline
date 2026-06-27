#!/usr/bin/env python3
"""Render the final cross-study PV precursor ON/OFF marker-expression panel."""

from __future__ import annotations

import argparse
from pathlib import Path

from mge_organoid_python.cross_study_marker_expression import (
    ALL_MARKER_GENES,
    PDF_EXPORT_DPI,
    PNG_EXPORT_DPI,
    PV_PRECURSOR_MARKER_PANEL,
    PV_PRECURSOR_TOP_GENE_SPANS,
    SVG_EXPORT_DPI,
    default_cross_study_marker_specs,
    load_marker_expression_tables,
    marker_expression_distribution_audit_table,
    plot_dir,
    plot_marker_umap_grid,
    table_dir,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", required=True)
    parser.add_argument("--run-label", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    project_root = Path(args.project_root)
    specs = default_cross_study_marker_specs(project_root)
    data = load_marker_expression_tables(specs, project_root, args.run_label, genes=ALL_MARKER_GENES)

    output_path = plot_dir(project_root, args.run_label) / "cross_study_marker_expression_pv_precursors_on_off_target.png"
    manifest = plot_marker_umap_grid(
        data=data,
        output_path=output_path,
        genes=PV_PRECURSOR_MARKER_PANEL.genes,
        specs=specs,
        title="Cross-study PV Precursors/OFF-target marker expression",
        on_genes_for_divider=PV_PRECURSOR_MARKER_PANEL.on_genes,
        gene_group_labels=PV_PRECURSOR_MARKER_PANEL.gene_group_labels(),
        top_gene_spans=PV_PRECURSOR_TOP_GENE_SPANS,
    )
    manifest.insert(0, "plot_token", "pv_precursors_on_off_target")
    manifest.insert(1, "plot_path", str(output_path))
    manifest.insert(2, "panel_id", PV_PRECURSOR_MARKER_PANEL.panel_id)
    manifest.insert(3, "panel_label", PV_PRECURSOR_MARKER_PANEL.panel_label)
    manifest["excluded_study_ids"] = ""
    manifest.to_csv(
        table_dir(project_root, args.run_label) / "cross_study_marker_expression_pv_precursors_on_off_target_plot_manifest.tsv",
        sep="\t",
        index=False,
    )

    audit = marker_expression_distribution_audit_table(
        data,
        genes=PV_PRECURSOR_MARKER_PANEL.genes,
        plot_token="pv_precursors_on_off_target",
        panel_id=PV_PRECURSOR_MARKER_PANEL.panel_id,
        panel_label=PV_PRECURSOR_MARKER_PANEL.panel_label,
        gene_group_labels=PV_PRECURSOR_MARKER_PANEL.gene_group_labels(),
    )
    audit.to_csv(
        table_dir(project_root, args.run_label) / "cross_study_marker_expression_pv_precursors_on_off_target_distribution_audit.tsv",
        sep="\t",
        index=False,
    )

    print(f"Rendered: {output_path}")
    print(f"PNG dpi: {PNG_EXPORT_DPI}")
    print(f"PDF dpi: {PDF_EXPORT_DPI}")
    print(f"SVG dpi: {SVG_EXPORT_DPI}")
    print(manifest[["study_label", "plot_filter_label", "n_cells_plotted"]].drop_duplicates().to_string(index=False))
    print(
        manifest[
            ["gene", "color_scale_min", "color_scale_min_rule", "color_scale_blue_start", "color_scale_max"]
        ]
        .drop_duplicates()
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
