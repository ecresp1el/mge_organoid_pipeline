#!/usr/bin/env python3
"""Render cleanly named DIV90 cortical and subpallial guidance feature UMAPs."""

from __future__ import annotations

import argparse
from datetime import datetime
import hashlib
from pathlib import Path
import shutil

import pandas as pd

import render_div90_loupe_subcluster_guidance_expression as workflow


PROJECT_ROOT_DEFAULT = Path("/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder")
OUTPUT_DIR_DEFAULT = Path("/home/elcrespo/Desktop/new_figures")
STAGING_DIR_DEFAULT = Path(
    "/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline/final_figures/"
    "div90_guidance_feature_umaps_clean_staging"
)

CLEAN_SET_TITLES = {
    "cortical_only": "Cortical lineage subclusters",
    "subpallial_only": "Subpallial lineage subclusters",
}
CLEAN_STEMS = {
    "cortical_only": "div90_cortical_guidance_feature_umaps",
    "subpallial_only": "div90_subpallial_guidance_feature_umaps",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT_DEFAULT)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR_DEFAULT)
    parser.add_argument("--staging-dir", type=Path, default=STAGING_DIR_DEFAULT)
    parser.add_argument("--csv-dir", type=Path, default=None)
    parser.add_argument("--membership", type=Path, default=None)
    parser.add_argument("--dpi", type=int, default=600)
    return parser.parse_args()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    args = parse_args()
    project_root = args.project_root.resolve()
    output_dir = args.output_dir.resolve()
    staging_dir = args.staging_dir.resolve()
    csv_dir = (args.csv_dir or project_root / workflow.CSV_DIR_RELATIVE).resolve()
    membership = (args.membership or project_root / workflow.MEMBERSHIP_RELATIVE).resolve()

    output_dir.mkdir(parents=True, exist_ok=True)
    workflow.setup_dirs(staging_dir)
    workflow.SET_TITLES.update(CLEAN_SET_TITLES)
    workflow.ATLAS_STEMS.update(CLEAN_STEMS)

    data, gene_matches, mapping_qc = workflow.load_joined_data(
        project_root, staging_dir, csv_dir, membership
    )
    vmax = workflow.common_vmax(data)

    final_outputs: list[Path] = []
    for set_id in workflow.SETS:
        generated = workflow.render_feature_atlas(
            data.loc[data["set_id"] == set_id],
            set_id,
            staging_dir,
            vmax,
            args.dpi,
        )
        for source in generated:
            destination = output_dir / source.name
            shutil.copy2(source, destination)
            final_outputs.append(destination)

    gene_matches_path = output_dir / "div90_guidance_feature_umap_gene_matches.tsv"
    mapping_qc_path = output_dir / "div90_guidance_feature_umap_mapping_qc.tsv"
    gene_matches.to_csv(gene_matches_path, sep="\t", index=False)
    mapping_qc.to_csv(mapping_qc_path, sep="\t", index=False)
    final_outputs.extend([gene_matches_path, mapping_qc_path])

    readme_path = output_dir / "DIV90_GUIDANCE_FEATURE_UMAPS_README.md"
    readme_path.write_text(
        "# DIV90 cortical and subpallial guidance feature UMAPs\n\n"
        "The two feature atlases use the existing local subcluster coordinates, "
        "subcluster annotations, and unchanged DIV90 log1p(CP10K) expression values.\n\n"
        f"All feature panels share the previously established 0 to {vmax:g} "
        "log1p(CP10K) display range. No reclustering, coordinate changes, or "
        "expression recalculation was performed.\n"
    )
    final_outputs.append(readme_path)

    source_script = Path(__file__).resolve()
    copied_script = output_dir / source_script.name
    shutil.copy2(source_script, copied_script)
    final_outputs.append(copied_script)

    render_manifest = pd.DataFrame(
        [
            ("rendered_at_local", datetime.now().isoformat(timespec="seconds")),
            ("project_root", str(project_root)),
            ("source_h5ad", str(workflow.div90_spec(project_root).resolved_h5ad_path(project_root))),
            ("source_coordinate_dir", str(csv_dir)),
            ("source_membership", str(membership)),
            ("genes", ",".join(workflow.GENES)),
            ("feature_range", f"0 to {vmax:g} log1p(CP10K)"),
            ("cortical_panel_title", CLEAN_SET_TITLES["cortical_only"]),
            ("subpallial_panel_title", CLEAN_SET_TITLES["subpallial_only"]),
            ("reclustering", "none; plot-only regeneration"),
        ],
        columns=["key", "value"],
    )
    render_manifest_path = output_dir / "div90_guidance_feature_umap_render_manifest.tsv"
    render_manifest.to_csv(render_manifest_path, sep="\t", index=False)
    final_outputs.append(render_manifest_path)

    checksum_path = output_dir / "div90_guidance_feature_umaps_sha256_manifest.txt"
    checksum_path.write_text(
        "".join(f"{sha256(path)}  {path.name}\n" for path in final_outputs)
    )

    print(f"Output folder: {output_dir}")
    print(f"Shared feature range: 0 to {vmax:g} log1p(CP10K)")
    for set_id in workflow.SETS:
        print(output_dir / f"{CLEAN_STEMS[set_id]}.png")


if __name__ == "__main__":
    main()
