#!/usr/bin/env python3
"""Extract CXCR4 and PLXNA2 from the canonical Varela AnnData caches."""

from __future__ import annotations

import argparse
from pathlib import Path

from mge_organoid_python.cross_study_marker_expression import (
    default_cross_study_marker_specs,
    extract_marker_expression_from_h5ad,
)


PROJECT_ROOT_DEFAULT = Path("/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder")
GENES = ["CXCR4", "PLXNA2"]
STUDIES = ["varela_div30", "varela_div90"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT_DEFAULT)
    parser.add_argument("--outdir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    specs = {spec.study_id: spec for spec in default_cross_study_marker_specs(args.project_root)}
    args.outdir.mkdir(parents=True, exist_ok=True)
    for study_id in STUDIES:
        spec = specs[study_id]
        output = args.outdir / f"{study_id}_marker_expression.tsv.gz"
        match_output = args.outdir / f"{study_id}_gene_matches.tsv"
        table, matches = extract_marker_expression_from_h5ad(
            spec,
            output,
            project_root=args.project_root,
            genes=GENES,
        )
        matches.to_csv(match_output, sep="\t", index=False)
        print(f"study_id={study_id} n_cells={len(table)} output={output}")


if __name__ == "__main__":
    main()
