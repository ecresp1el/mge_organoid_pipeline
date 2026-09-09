"""Regenerate the three existing DIV90 guidance figures with SST expression added.

Run with PYTHONPATH=python_notebooks/src using the mge-organoid-python environment.
Original coordinates, membership and selection criteria are reused unchanged.
"""
from pathlib import Path
import json
import gc
import sys
from zipfile import ZipFile, ZIP_DEFLATED

import numpy as np
import pandas as pd
import render_div90_loupe_subcluster_guidance_expression as atlas
import render_div90_loupe_pooled_subpopulation_violins as violins
from div90_hypergate_paths import PROJECT_ROOT, REPO_ROOT


def save_png(fig, outdir, stem, dpi):
    """Export only the requested raster figures and release rendering buffers."""
    path = outdir / "figures/png" / f"{stem}.png"
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    atlas.plt.close(fig)
    fig.clear()
    gc.collect()
    return [path]


def main():
    root = REPO_ROOT
    source = PROJECT_ROOT
    out = PROJECT_ROOT / "final_figures/div90_guidance_with_sst"
    atlas.setup_dirs(out)
    atlas.save_figure = save_png
    original_genes = list(atlas.GENES)
    atlas.GENES = [*original_genes, "SST"]
    violins.GENES = [*violins.GENES, "SST"]
    data = pd.read_csv(source / violins.INPUT_RELATIVE, sep="\t", dtype={"cell_id": str})
    expression, matches = atlas.extract_marker_expression_from_h5ad(
        atlas.div90_spec(source), out / "tables/expression_with_sst.tsv.gz",
        project_root=source, genes=atlas.GENES,
    )
    assert matches["matched"].all(), matches.to_string()
    assert expression["cell_id"].is_unique
    # Confirm that the source matrix reproduces the expression used previously.
    check = data.loc[data.expression_available].merge(expression[["cell_id", *original_genes]], on="cell_id",
                       validate="one_to_one", suffixes=("", "_source"))
    assert len(check) == int(data.expression_available.sum())
    for gene in original_genes:
        np.testing.assert_allclose(check[gene], check[gene + "_source"], atol=1e-6)
    data = data.merge(expression[["cell_id", "SST"]], on="cell_id", how="left", validate="one_to_one")
    assert np.isfinite(data.loc[data.expression_available, atlas.GENES].to_numpy()).all()
    joined = out / "tables/joined_with_sst.tsv.gz"
    data.to_csv(joined, sep="\t", index=False)
    matches.to_csv(out / "tables/gene_matches.tsv", sep="\t", index=False)
    print("Source expression matches all original genes; SST joined successfully.", flush=True)
    # Preserve the established feature scale for comparability with the originals.
    for set_id in atlas.SETS:
        atlas.ATLAS_STEMS[set_id] += "_with_sst"
        atlas.render_feature_atlas(data.loc[data.set_id == set_id], set_id, out, 4.0, 600)
        print(f"Rendered {set_id} atlas", flush=True)
    original, _ = violins.load_original_umap(source / violins.H5AD_RELATIVE)
    keep = original.loc[~original.cluster_id.isin(violins.EXCLUDED_DIV90_CLUSTERS), "cell_id"]
    filtered = data.loc[data.cell_id.isin(keep) & data.expression_available]
    cortical, subpallial, _, _ = violins.select_populations(filtered)
    assert (len(cortical), len(subpallial)) == (1071, 526)
    ymax = max(4, int(np.ceil(max(cortical[violins.GENES].max().max(),
                                subpallial[violins.GENES].max().max()))))
    output = out / "figures/png/div90_guidance_receptor_subpopulations_with_sst.png"
    sys.argv = [__file__, "--input", str(joined), "--output", str(output),
                "--violin-ymax", str(ymax), "--dpi", "600"]
    violins.main()
    (out / "provenance/validation.json").write_text(json.dumps({
        "source_h5ad": str(source / violins.H5AD_RELATIVE),
        "original_gene_expression_matches": True,
        "cortical_selected_n": len(cortical), "subpallial_selected_n": len(subpallial),
        "feature_display_range": [0, 4], "violin_ymax": ymax,
        "sst_detected_cells_by_lineage": data.groupby("set_id").SST.apply(lambda x: int((x > 0).sum())).to_dict(),
        "selection_criteria": "Unchanged; SST added as expression only",
    }, indent=2))
    archive = PROJECT_ROOT / "results/div90_cortical_subpallial_figures_with_SST.zip"
    pngs = sorted((out / "figures/png").glob("*.png"))
    assert len(pngs) == 3
    with ZipFile(archive, "w", ZIP_DEFLATED) as z:
        for path in pngs:
            z.write(path, path.name)
    with ZipFile(archive) as z:
        assert z.testzip() is None
    print(f"ZIP: {archive}", flush=True)


if __name__ == "__main__":
    main()
