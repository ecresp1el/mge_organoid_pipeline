"""Render individual fine-cluster locations from a completed hierarchy AnnData.

This read-only review supplement avoids overlapping labels and repeated colors
in a global many-cluster UMAP. It does not recompute embeddings or clustering.
The original scientific package remains unchanged; use a new output directory.
"""
from pathlib import Path
import argparse
import hashlib
import json
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def render_locations(input_path, output_dir):
    """Save paginated highlight views for every baseline fine cluster.

    Parameters
    ----------
    input_path : pathlib.Path
        Completed coarse/fine H5AD. Read backed; expression is not loaded.
    output_dir : pathlib.Path
        New review directory, separate from immutable scientific outputs.

    Returns
    -------
    pathlib.Path
        PDF with eight individual fine-cluster views per page. Coordinates,
        memberships, cell counts, page index and input hash are saved alongside.
    """
    output_dir=Path(output_dir)
    output_dir.mkdir(parents=True)
    obj=ad.read_h5ad(input_path,backed='r')
    try:
        xy=np.asarray(obj.obsm['X_umap_step06_display'])
        fine=obj.obs.hicat_fine_baseline.astype(str).to_numpy()
        coarse=obj.obs.hicat_coarse_baseline.astype(str).to_numpy()
        cells=obj.obs_names.to_numpy()
    finally:
        obj.file.close()
    clusters=sorted(set(fine))
    records=[]
    pdf=output_dir/'fine_cluster_locations.pdf'
    with PdfPages(pdf) as book:
        for start in range(0,len(clusters),8):
            fig,axes=plt.subplots(2,4,figsize=(16,9),sharex=True,sharey=True)
            for ax,label in zip(axes.flat,clusters[start:start+8]):
                selected=fine==label
                ax.scatter(xy[:,0],xy[:,1],s=1,color='#d4d4d4',rasterized=True)
                ax.scatter(xy[selected,0],xy[selected,1],s=4,color='#b2182b',rasterized=True)
                parent=np.unique(coarse[selected])
                if len(parent)!=1:raise ValueError('Invalid fine-parent nesting')
                ax.set_title(label+' | n='+str(selected.sum()),fontsize=11)
                ax.set_xticks([]);ax.set_yticks([])
                records.append(dict(page=start//8+1,fine_cluster=label,coarse_parent=parent[0],cells=int(selected.sum())))
            for ax in list(axes.flat)[len(clusters[start:start+8]):]:ax.axis('off')
            fig.suptitle('Fine-cluster locations | red = selected cluster, gray = all 12,000 pilot cells\nExisting Step 06 UMAP; display only',fontsize=14)
            fig.tight_layout(rect=(0,0,1,.93))
            book.savefig(fig)
            fig.savefig(output_dir/('locations_%02d.png'%(start//8+1)),dpi=160)
            plt.close(fig)
    pd.DataFrame(records).to_csv(output_dir/'cluster_page_index.tsv',sep='\t',index=False)
    pd.DataFrame(dict(cell_id=cells,coarse_cluster=coarse,fine_cluster=fine,
                      umap_1=xy[:,0],umap_2=xy[:,1])).to_csv(output_dir/'plotted_coordinates.tsv.gz',sep='\t',index=False)
    digest=hashlib.sha256()
    with Path(input_path).open('rb') as handle:
        for block in iter(lambda:handle.read(8*1024**2),b''):digest.update(block)
    (output_dir/'provenance.json').write_text(json.dumps(dict(input=str(input_path),input_sha256=digest.hexdigest(),
        annotation_performed=False,clustering_recomputed=False,embedding_recomputed=False,
        reason='Individual cluster highlights to avoid overlapping global labels',clusters=len(clusters),status='IN_REVIEW'),indent=2)+'\n')
    return pdf


def main():
    """Render only a completed object's saved baseline labels and coordinates."""
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args()
    print(render_locations(args.input,args.output_dir))


if __name__=='__main__':
    main()
