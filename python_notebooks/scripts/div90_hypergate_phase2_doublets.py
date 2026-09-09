#!/usr/bin/env python3
"""Exploratory per-sample Scrublet audit on all 22,338 available DIV90 cells.

Scores are an orthogonal diagnostic, not exclusion criteria. The original object
is already filtered, and the actual expected doublet rate is not documented.
"""
from pathlib import Path
import sys
import os
import json
_runtime=Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_hypergate_sst_pv_phase2')
sys.dont_write_bytecode=True
os.environ['MPLCONFIGDIR']=str(_runtime/'cache/matplotlib')
os.environ['NUMBA_CACHE_DIR']=str(_runtime/'cache/numba')
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_div90_hypergate_phase2 import OUT, PROJECT_ROOT, write_json


def run():
    sys.path.insert(0,str(OUT/'cache/deps'))
    source=ad.read_h5ad(PROJECT_ROOT/'results/python_anndata/varela_div90.h5ad',backed='r')
    # Backed AnnData layers are loaded by anndata; this is the unnormalised counts layer.
    all_counts=source.layers['counts']
    rows=[]; summaries=[]
    fig,axes=plt.subplots(2,3,figsize=(10,5.6),constrained_layout=True)
    for ax,(sample,ix) in zip(axes.flat,source.obs.groupby('orig.ident',observed=True).indices.items()):
        path=OUT/'tables'/f'doublets_{sample}.tsv.gz'
        simpath=OUT/'tables'/f'doublets_simulated_{sample}.tsv.gz'
        metapath=OUT/'provenance'/f'doublets_{sample}.json'
        if path.exists() and simpath.exists() and metapath.exists():
            d=pd.read_csv(path,sep='\t');sim=pd.read_csv(simpath,sep='\t').score.to_numpy();meta=json.loads(metapath.read_text())
        else:
            a=ad.AnnData(X=all_counts[ix,:].copy(),obs=source.obs.iloc[ix][['orig.ident']].copy(),var=source.var.copy())
            assert a.X.min()>=0 and np.allclose(a.X.data,np.round(a.X.data))
            sc.pp.scrublet(a,expected_doublet_rate=.05,sim_doublet_ratio=2.,random_state=90,use_approx_neighbors=False,n_prin_comps=30,verbose=True)
            d=a.obs[['doublet_score','predicted_doublet']].copy().reset_index(names='cell_id');d['sample']=str(sample)
            sim=np.asarray(a.uns['scrublet']['doublet_scores_sim'])
            threshold=float(a.uns['scrublet'].get('threshold',np.nan))
            meta=dict(sample=str(sample),n_input=len(ix),threshold=threshold,expected_doublet_rate_assumption=.05,predicted_fraction=float(d.predicted_doublet.mean()),scanpy_version=sc.__version__,random_seed=90,n_components=30)
            d.to_csv(path,sep='\t',index=False);pd.DataFrame({'score':sim}).to_csv(simpath,sep='\t',index=False);write_json(metapath,meta)
        rows.append(d);summaries.append(meta)
        bins=np.linspace(0,max(sim.max(),d.doublet_score.max()),45)
        ax.hist(d.doublet_score,bins=bins,density=True,histtype='step',color='#276a95',label='Observed cells')
        ax.hist(sim,bins=bins,density=True,histtype='step',color='#b85f28',label='Simulated doublets')
        threshold=meta.get('threshold')
        if threshold is not None and np.isfinite(threshold):ax.axvline(threshold,color='.35',ls='--',lw=.8)
        ax.set(title=str(sample),xlabel='Exploratory Scrublet score',ylabel='Density');ax.spines[['top','right']].set_visible(False)
        print(json.dumps(meta),flush=True)
    axes.flat[0].legend(frameon=False,fontsize=7)
    fig.suptitle('Doublet sensitivity: original full sample context; no cells excluded',fontsize=12)
    for ext in ['pdf','svg','png']:fig.savefig(OUT/f'figures/{ext}/supplement_doublet_score_distributions.{ext}',dpi=600)
    plt.close(fig)
    scores=pd.concat(rows,ignore_index=True);assert scores.cell_id.is_unique and len(scores)==source.n_obs
    scores.to_csv(OUT/'tables/doublet_scores_all_div90.tsv.gz',sep='\t',index=False)
    pd.DataFrame(summaries).to_csv(OUT/'tables/doublet_sample_summary.tsv',sep='\t',index=False)
    write_json(OUT/'provenance/doublet_method.json',dict(method='Scanpy Scrublet, simulated doublets, per original sample across all available cell types',reference='https://pmc.ncbi.nlm.nih.gov/articles/PMC6625319/',n_input_cells=len(scores),expected_doublet_rate=.05,expected_rate_is_assumed=True,seed=90,limitation='Available input is already QC-filtered. Continuous score and automated threshold are exploratory, not validated doublet truth; homotypic doublets can escape detection. No removal or alteration of primary states.',dependencies_location=str(OUT/'cache/deps')))
    source.file.close()


if __name__=='__main__':run()
