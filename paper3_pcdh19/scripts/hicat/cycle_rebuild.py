"""Minimum-gene-filtered HiCAT rerun with S/G2M residual clustering space.

Retains the prior Allen feature eligibility, recursive fitter, graph, DE,
80% sampling and consensus rules. Only the clustering/mapping representation
is globally S/G2M-regressed, SD-scaled and clipped at +10. Expression DE and
feature eligibility use original ln(1+CPM). No dense all-gene residual file.
"""
import argparse
from itertools import combinations
import gc
import json
import os
from pathlib import Path
import shutil
import tarfile
import tempfile
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
import transcriptomic_clustering as tc
from .fullscale_engine import FullScaleEngine, allocation_adapters
from .consensus_core import Membership
from .provenance import Progress, write_json, sha256

RAW=Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/primary_processing/02_qc_filtering/02_qc_filtering_20260830_124611_97e1bb5/objects/pcdh19_step02_qc_filtered.h5ad')


def fit_cycle_model(x, design):
    """OLS global coefficients and sample residual SD without dense all-gene X."""
    design=np.asarray(design,dtype=np.float64)
    if np.linalg.matrix_rank(design)!=design.shape[1]:raise ValueError('Rank-deficient covariates')
    xtx=design.T@design
    xty=np.asarray(x.T@design).T
    beta=np.linalg.solve(xtx,xty)
    squares=np.zeros(x.shape[1],dtype=np.float64)
    for start in range(0,x.shape[0],4096):
        block=x[start:start+4096]
        squares+=np.asarray(block.multiply(block).sum(axis=0)).ravel() if sparse.issparse(block) else np.square(block).sum(axis=0)
    residual_ss=np.maximum(squares-np.sum(beta*(xtx@beta),axis=0),0)
    sd=np.sqrt(residual_ss/(x.shape[0]-1));sd[sd<1e-12]=1
    return beta,sd


def corrected_block(x, design, beta, sd, columns):
    y=x[:,columns].toarray() if sparse.issparse(x) else np.asarray(x[:,columns]).copy()
    y=np.asarray(y,dtype=np.float64)
    y-=design@beta[:,columns]
    y/=sd[columns]
    # Same upper clipping as Scanpy scale(max_value=10, zero_center=False).
    np.minimum(y,10,out=y)
    if not np.isfinite(y).all():raise ValueError('Nonfinite residual representation')
    return y


class CycleEngine(FullScaleEngine):
    """Prior recursive graph/DE fitter with an explicit corrected PCA input."""
    def __init__(self,*args,**kwargs):
        self.model=kwargs.pop('model')
        super().__init__(*args,**kwargs)

    def _representation(self,obj,node,marker_mask=None):
        cfg=self.cfg
        with self.progress.track('hicat.features',node=node.name,cell_cycle_residual_PCA=True):
            if marker_mask is None:
                means,variances,eligible=tc.get_means_vars_genes(obj,low_thresh=cfg['low_thresh'],min_cells=cfg['min_gene_cells'])
                eligible=np.asarray(eligible,dtype=bool)
                if eligible.sum()<3 or not np.isfinite(means).all() or not np.isfinite(variances).all():raise ValueError('Invalid eligible gene statistics')
                hvgs=tc.highly_variable_genes(obj,means,variances,eligible,max_genes=cfg['max_hvgs'])
                self._table(node/'eligible_gene_statistics.tsv',pd.DataFrame({'gene_id':obj.var_names[eligible],'linear_cpm_mean':means,'linear_cpm_variance':variances}),index=False)
            else:hvgs=marker_mask
            columns=np.flatnonzero(np.asarray(hvgs,dtype=bool))
            self._table(node/'feature_selection.tsv',pd.DataFrame({'gene_id':obj.var_names,'selected':np.asarray(hvgs,dtype=bool)}),index=False)
        with self.progress.track('hicat.corrected_pca',node=node.name,cells=obj.n_obs,genes=len(columns)):
            rows=obj.obs.global_row.to_numpy(dtype=int)
            y=corrected_block(obj.X,self.model['design'][rows],self.model['beta'],self.model['sd'],columns)
            working=ad.AnnData(y,obs=pd.DataFrame(index=obj.obs_names),var=obj.var.iloc[columns].copy(),dtype=y.dtype)
            components,ratios,variances,means=tc.pca(working,cell_select=np.arange(obj.n_obs),n_comps=cfg['pca_components'],svd_solver='randomized',random_state=self.candidate['seed'])
            self._table(node/'pca_components_all.tsv',components);self._table(node/'pca_center.tsv',means)
            self._table(node/'pca_variance.tsv',pd.DataFrame({'variance_ratio':ratios,'variance':variances}),index=False)
            selected=tc.dimension_reduction.filter_components(components,variances.copy(),ratios.copy(),method=cfg['pc_filter'],max_pcs=cfg['max_retained_pcs'])
            if not selected.shape[1]:raise ValueError('PC filter retained zero components')
            self._table(node/'pca_components_retained.tsv',selected)
            projected=tc.project(working,selected,means)
            np.savez_compressed(node/'projection.npz',coordinates=projected.X,cell_ids=np.asarray(obj.obs_names,dtype=str))
            return projected


def map_corrected(obj,rows,labels,columns,model,block_cells=4096):
    """Previous mean-prototype Pearson assignment in corrected marker space."""
    unique=np.unique(labels);all_labels=np.full(obj.n_obs,-1,dtype=np.int32);all_labels[rows]=labels
    scores=np.full(obj.n_obs,np.nan);heldout=np.flatnonzero(all_labels<0)
    if len(unique)==1:
        all_labels[heldout]=unique[0]
        return all_labels,scores,np.empty((1,0)),unique
    if len(columns)<2:raise ValueError('Insufficient mapping markers')
    prototypes=np.zeros((len(unique),len(columns)))
    def values(idx):return corrected_block(obj.X[idx],model['design'][idx],model['beta'],model['sd'],columns)
    for j,key in enumerate(unique):
        group=rows[labels==key]
        for start in range(0,len(group),block_cells):prototypes[j]+=values(group[start:start+block_cells]).sum(axis=0)
        prototypes[j]/=len(group)
    centered=prototypes-prototypes.mean(axis=1,keepdims=True);norm=np.linalg.norm(centered,axis=1)
    for start in range(0,len(heldout),block_cells):
        idx=heldout[start:start+block_cells];x=values(idx);x-=x.mean(axis=1,keepdims=True)
        denominator=np.linalg.norm(x,axis=1)[:,None]*norm[None,:]
        corr=np.divide(x@centered.T,denominator,out=np.zeros((len(idx),len(unique))),where=denominator!=0)
        corr[~np.isfinite(corr)]=0;choice=corr.argmax(axis=1)
        all_labels[idx]=unique[choice];scores[idx]=corr[np.arange(len(idx)),choice]
    return all_labels,scores,prototypes,unique


def read_count_source(path, limit=None):
    """Read only counts and needed metadata; skip newer AnnData uns encodings."""
    import h5py
    from .consensus_benchmark import decode
    with h5py.File(path,'r') as h:
        matrix=h['X'];shape=tuple(map(int,matrix.attrs['shape']));n=shape[0] if limit is None else min(limit,shape[0])
        ptr=matrix['indptr'][:n+1];last=int(ptr[-1])
        x=sparse.csr_matrix((matrix['data'][:last],matrix['indices'][:last],ptr),shape=(n,shape[1]))
        on=h['obs'];vn=h['var'];ids=decode(on[on.attrs['_index']])[:n];genes=decode(vn[vn.attrs['_index']])
        obs=pd.DataFrame(index=pd.Index(ids,name='cell_id'))
        for key in ['technical_sample_id','design_group','genotype','sex','total_counts','n_genes_by_counts','pct_counts_mt']:
            if key in on:obs[key]=decode(on[key])[:n]
        var=pd.DataFrame({'gene_symbol':decode(vn['gene_symbol'])},index=genes)
    return ad.AnnData(x,obs=obs,var=var,dtype=x.dtype)


def prepare(root):
    cfg=json.loads((root/'config.json').read_text());comparison=Path(cfg['comparison'])
    if not (comparison/'SUCCESS.json').exists():raise ValueError('Corrected comparison is incomplete')
    obj=read_count_source(RAW);keep=obj.obs.n_genes_by_counts.to_numpy()>=cfg['minimum_detected_genes'];obj=obj[keep].copy()
    scores=pd.read_csv(comparison/'cells_and_cycle_scores.tsv.gz',sep='\t',index_col=0)
    if not np.array_equal(obj.obs_names.astype(str),scores.index.astype(str)):raise ValueError('Comparison cell IDs differ')
    if obj.n_obs!=cfg['expected_cells']:raise ValueError('Unexpected filtered population')
    obj.obs['global_row']=np.arange(obj.n_obs)
    totals=np.asarray(obj.X.sum(axis=1)).ravel();obj.X=obj.X.astype(np.float64)
    obj.X.data*=np.repeat(1e6/totals,np.diff(obj.X.indptr));obj.X.data=np.log1p(obj.X.data)
    design=np.column_stack([np.ones(obj.n_obs),scores[['S_score','G2M_score']].to_numpy()])
    beta,sd=fit_cycle_model(obj.X,design)
    np.savez_compressed(root/'cycle_model.npz',design=design,beta=beta,sd=sd,cell_id=np.asarray(obj.obs_names,dtype=str),gene_id=np.asarray(obj.var_names,dtype=str))
    obj.uns.clear();obj.uns['processing']='ln(1+CPM), no regression in X; global S/G2M residual model used only for representation/mapping; original gene detection/DE retained'
    obj.write_h5ad(root/'full_log1p_cpm.h5ad',compression='lzf')
    obj.obs.to_csv(root/'cell_universe.tsv.gz',sep='\t',index_label='cell_id')
    write_json(root/'PREPARED.json',dict(cells=obj.n_obs,genes=obj.n_vars,minimum_genes=cfg['minimum_detected_genes'],additional_QC_exclusions=int((~keep).sum()),input_sha256=sha256(RAW),normalized_sha256=sha256(root/'full_log1p_cpm.h5ad'),model_sha256=sha256(root/'cycle_model.npz')))


def load(root):
    for relative,expected in json.loads((root/'FROZEN.json').read_text()).items():
        if sha256(root/relative)!=expected:raise ValueError('Frozen source changed: '+relative)
    obj=ad.read_h5ad(root/'full_log1p_cpm.h5ad')
    model=dict(np.load(root/'cycle_model.npz',allow_pickle=False))
    if not np.array_equal(obj.obs_names,model['cell_id']) or not np.array_equal(obj.var_names,model['gene_id']):raise ValueError('Model IDs mismatch')
    obj.obs=pd.DataFrame({'global_row':np.arange(obj.n_obs)},index=obj.obs_names);obj.uns.clear()
    return obj,model


def fit(root,index):
    cfg=json.loads((root/'config.json').read_text());destination=root/('fullfit' if index<0 else 'iterations/%03d'%index)
    destination.mkdir(parents=True,exist_ok=True)
    if (destination/'SUCCESS.json').exists():return
    seeds=pd.read_csv(root/'seeds.tsv',sep='\t');sampling_seed=20260908 if index<0 else int(seeds.iloc[index].sampling_seed);fit_seed=20260909 if index<0 else int(seeds.iloc[index].fit_seed)
    from annoy_seed64.worker import install
    install(root/'code/annoy_seed64')
    obj,model=load(root)
    rows=np.arange(obj.n_obs) if index<0 else np.sort(np.random.default_rng(sampling_seed).choice(obj.n_obs,round(.8*obj.n_obs),replace=False))
    selected=obj if index<0 else obj[rows].copy()
    # Node graphs/projections exist only on job-local scratch, avoiding the old
    # hundreds-of-GB accumulation of 100 full recursive trace directories.
    with tempfile.TemporaryDirectory(prefix='pcdh_hicat_',dir=os.environ.get('SLURM_TMPDIR',os.environ.get('TMPDIR','/tmp'))) as temporary:
        temporary=Path(temporary);progress=Progress(destination/'progress')
        engine=CycleEngine(cfg['fitting'],dict(name='fullfit' if index<0 else 'iteration_%03d'%index,seed=fit_seed,score_thresh=150),temporary/'engine',progress,model=model)
        labels,summary,evidence,means=engine.run(selected)
        _,labels=np.unique(labels,return_inverse=True);labels=labels.astype(np.int32)+1
        columns=obj.var_names.get_indexer(sorted(engine.markers))
        np.save(destination/'sampled_source_rows.npy',rows);np.save(destination/'sampled_labels.npy',labels);np.save(destination/'marker_columns.npy',columns)
        if index<0:all_labels=labels;scores=np.full(obj.n_obs,np.nan)
        else:
            del selected;gc.collect()
            all_labels,scores,prototypes,keys=map_corrected(obj,rows,labels,columns,model)
            np.savez_compressed(destination/'mapping_model.npz',prototype_means=prototypes,cluster_ids=keys,marker_columns=columns)
        assert np.all(all_labels>0) and np.array_equal(all_labels[rows],labels)
        np.save(destination/'all_cell_labels.npy',all_labels)
        pd.DataFrame({'cell_id':obj.obs_names,'cluster':all_labels,'heldout_correlation':scores}).to_csv(destination/'all_cell_assignments.tsv.gz',sep='\t',index=False)
        if index>=0:sparse.save_npz(destination/'membership_B.npz',Membership([all_labels]).b)
        # Retain numerical feature/PCA/DE audit tables and node summaries in one
        # compressed archive; omit regenerable large graphs and repeated row lists.
        with tarfile.open(destination/'engine_audit.tar.gz','w:gz') as archive:
            for p in sorted((temporary/'engine').rglob('*')):
                if p.is_file() and p.suffix!='.npz' and p.name not in ('annoy.index','cell_membership.tsv','merged_membership.json') and not p.name.endswith('_membership.json'):
                    archive.add(p,arcname=str(p.relative_to(temporary)))
        write_json(destination/'SUCCESS.json',dict(summary=summary,cells=obj.n_obs,sampled=len(rows),sampling_seed=sampling_seed,fit_seed=fit_seed,model_sha256=sha256(root/'cycle_model.npz'),all_labels_sha256=sha256(destination/'all_cell_labels.npy'),cell_order_exact=True))


def aggregate(root):
    cfg=json.loads((root/'config.json').read_text());out=root/'final';out.mkdir(exist_ok=True)
    paths=[root/'iterations'/('%03d'%i) for i in range(cfg['iterations'])]
    for p in paths+[root/'fullfit']:
        if not (p/'SUCCESS.json').exists():raise ValueError('Missing completed fit: '+str(p))
        assert sha256(p/'all_cell_labels.npy')==json.loads((p/'SUCCESS.json').read_text())['all_labels_sha256']
    original=np.load(paths[0]/'all_cell_labels.npy');counts=np.unique(original,return_counts=True)[1]
    if sum(int(x)**2 for x in counts)<10**9:raise ValueError('Smaller Allen consensus branch requires qualification; not silently substituted')
    progress=Progress(out/'progress')
    with progress.track('consensus.membership',iterations=len(paths)):
        membership=Membership.from_blocks([sparse.load_npz(p/'membership_B.npz') for p in paths],4096)
        sparse.save_npz(out/'membership_B.npz',membership.b)
    initial=np.load(root/'fullfit/all_cell_labels.npy')
    with progress.track('consensus.merge_by_co'):
        labels,actions=membership.merge_by_co(initial);write_json(out/'merge_actions.json',actions)
    with progress.track('consensus.refine'):
        labels,trace,cells,medians=membership.refine(labels,min_cells=20)
        np.save(out/'labels_after_refinement.npy',labels);write_json(out/'refinement_trace.json',trace)
        cells.to_csv(out/'cell_consensus_diagnostics.tsv.gz',sep='\t',index_label='source_row');medians.to_csv(out/'cluster_consensus_diagnostics.tsv',sep='\t')
    obj,model=load(root);columns=np.load(root/'fullfit/marker_columns.npy')
    engine=CycleEngine(cfg['fitting'],dict(name='final_DE',seed=20260908,score_thresh=150),out/'engine',progress,model=model)
    groups={int(k):np.flatnonzero(labels==k).tolist() for k in np.unique(labels)}
    with progress.track('consensus.final_DE_merge'):
        # Final DE priorities use corrected marker PCs; significance/detection
        # are evaluated on unregressed ln(1+CPM), as in each recursive fit.
        with allocation_adapters():
            projected=engine._representation(obj,engine.directory,pd.Series(obj.var_names.isin(obj.var_names[columns]),index=obj.var_names))
            groups=engine._merge(obj,projected,groups,labels,engine.directory)
            labels=engine.labels(list(groups.values()),obj.n_obs)
            groups={k:np.flatnonzero(labels==k) for k in sorted(set(labels))}
            means,present,variances=tc.get_cluster_means(obj,groups,labels,low_th=1)
    thresholds=dict(cfg['fitting']['thresholds'])
    for key in ['low_thresh','min_genes','score_thresh']:thresholds.pop(key)
    pairs=list(combinations(groups,2))
    evidence=tc.de_pairs_ebayes(pairs,means,variances,present,{k:len(v) for k,v in groups.items()},thresholds) if pairs else pd.DataFrame()
    if len(evidence):
        evidence['meets_separation']=(evidence.score>=150)&(evidence.num>=5)
        for col in ['up_genes','down_genes']:evidence[col]=evidence[col].map(lambda x:json.dumps(list(x)))
    evidence.to_csv(out/'pairwise_DE.tsv',sep='\t');means.to_csv(out/'cluster_mean_log1p_cpm.tsv.gz',sep='\t');present.to_csv(out/'cluster_detection.tsv.gz',sep='\t')
    pd.DataFrame({'cell_id':obj.obs_names,'cluster':labels}).to_csv(out/'all_cell_assignments.tsv.gz',sep='\t',index=False)
    pd.crosstab(labels,pd.read_csv(root/'cell_universe.tsv.gz',sep='\t').technical_sample_id).to_csv(out/'cluster_by_sample_counts.tsv',sep='\t')
    render_clusters(root,obj.obs_names,labels,cfg)
    write_json(out/'SUCCESS.json',dict(cells=obj.n_obs,clusters=len(groups),real_iterations=len(paths),annotations_locked=False,residual_DE_pairs=int((~evidence.meets_separation).sum()) if len(evidence) else 0,status='COMPLETED_PENDING_BIOLOGICAL_REVIEW'))


def render_clusters(root, ids, labels, cfg):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    coords=np.load(Path(cfg['comparison'])/'before_after_embeddings.npz',allow_pickle=False)
    if not np.array_equal(np.asarray(ids,dtype=str),coords['cell_id']):raise ValueError('Plot IDs differ')
    obs=pd.read_csv(root/'cell_universe.tsv.gz',sep='\t');sample=obs.technical_sample_id.astype(str).to_numpy()
    samples=sorted(set(sample),key=lambda v:int(v.rsplit('-',1)[-1]));rng=np.random.RandomState(20260903)
    draw=np.concatenate([rng.choice(np.flatnonzero(sample==s),min(6000,int((sample==s).sum())),replace=False) for s in samples]);rng.shuffle(draw)
    xy=coords['X_umap_after'];keys=sorted(set(labels));palette=plt.get_cmap('turbo');colors={k:palette(i/max(1,len(keys)-1)) for i,k in enumerate(keys)}
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,axes=plt.subplots(1,2,figsize=(10,4.6))
    axes[0].scatter(xy[draw,0],xy[draw,1],s=.25,c=[colors[k] for k in labels[draw]],linewidths=0,rasterized=True)
    for k in keys:
        center=np.median(xy[labels==k],axis=0);axes[0].text(*center,str(k),fontsize=7,ha='center',bbox={'facecolor':'white','alpha':.7,'edgecolor':'none','pad':1})
    axes[0].set_title('HiCAT consensus · anonymous cluster IDs')
    sample_colors={k:plt.get_cmap('tab20')(i) for i,k in enumerate(samples)}
    axes[1].scatter(xy[draw,0],xy[draw,1],s=.25,c=[sample_colors[k] for k in sample[draw]],linewidths=0,rasterized=True);axes[1].set_title('Same cells · sample overlay')
    from matplotlib.lines import Line2D
    fig.legend(handles=[Line2D([0],[0],marker='o',ls='',color=sample_colors[k],label='JZ-'+k.rsplit('-',1)[-1],markersize=4) for k in samples],loc='lower center',ncol=6,frameon=False)
    for ax in axes:ax.set_xticks([]);ax.set_yticks([]);ax.set_xlabel('UMAP 1');ax.set_ylabel('UMAP 2')
    fig.suptitle('E14.5 mouse MGE · {:,} cells · ≥500 genes · S/G2M-regressed representation'.format(len(ids)))
    fig.tight_layout(rect=[0,.13,1,.95])
    for ext in ['png','pdf','svg']:fig.savefig(root/'final'/('HiCAT_consensus_UMAP.'+ext),dpi=250)
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--run-dir',type=Path,required=True);parser.add_argument('--stage',choices=['prepare','fit','aggregate'],required=True);parser.add_argument('--iteration',type=int,default=-1)
    args=parser.parse_args()
    for relative,expected in json.loads((args.run_dir/'FROZEN.json').read_text()).items():
        if sha256(args.run_dir/relative)!=expected:raise ValueError('Frozen source changed: '+relative)
    for path,expected in json.loads((args.run_dir/'REFERENCE.json').read_text()).items():
        if sha256(path)!=expected:raise ValueError('Pinned Allen source changed: '+path)
    import subprocess,sys
    if subprocess.check_output([sys.executable,'-m','pip','freeze'],text=True)!=(args.run_dir/'environment.txt').read_text():raise ValueError('HiCAT environment changed')
    if args.stage=='prepare':prepare(args.run_dir)
    elif args.stage=='fit':fit(args.run_dir,args.iteration)
    else:aggregate(args.run_dir)


if __name__=='__main__':main()
