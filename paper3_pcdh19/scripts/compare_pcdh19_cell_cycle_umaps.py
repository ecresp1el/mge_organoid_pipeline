#!/usr/bin/env python3
"""Python before/after S+G2M regression, preserving measured counts and baseline IDs."""
from __future__ import annotations
import argparse,datetime,gc,hashlib,importlib.metadata,json,os,shutil,time
from pathlib import Path
import numpy as np

ROOT=Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19')
REPO=Path(__file__).resolve().parents[1]


def regress_blocks(matrix, covariates, block_genes=64):
    """OLS residuals with intercept, in-place float32; never use residuals as counts."""
    design=np.column_stack([np.ones(len(covariates)),np.asarray(covariates,dtype=np.float64)])
    if not np.isfinite(design).all() or np.linalg.matrix_rank(design)!=design.shape[1]:
        raise ValueError('Nonfinite or rank-deficient cycle covariates')
    q,r=np.linalg.qr(design,mode='reduced')
    coefficients=np.empty((design.shape[1],matrix.shape[1]),dtype=np.float64)
    for start in range(0,matrix.shape[1],block_genes):
        stop=min(start+block_genes,matrix.shape[1]);y=matrix[:,start:stop]
        qt=q.T@y
        coefficients[:,start:stop]=np.linalg.solve(r,qt)
        matrix[:,start:stop]=y-q@qt
    return coefficients


def resolve_cycle_genes(symbols,feature_ids,config):
    """Audit unique conserved mouse-symbol matches and existing named aliases."""
    lookup={}
    for symbol,feature in zip(symbols,feature_ids):lookup.setdefault(str(symbol).upper(),[]).append(str(feature))
    aliases=config['cycle_reference']['renamed_aliases'];resolved={};audit=[]
    for program in config['programs']:
        if program['name'] not in ['S_phase','G2M']:continue
        matched=[]
        for symbol in program['genes']:
            query=aliases.get(symbol,symbol);hits=lookup.get(query.upper(),[])
            status='matched' if len(hits)==1 else ('missing' if not hits else 'ambiguous')
            audit.append(dict(phase=program['name'],requested=symbol,lookup=query,status=status,features=hits))
            if len(hits)==1:matched.append(hits[0])
        resolved[program['name']]=list(dict.fromkeys(matched))
        if len(matched)<10:raise ValueError('Insufficient cell-cycle feature coverage')
    return resolved,audit


def run(output, min_genes=0):
    import anndata as ad
    import pandas as pd
    import scanpy as sc
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scipy import sparse
    from matplotlib.lines import Line2D
    output.mkdir(parents=True,exist_ok=False)
    started=time.monotonic()
    def event(stage,**kwargs):
        row=dict(stage=stage,UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),elapsed_seconds=time.monotonic()-started,**kwargs)
        with (output/'progress.jsonl').open('a') as handle:handle.write(json.dumps(row)+'\n')
        (output/'status.json').write_text(json.dumps(row,indent=2));print(json.dumps(row),flush=True)
    base=ROOT/'results/rebuild_20261103'
    if json.loads((base/'cleanup_status.json').read_text())['stage']!='complete':
        raise RuntimeError('Cleanup must finish before this rebuild starts')
    shutil.copy2(__file__,output/Path(__file__).name)
    config_path=REPO/'config/hicat_validation_programs.json'
    shutil.copy2(config_path,output/config_path.name)
    cfg=json.loads(config_path.read_text());seed=20260903;sc.settings.n_jobs=int(os.environ.get('SLURM_CPUS_PER_TASK','8'))
    event('load_raw')
    raw_path=ROOT/'results/primary_processing/02_qc_filtering/02_qc_filtering_20260830_124611_97e1bb5/objects/pcdh19_step02_qc_filtered.h5ad'
    a=ad.read_h5ad(raw_path)
    assert a.shape==(446349,19071)
    baseline=np.load(base/'baseline/original_step06_embeddings.npz',allow_pickle=False)
    ids=np.asarray(a.obs_names,dtype=str);assert np.array_equal(ids,baseline['cell_id'])
    original_cells=a.n_obs
    if min_genes:
        keep=a.obs['n_genes_by_counts'].to_numpy()>=min_genes
        a=a[keep].copy();ids=np.asarray(a.obs_names,dtype=str)
        event('minimum_detected_genes_filter',minimum=min_genes,excluded=int((~keep).sum()),retained=a.n_obs)
        baseline=None  # Refit both arms on the newly retained cell universe.
    oldvar=pd.read_csv(base/'baseline/original_step06_var.tsv.gz',sep='\t',index_col=0)
    assert np.array_equal(a.var_names.astype(str),oldvar.index.astype(str))
    if not min_genes and 'highly_variable' in oldvar:
        hvg=oldvar['highly_variable'].astype(str).str.lower().isin(['true','1']).to_numpy()
        hvg_source='saved Step06 variable-gene flags'
    else:
        # Published Step06 reloads raw counts and does not attach the working HVG flags.
        # Recreate exactly its raw-count, sample-aware feature-selection call.
        event('rederive_original_HVGs',flavor='seurat_v3',n_top_genes=3000,batch_key='technical_sample_id')
        sc.pp.highly_variable_genes(a,flavor='seurat_v3',n_top_genes=3000,batch_key='technical_sample_id',subset=False)
        hvg=a.var['highly_variable'].to_numpy(dtype=bool)
        hvg_source='rederived on retained raw counts using original Step06 settings'
    assert hvg.sum()==3000
    pd.DataFrame({'feature_id':a.var_names,'highly_variable':hvg}).to_csv(output/'selected_HVGs.tsv',sep='\t',index=False)
    if 'gene_symbol' not in a.var:raise ValueError('Gene symbols missing')
    genes,audit=resolve_cycle_genes(a.var.gene_symbol.to_numpy(),a.var_names,cfg)
    (output/'cell_cycle_gene_coverage.json').write_text(json.dumps(audit,indent=2))
    event('normalize_and_score',s_genes=len(genes['S_phase']),g2m_genes=len(genes['G2M']))
    # The published raw checkpoint stays read-only; only this in-memory object changes.
    sc.pp.normalize_total(a,target_sum=10000);sc.pp.log1p(a)
    sc.tl.score_genes_cell_cycle(a,s_genes=genes['S_phase'],g2m_genes=genes['G2M'],random_state=seed)
    obs=a.obs.copy();obs.index.name='cell_id'
    covariates=obs[['S_score','G2M_score']].to_numpy()
    diagnostic=ad.AnnData(X=np.asarray(a[:,hvg].X.toarray() if sparse.issparse(a.X) else a.X[:,hvg],dtype=np.float32),obs=obs.copy(),var=a.var.loc[hvg].copy())
    # Keep a small explicit marker context from original, unregressed expression.
    marker_symbols=['Mki67','Top2a','Sox2','Ascl1','Dlx2','Gad1','Gad2','Lhx6','Sox6','Snap25','Syt1','Rbfox3','Pcdh19']
    markers={}
    for gene in marker_symbols:
        matches=np.flatnonzero(a.var.gene_symbol.astype(str).str.upper().to_numpy()==gene.upper())
        if len(matches)==1:
            x=a.X[:,matches[0]];markers[gene]=np.asarray(x.toarray() if sparse.issparse(x) else x).ravel()
    del a;gc.collect()
    if min_genes:
        event('refit_unregressed_baseline')
        uncorrected=diagnostic.copy()
        sc.pp.scale(uncorrected,zero_center=False,max_value=10)
        sc.tl.pca(uncorrected,n_comps=50,zero_center=True,svd_solver='arpack',random_state=seed)
        baseline={'X_pca':uncorrected.obsm['X_pca'].copy()}
        del uncorrected;gc.collect()
        uncorrected=ad.AnnData(X=np.empty((len(obs),0),dtype=np.float32),obs=obs.copy())
        uncorrected.obsm['X_pca']=baseline['X_pca']
        sc.pp.neighbors(uncorrected,n_neighbors=30,n_pcs=50,use_rep='X_pca',metric='cosine',random_state=seed)
        sc.tl.umap(uncorrected,min_dist=.3,spread=1.,random_state=seed)
        baseline['X_umap']=uncorrected.obsm['X_umap'].copy()
        del uncorrected;gc.collect()
    event('OLS_regression',matrix_shape=list(diagnostic.shape),covariates=['S_score','G2M_score'])
    coefficients=regress_blocks(diagnostic.X,covariates)
    sc.pp.scale(diagnostic,zero_center=False,max_value=10)
    assert np.isfinite(diagnostic.X).all()
    event('PCA',components=50)
    sc.tl.pca(diagnostic,n_comps=50,zero_center=True,svd_solver='arpack',random_state=seed)
    # Release the dense analysis copy after PCA; counts are never duplicated on disk.
    model=dict(regression_coefficients=coefficients,feature_ids=np.asarray(diagnostic.var_names,dtype=str),
               scale_sd=diagnostic.var['std'].to_numpy(),PCs=diagnostic.varm['PCs'],
               pca_variance=diagnostic.uns['pca']['variance'],pca_variance_ratio=diagnostic.uns['pca']['variance_ratio'])
    np.savez_compressed(output/'regression_pca_model.npz',**model)
    corrected_pca=diagnostic.obsm['X_pca'].copy()
    np.savez_compressed(output/'corrected_PCA_checkpoint.npz',cell_id=ids,X_pca=corrected_pca)
    obs.to_csv(output/'cells_and_cycle_scores.tsv.gz',sep='\t',compression='gzip')
    del diagnostic;gc.collect()
    embedding=ad.AnnData(X=np.empty((len(obs),0),dtype=np.float32),obs=obs.copy())
    embedding.obsm['X_pca']=corrected_pca
    event('neighbors',neighbors=30,metric='cosine')
    sc.pp.neighbors(embedding,n_neighbors=30,n_pcs=50,use_rep='X_pca',metric='cosine',random_state=seed)
    event('UMAP',min_dist=.3,spread=1.,seed=seed)
    sc.tl.umap(embedding,min_dist=.3,spread=1.,random_state=seed)
    before=baseline['X_umap'];after=embedding.obsm['X_umap']
    np.savez_compressed(output/'before_after_embeddings.npz',cell_id=ids,X_umap_before=before,X_umap_after=after,
                        X_pca_before=baseline['X_pca'],X_pca_after=corrected_pca)
    sparse.save_npz(output/'corrected_connectivities.npz',embedding.obsp['connectivities'])
    sparse.save_npz(output/'corrected_distances.npz',embedding.obsp['distances'])
    obs.to_csv(output/'cells_and_cycle_scores.tsv.gz',sep='\t',compression='gzip')
    np.savez_compressed(output/'unregressed_marker_context.npz',cell_id=ids,**markers)
    sample_field='technical_sample_id';assert sample_field in obs
    sample=obs[sample_field].astype(str).to_numpy();samples=sorted(set(sample),key=lambda v:int(v.rsplit('-',1)[-1]))
    assert len(samples)==12
    rng=np.random.RandomState(seed)
    selected=np.concatenate([rng.choice(np.flatnonzero(sample==sn),min(6000,int((sample==sn).sum())),replace=False) for sn in samples])
    rng.shuffle(selected)
    palette=plt.get_cmap('tab20');colors={sn:palette(i) for i,sn in enumerate(samples)}
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42,'svg.fonttype':'none'})
    def points(ax,xy,mask=None):
        ax.scatter(xy[selected,0],xy[selected,1],s=.2,c='#DADDE0',linewidths=0,rasterized=True)
        draw=selected if mask is None else selected[mask[selected]]
        color=[colors[v] for v in sample[draw]] if mask is None else colors[sample[draw[0]]] if len(draw) else '#333333'
        ax.scatter(xy[draw,0],xy[draw,1],s=.3,c=color if mask is None else [color],linewidths=0,rasterized=True)
        ax.set_xticks([]);ax.set_yticks([]);ax.set_xlabel('UMAP 1');ax.set_ylabel('UMAP 2')
    event('render_comparison')
    fig,axes=plt.subplots(1,2,figsize=(11,5))
    for ax,xy,title in zip(axes,[before,after],['Before · normalized/scaled, no cell-cycle regression','After · S and G2M scores regressed']):
        points(ax,xy);ax.set_title(title,fontsize=10)
    handles=[Line2D([0],[0],marker='o',ls='',color=colors[sn],label=sn,markersize=4) for sn in samples]
    fig.legend(handles=handles,loc='lower center',ncol=6,frameon=False)
    fig.suptitle('E14.5 mouse MGE · same {:,} cells · sample overlays'.format(len(obs)),fontsize=12)
    fig.tight_layout(rect=[0,.12,1,.92])
    for ext in ['png','pdf']:fig.savefig(output/('UMAP_before_after_samples.'+ext),dpi=250)
    plt.close(fig)
    fig,axes=plt.subplots(6,4,figsize=(12,17))
    for i,sn in enumerate(samples):
        row=i//2;col=(i%2)*2
        for ax,xy,label in zip(axes[row,col:col+2],[before,after],['Before','After']):
            points(ax,xy,sample==sn);ax.set_title(sn+' · '+label,fontsize=8)
    fig.suptitle('One sample highlighted at a time · same displayed cells',fontsize=12)
    fig.tight_layout(rect=[0,0,1,.97])
    for ext in ['png','pdf']:fig.savefig(output/('UMAP_before_after_individual_samples.'+ext),dpi=200)
    plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(12,7))
    for col,name in enumerate(['S_score','G2M_score','phase']):
        vals=obs[name].to_numpy()
        for row,(xy,label) in enumerate([(before,'Before'),(after,'After')]):
            ax=axes[row,col]
            if name=='phase':
                phase_colors={'G1':'#B8BCC1','S':'#D89A43','G2M':'#8A63A8'}
                c=[phase_colors[str(v)] for v in vals[selected]]
                ax.scatter(xy[selected,0],xy[selected,1],s=.25,c=c,linewidths=0,rasterized=True)
                ax.legend(handles=[Line2D([0],[0],marker='o',ls='',color=v,label=k,markersize=3) for k,v in phase_colors.items()],frameon=False,fontsize=6)
            else:
                limits=np.percentile(vals,[1,99]);scatter=ax.scatter(xy[selected,0],xy[selected,1],s=.25,c=vals[selected],vmin=limits[0],vmax=limits[1],cmap='viridis',linewidths=0,rasterized=True)
                fig.colorbar(scatter,ax=ax,shrink=.65)
            ax.set_title(label+' · '+name);ax.set_xticks([]);ax.set_yticks([])
    fig.suptitle('Original cell-cycle scores on both embeddings · scores are not refitted after regression')
    fig.tight_layout()
    for ext in ['png','pdf']:fig.savefig(output/('UMAP_before_after_cell_cycle.'+ext),dpi=220)
    plt.close(fig)
    lineage=[gene for gene in ['Sox2','Ascl1','Gad1','Snap25'] if gene in markers]
    if lineage:
        fig,axes=plt.subplots(2,len(lineage),figsize=(3*len(lineage),7),squeeze=False)
        for col,gene in enumerate(lineage):
            vals=markers[gene];upper=max(float(np.percentile(vals,99)),.01)
            for row,(xy,label) in enumerate([(before,'Before'),(after,'After')]):
                ax=axes[row,col]
                scatter=ax.scatter(xy[selected,0],xy[selected,1],s=.25,c=vals[selected],vmin=0,vmax=upper,cmap='viridis',linewidths=0,rasterized=True)
                ax.set_title(label+' · '+gene);ax.set_xticks([]);ax.set_yticks([])
                fig.colorbar(scatter,ax=ax,shrink=.65)
        fig.suptitle('Lineage context · original expression on both embeddings')
        fig.tight_layout()
        for ext in ['png','pdf']:fig.savefig(output/('UMAP_before_after_lineage_markers.'+ext),dpi=220)
        plt.close(fig)
    correlations=[]
    for arm,pcs in [('before',baseline['X_pca']),('after',corrected_pca)]:
        for j in range(pcs.shape[1]):
            for score in ['S_score','G2M_score']:
                correlations.append(dict(arm=arm,PC=j+1,score=score,pearson_r=float(np.corrcoef(pcs[:,j],obs[score])[0,1])))
    pd.DataFrame(correlations).to_csv(output/'PC_cycle_correlations.tsv',sep='\t',index=False)
    pd.crosstab(obs[sample_field],obs['phase']).to_csv(output/'sample_phase_counts.tsv',sep='\t')
    parameters=dict(input=str(raw_path),cells=len(obs),genes=19071,selected_HVGs=3000,hvg_source=hvg_source,seed=seed,
                    minimum_detected_genes=min_genes,additional_QC_exclusions=original_cells-len(obs),before=('refitted on minimum-gene-filtered cells' if min_genes else 'saved original Step06 UMAP/PCA, all cell IDs verified'),after='HVG log1p(CP10K) OLS residuals for intercept + S + G2M; SD scale/clip10; centered ARPACK PCA50; cosine k30; UMAP min_dist .3 spread1',
                    cells_removed=0,raw_input_modified=False,regressed_covariates=['S_score','G2M_score'],
                    batch_integration=False,clustering='not run in this comparison; Python HiCAT consensus is next',
                    packages={p:importlib.metadata.version(p) for p in ['scanpy','anndata','numpy','scipy','umap-learn','scikit-learn']})
    (output/'parameters.json').write_text(json.dumps(parameters,indent=2))
    (output/'README.md').write_text('# Cell-cycle regression comparison\n\nBoth arms contain the identical retained cells. For a minimum-gene filter both embeddings are refitted. After refits the high-dimensional graph after S/G2M regression. Same HVGs, cell IDs, normalization target and graph/UMAP settings. Separate layouts need not have the same orientation. Sample colors identify composition, not proof of batch correction. Scores/markers are from unregressed expression. Regression can remove real developmental biology; inspect lineage-marker context before adopting it for Python HiCAT. No new clusters or Allen labels are claimed yet. Raw counts remain the input for future DE and MapMyCells.\n')
    assert np.isfinite(after).all() and after.shape==before.shape
    event('complete',output=str(output),cells=len(obs),coordinates_finite=True)
    (output/'SUCCESS.json').write_text(json.dumps(parameters,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'results/rebuild_20261103/cell_cycle_comparison')
    parser.add_argument('--min-genes',type=int,default=0)
    args=parser.parse_args();run(args.output,args.min_genes)
