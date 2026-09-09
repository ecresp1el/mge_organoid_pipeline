#!/usr/bin/env python3
"""Separate DIV30 transcriptional-antecedent analysis; frozen DIV90 reference.

No pooling, cluster targets, UMAP, antibody recommendations or lineage claims.
All runtime files belong to the new Turbo output; original inputs are read-only.
"""
from __future__ import annotations

import os
from pathlib import Path
import sys
OUT = Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_pv_antecedent_v1')
for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[name] = '4'
os.environ['MPLCONFIGDIR'] = str(OUT/'cache/matplotlib')
sys.dont_write_bytecode = True
import argparse
import hashlib
import json
import shutil
import h5py
import numpy as np
import pandas as pd
from scipy import sparse, stats
from sklearn.decomposition import PCA
from sklearn.linear_model import ElasticNet
from sklearn.metrics import r2_score
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

REPO = Path(__file__).resolve().parents[2]
PROJECT = OUT.parents[1]
REF = PROJECT/'results/div90_hypergate_sst_pv_phase2'
OLD = PROJECT/'results/div90_hypergate_sst_pv'
H30 = PROJECT/'results/python_anndata/varela_div30.h5ad'
H90 = PROJECT/'results/python_anndata/varela_div90.h5ad'
PV = ['MEF2C','KCNC1','KCNC2','TAC1']
SST = ['SST','SATB1','NR2F2','CDK14','CACNG3']
PANELS = {
    'mge':['NKX2-1','NKX6-2','LHX6','SOX6','DLX1','DLX2'],
    'progenitor':['SOX2','HES1','HES5','VIM','NES'],
    'maturation':['RBFOX3','SNAP25','SYT1','STMN2','TUBB3','MAP2','DCX'],
    'progression':['DLX5','DLX6','GAD1','GAD2','SLC32A1'],
    'cycle':['MKI67','TOP2A','PCNA','MCM2','MCM5'],
}
SENS = ['MAF','MAFB','ELMO1','IGFBP4']
CONTEXT = ['ERBB4','LHX6','NKX2-1','MAF','MAFB','MEF2C','DCX','PVALB','FGFR2','PTPRS','FAT3','PTPRM']
EXCLUDED = sorted(set(PV+SST+SENS+sum(PANELS.values(), [])))
METHODS = ['direct_pv_score','neighbor_pv_score','sparse_pv_score','rank_pv_score']
SOURCES = [
 {'url':'https://www.nature.com/articles/nature25999','use':'Frozen DIV90 MEF2C interpretation, inherited phase-two evidence; mouse developmental association.'},
 {'url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC2749660/','use':'Frozen KCNC1/2 maturation association; not early fate proof.'},
 {'url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC9316417/','use':'Frozen DIV90 TAC1 and SST module reference; shared maturation caveat.'},
 {'url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC4191718/','use':'Human GE progenitors SOX2/ASCL1 versus DCX/LHX6 migratory differentiation, NKX2-1 spatial identity.'},
 {'url':'https://www.nature.com/articles/s41593-021-00940-3','use':'Human fetal GE DLX2/GAD2 differentiation and regional progenitor identity.'},
 {'url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC5374502/','use':'Operational progenitor, cell-cycle and neuronal context genes in mouse MGE model; panel not a validated developmental clock.'},
 {'url':'https://www.nature.com/articles/s41593-020-00794-1','use':'Human early SOX2 progenitors and DCX newborn neurons; neuronal maturation is not subtype.'},
 {'url':'https://doi.org/10.1371/journal.pone.0121314','use':'CSPA catalogue annotation, independently filtered for DIV30 expression; no neuronal surface-protein validation.'},
]


def clean(x):
    if isinstance(x,dict): return {str(k):clean(v) for k,v in x.items()}
    if isinstance(x,(tuple,list,np.ndarray)): return [clean(v) for v in x]
    if isinstance(x,np.generic): return clean(x.item())
    if isinstance(x,float) and not np.isfinite(x): return None
    return x


def dump(path,x): path.write_text(json.dumps(clean(x),indent=2,allow_nan=False)+'\n')
def tab(df,name): df.to_csv(OUT/'tables'/name,sep='\t',index=False)
def md_table(df):
    """Small dependency-free Markdown table renderer for the frozen environment."""
    def fmt(v):
        if isinstance(v,(float,np.floating)):return '' if not np.isfinite(v) else f'{v:.5g}'
        return str(v).replace('|','\\|').replace('\n',' ')
    lines=['| '+' | '.join(map(str,df.columns))+' |','| '+' | '.join(['---']*len(df.columns))+' |']
    lines.extend('| '+' | '.join(fmt(v) for v in row)+' |' for row in df.itertuples(index=False,name=None))
    return '\n'.join(lines)
def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(2**20),b''): h.update(b)
    return h.hexdigest()


def hvalues(x):
    if isinstance(x,h5py.Group):
        cats=hvalues(x['categories']);codes=x['codes'][:]
        return np.array([cats[i] if i>=0 else None for i in codes],object)
    return x.asstr()[:] if h5py.check_string_dtype(x.dtype) is not None else x[:]


def hmeta(path):
    with h5py.File(path,'r') as f:
        obs=pd.DataFrame({k:hvalues(f['obs'][k]) for k in f['obs']})
        genes=list(hvalues(f['var'][f['var'].attrs['_index']]))
        return obs,genes


def iter_x(path,block=1200):
    with h5py.File(path,'r') as f:
        grp=f['X'];ip=grp['indptr'][:]; n,p=grp.attrs['shape']
        for start in range(0,n,block):
            stop=min(start+block,n);a,b=int(ip[start]),int(ip[stop])
            mat=sparse.csr_matrix((grp['data'][a:b].astype('float32'),grp['indices'][a:b],ip[start:stop+1]-a),shape=(stop-start,p))
            yield start,stop,mat


def init():
    for p in ['tables','cache','logs','provenance','figures/png','figures/pdf','figures/svg']:(OUT/p).mkdir(parents=True,exist_ok=True)
    for p in [*REPO.glob('python_notebooks/scripts/div30_antecedent_*.py'),REPO/'python_notebooks/DIV90_HYPERGATE_PHASE2_EVIDENCE.md']:
        shutil.copy2(p,OUT/'provenance'/p.name)
    inputs=[REF/'cells.tsv.gz',REF/'state_parameters.json',REF/'gate_summary.json',OLD/'provenance/cspa_annotation_audit.tsv',OLD/'allowed_surface_markers.tsv']
    manifest={str(p):{'sha256':sha(p),'size':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns} for p in inputs}
    for p in [H30,H90]: manifest[str(p)]={'size':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns}
    before=OUT/'provenance/source_before.json'
    if before.exists(): assert json.loads(before.read_text())==manifest,'Original inputs changed since starting new version.'
    else: dump(before,manifest)
    dump(OUT/'provenance/references.json',SOURCES)


def reference_matrix():
    obs,genes=hmeta(H90);ref=pd.read_csv(REF/'cells.tsv.gz',sep='\t')
    positions=pd.Index(obs.cell_id).get_indexer(ref.cell_id);assert (positions>=0).all()
    inverse=np.full(len(obs),-1,dtype=int);inverse[positions]=np.arange(len(ref))
    path=OUT/'cache/reference_expression.npy'
    if not path.exists():
        mm=np.lib.format.open_memmap(path,mode='w+',dtype='float32',shape=(len(ref),len(genes)))
        for start,stop,x in iter_x(H90):
            keep=inverse[start:stop]>=0
            if keep.any():mm[inverse[start:stop][keep]]=x[keep].toarray()
        mm.flush()
    return ref,genes,np.load(path,mmap_mode='r')


def choose_features(refx,genes):
    cspa=pd.read_csv(OLD/'provenance/cspa_annotation_audit.tsv',sep='\t')
    curated=pd.read_csv(OLD/'allowed_surface_markers.tsv',sep='\t')
    # Exclude every catalogued potential surface gene from mapping, including
    # lower-confidence and noneligible categories, more than the discovery set.
    surface=set(cspa.gene.dropna())|set(curated.gene)|set(CONTEXT)
    allowed=np.array([g not in surface and g not in EXCLUDED and not g.startswith(('MT-','RPL','RPS')) for g in genes])
    detection=np.mean(refx>0,axis=0);variance=np.var(refx,axis=0,dtype=np.float64)
    ix=np.flatnonzero(allowed & (detection>=.05)&(detection<=.95))
    selected=ix[np.argsort(variance[ix])[-1000:]]
    return sorted(surface),selected


def attach_metadata(cells):
    """Use the repository's documented 9583/9853 alias; keep source IDs intact."""
    meta=pd.read_csv(REPO/'metadata/div30_div90_sample_id_to_biolabel_map.tsv',sep='\t').query("DIV == 'DIV30'").set_index('run_sample_id')
    cond=pd.read_csv(REPO/'metadata/culture_condition_user_20260909.tsv',sep='\t').query("DIV == 'DIV30'").set_index('run_sample_id')
    cells['canonical_metadata_sample']=cells['sample'].str.replace(r'^9583-','9853-',regex=True)
    assert set(cells.canonical_metadata_sample)<=set(meta.index)
    cells['biological_label']=cells.canonical_metadata_sample.map(meta.biological_label)
    cells['cell_line']=cells.biological_label.str.replace(r'_rep\d+$','',regex=True)
    cells['genetic_line']=cells.cell_line.str.replace('_old','',regex=False)
    for name in ['condition','culture_operator','relative_glucose']:cells[name]=cells.canonical_metadata_sample.map(cond[name])
    assert cells[['biological_label','cell_line','condition']].notna().all().all()
    observed=cells.groupby('canonical_metadata_sample').size()
    audit=meta[['biological_label']].copy();audit['n_cells_in_cached_DIV30']=observed.reindex(audit.index,fill_value=0)
    audit['present_in_analysis']=audit.n_cells_in_cached_DIV30>0
    audit['source_sample_id']=audit.index.to_series().str.replace(r'^9853-','9583-',regex=True)
    for name in ['condition','culture_operator','relative_glucose']:audit[name]=cond[name]
    tab(audit.reset_index(),'sample_metadata_and_presence.tsv')
    dump(OUT/'provenance/sample_alias_reconciliation.json',dict(source_ids_preserved=True,alias_rule='9583-MW-N (cached Seurat/H5AD) ↔ 9853-MW-N (original raw sample map), same numeric sample suffix',existing_repository_precedent='scripts/29_div30_jia_progenitor_metadata_umap_grid.R:sample_key_variants (lines147–154)',independent_inventory='python_notebooks/HANDOFF_shi_reference_div30_label_transfer.md:1534',observed_source_samples=sorted(cells['sample'].unique()),observed_n_samples=cells['sample'].nunique(),absent_metadata_samples=sorted(set(meta.index)-set(observed.index)),limitation='Cached DIV30 contains the six core samples only; metadata rows7–9 (79B_old) are not represented and cannot be analyzed from this frozen H5AD.'))
    return cells


def prepare():
    ref,allgenes,refx=reference_matrix()
    surface,featureix=choose_features(refx,allgenes)
    featuregenes=[allgenes[i] for i in featureix]
    chosen=sorted(set(surface+EXCLUDED+CONTEXT+featuregenes)&set(allgenes))
    obs,genes30=hmeta(H30);assert allgenes==genes30
    cols=pd.Index(allgenes).get_indexer(chosen)
    dump(OUT/'cache/genes.json',chosen)
    dump(OUT/'provenance/scoring_config.json',dict(pv=PV,sst=SST,panels=PANELS,sensitivity=SENS,excluded_genes=EXCLUDED,mapping_features=featuregenes,surface_catalogue_genes=surface,scaling='PV/SST use frozen DIV90 population means/std and [-3,3] gene clipping. Context panels use DIV30 gene z scores, unweighted, clipped [-3,3]. No context genes enter early-PV target.',rank_scoring='Within-cell reverse expression ranks among all 18082 genes; ties average; undetected genes zero; truncated at rank1501; mean(1501-min(rank,1501))/1500.',reference='4768 frozen DIV90 cortical LHX6+/ERBB4+ cells; cannot represent all DIV30 lineages.',consensus='Mean percentiles of four methods, with ties kept; this descriptive rank exists even if all absolute signals weak. >=3/4 methods in top20% is separately labeled multi-method support; direct and rank are mathematically related.',target_rule='Top q uses score >= quantile(1-q); ties are kept and actual target prevalence is reported.',mapping='Reference-only variable nonsurface/non-scoring genes, reference StandardScaler+30 PCs; 30 nearest reference cells, no cluster target.',sparse_model='ElasticNet alpha=0.02,l1_ratio=0.9,max_iter=5000; no hyperparameter tuning. Training labels only frozen DIV90 PV scores.',validation='Leave one DIV90 sample or cell line out. Standardization/PCA fitted only on training cells. Feature selection is unsupervised reference-wide for primary mapping; held-out model validation additionally reselects variable genes on training cells.',developmental_order='Contextual continuous scores and quantile bins only; no new trajectory or chronological lineage ordering.',sources=SOURCES))
    path=OUT/'cache/expression.npy';rankpath=OUT/'cache/rank_scores.npy'
    if not path.exists() or not rankpath.exists():
        expr=np.lib.format.open_memmap(path,mode='w+',dtype='float32',shape=(len(obs),len(chosen)))
        ranks=np.zeros((len(obs),2),dtype='float32')
        moduleix=[pd.Index(allgenes).get_indexer(gs) for gs in [PV,SST]]
        for start,stop,x in iter_x(H30):
            expr[start:stop]=x[:,cols].toarray()
            x.sort_indices();rd=np.zeros(x.nnz,dtype='float32')
            for r in range(x.shape[0]):
                aa,bb=x.indptr[r:r+2];rd[aa:bb]=np.maximum(1501-stats.rankdata(-x.data[aa:bb],method='average'),0)/1500
            rankx=sparse.csr_matrix((rd,x.indices,x.indptr),shape=x.shape)
            for j,mi in enumerate(moduleix):ranks[start:stop,j]=np.asarray(rankx[:,mi].mean(axis=1)).ravel()
            if start%12000==0:print(f'Extracted DIV30 {stop:,}/{len(obs):,}',flush=True)
        expr.flush();np.save(rankpath,ranks)
    expr=np.load(path,mmap_mode='r');lookup={g:i for i,g in enumerate(chosen)}
    cells=pd.DataFrame({'cell_id':obs.cell_id,'sample':obs['orig.ident']})
    cells=attach_metadata(cells)
    for source,target in [('nCount_RNA','total_counts'),('nFeature_RNA','n_genes'),('percent.mt','percent_mt'),('S.Score','s_score'),('G2M.Score','g2m_score')]:cells[target]=pd.to_numeric(obs[source],errors='raise')
    cells['log_counts']=np.log1p(cells.total_counts)
    cells['cell_cycle_phase']=obs.Phase
    for g in sorted(set(CONTEXT+EXCLUDED)):cells[g]=expr[:,lookup[g]] if g in lookup else 0.
    pars=json.loads((REF/'state_parameters.json').read_text())['modules']
    for kind,gs in [('pv',PV),('sst',SST)]:
        z=np.column_stack([np.clip((cells[g]-pars[kind]['means'][g])/pars[kind]['std_population'][g],-3,3) for g in gs])
        cells[f'direct_{kind}_score']=z.mean(axis=1)
        cells[f'{kind}_n_detected']=(cells[gs]>0).sum(axis=1)
        for j,g in enumerate(gs):cells[f'z90_{g}']=z[:,j]
    refpars={}
    for name,gs in PANELS.items():
        x=cells[gs].to_numpy();mu=x.mean(axis=0);sd=x.std(axis=0);sd[sd==0]=1
        cells[name+'_score']=np.clip((x-mu)/sd,-3,3).mean(axis=1)
        refpars[name]={'genes':gs,'means':dict(zip(gs,mu)),'std_population':dict(zip(gs,sd))}
    dump(OUT/'provenance/context_score_parameters.json',refpars)
    ranks=np.load(rankpath)
    cells['rank_pv_score']=ranks[:,0];cells['rank_sst_score']=ranks[:,1]
    cells['cycling']=cells.cell_cycle_phase.isin(['S','G2M'])
    cells['progenitor']=(cells.progenitor_score>=cells.progenitor_score.median())&(cells.progenitor_score>cells.maturation_score)
    cells['postmitotic']=(~cells.cycling)&(cells.maturation_score>cells.progenitor_score)
    # Histories / identities remain metadata only, never model labels.
    cells.to_csv(OUT/'cache/prepared_cells.tsv.gz',sep='\t',index=False)
    return cells,expr,chosen,ref,refx,allgenes,featureix


def fit_mapping(train,query,ypv,ysst,strong):
    scaler=StandardScaler().fit(train)
    a=np.clip(scaler.transform(train),-5,5)
    pca=PCA(n_components=30,svd_solver='randomized',random_state=30).fit(a)
    rp=pca.transform(a)
    knn=NearestNeighbors(n_neighbors=min(30,len(train)),n_jobs=4).fit(rp)
    distances=np.empty((len(query),30),dtype='float32');neighbors=np.empty((len(query),30),dtype='int32')
    for start in range(0,len(query),4000):
        stop=min(start+4000,len(query));qp=pca.transform(np.clip(scaler.transform(query[start:stop]),-5,5))
        d,ii=knn.kneighbors(qp);distances[start:stop]=d;neighbors[start:stop]=ii
    return distances,neighbors,ypv[neighbors].mean(axis=1),ysst[neighbors].mean(axis=1),strong[neighbors].mean(axis=1)


def validate_reference(ref,refx,allgenes):
    path=OUT/'tables/reference_validation.tsv'
    if path.exists():return pd.read_csv(path,sep='\t'),pd.read_csv(OUT/'tables/reference_oof_predictions.tsv.gz',sep='\t')
    y=ref.pv_score.to_numpy();ss=ref.sst_score.to_numpy();strong=y>=np.quantile(y,.75)
    rows=[];predictions=[]
    for grouping in ['sample','cell_line']:
        for group in sorted(ref[grouping].unique()):
            train=ref[grouping].to_numpy()!=group;test=~train
            _,fi=choose_features(refx[train],allgenes)
            a=np.asarray(refx[train][:,fi]);b=np.asarray(refx[test][:,fi])
            scaler=StandardScaler().fit(a);xa=scaler.transform(a);xb=scaler.transform(b)
            model=ElasticNet(alpha=.02,l1_ratio=.9,max_iter=5000,random_state=30).fit(xa,y[train]);pred=model.predict(xb)
            dist,nn,nscore,_,_=fit_mapping(a,b,y[train],ss[train],strong[train])
            for name,p in [('sparse',pred),('neighbor',nscore)]:
                rows.append(dict(grouping=grouping,held_out=group,n_train=int(train.sum()),n_test=int(test.sum()),method=name,spearman=stats.spearmanr(p,y[test]).statistic,pearson=stats.pearsonr(p,y[test]).statistic,r2=r2_score(y[test],p),mean_prediction=float(p.mean()),mean_truth=float(y[test].mean()),nonzero_coefficients=int(np.count_nonzero(model.coef_)) if name=='sparse' else None))
            predictions.append(pd.DataFrame({'cell_id':ref.cell_id[test],'grouping':grouping,'held_out':group,'observed_pv_score':y[test],'sparse_prediction':pred,'neighbor_prediction':nscore,'neighbor_distance':dist.mean(axis=1)}))
            print(f'Validated DIV90 {grouping}={group}: n={test.sum()}',flush=True)
    result=pd.DataFrame(rows);oof=pd.concat(predictions,ignore_index=True)
    tab(result,'reference_validation.tsv');tab(oof,'reference_oof_predictions.tsv.gz')
    return result,oof


def score_map(cells,expr,genes,ref,refx,allgenes,featureix):
    validation,oof=validate_reference(ref,refx,allgenes)
    features=[allgenes[i] for i in featureix];fi=pd.Index(genes).get_indexer(features)
    a=np.asarray(refx[:,featureix]);b=np.asarray(expr[:,fi]);y=ref.pv_score.to_numpy();s=ref.sst_score.to_numpy();strong=y>=np.quantile(y,.75)
    scaler=StandardScaler().fit(a);model=ElasticNet(alpha=.02,l1_ratio=.9,max_iter=5000,random_state=30).fit(scaler.transform(a),y)
    cells['sparse_pv_score']=model.predict(scaler.transform(b))
    coef=pd.DataFrame({'gene':features,'coefficient':model.coef_,'reference_mean':scaler.mean_,'reference_scale':scaler.scale_})
    tab(coef,'sparse_model_coefficients.tsv')
    dist,nn,pv,ss,fr=fit_mapping(a,b,y,s,strong)
    cells['neighbor_pv_score']=pv;cells['neighbor_sst_score']=ss;cells['neighbor_strong_pv_fraction']=fr
    cells['mapping_distance']=dist.mean(axis=1)
    base=oof.query("grouping == 'sample'").neighbor_distance
    cells['mapping_confidence']=np.exp(-cells.mapping_distance/base.median())
    cells['mapping_out_of_reference']=cells.mapping_distance>base.quantile(.95)
    np.savez_compressed(OUT/'cache/reference_neighbors.npz',reference_indices=nn,distances=dist)
    ref[['cell_id','sample','cell_line','pv_score','sst_score','state']].to_csv(OUT/'tables/reference_cell_index.tsv',sep='\t',index=False)
    with h5py.File(OUT/'cache/mapping_model.h5','w') as f:
        f.create_dataset('features',data=np.array(features,dtype=h5py.string_dtype()))
        f.create_dataset('sparse_coefficient',data=model.coef_);f.attrs['sparse_intercept']=model.intercept_
        f.create_dataset('feature_means',data=scaler.mean_);f.create_dataset('feature_scales',data=scaler.scale_)
    agreement=[]
    for i,aa in enumerate(METHODS):
        for bb in METHODS[i+1:]:
            yy=cells[aa]>=cells[aa].quantile(.8);zz=cells[bb]>=cells[bb].quantile(.8)
            agreement.append(dict(method_a=aa,method_b=bb,spearman=stats.spearmanr(cells[aa],cells[bb]).statistic,top20_jaccard=(yy&zz).sum()/(yy|zz).sum(),a_n=int(yy.sum()),b_n=int(zz.sum()),intersection_n=int((yy&zz).sum())))
    tab(pd.DataFrame(agreement),'score_agreement.tsv')
    # Compute only after the independent agreement table has been written.
    cells['early_pv_score']=cells[METHODS].rank(pct=True,method='average').mean(axis=1)
    votes=np.column_stack([cells[m]>=cells[m].quantile(.8) for m in METHODS]).sum(axis=1)
    cells['top20_method_votes']=votes;cells['multi_method_supported']=votes>=3
    absrows=[]
    for kind in ['pv','sst']:
        scores=cells['direct_'+kind+'_score'];rv=ref[kind+'_score']
        for label,cut in [('DIV90_median',rv.median()),('DIV90_top25',rv.quantile(.75)),('DIV90_top10',rv.quantile(.9))]:
            absrows.append(dict(program=kind,reference_cut=label,threshold=cut,n_div30=int((scores>=cut).sum()),fraction_div30=float((scores>=cut).mean()),div30_median=scores.median(),div90_median=rv.median(),two_genes_detected_fraction=(cells[kind+'_n_detected']>=2).mean()))
    tab(pd.DataFrame(absrows),'absolute_reference_overlap.tsv')
    component=[]
    for program,gs in [('pv',PV),('sst',SST)]:
        for g in gs:component.append(dict(program=program,gene=g,n_detected=int((cells[g]>0).sum()),fraction_detected=(cells[g]>0).mean(),mean_expression=cells[g].mean(),reference_fraction_detected=(ref[g]>0).mean(),reference_mean_expression=ref[g].mean()))
    tab(pd.DataFrame(component),'module_detection.tsv')
    sensitivities=[]
    variants={'without_MEF2C':[g for g in PV if g!='MEF2C'],'without_KCNC':[g for g in PV if not g.startswith('KCNC')],'without_MAF_MAFB':PV,'broad_with_MAF_MAFB':PV+['MAF','MAFB']}
    pars=json.loads((REF/'state_parameters.json').read_text())['modules']['pv']
    for name,gs in variants.items():
        cols=[]
        for g in gs:
            mu=pars['means'].get(g,ref[g].mean());sd=pars['std_population'].get(g,ref[g].std(ddof=0));cols.append(np.clip((cells[g]-mu)/sd,-3,3))
        v=np.column_stack(cols).mean(axis=1);cells['pv_'+name]=v
        # Consensus with the direct component removed/replaced diagnoses its influence;
        # the mapping/supervised labels remain the frozen reference, explicitly.
        alt=(stats.rankdata(v,method='average')/len(v)+cells[METHODS[1:]].rank(pct=True).sum(axis=1))/4
        original=cells.early_pv_score>=cells.early_pv_score.quantile(.8);alternate=alt>=np.quantile(alt,.8)
        sensitivity=dict(variant=name,genes=';'.join(gs),direct_spearman=stats.spearmanr(v,cells.direct_pv_score).statistic,consensus_spearman=stats.spearmanr(alt,cells.early_pv_score).statistic,top20_jaccard=(original&alternate).sum()/(original|alternate).sum(),unchanged_primary=name=='without_MAF_MAFB',interpretation='MAF/MAFB absent from primary by design' if name=='without_MAF_MAFB' else 'Direct-transfer sensitivity; frozen mapping and supervised target unchanged')
        sensitivities.append(sensitivity)
    sst_without=cells[[f'z90_{g}' for g in SST if g!='SST']].mean(axis=1)
    sstrho=stats.spearmanr(sst_without,cells.direct_sst_score).statistic
    sst_stable=bool(sstrho>=.7 and (cells.sst_n_detected>=2).mean()>=.1)
    sensitivities.append(dict(variant='SST_without_SST',genes=';'.join(g for g in SST if g!='SST'),direct_spearman=sstrho,interpretation='SST stable only if leave-SST-out rho>=0.7 AND >=10% cells detect >=2 SST genes; analyst diagnostic, not validated subtype identity.'))
    tab(pd.DataFrame(sensitivities),'score_sensitivity.tsv')
    dump(OUT/'provenance/sst_support.json',dict(stable_for_low_sst_target=sst_stable,leave_sst_out_spearman=sstrho,fraction_at_least_two_genes=(cells.sst_n_detected>=2).mean()))
    cells.to_csv(OUT/'cells.tsv.gz',sep='\t',index=False)
    return cells,sst_stable


def development(cells):
    rows=[];group_rows=[]
    for q in [.1,.2,.3]:
        selected=cells.early_pv_score>=cells.early_pv_score.quantile(1-q)
        for marker in ['ERBB4','LHX6','NKX2-1','MAF','MAFB','MEF2C','DCX']:
            positive=cells[marker]>0
            rows.append(dict(target=f'top{round(q*100)}',nominal_fraction=q,actual_target_fraction=selected.mean(),target_n=int(selected.sum()),marker=marker,positive_n=int((selected&positive).sum()),positive_fraction=float(positive[selected].mean()),lost_n=int((selected&~positive).sum()),lost_fraction=float((~positive[selected]).mean()),progenitor_fraction=cells.progenitor[selected].mean(),cycling_fraction=cells.cycling[selected].mean(),postmitotic_fraction=cells.postmitotic[selected].mean(),multi_method_supported_fraction=cells.multi_method_supported[selected].mean(),mapping_out_of_reference_fraction=cells.mapping_out_of_reference[selected].mean()))
        for axis in ['sample','cell_line','condition']:
            for label,group in cells.groupby(axis):
                ii=group.index;tg=selected.loc[ii]
                group_rows.append(dict(target=f'top{round(q*100)}',axis=axis,group=label,n_cells=len(ii),target_n=int(tg.sum()),target_fraction=tg.mean(),mean_early_pv_score=group.early_pv_score.mean(),erbb4_loss_fraction=float((group.loc[tg,'ERBB4']<=0).mean()),lhx6_loss_fraction=float((group.loc[tg,'LHX6']<=0).mean())))
    tab(pd.DataFrame(rows),'entry_marker_loss.tsv');tab(pd.DataFrame(group_rows),'target_by_sample_line_condition.tsv')
    bins=[]
    axes=['progenitor_score','mge_score','NKX2-1','LHX6','ERBB4','MAF','MAFB','MEF2C','DCX','maturation_score','progression_score','cycle_score']
    for ax in axes:
        # Quantile bins preserve ties; zeros are not split into artificial bins.
        cut=pd.qcut(cells[ax],20,duplicates='drop')
        for label,gg in cells.groupby(cut,observed=True):
            bins.append(dict(axis=ax,bin=str(label),n=len(gg),x_mean=gg[ax].mean(),x_min=gg[ax].min(),x_max=gg[ax].max(),early_pv_mean=gg.early_pv_score.mean(),early_pv_median=gg.early_pv_score.median(),early_pv_q25=gg.early_pv_score.quantile(.25),early_pv_q75=gg.early_pv_score.quantile(.75),direct_pv_mean=gg.direct_pv_score.mean(),neighbor_pv_mean=gg.neighbor_pv_score.mean(),top20_fraction=(gg.early_pv_score>=cells.early_pv_score.quantile(.8)).mean(),erbb4_positive=(gg.ERBB4>0).mean(),lhx6_positive=(gg.LHX6>0).mean(),cycling_fraction=gg.cycling.mean()))
    tab(pd.DataFrame(bins),'developmental_bins.tsv')
    masks={
      'MGE_progenitor':cells.progenitor&(cells.mge_score>=cells.mge_score.median()),
      'cycling':cells.cycling,
      'postmitotic_ERBB4_negative':cells.postmitotic&(cells.ERBB4<=0),
      'LHX6_positive_ERBB4_negative':(cells.LHX6>0)&(cells.ERBB4<=0),
      'LHX6_positive_ERBB4_positive':(cells.LHX6>0)&(cells.ERBB4>0),
      'higher_maturation':cells.maturation_score>=cells.maturation_score.quantile(.75),
    }
    stages=[]
    for name,m in masks.items():
        stages.append(dict(context=name,n=int(m.sum()),population_fraction=m.mean(),mean_early_pv_score=cells.loc[m,'early_pv_score'].mean(),mean_direct_pv_score=cells.loc[m,'direct_pv_score'].mean(),mean_neighbor_pv_score=cells.loc[m,'neighbor_pv_score'].mean(),top20_fraction=(cells.loc[m,'early_pv_score']>=cells.early_pv_score.quantile(.8)).mean(),supported_fraction=cells.loc[m,'multi_method_supported'].mean()))
    tab(pd.DataFrame(stages),'developmental_contexts.tsv')
    correlations=[]
    for score in METHODS+['early_pv_score']:
        for ax in axes+['log_counts','n_genes']:
            correlations.append(dict(score=score,axis=ax,spearman=stats.spearmanr(cells[score],cells[ax]).statistic))
    tab(pd.DataFrame(correlations),'score_developmental_correlations.tsv')


def sensitivity_end_to_end(cells):
    """Drop a gene family from all four scoring methods, retaining primary fits."""
    ref,allgenes,refx=reference_matrix()
    rankgenes=PV+SST+['MAF','MAFB']
    path=OUT/'cache/per_gene_rank_scores.npy'
    if not path.exists() or not (OUT/'cache/per_gene_rank_genes.json').exists():
        geneix=pd.Index(allgenes).get_indexer(rankgenes)
        ranked=np.lib.format.open_memmap(path,mode='w+',dtype='float32',shape=(len(cells),len(rankgenes)))
        for start,stop,x in iter_x(H30):
            x.sort_indices();rd=np.empty(x.nnz,dtype='float32')
            for r in range(x.shape[0]):
                aa,bb=x.indptr[r:r+2];rd[aa:bb]=np.maximum(1501-stats.rankdata(-x.data[aa:bb],method='average'),0)/1500
            rankx=sparse.csr_matrix((rd,x.indices,x.indptr),shape=x.shape)
            ranked[start:stop]=rankx[:,geneix].toarray()
            if start%24000==0:print(f'Per-gene rank sensitivity {stop:,}/{len(cells):,}',flush=True)
        ranked.flush();dump(OUT/'cache/per_gene_rank_genes.json',rankgenes)
    ranked=np.load(path,mmap_mode='r')
    assert np.allclose(ranked[:,:4].mean(axis=1),cells.rank_pv_score,atol=1e-6)
    cfg=json.loads((OUT/'provenance/scoring_config.json').read_text())
    expr=np.load(OUT/'cache/expression.npy',mmap_mode='r');genes=json.loads((OUT/'cache/genes.json').read_text())
    features=cfg['mapping_features'];ia=pd.Index(allgenes).get_indexer(features);ib=pd.Index(genes).get_indexer(features)
    a=np.asarray(refx[:,ia]);b=np.asarray(expr[:,ib]);scaler=StandardScaler().fit(a)
    xa=scaler.transform(a);xb=scaler.transform(b)
    nn=np.load(OUT/'cache/reference_neighbors.npz')['reference_indices']
    pars=json.loads((REF/'state_parameters.json').read_text())['modules']['pv']
    variants={'without_MEF2C':[g for g in PV if g!='MEF2C'],'without_KCNC':[g for g in PV if not g.startswith('KCNC')],'without_MAF_MAFB':PV,'broad_with_MAF_MAFB':PV+['MAF','MAFB']}
    primary=cells.early_pv_score.to_numpy();original=primary>=np.quantile(primary,.8);rows=[];scores={}
    percell=cells[['cell_id']].copy()
    for name,gs in variants.items():
        zz=[];rr=[]
        for g in gs:
            mu=pars['means'].get(g,ref[g].mean());sd=pars['std_population'].get(g,ref[g].std(ddof=0))
            zz.append(np.clip((cells[g].to_numpy()-mu)/sd,-3,3));rr.append(np.clip((ref[g].to_numpy()-mu)/sd,-3,3))
        direct=np.column_stack(zz).mean(axis=1);reference=np.column_stack(rr).mean(axis=1)
        neighbor=reference[nn].mean(axis=1)
        model=ElasticNet(alpha=.02,l1_ratio=.9,max_iter=5000,random_state=30).fit(xa,reference)
        predicted=model.predict(xb)
        rankscore=ranked[:,[rankgenes.index(g) for g in gs]].mean(axis=1)
        components=np.column_stack([direct,neighbor,predicted,rankscore])
        consensus=stats.rankdata(components,axis=0,method='average').mean(axis=1)/len(cells)
        selected=consensus>=np.quantile(consensus,.8)
        rows.append(dict(variant=name,genes=';'.join(gs),consensus_spearman=stats.spearmanr(primary,consensus).statistic,top20_jaccard=(selected&original).sum()/(selected|original).sum(),original_target_n=int(original.sum()),alternative_target_n=int(selected.sum()),original_target_recovery=float(selected[original].mean()),direct_spearman=stats.spearmanr(cells.direct_pv_score,direct).statistic,neighbor_spearman=stats.spearmanr(cells.neighbor_pv_score,neighbor).statistic,sparse_spearman=stats.spearmanr(cells.sparse_pv_score,predicted).statistic,rank_spearman=stats.spearmanr(cells.rank_pv_score,rankscore).statistic,interpretation='All four components recalculated: direct transfer, reference-neighbor labels, sparse supervised labels, and within-cell ranks. Reference neighbor identities/features remain fixed independently of scoring genes; primary analysis unchanged.'))
        scores[name]=consensus
        for j,method in enumerate(METHODS):percell[name+'__'+method]=components[:,j]
        percell[name+'__consensus']=consensus
        print(f'End-to-end sensitivity {name}: rho={rows[-1]["consensus_spearman"]:.3f}, top20 Jaccard={rows[-1]["top20_jaccard"]:.3f}',flush=True)
    tab(pd.DataFrame(rows),'end_to_end_score_sensitivity.tsv')
    tab(percell,'end_to_end_cell_score_sensitivity.tsv.gz')
    np.savez_compressed(OUT/'cache/end_to_end_consensus_sensitivity.npz',**scores)
    return scores


def sensitivity_fixed_gates(cells):
    path=OUT/'tables/surface_gate_best_at_recovery_floors.tsv'
    if not path.exists() or not (OUT/'cache/end_to_end_consensus_sensitivity.npz').exists():return
    from div30_antecedent_surface import apply,metric
    best=pd.read_csv(path,sep='\t').query("target == 'top20'").drop_duplicates('gate_id')
    expr=np.load(OUT/'cache/expression.npy',mmap_mode='r');genes=json.loads((OUT/'cache/genes.json').read_text());lookup={g:i for i,g in enumerate(genes)}
    alternate=np.load(OUT/'cache/end_to_end_consensus_sensitivity.npz');rows=[]
    for gate in best.to_dict('records'):
        keep=apply(gate,expr,lookup)
        for name in alternate.files:
            score=alternate[name];target=score>=np.quantile(score,.8)
            row={'gate_id':gate['gate_id'],'gate_label':gate['gate_label'],'variant':name,'evaluation':'Fixed primary top20-selected gate; alternate target from end-to-end gene-drop sensitivity; no gate refitting'}
            row.update({k:float(v) for k,v in metric(int(target[keep].sum()),int(keep.sum()),int(target.sum()),len(target)).items()});rows.append(row)
    tab(pd.DataFrame(rows),'fixed_gate_end_to_end_sensitivity.tsv')
    # Generic neuronal/MGE/progenitor composition and independent SST association
    # keep surface enrichment of a rank target from being mistaken for subtype fate.
    composition=[]
    for gate in best.to_dict('records'):
        keep=apply(gate,expr,lookup)
        for fraction,mask in [('all_DIV30',np.ones(len(cells),bool)),('retained',keep),('removed',~keep)]:
            d=cells.loc[mask]
            row={'gate_id':gate['gate_id'],'gate_label':gate['gate_label'],'fraction':fraction,'n_cells':int(mask.sum())}
            for field in ['direct_pv_score','direct_sst_score','neighbor_pv_score','neighbor_sst_score','early_pv_score','mge_score','progenitor_score','maturation_score','progression_score','log_counts']:
                row['mean_'+field]=d[field].mean()
            for field in ['progenitor','postmitotic','cycling','multi_method_supported','mapping_out_of_reference']:
                row[field+'_fraction']=d[field].mean()
            for marker in ['NKX2-1','LHX6','ERBB4','GAD1','GAD2','DCX','PVALB']:
                row[marker+'_detected_fraction']=(d[marker]>0).mean()
            composition.append(row)
    tab(pd.DataFrame(composition),'fixed_gate_biological_composition.tsv')


def figures(cells):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scipy.ndimage import gaussian_filter1d
    plt.rcParams.update({'font.size':10,'axes.spines.right':False,'axes.spines.top':False,'svg.fonttype':'none','pdf.fonttype':42})
    def save(fig,name):
        for typ in ['png','pdf','svg']:fig.savefig(OUT/f'figures/{typ}/{name}.{typ}',dpi=300,bbox_inches='tight')
        plt.close(fig)
    ref=pd.read_csv(REF/'cells.tsv.gz',sep='\t')
    fig,axs=plt.subplots(2,3,figsize=(15,8),layout='constrained')
    corr=cells[METHODS].corr(method='spearman');im=axs[0,0].imshow(corr,vmin=-1,vmax=1,cmap='coolwarm')
    labs=['Transfer','Neighbors','Sparse model','Within-cell rank']
    axs[0,0].set_xticks(range(4),labs,rotation=45,ha='right');axs[0,0].set_yticks(range(4),labs)
    for i in range(4):
        for j in range(4):axs[0,0].text(j,i,f'{corr.iloc[i,j]:.2f}',ha='center',va='center')
    axs[0,0].set_title('Methods agree only partially')
    axs[0,1].hist(cells.direct_pv_score,bins=65,density=True,alpha=.6,label='DIV30')
    axs[0,1].hist(ref.pv_score,bins=65,density=True,alpha=.6,label='Frozen DIV90')
    axs[0,1].set(xlabel='PV score on fixed DIV90 scale',ylabel='Density',title='Absolute transferred program');axs[0,1].legend()
    axs[0,2].hexbin(cells.direct_pv_score,cells.neighbor_pv_score,gridsize=55,bins='log',mincnt=1,cmap='viridis')
    axs[0,2].set(xlabel='Direct PV transfer',ylabel='Mean reference-neighbor PV',title='Independent feature mapping')
    axs[1,0].hexbin(cells.maturation_score,cells.early_pv_score,gridsize=55,bins='log',mincnt=1,cmap='viridis')
    axs[1,0].set(xlabel='Neuronal maturation context',ylabel='Early-PV descriptive consensus',title='Generic maturation confound')
    vv=pd.read_csv(OUT/'tables/reference_validation.tsv',sep='\t').query("grouping=='cell_line'")
    for j,m in enumerate(['neighbor','sparse']):
        d=vv[vv.method==m];axs[1,1].bar(np.arange(len(d))+(j-.5)*.35,d.spearman,width=.35,label=m)
    axs[1,1].set_xticks(np.arange(3),sorted(ref.cell_line.unique()));axs[1,1].set(ylabel='Held-line Spearman rho',title='Reference validation');axs[1,1].legend()
    counts=cells.top20_method_votes.value_counts().sort_index();axs[1,2].bar(counts.index,counts.values,color='#46758b')
    axs[1,2].set(xlabel='Methods calling cell in top 20%',ylabel='DIV30 cells',title='Show disagreement and supported subset')
    save(fig,'01_program_transfer_and_agreement')
    fig,axs=plt.subplots(3,4,figsize=(16,10),layout='constrained')
    b=pd.read_csv(OUT/'tables/developmental_bins.tsv',sep='\t')
    for ax,name in zip(axs.flat,b.axis.unique()):
        d=b[b.axis==name].sort_values('x_mean');ax.fill_between(d.x_mean,d.early_pv_q25,d.early_pv_q75,color='#87a6b6',alpha=.3)
        ax.plot(d.x_mean,d.early_pv_mean,'o',ms=3,color='#255a73')
        if len(d)>2:ax.plot(d.x_mean,gaussian_filter1d(d.early_pv_mean.to_numpy(),.9),color='#255a73')
        ax.set(xlabel=name.replace('_',' '),ylabel='Early-PV rank consensus')
    fig.suptitle('Cross-sectional relationships; quantile bins and smoothed means, not lineage trajectories')
    save(fig,'02_developmental_context_curves')
    losses=pd.read_csv(OUT/'tables/entry_marker_loss.tsv',sep='\t')
    fig,axs=plt.subplots(1,3,figsize=(15,4.2),layout='constrained')
    for ax,target in zip(axs,['top10','top20','top30']):
        d=losses[losses.target==target];ax.bar(d.marker,d.positive_fraction,color=['#a24d58' if m in ['ERBB4','LHX6'] else '#46758b' for m in d.marker]);ax.set(ylim=(0,1),ylabel='Fraction of target with detected RNA',title=f'{target}: n={int(d.target_n.iloc[0]):,}');ax.tick_params(axis='x',rotation=45)
    save(fig,'03_early_candidate_entry_marker_loss')
    (OUT/'figures/FIGURE_CAPTIONS.md').write_text('''# DIV30 figure captions\n\n1. Four complementary PV-related scores and their pairwise correlations, absolute DIV90-reference overlap, cross-feature mapping, generic maturation dependence, held-line reference prediction, and method-vote counts. Direct and within-cell rank scoring share genes; agreement does not make them independent biological replicates.\n2. All DIV30 cells, ordered separately along each displayed contextual expression axis. Quantile bins keep ties, points are means, ribbons are cell-distribution interquartiles, lines smooth adjacent bin means. These cross-sectional curves do not establish chronological onset or lineage. Maturation includes DCX; progression is generic inhibitory differentiation, not proof of cortical fate.\n3. Raw RNA detection (>0) in top10/20/30% rank-consensus sets, keeping ties. ERBB4 or LHX6 non-detection quantifies the hypothetical RNA-entry loss, not validated surface protein negativity.\n''')


def report(cells):
    agreement=pd.read_csv(OUT/'tables/score_agreement.tsv',sep='\t')
    loss=pd.read_csv(OUT/'tables/entry_marker_loss.tsv',sep='\t')
    absolute=pd.read_csv(OUT/'tables/absolute_reference_overlap.tsv',sep='\t')
    detection=pd.read_csv(OUT/'tables/module_detection.tsv',sep='\t')
    context=pd.read_csv(OUT/'tables/developmental_contexts.tsv',sep='\t')
    correlations=pd.read_csv(OUT/'tables/score_developmental_correlations.tsv',sep='\t')
    sensitivity=pd.read_csv(OUT/'tables/score_sensitivity.tsv',sep='\t')
    fullsens_path=OUT/'tables/end_to_end_score_sensitivity.tsv'
    fullsens=pd.read_csv(fullsens_path,sep='\t') if fullsens_path.exists() else pd.DataFrame()
    summary={'n_cells':len(cells),'n_samples':cells['sample'].nunique(),'n_pvalb_detected':int((cells.PVALB>0).sum()),'absolute_reference_overlap':absolute.to_dict('records'),'score_agreement':agreement.to_dict('records'),'top20_entry_marker_loss':loss.query("target=='top20'").to_dict('records'),'multi_method_supported_n':int(cells.multi_method_supported.sum()),'multi_method_supported_fraction':cells.multi_method_supported.mean(),'mapping_out_of_reference_fraction':cells.mapping_out_of_reference.mean(),'early_pv_maturation_spearman':stats.spearmanr(cells.early_pv_score,cells.maturation_score).statistic,'sensitivity':sensitivity.to_dict('records'),'classification':'transcriptional antecedent / developmental resemblance; no confirmed future PV cells','program_reference':'frozen independent DIV90 PV and SST axes; all 90631 DIV30 profiles scored separately','protein_endpoint':'Later PV protein production is user-supplied system-level endpoint; not individually linked to these profiled cells.'}
    summary['end_to_end_sensitivity']=fullsens.to_dict('records')
    surfacepath=OUT/'surface_summary.json'
    if surfacepath.exists():summary['surface']=json.loads(surfacepath.read_text())
    d=absolute.query("program=='pv' and reference_cut=='DIV90_median'").iloc[0]
    er=loss.query("target=='top20' and marker=='ERBB4'").iloc[0];lh=loss.query("target=='top20' and marker=='LHX6'").iloc[0]
    rho=summary['early_pv_maturation_spearman'];sens=sensitivity.set_index('variant')
    gateanswers={}
    if 'surface' in summary:
        sf=summary['surface'];best=pd.DataFrame(sf['best_at_recovery_floors']);top=best.query("target=='top20'")
        primary=top[np.isclose(top.recovery_floor,.5)&(top.maximum_markers==3)].iloc[0]
        high=top[np.isclose(top.recovery_floor,.8)&(top.maximum_markers==3)].iloc[0]
        benchmark=pd.DataFrame(sf['frozen_DIV90_benchmark']).query("target=='top20'").iloc[0]
        val=pd.DataFrame(sf['held_sample_validation']);held=val[np.isclose(val.recovery_floor,.5)].iloc[0]
        markers=pd.read_csv(OUT/'tables/surface_marker_confound_models.tsv',sep='\t').set_index('gene')
        chosen_genes=[g['gene'] for g in json.loads(primary.rules)]
        classes='; '.join(f'{g}: {markers.loc[g,"classification"]}' for g in chosen_genes)
        candidate=(' '+primary.logic+' ').join(f'{g["gene"]} {g["op"]} {g["threshold"]:.17g}' for g in json.loads(primary.rules))
        summary['primary_top20_gate_at_50pct_recovery']=primary.to_dict()
        summary['primary_gate_marker_classifications']={g:markers.loc[g,'classification'] for g in chosen_genes}
        summary['frozen_div90_benchmark_top20']=benchmark.to_dict()
        gateanswers[6]=f'6. **Does FGFR2/PTPRS contain useful information?** Applying the exact frozen DIV90 rule changes top20 antecedent composition **{benchmark.starting_target_fraction:.1%} → {benchmark.post_gate_target_fraction:.1%}**, **{benchmark.fold_enrichment:.2f}× enrichment**, recovering **{benchmark.target_recovery:.1%}** of targets at **{benchmark.total_cell_yield:.1%}** total yield. This is **{benchmark.enrichment_interpretation}** by the requested descriptive bins, and is a diagnostic RNA-threshold transfer rather than replication.'
        gateanswers[7]=f'7. **Is there a stronger DIV30 surface gate?** Yes in numerical, in-sample score enrichment: the best evaluated <=3-marker top20 rule at >=50% recovery is **retain `{candidate}`**, giving **{primary.starting_target_fraction:.1%} → {primary.post_gate_target_fraction:.1%}** target composition (**{primary.fold_enrichment:.2f}×**), **{primary.target_recovery:.1%}** target recovery, and **{primary.total_cell_yield:.1%}** total cell yield. This enriches a candidate transcriptional state; PV fate specificity remains unestablished.'
        gateanswers[8]=f'8. **Are best markers specific or maturation/progenitor markers?** For this selected gate: **{classes}**. Classes are based on the declared effect-size heuristics and regressions adjusting measured progenitor, maturation, MGE and depth axes, with separate neighbor and sample checks. Residual association can nominate a subtype candidate, but does not establish causal or prospective PV specificity.'
        gateanswers[9]=f'9. **Strongest achievable enrichment?** At >=50% target recovery the evaluated discovery optimum is **{primary.post_gate_target_fraction:.1%} purity/{primary.fold_enrichment:.2f}× enrichment**; at >=80% recovery it is **{high.post_gate_target_fraction:.1%}/{high.fold_enrichment:.2f}×**, with actual recovery **{high.target_recovery:.1%}**. With complete training-sample gate reselection, pooled held-out top20 performance at the 50% training floor is **{held.post_gate_target_fraction:.1%} purity, {held.fold_enrichment:.2f}×, {held.target_recovery:.1%} actual recovery**. Different gates are chosen across folds. Lower 10/20/30% floors and unconstrained maxima are separately reported to expose yield tradeoffs; no global or protein-level optimum is proven.'
        gateanswers[10]=f'10. **Is DIV30 worth sorting experimentally?** The RNA-derived state is numerically enrichable, but this is not yet a validated PV-directed purification. The primary candidate spans these classes: **{classes}**. Together with maturation correlation **rho={rho:.2f}**, partial method agreement, and gene-drop membership changes, this supports at most a controlled exploratory fraction-comparison experiment measuring immediate neuronal/progenitor composition and later PV protein, not a routine PV-production sort or a predicted PV fate label.'
        gateanswers[11]=f'11. **If tested, which gate?** The concrete top20, >=50%-recovery RNA hypothesis is **retain `{candidate}`**. Preserve retained, removed and unsorted comparators and evaluate later PV protein under matched culture. The rule is a hypothesis for surface phenotype testing; these log-RNA thresholds are not fluorescence cutoffs and no antibody or live-cell protocol is established here.'
    paragraphs=[
      '# DIV30: candidate transcriptional antecedents of the later PV-associated state',
      f'All **{len(cells):,}** cached DIV30 cells were analyzed separately from the **4,768** frozen DIV90 reference cells. No expression pooling, new UMAP, cluster targets, trajectory fit, or original-result changes were used. PVALB RNA is detected in **{summary["n_pvalb_detected"]:,}** cells. The user reports that this MGEO system later produces bona fide PV protein; cross-sectional resemblance cannot identify a confirmed future PV cell.',
      '## Twelve requested answers',
      f'1. **Is the DIV90 PV-associated program detectable?** Yes, a partial transcriptional program is detectable: all four components are present in some cells (per-gene table below). **{d.fraction_div30:.1%}** of DIV30 cells meet the fixed DIV90 median PV score, and **{(cells.pv_n_detected>=2).mean():.1%}** detect at least two of its four genes. Population-relative rank alone cannot establish a program absent in absolute expression.',
      f'2. **Do methods support it?** Pairwise Spearman correlations range **{agreement.spearman.min():.2f} to {agreement.spearman.max():.2f}**. **{summary["multi_method_supported_n"]:,} ({summary["multi_method_supported_fraction"]:.1%})** cells are in the top20% under at least three methods. Transfer and rank share module genes; neighbor and sparse methods share the frozen reference labels and independent nonsurface features, so four methods are not four independent biological validations. The consensus is descriptive and disagreement remains explicit.',
      f'3. **Which developmental state carries it?** The top20 candidates are **{er.postmitotic_fraction:.1%} operational postmitotic**, **{er.cycling_fraction:.1%} cycling**, and **{er.progenitor_fraction:.1%} progenitor**. Consensus correlation with maturation is **rho={rho:.2f}**. Candidate resemblance concentrates after the progenitor-dominant context and is present among postmitotic ERBB4-RNA-negative profiles. These are operational, overlapping contextual masks and cross-sectional relationships, not demonstrated temporal compartments.',
      f'4. **Does it appear before ERBB4?** **{er.lost_fraction:.1%}** of the top20% consensus candidates have no detected ERBB4 RNA. This supports cross-sectional coexistence of a candidate signature with ERBB4 RNA non-detection; it does not prove chronological onset before ERBB4 or absent ERBB4 protein.',
      f'5. **Would ERBB4 lose candidates?** An ERBB4 RNA>0 entry rule loses **{int(er.lost_n):,}/{int(er.target_n):,} ({er.lost_fraction:.1%})** top20 candidates; LHX6 RNA>0 loses **{int(lh.lost_n):,} ({lh.lost_fraction:.1%})**. Do not impose the DIV90 entry gate at DIV30 based only on late-state biology.',
      '6. **Does FGFR2/PTPRS contain useful information?** The frozen full-precision DIV90 removal rule is evaluated only as a diagnostic benchmark in the surface tables, alongside gene-level effects, overlap, and confound adjustment. Same-direction association is not replication.',
      '7. **Is there a stronger DIV30 surface gate?** The independent surface search evaluates actual target composition under distinct target definitions and recovery floors. Its selected candidates, absolute improvements, conditional third-marker gains, and held-sample performance appear in the surface report below.',
      '8. **Are best markers specific or maturation/progenitor markers?** Surface regressions condition on progenitor, neuronal maturation, MGE identity, and depth; neighbor-score regression is a separate check. Sample-adjusted effects and held-sample metrics qualify classifications. Residual association is a candidate subtype signal, not causal specificity.',
      '9. **Strongest achievable enrichment?** The gate table records the best observed purity and fold enrichment at each recovery floor (30/50/70/80/90%), with starting prevalence, retained yield, target recovery, and cell counts. This is a bounded in-sample discovery ceiling; cross-sample validation and restricted triplet search prevent claiming a global or experimentally achievable optimum.',
      '10. **Is DIV30 worth sorting?** There is no validated prospective sort here. An exploratory sort is worth considering only if surface enrichment survives sample/line differences and maturation adjustment and a matching surface-protein phenotype is confirmed. A high recovery of a rank-defined state alone is insufficient.',
      '11. **Which gate should be tested?** Consult the candidate rule and matched-recovery alternatives in the surface report. These RNA thresholds cannot be translated into fluorescence thresholds; experimental endpoint evaluation must distinguish generic neuronal enrichment from later PV protein enrichment.',
      '12. **If no practical sort is supported, what is learned?** The ERBB4/LHX6 losses, contextual curves and score sensitivity define where a candidate transcriptional resemblance is present or uncertain. They can guide a later sampling/sorting window or culture perturbation hypothesis. No culture manipulation effect is established by this analysis.',
      '## Score definition, absolute support and sensitivity',
      'The transferred PV module remains **MEF2C, KCNC1, KCNC2, TAC1**. SST remains **SST, SATB1, NR2F2, CDK14, CACNG3**. Equal-weight bounded gene z scores use the unchanged DIV90 means and population standard deviations, not DIV30 rescaling. MAF/MAFB are absent from the primary module and tested only in a declared broad sensitivity.',
      md_table(detection),
      md_table(absolute),
      md_table(sensitivity),
      '## End-to-end score sensitivity',
      ('All four methods were recalculated after each gene drop: transferred scores, frozen-neighbor score labels, newly fitted sparse-regression labels, and within-cell gene-rank components. Primary mapping neighbors stay fixed because their features exclude every tested gene. The original analysis remains unchanged.\n\n'+md_table(fullsens.drop(columns=['interpretation'],errors='ignore'))) if len(fullsens) else 'End-to-end sensitivity not available.',
      ('Removing MEF2C leaves broad ranks stable but changes membership; removing both KCNC genes is more consequential. No single MEF2C dependence explains the entire consensus, but sensitivity of top candidate identity prevents treating these cells as committed future PV cells. The broad MAF/MAFB result is explicitly a sensitivity module, not a redefined primary target.'),
      '## Method agreement and developmental contexts',
      md_table(agreement),md_table(context),
      '## Available samples and metadata aliases',
      'The cached DIV30 H5AD contains six core samples, with preserved source identifiers `9583-MW-1` through `9583-MW-6`. The repository metadata uses `9853-MW-N`; the existing alias handling in `scripts/29_div30_jia_progenitor_metadata_umap_grid.R` establishes this join. `canonical_metadata_sample` records the alias without changing the source sample or any barcode. The six samples form H9, 79B, and 2E pairs under the user-provided CV/MW assignments. Metadata samples 7–9 (79B_old) are absent from this cached H5AD, so no result is claimed for them. All six present samples have complete line/condition labels; the presence audit is `tables/sample_metadata_and_presence.tsv`.',
      f'**{summary["mapping_out_of_reference_fraction"]:.1%}** of DIV30 cells exceed the 95th percentile of mean 30-neighbor distance from held-sample DIV90 mapping. This measures reference mismatch, not identity correctness. The DIV90 reference contains only selected cortical interneurons and will assign even unrelated/progenitor DIV30 profiles to their closest available cortical profile.',
      'Sparse ElasticNet has fixed alpha 0.02 and l1 ratio 0.9. Leave-sample and leave-line validation reselect reference-variable nonsurface features and fit scaling/PCA/models only in each training fold. Reference prediction success cannot validate backward developmental extrapolation. Mapping/supervised genes exclude the entire CSPA catalogue, all curated gate genes, every score/context gene, and benchmark context genes.',
      '## Sorting strategy models',
      '**Model A (same gate both ages):** requires actual transfer enrichment; the frozen gate benchmark alone is insufficient. **Model B (early gate → later cortical entry → subtype gate):** a developmental hypothesis if an adjusted, cross-sample DIV30 candidate is supported; no lineage is shown. **Model C (wait until DIV90):** remains appropriate if DIV30 methods/markers lack robust subtype information. **Model D (culture manipulation):** not directly tested; sample/operator/glucose correlations cannot establish intervention effects.',
      '## Biological/statistical limits and reproduction',
      'The strongest supported endpoint remains later PV protein at the system level, supplied by the user. None of these cells is longitudinally linked to that endpoint. RNA zeros may reflect dropout; no protein absence follows. Context panels are literature-informed, analyst-defined indicators, not validated human developmental clocks. The maturation panel includes DCX (newborn neuronal differentiation), and the progression panel reflects inhibitory differentiation rather than uniquely cortical fate. Gene/module drop sensitivities can remove much of a short signature; the direct-only table is separate from the end-to-end sensitivity that recalculates all methods. Surface-gate discovery is exploratory and does not establish independent statistical significance; cells do not replace biological replicates. CV/MW includes operator/culture and reported glucose differences, so glucose causation cannot be isolated. The 79B_old samples are present in metadata but absent from the frozen cached expression input.',
      f'Reproduce with `/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/python {REPO}/python_notebooks/scripts/div30_antecedent_analysis.py --stage all`. Runtime caches are reused; choose a new output version to change definitions. Code, exact gene lists/scales, reference neighbors, model coefficients, targets, fitted gates and checks are under `provenance/`, `cache/`, and `tables/`. Four BLAS threads are used and expression is streamed from HDF5.',
      '## Sources',
    ]
    for i,p in enumerate(paragraphs):
        for number,answer in gateanswers.items():
            if p.startswith(f'{number}. **'):paragraphs[i]=answer
    summary['sample_alias_reconciliation']=json.loads((OUT/'provenance/sample_alias_reconciliation.json').read_text())
    dump(OUT/'summary.json',summary)
    paragraphs.extend(f'- [{s["url"]}]({s["url"]}): {s["use"]}' for s in SOURCES)
    for candidate in ['SURFACE_REPORT.md','surface_report.md']:
        if (OUT/candidate).exists():paragraphs+=['## Separate surface-discovery report',(OUT/candidate).read_text()]
    (OUT/'REPORT.md').write_text('\n\n'.join(paragraphs)+'\n')
    # A tiny self-contained index leaves all results reviewable without an app.
    (OUT/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>DIV30 antecedent analysis</title><style>body{font:18px system-ui;max-width:1100px;margin:3em auto;color:#173c50}img{width:100%}a{color:#246a8b}</style><h1>DIV30 transcriptional antecedents</h1><p>Separate from DIV90. Developmental resemblance, not lineage tracing.</p><p><a href="REPORT.md">Complete report and twelve answers</a> · <a href="summary.json">Machine-readable summary</a> · <a href="tables/entry_marker_loss.tsv">ERBB4/LHX6 entry losses</a></p>'+''.join(f'<figure><img src="figures/png/{p.name}"><figcaption>{p.stem.replace("_"," ")}</figcaption></figure>' for p in sorted((OUT/'figures/png').glob('*.png'))))


def verify(cells):
    original=json.loads((OUT/'provenance/source_before.json').read_text());checks={}
    for path,v in original.items():
        p=Path(path);ok=p.stat().st_size==v['size'] and p.stat().st_mtime_ns==v['mtime_ns']
        if 'sha256' in v:ok=ok and sha(p)==v['sha256']
        checks[path]=ok
    assert all(checks.values())
    cfg=json.loads((OUT/'provenance/scoring_config.json').read_text())
    assert not set(cfg['mapping_features']) & (set(cfg['surface_catalogue_genes'])|set(EXCLUDED))
    assert cells.cell_id.is_unique and len(cells)==90631
    assert np.isfinite(cells[METHODS+['early_pv_score']].to_numpy()).all()
    assert cells['sample'].nunique()==6
    assert cells[['biological_label','cell_line','condition']].notna().all().all()
    assert set(cells.cell_line)=={'H9','79B','2E'}
    summary={'all_original_inputs_preserved':all(checks.values()),'inputs':checks,'unique_cell_ids':True,'n_cells':len(cells),'no_scoring_surface_mapping_feature_overlap':True,'finite_scores':True,'no_new_umap_or_clusters':True,'runtime_under_turbo':str(OUT)}
    dump(OUT/'provenance/validation.json',summary)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--stage',choices=['prepare','score','metadata','sensitivity','surface','figures','report','verify','all'],default='all');args=ap.parse_args()
    init()
    with threadpool_limits(limits=4):
        if args.stage in ['prepare','score','all']:
            prepared=prepare()
            if args.stage=='prepare':return
            cells,sst_stable=score_map(*prepared)
            development(cells)
        else:
            cells=pd.read_csv(OUT/'cells.tsv.gz',sep='\t');sst_stable=json.loads((OUT/'provenance/sst_support.json').read_text())['stable_for_low_sst_target']
        if args.stage=='metadata':
            cells=attach_metadata(cells)
            pending=OUT/'.cells.metadata.pending.tsv.gz';cells.to_csv(pending,sep='\t',index=False);pending.replace(OUT/'cells.tsv.gz')
            development(cells)
        if args.stage in ['sensitivity','all']:sensitivity_end_to_end(cells)
        if args.stage in ['surface','all']:
            from div30_antecedent_surface import run_surface
            cfg=json.loads((OUT/'provenance/scoring_config.json').read_text());expr=np.load(OUT/'cache/expression.npy',mmap_mode='r');genes=json.loads((OUT/'cache/genes.json').read_text())
            with threadpool_limits(limits=2):run_surface(OUT,cells,expr,genes,EXCLUDED+cfg['mapping_features'],sst_stable=sst_stable)
        if args.stage in ['sensitivity','surface','report','all']:sensitivity_fixed_gates(cells)
        if args.stage in ['figures','all']:figures(cells)
        if args.stage in ['report','all']:report(cells)
        if args.stage in ['verify','all']:verify(cells)
    print(f'Completed DIV30 {args.stage}: {OUT}',flush=True)


if __name__=='__main__':main()
