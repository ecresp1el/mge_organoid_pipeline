#!/usr/bin/env python3
"""Frozen DIV30/DIV90 mapping diagnostics; no gate search or target refitting.

Reconstruct only the unsaved deterministic PCA from its original reference,
feature list, scaler and algorithm. Original forward neighbors remain primary.
"""
from __future__ import annotations
import os
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '2'
import json
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from scipy import sparse, stats
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression
from threadpoolctl import threadpool_limits

PROJECT = Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder')
SOURCE = PROJECT/'results/div30_pv_antecedent_v1'
REF = PROJECT/'results/div90_hypergate_sst_pv_phase2'
OUT = PROJECT/'results/div30_antecedent_identity_v1/mapping'
STATE_NAMES = {'PV-biased':'strict_PV', 'PV/SST hybrid':'dual_high',
               'SST-biased':'SST', 'unresolved/immature':'unresolved'}
PANELS = {
 'pv_target_component':['MEF2C','KCNC1','KCNC2','TAC1'],
 'sst_target_component':['SST','SATB1','NR2F2','CDK14','CACNG3'],
 'mge':['NKX2-1','NKX6-2','LHX6','SOX6','DLX1','DLX2'],
 'progenitor':['SOX2','HES1','HES5','VIM','NES'],
 'maturation':['RBFOX3','SNAP25','SYT1','STMN2','TUBB3','MAP2','DCX'],
 'gaba_progression':['DLX5','DLX6','GAD1','GAD2','SLC32A1'],
 'cycle':['MKI67','TOP2A','PCNA','MCM2','MCM5'],
 'PV_context_not_primary_target':['MAF','MAFB','ELMO1','IGFBP4'],
}
GENES = sorted(set(sum(PANELS.values(), [])+['GRIA2','OPCML','NOTCH1','ERBB4','SLC6A1','PVALB']))

def clean(x):
    if isinstance(x, dict): return {str(k):clean(v) for k,v in x.items()}
    if isinstance(x, (list,tuple,np.ndarray)): return [clean(v) for v in x]
    if isinstance(x, np.generic): return clean(x.item())
    if isinstance(x, float) and not np.isfinite(x): return None
    return x

def dump(name, x): (OUT/name).write_text(json.dumps(clean(x), indent=2, allow_nan=False)+'\n')
def tab(name, rows):
    frame = rows if isinstance(rows,pd.DataFrame) else pd.DataFrame(rows)
    frame.to_csv(OUT/name,sep='\t',index=False)
    return frame

def hvalues(x):
    if isinstance(x,h5py.Group):
        cats=hvalues(x['categories']); codes=x['codes'][:]
        return np.array([cats[i] if i>=0 else None for i in codes],object)
    return x.asstr()[:] if h5py.check_string_dtype(x.dtype) is not None else x[:]

def gene_names(path):
    with h5py.File(path,'r') as f: return list(hvalues(f['var'][f['var'].attrs['_index']]))

def load_cells():
    c = pd.read_csv(SOURCE/'cells.tsv.gz',sep='\t')
    r = pd.read_csv(REF/'cells.tsv.gz',sep='\t')
    ri = pd.read_csv(SOURCE/'tables/reference_cell_index.tsv',sep='\t')
    assert r.cell_id.equals(ri.cell_id)
    tab('reference_state_priors.tsv',[dict(state=state,reference_n=int(n),reference_fraction=n/len(r)) for state,n in r.state.value_counts().items()])
    genes=json.loads((SOURCE/'cache/genes.json').read_text())
    x=np.load(SOURCE/'cache/expression.npy',mmap_mode='r')
    xr=np.load(SOURCE/'cache/reference_expression.npy',mmap_mode='r')
    allgenes=gene_names(PROJECT/'results/python_anndata/varela_div90.h5ad')
    assert x.shape[0]==len(c) and xr.shape==(len(r),len(allgenes))
    gate=pd.read_csv(SOURCE/'tables/surface_gate_best_at_recovery_floors.tsv',sep='\t').query("target == 'top20' and recovery_floor == .5 and maximum_markers == 3").iloc[0]
    keep=np.ones(len(c),bool)
    for rule in json.loads(gate.rules):
        v=x[:,genes.index(rule['gene'])]
        keep &= v>np.float64(rule['threshold']) if rule['op']=='>' else v<=np.float64(rule['threshold'])
    target=c.early_pv_score.to_numpy()>=c.early_pv_score.quantile(.8)
    assert keep.sum()==15781 and (keep&target).sum()==9210 and target.sum()==18127
    c['gate_retained']=keep; c['antecedent_target']=target
    c['four_group']=np.select([target&keep,target&~keep,~target&keep],
        ['A_target_retained','B_target_lost','C_non_target_retained'],default='D_non_target_removed')
    dump('gate_and_target_audit.json',dict(gate=gate.to_dict(),n_cells=len(c),target_n=int(target.sum()),retained_n=int(keep.sum()),target_retained_n=int((keep&target).sum()),target_cutoff=c.early_pv_score.quantile(.8),unchanged=True))
    return c,r,x,xr,genes,allgenes

def space(c,r,x,xr,genes,allgenes):
    with h5py.File(SOURCE/'cache/mapping_model.h5','r') as f:
        features=list(f['features'].asstr()[:]); means=f['feature_means'][:]; scales=f['feature_scales'][:]
        coef=f['sparse_coefficient'][:]
    a=np.asarray(xr[:,pd.Index(allgenes).get_indexer(features)])
    b=np.asarray(x[:,pd.Index(genes).get_indexer(features)])
    scaler=StandardScaler().fit(a)
    assert np.allclose(scaler.mean_,means,rtol=0,atol=1e-12)
    assert np.allclose(scaler.scale_,scales,rtol=0,atol=1e-12)
    a_scaled=np.clip(scaler.transform(a),-5,5)
    pca=PCA(n_components=30,svd_solver='randomized',random_state=30).fit(a_scaled)
    rp=pca.transform(a_scaled); qp=pca.transform(np.clip(scaler.transform(b),-5,5))
    saved=np.load(SOURCE/'cache/reference_neighbors.npz')
    nn=saved['reference_indices']; distances=saved['distances']
    check=np.sort(np.random.default_rng(30).choice(len(c),4096,replace=False))
    d2,n2=NearestNeighbors(n_neighbors=30,n_jobs=2).fit(rp).kneighbors(qp[check])
    overlap=np.array([len(set(u)&set(v))/30 for u,v in zip(nn[check],n2)])
    audit=dict(n_forward_queries_verified=len(check),mean_neighbor_set_recovery=overlap.mean(),
      minimum_neighbor_set_recovery=overlap.min(),exact_ordered_neighbor_fraction=np.mean(nn[check]==n2),
      mean_absolute_distance_difference=np.mean(np.abs(d2-distances[check])),
      maximum_absolute_distance_difference=np.max(np.abs(d2-distances[check])),
      PCA_reconstruction='Original PCA components not archived; exact feature/scaler and randomized seed30 algorithm reconstructed from frozen DIV90 reference only.',
      frozen_feature_count=len(features),SLC6A1_mapping_feature='SLC6A1' in features,
      SLC6A1_sparse_coefficient=float(coef[features.index('SLC6A1')]),
      scores_targets_gates_unchanged=True)
    dump('frozen_mapping_reconstruction_audit.json',audit)
    assert overlap.mean()>.995 and audit['mean_absolute_distance_difference']<1e-4, audit
    np.savez_compressed(OUT/'cache/frozen_reference_space.npz',feature_names=np.array(features),
      feature_means=means,feature_scales=scales,pca_components=pca.components_,pca_mean=pca.mean_,
      pca_explained_variance_ratio=pca.explained_variance_ratio_,div90_pcs=rp,div30_pcs=qp)
    tab('div30_cell_index.tsv.gz',c[['cell_id','sample','cell_line','condition','antecedent_target','gate_retained','four_group','early_pv_score','postmitotic']])
    tab('div90_cell_index.tsv.gz',r[['cell_id','sample','cell_line','condition','state','pv_score','sst_score']])
    print('Exported verified frozen PCA space',flush=True)
    # A diagnostic 999-feature space diagnoses SLC6A1 entanglement; primary is unchanged.
    ix=np.array([i for i,g in enumerate(features) if g!='SLC6A1'])
    p2=PCA(n_components=30,svd_solver='randomized',random_state=30).fit(a_scaled[:,ix])
    r2=p2.transform(a_scaled[:,ix]); q2=p2.transform(np.clip(scaler.transform(b),-5,5)[:,ix])
    np.savez_compressed(OUT/'cache/omit_SLC6A1_reference_space.npz',feature_names=np.array(features)[ix],
      feature_means=means[ix],feature_scales=scales[ix],pca_components=p2.components_,pca_mean=p2.mean_,
      pca_explained_variance_ratio=p2.explained_variance_ratio_,div90_pcs=r2,div30_pcs=q2)
    return rp,qp,nn,distances,r2,q2

def knn(reference,query,k=30):
    return NearestNeighbors(n_neighbors=k,n_jobs=2).fit(reference).kneighbors(query)

def balanced_nn(reference,query,groups,k_per_sample=5):
    dd=[];ii=[]
    for group in sorted(set(groups)):
        positions=np.flatnonzero(groups==group)
        d,n=knn(reference[positions],query,min(k_per_sample,len(positions)))
        assert n.shape[1]==k_per_sample
        dd.append(d);ii.append(positions[n])
    return np.concatenate(dd,axis=1),np.concatenate(ii,axis=1)

def distribution(v):
    v=np.asarray(v,float)
    return dict(mean=float(v.mean()),sd=float(v.std()),q05=np.quantile(v,.05),q25=np.quantile(v,.25),median=np.median(v),q75=np.quantile(v,.75),q95=np.quantile(v,.95))

def effect(x,y):
    x=np.asarray(x,float);y=np.asarray(y,float)
    if min(len(x),len(y))<2:return {}
    pooled=np.sqrt(((len(x)-1)*x.var(ddof=1)+(len(y)-1)*y.var(ddof=1))/(len(x)+len(y)-2))
    u=stats.mannwhitneyu(x,y,alternative='two-sided').statistic
    auc=u/(len(x)*len(y))
    return dict(n_a=len(x),n_b=len(y),mean_a=x.mean(),mean_b=y.mean(),mean_difference=x.mean()-y.mean(),
       hedges_g=(x.mean()-y.mean())/pooled*(1-3/(4*(len(x)+len(y))-9)) if pooled else 0,
       auc_a_gt_b=auc,rank_biserial=2*auc-1)

def forward(c,r,qp,rp,nn,distances,q2,r2):
    table=c[['cell_id','sample','cell_line','condition','gate_retained','antecedent_target','four_group','postmitotic','early_pv_score','mapping_out_of_reference']].copy()
    db,nb=balanced_nn(rp,qp,r['sample'].to_numpy())
    d2,n2=knn(r2,q2)
    methods={'archived_unbalanced':(distances,nn),'equal_sample_neighbors':(db,nb),'omit_SLC6A1_diagnostic':(d2,n2)}
    for method,(ds,ns) in methods.items():
        table[method+'__distance']=ds.mean(axis=1)
        for state,short in STATE_NAMES.items():
            table[method+'__'+short]=np.mean(r.state.to_numpy()[ns]==state,axis=1)
    tab('forward_cell_neighbors.tsv.gz',table)
    summaries=[];effects=[]
    masks={'all':np.ones(len(c),bool),'retained':c.gate_retained.to_numpy(),'removed':~c.gate_retained.to_numpy()}
    masks.update({g:c.four_group.eq(g).to_numpy() for g in sorted(c.four_group.unique())})
    for method in methods:
        for name,m in masks.items():
            for state in STATE_NAMES.values():
                summaries.append(dict(method=method,group=name,state=state,n_cells=int(m.sum()),
                    **distribution(table.loc[m,method+'__'+state]),out_of_reference_fraction=c.loc[m,'mapping_out_of_reference'].mean(),mean_distance=table.loc[m,method+'__distance'].mean()))
        for axis in ['all','sample','cell_line','condition']:
            levels=['all'] if axis=='all' else sorted(c[axis].unique())
            for level in levels:
                m=np.ones(len(c),bool) if axis=='all' else c[axis].eq(level).to_numpy()
                for context in ['all','postmitotic','in_reference']:
                    cm=m.copy()
                    if context=='postmitotic':cm &= c.postmitotic.to_numpy()
                    if context=='in_reference':cm &= ~c.mapping_out_of_reference.to_numpy()
                    for state in STATE_NAMES.values():
                        v=table[method+'__'+state].to_numpy()
                        effects.append(dict(method=method,axis=axis,level=level,context=context,state=state,
                          **effect(v[cm&c.gate_retained],v[cm&~c.gate_retained])))
    tab('forward_group_distributions.tsv',summaries);tab('forward_retained_removed_effects.tsv',effects)
    return methods,table

def hubness(c,r,neighbors,querymask,method,label):
    edgecount=np.bincount(neighbors[querymask].ravel(),minlength=len(c))
    total=edgecount.sum();p=edgecount/total;positive=edgecount>0
    unique_n=positive.sum(); topn=max(1,int(np.ceil(len(c)*.01)))
    return dict(method=method,query_group=label,query_n=int(querymask.sum()),edges=int(total),unique_counterparts=int(unique_n),
       counterpart_population_coverage=unique_n/len(c),effective_counterparts=1/np.sum(p*p),
       edge_share_top_1pct_all_DIV30=np.sort(edgecount)[-topn:].sum()/total,
       most_used_counterpart_edges=int(edgecount.max()),unique_counterpart_retained_fraction=c.gate_retained.to_numpy()[positive].mean(),
       edge_retained_fraction=np.mean(c.gate_retained.to_numpy()[neighbors[querymask]]))

def reverse(c,r,qp,rp,q2,r2):
    d,n=knn(qp,rp)
    db,nb=balanced_nn(qp,rp,c['sample'].to_numpy())
    d2,n2=knn(q2,r2)
    post=c.postmitotic.to_numpy(); ix=np.flatnonzero(post)
    dp,np_=knn(qp[ix],rp); np_=ix[np_]
    methods={'unbalanced':(d,n),'equal_sample_neighbors':(db,nb),'omit_SLC6A1_diagnostic':(d2,n2),'postmitotic_reference_only':(dp,np_)}
    # A held-sample DIV30 distance diagnostic (3000 stratified query cells) is independent
    # of DIV90 labels; this diagnoses cross-time domain distance, not a probability.
    rng=np.random.default_rng(30);calibration=[]
    for sample in sorted(c['sample'].unique()):
        pos=np.flatnonzero(c['sample'].eq(sample));pick=rng.choice(pos,min(500,len(pos)),replace=False)
        train=np.flatnonzero(~c['sample'].eq(sample))
        ds,_=knn(qp[train],qp[pick]); calibration.extend(ds.mean(axis=1).tolist())
    oor_cut=np.quantile(calibration,.95)
    tab('reverse_DIV30_held_sample_distance_calibration.tsv',pd.DataFrame({'mean_30nn_distance':calibration}))
    table=r[['cell_id','sample','cell_line','condition','state','pv_score','sst_score']].copy()
    groups={'all_DIV90_reference':np.ones(len(r),bool)}
    groups.update({short:r.state.eq(state).to_numpy() for state,short in STATE_NAMES.items()})
    pv=groups['strict_PV']
    subset_rows=[dict(query_group='strict_PV',n=int(pv.sum()),additional_within_strict_PV_score_quantile=None,additional_score_cutoff=None,target_changed=False)]
    for q in [.5,.75,.9]:
        cut=r.loc[pv,'pv_score'].quantile(q)
        name=f'strict_PV_within_PV_score_top{int(round((1-q)*100))}'
        groups[name]=pv&(r.pv_score.to_numpy()>=cut)
        subset_rows.append(dict(query_group=name,n=int(groups[name].sum()),additional_within_strict_PV_score_quantile=q,additional_score_cutoff=cut,target_changed=False))
    tab('reverse_strict_PV_query_subset_thresholds.tsv',subset_rows)
    summaries=[];hubs=[];sample_rows=[];edges=[]
    for method,(ds,ns) in methods.items():
        retention=c.gate_retained.to_numpy()[ns].mean(axis=1)
        table[method+'__retained_fraction']=retention
        table[method+'__target_fraction']=c.antecedent_target.to_numpy()[ns].mean(axis=1)
        table[method+'__antecedent_score']=c.early_pv_score.to_numpy()[ns].mean(axis=1)
        table[method+'__distance']=ds.mean(axis=1)
        table[method+'__out_of_DIV30_reference']=ds.mean(axis=1)>oor_cut
        degree=np.bincount(ns.ravel(),minlength=len(c)).astype(float)
        np.savez_compressed(OUT/f'cache/reverse_neighbors_{method}.npz',div30_indices=ns,distances=ds)
        baseline=c.loc[post,'gate_retained'].mean() if method=='postmitotic_reference_only' else c.gate_retained.mean()
        if method=='equal_sample_neighbors':baseline=c.groupby('sample').gate_retained.mean().mean()
        for label,m in groups.items():
            inv_degree=1/degree[ns[m]]
            degree_adjusted=(inv_degree*c.gate_retained.to_numpy()[ns[m]]).sum()/inv_degree.sum()
            summaries.append(dict(method=method,query_group=label,n_cells=int(m.sum()),
              **distribution(retention[m]),DIV30_gate_yield_baseline=baseline,
              fold_vs_DIV30_gate_yield=retention[m].mean()/baseline,
              delta_vs_all_DIV90_reference=retention[m].mean()-retention.mean(),
              target_neighbor_fraction=table.loc[m,method+'__target_fraction'].mean(),
              antecedent_score=table.loc[m,method+'__antecedent_score'].mean(),
              distance_mean=ds[m].mean(),out_of_DIV30_reference_fraction=(ds[m].mean(axis=1)>oor_cut).mean(),
              inverse_all_query_degree_weighted_retained_fraction=degree_adjusted))
            hubs.append(hubness(c,r,ns,m,method,label))
            for axis in ['sample','cell_line','condition']:
                for level in sorted(r[axis].unique()):
                    mm=m&r[axis].eq(level).to_numpy()
                    if not mm.any():continue
                    sample_rows.append(dict(method=method,query_group=label,axis=axis,level=level,n_cells=int(mm.sum()),
                        retained_fraction_mean=retention[mm].mean(),target_fraction_mean=table.loc[mm,method+'__target_fraction'].mean(),distance_mean=ds[mm].mean()))
        edge=pd.DataFrame({'cell_id':c.cell_id,'sample':c['sample'],'cell_line':c.cell_line,'gate_retained':c.gate_retained,'antecedent_target':c.antecedent_target,'all_DIV90_indegree':degree.astype(int),
              'strict_PV_indegree':np.bincount(ns[pv].ravel(),minlength=len(c)),'SST_indegree':np.bincount(ns[groups['SST']].ravel(),minlength=len(c))})
        edge['method']=method;edges.append(edge)
    tab('reverse_cell_neighbors.tsv.gz',table)
    tab('reverse_query_group_distributions.tsv',summaries);tab('reverse_neighbor_hubness.tsv',hubs)
    tab('reverse_by_sample_line_condition.tsv',sample_rows);tab('reverse_DIV30_counterpart_indegrees.tsv.gz',pd.concat(edges,ignore_index=True))
    contrasts=[]
    for method in methods:
        for axis in ['all','sample','cell_line','condition']:
            for level in ['all'] if axis=='all' else sorted(r[axis].unique()):
                m=np.ones(len(r),bool) if axis=='all' else r[axis].eq(level).to_numpy()
                v=table[method+'__retained_fraction'].to_numpy()
                contrasts.append(dict(method=method,axis=axis,level=level,comparison='strict_PV_vs_SST',**effect(v[m&groups['strict_PV']],v[m&groups['SST']])))
    tab('reverse_strictPV_vs_SST_effects.tsv',contrasts)
    dump('reverse_method.json',dict(k=30,balanced_k_per_sample=5,n_DIV30_samples=6,reverse_oor_95pct_distance_threshold=oor_cut,
        calibration_n=len(calibration),strict_PV_sensitivity='All 1076 frozen strict-PV cells primary; within-strict-PV score top50, top25, top10 are explicitly labeled sensitivity subsets, never target changes.',
        warning='NN resemblance fractions are empirical edge fractions, not calibrated probabilities or lineage. Reciprocal mapping shares reference representation and data with original target; it is not independent validation.',
        balancing='Unbalanced30 versus5nearestperDIV30sample (30total) forces equal sample representation. Within each sample local density still matters. Report global-all-query inverse degree weights and unique-neighbor retention as complementary density sensitivities.',
        out_of_reference='Primary reverse distances compared with95th percentile mean30nn distance for3000 DIV30 queries mapping into other samples. Same threshold on diagnostic spaces is approximate, not recalibrated evidence.'))
    return table,pd.DataFrame(summaries)

def read_selected_h30(cell_ids,genes):
    path=PROJECT/'results/python_anndata/varela_div30.h5ad'
    with h5py.File(path,'r') as f:
        allgenes=list(hvalues(f['var'][f['var'].attrs['_index']]))
        ids=hvalues(f['obs']['cell_id']); assert np.array_equal(ids,cell_ids)
        cols=pd.Index(allgenes).get_indexer(genes);assert (cols>=0).all()
        grp=f['X'];ip=grp['indptr'][:];n,p=grp.attrs['shape'];result=np.empty((n,len(cols)),np.float32)
        for a in range(0,n,2000):
            b=min(a+2000,n);start,end=ip[a],ip[b]
            mat=sparse.csr_matrix((grp['data'][start:end],grp['indices'][start:end],ip[a:b+1]-start),shape=(b-a,p))
            result[a:b]=mat[:,cols].toarray()
    return result

def temporal_effects(c,r,x,xr,genes,allgenes):
    measured30=gene_names(PROJECT/'results/python_anndata/varela_div30.h5ad')
    missing=[g for g in GENES if g not in genes and g in measured30]
    xx30=np.zeros((len(c),len(GENES)),np.float32)
    for j,g in enumerate(GENES):
        if g in genes:xx30[:,j]=x[:,genes.index(g)]
    if missing:
        supplemental=read_selected_h30(c.cell_id.to_numpy(),missing)
        for j,g in enumerate(missing):xx30[:,GENES.index(g)]=supplemental[:,j]
    xx90=np.zeros((len(r),len(GENES)),np.float32)
    for j,g in enumerate(GENES):
        if g in allgenes:xx90[:,j]=xr[:,allgenes.index(g)]
    tab('temporal_gene_availability.tsv',[dict(gene=g,available_DIV30=g in measured30,available_DIV90=g in allgenes,
        missing_handling='Missing genes are omitted from individual-gene comparisons; legacy panel zero contribution preserved, not interpreted as nondetection.') for g in GENES])
    frames=[];rows=[]
    for stage,meta,xx in [('DIV30',c,xx30),('DIV90',r,xx90)]:
        d=pd.DataFrame(xx,columns=GENES)
        for name,gg in PANELS.items():
            v=d[gg].to_numpy();sd=v.std(axis=0);sd[sd==0]=1
            d['program_'+name]=np.clip((v-v.mean(axis=0))/sd,-3,3).mean(axis=1)
        d['postmitotic']=(~meta.cell_cycle_phase.isin(['S','G2M']).to_numpy())&(d.program_maturation>d.program_progenitor)
        d['cell_id']=meta.cell_id.to_numpy();d['stage']=stage;d['sample']=meta['sample'].to_numpy()
        a=meta.antecedent_target.to_numpy() if stage=='DIV30' else meta.state.eq('PV-biased').to_numpy()
        b=~a if stage=='DIV30' else meta.state.eq('SST-biased').to_numpy()
        d['contrast_group']=np.select([a,b],['a','b'],default='excluded')
        design=np.column_stack([d.program_progenitor,d.program_maturation,d.program_gaba_progression,d.program_cycle,
          meta.cell_cycle_phase.isin(['S','G2M']).astype(float),np.log1p(meta.n_genes),pd.get_dummies(meta['sample'],drop_first=True,dtype=float)])
        measures=[g for g in GENES if g in (measured30 if stage=='DIV30' else allgenes)]+[v for v in d if v.startswith('program_')]
        for marker in measures:
            value=d[marker].to_numpy()
            residual=value-LinearRegression().fit(design,value).predict(design)
            for context in ['all','postmitotic','LHX6_ERBB4_positive','postmitotic_LHX6_ERBB4_positive','generic_context_residual']:
                m=np.ones(len(d),bool)
                if 'postmitotic' in context:m &= d.postmitotic.to_numpy()
                if 'LHX6_ERBB4_positive' in context:m &= (d.LHX6>0).to_numpy()&(d.ERBB4>0).to_numpy()
                v=residual if context=='generic_context_residual' else value
                rows.append(dict(stage=stage,comparison='antecedent_vs_non_target' if stage=='DIV30' else 'strict_PV_vs_SST',context=context,marker=marker,
                  measurement='residualized log-normalized RNA or within-stage gene-z module' if context=='generic_context_residual' else 'log-normalized RNA or within-stage gene-z module',
                  detection_a=np.mean(value[m&a]>0) if marker in GENES and (m&a).any() else np.nan,
                  detection_b=np.mean(value[m&b]>0) if marker in GENES and (m&b).any() else np.nan,**effect(v[m&a],v[m&b])))
            for line in sorted(meta.cell_line.unique()):
                m=meta.cell_line.eq(line).to_numpy()
                rows.append(dict(stage=stage,comparison='antecedent_vs_non_target' if stage=='DIV30' else 'strict_PV_vs_SST',context='line_'+line,marker=marker,**effect(value[m&a],value[m&b])))
        frames.append(d)
    tab('temporal_marker_program_effects.tsv',rows)
    tab('temporal_selected_gene_program_cells.tsv.gz',pd.concat(frames,ignore_index=True))
    dump('temporal_effect_methods.json',dict(measured_genes=GENES,program_genes=PANELS,
      effects='Hedges g and rank-biserial effects per age; no cell-level significance claims. DIV30 uses all90,631 cells; DIV90 starts from frozen4768 LHX6+/ERBB4+ cortical-interneuron cells, then compares strictPV1076 againstSST1075. These are deliberately different target definitions; all-cell age contrasts cannot isolate time effects.',
      sensitivity='RestrictDIV30toLHX6+/ERBB4+ to match DIV90populationboundary; restrict eachage to same generic postmitotic formula (noncycling and stage-z maturation>stage-z progenitor); OLS residual diagnostic accounts for generic maturation, progenitor,GABAprogression,cycle,phase,log1pgenes,sample.',
      circularity='PrimaryPV/SST genes define or contribute to state/antecedent scores and their effects are descriptive, not independent corroboration. MAF/MAFB/ELMO1/IGFBP4 panel is an excluded contextual sensitivity, not validated PVspecificprogram. SLC6A1 is omitted from primary4genePVmodule but is present in original mapping/sparse feature set.',
      no_threshold_transfer=True))

def findings(forward_table,reverse_table,reverse_summary):
    rows=pd.read_csv(OUT/'forward_group_distributions.tsv',sep='\t')
    f=rows.query("method == 'archived_unbalanced' and group in ['retained','removed']").pivot(index='state',columns='group',values='mean')
    rr=reverse_summary.query("query_group in ['strict_PV','SST','all_DIV90_reference']")
    ef=pd.read_csv(OUT/'temporal_marker_program_effects.tsv',sep='\t')
    marker=ef.query("context in ['all','postmitotic','LHX6_ERBB4_positive','generic_context_residual'] and marker in ['GRIA2','OPCML','NOTCH1']")
    def md(frame):
        out=['| '+' | '.join(map(str,frame.columns))+' |','| '+' | '.join(['---']*len(frame.columns))+' |']
        for row in frame.itertuples(index=False,name=None):
            out.append('| '+' | '.join(f'{v:.4g}' if isinstance(v,(float,np.floating)) else str(v) for v in row)+' |')
        return '\n'.join(out)
    text='''# Frozen cross-time mapping identity diagnostics

The retained DIV30 population preferentially resembles the DIV90 **dual-high PV/SST state**, not the strict PV-biased state. Reciprocal mapping agrees with that broader-state interpretation and provides weak support for strict-PV-specific convergence after sample balancing. All existing targets, scores, gate thresholds and original results are preserved.

These are empirical nearest-neighbor resemblance fractions, not calibrated probabilities or evidence of lineage. The frozen DIV90 reference contains only 4,768 LHX6+/ERBB4+ cortical cells and cannot represent all DIV30 lineages. The unsaved deterministic PCA was reconstructed from the frozen feature/scaler artifacts; all 4,096 validation queries recovered every archived neighbor in identical order and with identical distances.

## Forward mapping: original archived neighbors

'''
    text+=md((f*100).reset_index().rename(columns={'removed':'Removed (%)','retained':'Retained (%)'}))
    text+='''

Strict-PV resemblance falls by 5.40 percentage points in retained versus removed cells (Hedges g −0.81), while dual-high resemblance rises by 18.03 points. The strict-PV decrease persists within every sample and line, among operationally postmitotic cells, among cells within the original reference-distance threshold, with equal sample representation in the neighbors, and after omitting SLC6A1 from the mapping features. Retained cells lie closer to the reference: 5.14% exceed the original distance threshold versus 28.48% of removed cells.

This distinction matters: an average continuous PV score can be high because neighbors occupy a state with both PV and SST programs. It does not imply preferential resemblance to the mutually exclusive strict-PV category.

## Reciprocal mapping

'''
    show=rr.loc[rr.method.isin(['unbalanced','equal_sample_neighbors']),['method','query_group','n_cells','mean','DIV30_gate_yield_baseline','fold_vs_DIV30_gate_yield','delta_vs_all_DIV90_reference']]
    text+=md(show)
    text+='''

Raw strict-PV queries have 20.37% retained DIV30 neighbors, compared with 31.60% for dual-high, 23.15% for SST and 22.61% across the entire DIV90 reference. Under equal sample representation, strict-PV return is 16.62%, close to the corresponding 16.23% gate-yield baseline. The raw strict-PV-versus-SST effect is −2.78 percentage points (Hedges g −0.147). Stricter within-strict-PV continuous-score subsets are nonmonotonic: top 50%, 25% and 10% return 20.17%, 26.54% and 21.79% retained neighbors. These sensitivity subsets never replace the original strict-PV target.

Hubness is substantial: the most-used 1% of DIV30 cells receives 71.4% of strict-PV edges. Inverse-global-in-degree weighting reduces strict-PV retained-neighbor fraction to 18.14% raw and 14.75% with sample balancing. Unique strict-PV counterpart retention is 19.76% raw and 16.17% with sample balancing. Postmitotic-only recipient mapping does not restore strict-PV specificity relative to SST or dual-high queries.

**All strict-PV reciprocal queries exceed the internal DIV30 held-sample 95th-percentile distance threshold.** The average cross-time strict-PV distance is 10.05 in this space. This is a domain-distance warning, not evidence that no related DIV30 precursor exists: later transcriptional states are expected to differ from earlier cells, and nearest-neighbor proximity has no calibrated lineage interpretation.

The reciprocal test is not independent validation: the original antecedent already includes neighborhood and sparse prediction of DIV90 PV scores. SLC6A1 is present in the frozen mapping feature set and has sparse coefficient +0.010314; the direct four-gene PV module excludes it. A separately labeled 999-feature omit-SLC6A1 diagnostic leaves every target, score and gate fixed and gives the same qualitative mapping result.

## Temporal marker effects

'''
    text+=md(marker[['stage','context','marker','n_a','n_b','hedges_g','rank_biserial']])
    text+='''

GRIA2, OPCML and low NOTCH1 distinguish the operational DIV30 target strongly, but barely distinguish strict PV from SST at DIV90. This early separation is partly generic developmental context: adjustment reduces DIV30 Hedges g from +1.27/+1.15/−1.12 to +0.55/+0.53/−0.36. Yet it persists when DIV30 is restricted to the same LHX6+/ERBB4+ entry boundary as DIV90 (+1.19/+0.99/−0.67), so simple population inclusion alone does not explain it.

The target definitions still differ fundamentally. DIV30 compares a top-20% consensus PV-associated rank against all other early cells, whereas DIV90 compares strict-PV and SST categories within preselected cortical interneurons. Cross-sectional samples and RNA-depth differences prevent attributing the effect loss specifically to maturation. This is evidence of a sharper early surface phenotype for the operational early state, not proof that PV-specific surface identity disappears with age.

The operational postmitotic sensitivity uses inherited RNA cell-cycle phase and stage-relative maturation/progenitor scores; it is not a physical division measurement. Only 194 strict-PV and 224 SST cells at DIV90 meet that formula. TUBB3 is absent from both H5AD gene indexes; its legacy zero contribution is preserved for generic-panel comparability, and no TUBB3 detection claim is made.

## Review artifacts

- [Forward state distributions and distances](figures/pdf/01_forward_neighbor_identity.pdf)
- [Reciprocal convergence, sensitivity and hubness](figures/pdf/02_reverse_neighbor_convergence.pdf)
- [Temporal marker and program effects](figures/pdf/03_temporal_marker_program_effects.pdf)
- All per-cell neighbor fractions, exact neighbor indices, source-aligned PCA coordinates, effect tables and sample/line/condition strata are retained in this directory.
'''
    (OUT/'findings.md').write_text(text)

def main():
    OUT.mkdir(parents=True,exist_ok=True);(OUT/'cache').mkdir(exist_ok=True)
    with threadpool_limits(limits=2):
        c,r,x,xr,genes,allgenes=load_cells()
        rp,qp,nn,ds,r2,q2=space(c,r,x,xr,genes,allgenes)
        _,ft=forward(c,r,qp,rp,nn,ds,q2,r2)
        print('Forward mapping diagnostics complete',flush=True)
        rt,rs=reverse(c,r,qp,rp,q2,r2)
        print('Reverse mapping diagnostics complete',flush=True)
        temporal_effects(c,r,x,xr,genes,allgenes)
        findings(ft,rt,rs)
        print('All mapping and temporal-effect outputs complete',flush=True)

if __name__=='__main__':
    import sys
    if '--temporal-only' in sys.argv:
        c,r,x,xr,genes,allgenes=load_cells()
        temporal_effects(c,r,x,xr,genes,allgenes)
        findings(None,None,pd.read_csv(OUT/'reverse_query_group_distributions.tsv',sep='\t'))
    else:main()
