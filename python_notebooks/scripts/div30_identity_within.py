#!/usr/bin/env python3
"""Describe the frozen DIV30 antecedent and fixed gate; never search or relabel.

All original artifacts are opened read-only. Supplementary programs are clearly
named and cannot change the frozen target, gate, or reference mapping.
"""
from __future__ import annotations
import os
import sys
for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_key] = '2'
import json
import itertools
import shutil
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from scipy import sparse, stats
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.preprocessing import StandardScaler, SplineTransformer
from sklearn.metrics import roc_auc_score, average_precision_score, log_loss, brier_score_loss, r2_score, mean_squared_error
from threadpoolctl import threadpool_limits

ROOT=Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results')
OLD=ROOT/'div30_pv_antecedent_v1'
OUT=ROOT/'div30_antecedent_identity_v1/within'
H30=ROOT/'python_anndata/varela_div30.h5ad'
MARKERS=['GRIA2','OPCML','NOTCH1']
REQUESTED=['NKX2-1','LHX6','SOX6','DLX1','DLX2','DLX5','DLX6','MAF','MAFB','MEF2C','ERBB4','KCNC1','KCNC2','SLC6A1','SST','SATB1','TAC1','PVALB']
SUPPLEMENTARY={
    'gaba':['GAD1','GAD2','SLC32A1','ABAT','GABRA1','GABRB2','GABRG2'],
    'stress':['FOS','JUN','JUNB','ATF3','HSPA1A','HSPA1B','DDIT3'],
    'glutamatergic_context':['SLC17A6','SLC17A7','TBR1','SATB2'],
    'glial_context':['GFAP','AQP4','S100B','OLIG1','OLIG2','PDGFRA','SOX10'],
    'cholinergic_context':['CHAT','SLC18A3','LHX8','ISL1'],
    'nonmodule_pv_context':['MAF','MAFB','ERBB4','LHX6','SOX6'],
}
PROGRAMS=['early_pv_score','direct_pv_score','direct_sst_score','mge_score','progression_score','progenitor_score','maturation_score','cycle_score','postmitotic','cycling','progenitor','n_genes','log_counts','percent_mt','neighbor_pv_score','neighbor_sst_score','mapping_distance','mapping_out_of_reference']+[f'supplementary_{x}_score' for x in SUPPLEMENTARY]
GROUPS={'A':'target + retained','B':'target + lost','C':'non-target + retained','D':'non-target + removed'}

def clean(v):
    if isinstance(v,dict):return {str(k):clean(x) for k,x in v.items()}
    if isinstance(v,(list,tuple,np.ndarray)):return [clean(x) for x in v]
    if isinstance(v,np.generic):return clean(v.item())
    if isinstance(v,float) and not np.isfinite(v):return None
    return v
def dump(name,v): (OUT/name).write_text(json.dumps(clean(v),indent=2,allow_nan=False)+'\n')
def tab(name,df):df.to_csv(OUT/'tables'/name,sep='\t',index=False)
def progress(msg):print(msg,flush=True)
def hvalues(x):
    if isinstance(x,h5py.Group):
        cats=hvalues(x['categories']);return np.array([cats[i] if i>=0 else None for i in x['codes'][:]],object)
    return x.asstr()[:] if h5py.check_string_dtype(x.dtype) is not None else x[:]

def metrics(y,keep):
    n=len(y); nt=int(y.sum()); nr=int(keep.sum()); tp=int((y&keep).sum())
    return dict(n_cells=n,target_n=nt,retained_n=nr,target_retained_n=tp,starting_target_fraction=nt/n if n else np.nan,post_gate_target_fraction=tp/nr if nr else np.nan,target_recovery=tp/nt if nt else np.nan,fold_enrichment=(tp/nr)/(nt/n) if nr and nt else np.nan,total_cell_yield=nr/n if n else np.nan)

def extract(cells,genes):
    frozen_columns=set(cells.columns)
    existing=json.loads((OLD/'cache/genes.json').read_text());ix={g:i for i,g in enumerate(existing)}
    expr=np.load(OLD/'cache/expression.npy',mmap_mode='r')
    missing=[]
    for gene in genes:
        if gene in cells:continue
        if gene in ix:cells[gene]=np.asarray(expr[:,ix[gene]])
        else:missing.append(gene)
    audit=[]
    with h5py.File(H30,'r') as f:
        allgenes=list(hvalues(f['var'][f['var'].attrs['_index']]))
        ids=hvalues(f['obs']['cell_id']);assert np.array_equal(ids,cells.cell_id.to_numpy())
        present=[g for g in missing if g in allgenes]
        cols=pd.Index(allgenes).get_indexer(present)
        if present:
            values=np.empty((len(cells),len(present)),np.float32)
            group=f['X'];ip=group['indptr'][:];n,p=group.attrs['shape']
            for start in range(0,n,2000):
                stop=min(start+2000,n);a,b=int(ip[start]),int(ip[stop])
                block=sparse.csr_matrix((group['data'][a:b],group['indices'][a:b],ip[start:stop+1]-a),shape=(stop-start,p))
                values[start:stop]=block[:,cols].toarray()
            cells=pd.concat([cells,pd.DataFrame(values,columns=present,index=cells.index)],axis=1)
        for g in genes:
            assayed=g in allgenes
            audit.append(dict(gene=g,available=g in cells,assayed_in_H5AD=assayed,frozen_unassayed_zero_placeholder=(g in cells and not assayed and np.all(cells[g]==0)),source='frozen_unassayed_zero_placeholder_preserved' if g in cells and not assayed else 'frozen_cells' if g in frozen_columns else 'frozen_expression_cache' if g in ix else 'H5AD_X_bounded_extraction' if assayed else 'absent_from_H5AD'))
    tab('gene_availability.tsv',pd.DataFrame(audit))
    return cells

def prepare():
    for p in [OUT,OUT/'tables',OUT/'provenance',OUT/'logs']:p.mkdir(parents=True,exist_ok=True)
    cells=pd.read_csv(OLD/'cells.tsv.gz',sep='\t')
    cfg=json.loads((OLD/'provenance/scoring_config.json').read_text())
    genes=sorted(set(REQUESTED+MARKERS+sum(SUPPLEMENTARY.values(),[])+sum(cfg['panels'].values(),[])+cfg['sst']))
    cells=extract(cells,genes)
    chosen=pd.read_csv(OLD/'tables/surface_gate_best_at_recovery_floors.tsv',sep='\t')
    chosen=chosen[(chosen.target=='top20')&(chosen.recovery_floor==.5)&(chosen.maximum_markers==3)].iloc[0]
    rules=json.loads(chosen.rules);assert chosen.logic=='AND';assert [c['gene'] for c in rules]==MARKERS
    clauses={c['gene']:(cells[c['gene']].to_numpy()>np.float64(c['threshold'])) if c['op']=='>' else (cells[c['gene']].to_numpy()<=np.float64(c['threshold'])) for c in rules}
    cells['retained']=np.logical_and.reduce(list(clauses.values()))
    threshold=float(cells.early_pv_score.quantile(.8));cells['target']=cells.early_pv_score>=threshold
    cells['group']=np.select([cells.target&cells.retained,cells.target&~cells.retained,~cells.target&cells.retained],['A','B','C'],default='D')
    cells['antecedent_percentile']=cells.early_pv_score.rank(method='average',pct=True)
    for m in MARKERS:cells['passes_'+m]=clauses[m]
    supp=[]
    for name,panel in SUPPLEMENTARY.items():
        available=[g for g in panel if g in cells];x=cells[available].to_numpy(float);mu=x.mean(0);sd=x.std(0);sd[sd==0]=1
        cells[f'supplementary_{name}_score']=np.clip((x-mu)/sd,-3,3).mean(1)
        supp.append(dict(program=name,requested_genes=panel,available_genes=available,means=mu,std_population=sd,definition='DIV30 gene z-scores clipped [-3,3], unweighted mean; supplementary context only, never target or gate',interpretation='Shared developmental context, not subtype-specific' if name=='nonmodule_pv_context' else 'Context only, not validated lineage classification'))
    metric=metrics(cells.target.to_numpy(),cells.retained.to_numpy())
    assert metric['n_cells']==90631 and metric['target_n']==18127 and metric['retained_n']==15781 and metric['target_retained_n']==9210
    assert np.isclose(metric['post_gate_target_fraction'],chosen.post_gate_target_fraction,atol=1e-9)
    dump('provenance/frozen_gate.json',dict(gate_id=chosen.gate_id,rules=rules,logic=chosen.logic,target='frozen early_pv_score >= original 80th percentile',target_threshold=threshold,metrics=metric,threshold_precision='np.float64 archive thresholds applied against expression vectors'))
    dump('provenance/supplementary_programs.json',supp)
    coef=pd.read_csv(OLD/'tables/sparse_model_coefficients.tsv',sep='\t').set_index('gene')
    independence=[]
    for g in genes:
        independence.append(dict(gene=g,in_direct_pv=g in cfg['pv'],in_direct_sst=g in cfg['sst'],in_reference_mapping_features=g in cfg['mapping_features'],sparse_coefficient=float(coef.loc[g,'coefficient']) if g in coef.index else 0.,in_explicit_context_panel=any(g in v for v in cfg['panels'].values()),in_fixed_gate=g in MARKERS))
    tab('gene_independence_audit.tsv',pd.DataFrame(independence))
    cells.to_csv(OUT/'cells.tsv.gz',sep='\t',index=False)
    shutil.copy2(__file__,OUT/'provenance'/Path(__file__).name)
    progress(f'Prepared {len(cells):,} cells; group counts {cells.group.value_counts().to_dict()}')
    return cells,genes,rules,metric

def stats_row(x):
    x=np.asarray(x,float);q=np.quantile(x,[.05,.1,.25,.5,.75,.9,.95])
    return dict(n=len(x),mean=float(x.mean()),std=float(x.std(ddof=1)),q05=q[0],q10=q[1],q25=q[2],median=q[3],q75=q[4],q90=q[5],q95=q[6],positive_fraction=float((x>0).mean()))

def contrast(x,y):
    x=np.asarray(x,float);y=np.asarray(y,float)
    u=stats.mannwhitneyu(x,y,alternative='two-sided',method='asymptotic')
    pooled=np.sqrt(((len(x)-1)*x.var(ddof=1)+(len(y)-1)*y.var(ddof=1))/(len(x)+len(y)-2))
    return dict(n_a=len(x),n_b=len(y),mean_a=x.mean(),mean_b=y.mean(),mean_difference=x.mean()-y.mean(),median_a=np.median(x),median_b=np.median(y),standardized_mean_difference=(x.mean()-y.mean())/pooled if pooled else 0.,probability_a_greater_with_half_ties=u.statistic/(len(x)*len(y)),rank_biserial=2*u.statistic/(len(x)*len(y))-1,positive_fraction_a=(x>0).mean(),positive_fraction_b=(y>0).mean(),positive_difference_pp=100*((x>0).mean()-(y>0).mean()),cell_level_mw_p_descriptive_only=u.pvalue)

def summarize_groups(cells,genes):
    masks={**{g:cells.group.eq(g).to_numpy() for g in GROUPS},'all':np.ones(len(cells),bool),'target':cells.target.to_numpy(),'non_target':~cells.target.to_numpy(),'retained':cells.retained.to_numpy(),'removed':~cells.retained.to_numpy()}
    programs=[p for p in PROGRAMS if p in cells]
    group_rows=[];gene_rows=[];program_rows=[]
    for label,mask in masks.items():
        part=cells.loc[mask]
        row=dict(group=label,group_label=GROUPS.get(label,label),n=len(part),target_fraction=part.target.mean(),fraction_all=mask.mean(),retained_fraction=part.retained.mean())
        for p in programs:row[p+'_mean']=part[p].mean()
        row['SLC6A1_positive_fraction']=(part.SLC6A1>0).mean();row['DCX_positive_fraction']=(part.DCX>0).mean()
        group_rows.append(row)
        for g in genes:
            if g in part:gene_rows.append(dict(group=label,gene=g,**stats_row(part[g])))
        for p in programs:program_rows.append(dict(group=label,program=p,**stats_row(part[p])))
    availability=pd.read_csv(OUT/'tables/gene_availability.tsv',sep='\t')
    tab('group_summary.tsv',pd.DataFrame(group_rows));tab('group_gene_summary.tsv',pd.DataFrame(gene_rows).merge(availability,on='gene',how='left'));tab('group_program_summary.tsv',pd.DataFrame(program_rows))
    contrasts=[('A','B'),('A','C'),('C','D'),('target','non_target'),('retained','removed')]
    for features,name in [(genes,'gene_contrasts.tsv'),(programs,'program_contrasts.tsv')]:
        rows=[]
        for a,b in contrasts:
            for feature in features:
                if feature in cells:rows.append(dict(comparison=f'{a}_vs_{b}',feature=feature,**contrast(cells.loc[masks[a],feature],cells.loc[masks[b],feature])))
        df=pd.DataFrame(rows)
        # Descriptive cell-level FDR, explicitly not independent biological evidence.
        for _,ix in df.groupby('comparison').groups.items():
            p=df.loc[ix,'cell_level_mw_p_descriptive_only'].to_numpy();order=np.argsort(p);adj=np.minimum.accumulate((p[order]*len(p)/(np.arange(len(p))+1))[::-1])[::-1];vals=np.empty(len(p));vals[order]=np.minimum(adj,1);df.loc[ix,'cell_level_bh_q_descriptive_only']=vals
        if name=='gene_contrasts.tsv':df=df.merge(availability.rename(columns={'gene':'feature'}),on='feature',how='left')
        tab(name,df)
    sample_rows=[]
    for (sample,group),part in cells.groupby(['sample','group']):
        row=dict(sample=sample,group=group,n=len(part),cell_line=part.cell_line.iloc[0],condition=part.condition.iloc[0])
        for p in programs+REQUESTED:row[p+'_mean']=part[p].mean()
        row['SLC6A1_positive_fraction']=(part.SLC6A1>0).mean();sample_rows.append(row)
    tab('group_by_sample.tsv',pd.DataFrame(sample_rows))
    # Same direction among six samples is descriptive: one sample per line/condition.
    rows=[]
    for a,b in contrasts:
        for feature in programs+REQUESTED:
            vals=[]
            for sample,ix in cells.groupby('sample').groups.items():
                ia=np.asarray(ix)[masks[a][ix]];ib=np.asarray(ix)[masks[b][ix]]
                if len(ia) and len(ib):vals.append((sample,cells.loc[ia,feature].mean()-cells.loc[ib,feature].mean()))
            rows.append(dict(comparison=f'{a}_vs_{b}',feature=feature,n_samples=len(vals),samples_positive=sum(v>0 for _,v in vals),sample_mean_difference_mean=np.mean([v for _,v in vals]),sample_mean_difference_min=min(v for _,v in vals),sample_mean_difference_max=max(v for _,v in vals),sample_differences=json.dumps(clean(dict(vals)))))
    tab('sample_effect_direction.tsv',pd.DataFrame(rows))
    return masks

def gate_dissection(cells,rules):
    masks={g:cells['passes_'+g].to_numpy() for g in MARKERS};target=cells.target.to_numpy();abl=[];index={}
    for n in range(0,4):
        for subset in itertools.combinations(MARKERS,n):
            keep=np.logical_and.reduce([masks[g] for g in subset]) if subset else np.ones(len(cells),bool)
            row=dict(markers='+'.join(subset) or 'none',n_markers=n,**metrics(target,keep));abl.append(row);index[frozenset(subset)]=row
    increments=[]
    for base in index:
        for marker in set(MARKERS)-set(base):
            before=index[base];after=index[base|{marker}]
            increments.append(dict(base_markers=before['markers'],added_marker=marker,result_markers=after['markers'],purity_gain_pp=100*(after['post_gate_target_fraction']-before['post_gate_target_fraction']),recovery_change_pp=100*(after['target_recovery']-before['target_recovery']),additional_target_lost=before['target_retained_n']-after['target_retained_n'],additional_non_target_removed=(before['retained_n']-before['target_retained_n'])-(after['retained_n']-after['target_retained_n']),retained_n_before=before['retained_n'],retained_n_after=after['retained_n']))
    tab('gate_ablations.tsv',pd.DataFrame(abl));tab('gate_incremental.tsv',pd.DataFrame(increments))
    def branchrow(label,keep):
        p=cells.loc[keep];return dict(branch=label,n=len(p),fraction_all=len(p)/len(cells),antecedent_fraction=p.target.mean(),pv_score_mean=p.direct_pv_score.mean(),sst_score_mean=p.direct_sst_score.mean(),antecedent_score_mean=p.early_pv_score.mean(),SLC6A1_positive_fraction=(p.SLC6A1>0).mean(),SLC6A1_mean=p.SLC6A1.mean(),mge_score_mean=p.mge_score.mean(),postmitotic_fraction=p.postmitotic.mean(),maturation_score_mean=p.maturation_score.mean(),cycling_fraction=p.cycling.mean())
    truth=[]
    for pattern in itertools.product([False,True],repeat=3):
        keep=np.logical_and.reduce([masks[g] if b else ~masks[g] for g,b in zip(MARKERS,pattern)])
        truth.append(dict(**{f'passes_{g}':b for g,b in zip(MARKERS,pattern)},**branchrow(' '.join(f'{g}={b}' for g,b in zip(MARKERS,pattern)),keep)))
    tab('gate_truth_patterns.tsv',pd.DataFrame(truth))
    branches=[('GRIA2 fails',~masks['GRIA2']),('GRIA2 passes / OPCML fails',masks['GRIA2']&~masks['OPCML']),('GRIA2 passes / OPCML passes / NOTCH1 fails',masks['GRIA2']&masks['OPCML']&~masks['NOTCH1']),('GRIA2 passes / OPCML passes / NOTCH1 passes (retained)',cells.retained.to_numpy())]
    tab('gate_decision_branches.tsv',pd.DataFrame([branchrow(l,m) for l,m in branches]))
    assert sum(m.sum() for _,m in branches)==len(cells)
    conditional=[]
    for m1,m2 in itertools.combinations(MARKERS,2):
        for b1,b2 in itertools.product([False,True],repeat=2):
            keep=(masks[m1] if b1 else ~masks[m1])&(masks[m2] if b2 else ~masks[m2]);row=branchrow(f'{m1}={b1} / {m2}={b2}',keep);row.update(marker_a=m1,passes_a=b1,marker_b=m2,passes_b=b2)
            for m in MARKERS:row[m+'_mean']=cells.loc[keep,m].mean()
            conditional.append(row)
    tab('marker_conditional.tsv',pd.DataFrame(conditional))
    perf=[]
    for axes in [['sample'],['cell_line'],['condition'],['cell_line','condition']]:
        for labels,part in cells.groupby(axes,observed=True):
            if not isinstance(labels,tuple):labels=(labels,)
            row=dict(grouping='+'.join(axes),group=' | '.join(str(x) for x in labels),**dict(zip(axes,labels)),**metrics(part.target.to_numpy(),part.retained.to_numpy()),fraction_all_retained=part.retained.sum()/cells.retained.sum(),fraction_all_retained_targets=(part.target&part.retained).sum()/(cells.target&cells.retained).sum())
            perf.append(row)
    tab('line_condition_performance.tsv',pd.DataFrame(perf))
    omission=[]
    for sample in sorted(cells['sample'].unique()):
        p=cells.loc[cells['sample']!=sample]
        omission.append(dict(excluded_sample=sample,**metrics(p.target.to_numpy(),p.retained.to_numpy()),interpretation='Descriptive removal of one sample with target and gate frozen; NOT refitting, gate selection, or independent validation'))
    tab('leave_sample_out_fixed_gate_composition.tsv',pd.DataFrame(omission))

def score_and_slc(cells,masks):
    summary=[];tops=[]
    for label,mask in masks.items():summary.append(dict(group=label,**stats_row(cells.loc[mask,'early_pv_score'])))
    tab('continuous_score_summary.tsv',pd.DataFrame(summary))
    for q in [.05,.1,.2,.3,.4]:
        threshold=float(cells.early_pv_score.quantile(1-q));top=cells.early_pv_score.to_numpy()>=threshold
        tops.append(dict(nominal_top_fraction=q,threshold=threshold,**metrics(top,cells.retained.to_numpy())))
    tab('score_top_fraction_recovery.tsv',pd.DataFrame(tops))
    slc=[]
    for label,mask in masks.items():slc.append(dict(group=label,**stats_row(cells.loc[mask,'SLC6A1'])))
    tab('slc6a1_group_summary.tsv',pd.DataFrame(slc))
    rows=[]
    axes=['early_pv_score','direct_pv_score','direct_sst_score','neighbor_pv_score','sparse_pv_score','mge_score','maturation_score','progression_score','cycle_score','postmitotic','cycling','progenitor_score','supplementary_gaba_score','SLC6A1','ERBB4','LHX6','DCX','n_genes','log_counts']
    for marker in ['SLC6A1']+MARKERS:
        for axis in axes:
            if marker==axis:continue
            for subset,mask in [('all',masks['all']),('postmitotic',cells.postmitotic.to_numpy()),('target',masks['target']),('non_target',masks['non_target'])]:
                a,b=cells.loc[mask,marker],cells.loc[mask,axis]
                rho,p=stats.spearmanr(a,b) if a.nunique()>1 and b.nunique()>1 else (np.nan,np.nan)
                rows.append(dict(marker=marker,axis=axis,subset=subset,n=int(mask.sum()),spearman=rho,cell_level_p_descriptive_only=p))
    tab('axis_correlations.tsv',pd.DataFrame(rows))
    # Quartiles keep ties; zeros are retained as their own meaningful group.
    positive=cells.SLC6A1>0;upper=float(cells.loc[positive,'SLC6A1'].quantile(.75))
    bins={'undetected':~positive,'detected_below_positive_q75':positive&(cells.SLC6A1<upper),'detected_at_or_above_positive_q75':positive&(cells.SLC6A1>=upper)}
    slcrows=[]
    for label,mask in bins.items():
        p=cells.loc[mask];slcrows.append(dict(slc6a1_stratum=label,n=len(p),fraction_all=len(p)/len(cells),target_fraction=p.target.mean(),retained_fraction=p.retained.mean(),antecedent_mean=p.early_pv_score.mean(),direct_pv_mean=p.direct_pv_score.mean(),sst_mean=p.direct_sst_score.mean(),maturation_mean=p.maturation_score.mean(),postmitotic_fraction=p.postmitotic.mean(),threshold_positive_q75=upper))
    tab('slc6a1_strata.tsv',pd.DataFrame(slcrows))
    # Operational evidence flags overlap; do not assert exclusive lineages.
    fp=masks['C'];fn=masks['B'];evidence={
       'score_percentile_70_to_80':(cells.antecedent_percentile>=.7)&~cells.target,
       'score_percentile_60_to_80':(cells.antecedent_percentile>=.6)&~cells.target,
       'score_percentile_below_50':cells.antecedent_percentile<.5,
       'postmitotic':cells.postmitotic,
       'operational_progenitor':cells.progenitor,
       'cycling':cells.cycling,
       'sst_program_above_all_cell_q75':cells.direct_sst_score>=cells.direct_sst_score.quantile(.75),
       'SST_detected':cells.SST>0,
       'LHX6_or_ERBB4_detected':(cells.LHX6>0)|(cells.ERBB4>0),
       'NKX2-1_detected':cells['NKX2-1']>0,
       'GAD1_or_GAD2_detected':(cells.GAD1>0)|(cells.GAD2>0),
       'glutamate_transporter_detected':(cells.SLC17A6>0)|(cells.SLC17A7>0),
       'astrocyte_context_2_of_GFAP_AQP4_S100B':(cells[['GFAP','AQP4','S100B']]>0).sum(axis=1)>=2,
       'oligodendroglial_context_2_of_OLIG1_OLIG2_PDGFRA_SOX10':(cells[['OLIG1','OLIG2','PDGFRA','SOX10']]>0).sum(axis=1)>=2,
       'CHAT_or_SLC18A3_detected':(cells.CHAT>0)|(cells.SLC18A3>0),
       'out_of_DIV90_reference':cells.mapping_out_of_reference,
    }
    flags=[]
    for label,mask in masks.items():
        for criterion,yes in evidence.items():flags.append(dict(group=label,criterion=criterion,n_group=int(mask.sum()),n_flagged=int((mask&yes).sum()),fraction=float(yes[mask].mean()),interpretation='Overlapping descriptive evidence; RNA nondetection is not lineage absence; not validated subtype annotation'))
    tab('false_positive_negative_evidence.tsv',pd.DataFrame(flags))

def evaluate_binary(y,p):
    return dict(roc_auc=roc_auc_score(y,p) if len(np.unique(y))>1 else np.nan,average_precision=average_precision_score(y,p),log_loss=log_loss(y,p,labels=[0,1]),brier_score=brier_score_loss(y,p),observed_fraction=y.mean(),mean_prediction=p.mean())

def additional_coupling(cells):
    """Frozen-score SLC6A1 coupling; conditional ranks do not redefine scores."""
    axes=['early_pv_score','direct_pv_score','rank_pv_score','neighbor_pv_score','sparse_pv_score','direct_sst_score','maturation_score','supplementary_gaba_score']
    controls=['progenitor_score','maturation_score','postmitotic','DCX','cycle_score','cycling','supplementary_gaba_score','log_counts','n_genes','percent_mt']
    rows=[]
    groups=[('all','all',cells.index)]
    for grouping in ['sample','cell_line','condition']:
        groups.extend((grouping,label,ix) for label,ix in cells.groupby(grouping).groups.items())
    for grouping,label,ix in groups:
        p=cells.loc[ix];x=stats.rankdata(p.SLC6A1)
        for axis in axes:
            y=stats.rankdata(p[axis]);used=[c for c in controls if c!=axis]
            design=np.column_stack([stats.rankdata(p[c]) for c in used]);design=StandardScaler().fit_transform(design)
            # Residual Spearman is an adjusted association diagnostic, not a
            # causal effect or independent development/fate endpoint.
            design=np.column_stack([np.ones(len(p)),design]);rx=x-design@np.linalg.lstsq(design,x,rcond=None)[0];ry=y-design@np.linalg.lstsq(design,y,rcond=None)[0]
            row=dict(grouping=grouping,group=label,axis=axis,n=len(p),SLC6A1_positive_fraction=(p.SLC6A1>0).mean(),SLC6A1_mean=p.SLC6A1.mean(),spearman=stats.spearmanr(x,y).statistic,partial_spearman_linear_rank_residuals=stats.pearsonr(rx,ry).statistic,controls=';'.join(used))
            rows.append(row)
    tab('slc6a1_adjusted_by_sample_line_condition.tsv',pd.DataFrame(rows))
    rows=[]
    for a,b in itertools.combinations(MARKERS,2):
        for label,part in [('all',cells),('target',cells[cells.target]),('non_target',cells[~cells.target]),('postmitotic',cells[cells.postmitotic])]:rows.append(dict(marker_a=a,marker_b=b,subset=label,n=len(part),spearman=stats.spearmanr(part[a],part[b]).statistic,both_detected_fraction=((part[a]>0)&(part[b]>0)).mean()))
    tab('surface_marker_pairwise.tsv',pd.DataFrame(rows))

def models(cells):
    core=['progenitor_score','maturation_score','postmitotic','DCX']
    extended=core+['cycle_score','s_score','g2m_score','cycling','supplementary_gaba_score','log_counts','n_genes','percent_mt']
    configs={
      'A_core':(core,False),
      'B_core_plus_PV':(core+['direct_pv_score'],False),
      'A_extended_generic':(extended,False),
      'B_extended_plus_PV':(extended+['direct_pv_score'],False),
      'A_flexible_generic':(extended,True),
      'B_flexible_plus_PV':(extended+['direct_pv_score'],True),
      'B_extended_nonmodule_context':(extended+['supplementary_nonmodule_pv_context_score'],False),
    }
    ydict={'antecedent_target':cells.target.to_numpy(int),'fixed_gate_retained':cells.retained.to_numpy(int)}
    sample=cells['sample'].to_numpy();folds=sorted(set(sample));fold_rows=[];coef_rows=[];preds=[];overall=[];calibration=[]
    for outcome,y in ydict.items():
        for name,(features,flexible) in configs.items():
            x=cells[features].to_numpy(float);oof=np.zeros(len(cells));null=np.zeros(len(cells))
            for fold in folds:
                train=sample!=fold;test=~train;scaler=StandardScaler().fit(x[train]);a=scaler.transform(x[train]);b=scaler.transform(x[test])
                if flexible:
                    # Univariate cubic spline effects provide a flexible, additive
                    # generic-maturation sensitivity; not a new optimized gate.
                    spline=SplineTransformer(n_knots=4,degree=3,include_bias=False,knots='quantile',extrapolation='linear').fit(a)
                    a=spline.transform(a);b=spline.transform(b)
                model=LogisticRegression(C=1.0,max_iter=800,tol=1e-6,solver='lbfgs').fit(a,y[train]);p=model.predict_proba(b)[:,1]
                oof[test]=p;null[test]=y[train].mean()
                row=dict(outcome=outcome,model=name,held_sample=fold,n_train=int(train.sum()),n_test=int(test.sum()),n_iter=int(model.n_iter_[0]),**evaluate_binary(y[test],p));fold_rows.append(row)
                if not flexible:
                    for feature,coefficient in zip(features,model.coef_[0]):coef_rows.append(dict(outcome=outcome,model=name,held_sample=fold,feature=feature,coefficient_per_training_sd=coefficient,odds_ratio_per_training_sd=np.exp(coefficient)))
            met=evaluate_binary(y,oof);null_loss=log_loss(y,null,labels=[0,1]);met.update(outcome=outcome,model=name,features=';'.join(features),flexible_additive_splines=flexible,null_log_loss=null_loss,predictive_log_loss_reduction_fraction=1-met['log_loss']/null_loss)
            overall.append(met);preds.append(pd.DataFrame({'cell_id':cells.cell_id,'sample':sample,'outcome':outcome,'model':name,'observed':y,'oof_probability':oof}))
            bins=pd.qcut(oof,10,duplicates='drop')
            for interval,ix in pd.Series(np.arange(len(y))).groupby(bins,observed=True):
                ii=ix.to_numpy();calibration.append(dict(outcome=outcome,model=name,bin=str(interval),n=len(ii),mean_probability=oof[ii].mean(),observed_fraction=y[ii].mean()))
            progress(f'Model {outcome} {name}: AUC={met["roc_auc"]:.3f}, AP={met["average_precision"]:.3f}, logloss={met["log_loss"]:.3f}')
    tab('model_oof_metrics.tsv',pd.DataFrame(overall));tab('model_fold_metrics.tsv',pd.DataFrame(fold_rows));tab('model_coefficients.tsv',pd.DataFrame(coef_rows));tab('model_predictions.tsv.gz',pd.concat(preds,ignore_index=True));tab('model_calibration.tsv',pd.DataFrame(calibration))
    contrasts=[('A_core','B_core_plus_PV'),('A_extended_generic','B_extended_plus_PV'),('A_flexible_generic','B_flexible_plus_PV'),('A_extended_generic','B_extended_nonmodule_context')]
    rows=[]
    for granularity,table,keys in [('pooled',pd.DataFrame(overall),['outcome']),('sample',pd.DataFrame(fold_rows),['outcome','held_sample'])]:
        for labels,part in table.groupby(keys):
            labels=labels if isinstance(labels,tuple) else (labels,);d=part.set_index('model')
            for aa,bb in contrasts:
                a,b=d.loc[aa],d.loc[bb];rows.append(dict(granularity=granularity,**dict(zip(keys,labels)),base_model=aa,augmented_model=bb,roc_auc_improvement=b.roc_auc-a.roc_auc,average_precision_improvement=b.average_precision-a.average_precision,log_loss_reduction=a.log_loss-b.log_loss,brier_reduction=a.brier_score-b.brier_score,additional_information_bits_per_cell=(a.log_loss-b.log_loss)/np.log(2)))
    tab('model_incremental.tsv',pd.DataFrame(rows))
    neighbor_models_and_design(cells,configs)

def neighbor_models_and_design(cells,configs):
    sample=cells['sample'].to_numpy();folds=sorted(set(sample))
    # Predict neighbor resemblance from disjoint query scoring genes as context;
    # the neighbor labels still come from PV-scored DIV90 reference cells.
    reg=[];regpred=[]
    y=cells.neighbor_pv_score.to_numpy()
    for name in ['A_extended_generic','B_extended_plus_PV','A_flexible_generic','B_flexible_plus_PV']:
        features,flexible=configs[name];x=cells[features].to_numpy(float);pred=np.zeros(len(cells))
        for fold in folds:
            train=sample!=fold;test=~train;scaler=StandardScaler().fit(x[train]);a=scaler.transform(x[train]);b=scaler.transform(x[test])
            if flexible:
                spline=SplineTransformer(n_knots=4,degree=3,include_bias=False,knots='quantile',extrapolation='linear').fit(a);a=spline.transform(a);b=spline.transform(b)
            model=Ridge(alpha=10).fit(a,y[train]);pred[test]=model.predict(b)
        reg.append(dict(outcome='frozen_neighbor_pv_score',model=name,oof_r2=r2_score(y,pred),oof_mse=mean_squared_error(y,pred),oof_spearman=stats.spearmanr(y,pred).statistic));regpred.append(pd.DataFrame({'cell_id':cells.cell_id,'model':name,'observed_neighbor_pv':y,'oof_prediction':pred}))
    tab('neighbor_resemblance_model_metrics.tsv',pd.DataFrame(reg));tab('neighbor_resemblance_model_predictions.tsv.gz',pd.concat(regpred,ignore_index=True))
    dump('provenance/model_design.json',dict(configs=configs,folds=folds,validation='leave one original DIV30 sample out; scaler, spline knots and model fitted on training cells only; frozen target is never recomputed within folds',logistic='L2 C=1, lbfgs max_iter800 tol1e-6, no hyperparameter tuning',calibration='held-sample out-of-fold probability deciles; intercepts and slopes not recalibrated on test cells',important_limit='PV direct-score improvement for antecedent classification is partly mathematical circularity: direct and rank PV gene scores are two of four target components. Fixed gate was selected against full-data target, so held-sample gate prediction is descriptive biological explanation, NOT new prospective validation or a replacement of the preserved held-sample gate selection result.',context_program='MAF/MAFB/ERBB4/LHX6/SOX6 are excluded from original mapping features and explicit PV module but are broad developmental/MGE context, not a PV-specific independent endpoint',neighbor_outcome='Neighbor resemblance uses frozen nonsurface mapping features excluding direct PV genes, but DIV90 neighbor labels still derive from the same PV program; association does not validate cell fate.'))

def requested_model_audit(cells):
    """Exact user formulas, with DCX-augmented core retained as sensitivity."""
    configs={'A_requested':['progenitor_score','maturation_score','postmitotic'],'B_requested_plus_PV':['progenitor_score','maturation_score','postmitotic','direct_pv_score']}
    sample=cells['sample'].to_numpy();folds=sorted(set(sample));fold_rows=[];coef_rows=[];preds=[];overall=[]
    for outcome,y in {'antecedent_target':cells.target.to_numpy(int),'fixed_gate_retained':cells.retained.to_numpy(int)}.items():
        for name,features in configs.items():
            x=cells[features].to_numpy(float);oof=np.zeros(len(cells));null=np.zeros(len(cells))
            for fold in folds:
                train=sample!=fold;test=~train;scaler=StandardScaler().fit(x[train]);a=scaler.transform(x[train]);b=scaler.transform(x[test])
                model=LogisticRegression(C=1,max_iter=800,tol=1e-6,solver='lbfgs').fit(a,y[train]);p=model.predict_proba(b)[:,1];oof[test]=p;null[test]=y[train].mean()
                fold_rows.append(dict(outcome=outcome,model=name,held_sample=fold,n_train=int(train.sum()),n_test=int(test.sum()),n_iter=int(model.n_iter_[0]),**evaluate_binary(y[test],p)))
                for f,co in zip(features,model.coef_[0]):coef_rows.append(dict(outcome=outcome,model=name,held_sample=fold,feature=f,coefficient_per_training_sd=co,odds_ratio_per_training_sd=np.exp(co)))
            met=evaluate_binary(y,oof);null_loss=log_loss(y,null,labels=[0,1]);met.update(outcome=outcome,model=name,features=';'.join(features),flexible_additive_splines=False,null_log_loss=null_loss,predictive_log_loss_reduction_fraction=1-met['log_loss']/null_loss);overall.append(met)
            preds.append(pd.DataFrame({'cell_id':cells.cell_id,'sample':sample,'outcome':outcome,'model':name,'observed':y,'oof_probability':oof}))
            progress(f'Exact requested {outcome} {name}: AUC={met["roc_auc"]:.4f}, AP={met["average_precision"]:.4f}')
    for name,rows in [('model_oof_metrics.tsv',overall),('model_fold_metrics.tsv',fold_rows),('model_coefficients.tsv',coef_rows)]:
        old=pd.read_csv(OUT/'tables'/name,sep='\t');old=old[~old.model.isin(configs)];tab(name,pd.concat([old,pd.DataFrame(rows)],ignore_index=True))
    requested=pd.concat(preds,ignore_index=True);tab('requested_model_predictions.tsv.gz',requested)
    cal=[]
    for (outcome,model),p in requested.groupby(['outcome','model']):
        for interval,part in p.groupby(pd.qcut(p.oof_probability,10,duplicates='drop'),observed=True):cal.append(dict(outcome=outcome,model=model,bin=str(interval),n=len(part),mean_probability=part.oof_probability.mean(),observed_fraction=part.observed.mean()))
    oldcal=pd.read_csv(OUT/'tables/model_calibration.tsv',sep='\t');tab('model_calibration.tsv',pd.concat([oldcal[~oldcal.model.isin(configs)],pd.DataFrame(cal)],ignore_index=True))
    rows=[]
    for granularity,table,keys in [('pooled',pd.DataFrame(overall),['outcome']),('sample',pd.DataFrame(fold_rows),['outcome','held_sample'])]:
        for labels,part in table.groupby(keys):
            labels=labels if isinstance(labels,tuple) else (labels,);d=part.set_index('model');a,b=d.loc['A_requested'],d.loc['B_requested_plus_PV']
            rows.append(dict(granularity=granularity,**dict(zip(keys,labels)),base_model='A_requested',augmented_model='B_requested_plus_PV',roc_auc_improvement=b.roc_auc-a.roc_auc,average_precision_improvement=b.average_precision-a.average_precision,log_loss_reduction=a.log_loss-b.log_loss,brier_reduction=a.brier_score-b.brier_score,additional_information_bits_per_cell=(a.log_loss-b.log_loss)/np.log(2)))
    old=pd.read_csv(OUT/'tables/model_incremental.tsv',sep='\t');tab('model_incremental.tsv',pd.concat([old[old.base_model!='A_requested'],pd.DataFrame(rows)],ignore_index=True))
    design=json.loads((OUT/'provenance/model_design.json').read_text());design['exact_requested_formulas']=configs;dump('provenance/model_design.json',design)
    # Every subgroup prediction was made with its original sample held out;
    # splitting the predictions by line is not a new leave-line-out validation.
    pred=pd.concat([pd.read_csv(OUT/'tables/model_predictions.tsv.gz',sep='\t'),requested],ignore_index=True)
    pred=pred.merge(cells[['cell_id','cell_line','condition']],on='cell_id',how='left',validate='many_to_one')
    subgroup=[]
    for grouping in ['cell_line','condition']:
        for (outcome,model,group),p in pred.groupby(['outcome','model',grouping]):
            subgroup.append(dict(grouping=grouping,group=group,outcome=outcome,model=model,n=len(p),**evaluate_binary(p.observed.to_numpy(),p.oof_probability.to_numpy())))
    sub=pd.DataFrame(subgroup);tab('model_subgroup_metrics.tsv',sub)
    rows=[]
    pairs=[('A_requested','B_requested_plus_PV'),('A_core','B_core_plus_PV'),('A_extended_generic','B_extended_plus_PV'),('A_flexible_generic','B_flexible_plus_PV')]
    for (grouping,group,outcome),part in sub.groupby(['grouping','group','outcome']):
        d=part.set_index('model')
        for aa,bb in pairs:
            a,b=d.loc[aa],d.loc[bb];rows.append(dict(grouping=grouping,group=group,outcome=outcome,base_model=aa,augmented_model=bb,roc_auc_improvement=b.roc_auc-a.roc_auc,average_precision_improvement=b.average_precision-a.average_precision,log_loss_reduction=a.log_loss-b.log_loss,brier_reduction=a.brier_score-b.brier_score))
    tab('model_subgroup_incremental.tsv',pd.DataFrame(rows))

def findings(cells,metric):
    group=pd.read_csv(OUT/'tables/group_summary.tsv',sep='\t').set_index('group')
    inc=pd.read_csv(OUT/'tables/model_incremental.tsv',sep='\t');mods=pd.read_csv(OUT/'tables/model_oof_metrics.tsv',sep='\t')
    flags=pd.read_csv(OUT/'tables/false_positive_negative_evidence.tsv',sep='\t')
    slc=pd.read_csv(OUT/'tables/slc6a1_group_summary.tsv',sep='\t').set_index('group')
    cor=pd.read_csv(OUT/'tables/axis_correlations.tsv',sep='\t');line=pd.read_csv(OUT/'tables/line_condition_performance.tsv',sep='\t')
    lines=['# Frozen DIV30 within-time biological identity findings','',f'Exact four groups: A {int(group.loc["A","n"]):,}; B {int(group.loc["B","n"]):,}; C {int(group.loc["C","n"]):,}; D {int(group.loc["D","n"]):,}. Original target and fixed thresholds are unchanged.','',
    '**Independence correction:** SLC6A1 is absent from the explicit MEF2C/KCNC1/KCNC2/TAC1 module, but is present in the frozen 1,000-gene mapping and has a nonzero sparse coefficient (+0.010314168). Its association with the frozen consensus therefore cannot be called wholly independent. No score was changed to remove it. SLC6A1 is excluded from every new supplementary GABA and PV-context program.','',
    '**Models:** Leave-one-sample-out models explain the existing labels. They do not validate future fate or redo gate discovery. The direct PV program is already part of the antecedent definition; adding it to antecedent models partly tests construction, not independent specificity. Six samples include one sample for each of three lines × two conditions.','']
    for outcome in ['antecedent_target','fixed_gate_retained']:
        pairs=[('A_core','B_core_plus_PV'),('A_extended_generic','B_extended_plus_PV'),('A_flexible_generic','B_flexible_plus_PV')]
        if (mods.model=='A_requested').any():pairs=[('A_requested','B_requested_plus_PV')]+pairs
        for aa,bb in pairs:
            a=mods[(mods.outcome==outcome)&(mods.model==aa)].iloc[0];b=mods[(mods.outcome==outcome)&(mods.model==bb)].iloc[0]
            lines.append(f'- {outcome}, {aa} → {bb}: AUC {a.roc_auc:.4f} → {b.roc_auc:.4f}; AP {a.average_precision:.4f} → {b.average_precision:.4f}; logloss {a.log_loss:.4f} → {b.log_loss:.4f}; generic model logloss reduction vs held-sample null {a.predictive_log_loss_reduction_fraction:.1%}.')
    lines+=['','SLC6A1 detection: '+', '.join(f'{label} {slc.loc[label,"positive_fraction"]:.1%}' for label in ['all','target','non_target','retained','removed','A','B','C','D'])+'.','']
    for axis in ['early_pv_score','direct_pv_score','neighbor_pv_score','maturation_score']:
        row=cor[(cor.marker=='SLC6A1')&(cor.axis==axis)&(cor.subset=='all')].iloc[0];lines.append(f'- SLC6A1 vs {axis}: Spearman {row.spearman:.3f}.')
    lines+=['','Retained non-target evidence (overlapping flags; no imposed lineage labels):']
    for flag in ['score_percentile_70_to_80','score_percentile_60_to_80','score_percentile_below_50','postmitotic','cycling','operational_progenitor','SST_detected','sst_program_above_all_cell_q75','LHX6_or_ERBB4_detected','GAD1_or_GAD2_detected','glutamate_transporter_detected','out_of_DIV90_reference']:
        row=flags[(flags.group=='C')&(flags.criterion==flag)].iloc[0];lines.append(f'- {flag}: {row.fraction:.1%} ({row.n_flagged:,}/{row.n_group:,}).')
    lines+=['','Line contributions to retained targets:']
    for _,r in line[line.grouping=='cell_line'].iterrows():lines.append(f'- {r["group"]}: {r.fraction_all_retained_targets:.1%} of retained targets; starting {r.starting_target_fraction:.1%}, post-gate {r.post_gate_target_fraction:.1%}, recovery {r.target_recovery:.1%}.')
    lines+=['','Every gene and program comparison reports detection, continuous effect sizes, ranks and six-sample effect directions. Cell-level Mann–Whitney and BH columns are descriptive only; their small values are not independent biological-replicate significance. New supplementary stress/GABA/lineage-context programs do not redefine targets. Absence of RNA is not absence of protein or lineage.','', 'TUBB3 is not assayed in the source H5AD feature catalogue; its inherited all-zero column and seven-gene maturation-score denominator are preserved exactly. Gene tables flag this placeholder rather than interpreting it as measured nondetection. All explicitly requested identity genes are assayed.','', 'Missing metadata-listed samples 7–9 remain missing; the frozen six-sample coverage and prior held-sample gate-selection results are preserved.']
    (OUT/'findings.md').write_text('\n'.join(lines)+'\n')
    dump('summary.json',dict(metrics=metric,groups=group.reset_index().to_dict('records'),slc6a1=slc.reset_index().to_dict('records'),models=mods.to_dict('records'),model_incremental=inc.query("granularity == 'pooled'").to_dict('records'),prior_results_read_only=True,large_gate_search_rerun=False,target_changed=False))

def validate_outputs(cells):
    truth=pd.read_csv(OUT/'tables/gate_truth_patterns.tsv',sep='\t');branches=pd.read_csv(OUT/'tables/gate_decision_branches.tsv',sep='\t')
    assert truth.n.sum()==len(cells) and branches.n.sum()==len(cells)
    assert cells.group.eq('A').sum()==9210 and cells.group.eq('B').sum()==8917 and cells.group.eq('C').sum()==6571 and cells.group.eq('D').sum()==65933
    audit=pd.read_csv(OUT/'tables/gene_availability.tsv',sep='\t').set_index('gene');assert audit.loc[REQUESTED,'assayed_in_H5AD'].all()
    assert not audit.loc['TUBB3','assayed_in_H5AD'] and audit.loc['TUBB3','frozen_unassayed_zero_placeholder']
    top=pd.read_csv(OUT/'tables/score_top_fraction_recovery.tsv',sep='\t');assert np.all(np.diff(top.target_retained_n)>=0)
    folds=pd.read_csv(OUT/'tables/model_fold_metrics.tsv',sep='\t');assert folds.groupby(['outcome','model']).n_test.sum().eq(len(cells)).all()
    assert folds.n_iter.max()<800
    assert pd.read_csv(OUT/'tables/model_oof_metrics.tsv',sep='\t').set_index(['outcome','model']).index.is_unique
    dump('provenance/validation.json',dict(status='passed',n_cells=len(cells),four_groups_exact=True,eight_marker_truth_patterns_partition_all_cells=True,decision_branches_partition_all_cells=True,all_named_identity_genes_assayed=True,TUBB3_unassayed_placeholder_flagged=True,top_fraction_retained_counts_monotone=True,every_model_predicts_each_cell_once_in_held_sample=True,logistic_iterations_below_cap=True,large_gate_search_rerun=False,target_definitions_changed=False))

def main():
    if '--prepared' in sys.argv or '--finish' in sys.argv:
        cells=pd.read_csv(OUT/'cells.tsv.gz',sep='\t');genes=pd.read_csv(OUT/'tables/gene_availability.tsv',sep='\t').gene.tolist()
        frozen=json.loads((OUT/'provenance/frozen_gate.json').read_text());rules=frozen['rules'];metric=frozen['metrics']
        assert metrics(cells.target.to_numpy(),cells.retained.to_numpy())==metric
        shutil.copy2(__file__,OUT/'provenance'/Path(__file__).name)
    else:cells,genes,rules,metric=prepare()
    if '--finish' not in sys.argv:
        masks=summarize_groups(cells,genes);progress('Group and gene summaries complete')
        gate_dissection(cells,rules);score_and_slc(cells,masks);additional_coupling(cells);progress('Fixed-marker ablations, continuous scores, SLC6A1 and lineage context complete')
    with threadpool_limits(limits=2):
        if '--finish' in sys.argv:
            rows=pd.read_csv(OUT/'tables/model_oof_metrics.tsv',sep='\t').drop_duplicates('model')
            configs={r.model:(r.features.split(';'),bool(r.flexible_additive_splines)) for _,r in rows.iterrows()}
            neighbor_models_and_design(cells,configs)
        else:models(cells)
        requested_model_audit(cells)
    findings(cells,metric);validate_outputs(cells);progress(f'Complete: {OUT}')

if __name__=='__main__':main()
