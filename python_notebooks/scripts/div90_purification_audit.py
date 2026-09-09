#!/usr/bin/env python3
"""Purification-first, descriptive DIV90 audit of frozen phase-two cells.

No biological fitting, label replacement, or source-output writes.  Every old
candidate is reconstructed; additional signed grid pairs and bounded hierarchical
triplets extend the old preservation-oriented search. Run with --stage all, or
--stage summarize to reproduce summaries from the completed audit HDF5.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import time

import h5py
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr
from threadpoolctl import threadpool_limits

ROOT = Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results')
SOURCE = ROOT / 'div90_hypergate_sst_pv_phase2'
OUT = ROOT / 'div90_purification_audit_v1'
SOURCE_FILES = ['depletion_all_candidate_metrics.h5', 'depletion_supplemental_candidate_metrics.h5']
STATES = ['PV-biased', 'PV/SST hybrid', 'SST-biased', 'unresolved/immature']
COUNT_NAMES = ['n', 'pv', 'dual', 'sst', 'unresolved', 'hc_dual', 'hc_loose',
               'hc_strict', 'hc_no_qc', 'continuous10', 'continuous20', 'continuous30',
               'rankdiff10', 'rankdiff20', 'rankdiff30']
RECOVERIES = [0., .1, .2, .3, .4, .5, .6, .7, .8]


def log(s):
    print(time.strftime('%H:%M:%S'), s, flush=True)


def clean(x):
    if isinstance(x, dict): return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple, np.ndarray)): return [clean(v) for v in x]
    if isinstance(x, (np.integer,)): return int(x)
    if isinstance(x, (np.bool_,)): return bool(x)
    if isinstance(x, (float, np.floating)): return float(x) if np.isfinite(x) else None
    return x


def dump(path, x):
    path.write_text(json.dumps(clean(x), indent=2, allow_nan=False) + '\n')


def div(a, b, zero=np.nan):
    a, b = np.broadcast_arrays(np.asarray(a, float), np.asarray(b, float))
    return np.divide(a, b, out=np.full(a.shape, zero), where=b != 0)


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8 << 20), b''): h.update(block)
    return h.hexdigest()


def boolcol(c, key):
    return c[key].astype(str).str.lower().isin(['true', '1']).to_numpy()


def prepare_targets(c):
    state = c.state.to_numpy()
    n = len(c)
    pv = c.pv_score.to_numpy(float)
    sst = c.sst_score.to_numpy(float)
    pr = (rankdata(pv, method='average') - .5) / n
    sr = (rankdata(sst, method='average') - .5) / n
    preference = np.minimum(pr, 1 - sr)
    rankdiff = pr - sr
    qc = ~(boolcol(c, 'predicted_doublet') | boolcol(c, 'high_complexity_tail'))
    robust = pv >= np.quantile(pv, .65)
    multigene = (c.pv_n_detected.to_numpy() >= 2) & (c.sst_n_detected_without_SST.to_numpy() >= 2)
    dual = state == STATES[1]
    hc = dual & multigene & robust & qc
    loose = dual & (c.pv_n_detected.to_numpy() >= 2) & (c.sst_n_detected_without_SST.to_numpy() >= 1) & qc
    strict = dual & (c.pv_n_detected.to_numpy() >= 3) & (c.sst_n_detected_without_SST.to_numpy() >= 2) & (pv >= np.quantile(pv, .75)) & qc
    masks = [np.ones(n, bool)] + [state == s for s in STATES] + [hc, loose, strict, dual & multigene & robust]
    extra = c[['cell_id','sample','cell_line','condition','state','pv_score','sst_score',
               'pv_n_detected','sst_n_detected_without_SST','predicted_doublet','high_complexity_tail']].copy()
    extra['pv_rank_percentile'], extra['sst_rank_percentile'] = pr, sr
    extra['continuous_preference'] = preference
    extra['rank_difference_sensitivity'] = rankdiff
    cuts = {}
    for label, v in [('continuous', preference), ('rankdiff', rankdiff)]:
        for pct in [10,20,30]:
            cutoff = np.quantile(v, 1-pct/100)
            m = v >= cutoff  # All ties at the cutoff are included and disclosed.
            masks.append(m)
            cuts[label+str(pct)] = {'cutoff':cutoff, 'n':int(m.sum()), 'fraction':m.mean(), 'requested_fraction':pct/100}
    W = np.asarray(masks, np.float32).T
    for j, name in enumerate(COUNT_NAMES[5:]): extra[name] = W[:,j+5].astype(bool)
    extra.to_csv(OUT/'tables/target_cell_assignments.tsv.gz', sep='\t', index=False)
    method = {
        'primary_hc_dual_definition':'Frozen dual-high AND >=2 of MEF2C/KCNC1/KCNC2/TAC1 detected AND >=2 of SATB1/NR2F2/CDK14/CACNG3 detected AND PV score >= pooled 65th percentile AND no saved predicted_doublet AND no saved within-sample top-5%-UMI high_complexity_tail.',
        'hc_label_caution':'Operational high-confidence PV-supported dual-high; not independently validated biological hybrids. Low/negative cells remain in all gate denominators.',
        'hc_pv_score_cutoff':np.quantile(pv,.65), 'hc_strict_pv_score_cutoff':np.quantile(pv,.75),
        'hc_loose_definition':'Frozen dual-high, >=2 PV genes, >=1 non-SST SST-module gene, same QC exclusions; original median-high PV condition only.',
        'hc_strict_definition':'Frozen dual-high, >=3 PV genes, >=2 non-SST SST-module genes, PV score >=75th pooled percentile, same QC exclusions.',
        'hc_no_qc_definition':'Primary HC rule before QC exclusions.',
        'continuous_definition':'min(midrank_percentile(PV score), 1-midrank_percentile(SST score)); a maximin joint high-PV/low-SST preference without assuming commensurate raw score scales.',
        'continuous_ties':'Include all cells at each quantile cutoff; achieved starting fraction may differ from nominal percentage.',
        'continuous_sensitivity':'midrank_percentile(PV score)-midrank_percentile(SST score); corresponding top 10/20/30% tie-inclusive subsets.',
        'score_scaling':{'pv_sd':pv.std(), 'sst_sd':sst.std(), 'pv_iqr':np.subtract(*np.quantile(pv,[.75,.25])), 'sst_iqr':np.subtract(*np.quantile(sst,[.75,.25])), 'raw_difference_vs_rank_difference_spearman':spearmanr(pv-sst,rankdiff).statistic},
        'continuous_targets':cuts, 'baseline_counts':dict(zip(COUNT_NAMES,W.sum(0).astype(int))),
        'biological_scope':'Frozen 4768 DIV90 entry cells and original independent scores; no identity, UMAP, trajectory, or gene-score recomputation.',
    }
    dump(OUT/'target_definitions.json',method)
    return W, np.column_stack([pv,sst]), method


def expr_from_old(rule):
    terms = [{'gene':x['gene'],'op':x['op'],'threshold':float(x['threshold'])} for x in rule['rules']]
    e = {'logic':rule['logic'], 'terms':terms}
    return {'not':e} if rule['action']=='remove' else e


def evaluate(c, e):
    if 'not' in e: return ~evaluate(c,e['not'])
    if 'gene' in e:
        x = c[e['gene']].to_numpy(float)
        return {'>':np.greater,'<=':np.less_equal,'>=':np.greater_equal,'<':np.less}[e['op']](x,e['threshold'])
    terms = [evaluate(c,t) for t in e['terms']]
    if not terms: return np.ones(len(c),bool) if e['logic']=='AND' else np.zeros(len(c),bool)
    return np.logical_and.reduce(terms) if e['logic']=='AND' else np.logical_or.reduce(terms)


def label(e):
    if 'not' in e: return 'NOT ('+label(e['not'])+')'
    if 'gene' in e: return f"{e['gene']} {e['op']} {e['threshold']:.17g}"
    return '('+(' '+e['logic']+' ').join(label(t) for t in e['terms'])+')' if e['terms'] else 'ALL'


def genes_in(e):
    if 'not' in e:return genes_in(e['not'])
    if 'gene' in e:return {e['gene']}
    return set().union(*(genes_in(x) for x in e['terms'])) if e['terms'] else set()


def get_rule(h, group, i):
    g=h[group]
    kind=g.attrs['kind']
    atoms=json.loads(h.attrs['atoms'])
    if kind=='rules':return json.loads(g['rules_json'].asstr()[i])
    if kind=='single':
        return {'gene':json.loads(h.attrs['genes'])[int(g['gene_index'][i])], 'op':str(g.attrs['op']), 'threshold':float(g['threshold'][i])}
    if kind=='pair':
        a,b=atoms[int(g['a'][i])].copy(),atoms[int(g['b'][i])].copy()
        a['op']='>' if int(g.attrs['sign_a']) else '<='
        b['op']='>' if int(g.attrs['sign_b']) else '<='
        e={'logic':'AND','terms':[a,b]}
        return {'not':e} if int(g.attrs['complement']) else e
    if kind=='triplet':
        p=json.loads(g.attrs['parents'])[int(g['parent'][i])]
        a=atoms[int(g['atom'][i])].copy()
        a['op']='>' if int(g['sign'][i]) else '<='
        return {'logic':'OR' if int(g['union'][i]) else 'AND','terms':[p,a]}
    raise ValueError(kind)


def derived(counts, scores, total):
    x=np.asarray(counts,dtype=float)
    d={'retained_n':x[:,0], 'total_yield':x[:,0]/total[0]}
    for j,name in enumerate(COUNT_NAMES[1:],1):
        d[name+'_fraction']=div(x[:,j],x[:,0])
        d[name+'_recovery']=div(x[:,j],total[j])
    d['pv_fold_enrichment']=div(d['pv_fraction'],total[1]/total[0])
    d['pv_all_dual_fraction']=div(x[:,1]+x[:,2],x[:,0])
    d['pv_all_dual_recovery']=div(x[:,1]+x[:,2],total[1]+total[2])
    d['pv_all_dual_fold_enrichment']=div(d['pv_all_dual_fraction'],(total[1]+total[2])/total[0])
    d['supported_fraction']=div(x[:,1]+x[:,5],x[:,0])
    d['supported_recovery']=div(x[:,1]+x[:,5],total[1]+total[5])
    d['supported_fold_enrichment']=div(d['supported_fraction'],(total[1]+total[5])/total[0])
    for name in ['loose','strict','no_qc']:
        j=COUNT_NAMES.index('hc_'+name)
        d['supported_'+name+'_fraction']=div(x[:,1]+x[:,j],x[:,0])
        d['supported_'+name+'_recovery']=div(x[:,1]+x[:,j],total[1]+total[j])
        d['supported_'+name+'_fold_enrichment']=div(d['supported_'+name+'_fraction'],(total[1]+total[j])/total[0])
    d['sst_fold_depletion']=div(total[3]/total[0],d['sst_fraction'],zero=np.inf)
    d['pv_sst_ratio']=div(x[:,1],x[:,3],zero=np.inf)
    d['pv_sst_ratio'][(x[:,1]==0)&(x[:,3]==0)]=np.nan
    d['sst_removal']=1-d['sst_recovery']
    d['pv_loss']=1-d['pv_recovery']
    d['hc_dual_loss']=1-d['hc_dual_recovery']
    d['pv_score_mean']=div(scores[:,0],x[:,0])
    d['sst_score_mean']=div(scores[:,1],x[:,0])
    for name in ['continuous10','continuous20','continuous30','rankdiff10','rankdiff20','rankdiff30']:
        d[name+'_fold_enrichment']=div(d[name+'_fraction'],total[COUNT_NAMES.index(name)]/total[0])
    return d


def write_group(h,name,counts,scores,**defs):
    counts=np.rint(counts).astype(np.int32)
    totals=np.asarray(json.loads(h.attrs['totals']))
    assert np.all(counts>=0) and np.all(counts<=totals)
    assert np.array_equal(counts[:,0],counts[:,1:5].sum(1))
    assert np.all(counts[:,5]<=counts[:,2])
    g=h.create_group(name)
    g.create_dataset('counts',data=counts.astype(np.uint16),compression='gzip',compression_opts=2,shuffle=True)
    g.create_dataset('score_sums',data=scores,compression='gzip',compression_opts=2,shuffle=True)
    for k,v in defs.items():
        if np.isscalar(v):g.attrs[k]=v
        else:g.create_dataset(k,data=np.asarray(v),compression='gzip',compression_opts=2,shuffle=True)
    # Every candidate gets persisted composition/recovery/fold metrics, not merely shortlisted rows.
    keys=list(derived(counts[:1],scores[:1],totals))
    g.attrs['metric_columns']=json.dumps(keys)
    ds=g.create_dataset('metrics',shape=(len(counts),len(keys)),dtype='f8',chunks=(min(8192,len(counts)),len(keys)),compression='gzip',compression_opts=2,shuffle=True)
    for start in range(0,len(counts),65536):
        d=derived(counts[start:start+65536],scores[start:start+65536],totals)
        ds[start:start+65536]=np.column_stack([d[k] for k in keys])
    h.flush()
    log(f'{name}: {len(counts):,} candidates')
    return g


def iter_count_chunks(h, seed_only=False):
    """Bound ranking memory independently of total search size."""
    for gn,g in h.items():
        if not isinstance(g,h5py.Group) or 'counts' not in g:continue
        if seed_only and int(g.attrs.get('n_markers',0))!=2:continue
        for start in range(0,len(g['counts']),65536):
            yield gn,g,start,g['counts'][start:start+65536].astype(np.int32),g['score_sums'][start:start+65536]


def summarize(h, c, W, S, seed_only=False):
    total=W.sum(0).astype(int)
    # Objectives use descriptive fractions and counts. No F1 or significance ranking.
    objectives=[]
    for target in ['pv','supported','supported_loose','supported_strict','supported_no_qc','continuous10','continuous20','continuous30','rankdiff10','rankdiff20','rankdiff30']:
        for r in RECOVERIES:
            objectives.append((f'{target}_purity_recovery_{r:.1f}',target+'_fraction',target+'_recovery',r,None))
    for r in RECOVERIES: objectives.append((f'pv_ratio_recovery_{r:.1f}','pv_sst_ratio','pv_recovery',r,None))
    for j in [1,2,3]:
        for r in RECOVERIES:objectives.append((f'markers_{j}_pv_recovery_{r:.1f}','pv_fraction','pv_recovery',r,j))
    for r in [.2,.3,.4,.5,.6,.7,.8,1.]:objectives.append((f'sst_removal_{r:.1f}','pv_recovery','sst_removal',r,None))
    objectives += [('minimum_sst_fraction','negative_sst_fraction','pv_recovery',0.,None)]
    thresholds=[]
    for target in ['pv','supported','pv_all_dual']:
        for v in [1.5,2.]:thresholds.append((f'{target}_enrichment_{v}',target+'_fold_enrichment',v,target if target!='pv_all_dual' else 'pv'))
        for v in [.65,.7,.75]:thresholds.append((f'{target}_composition_{v}',target+'_fraction',v,target if target!='pv_all_dual' else 'pv'))
    for v in [.2,.15,.1]:thresholds.append((f'sst_below_{v}','negative_sst_fraction',-v,'pv'))
    best={}
    fronts={k:{} for k in ['pv','supported','continuous10','continuous20','continuous30','sst_removal']}
    group_inventory=[]
    for gn,g,start,co,sc in iter_count_chunks(h,seed_only):
        nmark=int(g.attrs.get('n_markers',0))
        d=derived(co,sc,total)
        d['negative_sst_fraction']=-d['sst_fraction']
        valid=co[:,0]>0
        if start==0:group_inventory.append({'group':gn,'n':len(g['counts']),'n_markers':nmark,'scope':str(g.attrs.get('scope',''))})
        for name,obj,con,minimum,markers in objectives:
            if markers and nmark!=markers:continue
            # Subtraction-defined removal rates can lie one ULP below an exact
            # rational boundary, e.g. 1 - 860/1075 for 20% SST removal.
            ok=valid & (d[con]>=minimum-1e-12) & ~np.isnan(d[obj])
            ix=np.flatnonzero(ok)
            if not len(ix):continue
            # Objective, PV recovery, lower SST, larger capture; final exact row ID tie-break.
            vals=d[obj][ix]; ix=ix[vals==np.max(vals)]
            idx=ix[np.lexsort((ix,-co[ix,0],d['sst_fraction'][ix],-d['pv_recovery'][ix]))[0]]
            key=(float(d[obj][idx]),float(d['pv_recovery'][idx]),-float(d['sst_fraction'][idx]),int(co[idx,0]),-nmark)
            if name not in best or key>best[name][0]:best[name]=(key,gn,start+int(idx))
        for name,obj,value,target in ([] if seed_only else thresholds):
            ok=valid & ((d[obj]>value) if obj=='negative_sst_fraction' else (d[obj]>=value))
            ix=np.flatnonzero(ok)
            if not len(ix):continue
            rec=target+'_recovery'; v=d[rec][ix];ix=ix[v==v.max()]
            idx=ix[np.lexsort((ix,-co[ix,0],-d[obj][ix]))[0]]
            key=(float(d[rec][idx]),float(d[obj][idx]),int(co[idx,0]),-nmark)
            if name not in best or key>best[name][0]:best[name]=(key,gn,start+int(idx))
        for target,positions in fronts.items():
            if target=='sst_removal':
                xx=total[3]-co[:,3]; yy=d['pv_recovery']
            else:
                xx=co[:,1]+co[:,5] if target=='supported' else co[:,COUNT_NAMES.index(target)]
                yy=d[target+'_fraction']
            ix=np.flatnonzero(valid & np.isfinite(yy))
            # One best purity/loss point per exact recovered-count value.
            order=ix[np.lexsort((-d['pv_recovery'][ix],-yy[ix],xx[ix]))]
            _,first=np.unique(xx[order],return_index=True)
            for i in order[first]:
                x=int(xx[i]);key=(float(yy[i]),float(d['pv_recovery'][i]),int(co[i,0]),-nmark)
                if x not in positions or key>positions[x][0]:positions[x]=(key,gn,start+int(i))
    if seed_only:
        # Purity-oriented pair seeds across every requested recovery and state target,
        # plus ratio and removal objectives; exact phenotype duplicates removed.
        candidates=[];seen=set()
        for name,(_,gn,i) in best.items():
            if name.startswith('rankdiff') or name.startswith('markers_'):continue
            e=get_rule(h,gn,i);key=json.dumps(e,sort_keys=True)
            if len(genes_in(e))==2 and key not in seen:
                seen.add(key);candidates.append({'objective':name,'group':gn,'row':i,'rule':e})
        for target,points in fronts.items():
            frontier=[];best_y=-np.inf
            for x,(key,gn,i) in sorted(points.items(),reverse=True):
                if key[0]>best_y:
                    frontier.append((gn,i));best_y=key[0]
            for ix in np.unique(np.linspace(0,len(frontier)-1,min(8,len(frontier))).astype(int)):
                gn,i=frontier[ix];e=get_rule(h,gn,i);key=json.dumps(e,sort_keys=True)
                if len(genes_in(e))==2 and key not in seen:
                    seen.add(key);candidates.append({'objective':'bounded_full_pareto_'+target,'group':gn,'row':i,'rule':e})
        return candidates
    def record(gn,i,objective=''):
        co=h[gn]['counts'][i:i+1].astype(int); sc=h[gn]['score_sums'][i:i+1]
        e=get_rule(h,gn,i)
        r={'objective':objective,'group':gn,'source_row':i,'n_markers':len(genes_in(e)), 'exact_retention_rule':label(e),'rule_json':json.dumps(e,separators=(',',':'))}
        for j,k in enumerate(COUNT_NAMES):
            r['retained_'+k+'_n']=int(co[0,j]);r['baseline_'+k+'_n']=int(total[j]);r['baseline_'+k+'_fraction']=total[j]/total[0]
        r.update({k:float(v[0]) for k,v in derived(co,sc,total).items()})
        r['baseline_supported_fraction']=(total[1]+total[5])/total[0]
        r['baseline_pv_all_dual_fraction']=(total[1]+total[2])/total[0]
        r['baseline_pv_score_mean']=S[:,0].mean();r['baseline_sst_score_mean']=S[:,1].mean()
        r['pv_enrichment_class']='WEAK' if r['pv_fold_enrichment']<1.5 else '>=2x' if r['pv_fold_enrichment']>=2 else '>=1.5x'
        r['supported_enrichment_class']='WEAK' if r['supported_fold_enrichment']<1.5 else '>=2x' if r['supported_fold_enrichment']>=2 else '>=1.5x'
        return r
    rows=[record(gn,i,name) for name,(_,gn,i) in best.items()]
    pd.DataFrame(rows).to_csv(OUT/'tables/recovery_sweeps.tsv',sep='\t',index=False)
    pd.DataFrame(group_inventory).to_csv(OUT/'tables/search_inventory.tsv',sep='\t',index=False)
    critical_map={'A':'pv_purity_recovery_0.0','B':'pv_purity_recovery_0.5','C':'pv_purity_recovery_0.7','D':'supported_purity_recovery_0.0','E':'minimum_sst_fraction','F':'pv_ratio_recovery_0.0','G':'markers_1_pv_recovery_0.1','H':'markers_2_pv_recovery_0.1','I':'markers_3_pv_recovery_0.1'}
    critical=[]
    for letter,name in critical_map.items():
        if name in best:
            _,gn,i=best[name];critical.append({'criterion':letter,**record(gn,i,name)})
    pd.DataFrame(critical).to_csv(OUT/'tables/critical_best_gates.tsv',sep='\t',index=False)
    trows=[]
    for name,obj,value,target in thresholds:
        if name in best:
            _,gn,i=best[name];r=record(gn,i,name);r['attained']=True
        else:r={'objective':name,'attained':False}
        r.update(threshold_metric=obj,threshold=value,recovery_optimized=target, useful_recovery_definition='>=10% of specified target',attained_at_10pct_recovery=bool(r.get(target+'_recovery',0)>=.1))
        trows.append(r)
    pd.DataFrame(trows).to_csv(OUT/'tables/threshold_attainment.tsv',sep='\t',index=False)
    frontier_rows=[]
    for target,points in fronts.items():
        frontier=[];best_y=-np.inf
        for x,(key,gn,i) in sorted(points.items(),reverse=True):
            if key[0]>best_y:
                frontier.append(record(gn,i,target));best_y=key[0]
        pd.DataFrame(frontier).to_csv(OUT/f'tables/pareto_{target}.tsv',sep='\t',index=False)
        frontier_rows.extend(r for r in frontier if r['pv_recovery']>=.1)
    # Fixed selected gates: robust counts, target-definition sensitivity, sample/condition composition.
    selected={(r['group'],r['source_row']) for r in rows+critical+frontier_rows}
    bygroup=[];validation=[];registry=[];selected_masks=[]
    for gn,i in sorted(selected):
        e=get_rule(h,gn,i);mask=evaluate(c,e)
        actual=np.rint(mask.astype(float)@W).astype(int)
        expected=h[gn]['counts'][i].astype(int)
        if not np.array_equal(actual,expected):raise AssertionError((gn,i,actual,expected))
        if not np.allclose(mask@S,h[gn]['score_sums'][i],atol=1e-8):raise AssertionError('score mismatch')
        validation.append({'group':gn,'row':i,'counts_match':True,'score_sums_match':True})
        registry.append({'group':gn,'row':i,'retention_rule':e})
        selected_masks.append(mask)
        for groupby in ['sample','cell_line','condition']:
            for value,ix in c.groupby(groupby,sort=False).indices.items():
                t=W[ix].sum(0).astype(int);cc=np.rint(mask[ix]@W[ix]).astype(int)[None,:]
                ds=derived(cc,(mask[ix]@S[ix])[None,:],t)
                r={'group':gn,'row':i,'groupby':groupby,'stratum':str(value),'baseline_n':len(ix),'retained_n':int(mask[ix].sum())}
                r.update({k:float(v[0]) for k,v in ds.items()})
                bygroup.append(r)
    pd.DataFrame(bygroup).to_csv(OUT/'tables/selected_gate_fixed_sample_evaluations.tsv',sep='\t',index=False)
    dump(OUT/'selected_rules.json',registry)
    np.savez_compressed(OUT/'tables/selected_gate_retention_masks.npz',cell_id=c.cell_id.to_numpy(dtype=str),masks=np.asarray(selected_masks),rule_group=np.asarray([r['group'] for r in registry]),rule_row=np.asarray([r['row'] for r in registry]))
    dump(OUT/'provenance/selected_mask_validation.json',validation)
    summary={'total_evaluated_rows':sum(x['n'] for x in group_inventory),'legacy_evaluated_rows':sum(x['n'] for x in group_inventory if x['group'].startswith('legacy_')),'critical':critical,'thresholds':trows,'selected_masks_independently_reproduced':len(validation),'baseline_counts':dict(zip(COUNT_NAMES,total)), 'critical_G_H_I_definition':'Maximum PV purity at >=10% PV recovery for each exact marker count. Absolute maxima and all recovery thresholds are also in recovery_sweeps.tsv.', 'scope_caution':'Observed best in frozen surface feature set and explicitly evaluated boundaries; expanded pairs use the frozen grid, triplets are bounded purity/Pareto-seeded extensions. No theoretical global purification bound, held-out significance, fate validation, or protein sorting claim.'}
    dump(OUT/'audit_summary.json',summary)
    log(f'Summaries complete: {summary["total_evaluated_rows"]:,} rows, {len(validation)} selected masks verified')
    return summary


def extend_triplets(h,c,W,S):
    atoms=json.loads(h.attrs['atoms']);totals=W.sum(0).astype(int);scoretot=S.sum(0)
    B=np.asarray([evaluate(c,a) for a in atoms],np.float32)
    singles=np.rint(B@W).astype(int);single_scores=B@S
    seeds=summarize(h,c,W,S,seed_only=True)
    dump(OUT/'triplet_seed_rules.json',seeds)
    log(f'Extending {len(seeds)} distinct purity/recovery/removal pair seeds with all remaining marker atoms')
    allco=[];allsc=[];pi=[];ai=[];si=[];ui=[]
    for parent_idx,seed in enumerate(seeds):
        p=evaluate(c,seed['rule']);used=genes_in(seed['rule'])
        ix=np.asarray([j for j,a in enumerate(atoms) if a['gene'] not in used],int)
        pm=p.astype(float)@W;ps=p.astype(float)@S
        cross=np.rint((B[ix]*p)@W).astype(int);cross_score=(B[ix]*p)@S
        for sign in [0,1]:
            ac=singles[ix] if sign else totals-singles[ix];asc=single_scores[ix] if sign else scoretot-single_scores[ix]
            ic=cross if sign else pm-cross;iss=cross_score if sign else ps-cross_score
            for union in [0,1]:
                allco.append((pm+ac-ic if union else ic).astype(np.uint16));allsc.append(ps+asc-iss if union else iss)
                pi.extend([parent_idx]*len(ix));ai.extend(ix);si.extend([sign]*len(ix));ui.extend([union]*len(ix))
    write_group(h,'expanded_hierarchical_triplets',np.concatenate(allco),np.concatenate(allsc),kind='triplet',parents=json.dumps([s['rule'] for s in seeds]),parent=np.asarray(pi,np.uint16),atom=np.asarray(ai,np.uint16),sign=np.asarray(si,np.uint8),union=np.asarray(ui,np.uint8),n_markers=3,scope='bounded purity/recovery/removal-seeded AND/OR extension with either sign of every third surface marker on frozen grid')


def run():
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['all','summarize','extend'],default='all');args=parser.parse_args()
    for x in ['tables','provenance','logs','cache']: (OUT/x).mkdir(parents=True,exist_ok=True)
    c=pd.read_csv(SOURCE/'cells.tsv.gz',sep='\t',float_precision='round_trip')
    assert len(c)==4768 and not c.cell_id.duplicated().any()
    W,S,method=prepare_targets(c)
    dest=OUT/'tables/all_purification_candidates.h5'
    if args.stage=='summarize':
        with h5py.File(dest,'r') as h:summarize(h,c,W,S)
        return
    if args.stage=='extend':
        with h5py.File(dest,'r+') as h:
            assert 'expanded_hierarchical_triplets' not in h
            extend_triplets(h,c,W,S)
            summarize(h,c,W,S)
        before=json.loads((OUT/'provenance/source_hashes_before.json').read_text())
        after={p:sha(Path(p)) for p in before}
        assert before==after
        dump(OUT/'provenance/source_preservation.json',{'all_inputs_unchanged':True,'before':before,'after':after})
        return
    if dest.exists():raise RuntimeError(f'Preserving existing audit HDF5: {dest}; use --stage summarize or a new version.')
    sources=[SOURCE/'tables'/f for f in SOURCE_FILES]+[SOURCE/'cells.tsv.gz',SOURCE/'state_parameters.json',SOURCE/'hybrid_validation.json',SOURCE/'gate_summary.json']
    before={str(p):sha(p) for p in sources};dump(OUT/'provenance/source_hashes_before.json',before)
    with h5py.File(sources[0],'r') as old:
        atoms=json.loads(old.attrs['pair_threshold_definitions']);genes=json.loads(old.attrs['genes'])
        aa=old['all_pairs_and/threshold_a_index'][:];bb=old['all_pairs_and/threshold_b_index'][:]
    B=np.asarray([evaluate(c,a) for a in atoms],np.float32)
    totals=W.sum(0).astype(int);scoretot=S.sum(0)
    singles=np.rint(B@W).astype(int);single_scores=B@S
    log(f'{len(genes)} surface genes, {len(atoms)} frozen grid atoms; HC dual n={totals[5]}')
    inter=[]
    for j in range(W.shape[1]):
        sub=B[:,W[:,j].astype(bool)]
        mm=np.rint(sub@sub.T).astype(np.int32)
        inter.append(mm[aa,bb])
    pair_inter=np.column_stack(inter);del inter
    pair_score=np.column_stack([((B*S[:,j])@B.T)[aa,bb] for j in range(S.shape[1])])
    log('Exact pair intersections reconstructed for all count and continuous channels')
    with h5py.File(dest,'w') as h:
        h.attrs['count_columns']=json.dumps(COUNT_NAMES);h.attrs['totals']=json.dumps(totals.tolist());h.attrs['genes']=json.dumps(genes);h.attrs['atoms']=json.dumps(atoms);h.attrs['score_columns']=json.dumps(['pv_score','sst_score']);h.attrs['target_definitions']=json.dumps(clean(method))
        audits=[]
        for source_idx,path in enumerate(sources[:2]):
            with h5py.File(path,'r') as old:
                for key,g in old.items():
                    names=json.loads(g.attrs['metric_columns']);nm=1 if key=='exact_single' else 3 if 'third' in key else 2
                    defs={'scope':'legacy','n_markers':nm,'source_file':str(path),'source_group':key}
                    if key.startswith('all_pairs'):
                        if key.endswith('_and'):
                            co=totals-pair_inter;sc=scoretot-pair_score
                            defs.update(kind='pair',a=aa,b=bb,sign_a=1,sign_b=1,complement=1)
                        else:
                            co=totals-singles[aa]-singles[bb]+pair_inter;sc=scoretot-single_scores[aa]-single_scores[bb]+pair_score
                            defs.update(kind='pair',a=aa,b=bb,sign_a=0,sign_b=0,complement=0)
                    elif key=='exact_single':
                        gi=g['gene_index'][:];th=g['threshold'][:]
                        co=np.empty((len(gi),W.shape[1]),int);sc=np.empty((len(gi),2),float)
                        for gene_idx,gene in enumerate(genes):
                            ix=np.flatnonzero(gi==gene_idx);xx=c[gene].to_numpy(float);order=np.argsort(xx,kind='stable');pos=np.searchsorted(xx[order],th[ix],side='right')
                            co[ix]=np.r_[np.zeros((1,W.shape[1])),W[order].cumsum(0)][pos]
                            sc[ix]=np.r_[np.zeros((1,2)),S[order].cumsum(0)][pos]
                        defs.update(kind='single',gene_index=gi,threshold=th,op='<=')
                    else:
                        rr=[expr_from_old(json.loads(s)) for s in g['rules_json'].asstr()[:]]
                        co=np.empty((len(rr),W.shape[1]),int);sc=np.empty((len(rr),2),float)
                        for start in range(0,len(rr),2048):
                            masks=np.asarray([evaluate(c,e) for e in rr[start:start+2048]],np.float32)
                            co[start:start+2048]=np.rint(masks@W).astype(int);sc[start:start+2048]=masks@S
                        defs.update(kind='rules',rules_json=np.asarray([json.dumps(r,separators=(',',':')) for r in rr],dtype=h5py.string_dtype()))
                    # Compare reconstructed total and all four categorical counts for EVERY legacy row.
                    countix=[names.index('retained_n'),names.index('retained_pv_n'),names.index('retained_hybrid_n'),names.index('retained_sst_n'),names.index('retained_unresolved_n')]
                    for start in range(0,len(co),32768):
                        expected=np.asarray(g['metrics'][start:start+32768])[:,countix].astype(int)
                        if not np.array_equal(co[start:start+32768,:5],expected):raise AssertionError(f'Legacy reconstruction mismatch {path}/{key}, starting row {start}')
                    name=f'legacy_{source_idx}_{key}'
                    write_group(h,name,co,sc,**defs)
                    audits.append({'source':str(path),'group':key,'rows':len(co),'all_categorical_counts_match':True})
        dump(OUT/'provenance/all_legacy_counts_validation.json',{'groups':audits,'total':sum(x['rows'] for x in audits)})
        assert sum(x['rows'] for x in audits)==2885531
        log('ALL 2,885,531 legacy rows exactly reproduced')
        # Signed single complements add positive enrichment at every observed boundary.
        g=h['legacy_0_exact_single']
        write_group(h,'expanded_positive_singles',totals-g['counts'][:],scoretot-g['score_sums'][:],kind='single',gene_index=g['gene_index'][:],threshold=g['threshold'][:],op='>',n_markers=1,scope='all exact single complements')
        # Add the six missing truth patterns; original rows cover the other two.
        for sa,sb in [(1,1),(1,0),(0,1),(0,0)]:
            if sa and sb:pc,ps=pair_inter,pair_score
            elif sa:pc,ps=singles[aa]-pair_inter,single_scores[aa]-pair_score
            elif sb:pc,ps=singles[bb]-pair_inter,single_scores[bb]-pair_score
            else:pc,ps=totals-singles[aa]-singles[bb]+pair_inter,scoretot-single_scores[aa]-single_scores[bb]+pair_score
            for comp in [0,1]:
                if (sa,sb,comp) in [(1,1,1),(0,0,0)]:continue
                write_group(h,f'expanded_pairs_{sa}{sb}_complement{comp}',totals-pc if comp else pc,scoretot-ps if comp else ps,kind='pair',a=aa,b=bb,sign_a=sa,sign_b=sb,complement=comp,n_markers=2,scope='all signed pairs on frozen 5-threshold per-gene grid')
        # Preserve original R fits and frozen comparator explicitly, without refitting.
        oldsum=json.loads((SOURCE/'gate_summary.json').read_text())
        bench=[]
        for key in ['phase1_reference','experimental_depletion']:
            gid=oldsum[key]['gate_id'];bench.append({'name':key,'rule':expr_from_old(oldsum['rules'][gid])})
        rfits=pd.read_csv(SOURCE/'tables/hypergate_four_state_evaluations.tsv',sep='\t')
        for r in rfits.itertuples():bench.append({'name':r.job_id,'rule':expr_from_old(dict(action=r.action,logic=r.logic,rules=json.loads(r.rules)))})
        for nm in [1,2,3]:
            rr=[r for r in bench if len(genes_in(r['rule']))==nm]
            if not rr:continue
            masks=np.asarray([evaluate(c,r['rule']) for r in rr],np.float32)
            write_group(h,f'frozen_benchmarks_{nm}',masks@W,masks@S,kind='rules',rules_json=np.asarray([json.dumps(r['rule'],separators=(',',':')) for r in rr],dtype=h5py.string_dtype()),n_markers=nm,scope='frozen actual R fits and historical comparator rules')
        dump(OUT/'frozen_comparator_rules.json',bench[:2])
        del pair_inter,pair_score,B,co,sc,pc,ps
        extend_triplets(h,c,W,S)
        summarize(h,c,W,S)
    after={str(p):sha(p) for p in sources}
    assert before==after
    dump(OUT/'provenance/source_preservation.json',{'all_inputs_unchanged':before==after,'before':before,'after':after})
    log('Completed; source hashes unchanged')


if __name__=='__main__':
    with threadpool_limits(limits=4):run()
