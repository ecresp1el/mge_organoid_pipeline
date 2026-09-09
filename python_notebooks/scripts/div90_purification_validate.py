#!/usr/bin/env python3
"""Direct-mask checks and comparator metrics after the completed DIV90 audit."""
import json
import sys

import h5py
import numpy as np
import pandas as pd

from div90_purification_audit import OUT, SOURCE, COUNT_NAMES, get_rule, evaluate, derived, dump, label


def main():
    c=pd.read_csv(SOURCE/'cells.tsv.gz',sep='\t',float_precision='round_trip')
    targets=pd.read_csv(OUT/'tables/target_cell_assignments.tsv.gz',sep='\t',float_precision='round_trip')
    assert c.cell_id.equals(targets.cell_id)
    state=c.state.to_numpy()
    columns=[np.ones(len(c),bool)]+[state==s for s in ['PV-biased','PV/SST hybrid','SST-biased','unresolved/immature']]
    columns.extend(targets[k].astype(str).str.lower().isin(['true','1']).to_numpy() for k in COUNT_NAMES[5:])
    W=np.column_stack(columns).astype(np.int32)
    S=c[['pv_score','sst_score']].to_numpy(float)
    totals=W.sum(0)
    rng=np.random.default_rng(900209)
    sample_rows=[]
    if '--fix-zero-pv-ratios' in sys.argv:
        fixes=[]
        # The first running audit process predated the final 0/POSITIVE -> 0
        # correction. This changes only that arithmetic field in NEW outputs;
        # boundaries, masks, counts, all other metrics, and upstream files stay fixed.
        with h5py.File(OUT/'tables/all_purification_candidates.h5','r+') as h:
            for gn,g in h.items():
                if not isinstance(g,h5py.Group) or 'counts' not in g:continue
                col=json.loads(g.attrs['metric_columns']).index('pv_sst_ratio')
                fixed=0
                for start in range(0,len(g['counts']),8192):
                    co=g['counts'][start:start+8192]
                    ix=(co[:,1]==0)&(co[:,3]>0)
                    if not ix.any():continue
                    met=g['metrics'][start:start+8192]
                    wrong=ix & ~np.equal(met[:,col],0)
                    if wrong.any():
                        met[wrong,col]=0.;g['metrics'][start:start+8192]=met;fixed+=int(wrong.sum())
                fixes.append({'group':gn,'zero_pv_positive_sst_ratios_set_to_zero':fixed})
        dump(OUT/'provenance/zero_pv_ratio_arithmetic_correction.json',fixes)
    with h5py.File(OUT/'tables/all_purification_candidates.h5','r') as h:
        for gn,g in h.items():
            if not isinstance(g,h5py.Group) or 'counts' not in g:continue
            n=len(g['counts']);ix=np.unique(np.r_[0,n-1,rng.choice(n,min(n,64),replace=False)])
            for i in ix:
                rule=get_rule(h,gn,int(i));mask=evaluate(c,rule)
                obs=mask@W;expected=g['counts'][i].astype(int)
                assert np.array_equal(obs,expected),(gn,int(i),obs,expected)
                score=mask@S
                assert np.allclose(score,g['score_sums'][i],atol=1e-8),(gn,int(i),'scores')
                sample_rows.append({'group':gn,'row':int(i),'retained_n':int(mask.sum()),'all_count_channels_match':True,'score_sums_match':True})
    dump(OUT/'provenance/sampled_direct_mask_validation.json',{'seed':900209,'n_checked':len(sample_rows),'all_passed':True,'checks':sample_rows})
    benchmarks=json.loads((OUT/'frozen_comparator_rules.json').read_text())
    rows=[]
    masks=[]
    for b in benchmarks:
        e=b['rule'];m=evaluate(c,e);co=(m@W)[None,:];sc=(m@S)[None,:]
        r={'gate_name':b['name'],'exact_retention_rule':label(e),'rule_json':json.dumps(e,separators=(',',':'))}
        r.update({k:float(v[0]) for k,v in derived(co,sc,totals).items()})
        for j,name in enumerate(COUNT_NAMES):
            r['retained_'+name+'_n']=int(co[0,j]);r['baseline_'+name+'_n']=int(totals[j]);r['baseline_'+name+'_fraction']=totals[j]/totals[0]
        r['baseline_supported_fraction']=(totals[1]+totals[5])/totals[0]
        r['baseline_pv_all_dual_fraction']=(totals[1]+totals[2])/totals[0]
        r['baseline_sst_score_mean']=S[:,1].mean()
        r['pv_enrichment_class']='WEAK' if r['pv_fold_enrichment']<1.5 else '>=1.5x'
        r['supported_enrichment_class']='WEAK' if r['supported_fold_enrichment']<1.5 else '>=1.5x'
        rows.append(r);masks.append(m)
    pd.DataFrame(rows).to_csv(OUT/'tables/frozen_benchmark_metrics.tsv',sep='\t',index=False)
    np.savez_compressed(OUT/'tables/frozen_benchmark_masks.npz',cell_id=c.cell_id.to_numpy(str),masks=np.asarray(masks),gate_name=np.asarray([b['name'] for b in benchmarks]))
    print(json.dumps({'direct_masks_checked':len(sample_rows),'benchmark_rows':len(rows),'all_passed':True}),flush=True)


if __name__=='__main__':main()
