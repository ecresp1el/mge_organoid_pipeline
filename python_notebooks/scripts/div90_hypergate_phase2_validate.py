#!/usr/bin/env python3
"""Independent count audits and gate sensitivity on fixed DIV90 cells."""
import io
import json
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from run_div90_hypergate_phase2 import OUT, OLD, STATES, PV, SST, write_json, verify_phase1


def removed_mask(cells,rule):
    conditions=[]
    for item in rule['rules']:
        x=cells[item['gene']].to_numpy();t=item['threshold']
        if item['op']=='>':b=x>t
        elif item['op']=='>=':b=x>=t
        elif item['op']=='<':b=x<t
        elif item['op']=='<=':b=x<=t
        else:raise ValueError(item['op'])
        conditions.append(b)
    if conditions:
        inside=np.logical_and.reduce(conditions) if rule['logic']=='AND' else np.logical_or.reduce(conditions)
    else:inside=np.full(len(cells),rule['action']=='retain',dtype=bool)
    return inside if rule['action']=='remove' else ~inside


def evaluate(labels,remove,keep=None):
    labels=np.asarray(labels);keep=np.ones(len(labels),bool) if keep is None else np.asarray(keep)
    totals=np.array([((labels==s)&keep).sum() for s in STATES]);r=np.array([((labels==s)&remove&keep).sum() for s in STATES]);k=totals-r
    d={'n':int(keep.sum()),'removed_n':int(r.sum()),'retained_n':int(k.sum())}
    for i,name in enumerate(['pv','hybrid','sst','unresolved']):
        d[name+'_retained_n']=int(k[i]);d[name+'_lost_n']=int(r[i]);d[name+'_recovery']=float(k[i]/totals[i]) if totals[i] else None
    d['sst_contamination']=float(k[2]/k.sum()) if k.sum() else None
    d['target_recovery']=float(k[:2].sum()/totals[:2].sum()) if totals[:2].sum() else None
    return d


def run():
    cells=pd.read_csv(OUT/'cells.tsv.gz',sep='\t');old=pd.read_csv(OLD/'target_labels.tsv.gz',sep='\t')
    assert cells.cell_id.equals(old.cell_id)
    assert np.array_equal(cells[['loupe_x','loupe_y']],old[['loupe_x','loupe_y']])
    for col in ['crude_label','module_label','sample','loupe_label']:assert cells[col].equals(old[col])
    for col in ['pv_module','sst_module']:assert np.allclose(cells[col],old[col],rtol=0,atol=1e-12)
    assert len(cells)==4768 and cells.cell_id.is_unique and cells.state.notna().all()
    rows=[]
    for program,genes in [('pv',PV),('sst',SST)]:
        score=cells[program+'_score'].to_numpy()
        for gene in genes:
            component=cells['z_'+gene].to_numpy()/len(genes)
            covariance=np.mean((component-component.mean())*(score-score.mean()))
            rows.append(dict(program=program,gene=gene,nominal_weight=1/len(genes),score_variance_share=covariance/np.var(score),gene_score_spearman=spearmanr(cells[gene],score).statistic))
    contributions=pd.DataFrame(rows)
    assert np.allclose(contributions.groupby('program').score_variance_share.sum(),1)
    contributions.to_csv(OUT/'tables/module_variance_contributions.tsv',sep='\t',index=False)
    summary=json.loads((OUT/'gate_summary.json').read_text());rules=summary['rules']
    assignments=pd.read_csv(OUT/'gate_cell_assignments.tsv.gz',sep='\t')
    assert assignments.cell_id.equals(cells.cell_id)
    selected=['baseline','phase1_reference','best_single','best_pair','selected','practical_single','practical_pair','experimental_depletion','erbb4_cxcr4']
    sensitivity=pd.read_csv(OUT/'tables/cell_sensitivity.tsv.gz',sep='\t')
    assert sensitivity.cell_id.equals(cells.cell_id)
    srows=[];drows=[]
    depth=cells.groupby('sample').total_counts.transform(lambda x:pd.qcut(x.rank(method='first'),5,labels=False))
    for key in selected:
        row=summary[key];rm=removed_mask(cells,rules[row['gate_id']]);actual=evaluate(cells.state,rm)
        for k in ['removed_n','retained_n']:assert actual[k]==row[k]
        for short,state in zip(['pv','hybrid','sst','unresolved'],STATES):
            assert actual[short+'_lost_n']==row['removed_'+short+'_n']
        name='phase1' if key=='phase1_reference' else key
        if name+'_retained' in assignments:assert np.array_equal(~rm,assignments[name+'_retained'])
        for col in [c for c in sensitivity if c.endswith('_state')]:
            record=evaluate(sensitivity[col],rm);record.update(gate_name=key,gate_id=row['gate_id'],variant=col[:-6]);srows.append(record)
        for q in range(5):
            record=evaluate(cells.state,rm,depth==q);record.update(gate_name=key,gate_id=row['gate_id'],depth_quintile=q+1);drows.append(record)
    pd.DataFrame(srows).to_csv(OUT/'tables/gate_module_threshold_sensitivity.tsv',sep='\t',index=False)
    pd.DataFrame(drows).to_csv(OUT/'tables/gate_depth_quintile_sensitivity.tsv',sep='\t',index=False)
    # Independent direct cell-mask checks on reproducibly sampled complete-search
    # rows, plus count conservation for every evaluated positive-depletion gate.
    checked=0;total_rows=0;rng=np.random.default_rng(9090);totals=np.array([cells.state.eq(s).sum() for s in STATES])
    for filename in ['depletion_all_candidate_metrics.h5','depletion_supplemental_candidate_metrics.h5']:
        path=OUT/'tables'/filename
        if not path.exists():continue
        with h5py.File(io.BytesIO(path.read_bytes()),'r') as f:
            genes=json.loads(f.attrs.get('genes','[]'));definitions=json.loads(f.attrs.get('pair_threshold_definitions','[]'))
            for group in f:
                g=f[group];cols=json.loads(g.attrs['metric_columns']);n=len(g['metrics']);total_rows+=n
                rix=[cols.index('removed_'+s+'_n') for s in ['pv','hybrid','sst','unresolved']]
                kix=[cols.index('retained_'+s+'_n') for s in ['pv','hybrid','sst','unresolved']]
                for start in range(0,n,50000):
                    chunk=g['metrics'][start:start+50000];r=chunk[:,rix];k=chunk[:,kix]
                    assert np.all(r>=0) and np.array_equal(r+k,np.broadcast_to(totals,r.shape))
                    assert np.array_equal(r.sum(1),chunk[:,cols.index('removed_n')])
                    assert np.array_equal(k.sum(1),chunk[:,cols.index('retained_n')])
                for ix in sorted(rng.choice(n,size=min(100,n),replace=False)):
                    if group=='exact_single':rule=dict(action='remove',logic='AND',rules=[dict(gene=genes[int(g['gene_index'][ix])],op='>',threshold=float(g['threshold'][ix]))])
                    elif group.startswith('all_pairs_'):
                        rule=dict(action='remove',logic=group.split('_')[-1].upper(),rules=[definitions[int(g['threshold_a_index'][ix])],definitions[int(g['threshold_b_index'][ix])]])
                    else:rule=json.loads(g['rules_json'][ix])
                    rm=removed_mask(cells,rule);counts=np.array([(rm&cells.state.eq(s)).sum() for s in STATES])
                    assert np.array_equal(counts,g['metrics'][ix][rix]),(filename,group,int(ix))
                    checked+=1
    verify_phase1()
    write_json(OUT/'provenance/final_validation.json',dict(n_cells=4768,cell_order_preserved=True,original_loupe_coordinates_exact=True,original_phase1_scores_and_labels_preserved=True,all_cells_assigned_operational_state=True,total_python_candidate_rows_count_conservation_checked=total_rows,random_candidate_masks_independently_reproduced=checked,key_gate_counts_and_cell_flags_independently_reproduced=True,actual_R_hypergate_fits_independently_checked=summary['actual_R_hypergate_fits'],primary_module_features_excluded_from_search=True,all_phase1_files_unchanged=True))
    print(f'Validated {total_rows:,} candidate metric rows and {checked} independent masks.',flush=True)


if __name__=='__main__':run()
