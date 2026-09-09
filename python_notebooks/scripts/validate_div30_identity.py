#!/usr/bin/env python3
"""Independent source/order check for the frozen target and archived neighbors."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
ROOT=Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results')
OLD=ROOT/'div30_pv_antecedent_v1'
OUT=ROOT/'div30_antecedent_identity_v1'

def read(p):return pd.read_csv(p,sep='\t',float_precision='round_trip')

def main():
    old=read(OLD/'cells.tsv.gz');new=read(OUT/'within/cells.tsv.gz')
    assert len(old)==len(new)==90631 and old.cell_id.equals(new.cell_id)
    for col in old:
        assert col in new,col
        if pd.api.types.is_numeric_dtype(old[col]):
            assert np.allclose(old[col],new[col],rtol=0,atol=1e-12,equal_nan=True),col
        else:
            assert old[col].fillna('').equals(new[col].fillna('')),col
    # Reapply the frozen rule from the old selected archive and old expression.
    selected=read(OLD/'tables/surface_gate_best_at_recovery_floors.tsv')
    chosen=selected[(selected.target=='top20')&(selected.recovery_floor==.5)&(selected.maximum_markers==3)].iloc[0]
    rules=json.loads(chosen.rules);genes=json.loads((OLD/'cache/genes.json').read_text())
    expr=np.load(OLD/'cache/expression.npy',mmap_mode='r')
    clauses={}
    for r in rules:
        values=np.asarray(expr[:,genes.index(r['gene'])],dtype='float64')
        cut=np.float64(r['threshold'])
        clauses[r['gene']]=values>cut if r['op']=='>' else values<=cut
    keep=np.logical_and.reduce(list(clauses.values()))
    target=old.early_pv_score>=old.early_pv_score.quantile(.8)
    assert np.array_equal(keep,new.retained) and np.array_equal(target,new.target)
    groups=np.select([target&keep,target&~keep,~target&keep],['A','B','C'],default='D')
    assert np.array_equal(groups,new.group)
    assert [int((groups==g).sum()) for g in 'ABCD']==[9210,8917,6571,65933]
    ab=read(OUT/'within/tables/gate_ablations.tsv')
    for r in ab.itertuples():
        sub=np.logical_and.reduce([clauses[g] for g in r.markers.split('+')]) if r.markers!='none' else np.ones(len(old),bool)
        assert sub.sum()==r.retained_n and (sub&target).sum()==r.target_retained_n
    branches=read(OUT/'within/tables/gate_decision_branches.tsv')
    truth=read(OUT/'within/tables/gate_truth_patterns.tsv')
    assert branches.n.sum()==truth.n.sum()==len(old)
    # Recompute every archived empirical state fraction, independently of the
    # new PCA reconstruction, new summary tables and reverse-map implementation.
    neighbors=np.load(OLD/'cache/reference_neighbors.npz')['reference_indices']
    ref=read(OLD/'tables/reference_cell_index.tsv')
    forward=read(OUT/'mapping/forward_cell_neighbors.tsv.gz')
    assert forward.cell_id.equals(old.cell_id)
    state_names={'strict_PV':'PV-biased','dual_high':'PV/SST hybrid','SST':'SST-biased','unresolved':'unresolved/immature'}
    assert set(ref.state)==set(state_names.values()),set(ref.state)
    totals=np.zeros(len(old));means={}
    for name,label in state_names.items():
        frac=(ref.state.to_numpy()[neighbors]==label).mean(axis=1)
        assert np.allclose(frac,forward['archived_unbalanced__'+name],rtol=0,atol=1e-15)
        totals+=frac;means[name]=dict(retained=float(frac[keep].mean()),removed=float(frac[~keep].mean()))
    assert np.allclose(totals,1)
    result=dict(passed=True,n_cells=len(old),all_original_cell_columns_preserved=list(old.columns),
        exact_primary_retained_n=int(keep.sum()),exact_target_n=int(target.sum()),
        group_counts={g:int((groups==g).sum()) for g in 'ABCD'},
        fixed_ablation_subsets_verified=len(ab),full_forward_cell_state_fractions_verified=len(old)*4,
        forward_means=means,large_gate_search_rerun=False,target_changed=False)
    (OUT/'provenance/independent_identity_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='all_original_cell_columns_preserved'},indent=2))

if __name__=='__main__':main()
