#!/usr/bin/env python3
"""Independently verify saved DIV30 gate counts against direct cell masks.

This audit imports no optimization or metric code from the analysis helper.
It only reads results and writes a new provenance JSON after all checks pass.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True

import h5py
import numpy as np
import pandas as pd


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',type=Path,default=Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_pv_antecedent_v1'))
    args=parser.parse_args();out=args.out
    result_path=out/'provenance/surface_independent_validation.json'
    if result_path.exists(): raise FileExistsError(f'Preserving existing audit: {result_path}')
    summary_path=out/'provenance/surface_validation.json'
    if summary_path.exists(): raise FileExistsError(f'Preserving existing audit: {summary_path}')
    cells=pd.read_csv(out/'cells.tsv.gz',sep='\t')
    expression=np.load(out/'cache/expression.npy',mmap_mode='r')
    genes=json.loads((out/'cache/genes.json').read_text());lookup={g:i for i,g in enumerate(genes)}
    assert cells.cell_id.is_unique and len(cells)==expression.shape[0]
    assert len(genes)==expression.shape[1] and len(genes)==len(set(genes))
    score=cells.early_pv_score.to_numpy(float)
    targets={f'top{q}':score>=np.quantile(score,1-q/100) for q in (10,20,30)}
    targets['top20_strong_MGE']=targets['top20']&(cells.mge_score.to_numpy(float)>=cells.mge_score.median())
    targets['top20_postmitotic']=targets['top20']&cells.postmitotic.to_numpy(bool)
    support=json.loads((out/'provenance/sst_support.json').read_text())
    if support['stable_for_low_sst_target']:
        targets['top20_lowSST']=targets['top20']&(cells.direct_sst_score.to_numpy(float)<cells.direct_sst_score.median())
    design=json.loads((out/'provenance/DIV30_surface_search_design.json').read_text())
    eligible=set(design['eligible_genes']);excluded=set(design['excluded_score_mapping_supervised_genes'])
    assert not eligible&excluded
    best=pd.read_csv(out/'tables/surface_gate_best_at_recovery_floors.tsv',sep='\t')

    clause_masks={}
    def mask(clauses,logic):
        values=[]
        for c in clauses:
            assert c['op'] in ('<=','>')
            key=(c['gene'],c['op'],float(c['threshold']))
            if key not in clause_masks:
                x=expression[:,lookup[c['gene']]]
                clause_masks[key]=x<=np.float64(c['threshold']) if c['op']=='<=' else x>np.float64(c['threshold'])
            values.append(clause_masks[key])
        if not values: return np.ones(len(cells),bool)
        assert logic in ('AND','OR')
        return np.logical_and.reduce(values) if logic=='AND' else np.logical_or.reduce(values)

    checked=[]
    for row in best.to_dict('records'):
        y=targets[row['target']];clauses=json.loads(row['rules']);keep=mask(clauses,row['logic'])
        assert {c['gene'] for c in clauses}<=eligible
        n,tp=int(keep.sum()),int((keep&y).sum());N,T=len(y),int(y.sum())
        assert n>0 and T>0
        expected=dict(retained_n=n,target_retained_n=tp,starting_target_fraction=T/N,
                      post_gate_target_fraction=tp/n,fold_enrichment=(tp/n)/(T/N),
                      target_recovery=tp/T,total_cell_yield=n/N,composition_change_pp=100*(tp/n-T/N))
        for key,value in expected.items():
            assert np.isclose(row[key],value,rtol=1e-11,atol=1e-11),(row['target'],row['gate_id'],key,row[key],value)
        assert expected['target_recovery']>=row['recovery_floor']-1e-12
        assert len({c['gene'] for c in clauses})==row['n_markers']<=row['maximum_markers']
        checked.append(dict(target=row['target'],gate_id=row['gate_id'],recovery_floor=row['recovery_floor'],maximum_markers=row['maximum_markers'],retained_n=n,target_retained_n=tp))

    rng=np.random.default_rng(20260909);single_checks=pair_checks=triplet_checks=0
    bounded_single_rules=bounded_pair_rules=0;frontiers={}
    with h5py.File(out/'tables/surface_gate_complete_counts.h5','r') as h5:
        for target in h5:
            group=h5[target];y=targets[target]
            assert int(group.attrs['n_total'])==len(y) and int(group.attrs['n_target'])==int(y.sum())
            marker_names=list(group['exact_single'])
            for gene in marker_names:
                sg=group['exact_single'][gene];thresholds=sg['threshold'][:]
                n,tp=sg['retained_n_le'][:],sg['target_n_le'][:]
                assert np.all(np.isfinite(thresholds)) and np.all(np.diff(thresholds)>0)
                assert np.all(np.diff(n)>0) and np.all(np.diff(tp)>=0)
                assert np.all(np.diff(tp)<=np.diff(n))
                assert np.all((0<=tp)&(tp<=n)&(n<=len(y))&(tp<=int(y.sum())))
                assert np.all((int(y.sum())-tp)<=(len(y)-n))
                assert n[-1]==len(y) and tp[-1]==int(y.sum())
                bounded_single_rules+=2*len(n)
            for _ in range(64):
                gene=marker_names[int(rng.integers(len(marker_names)))];g=group['exact_single'][gene]
                j=int(rng.integers(len(g['threshold'])));threshold=float(g['threshold'][j])
                n,tp=int(g['retained_n_le'][j]),int(g['target_n_le'][j])
                for op,expected_n,expected_tp in [('<=',n,tp),('>',len(y)-n,int(y.sum())-tp)]:
                    keep=mask([dict(gene=gene,op=op,threshold=threshold)],'AND')
                    assert int(keep.sum())==expected_n and int((keep&y).sum())==expected_tp
                    single_checks+=1
            pair_names=list(group['screened_pairs'])
            for pair in pair_names:
                pg=group['screened_pairs'][pair]
                for rule in pg:
                    if rule.startswith('threshold_'): continue
                    counts=pg[rule][:];n,tp=counts[:,0],counts[:,1]
                    assert len(counts)==len(pg['threshold_a'])*len(pg['threshold_b'])
                    assert np.all((0<=tp)&(tp<=n)&(n<=len(y))&(tp<=int(y.sum())))
                    assert np.all((int(y.sum())-tp)<=(len(y)-n))
                    bounded_pair_rules+=len(counts)
            for _ in range(64):
                name=pair_names[int(rng.integers(len(pair_names)))];g=group['screened_pairs'][name]
                rules=[k for k in g if k not in ('threshold_a','threshold_b')]
                rule=rules[int(rng.integers(len(rules)))];a,b,logic=rule.split('_')
                ia=int(rng.integers(len(g['threshold_a'])));ib=int(rng.integers(len(g['threshold_b'])))
                clauses=[dict(gene=g.attrs['gene_a'],op='<=' if a=='le' else '>',threshold=float(g['threshold_a'][ia])),
                         dict(gene=g.attrs['gene_b'],op='<=' if b=='le' else '>',threshold=float(g['threshold_b'][ib]))]
                keep=mask(clauses,logic);j=ia*len(g['threshold_b'])+ib
                expected=g[rule][j]
                assert np.array_equal(expected,[int(keep.sum()),int((keep&y).sum())]),(target,name,rule,ia,ib)
                pair_checks+=1
            frontier=group['best_single_or_pair_purity_at_minimum_target_count'][:]
            frontiers[target]=frontier
            assert np.all(np.diff(frontier)<=1e-12)
            for row in best[(best.target==target)&(best.n_markers==3)].to_dict('records'):
                gain=row['post_gate_target_fraction']-frontier[int(row['target_retained_n'])]
                assert gain>=.02-1e-12,(target,row['gate_id'],gain)
                triplet_checks+=1
            print(f'Complete archive bounds/conservation and sampled direct masks passed: {target}',flush=True)
    triples=pd.read_csv(out/'tables/surface_triplet_materiality_audit.tsv.gz',sep='\t')
    for row in triples.to_dict('records'):
        y=targets[row['target']];keep=mask(json.loads(row['rules']),row['logic'])
        n,tp=int(keep.sum()),int((keep&y).sum())
        assert n==row['retained_n'] and tp==row['target_retained_n'],('triplet',row['target'],row['gate_id'],n,tp,row['retained_n'],row['target_retained_n'])
        gain=tp/n-frontiers[row['target']][tp]
        assert np.isclose(gain,row['purity_gain_vs_1_or_2_markers_at_equal_or_higher_recovery'],rtol=1e-11,atol=1e-11)
        assert bool(row['material_triplet_improvement'])==(gain>=.02)
    print(f'All {len(triples)} triplet direct masks and materiality flags passed',flush=True)
    for (_,floor),df in best.groupby(['target','recovery_floor']):
        assert np.all(np.diff(df.sort_values('maximum_markers').post_gate_target_fraction)>=-1e-12)
    benchmark=pd.read_csv(out/'tables/frozen_DIV90_gate_on_DIV30.tsv',sep='\t')
    frozen=json.loads((out/'provenance/frozen_DIV90_benchmark_rule.json').read_text())
    assert frozen['source_key']=='experimental_depletion'
    original=json.loads(Path(frozen['source']).read_text())['experimental_depletion']
    assert frozen['source_gate_id']==original['gate_id']
    original_clauses=json.loads(original['rules'])
    assert frozen['original_rule']==original_clauses
    for row in benchmark.to_dict('records'):
        y=targets[row['target']];keep=mask(json.loads(row['rules']),row['logic'])
        assert int(keep.sum())==row['retained_n'] and int((keep&y).sum())==row['target_retained_n']
    validation=pd.read_csv(out/'tables/held_sample_surface_gate_validation.tsv',sep='\t')
    for row in validation.to_dict('records'):
        test=cells['sample'].astype(str).eq(str(row['held_out_sample'])).to_numpy()
        keep=mask(json.loads(row['rules']),row['logic'])&test;y=targets[row['target']]
        n,tp=int(keep.sum()),int((keep&y).sum());T=int(y[test].sum());N=int(test.sum())
        assert n==row['retained_n'] and tp==row['target_retained_n'],('held_sample',row['held_out_sample'],row['gate_id'],n,tp,row['retained_n'],row['target_retained_n'])
        assert np.isclose(row['post_gate_target_fraction'],tp/n if n else np.nan,atol=1e-11,rtol=1e-11,equal_nan=True)
        assert np.isclose(row['target_recovery'],tp/T,atol=1e-11,rtol=1e-11)
        assert np.isclose(row['total_cell_yield'],n/N,atol=1e-11,rtol=1e-11)
    result=dict(passed=True,n_DIV30_cells=len(cells),best_gate_rows_checked=len(checked),best_gate_direct_counts=checked,
                random_exact_single_direction_checks=single_checks,random_signed_pair_checks=pair_checks,
                material_best_triplet_checks=triplet_checks,frozen_DIV90_benchmark_rows_checked=len(benchmark),
                held_sample_direct_mask_rows_checked=len(validation),
                all_exact_single_rules_bounds_checked=bounded_single_rules,all_screened_pair_rules_bounds_checked=bounded_pair_rules,
                all_triplet_direct_mask_rows_checked=len(triples),
                random_seed=20260909,score_mapping_gene_overlap=len(eligible&excluded),
                checks=['all best-gate counts and derived composition/enrichment/recovery/yield match independent direct masks',
                        'random complete HDF5 single cut points and both complements match direct masks',
                        'random complete HDF5 signed AND/OR pair counts match direct masks',
                        'triplet improvements checked against full evaluated single/pair frontier',
                        'fixed recovery floors and monotonicity by maximum marker count verified',
                        'frozen benchmark matches experimental_depletion provenance and exact thresholds'],
                scope='Numerical implementation audit only; does not establish biological significance or external prospective validity.',
                validator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    result_path.write_text(json.dumps(result,indent=2)+'\n')
    summary_path.write_text(json.dumps({k:v for k,v in result.items() if k!='best_gate_direct_counts'},indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='best_gate_direct_counts'},indent=2))


if __name__=='__main__': main()
