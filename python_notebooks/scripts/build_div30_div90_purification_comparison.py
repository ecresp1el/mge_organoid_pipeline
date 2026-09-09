#!/usr/bin/env python3
"""Preserve frozen inputs and assemble independently solved DIV30/DIV90 work.

Scientific artifacts and caches belong on Turbo; this checkout holds code only.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

from div90_hypergate_paths import PROJECT_ROOT, REPO_ROOT

OUT = PROJECT_ROOT / 'results/div30_div90_purification_comparison_v1'
DIV90 = PROJECT_ROOT / 'results/div90_purification_audit_v1'
DIV30 = PROJECT_ROOT / 'results/div30_pv_antecedent_v1'
FROZEN = [PROJECT_ROOT / 'results/div90_hypergate_sst_pv',
          PROJECT_ROOT / 'results/div90_hypergate_sst_pv_phase2']
REQUESTS = {
    'DIV90_user_instructions.txt': Path('/home/elcrespo/.codex/attachments/9babca15-b519-4690-b092-e70d6ee5eb6e/pasted-text.txt'),
    'DIV30_user_instructions.txt': Path('/home/elcrespo/.codex/attachments/9c2190da-6499-47a4-97d0-df85d8826e5e/pasted-text.txt'),
}


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')


def inventory():
    result = {}
    for root in FROZEN:
        for p in sorted(root.rglob('*')):
            if p.is_file():
                rel = p.relative_to(root)
                # Runtime dependency and library cache bytes are not scientific
                # results. Record their metadata; hash all analysis/PI artifacts.
                auxiliary = rel.parts[0] in {'cache', 'runtime_deps'}
                info = {'size': p.stat().st_size, 'mtime_ns': p.stat().st_mtime_ns}
                if not auxiliary:
                    info['sha256'] = sha(p)
                result[str(p)] = info
    return result


def initialize():
    (OUT / 'provenance').mkdir(parents=True, exist_ok=True)
    for name, source in REQUESTS.items():
        dest = OUT / 'provenance' / name
        if not dest.exists():
            shutil.copy2(source, dest)
        assert sha(dest) == sha(source)
    dest = OUT / 'provenance/frozen_inputs_before.json'
    if dest.exists():
        print('Preserving original before-manifest', flush=True)
    else:
        write_json(dest, inventory())
        print('Frozen phase-one, phase-two and PI inventory recorded', flush=True)
    write_json(OUT / 'analysis_contract.json', {
        'div90': str(DIV90), 'div30': str(DIV30),
        'source_package': str(FROZEN[1] / 'pi_figure_package_v1'),
        'starting_conclusion': 'FGFR2/PTPRS preserves about 89% of PV-biased and dual-high cells; SST-biased composition changes 22.5% to 21.3%. Weak purification; independent statistical significance is not established.',
        'separate_optimization': True, 'rerun_full_biological_pipeline': False,
        'principal_objective': 'Actual retained purity at explicit recovery thresholds; report enrichment, contamination and cell yield.',
        'scope_limit': 'Observed RNA gate search; three-marker searches are bounded unless explicitly proven exhaustive. No claim of a universal theoretical surfaceome maximum.',
        'no_antibody_recommendations': True,
        'cross_sectional_interpretation': 'Transcriptional antecedent and developmental resemblance, never lineage tracing.',
    })


def verify():
    before = json.loads((OUT / 'provenance/frozen_inputs_before.json').read_text())
    after = inventory()
    changed = [p for p in before if after.get(p) != before[p]]
    added = sorted(set(after) - set(before))
    result = dict(all_existing_results_preserved=not changed,
                  files_checked=len(before), scientific_files_sha256_checked=sum('sha256' in v for v in before.values()),
                  changed_or_missing=changed, added_within_frozen_directories=added)
    write_json(OUT / 'provenance/source_preservation.json', result)
    print(json.dumps(result), flush=True)
    assert not changed and not added, 'Frozen result directories changed'


def legacy_check():
    """Independent count-derived optima, without importing the discovery engine."""
    import h5py
    import numpy as np
    import pandas as pd

    source = FROZEN[1]
    cells = pd.read_csv(source / 'cells.tsv.gz', sep='\t')
    totals = cells.state.value_counts()
    npv = int(totals['PV-biased'])
    rows, checked = [], 0
    for fn in ['depletion_all_candidate_metrics.h5', 'depletion_supplemental_candidate_metrics.h5']:
        with h5py.File(source / 'tables' / fn, 'r') as h:
            for name, group in h.items():
                a = group['metrics'][:]
                columns = json.loads(group.attrs['metric_columns'])
                get = lambda col: a[:, columns.index(col)]
                n, pv, sst, dual = [get(k) for k in ['retained_n', 'retained_pv_n', 'retained_sst_n', 'retained_hybrid_n']]
                assert np.array_equal(n, pv + sst + dual + get('retained_unresolved_n'))
                assert np.array_equal(n + get('removed_n'), np.full(len(n), len(cells)))
                checked += len(n)
                purity = np.divide(pv, n, out=np.zeros_like(pv), where=n > 0)
                for floor in [0, .1, .2, .3, .4, .5, .6, .7, .8]:
                    ids = np.flatnonzero((pv >= npv * floor) & (n > 0))
                    if not len(ids):
                        continue
                    i = int(ids[np.lexsort((-n[ids], -pv[ids], -purity[ids]))[0]])
                    if 'rules_json' in group:
                        rule = json.loads(group['rules_json'][i])
                    elif name == 'exact_single':
                        genes = json.loads(h.attrs['genes'])
                        rule = dict(action='remove', logic='AND', rules=[dict(gene=genes[int(group['gene_index'][i])], op='>', threshold=float(group['threshold'][i]))])
                    else:
                        defs = json.loads(h.attrs['pair_threshold_definitions'])
                        rule = dict(action='remove', logic=name.rsplit('_', 1)[1].upper(), rules=[defs[int(group[k][i])] for k in ['threshold_a_index', 'threshold_b_index']])
                    rows.append(dict(file=fn, group=name, row=i, min_recovery=floor,
                                     retained_n=int(n[i]), pv_n=int(pv[i]), sst_n=int(sst[i]), dual_n=int(dual[i]),
                                     pv_purity=float(purity[i]), pv_recovery=float(pv[i] / npv), sst_fraction=float(sst[i] / n[i]),
                                     rule_json=json.dumps(rule, separators=(',', ':'))))
    d = pd.DataFrame(rows).sort_values(['min_recovery', 'pv_purity', 'pv_recovery'], ascending=[True, False, False]).groupby('min_recovery').head(1)
    for r in d.itertuples():
        rule = json.loads(r.rule_json)
        masks = []
        for term in rule['rules']:
            x = cells[term['gene']].to_numpy()
            op = {'>': np.greater, '>=': np.greater_equal, '<': np.less, '<=': np.less_equal}[term['op']]
            masks.append(op(x, term['threshold']))
        selected = np.logical_or.reduce(masks) if rule['logic'] == 'OR' else np.logical_and.reduce(masks)
        retained = ~selected if rule['action'] == 'remove' else selected
        assert retained.sum() == r.retained_n
        assert (retained & cells.state.eq('PV-biased')).sum() == r.pv_n
        assert (retained & cells.state.eq('SST-biased')).sum() == r.sst_n
    d.to_csv(OUT / 'provenance/independent_legacy_purity_optima.tsv', sep='\t', index=False)
    write_json(OUT / 'provenance/independent_legacy_validation.json', dict(rows_checked=checked, count_conservation=True, direct_masks_verified=len(d), full_biological_rerun=False))
    assert checked == 2885531
    print(d.drop(columns='rule_json').to_string(index=False), flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--stage', choices=['initialize', 'verify', 'legacy-check'], required=True)
    args = ap.parse_args()
    {'initialize': initialize, 'verify': verify, 'legacy-check': legacy_check}[args.stage]()
