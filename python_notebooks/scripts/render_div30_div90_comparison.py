#!/usr/bin/env python3
"""Link completed, independently optimized ages without fitting a pooled gate."""
from __future__ import annotations

import json
import os
from pathlib import Path
import textwrap

ROOT = Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results')
OUT = ROOT / 'div30_div90_purification_comparison_v1'
D30 = ROOT / 'div30_pv_antecedent_v1'
D90 = ROOT / 'div90_purification_audit_v1'
OLD90 = ROOT / 'div90_hypergate_sst_pv_phase2'
os.environ['MPLCONFIGDIR'] = str(OUT / 'cache/matplotlib')

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

from render_div90_purification_report import mask_from_json


def read(base, name):
    return pd.read_csv(base / 'tables' / name, sep='\t', float_precision='round_trip')


def markdown_table(frame):
    def cell(value):
        if pd.isna(value): return ''
        if isinstance(value, (float, np.floating)): return f'{value:.3f}'
        return str(value).replace('|', '\\|').replace('\n', '<br>')
    lines = ['| ' + ' | '.join(map(str, frame.columns)) + ' |',
             '| ' + ' | '.join(['---'] * len(frame.columns)) + ' |']
    lines.extend('| ' + ' | '.join(cell(v) for v in row) + ' |' for row in frame.itertuples(index=False, name=None))
    return '\n'.join(lines)


def get_genes(e):
    if 'gene' in e: return {e['gene']}
    if 'not' in e: return get_genes(e['not'])
    return set().union(*(get_genes(x) for x in e['terms']))


def add_div90_expression(cells, genes):
    missing = sorted(set(genes) - set(cells.columns))
    if not missing: return cells
    with h5py.File(ROOT / 'python_anndata/varela_div90.h5ad', 'r') as f:
        names = f['var'][f['var'].attrs['_index']].asstr()[:]
        cols = pd.Index(names).get_indexer(missing)
        if (cols < 0).any(): raise ValueError('Cross-age gate gene missing from DIV90')
        ids = f['obs/cell_id'].asstr()[:]
        positions = pd.Index(ids).get_indexer(cells.cell_id)
        assert (positions >= 0).all()
        inverse = np.full(len(ids), -1, int); inverse[positions] = np.arange(len(cells))
        x = f['X']; ip = x['indptr'][:]
        vals = np.zeros((len(cells), len(missing)))
        for start in range(0, len(ids), 1200):
            stop = min(len(ids), start + 1200)
            selected = inverse[start:stop] >= 0
            if not selected.any(): continue
            a, b = int(ip[start]), int(ip[stop])
            mat = sparse.csr_matrix((x['data'][a:b], x['indices'][a:b], ip[start:stop+1] - a), shape=(stop-start, len(names)))
            vals[inverse[start:stop][selected]] = mat[selected][:, cols].toarray()
        for j, g in enumerate(missing): cells[g] = vals[:, j]
    return cells


def cross_age(c30, c90, best30, best90):
    e30 = {'logic': best30['logic'], 'terms': json.loads(best30['rules'])}
    e90 = json.loads(best90['rule_json'])
    genes = get_genes(e30) | get_genes(e90)
    expression = np.load(D30 / 'cache/expression.npy', mmap_mode='r')
    gene_index = {g: i for i, g in enumerate(json.loads((D30 / 'cache/genes.json').read_text()))}
    d30 = c30.copy()
    for g in genes:
        if g not in gene_index: raise ValueError(f'DIV30 cross-age expression unavailable: {g}')
        d30[g] = expression[:, gene_index[g]]
    c90 = add_div90_expression(c90.copy(), genes)
    targets = {'DIV30': (d30, d30.early_pv_score >= d30.early_pv_score.quantile(.8), 'DIV30 top20 early-PV rank consensus'),
               'DIV90': (c90, c90.state.eq('PV-biased'), 'DIV90 strict PV-biased')}
    rows = []
    for rule_age, e in [('DIV30', e30), ('DIV90', e90)]:
        for data_age, (c, y, description) in targets.items():
            m = mask_from_json(c, e); n = int(m.sum()); tp = int((m & y).sum())
            purity = tp / n if n else np.nan
            rows.append(dict(rule_selected_at=rule_age, evaluated_at=data_age, target=description,
                             n_total=len(c), n_target=int(y.sum()), retained_n=n, retained_target_n=tp,
                             starting_target_fraction=float(y.mean()), post_gate_target_fraction=purity,
                             fold_enrichment=purity / y.mean(), target_recovery=tp / y.sum(), total_yield=n / len(c),
                             rule_json=json.dumps(e, separators=(',', ':')),
                             interpretation='Exact RNA-threshold diagnostic only; independent age-specific targets and denominators, no pooled fitting or protein threshold equivalence.'))
    table = pd.DataFrame(rows)
    own30 = table[(table.rule_selected_at == 'DIV30') & (table.evaluated_at == 'DIV30')].iloc[0]
    own90 = table[(table.rule_selected_at == 'DIV90') & (table.evaluated_at == 'DIV90')].iloc[0]
    assert own30.retained_n == best30['retained_n']
    assert own30.retained_target_n == best30['target_retained_n']
    assert own90.retained_n == best90['retained_n']
    assert own90.retained_target_n == best90['retained_pv_n']
    table.to_csv(OUT / 'tables/cross_age_gate_benchmarks.tsv', sep='\t', index=False)
    return table, genes


def schematic(s30, b30, b90, loss):
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'svg.fonttype': 'none', 'pdf.fonttype': 42})
    fig, ax = plt.subplots(figsize=(13.5, 6))
    ax.set(xlim=(0, 1), ylim=(0, 1)); ax.axis('off')
    xs = [.02, .37, .72]; width = .26
    colors = ['#ECF4F2', '#EAF1F7', '#F7EFE4']
    headers = ['DIV30', 'DIV90', '3–4+ months']
    bodies = [
        f"MGE / precursor context\nCandidates mainly postmitotic\n\nGRIA2-high AND OPCML-high\nAND NOTCH1-low\n\n{b30['fold_enrichment']:.2f}× top20 enrichment\nat {b30['target_recovery']:.1%} recovery\n\nERBB4 RNA entry would lose\n{loss:.1%} of top20 candidates",
        f"Frozen cortical\nLHX6+/ERBB4+ entry\n\n(ERBB4-high AND LRRC4B-low)\nOR ACKR3-high\n\n{b90['pv_fraction']:.1%} strict PV purity\n({b90['pv_fold_enrichment']:.2f}× baseline)\nat {b90['pv_recovery']:.1%} recovery\n\nProtein gate unvalidated",
        "Later PVALB protein\n\nSUPPORTED\nUser-reported system endpoint\n\nLink from either profiled\ncell population to this endpoint:\nTO BE VALIDATED",
    ]
    for x, color, header, body in zip(xs, colors, headers, bodies):
        ax.add_patch(FancyBboxPatch((x, .19), width, .65, boxstyle='round,pad=.015', facecolor=color, edgecolor='#8A9AA5', linewidth=1))
        ax.text(x + width/2, .785, header, ha='center', va='center', fontsize=18, fontweight='bold', color='#233B4A')
        ax.text(x + width/2, .65, body, ha='center', va='top', fontsize=10, linespacing=1.45)
    for start, end, label in [(xs[0] + width + .012, xs[1] - .012, 'TO BE\nVALIDATED'), (xs[1] + width + .012, xs[2] - .012, 'TO BE\nVALIDATED')]:
        ax.annotate('', xy=(end, .51), xytext=(start, .51), arrowprops=dict(arrowstyle='->', color='#687D88', lw=1.6, linestyle='--'))
        ax.text((start+end)/2, .565, label, ha='center', va='bottom', fontsize=8, color='#536875')
    ax.text(.5, .955, 'Linked biological question • independent age-specific analyses', ha='center', fontsize=17, color='#233B4A')
    ax.text(.5, .115, 'Dashed connections are developmental hypotheses. Cross-sectional scRNA-seq does not trace lineage.', ha='center', fontsize=10, color='#536875')
    ax.text(.5, .06, 'Gates shown use a minimum 50% recovery exploration constraint; later PV protein is not predicted cell by cell.', ha='center', fontsize=9, color='#536875')
    for ext in ['pdf', 'svg', 'png']:
        fig.savefig(OUT / 'figures' / ('DIV30_DIV90_developmental_schematic.' + ext), dpi=600, bbox_inches='tight', facecolor='white')
    plt.close(fig)


def main():
    for d in ['tables', 'figures', 'provenance', 'cache']: (OUT / d).mkdir(parents=True, exist_ok=True)
    s30 = json.loads((D30 / 'summary.json').read_text())
    s90 = json.loads((D90 / 'audit_summary.json').read_text())
    c30 = pd.read_csv(D30 / 'cells.tsv.gz', sep='\t', float_precision='round_trip')
    c90 = pd.read_csv(OLD90 / 'cells.tsv.gz', sep='\t', float_precision='round_trip')
    best30s = read(D30, 'surface_gate_best_at_recovery_floors.tsv')
    b30 = best30s[(best30s.target == 'top20') & np.isclose(best30s.recovery_floor, .5) & (best30s.maximum_markers == 3)].iloc[0].to_dict()
    sweeps90 = read(D90, 'recovery_sweeps.tsv').set_index('objective')
    b90 = sweeps90.loc['pv_purity_recovery_0.5'].to_dict()
    cross, genes = cross_age(c30, c90, b30, b90)
    transferred30 = cross[(cross.rule_selected_at == 'DIV30') & (cross.evaluated_at == 'DIV90')].iloc[0]
    transferred90 = cross[(cross.rule_selected_at == 'DIV90') & (cross.evaluated_at == 'DIV30')].iloc[0]
    gate_temporal = pd.DataFrame([
        dict(rule_selected_at='DIV30', own_age_fold=b30['fold_enrichment'], other_age_fold=transferred30.fold_enrichment,
             temporal_classification='EARLY-TRANSIENT (fixed RNA rule diagnostic)',
             limitation='Distinct age-specific targets; not lineage tracing or protein panel validation.'),
        dict(rule_selected_at='DIV90', own_age_fold=b90['pv_fold_enrichment'], other_age_fold=transferred90.fold_enrichment,
             temporal_classification='LATE-EMERGING (modest fixed RNA rule diagnostic)',
             limitation='Distinct age-specific targets; own-age enrichment remains below1.5x at this recovery.'),
    ])
    gate_temporal.to_csv(OUT / 'tables/selected_gate_temporal_classification.tsv', sep='\t', index=False)
    pars = json.loads((OLD90 / 'state_parameters.json').read_text())['modules']
    score_errors = {}
    for name, p in pars.items():
        z = np.column_stack([np.clip((c30[g] - p['means'][g]) / p['std_population'][g], -3, 3) for g in p['genes']])
        error = float(np.abs(z.mean(axis=1) - c30['direct_' + name + '_score']).max())
        assert error < 2e-6, (name, error)
        score_errors[name] = error
    (OUT / 'provenance/independent_comparison_validation.json').write_text(json.dumps({'div30_cells': len(c30), 'div90_cells': len(c90), 'own_age_selected_gate_counts_match': True, 'cross_age_rules_evaluated_without_refitting': 4, 'frozen_scale_direct_score_max_abs_errors': score_errors}, indent=2) + '\n')
    loss = float((c30.loc[c30.early_pv_score >= c30.early_pv_score.quantile(.8), 'ERBB4'] <= 0).mean())
    temporal = read(D30, 'temporal_surface_marker_classification.tsv')
    marker = read(D30, 'surface_marker_evidence.tsv')
    important = set(genes) | {'FGFR2', 'PTPRS', 'FAT3', 'PTPRM', 'ERBB4'}
    for r in read(D90, 'critical_best_gates.tsv').rule_json: important |= get_genes(json.loads(r))
    for r in best30s[(best30s.target.isin(['top10', 'top20', 'top30'])) & np.isclose(best30s.recovery_floor, .5)].rules:
        important |= {g['gene'] for g in json.loads(r)}
    selected = temporal[temporal.gene.isin(important)].merge(marker[['gene', 'partial_r_early_pv', 'partial_r_neighbor_pv', 'partial_r_early_pv_plus_sample']], on='gene', how='left')
    selected.to_csv(OUT / 'tables/important_marker_temporal_classification.tsv', sep='\t', index=False)
    schematic(s30, b30, b90, loss)
    v30 = read(D30, 'held_sample_surface_gate_summary.tsv')
    held = v30[np.isclose(v30.recovery_floor, .5)].iloc[0]
    old30 = read(D30, 'frozen_DIV90_gate_on_DIV30.tsv').query("target == 'top20'").iloc[0]
    m30 = read(D30, 'surface_gate_confound_models.tsv').set_index('gate_id').loc[b30['gate_id']]
    if isinstance(m30, pd.DataFrame): m30 = m30.iloc[0]
    rows = [
        '# Linked DIV30 / DIV90 purification analysis', '',
        'DIV90 and DIV30 were solved independently. DIV90 asks about purification within a frozen cortical interneuron entry population; DIV30 asks about a possible earlier transcriptional antecedent among all available cells. The analyses do not optimize a pooled gate, and their target percentages have different biological meanings.', '',
        f"The available DIV30 cache contains {len(c30):,} cells from {c30['sample'].nunique()} samples. Three additional samples listed in project metadata are absent from this cache. The historical 9583/9853 sample-prefix alias is reconciled for line and culture summaries, while original barcodes and sample IDs remain unchanged.", '',
        f"**DIV90:** at a minimum 50% strict-PV recovery, the best observed gate yields {b90['pv_fraction']:.1%} PV-biased composition ({b90['pv_fold_enrichment']:.2f}×), retains {int(b90['retained_n']):,} cells and recovers {b90['pv_recovery']:.1%} of PV-biased cells. The old FGFR2/PTPRS gate remains a weak benchmark: SST-biased composition changes only 22.5% → 21.3%; independent statistical significance is not established.", '',
        f"**DIV30:** at a minimum 50% recovery of its own top20 consensus target, the best observed gate yields {b30['post_gate_target_fraction']:.1%} target composition ({b30['fold_enrichment']:.2f}×), with {b30['target_recovery']:.1%} target recovery and {int(b30['retained_n']):,} cells. Pooled held-sample surface selection yields {held.post_gate_target_fraction:.1%} composition ({held.fold_enrichment:.2f}×) and {held.target_recovery:.1%} recovery at the same training constraint, using different gates across folds. These are targets defined by transcriptional resemblance, not confirmed future PV cells.", '',
        f"The DIV30 consensus has Spearman correlation {s30['early_pv_maturation_spearman']:.2f} with neuronal maturation. The selected DIV30 gate's residual correlations are {m30.partial_r_early_pv:.2f} with consensus after measured developmental/depth adjustment, {m30.partial_r_neighbor_pv:.2f} with the independent-feature neighbor score, and {m30.partial_r_early_pv_plus_sample:.2f} with consensus after also adjusting sample. These diagnostics qualify subtype interpretation even when numerical enrichment is strong.", '',
        f"The frozen FGFR2/PTPRS benchmark is more informative for the different DIV30 target: {old30.starting_target_fraction:.1%} → {old30.post_gate_target_fraction:.1%}, {old30.fold_enrichment:.2f}× enrichment and {old30.target_recovery:.1%} recovery. Its DIV30 marker diagnostics primarily reflect progenitor/maturation relationships. This is age-specific diagnostic information, not replication of a validated subtype gate.", '',
        '## Age-specific rule details', '',
        '**DIV30 retention rule** (full precision in the linked table):', '', '```json', json.dumps({'logic': b30['logic'], 'terms': json.loads(b30['rules'])}, indent=2), '```', '',
        '**DIV90 retention rule** (full precision in the linked table):', '', '```json', json.dumps(json.loads(b90['rule_json']), indent=2), '```', '',
        'Each rule was additionally applied, without refitting, to both ages as a diagnostic. The retained target remains age specific. This table tests threshold transfer; it neither assumes common thresholds nor establishes replication from a shared direction.', '',
        markdown_table(cross.drop(columns=['rule_json', 'interpretation'])), '',
        f"Neither selected rule transfers usefully under these fixed RNA thresholds: the DIV30 rule gives {transferred30.fold_enrichment:.2f}× strict-PV enrichment at DIV90, and the DIV90 rule gives {transferred90.fold_enrichment:.2f}× early-consensus enrichment at DIV30. Both are below their receiving-age baseline. The selected DIV30 panel is therefore EARLY-TRANSIENT in this diagnostic; the selected DIV90 panel is LATE-EMERGING with modest own-age purification. These labels describe the two tested problems, not a tracked fate transition.", '',
        '## Temporal marker classification', '',
        'The [important-marker table](tables/important_marker_temporal_classification.tsv) joins independent DIV30 and DIV90 effect sizes with measured developmental confound diagnostics. Labels EARLY-STABLE, EARLY-TRANSIENT, LATE-EMERGING, MATURATION-RELATED and SURFACE-INFORMATIVE BUT NOT SUBTYPE-SPECIFIC are descriptive screening categories. The declared effect-size threshold is |d|≥0.3, with maturation/progenitor evidence taking precedence; it is not a test of developmental lineage or a formal age interaction. A marker unavailable in the old DIV90 surface screen has an unestablished late classification.', '',
        markdown_table(selected[['gene', 'DIV30_top20_cohens_d', 'DIV90_PV_biased_cohens_d', 'DIV30_adjusted_classification', 'temporal_classification']]), '',
        '## Sorting models A–D', '',
        '- **A — same surface gate at both ages:** unsupported for the two selected fixed RNA rules, because both de-enrich the receiving-age target. This does not rule out every other common panel or separately calibrated thresholds.',
        '- **B — early DIV30 gate, later cortical entry, then later subtype gate:** the strongest working hypothesis among these sorting models. Independent DIV30 enrichment, measured-confound residual associations and gene-drop sensitivity support testing an early candidate fraction; the poor rule transfer supports distinct age-specific surface strategies. The later subtype gate still has a substantial recovery/purity tradeoff, and none of the transitions is validated.',
        '- **C — wait until DIV90:** not compelled by the available RNA data, since DIV30 already has an enrichable candidate transcriptional state. Whether waiting improves actual PV production cannot be decided by comparing different age-specific score targets.',
        '- **D — manipulate DIV30 culture instead of sorting:** these cross-sectional data do not compare perturbations and therefore cannot establish that manipulation is more useful. Developmental context can motivate a future test.', '',
        'No model is established as an experimentally useful PV sorting strategy by these RNA computations alone. The separate reports identify which computational patterns support further assessment; no antibody recommendations are made.', '',
        '## Developmental schematic and complete results', '',
        '[Developmental schematic, PDF](figures/DIV30_DIV90_developmental_schematic.pdf) · [editable SVG](figures/DIV30_DIV90_developmental_schematic.svg) · [600-dpi PNG](figures/DIV30_DIV90_developmental_schematic.png). Both connections are labeled TO BE VALIDATED. The later PVALB protein endpoint is the user-reported system observation; individual profiled cells are not traced to it.', '',
        '[Complete DIV90 report](../div90_purification_audit_v1/REPORT.md) · [Complete DIV30 report](../div30_pv_antecedent_v1/REPORT.md). Original requests and SHA256/source preservation checks are in `provenance/`. The original PI package remains unchanged.', '',
    ]
    (OUT / 'REPORT.md').write_text('\n'.join(rows))
    (OUT / 'comparison_summary.json').write_text(json.dumps({'separate_optimization': True, 'div30_top20_erbb4_rna_loss': loss, 'div90_selected_objective': 'strict PV purity at >=50% PV recovery', 'div30_selected_objective': 'top20 consensus purity at >=50% target recovery', 'cross_age_rule_tests': len(cross), 'important_markers': len(selected), 'lineage_established': False, 'antibodies_recommended': False}, indent=2) + '\n')
    print('Linked report, temporal comparison, four cross-age diagnostic masks and schematic complete', flush=True)


if __name__ == '__main__': main()
