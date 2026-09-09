#!/usr/bin/env python3
"""Render and independently verify the completed purification-first DIV90 audit."""
from __future__ import annotations

import json
import os
from pathlib import Path
import textwrap

ROOT = Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results')
OUT = ROOT / 'div90_purification_audit_v1'
SOURCE = ROOT / 'div90_hypergate_sst_pv_phase2'
os.environ['MPLCONFIGDIR'] = str(OUT / 'cache/matplotlib')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np
import pandas as pd

COLORS = ['#176B9B', '#A879BB', '#D1B9DC', '#DD7250', '#B6BFC5']
LABELS = ['PV-biased', 'PV-supported dual-high', 'Other dual-high', 'SST-biased', 'Unresolved']


def read(name):
    return pd.read_csv(OUT / 'tables' / name, sep='\t', float_precision='round_trip')


def save(fig, name):
    for suffix, kwargs in [('pdf', {}), ('svg', {}), ('png', {'dpi': 600})]:
        fig.savefig(OUT / 'figures' / (name + '.' + suffix), bbox_inches='tight', facecolor='white', **kwargs)
    plt.close(fig)


def mask_from_json(cells, expression):
    """Independent recursive rule evaluator: True always means retained."""
    if 'gene' in expression:
        v = cells[expression['gene']].to_numpy(dtype=np.float64)
        threshold = np.float64(expression['threshold'])
        if expression['op'] == '>': return v > threshold
        if expression['op'] == '>=': return v >= threshold
        if expression['op'] == '<': return v < threshold
        if expression['op'] == '<=': return v <= threshold
        raise ValueError(expression['op'])
    if 'not' in expression:
        return np.logical_not(mask_from_json(cells, expression['not']))
    result = np.full(len(cells), expression['logic'] == 'AND', dtype=bool)
    for child in expression['terms']:
        child_mask = mask_from_json(cells, child)
        result = result & child_mask if expression['logic'] == 'AND' else result | child_mask
    return result


def baseline_and_benchmarks(cells, targets):
    summary = json.loads((SOURCE / 'gate_summary.json').read_text())
    masks = {'Starting population': np.ones(len(cells), bool)}
    for name, key in [('FGFR2/PTPRS benchmark', 'experimental_depletion'), ('FAT3/PTPRM benchmark', 'phase1_reference')]:
        r = summary['rules'][summary[key]['gate_id']]
        e = {'logic': r['logic'], 'terms': r['rules']}
        if r['action'] == 'remove': e = {'not': e}
        masks[name] = mask_from_json(cells, e)
    rows = []
    for name, m in masks.items():
        row = {'objective': name, 'retained_n': int(m.sum())}
        for short, label in [('pv', 'PV-biased'), ('dual', 'PV/SST hybrid'), ('sst', 'SST-biased'), ('unresolved', 'unresolved/immature')]:
            y = cells.state.eq(label).to_numpy()
            row[short + '_fraction'] = float(y[m].mean())
            row[short + '_recovery'] = float((y & m).sum() / y.sum())
        hc = targets.hc_dual.to_numpy(bool)
        row['hc_dual_fraction'] = float(hc[m].mean())
        row['hc_dual_recovery'] = float((hc & m).sum() / hc.sum())
        row['supported_fraction'] = row['pv_fraction'] + row['hc_dual_fraction']
        row['pv_fold_enrichment'] = row['pv_fraction'] / cells.state.eq('PV-biased').mean()
        rows.append(row)
    d = pd.DataFrame(rows)
    d.to_csv(OUT / 'tables/independently_verified_benchmarks.tsv', sep='\t', index=False)
    assert int(d.iloc[1].retained_n) == 4024
    assert abs(d.iloc[1].sst_fraction - 856 / 4024) < 1e-12
    return d


def independent_validation(cells, targets, sweep):
    rows = []
    for row in sweep.drop_duplicates(['group', 'source_row']).itertuples():
        mask = mask_from_json(cells, json.loads(row.rule_json))
        assert mask.sum() == row.retained_n
        for short, label in [('pv', 'PV-biased'), ('dual', 'PV/SST hybrid'), ('sst', 'SST-biased'), ('unresolved', 'unresolved/immature')]:
            assert (mask & cells.state.eq(label)).sum() == getattr(row, 'retained_' + short + '_n')
        for short in ['hc_dual', 'hc_loose', 'hc_strict', 'continuous10', 'continuous20', 'continuous30']:
            assert (mask & targets[short].to_numpy(bool)).sum() == getattr(row, 'retained_' + short + '_n')
        assert abs(cells.loc[mask, 'sst_score'].mean() - row.sst_score_mean) < 1e-10
        rows.append({'group': row.group, 'row': int(row.source_row), 'retained_n': int(mask.sum()), 'all_counts_and_sst_mean_match': True})
    result = {'independent_implementation': True, 'rules_checked': len(rows), 'results': rows}
    (OUT / 'provenance/report_independent_validation.json').write_text(json.dumps(result, indent=2) + '\n')


def figures(critical, sweep, bench):
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                         'svg.fonttype': 'none', 'pdf.fonttype': 42})
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    front = read('pareto_pv.tsv')
    ax = axes[0]
    s = ax.scatter(front.pv_recovery, front.pv_fraction, c=front.sst_fraction, s=15 + front.retained_n / 35, cmap='magma_r', vmin=0, vmax=.35, edgecolors='white', linewidths=.4, zorder=3)
    ax.plot(front.pv_recovery, front.pv_fraction, color='#9AA7B1', lw=.8, zorder=1)
    for b in bench.itertuples():
        marker = 'D' if b.objective == 'Starting population' else '*' if b.objective.startswith('FGFR2') else 's'
        ax.scatter(b.pv_recovery, b.pv_fraction, marker=marker, s=100, facecolors='none', edgecolors='#172B3A', label=b.objective, zorder=4)
    ax.set(xlabel='PV-biased recovery', ylabel='PV-biased retained purity', title='Full observed PV purity / recovery frontier', xlim=(-.015, 1.025), ylim=(0, 1.03))
    ax.xaxis.set_major_formatter(PercentFormatter(1)); ax.yaxis.set_major_formatter(PercentFormatter(1))
    ax.legend(loc='upper right', fontsize=8, frameon=False)
    fig.colorbar(s, ax=ax, label='Retained SST-biased fraction', shrink=.65)
    ax.text(.02, .03, 'Point area increases with total retained cells', transform=ax.transAxes, fontsize=8)
    ax = axes[1]
    for nm, color in [(1, '#6999B5'), (2, '#A57AB6'), (3, '#C47B46')]:
        d = sweep[sweep.objective.str.startswith(f'markers_{nm}_pv_recovery_')].copy()
        d['floor'] = d.objective.str.rsplit('_', n=1).str[-1].astype(float)
        d = d[d.floor >= .1].sort_values('floor')
        ax.plot(d.floor, d.pv_fraction, 'o-', label=f'{nm} marker' + ('s' if nm > 1 else ''), color=color)
    ax.axhline(bench.iloc[0].pv_fraction, color='#697780', linestyle=':', label='Starting purity')
    ax.set(xlabel='Minimum PV-biased recovery required', ylabel='Maximum observed PV-biased purity', title='Does a third marker improve useful recovery?')
    ax.xaxis.set_major_formatter(PercentFormatter(1)); ax.yaxis.set_major_formatter(PercentFormatter(1)); ax.legend(frameon=False)
    save(fig, '01_pv_purification_frontier')

    comp = pd.concat([bench.iloc[:2], critical], ignore_index=True)
    objective_labels = {
        'A': 'Highest PV purity', 'B': 'Highest PV purity; ≥50% PV recovery',
        'C': 'Highest PV purity; ≥70% PV recovery', 'D': 'Highest PV-supported purity',
        'E': 'Lowest SST contamination', 'F': 'Highest PV:SST ratio',
        'G': 'Best 1-marker; ≥10% PV recovery', 'H': 'Best 2-marker; ≥10% PV recovery',
        'I': 'Best 3-marker; ≥10% PV recovery',
    }
    labels = ['Starting population', 'FGFR2/PTPRS'] + [f'{r.criterion}: {objective_labels[r.criterion]}' for r in critical.itertuples()]
    fig, ax = plt.subplots(figsize=(13, 7), constrained_layout=True)
    left = np.zeros(len(comp))
    vals = np.column_stack([comp.pv_fraction, comp.hc_dual_fraction, comp.dual_fraction - comp.hc_dual_fraction, comp.sst_fraction, comp.unresolved_fraction])
    for j, (name, color) in enumerate(zip(LABELS, COLORS)):
        ax.barh(np.arange(len(comp)), vals[:, j], left=left, color=color, label=name, height=.7)
        left += vals[:, j]
    ax.set_yticks(np.arange(len(comp)), [textwrap.fill(s, 40) for s in labels], fontsize=8)
    ax.invert_yaxis(); ax.set(xlim=(0, 1.48), xlabel='Retained composition (all retained cells in denominator)', title='Best gate for each objective: composition and recovery')
    ax.set_xticks([0, .25, .5, .75, 1]); ax.xaxis.set_major_formatter(PercentFormatter(1))
    for i, r in enumerate(comp.itertuples()):
        ax.text(1.015, i, f'n={int(r.retained_n):,}  |  PV recovery {r.pv_recovery:.1%}', va='center', fontsize=9)
    ax.legend(loc='upper center', bbox_to_anchor=(.45, -.09), ncol=3, frameon=False, fontsize=9)
    save(fig, '02_best_objective_compositions')

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    ax = axes[0]
    for pct, color in [(10, '#286F8E'), (20, '#A05AA5'), (30, '#C9803F')]:
        d = read(f'pareto_continuous{pct}.tsv')
        ax.plot(d[f'continuous{pct}_recovery'], d[f'continuous{pct}_fold_enrichment'], label=f'Top {pct}% preference', color=color)
    ax.axhline(1.5, color='#999999', lw=.7, linestyle=':'); ax.axhline(2, color='#777777', lw=.7, linestyle='--')
    ax.set(xlabel='Continuous-target recovery', ylabel='Fold enrichment of target', title='Scaled high-PV / low-SST continuous targets', xlim=(0, 1.01))
    ax.xaxis.set_major_formatter(PercentFormatter(1)); ax.legend(frameon=False)
    ax = axes[1]
    d = read('pareto_sst_removal.tsv')
    s = ax.scatter(d.sst_removal, d.pv_loss, c=d.hc_dual_loss, s=25, cmap='viridis', vmin=0, vmax=1)
    fig.colorbar(s, ax=ax, label='PV-supported dual-high fraction lost', shrink=.65)
    ax.plot(d.sst_removal, d.pv_loss, color='#467A8C', lw=.8)
    ax.plot([0, 1], [0, 1], linestyle=':', color='#999999', label='Equal proportional loss')
    ax.set(xlabel='SST-biased fraction removed', ylabel='PV-biased fraction lost', title='Full observed SST removal / PV loss frontier', xlim=(0, 1), ylim=(0, 1))
    ax.xaxis.set_major_formatter(PercentFormatter(1)); ax.yaxis.set_major_formatter(PercentFormatter(1)); ax.legend(frameon=False)
    save(fig, '03_continuous_targets_and_sst_depletion')


def report(critical, sweep, bench, summary):
    ix = sweep.set_index('objective')
    pick = lambda name: ix.loc[name]
    pv10, pv50, pv70 = [pick(f'pv_purity_recovery_{v}') for v in ['0.1', '0.5', '0.7']]
    sup10 = pick('supported_purity_recovery_0.1')
    absolute = pick('pv_purity_recovery_0.0')
    pair = pick('markers_2_pv_recovery_0.1'); triple = pick('markers_3_pv_recovery_0.1')
    enrich15, enrich2 = pick('pv_enrichment_1.5'), pick('pv_enrichment_2.0')
    below15, below10 = pick('sst_below_0.15'), pick('sst_below_0.1')
    targets = json.loads((OUT / 'target_definitions.json').read_text())
    lines = [
        '# DIV90 purification audit', '',
        'This audit resumes the completed Phase 2 PI package. It preserves the frozen 4,768 cells, their original scores, labels and coordinates, and all previous results. It ranks actual retained composition at explicit recovery floors rather than F1 or preservation-heavy scores.', '',
        f"**Main numerical result:** the best observed strict PV-biased purity is {pv10.pv_fraction:.1%} at a minimum 10% PV recovery ({pv10.pv_recovery:.1%} achieved; {int(pv10.retained_n):,} retained cells). At minimum 50% and 70% recovery, the optima are {pv50.pv_fraction:.1%} and {pv70.pv_fraction:.1%}. Starting strict PV-biased composition is {bench.iloc[0].pv_fraction:.1%}.", '',
        '**Exact retention rule for the ≥10% recovery optimum:**', '', '```text', pv10.exact_retention_rule, '```', '',
        f"FGFR2/PTPRS remains a weak benchmark: PV-biased composition {bench.iloc[0].pv_fraction:.1%} → {bench.iloc[1].pv_fraction:.1%}, SST-biased composition {bench.iloc[0].sst_fraction:.1%} → {bench.iloc[1].sst_fraction:.1%}, PV recovery {bench.iloc[1].pv_recovery:.1%}. Independent statistical significance is not established.", '',
        '## Search coverage and interpretation', '',
        f"Every one of the **{summary['legacy_evaluated_rows']:,} existing rows** was reconstructed. The expanded audit evaluates **{summary['total_evaluated_rows']:,} rows**, including repeated phenotypes and frozen R-fit benchmarks. Exact signed single thresholds, signed grid pairs and bounded hierarchical three-marker extensions are separately indexed in `tables/search_inventory.tsv`. Every row has count, composition, recovery, yield, fold and score metrics in `tables/all_purification_candidates.h5`.", '',
        'Absolute maxima and useful-recovery screening are distinct. A 10% recovery floor is an explicit exploration convention, not evidence that a sort is biologically or practically useful. All results are selected in the same data used to evaluate them; fixed sample/line/condition summaries are descriptive and do not provide independent validation.', '',
        f"The operational high-confidence PV-supported dual-high subset contains **{targets['baseline_counts']['hc_dual']} of 1,309** original dual-high cells. It requires ≥2 detected PV genes, ≥2 non-SST SST-module genes, PV score ≥ the pooled 65th percentile, and absence of the saved exploratory doublet/high-RNA flags. Looser, stricter and no-QC versions are retained as sensitivity channels. These are computational support filters, not independently validated hybrid identities.", '',
        'Continuous preference is `min(midrank_percentile(PV), 1 − midrank_percentile(SST))`. This requires both high PV and low SST on comparable percentile scales; raw PV−SST is not used as the primary ranking. Whole quantile ties are retained. Rank-difference targets provide a scaling sensitivity.', '',
        '## Critical objective table', '',
        'The complete requested table is [critical_best_gates.tsv](tables/critical_best_gates.tsv), including the exact full-precision retention expression, retained cells, PV/HC-dual/all-dual/SST/unresolved fractions, baseline fractions beside every row, recoveries and fold enrichment. The abbreviated view below shows its main outcomes. HC dual-high is a subset of all dual-high and must not be added to it twice. G/H/I maximize strict PV purity at ≥10% PV recovery; all other recovery floors and unconstrained maxima are in [recovery_sweeps.tsv](tables/recovery_sweeps.tsv).', '',
        '| Criterion | Retained cells | PV purity | HC dual | All dual | SST | Unresolved | PV recovery | HC recovery | SST recovery | PV fold |',
        '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|',
    ]
    for r in critical.itertuples():
        lines.append(f'| {r.criterion} | {int(r.retained_n):,} | {r.pv_fraction:.1%} | {r.hc_dual_fraction:.1%} | {r.dual_fraction:.1%} | {r.sst_fraction:.1%} | {r.unresolved_fraction:.1%} | {r.pv_recovery:.1%} | {r.hc_dual_recovery:.1%} | {r.sst_recovery:.1%} | {r.pv_fold_enrichment:.2f}× |')
    lines += ['', 'A highest PV purity; B highest PV purity at ≥50% recovery; C at ≥70%; D highest PV plus HC-dual purity; E lowest retained SST fraction (PV recovery breaks ties); F highest PV:SST ratio; G/H/I best one/two/three-marker purity at ≥10% PV recovery.', '',
              '## Continuous PV-preference target results', '',
              'Representative purity optima for each independent high-PV/low-SST quantile. Actual starting fractions include whole ties. Scores below are the unchanged independent SST module score, averaged over the starting or retained population. Full exact rules and every recovery floor remain in `recovery_sweeps.tsv`.', '',
              '| Target | Recovery floor | Starting target | Retained target | Fold | Target recovery | SST score before → after | Total yield | Retained cells |',
              '|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for pct in [10, 20, 30]:
        for floor in [.1, .5]:
            target = f'continuous{pct}'
            r = pick(f'{target}_purity_recovery_{floor:.1f}')
            lines.append(f"| Top {pct}% | {floor:.0%} | {r['baseline_'+target+'_fraction']:.2%} | {r[target+'_fraction']:.1%} | {r[target+'_fold_enrichment']:.2f}× | {r[target+'_recovery']:.1%} | {r.baseline_sst_score_mean:.3f} → {r.sst_score_mean:.3f} | {r.total_yield:.1%} | {int(r.retained_n):,} |")
    lines += ['',
              '## Answers to the six DIV90 questions', '',
              '1. **Did preservation weighting hide better purification?** Yes, it prioritized retaining PV-biased and all dual-high cells. The independent legacy-table scan alone found 28.56% strict PV purity at ≥50% PV recovery, versus 23.71% with FGFR2/PTPRS. The stronger low-recovery results additionally require signed and hierarchical logic that was absent or restricted in the earlier search. Both the ranking objective and the explored rule set matter.', '',
              f'2. **Is stronger strict PV purification possible at lower recovery?** The expanded search reaches {pv10.pv_fraction:.1%} purity ({pv10.pv_fold_enrichment:.2f}×) at {pv10.pv_recovery:.1%} PV recovery. The absolute maximum is {absolute.pv_fraction:.1%}, but captures only {int(absolute.retained_n)} cells and {absolute.pv_recovery:.2%} of PV-biased cells. Judge these together rather than reporting the maximum alone.', '',
              f'3. **Are the strongest gates enriching strict PV rather than ambiguous dual-high?** The ≥10%-recovery strict-PV optimum contains {pv10.pv_fraction:.1%} strict PV, {pv10.dual_fraction:.1%} all dual-high and {pv10.hc_dual_fraction:.1%} conservative HC dual-high. The separate supported-state optimum at ≥10% supported recovery reaches {sup10.supported_fraction:.1%} supported purity, comprising {sup10.pv_fraction:.1%} strict PV plus {sup10.hc_dual_fraction:.1%} HC dual-high. The original quadrant labels themselves remain proxies.', '',
              '4. **Can SST removal avoid catastrophic PV loss?** The full removal/loss Pareto frontier and fixed removal thresholds are reported below; no preservation-first objective is used to choose the removal amount.', '']
    lines += ['| SST removal requested | SST removed | PV lost | HC dual lost | Retained SST composition | Cells retained |', '|---:|---:|---:|---:|---:|---:|']
    for v in [.2, .3, .4, .5, .6, .7, .8]:
        r = pick(f'sst_removal_{v:.1f}')
        lines.append(f'| {v:.0%} | {r.sst_removal:.1%} | {r.pv_loss:.1%} | {r.hc_dual_loss:.1%} | {r.sst_fraction:.1%} | {int(r.retained_n):,} |')
    lines += ['', f'5. **Is there a useful three-marker improvement?** At the same ≥10% PV recovery constraint, the best two-marker purity is {pair.pv_fraction:.1%}, and the best three-marker purity is {triple.pv_fraction:.1%}: an increase of {100*(triple.pv_fraction-pair.pv_fraction):.1f} percentage points. Compare the complete marker-count curves at higher recovery as well. The triplet search is a bounded extension of purity/recovery/removal-selected pairs, so an untested triplet cannot be ruled out.', '',
              '6. **What is the theoretical best purification?** These results establish an observed optimum over the evaluated RNA rules, not a mathematical bound on every threshold, every three-marker expression or the entire protein surfaceome. Exact single boundaries are exhaustive for the supplied genes; pairs use a threshold grid and triplets are bounded. Prospective protein sorting and later PV phenotype remain unvalidated.', '',
              '## Translational exploration thresholds', '',
              f"The greatest strict-PV recovery reaching ≥1.5× enrichment is {enrich15.pv_recovery:.1%} ({enrich15.pv_fraction:.1%} purity). A ≥2× rule reaches {enrich2.pv_fraction:.1%} purity at {enrich2.pv_recovery:.2%} PV recovery—approximately 10%, although just below the table's exact 10% screen. These exploration cutoffs should not turn a small numerical threshold difference into a biological success/failure distinction.", '',
              f"Reducing retained SST-biased composition below 15% can preserve at most {below15.pv_recovery:.1%} of PV-biased cells among evaluated rules; below 10% can preserve {below10.pv_recovery:.1%}. The corresponding PV losses are {1-below15.pv_recovery:.1%} and {1-below10.pv_recovery:.1%}. This distinguishes removing many SST cells from obtaining a low-SST retained population.", '',
              '[threshold_attainment.tsv](tables/threshold_attainment.tsv) also tests 65%/70%/75% composition, reporting the greatest achievable recovery and whether it reaches the stated 10% exploration floor. Small composition shifts remain weak regardless of significance.', '',
              'For the historical PV-plus-all-dual target, starting composition is 50.02097%, so its mathematical maximum fold enrichment is 1.99916×. An exact ≥2.0× flag is therefore impossible even for a 100%-target capture. Interpret the approximate biological threshold using purity and retained counts rather than that rounding boundary.', '',
              '## Figures and reproducibility', '',
              '- [PV purity/recovery Pareto frontier and marker complexity](figures/01_pv_purification_frontier.pdf): SST contamination is color; total retained cells determine point size.',
              '- [Composition for the nine requested objectives](figures/02_best_objective_compositions.pdf).',
              '- [Continuous targets and SST removal/PV loss frontier](figures/03_continuous_targets_and_sst_depletion.pdf).', '',
              'Each figure is also supplied as editable SVG and 600-dpi PNG. The entire numeric Pareto frontier for every primary objective is in `tables/pareto_*.tsv`. Fixed sample, cell-line and culture-condition results are in `tables/selected_gate_fixed_sample_evaluations.tsv`; CV/MW remains confounded with other culture differences, and 2E/CV has only 37 entry cells.', '',
              'The discovery script verifies all legacy categorical counts and directly reconstructs selected masks. A separate report implementation checks selected rule counts, HC/continuous targets and SST means in `provenance/report_independent_validation.json`. The comparison package verifies preservation of the original PI package and all earlier scientific results. No antibodies are recommended by this audit.', '']
    (OUT / 'REPORT.md').write_text('\n'.join(lines))
    (OUT / 'figures/FIGURE_CAPTIONS.md').write_text('\n'.join([
        '# DIV90 purification figures', '',
        'Figure 1. All nondominated strict-PV purity/recovery points among evaluated rules. Color is retained SST-biased contamination; area increases linearly with total retained cells (with a small visibility offset). Hollow points mark frozen benchmarks. Marker-complexity curves optimize purity at identical requested recovery floors. None are held-out estimates.', '',
        'Figure 2. Retained composition for baseline, FGFR2/PTPRS and criteria A–I. Conservative PV-supported dual-high cells are separated from other dual-high cells so stacked categories sum to 100%. Retained cell numbers and PV recoveries expose tiny-capture optima. G/H/I use a 10% recovery exploration floor.', '',
        'Figure 3. Left: full observed enrichment/recovery frontiers for three continuous PV-preference quantiles. Right: full observed SST-removal/PV-loss frontier; color denotes conservative HC-dual loss from 0 (dark) to 1 (bright). The diagonal is equal proportional PV and SST loss. All comparisons are cross-sectional, in-sample RNA calculations.', '']))


def main():
    (OUT / 'figures').mkdir(parents=True, exist_ok=True)
    summary = json.loads((OUT / 'audit_summary.json').read_text())
    critical, sweep = read('critical_best_gates.tsv'), read('recovery_sweeps.tsv')
    cells = pd.read_csv(SOURCE / 'cells.tsv.gz', sep='\t', float_precision='round_trip')
    targets = read('target_cell_assignments.tsv.gz')
    assert cells.cell_id.equals(targets.cell_id)
    independent_validation(cells, targets, sweep)
    bench = baseline_and_benchmarks(cells, targets)
    figures(critical, sweep, bench)
    report(critical, sweep, bench, summary)
    print('DIV90 report, three figure sets and independent selected-rule checks complete', flush=True)


if __name__ == '__main__':
    main()
