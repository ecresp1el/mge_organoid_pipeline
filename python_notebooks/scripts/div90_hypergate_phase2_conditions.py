#!/usr/bin/env python3
"""Descriptive CV/MW comparisons preserving cell lines, state mixture and gates.

The user-supplied glucose/operator mapping supplements technical sample IDs.
All outputs are restricted to the separate Turbo phase-two runtime directory.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from div90_hypergate_paths import PROJECT_ROOT, REPO_ROOT
from div90_hypergate_phase2_gates import STATES, SHORT, apply_rule, clean_json, metrics

OUT = PROJECT_ROOT / 'results/div90_hypergate_sst_pv_phase2'
REFERENCE_GENES = ['ERBB4', 'FAT3', 'PTPRM', 'CXCR4', 'ACKR3', 'NRP2']
GROUPS = [[], ['culture_operator'], ['cell_line'], ['culture_operator', 'cell_line'],
          ['culture_operator', 'cell_line', 'sample']]
GATE_NAMES = ['baseline', 'phase1_reference', 'best_single', 'best_pair', 'best_three',
              'selected', 'erbb4_cxcr4', 'practical_single', 'practical_pair', 'experimental_depletion']


def write_json(path, data):
    path.write_text(json.dumps(clean_json(data), indent=2, allow_nan=False) + '\n')


def groups(cells):
    for keys in GROUPS:
        iterator = [((), cells)] if not keys else cells.groupby(keys, observed=True, sort=True)
        for values, d in iterator:
            values = values if isinstance(values, tuple) else (values,)
            meta = dict(grouping='all' if not keys else '|'.join(keys),
                        culture_operator='all', cell_line='all', sample='all')
            meta.update(zip(keys, values))
            yield meta, d


def counts(d):
    return np.asarray([(d.state == s).sum() for s in STATES], dtype=int)


def verify_metadata(cells, out):
    path = REPO_ROOT / 'metadata/div30_div90_sample_id_to_biolabel_map.tsv'
    source = pd.read_csv(path, sep='\t').query("DIV == 'DIV90'").set_index('run_sample_id')
    expected = {f'10496-MW-{i}': ('CV' if i % 2 else 'MW') for i in range(1, 7)}
    if len(cells) != 4768 or not cells.cell_id.is_unique:
        raise AssertionError('Condition analysis requires the frozen 4,768 unique cells')
    if set(cells['sample']) != set(expected) or set(source.index) != set(expected):
        raise AssertionError('User condition mapping does not match recorded six sample IDs')
    observed = cells[['sample', 'biological_label', 'cell_line', 'culture_operator', 'condition']].drop_duplicates()
    if len(observed) != 6:
        raise AssertionError('Expected one metadata assignment per sample')
    for row in observed.itertuples():
        if row.culture_operator != expected[row.sample]:
            raise AssertionError(f'Culture operator mismatch for {row.sample}')
        label = source.loc[row.sample, 'biological_label']
        if row.biological_label != label or row.cell_line != label.rsplit('_rep', 1)[0]:
            raise AssertionError(f'Cell-line mapping changed for {row.sample}')
    observed = observed.sort_values('sample')
    observed['relative_glucose_approx'] = observed.culture_operator.map({'CV': 2., 'MW': 1.})
    observed['condition_source'] = 'User experimental metadata supplied in this conversation, 2026-09-09'
    observed.to_csv(out / 'tables/condition_metadata_verified.tsv', sep='\t', index=False)
    return dict(sample_mapping=observed.to_dict('records'), source_map=str(path),
                source_map_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                all_six_ids_match=True, no_phase1_metadata_changed=True,
                technical_id_note='MW in every original run ID is not culture-operator metadata.',
                interpretation='CV higher glucose (~2x) versus MW lower glucose (~1x); operator and glucose are confounded.',
                replication='One CV and one MW sample per cell line; three paired cell-line contrasts, descriptive only.')


def summarize_baseline(cells, genes, out):
    rows = []
    for meta, d in groups(cells):
        row = dict(**meta, n=len(d), sst_rna_detected_n=int(d.SST.gt(0).sum()),
                   sst_rna_detected_fraction=d.SST.gt(0).mean(),
                   pv_score_median=d.pv_score.median(), sst_score_median=d.sst_score.median(),
                   pv_score_mean=d.pv_score.mean(), sst_score_mean=d.sst_score.mean())
        for program in ['pv', 'sst']:
            score = d[program + '_score']
            row.update({program + '_score_min': score.min(), program + '_score_max': score.max(),
                        program + '_score_q10': score.quantile(.1), program + '_score_q90': score.quantile(.9)})
        for key in ['total_counts', 'n_genes', 'mito_fraction', 'stress_score']:
            if key in d:
                row[key + '_median'] = d[key].median()
        for state, short in zip(STATES, SHORT):
            row[short + '_n'] = int(d.state.eq(state).sum())
            row[short + '_fraction'] = d.state.eq(state).mean()
        for gene in genes:
            row[gene + '_mean_logexpr'] = d[gene].mean()
            row[gene + '_detected_fraction'] = d[gene].gt(0).mean()
        rows.append(row)
    table = pd.DataFrame(rows)
    table.to_csv(out / 'tables/condition_baseline.tsv', sep='\t', index=False)
    return table


def summarize_gates(cells, summary, out):
    rows = []
    for name in GATE_NAMES:
        gate = summary.get(name)
        if not isinstance(gate, dict):
            continue
        rule = summary['rules'][gate['gate_id']]
        removed = pd.Series(apply_rule(cells, rule), index=cells.index)
        for meta, d in groups(cells):
            totals = counts(d)
            rm = removed.loc[d.index].to_numpy()
            removed_counts = counts(d.loc[rm])
            if not np.array_equal(removed_counts + counts(d.loc[~rm]), totals):
                raise AssertionError('Gate fractions do not reconstruct their baseline')
            row = dict(**meta, gate_name=name, gate_id=gate['gate_id'],
                       gate_label=gate['gate_label'], baseline_n=len(d),
                       **{k: float(v[0]) for k, v in metrics(removed_counts, totals).items()})
            for state, short in zip(STATES, SHORT):
                row['baseline_' + short + '_fraction'] = d.state.eq(state).mean()
            for fraction, subset in [('removed', d.loc[rm]), ('retained', d.loc[~rm])]:
                row[fraction + '_pv_score_median'] = subset.pv_score.median()
                row[fraction + '_sst_score_median'] = subset.sst_score.median()
                row[fraction + '_sst_rna_detected_fraction'] = subset.SST.gt(0).mean() if len(subset) else np.nan
            rows.append(row)
    table = pd.DataFrame(rows)
    table.to_csv(out / 'tables/condition_gate_metrics.tsv', sep='\t', index=False)
    return table


def paired_table(table, index_keys, excluded, outpath):
    """CV minus MW by line; never a cell-level significance test."""
    d = table[table.grouping == 'culture_operator|cell_line']
    measures = [x for x in d.select_dtypes(include=np.number).columns if x not in excluded]
    rows = []
    for key, group in d.groupby(index_keys, observed=True, sort=True):
        key = key if isinstance(key, tuple) else (key,)
        cv, mw = group[group.culture_operator == 'CV'], group[group.culture_operator == 'MW']
        if len(cv) != 1 or len(mw) != 1:
            raise AssertionError('A paired contrast must have exactly one CV and one MW value')
        for measure in measures:
            a, b = float(cv.iloc[0][measure]), float(mw.iloc[0][measure])
            rows.append(dict(zip(index_keys, key), metric=measure, CV=a, MW=b, CV_minus_MW=a-b))
    result = pd.DataFrame(rows)
    result.to_csv(outpath, sep='\t', index=False)
    return result


def marker_analysis(cells, genes, out):
    """Within-state expression and state associations, pooled and in each line."""
    effects, standardized, shifts = [], [], []
    scopes = [('all', cells)] + list(cells.groupby('cell_line', observed=True, sort=True))
    for line, scope in scopes:
        state_mean, state_det, state_prop = {}, {}, {}
        for operator, d in scope.groupby('culture_operator', observed=True, sort=True):
            state_mean[operator] = d.groupby('state', observed=True)[genes].mean().reindex(STATES)
            state_det[operator] = d[genes].gt(0).groupby(d.state, observed=True).mean().reindex(STATES)
            state_prop[operator] = d.state.value_counts(normalize=True).reindex(STATES, fill_value=0)
            supported = d.state.value_counts().reindex(STATES, fill_value=0)
            pvr = d[genes].corrwith(d.pv_score, method='spearman')
            sstr = d[genes].corrwith(d.sst_score, method='spearman')
            for gene in genes:
                means, detection = state_mean[operator][gene], state_det[operator][gene]
                complete = bool(means.notna().all())
                standardized.append(dict(cell_line=line, culture_operator=operator, gene=gene, n=len(d),
                    n_supported_states=int(means.notna().sum()), minimum_state_n=int(supported.min()),
                    raw_mean_logexpr=d[gene].mean(), raw_detected_fraction=d[gene].gt(0).mean(),
                    equal_state_mean_logexpr=means.mean() if complete else np.nan,
                    equal_state_detected_fraction=detection.mean() if complete else np.nan,
                    rho_pv_score=pvr[gene], rho_sst_score=sstr[gene]))
                for state in STATES:
                    a, b = d.loc[d.state.eq(state), gene], d.loc[~d.state.eq(state), gene]
                    other = [s for s in STATES if s != state]
                    complete_rest = means.loc[other].notna().all()
                    effects.append(dict(cell_line=line, culture_operator=operator, gene=gene, state=state,
                        n_state=len(a), n_other=len(b), mean_logexpr_state=a.mean(), mean_logexpr_other=b.mean(),
                        mean_difference=a.mean()-b.mean(), detection_state=a.gt(0).mean() if len(a) else np.nan,
                        detection_other=b.gt(0).mean() if len(b) else np.nan,
                        detection_difference=(a.gt(0).mean()-b.gt(0).mean()) if len(a) and len(b) else np.nan,
                        standardized_other_mean=means.loc[other].mean() if complete_rest else np.nan,
                        standardized_mean_difference=a.mean()-means.loc[other].mean() if complete_rest else np.nan,
                        standardized_detection_difference=detection.loc[state]-detection.loc[other].mean() if complete_rest else np.nan))
        for gene in genes:
            for endpoint, value in [('mean_logexpr', state_mean), ('detected_fraction', state_det)]:
                cv, mw = value['CV'][gene], value['MW'][gene]
                pcv, pmw = state_prop['CV'], state_prop['MW']
                complete = bool(cv.notna().all() and mw.notna().all())
                xcv, xmw = scope.loc[scope.culture_operator.eq('CV'), gene], scope.loc[scope.culture_operator.eq('MW'), gene]
                if endpoint == 'detected_fraction':
                    xcv, xmw = xcv.gt(0), xmw.gt(0)
                raw = xcv.mean()-xmw.mean()
                mix = float(((pcv-pmw)*(cv+mw)/2).sum()) if complete else np.nan
                within = float(((cv-mw)*(pcv+pmw)/2).sum()) if complete else np.nan
                if complete and not np.isclose(raw, mix+within, atol=1e-10):
                    raise AssertionError('Mixture decomposition did not reconstruct raw difference')
                shifts.append(dict(cell_line=line, gene=gene, endpoint=endpoint,
                    complete_four_state_support=complete, raw_CV_minus_MW=raw,
                    equal_state_CV_minus_MW=(cv-mw).mean() if complete else np.nan,
                    composition_component=mix, within_state_component=within,
                    reconstruction_error=raw-mix-within if complete else np.nan))
    effects = pd.DataFrame(effects)
    standard = pd.DataFrame(standardized)
    shifts = pd.DataFrame(shifts)
    effects.to_csv(out / 'tables/condition_marker_state_effects.tsv.gz', sep='\t', index=False)
    standard.to_csv(out / 'tables/condition_marker_state_standardized.tsv', sep='\t', index=False)
    shifts.to_csv(out / 'tables/condition_marker_shift_decomposition.tsv', sep='\t', index=False)
    effect_measures = ['mean_logexpr_state', 'detection_state', 'mean_difference', 'detection_difference',
                       'standardized_mean_difference', 'standardized_detection_difference']
    cv = effects[effects.culture_operator.eq('CV')].set_index(['cell_line', 'gene', 'state'])
    mw = effects[effects.culture_operator.eq('MW')].set_index(['cell_line', 'gene', 'state'])
    delta = cv[effect_measures]-mw[effect_measures]
    delta.columns = [x+'_CV_minus_MW' for x in delta.columns]
    for operator, frame in [('CV', cv), ('MW', mw)]:
        delta[operator+'_n_state'] = frame.n_state
        delta[operator+'_n_other'] = frame.n_other
    delta.reset_index().to_csv(out / 'tables/condition_marker_state_interactions.tsv.gz', sep='\t', index=False)
    return effects, standard, shifts


def fmt(value, percent=False):
    if pd.isna(value):
        return 'undefined'
    return f'{100*value:.1f}%' if percent else f'{value:.3f}'


def make_report(out, baseline, paired, gate_pairs, shifts, genes, provenance):
    lines = ['# DIV90 CV/MW culture-condition comparison', '',
        'CV denotes the user-supplied higher-glucose condition (~2x); MW denotes lower glucose (~1x). '
        'All technical sample IDs contain MW. Original IDs and replicate labels are preserved.', '',
        'Condition is confounded with operator and other culture differences. Each cell line contributes one sample '
        'per condition. The comparisons below are descriptive paired line contrasts; cells are not independent '
        'biological replicates and no glucose-specific causal effect is claimed.', '',
        '## Starting population', '',
        '| Sample | Line | Condition | N | SST RNA detected | PV-biased | Hybrid | SST-biased | Unresolved | Median PV score | Median SST score |',
        '|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    d = baseline[baseline.grouping.eq('culture_operator|cell_line|sample')].sort_values('sample')
    for r in d.itertuples():
        lines.append(f'| {r.sample} | {r.cell_line} | {r.culture_operator} | {r.n} | {fmt(r.sst_rna_detected_fraction,True)} | '
            f'{fmt(r.pv_fraction,True)} | {fmt(r.hybrid_fraction,True)} | {fmt(r.sst_fraction,True)} | '
            f'{fmt(r.unresolved_fraction,True)} | {fmt(r.pv_score_median)} | {fmt(r.sst_score_median)} |')
    pooled = baseline[baseline.grouping.eq('culture_operator')].sort_values('culture_operator')
    lines += ['', 'Pooled conditions have unequal cell-line weights; use the paired line contrasts below when assessing consistency.', '',
        '| Condition | N | SST RNA detected | PV-biased | Hybrid | SST-biased | Unresolved | PV median (10–90%) | SST median (10–90%) |',
        '|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in pooled.itertuples():
        lines.append(f'| {r.culture_operator} | {r.n} | {fmt(r.sst_rna_detected_fraction,True)} | '
            f'{fmt(r.pv_fraction,True)} | {fmt(r.hybrid_fraction,True)} | {fmt(r.sst_fraction,True)} | '
            f'{fmt(r.unresolved_fraction,True)} | {fmt(r.pv_score_median)} ({fmt(r.pv_score_q10)}–{fmt(r.pv_score_q90)}) | '
            f'{fmt(r.sst_score_median)} ({fmt(r.sst_score_q10)}–{fmt(r.sst_score_q90)}) |')
    lines += ['', '| Condition | FAT3 mean RNA | FAT3 detected | PTPRM mean RNA | PTPRM detected |',
              '|---|---:|---:|---:|---:|']
    for r in pooled.itertuples():
        lines.append(f'| {r.culture_operator} | {fmt(r.FAT3_mean_logexpr)} | {fmt(r.FAT3_detected_fraction,True)} | '
            f'{fmt(r.PTPRM_mean_logexpr)} | {fmt(r.PTPRM_detected_fraction,True)} |')
    lines += ['', '## Paired condition differences', '', 'Every value below is CV minus MW within a cell line. Fractions are percentage-point differences.', '',
        '| Line | Δ median PV | Δ median SST | Δ PV-biased | Δ hybrid | Δ SST-biased | Δ SST RNA detected |',
        '|---|---:|---:|---:|---:|---:|---:|']
    for line, group in paired.groupby('cell_line', sort=True):
        values = group.set_index('metric').CV_minus_MW
        keys = ['pv_score_median', 'sst_score_median', 'pv_fraction', 'hybrid_fraction', 'sst_fraction', 'sst_rna_detected_fraction']
        display = [fmt(values[k]) if i < 2 else f'{100*values[k]:+.1f}' for i,k in enumerate(keys)]
        lines.append('| '+str(line)+' | '+' | '.join(display)+' |')
    lines += ['', 'H9 and 79B have higher PV scores, lower SST scores and fewer SST-biased cells in CV; 2E shows the opposite direction. '
        'The 2E/CV sample contains only 37 entry cells. These results support line-dependent culture associations, not one uniform glucose effect. '
        'FAT3 and PTPRM mean RNA are lower in CV in all three paired lines.']
    lines += ['', '## State mixture versus marker expression within states', '',
        'The following differences use the same four operational states in both conditions. Raw mean differences '
        'include state-mixture changes. Equal-state differences assign 25% to each state in each condition. '
        'The symmetric decomposition partitions the raw difference into state-proportion and within-state-expression '
        'components and verifies they sum to the raw difference. These are descriptive standardizations, not causal '
        'mediation or a regression of culture effects. Values are undefined if any state lacks support in either condition.', '',
        '| Scope | Marker | Raw mean Δ | Equal-state mean Δ | Composition component | Within-state component |',
        '|---|---|---:|---:|---:|---:|']
    for r in shifts[shifts.gene.isin(genes)&shifts.endpoint.eq('mean_logexpr')].sort_values(['gene','cell_line']).itertuples():
        lines.append(f'| {r.cell_line} | {r.gene} | {fmt(r.raw_CV_minus_MW)} | {fmt(r.equal_state_CV_minus_MW)} | '
            f'{fmt(r.composition_component)} | {fmt(r.within_state_component)} |')
    marker_shifts = shifts[shifts.cell_line.eq('all') & shifts.endpoint.eq('mean_logexpr')].set_index('gene')
    lines += ['', 'The FAT3/PTPRM shifts persist within operational states: '
        + '; '.join(f'{g} raw CV−MW {fmt(marker_shifts.loc[g,"raw_CV_minus_MW"])} versus within-state component '
                    f'{fmt(marker_shifts.loc[g,"within_state_component"])}' for g in ['FAT3','PTPRM'])
        + '. Thus changing four-state composition alone does not explain these marker differences. '
        'The broad states still contain continuous biological and technical variation; this is not proof of direct glucose regulation.']
    lines += ['', 'The marker-state interaction table reports CV-minus-MW differences in state-versus-rest effects '
        'within each line and pooled. A second contrast gives equal weight to the other three states so changing '
        'rest-population composition is not mistaken for a change in the association. Review support counts for '
        'each contrast; the 37-cell 2E/CV sample has limited resolution.']
    effects = pd.read_csv(out / 'tables/condition_marker_state_effects.tsv.gz', sep='\t')
    focus = effects[effects.cell_line.eq('all') & effects.gene.isin(['FAT3','PTPRM','FGFR2','PTPRS'])
                    & effects.state.isin(STATES[:3])]
    lines += ['', 'Pooled state-versus-equal-weight-rest mean RNA effects; positive values mean higher marker expression in that state:', '',
        '| Marker | State | CV effect | MW effect | CV minus MW |', '|---|---|---:|---:|---:|']
    for (gene, state), d in focus.groupby(['gene','state'], sort=True):
        v = d.set_index('culture_operator').standardized_mean_difference
        lines.append(f'| {gene} | {state} | {fmt(v.CV)} | {fmt(v.MW)} | {fmt(v.CV-v.MW)} |')
    lines += ['', 'FAT3 and PTPRM associate positively with hybrid and SST-biased states and negatively with PV-biased states '
        'in both pooled conditions, with stronger hybrid contrasts in MW. Their association strength changes even after fixing '
        'the rest-state weights; their major direction is preserved. This explains why retaining low/undetected expression can '
        'raise PV-biased purity while losing many hybrid cells. PTPRS remains higher in SST-biased and lower in hybrid cells in '
        'both conditions; FGFR2 has a small positive SST-biased contrast in both, but its weak hybrid contrast changes sign. '
        'These are descriptive associations within the frozen entry population, not proof of marker specificity or stable neuronal protein staining.', '',
        '## Fixed global gates', '',
        'All gate thresholds are frozen from pooled discovery or phase one. Applying them within condition and '
        'cell line is a biological/sample comparison, not held-out validation or separate threshold optimization.', '',
        '| Line | Gate | Δ PV recovery (pp) | Δ hybrid recovery (pp) | Δ retained SST contamination (pp) | Δ total recovery (pp) |',
        '|---|---|---:|---:|---:|---:|']
    for (line, gate), group in gate_pairs[gate_pairs.gate_name.isin(['phase1_reference','selected','experimental_depletion'])].groupby(['cell_line','gate_name']):
        v = group.set_index('metric').CV_minus_MW
        display = ['undefined' if pd.isna(v[k]) else f'{100*v[k]:+.1f}' for k in ['pv_recovery','hybrid_recovery','sst_contamination','total_recovery']]
        lines.append('| '+line+' | '+gate+' | '+' | '.join(display)+' |')
    gates = pd.read_csv(out / 'tables/condition_gate_metrics.tsv', sep='\t')
    fixed = gates[gates.grouping.eq('culture_operator') & gates.gate_name.isin(['phase1_reference','experimental_depletion'])]
    lines += ['', 'Absolute pooled condition metrics for the same fixed gates:', '',
        '| Condition | Gate | Retained N | PV recovery | Hybrid recovery | Retained PV-biased | Retained hybrid | Retained SST-biased | Retained unresolved |',
        '|---|---|---:|---:|---:|---:|---:|---:|---:|']
    for r in fixed.itertuples():
        lines.append(f'| {r.culture_operator} | {r.gate_name} | {int(r.retained_n)} | {fmt(r.pv_recovery,True)} | '
            f'{fmt(r.hybrid_recovery,True)} | {fmt(r.retained_pv_fraction,True)} | {fmt(r.retained_hybrid_fraction,True)} | '
            f'{fmt(r.retained_sst_fraction,True)} | {fmt(r.retained_unresolved_fraction,True)} |')
    lines += ['', 'The old FAT3/PTPRM gate enriches the PV-biased fraction in both pooled conditions, but reduces the hybrid fraction in both. '
        'It is therefore not a joint PV/hybrid enrichment gate. Its hybrid recovery is substantially worse in MW than CV, '
        'and the paired line table shows that gate behavior is not invariant across samples. '
        'The final FGFR2/PTPRS candidate retains more than 87% of hybrid cells in each pooled condition and slightly increases both '
        'PV-biased and hybrid fractions. Its SST depletion remains modest, and the small 2E/CV group gives an unstable PV-recovery estimate. '
        'The absolute per-sample gate metrics remain available in the table files; pooled comparisons must not hide line imbalance.']
    lines += ['', '## Files', '',
        '- `condition_baseline.tsv`: counts, SST RNA prevalence, program scores, state composition, selected-marker means and detection.',
        '- `condition_paired_line_differences.tsv`: all descriptive CV-minus-MW baseline differences within the three cell lines.',
        '- `condition_gate_metrics.tsv`: both removed/retained fractions for every key gate, including each condition × line × sample.',
        '- `condition_gate_paired_line_differences.tsv`: line-paired differences for fixed gate metrics.',
        '- `condition_marker_state_effects.tsv.gz`: all allowed surface-marker state effects separately by condition and line.',
        '- `condition_marker_state_interactions.tsv.gz`: differences in marker-state effects between conditions, with state support counts.',
        '- `condition_marker_state_standardized.tsv`: raw/equal-state marker means and detection plus continuous-score associations.',
        '- `condition_marker_shift_decomposition.tsv`: state-mixture and within-state components, with reconstruction errors.',
        '- `condition_metadata_verified.tsv`: user-supplied condition metadata verified against original sample and cell-line identities.', '',
        'Selected surface genes remain separate from developmental state definitions. All expression thresholds and means '
        'are RNA log1p(CP10K), not measured cell-surface protein or fluorescence.']
    (out / 'CONDITION_REPORT.md').write_text('\n'.join(lines)+'\n')


def run(cells: pd.DataFrame, out: Path | str, summary: dict):
    out = Path(out).resolve()
    if out != OUT.resolve() or REPO_ROOT.resolve() in out.parents:
        raise ValueError(f'Condition runtime must be {OUT}')
    (out / 'tables').mkdir(exist_ok=True)
    (out / 'provenance').mkdir(exist_ok=True)
    cells = cells.reset_index(drop=True)
    provenance = verify_metadata(cells, out)
    selected = set(REFERENCE_GENES)
    for name in GATE_NAMES:
        if isinstance(summary.get(name), dict):
            selected.update(x['gene'] for x in summary['rules'][summary[name]['gate_id']]['rules'])
    selected = sorted(g for g in selected if g in cells)
    annotations = pd.read_csv(PROJECT_ROOT / 'results/div90_hypergate_sst_pv/allowed_surface_markers.tsv', sep='\t')
    excluded = set(summary.get('excluded_module_features', []))
    allgenes = [g for g in annotations.gene if g in cells and g not in excluded]
    baseline = summarize_baseline(cells, selected, out)
    gates = summarize_gates(cells, summary, out)
    paired = paired_table(baseline, ['cell_line'], [], out / 'tables/condition_paired_line_differences.tsv')
    gate_pairs = paired_table(gates, ['cell_line','gate_name'], [], out / 'tables/condition_gate_paired_line_differences.tsv')
    effects, standard, shifts = marker_analysis(cells, allgenes, out)
    provenance.update(selected_surface_genes=selected, analyzed_surface_genes=len(allgenes),
        equal_state_weights={s: .25 for s in STATES},
        missing_state_policy='Standardization/decomposition undefined when any state lacks cells in a condition.',
        decomposition='CV-MW = sum_s[(pCV-pMW)*(mCV+mMW)/2] + sum_s[(mCV-mMW)*(pCV+pMW)/2]',
        max_decomposition_error=float(shifts.reconstruction_error.abs().max()),
        no_cell_level_p_values=True, no_condition_regression=True,
        gate_thresholds='Frozen global gates; no per-condition reoptimization.')
    write_json(out / 'provenance/condition_analysis.json', provenance)
    make_report(out, baseline, paired, gate_pairs, shifts, selected, provenance)
    result = dict(condition_analysis_complete=True, selected_surface_genes=selected,
                  analyzed_surface_genes=len(allgenes), sample_n=6, paired_cell_lines=['H9','79B','2E'],
                  report=str(out / 'CONDITION_REPORT.md'))
    write_json(out / 'condition_summary.json', result)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, default=OUT)
    args = parser.parse_args()
    run(pd.read_csv(args.out / 'cells.tsv.gz', sep='\t'), args.out,
        json.loads((args.out / 'gate_summary.json').read_text()))
