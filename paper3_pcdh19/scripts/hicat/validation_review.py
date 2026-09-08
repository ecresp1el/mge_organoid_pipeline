"""Conservative, reproducible review of the existing E14.5 mouse MGE hierarchy.

This module does not fit clusters, change labels, or remove cells. User-supplied
identities are review hypotheses. Independent canonical programs, representation
and stability metrics determine operational review categories, not annotations.
Every threshold is saved with the results; descriptive association is not proof
that a state or technical variable caused a boundary.
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse

DEVELOPMENT = ['rg_stemness', 'ipc_neurogenic', 'neuroblast', 'immature_inhibitory',
               'early_neuronal_maturation', 'later_neuronal_maturation', 'endothelial', 'pericyte', 'erythroid']
REGION = ['MGE', 'MGE_interneuron', 'LGE', 'CGE', 'POA', 'basal_forebrain_alternative']
EXPECTED = {'RG': 'rg_stemness', 'IPC': 'ipc_neurogenic', 'neuroblast': 'neuroblast',
            'immature inhibitory': 'immature_inhibitory', 'maturation': 'later_neuronal_maturation',
            'endothelial': 'endothelial', 'erythroid': 'erythroid'}
DEFAULT_THRESHOLDS = dict(program_mean_positive=0.10, marker_detection_fraction=0.10,
    canonical_minimum_markers=2, strong_hypothesis_minimum_markers=3,
    regional_anchor_minimum_markers=2, MGE_generic_genes_excluded_from_regional_support=['Gad1','Gad2','Arx'],
    cycle_fraction=0.70, cycle_phase_tv=0.40, maturation_z=0.70,
    technical_robust_z=3.0, minimum_review_cells=30,
    readiness_minimum_nonlow_seed_fraction=0.25, readiness_median_seed_jaccard=0.50,
    readiness_minimum_surviving_allen_fraction=0.25, readiness_major_parent_fraction=0.05,
    maturation_max_regional_effect=0.30, maturation_max_phase_tv=0.20)


def _number(value, default=np.nan):
    """Convert numeric table values while preserving unavailable evidence as NaN."""
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _get(row, names, default=np.nan):
    """Read supported metric aliases without interpreting absent values as zero."""
    for name in names:
        if name in row:
            return row[name]
    return default


def _indexed(frame):
    """Normalize an explicit cluster column to an index for metric joins."""
    frame = frame.copy()
    for key in ['cluster', 'baseline_cluster', 'current_cluster', 'fine_cluster']:
        if key in frame.columns:
            return frame.set_index(key)
    return frame


def _top_program(row, names, threshold):
    """Name a positive leading program, or retain the unresolved state."""
    values = row.reindex(names).dropna()
    if values.empty or values.max() < threshold:
        return 'unresolved'
    return str(values.idxmax())


def expression_summaries(pilot):
    """Return coarse/fine all-assayed-gene means and detection fractions.

    Duplicate symbols are aggregated by maximum mean/detection for display only;
    gene IDs remain available in the saved pilot and score provenance. No DE
    genes are used to select the canonical programs or fit any model.
    """
    matrix = sparse.csr_matrix(pilot.layers['log1p_cpm'])
    symbols = pilot.var['gene_symbol'].astype(str).to_numpy()
    means, fractions = {}, {}
    for level in ['coarse', 'fine']:
        labels = pilot.obs['hicat_' + level + '_baseline'].astype(str).to_numpy()
        groups = sorted(set(labels))
        mean = pd.DataFrame([np.asarray(matrix[labels == c].mean(axis=0)).ravel() for c in groups], index=groups, columns=symbols)
        frac = pd.DataFrame([np.asarray((matrix[labels == c] > 0).mean(axis=0)).ravel() for c in groups], index=groups, columns=symbols)
        means[level] = mean.T.groupby(level=0, sort=False).max().T
        fractions[level] = frac.T.groupby(level=0, sort=False).max().T
    return means, fractions


def _phase_counts(pilot, level):
    """Calculate observed cycle composition, retaining all standard phases."""
    phase = pilot.obs['validation_phase'].astype(str).replace({'G1': 'G1-like', 'G2/M': 'G2M'})
    counts = pd.crosstab(pilot.obs['hicat_' + level + '_baseline'].astype(str), phase)
    return counts.reindex(columns=['G1-like', 'S', 'G2M'], fill_value=0)


def _pair_concerns(pair_table, cluster, cycle_tv_threshold):
    """Identify phase-associated sibling comparisons without asserting causality."""
    if pair_table is None or pair_table.empty:
        return 0, 0, 'not_available'
    a = next((v for v in ['cluster_a', 'cluster1', 'cluster_1', 'left'] if v in pair_table), None)
    b = next((v for v in ['cluster_b', 'cluster2', 'cluster_2', 'right'] if v in pair_table), None)
    if not a or not b:
        return 0, 0, 'not_available'
    selected = pair_table[(pair_table[a] == cluster) | (pair_table[b] == cluster)]
    status_col = next((v for v in ['identity_difference_persists_after_cell_cycle_stratification', 'identity_persists', 'persists'] if v in selected), None)
    tv_col = next((v for v in ['cycle_composition_tv', 'phase_total_variation', 'phase_composition_tv', 'cycle_tv', 'phase_tv_distance'] if v in selected), None)
    if not status_col:
        return len(selected), 0, 'not_available'
    status = selected[status_col].astype(str)
    associated = (pd.to_numeric(selected[tv_col], errors='coerce') >= cycle_tv_threshold) if tv_col else pd.Series(False, index=selected.index)
    unresolved = associated & (status == 'no')
    counts = status.value_counts()
    description = '; '.join('%s=%d' % (key, value) for key, value in counts.items())
    return len(selected), int(unresolved.sum()), description


def maturation_boundaries(pilot, scores, pairs, thresholds):
    """Describe sibling maturation contrasts separately from region and phase.

    Effects use full-pilot score standard deviations. A candidate maturation
    boundary requires a large maturation-axis contrast, small regional contrast,
    similar phase composition, and at least one adequately populated matched
    phase. This is association, not a fitted trajectory or causal explanation.
    """
    labels=pilot.obs.hicat_fine_baseline.astype(str)
    axis=scores['later_neuronal_maturation']-scores['rg_stemness']
    means=scores.groupby(labels,observed=True).mean()
    axis_mean=axis.groupby(labels,observed=True).mean()
    axis_sd=float(axis.std(ddof=0))
    regional_sd=scores.reindex(columns=REGION).std(ddof=0).replace(0,np.nan)
    fractions=_phase_counts(pilot,'fine')
    fractions=fractions.div(fractions.sum(axis=1),axis=0)
    records=[]
    for parent in sorted(set(c.split('.')[0] for c in means.index)):
        groups=[c for c in means.index if c.split('.')[0]==parent]
        for i,a in enumerate(groups):
            for b in groups[i+1:]:
                effect=abs(float(axis_mean[a]-axis_mean[b]))/axis_sd if axis_sd else np.nan
                regional=(means.loc[a,REGION]-means.loc[b,REGION]).abs().div(regional_sd)
                region_max=float(regional.max())
                tv=float((fractions.loc[a]-fractions.loc[b]).abs().sum()/2)
                eligible=0
                if pairs is not None and not pairs.empty and 'cluster_a' in pairs and 'cluster_b' in pairs:
                    match=pairs[((pairs.cluster_a==a)&(pairs.cluster_b==b))|((pairs.cluster_a==b)&(pairs.cluster_b==a))]
                    if len(match) and 'n_eligible_phases' in match:eligible=int(match.iloc[0].n_eligible_phases)
                flag=bool(effect>=thresholds['maturation_z'] and region_max<=thresholds['maturation_max_regional_effect'] and tv<=thresholds['maturation_max_phase_tv'] and eligible>0)
                records.append(dict(parent=parent,cluster_a=a,cluster_b=b,maturation_axis_standardized_difference=effect,
                    maximum_regional_standardized_difference=region_max,cell_cycle_total_variation=tv,n_eligible_phases=eligible,
                    maturation_boundary_review=flag,interpretation='Descriptive maturation-associated boundary candidate; no trajectory, regression or causal claim'))
    return pd.DataFrame(records)


def review_clusters(pilot, score_frame, metrics, hypotheses_path, run_root, config,
                    program_metadata=None):
    """Build cluster, parent and hypothesis decision tables from saved evidence.

    Parameters
    ----------
    pilot : anndata.AnnData
        Existing 12,000-cell object with unchanged baseline/repeat memberships,
        saved natural log1p(CPM), and newly computed validation_phase.
    score_frame : pandas.DataFrame
        Cell-indexed canonical scores; mean expression minus fixed controls.
    metrics : mapping
        ``seed_coarse``, ``seed_fine``, ``allen_fine``, ``sample_coarse``,
        ``sample_fine`` and ``phase_pairs`` dictionaries from validation_metrics.
    hypotheses_path : path
        Review-only user reference. It is read after all clustering/metrics.
    run_root : path
        Working run directory; outputs are written beneath outputs/.
    config : mapping
        ``programs_path`` plus optional ``review_thresholds`` overrides.
    program_metadata : mapping, optional
        Score model/coverage metadata retained by the scorer.

    Returns
    -------
    dict
        DataFrames fine/coarse/parent/hypotheses, source expression/program
        summaries, and a serializable readiness summary. No adopted labels.
    """
    out = Path(run_root) / 'outputs'
    review_dir = out / 'annotation_review'
    review_dir.mkdir(parents=True, exist_ok=True)
    thresholds = dict(DEFAULT_THRESHOLDS)
    thresholds.update(config.get('review_thresholds', {}))
    (review_dir / 'review_thresholds.json').write_text(json.dumps(thresholds, indent=2) + '\n')
    reference = json.loads(Path(hypotheses_path).read_text())
    hypotheses = {row['cluster']: row for row in reference['clusters']}
    program_cfg = json.loads(Path(config['programs_path']).read_text())
    programs = {p['name']: p for p in program_cfg['programs']}
    gene_means, gene_fractions = expression_summaries(pilot)
    all_review, program_means, phase_fractions = {}, {}, {}
    pair_table = metrics.get('phase_pairs', {}).get('pairs', pd.DataFrame())
    maturation_pairs=maturation_boundaries(pilot,score_frame,pair_table,thresholds)
    maturation_pairs.to_csv(review_dir/'maturation_boundary_review.tsv',sep='\t',index=False)
    qc_columns = [key for key in ['total_counts', 'n_genes_by_counts', 'pct_counts_mt'] if key in pilot.obs]
    # Robust z-scores are descriptive QC contrasts across cells, never filters.
    qc = pilot.obs[qc_columns].astype(float).copy()
    for key in ['total_counts', 'n_genes_by_counts']:
        if key in qc:
            qc[key] = np.log1p(qc[key])
    qc_center, qc_mad = qc.median(), (qc - qc.median()).abs().median() * 1.4826
    qc_z = (qc - qc_center) / qc_mad.replace(0, np.nan)
    for level in ['coarse', 'fine']:
        labels = pilot.obs['hicat_' + level + '_baseline'].astype(str)
        program_means[level] = score_frame.groupby(labels, observed=True).mean()
        phase_count = _phase_counts(pilot, level)
        phase_fractions[level] = phase_count.div(phase_count.sum(axis=1), axis=0)
        sample = _indexed(metrics['sample_' + level]['clusters'])
        seed = _indexed(metrics['seed_' + level]['clusters'])
        allen = _indexed(metrics['allen_fine']['clusters']) if level == 'fine' else pd.DataFrame()
        standard = program_means[level].sub(score_frame.mean()).div(score_frame.std(ddof=0).replace(0, np.nan))
        records = []
        for cluster, n_cells in labels.value_counts().sort_index().items():
            h = hypotheses.get(cluster, {})
            row = program_means[level].loc[cluster]
            detect = gene_fractions[level].loc[cluster]
            sm = sample.loc[cluster]
            sd = seed.loc[cluster]
            al = allen.loc[cluster] if cluster in allen.index else pd.Series(dtype=object)
            phase = phase_fractions[level].loc[cluster]
            developmental = _top_program(row, DEVELOPMENT, thresholds['program_mean_positive'])
            regional = _top_program(row, REGION, thresholds['program_mean_positive'])
            expected_dev = EXPECTED.get(h.get('developmental_hypothesis', ''))
            expected_region = h.get('regional_hypothesis', 'unresolved')
            # Check user hypotheses against predeclared canonical signatures,
            # not against the user's DEG-derived supporting-gene list.
            expected_programs = [expected_dev] if expected_dev else []
            if expected_region in REGION:
                expected_programs.append(expected_region)
            if expected_region == 'MGE':
                expected_programs = list(dict.fromkeys(expected_programs + ['MGE', 'MGE_interneuron']))
            independent_genes = list(dict.fromkeys(g for name in expected_programs for g in programs.get(name, {}).get('genes', [])))
            supporting = [g for g in independent_genes if g in detect and detect[g] >= thresholds['marker_detection_fraction']]
            positive_programs = [name for name in expected_programs if _number(row.get(name)) >= thresholds['program_mean_positive']]
            dev_supported = expected_dev in positive_programs
            regional_signature = (programs.get(expected_region, {}).get('genes', []) +
                (programs.get('MGE_interneuron', {}).get('genes', []) if expected_region == 'MGE' else []))
            regional_anchor_signature=[g for g in regional_signature if expected_region!='MGE' or g not in thresholds['MGE_generic_genes_excluded_from_regional_support']]
            regional_detected = sorted(g for g in set(regional_anchor_signature) if g in detect and detect[g] >= thresholds['marker_detection_fraction'])
            region_supported = (expected_region not in REGION or
                ((expected_region in positive_programs or (expected_region == 'MGE' and 'MGE_interneuron' in positive_programs)) and
                 len(regional_detected) >= thresholds['regional_anchor_minimum_markers']))
            conflicts = []
            if expected_dev and _number(row.get(expected_dev)) < -thresholds['program_mean_positive']:
                conflicts.append('Expected %s mean score is negative (%.3f)' % (expected_dev, row[expected_dev]))
            other_region = 'LGE' if expected_region == 'MGE' else 'MGE_interneuron' if expected_region == 'LGE' else None
            counter_genes = []
            if other_region and _number(row.get(other_region)) >= thresholds['program_mean_positive']:
                counter_genes = [g for g in programs[other_region]['genes'] if g in detect and detect[g] >= thresholds['marker_detection_fraction']]
                counter_anchors=[g for g in counter_genes if other_region!='MGE_interneuron' or g not in thresholds['MGE_generic_genes_excluded_from_regional_support']]
                if len(counter_anchors) >= thresholds['regional_anchor_minimum_markers']:
                    conflicts.append('Alternative %s program is also positive; regional mixture/shared expression remains possible' % other_region)
            if expected_region == 'POA':
                conflicts.append('Nkx2-1 is shared across MGE/POA; POA hypothesis requires multiple independent markers')
            cycle_frac = float(phase.get('S', 0) + phase.get('G2M', 0))
            if h.get('cell_cycle_hypothesis') == 'postmitotic' and cycle_frac >= thresholds['cycle_fraction']:
                conflicts.append('Postmitotic hypothesis conflicts with %.1f%% assigned S/G2M cells' % (100 * cycle_frac))
            if expected_dev is None:
                hypothesis_support = 'PRESENTLY_UNRESOLVABLE'
            elif dev_supported and region_supported and len(supporting) >= thresholds['strong_hypothesis_minimum_markers'] and not conflicts:
                hypothesis_support = 'STRONGLY_SUPPORTED'
            elif dev_supported or len(positive_programs):
                hypothesis_support = 'WEAKLY_SUPPORTED'
            elif _number(row.get(expected_dev)) < -thresholds['program_mean_positive'] and developmental != 'unresolved':
                hypothesis_support = 'CONTRADICTED'
            else:
                hypothesis_support = 'PRESENTLY_UNRESOLVABLE'
            top_fraction = _number(_get(sm, ['top_sample_fraction', 'maximum_sample_fraction', 'max_sample_fraction']))
            entropy = _number(_get(sm, ['sample_entropy', 'entropy']))
            effective = _number(_get(sm, ['effective_samples', 'effective_number_of_samples', 'effective_sample_number']))
            sample_flag = _get(sm, ['sample_concern', 'sample_flag', 'flag_sample_associated'], False)
            if isinstance(sample_flag, str):
                sample_flag = sample_flag.lower() not in ['false', 'none', 'no', 'none_flagged', 'not_flagged', '']
            sample_flag = bool(sample_flag) if pd.notna(sample_flag) else False
            pair_n, cycle_pair_n, pair_summary = _pair_concerns(pair_table, cluster, thresholds['cycle_phase_tv']) if level == 'fine' else (0, 0, 'see fine sibling comparisons')
            cycle_concern = cycle_pair_n > 0
            stress_z = _number(standard.loc[cluster].get('stress'))
            local_qc_z = qc_z.loc[labels == cluster].median()
            qc_flags = []
            for key in ['total_counts', 'n_genes_by_counts']:
                if key in local_qc_z and local_qc_z[key] <= -thresholds['technical_robust_z']:
                    qc_flags.append(key + '_low')
            if 'pct_counts_mt' in local_qc_z and local_qc_z['pct_counts_mt'] >= thresholds['technical_robust_z']:
                qc_flags.append('pct_counts_mt_high')
            if stress_z >= thresholds['technical_robust_z'] and _number(row.get('stress')) > thresholds['program_mean_positive']:
                qc_flags.append('stress_program_high')
            technical_concern = bool(qc_flags)
            maturity_selected=maturation_pairs[(maturation_pairs.cluster_a==cluster)|(maturation_pairs.cluster_b==cluster)] if level=='fine' and len(maturation_pairs) else pd.DataFrame()
            maturation_concern=bool(maturity_selected.maturation_boundary_review.any()) if len(maturity_selected) else False
            seed_stability = str(_get(sd, ['stability', 'seed_stability'], 'unavailable'))
            allen_survives = str(_get(al, ['boundary_survives', 'Allen_boundary_survival'], 'not_applicable_fixed_parent'))
            stable = seed_stability == 'HIGH_STABILITY'
            allen_strong = allen_survives == 'strong'
            confidence = ('Moderate' if hypothesis_support == 'STRONGLY_SUPPORTED' else
                          'Low' if hypothesis_support == 'WEAKLY_SUPPORTED' else 'Insufficient')
            if n_cells < thresholds['minimum_review_cells']:
                confidence += '; small cluster'
            # Labels are operational and conservative: cycling alone does not
            # imply that cycling caused the cluster boundary.
            if technical_concern:
                status, action = 'TECHNICAL_OR_STRESS_ASSOCIATED', 'retain_but_review'
            elif sample_flag:
                status, action = 'SAMPLE_ASSOCIATED', 'sample_effect_review'
            elif cycle_concern:
                status, action = 'CELL_CYCLE_DOMINATED', 'cell_cycle_boundary_review'
            elif maturation_concern:
                status, action = 'MATURATION_DOMINATED', 'retain_but_review'
            elif hypothesis_support in ['PRESENTLY_UNRESOLVABLE', 'CONTRADICTED']:
                status, action = 'INSUFFICIENT_EVIDENCE', 'insufficient_evidence'
            elif stable and (level == 'coarse' or allen_strong) and hypothesis_support == 'STRONGLY_SUPPORTED':
                status, action = 'SUPPORTED', 'retain_candidate_cluster'
            else:
                status, action = 'PLAUSIBLE_BUT_UNSTABLE', 'retain_but_review'
            record = dict(cluster=cluster, fine_cluster=cluster if level == 'fine' else '', coarse_parent=cluster.split('.')[0],
                level=level, n_cells=int(n_cells), fraction_of_pilot=n_cells/pilot.n_obs,
                provisional_identity=h.get('provisional_identity', 'No supplied hypothesis'), identity_confidence=confidence,
                user_confidence=h.get('user_confidence', 'not supplied'), hypothesis_support=hypothesis_support,
                dominant_developmental_program=developmental, dominant_regional_program=regional,
                RG_score=row.get('rg_stemness'), IPC_score=row.get('ipc_neurogenic'), neuroblast_score=row.get('neuroblast'),
                MGE_score=row.get('MGE'), LGE_score=row.get('LGE'), maturation_score=row.get('later_neuronal_maturation'),
                S_score=row.get('S_phase'), G2M_score=row.get('G2M'), cell_cycle_class=str(phase.idxmax()),
                S_fraction=phase.get('S', 0), G2M_fraction=phase.get('G2M', 0), G1_like_fraction=phase.get('G1-like', 0),
                cell_cycle_concern=cycle_concern, maturation_concern=maturation_concern, cycle_associated_unresolved_sibling_pairs=cycle_pair_n,
                phase_stratified_pair_summary=pair_summary, sibling_pairs=pair_n,
                top_sample=str(_get(sm, ['top_sample', 'largest_sample'], 'unavailable')), top_sample_fraction=top_fraction,
                sample_entropy=entropy, effective_samples=effective, sample_concern=sample_flag,
                seed_best_match=str(_get(sd, ['best_match'], 'unavailable')), seed_jaccard=_number(_get(sd, ['jaccard'])),
                seed_stability=seed_stability, Allen_best_match=str(_get(al, ['best_match'], 'not_applicable')),
                Allen_jaccard=_number(_get(al, ['jaccard'])), Allen_boundary_survival=allen_survives,
                Allen_merges_with=str(_get(al, ['merges_with'], '[]')), Allen_splits_into=str(_get(al, ['splits_into', 'split_targets'], '[]')),
                stress_score=row.get('stress'), technical_concern=technical_concern, technical_flags='; '.join(qc_flags) or 'none_at_review_threshold',
                canonical_supporting_markers='; '.join(supporting), supporting_programs='; '.join(positive_programs),
                regional_anchor_markers='; '.join(regional_detected), regional_anchor_count=len(regional_detected),
                conflicting_markers='; '.join(counter_genes), conflicting_evidence='; '.join(conflicts) or 'No contradiction at these descriptive thresholds; absence is not proof',
                alternative_interpretation=h.get('alternative_interpretation', 'Unresolved'),
                additional_markers='; '.join(h.get('additional_markers', [])),
                overall_validation_status=status, recommended_action=action)
            records.append(record)
        all_review[level] = pd.DataFrame(records).set_index('cluster', drop=False)
        program_means[level].to_csv(out/'marker_programs'/('%s_cluster_program_means.tsv' % level), sep='\t', index_label='cluster')
        score_frame.groupby(labels, observed=True).median().to_csv(out/'marker_programs'/('%s_cluster_program_medians.tsv' % level), sep='\t', index_label='cluster')
        phase_count.to_csv(out/'cell_cycle'/('%s_phase_counts.tsv' % level), sep='\t', index_label='cluster')
        phase_fractions[level].to_csv(out/'cell_cycle'/('%s_phase_fractions.tsv' % level), sep='\t', index_label='cluster')
    fine, coarse = all_review['fine'], all_review['coarse']
    parent_records = []
    for cluster, parent in coarse.iterrows():
        children = fine[fine.coarse_parent == cluster]
        mask = pilot.obs.hicat_coarse_baseline.astype(str) == cluster
        repeat_coarse = parent.seed_best_match
        repeat_mask = pilot.obs.hicat_coarse_seed_repeat.astype(str) == repeat_coarse
        n_allen = len(set(pilot.obs.loc[mask, 'hicat_fine_allen_reference'].astype(str))) if 'hicat_fine_allen_reference' in pilot.obs else np.nan
        n_repeat = pilot.obs.loc[repeat_mask, 'hicat_fine_seed_repeat'].astype(str).nunique()
        profiles=children[['S_fraction', 'G2M_fraction', 'G1_like_fraction']].to_numpy(dtype=float)
        cycle_spread = float(max((np.abs(a-b).sum()/2 for a in profiles for b in profiles), default=0))
        cycling_parent = (parent.S_fraction + parent.G2M_fraction) >= thresholds['cycle_fraction']
        if parent.sample_concern:
            coherence = 'sample-associated'
        elif cycling_parent and cycle_spread >= thresholds['cycle_phase_tv']:
            coherence = 'primarily cell-cycle-defined'
        elif children.maturation_concern.any() and not children.cell_cycle_concern.any():
            coherence = 'primarily maturation-defined'
        elif parent.hypothesis_support == 'PRESENTLY_UNRESOLVABLE':
            coherence = 'poorly resolved'
        elif children.dominant_developmental_program.nunique() > 2:
            coherence = 'mixed but interpretable'
        elif parent.hypothesis_support == 'STRONGLY_SUPPORTED':
            coherence = 'biologically coherent'
        else:
            coherence = 'mixed but interpretable'
        parent_records.append(dict(parent=cluster, n_cells=parent.n_cells, n_fine_current=len(children),
            n_fine_repeat=n_repeat, repeat_parent_best_match=repeat_coarse,
            repeat_count_definition='children of best-matching repeat coarse parent; not assumed identical parent',
            n_fine_repeat_intersecting_current_parent=pilot.obs.loc[mask,'hicat_fine_seed_repeat'].astype(str).nunique(),
            n_fine_Allen_reference=n_allen, dominant_identity=parent.dominant_regional_program,
            dominant_developmental_state=parent.dominant_developmental_program,
            cell_cycle_structure='G1-like %.1f%%; S %.1f%%; G2M %.1f%%' % tuple(100 * parent[k] for k in ['G1_like_fraction','S_fraction','G2M_fraction']),
            sample_structure='effective samples %.2f; largest sample %.1f%%' % (parent.effective_samples,100*parent.top_sample_fraction),
            fine_boundary_stability='%d/%d HIGH seed; %d/%d strong Allen' % ((children.seed_stability=='HIGH_STABILITY').sum(),len(children),(children.Allen_boundary_survival=='strong').sum(),len(children)),
            cross_parent_DE_issue='See hierarchy_relationship_flags.tsv and original final-pair audits; no cross-parent merging',
            strongest_canonical_markers=parent.canonical_supporting_markers,
            strongest_conflicting_markers=parent.conflicting_markers, provisional_interpretation=parent.provisional_identity,
            confidence=parent.identity_confidence, major_concern=parent.overall_validation_status,
            recommended_status=coherence))
    parents = pd.DataFrame(parent_records)
    fine.to_csv(review_dir/'cluster_validation_summary.tsv',sep='\t',index=False)
    coarse.to_csv(review_dir/'coarse_cluster_validation_summary.tsv',sep='\t',index=False)
    parents.to_csv(review_dir/'parent_validation_summary.tsv',sep='\t',index=False)
    hypothesis_table = pd.concat([coarse, fine], ignore_index=True)
    hypothesis_table.to_csv(review_dir/'provisional_hypothesis_review.tsv',sep='\t',index=False)
    high = int((fine.seed_stability == 'HIGH_STABILITY').sum())
    strong = int((fine.Allen_boundary_survival == 'strong').sum())
    nonlow = int(fine.seed_stability.isin(['HIGH_STABILITY','MODERATE_STABILITY']).sum())
    surviving = int(fine.Allen_boundary_survival.isin(['strong','partial']).sum())
    median_jaccard = float(fine.seed_jaccard.median())
    readiness_reasons = []
    if nonlow / len(fine) < thresholds['readiness_minimum_nonlow_seed_fraction'] and median_jaccard < thresholds['readiness_median_seed_jaccard']:
        readiness_reasons.append('Only %d/%d fine clusters have HIGH/MODERATE second-seed stability; median best-match Jaccard %.3f' % (nonlow,len(fine),median_jaccard))
    if surviving / len(fine) < thresholds['readiness_minimum_surviving_allen_fraction']:
        readiness_reasons.append('Only %d/%d fine boundaries have strong/partial survival in the controlled Allen condition' % (surviving,len(fine)))
    major_unresolved = parents[(parents.recommended_status == 'poorly resolved') & (parents.n_cells / pilot.n_obs >= thresholds['readiness_major_parent_fraction'])]
    if len(major_unresolved):
        readiness_reasons.append('Major coarse parent(s) remain poorly resolved: '+', '.join(major_unresolved.parent))
    readiness = 'NOT_READY' if readiness_reasons else 'READY_WITH_CAVEATS'
    summary = dict(readiness=readiness, readiness_reasons=readiness_reasons,
        fine_clusters=len(fine), coarse_clusters=len(coarse), high_seed_clusters=high,
        strong_allen_clusters=strong, nonlow_seed_clusters=nonlow, surviving_allen_clusters=surviving, median_seed_jaccard=median_jaccard, fine_status_counts=fine.overall_validation_status.value_counts().to_dict(),
        coarse_status_counts=coarse.overall_validation_status.value_counts().to_dict(),
        hypothesis_support_counts=hypothesis_table.hypothesis_support.value_counts().to_dict(),
        annotations_locked=False, full_data_clustered=False, review_thresholds=thresholds,
        interpretation='Operational review categories; no automatic biological labels or cluster merges. Readiness refers to scaling this hierarchy, not feasibility of a future revised approach.')
    (review_dir/'review_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return dict(fine=fine,coarse=coarse,parent=parents,hypotheses=hypothesis_table,summary=summary,
                program_means=program_means,gene_means=gene_means,gene_fractions=gene_fractions,
                phase_fractions=phase_fractions)
