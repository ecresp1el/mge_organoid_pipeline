"""Label-blind diagnostics for the E14.5 dissected mouse MGE HiCAT pilot.

This module compares already fitted partitions; it never fits, renames, removes,
merges, annotates, or transfers clusters. Inputs are pandas objects indexed by
unique cell IDs. All secondary inputs are explicitly aligned to those IDs, and
missing/extra IDs raise an error. Results are ordinary tables plus JSON-safe
summaries so callers can save every plotted number and the exact decision rules.

Public entry points
-------------------
``partition_overlap`` compares two partitions of exactly the same cells.
``sample_composition`` describes each cluster's representation across samples.
``phase_matched_identity`` asks whether canonical identity-score differences
remain within observed cell-cycle strata. This is descriptive stratification,
not regression, formal DE, causal adjustment, or proof of a biological identity.
"""
from itertools import combinations
import json

import numpy as np
import pandas as pd
from scipy.stats import chisquare
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


def _aligned(reference, other, name):
    """Check unique, identical cell-ID sets and return an aligned copy of other."""
    if not isinstance(reference, (pd.Series, pd.DataFrame)):
        raise TypeError('Reference must be a pandas object indexed by cell ID.')
    if not isinstance(other, (pd.Series, pd.DataFrame)):
        raise TypeError('%s must be a pandas object indexed by cell ID.' % name)
    if not reference.index.is_unique or not other.index.is_unique:
        raise ValueError('Cell IDs must be unique in every input.')
    if reference.empty or reference.index.hasnans or other.index.hasnans:
        raise ValueError('Cell IDs must be nonempty and nonmissing.')
    missing = reference.index.difference(other.index)
    extra = other.index.difference(reference.index)
    if len(missing) or len(extra):
        raise ValueError('%s cell IDs differ: %d missing, %d extra.' %
                         (name, len(missing), len(extra)))
    return other.reindex(reference.index).copy()


def _labels(series, name):
    """Validate a nonmissing label Series without silently dropping any cells."""
    if not isinstance(series, pd.Series):
        raise TypeError('%s must be a pandas Series.' % name)
    _aligned(series, series, name)
    if series.isna().any():
        raise ValueError('%s contains missing labels.' % name)
    # Preserve ID strings in files rather than categorical positions or codes.
    return series.astype(str)


def _fraction(value, name):
    """Reject a malformed proportion threshold before producing a decision."""
    if not np.isfinite(value) or not 0 <= value <= 1:
        raise ValueError('%s must be between zero and one.' % name)


def _list_json(values):
    """Encode an ordered label list unambiguously inside a TSV cell."""
    return json.dumps(list(values), ensure_ascii=False)


def partition_overlap(baseline, comparison, *, high_jaccard=.75,
                      high_precision=.8, high_recall=.8,
                      moderate_jaccard=.5, split_fraction=.1,
                      partial_jaccard=.35, partial_recall=.5):
    """Compare two labelings of exactly the same cells, without matching IDs.

    Parameters
    ----------
    baseline, comparison : pandas.Series
        Cluster IDs indexed by unique cell IDs. Comparison order may differ;
        its cell-ID set must be identical. No cells are intersected or dropped.
    high_jaccard, high_precision, high_recall : float
        All three inclusive thresholds must pass for HIGH_STABILITY and strong
        Allen boundary survival. Defaults are 0.75, 0.8 and 0.8 respectively.
    moderate_jaccard : float
        Minimum Jaccard for MODERATE_STABILITY when the high rule fails.
    split_fraction : float
        Strictly greater than this fraction marks a substantial overlap.
        Row fractions identify split targets; column fractions identify merging
        baseline contributors. The default implements the requested >10% rule.
    partial_jaccard, partial_recall : float
        Inclusive thresholds for partial Allen boundary survival when strong
        fails. These describe cluster membership retention, not an independent
        pairwise DE test or adoption of the alternative clustering.

    Returns
    -------
    dict
        ``counts``, ``row_fractions`` and ``column_fractions`` have baseline
        clusters as rows and comparison clusters as columns. ``clusters`` has
        one row per baseline cluster, indexed by ``cluster``. ``summary`` holds
        ARI, arithmetic-mean NMI, sample size, thresholds, and display ordering.
        Best match maximizes Jaccard, then overlap, then the literal target ID.
        Precision = intersection / comparison size; recall = intersection /
        baseline size. ``fraction_retained`` is recall for that best match.
        ``max_overlap_fraction`` separately records the largest row fraction.

    Notes
    -----
    ``merges_with`` lists other baseline clusters each contributing >10% of the
    best comparison match, only if the focal baseline also contributes >10%.
    ``merge_targets`` examines all comparison targets under that same column
    rule. JSON strings preserve label lists inside TSV columns. Actual labels
    and coarse-parent blocks are retained in all suggested plot orders.
    """
    for name, value in [('high_jaccard', high_jaccard),
                        ('high_precision', high_precision),
                        ('high_recall', high_recall),
                        ('moderate_jaccard', moderate_jaccard),
                        ('split_fraction', split_fraction),
                        ('partial_jaccard', partial_jaccard),
                        ('partial_recall', partial_recall)]:
        _fraction(value, name)
    left = _labels(baseline, 'baseline')
    right = _labels(_aligned(left, comparison, 'comparison'), 'comparison')
    counts = pd.crosstab(left.rename('cluster'), right.rename('comparison_cluster'))
    counts = counts.reindex(index=sorted(counts.index), columns=sorted(counts.columns))
    row_n, col_n = counts.sum(axis=1), counts.sum(axis=0)
    row_fraction = counts.div(row_n, axis=0)
    col_fraction = counts.div(col_n, axis=1)
    rows = []
    for cluster in counts.index:
        overlap = counts.loc[cluster]
        jaccard = overlap / (row_n.loc[cluster] + col_n - overlap)
        best = sorted(counts.columns,
                      key=lambda target: (-jaccard[target], -overlap[target], target))[0]
        shared = int(overlap[best])
        precision = float(col_fraction.loc[cluster, best])
        recall = float(row_fraction.loc[cluster, best])
        high = (jaccard[best] >= high_jaccard and precision >= high_precision
                and recall >= high_recall)
        stability = ('HIGH_STABILITY' if high else 'MODERATE_STABILITY'
                     if jaccard[best] >= moderate_jaccard else 'LOW_STABILITY')
        survives = ('strong' if high else 'partial'
                    if jaccard[best] >= partial_jaccard and recall >= partial_recall
                    else 'no')
        split_targets = row_fraction.columns[row_fraction.loc[cluster] > split_fraction].tolist()
        contributors = col_fraction.index[col_fraction[best] > split_fraction].tolist()
        merges_with = [c for c in contributors if c != cluster] if cluster in contributors else []
        merge_targets = [target for target in counts.columns
                         if col_fraction.loc[cluster, target] > split_fraction
                         and int((col_fraction[target] > split_fraction).sum()) > 1]
        rows.append(dict(cluster=cluster, n_cells=int(row_n.loc[cluster]),
                         best_match=best, best_match_n_cells=int(col_n.loc[best]),
                         intersection=shared, jaccard=float(jaccard[best]),
                         precision=precision, recall=recall,
                         F1=float(2 * shared / (row_n.loc[cluster] + col_n.loc[best])),
                         fraction_retained=recall,
                         fraction_of_cells_retained_together=recall,
                         max_overlap_fraction=float(row_fraction.loc[cluster].max()),
                         n_targets_over_split_fraction=len(split_targets),
                         splits=len(split_targets) > 1,
                         split_targets=_list_json(split_targets),
                         splits_into=_list_json(split_targets if len(split_targets) > 1 else []),
                         merges=len(merges_with) > 0, merges_with=_list_json(merges_with),
                         merge_targets=_list_json(merge_targets),
                         stability=stability, boundary_survives=survives))
    clusters = pd.DataFrame(rows).set_index('cluster')
    # A repeat's fine IDs stay in their own parent block. Only sibling display
    # order follows the baseline row receiving the largest target contribution.
    row_position = dict(zip(counts.index, range(len(counts))))
    comparison_order = sorted(counts.columns, key=lambda c: (
        c.split('.')[0] if '.' in c else '',
        row_position[counts[c].idxmax()], c))
    summary = dict(n_cells=len(left), n_baseline_clusters=len(counts),
                   n_comparison_clusters=len(counts.columns),
                   ARI=float(adjusted_rand_score(left, right)),
                   NMI=float(normalized_mutual_info_score(left, right, average_method='arithmetic')),
                   comparison_is_same_cells=True, bootstrap_stability=False,
                   baseline_order=counts.index.tolist(), comparison_order=comparison_order,
                   thresholds=dict(high_jaccard=high_jaccard, high_precision=high_precision,
                                   high_recall=high_recall, moderate_jaccard=moderate_jaccard,
                                   split_fraction=split_fraction, split_rule='strictly_greater_than',
                                   partial_jaccard=partial_jaccard, partial_recall=partial_recall),
                   stability_counts=clusters.stability.value_counts().to_dict(),
                   boundary_survival_counts=clusters.boundary_survives.value_counts().to_dict())
    return dict(counts=counts, row_fractions=row_fraction, column_fractions=col_fraction,
                clusters=clusters, summary=summary)


def sample_composition(labels, samples, *, expected_fractions=None,
                       max_sample_fraction=.25, min_effective_samples=6,
                       min_cells_for_flag=30):
    """Describe cluster/sample representation and flag substantial imbalance.

    ``labels`` and ``samples`` are cell-indexed Series. By default expected sample
    fractions equal their frequencies in the complete supplied pilot; this is
    exactly 1/12 for the balanced 12-sample, 12,000-cell design. An optional Series
    indexed by sample ID may instead supply strictly positive fractions summing
    to one. It must include all observed samples; additional expected samples
    are retained as zero observed columns, never silently omitted.

    Returns counts, within-cluster ``cluster_fractions``, reciprocal
    ``sample_fractions`` (cluster fraction of each sample), observed/expected
    ``enrichment``, cluster statistics, and a summary with expected fractions and
    rules. Entropy uses natural logarithms; effective samples = exp(entropy).
    Zero contributions add zero entropy. Fractions for completely absent
    expected samples in the reciprocal view are NaN, not invented zeros.

    A descriptive sample_concern requires at least ``min_cells_for_flag`` cells
    and either maximum fraction > ``max_sample_fraction`` or effective samples
    < min(``min_effective_samples``, number of expected samples). Small clusters
    retain every statistic but are marked limited_by_small_cluster. Chi-square
    p-values are descriptive approximations; low expected counts are recorded,
    and no p-value participates in the concern rule. Neither imbalance nor a
    concern flag distinguishes technical batch effects from real biology.
    """
    _fraction(max_sample_fraction, 'max_sample_fraction')
    if min_effective_samples <= 0 or min_cells_for_flag < 1:
        raise ValueError('Sample-count thresholds must be positive.')
    labels = _labels(labels, 'labels')
    samples = _labels(_aligned(labels, samples, 'samples'), 'samples')
    if expected_fractions is None:
        expected = samples.value_counts(sort=False).astype(float) / len(samples)
    else:
        expected = pd.Series(expected_fractions, dtype=float).copy()
        expected.index = expected.index.astype(str)
        if not expected.index.is_unique or expected.isna().any() or (expected <= 0).any():
            raise ValueError('Expected sample fractions must have unique IDs and positive finite values.')
        if not np.isfinite(expected).all() or not np.isclose(expected.sum(), 1., atol=1e-10):
            raise ValueError('Expected sample fractions must sum to one.')
        if not set(samples).issubset(expected.index):
            raise ValueError('Expected fractions omit an observed sample.')
    expected = expected.sort_index()
    counts = pd.crosstab(labels.rename('cluster'), samples.rename('sample'))
    counts = counts.reindex(index=sorted(counts.index), columns=expected.index, fill_value=0)
    counts.columns.name = 'sample'
    n_cells = counts.sum(axis=1)
    cluster_fraction = counts.div(n_cells, axis=0)
    sample_fraction = counts.div(counts.sum(axis=0).replace(0, np.nan), axis=1)
    enrichment = cluster_fraction.div(expected, axis=1)
    rows = []
    for cluster in counts.index:
        fractions = cluster_fraction.loc[cluster]
        nonzero = fractions[fractions > 0]
        entropy = float(-(nonzero * np.log(nonzero)).sum())
        effective = float(np.exp(entropy))
        top_sample = fractions.idxmax()
        n = int(n_cells[cluster])
        expected_counts = n * expected
        chi = chisquare(counts.loc[cluster].to_numpy(), f_exp=expected_counts.to_numpy())
        max_flag = bool(fractions.max() > max_sample_fraction)
        few_flag = bool(effective + 1e-10 < min(min_effective_samples, len(expected)))
        adequate = n >= min_cells_for_flag
        rows.append(dict(cluster=cluster, n_cells=n, top_sample=top_sample,
                         top_sample_fraction=float(fractions.max()),
                         max_sample_fraction=float(fractions.max()),
                         top_sample_enrichment=float(enrichment.loc[cluster, top_sample]),
                         max_sample_enrichment=float(enrichment.loc[cluster].max()),
                         sample_entropy=entropy,
                         normalized_sample_entropy=float(entropy / np.log(len(expected)))
                         if len(expected) > 1 else np.nan,
                         effective_samples=effective, n_samples_present=int((fractions > 0).sum()),
                         n_samples_with_at_least_3_cells=int((counts.loc[cluster] >= 3).sum()),
                         n_samples_absent=int((counts.loc[cluster] == 0).sum()),
                         chi_square=float(chi.statistic), chi_square_p_descriptive=float(chi.pvalue),
                         minimum_expected_count=float(expected_counts.min()),
                         chi_square_small_expected_count=bool((expected_counts < 5).any()),
                         disproportionate_sample_fraction=max_flag,
                         few_effective_samples=few_flag,
                         limited_by_small_cluster=not adequate,
                         sample_concern=adequate and (max_flag or few_flag)))
    clusters = pd.DataFrame(rows).set_index('cluster')
    summary = dict(n_cells=len(labels), n_clusters=len(counts), n_samples=len(expected),
                   expected_fractions=expected.to_dict(),
                   expected_source='observed_pilot_sample_frequencies' if expected_fractions is None
                   else 'explicit_expected_fractions',
                   observed_sample_counts=counts.sum(axis=0).to_dict(),
                   exactly_balanced_observed_samples=bool(counts.sum(axis=0).nunique() == 1),
                   n_sample_concerns=int(clusters.sample_concern.sum()),
                   thresholds=dict(max_sample_fraction=max_sample_fraction,
                                   min_effective_samples=min_effective_samples,
                                   min_cells_for_flag=min_cells_for_flag),
                   interpretation='Descriptive sample association; no causal or artifact assignment; no removal.')
    return dict(counts=counts, cluster_fractions=cluster_fraction,
                sample_fractions=sample_fraction, enrichment=enrichment,
                clusters=clusters, summary=summary)


def phase_matched_identity(labels, parents, phases, program_scores, *,
                           min_cells_per_phase=20, effect_threshold=.5,
                           cycle_tv_threshold=.4):
    """Reassess every within-parent fine pair in matching cell-cycle strata.

    Inputs are Series ``labels`` (fine ID), ``parents`` (coarse ID), ``phases``
    (observed cell-cycle state), and a numeric DataFrame ``program_scores``
    containing independently defined identity/development programs. Callers
    must exclude cell-cycle genes/programs from this identity-score input and
    record that independent feature definition; this function cannot infer
    gene membership from a score column name. Regional and developmental
    programs can be assessed together without fitting a trajectory.

    Effects are differences of means divided by each program's standard
    deviation across *all supplied pilot cells*. This fixed common scale is
    retained in every stratum; within-phase variance shrinkage cannot by itself
    inflate a standardized effect. No expression values are modified. Constant
    or unavailable programs are retained as NaN effects. Each program/phase
    comparison requires ``min_cells_per_phase`` finite scores in each cluster.

    A baseline distinguishing program has abs(unstratified effect) >=
    ``effect_threshold``. A phase retains a difference if at least one such
    program has the same sign and reaches that threshold in the phase. Every
    qualifying program effect, including effects that emerge only within a
    phase, is returned for inspection. The overall persistence status is:

    - insufficient_cells: no phase can assess a baseline distinguishing program
      (or no phase has adequate scores at all);
    - yes: at least two eligible phases, and every eligible phase retains a
      baseline difference;
    - partial: at least one eligible phase retains a baseline difference, but
      only one phase is eligible or other eligible phases do not retain it;
    - no: adequate phase data exist, but none retain a baseline difference, or
      no unstratified identity difference reached the threshold to begin with.

    Returned dict contains ``pairs`` (one row per within-parent pair),
    ``phase_effects`` (one row per pair/phase/program, including an ``all``
    unstratified row), ``phase_composition`` (cluster fractions), ``score_scale``
    (the numerical standardization reference), and ``summary``. Phase-composition
    total variation is half the L1 distance. A ``cell_cycle_boundary_review``
    flag requires TV >= ``cycle_tv_threshold`` AND status no/partial; it is an
    operational review flag, never proof that cell cycle caused the boundary.
    """
    if min_cells_per_phase < 2 or effect_threshold <= 0:
        raise ValueError('At least two cells per phase and a positive effect threshold are required.')
    _fraction(cycle_tv_threshold, 'cycle_tv_threshold')
    labels = _labels(labels, 'labels')
    parents = _labels(_aligned(labels, parents, 'parents'), 'parents')
    phases = _labels(_aligned(labels, phases, 'phases'), 'phases')
    scores = _aligned(labels, program_scores, 'program_scores')
    if not isinstance(scores, pd.DataFrame) or not scores.columns.is_unique or not len(scores.columns):
        raise ValueError('program_scores must have unique, nonempty program columns.')
    scores = scores.astype(float)
    if np.isinf(scores.to_numpy()).any():
        raise ValueError('Program scores may be unavailable (NaN), but cannot be infinite.')
    nesting = pd.DataFrame({'cluster': labels, 'parent': parents}).groupby('cluster').parent.nunique()
    if (nesting != 1).any():
        raise ValueError('A fine cluster spans multiple coarse parents.')
    parent_map = pd.DataFrame({'cluster': labels, 'parent': parents}).drop_duplicates().set_index('cluster').parent
    composition = pd.crosstab(labels.rename('cluster'), phases.rename('phase'))
    composition = composition.div(composition.sum(axis=1), axis=0)
    sd = scores.std(ddof=1).where(lambda x: x > np.finfo(float).eps)
    scale = pd.DataFrame({'pilot_mean': scores.mean(), 'pilot_sd': sd,
                          'finite_cells': scores.notna().sum()})
    scale.index.name = 'program'
    # Cache group means and finite-score counts, avoiding repeated cell scans
    # for the hundreds of sibling-pair and phase combinations.
    means_all = scores.groupby(labels).mean()
    counts_all = scores.groupby(labels).count()
    means_phase = scores.groupby([labels.rename('cluster'), phases.rename('phase')]).mean()
    counts_phase = scores.groupby([labels.rename('cluster'), phases.rename('phase')]).count()
    phase_counts = pd.crosstab(labels, phases)
    effects, pairs = [], []
    for parent in sorted(parent_map.unique()):
        children = sorted(parent_map[parent_map == parent].index)
        for a, b in combinations(children, 2):
            delta = (means_all.loc[a] - means_all.loc[b]) / sd
            all_eligible = (counts_all.loc[a] >= min_cells_per_phase) & (counts_all.loc[b] >= min_cells_per_phase) & sd.notna()
            delta = delta.where(all_eligible)
            original = delta.abs() >= effect_threshold
            for program in scores.columns:
                effects.append(dict(parent=parent, cluster_a=a, cluster_b=b, phase='all',
                                    program=program, n_a=int(counts_all.loc[a, program]),
                                    n_b=int(counts_all.loc[b, program]),
                                    mean_a=float(means_all.loc[a, program]), mean_b=float(means_all.loc[b, program]),
                                    pilot_sd=float(sd[program]), standardized_effect=float(delta[program]),
                                    eligible=bool(all_eligible[program]),
                                    unstratified_distinguishing=bool(original[program]),
                                    same_direction_as_unstratified=bool(original[program]),
                                    identity_difference_retained=bool(original[program])))
            eligible_phases, retained_phases, adequate_any = [], [], []
            retained_programs = set()
            for phase in composition.columns:
                a_present = int(phase_counts.loc[a, phase])
                b_present = int(phase_counts.loc[b, phase])
                eligible = pd.Series(False, index=scores.columns)
                phase_delta = pd.Series(np.nan, index=scores.columns)
                if a_present and b_present:
                    eligible = (counts_phase.loc[(a, phase)] >= min_cells_per_phase) & (counts_phase.loc[(b, phase)] >= min_cells_per_phase) & sd.notna()
                    phase_delta = ((means_phase.loc[(a, phase)] - means_phase.loc[(b, phase)]) / sd).where(eligible)
                if eligible.any():
                    adequate_any.append(phase)
                if (eligible & original).any():
                    eligible_phases.append(phase)
                retained = original & eligible & (phase_delta.abs() >= effect_threshold) & (np.sign(phase_delta) == np.sign(delta))
                if retained.any():
                    retained_phases.append(phase)
                    retained_programs.update(scores.columns[retained])
                for program in scores.columns:
                    effects.append(dict(parent=parent, cluster_a=a, cluster_b=b, phase=phase,
                                        program=program,
                                        n_a=int(counts_phase.loc[(a, phase), program]) if a_present else 0,
                                        n_b=int(counts_phase.loc[(b, phase), program]) if b_present else 0,
                                        mean_a=float(means_phase.loc[(a, phase), program]) if a_present else np.nan,
                                        mean_b=float(means_phase.loc[(b, phase), program]) if b_present else np.nan,
                                        pilot_sd=float(sd[program]), standardized_effect=float(phase_delta[program]),
                                        eligible=bool(eligible[program]),
                                        unstratified_distinguishing=bool(original[program]),
                                        same_direction_as_unstratified=bool(eligible[program] and original[program]
                                                                          and np.sign(phase_delta[program]) == np.sign(delta[program])),
                                        identity_difference_retained=bool(retained[program])))
            if not adequate_any:
                status, reason = 'insufficient_cells', 'no_phase_with_sufficient_finite_identity_scores'
            elif not original.any():
                status, reason = 'no', 'no_unstratified_identity_effect_at_threshold'
            elif not eligible_phases:
                status, reason = 'insufficient_cells', 'baseline_distinguishing_programs_not_assessable_in_any_phase'
            elif len(eligible_phases) >= 2 and len(retained_phases) == len(eligible_phases):
                status, reason = 'yes', 'difference_retained_in_all_of_at_least_two_eligible_phases'
            elif retained_phases:
                status, reason = 'partial', 'difference_retained_in_only_one_or_some_eligible_phases'
            else:
                status, reason = 'no', 'baseline_identity_difference_not_retained_in_eligible_phases'
            tv = float(.5 * np.abs(composition.loc[a] - composition.loc[b]).sum())
            pairs.append(dict(parent=parent, cluster_a=a, cluster_b=b,
                              n_a=int((labels == a).sum()), n_b=int((labels == b).sum()),
                              cell_cycle_total_variation=tv, cycle_composition_tv=tv,
                              strong_cell_cycle_composition_difference=tv >= cycle_tv_threshold,
                              n_unstratified_distinguishing_programs=int(original.sum()),
                              unstratified_distinguishing_programs=_list_json(scores.columns[original]),
                              max_abs_unstratified_identity_effect=float(delta.abs().max()) if delta.notna().any() else np.nan,
                              n_eligible_phases=len(eligible_phases), eligible_phases=_list_json(eligible_phases),
                              n_retained_phases=len(retained_phases), retained_phases=_list_json(retained_phases),
                              retained_programs=_list_json(sorted(retained_programs)),
                              identity_difference_persists_after_cell_cycle_stratification=status,
                              interpretation_reason=reason,
                              cell_cycle_boundary_review=tv >= cycle_tv_threshold and status in ('no', 'partial')))
    pair_columns = ['parent', 'cluster_a', 'cluster_b', 'n_a', 'n_b',
                    'cell_cycle_total_variation', 'cycle_composition_tv', 'strong_cell_cycle_composition_difference',
                    'n_unstratified_distinguishing_programs', 'unstratified_distinguishing_programs',
                    'max_abs_unstratified_identity_effect', 'n_eligible_phases', 'eligible_phases',
                    'n_retained_phases', 'retained_phases', 'retained_programs',
                    'identity_difference_persists_after_cell_cycle_stratification',
                    'interpretation_reason', 'cell_cycle_boundary_review']
    pair_table = pd.DataFrame(pairs, columns=pair_columns)
    status_column = 'identity_difference_persists_after_cell_cycle_stratification'
    summary = dict(n_cells=len(labels), n_pairs=len(pair_table),
                   programs=scores.columns.tolist(),
                   thresholds=dict(min_cells_per_phase=min_cells_per_phase,
                                   effect_threshold=effect_threshold, cycle_tv_threshold=cycle_tv_threshold),
                   effect_scale='fixed_program_standard_deviation_across_all_pilot_cells_ddof_1',
                   status_counts=pair_table[status_column].value_counts().to_dict(),
                   n_cell_cycle_boundary_review=int(pair_table.cell_cycle_boundary_review.sum()),
                   interpretation='Descriptive stratification; no regression, DE, annotation fitting, or causal conclusion.')
    return dict(pairs=pair_table, phase_effects=pd.DataFrame(effects),
                phase_composition=composition, score_scale=scale, summary=summary)
