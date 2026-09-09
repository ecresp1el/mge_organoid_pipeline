#!/usr/bin/env python3
"""Build the phase-two scientific report and browsable Turbo results index."""
from __future__ import annotations

import argparse
import html
import json
from pathlib import Path

import numpy as np
import pandas as pd

from div90_hypergate_paths import PROJECT_ROOT, REPO_ROOT
from div90_hypergate_phase2_gates import STATES, SHORT, apply_rule, clean_json

OUT = PROJECT_ROOT / 'results/div90_hypergate_sst_pv_phase2'
FIGURES = [
    ('figure_01_pv_sst_landscape', '1. Independent PV × SST program landscape'),
    ('figure_02_sample_facets', '2. Cell-line and culture-condition landscapes'),
    ('figure_03_hybrid_validation', '3. Gene-level and technical assessment of dual-high cells'),
    ('figure_04_surface_marker_discovery', '4. State-associated surface markers'),
    ('figure_05_fat3_ptprm_reinterpretation', '5. FAT3/PTPRM capture in the four-state model'),
    ('figure_06_best_sst_depletion', '6. Exploratory SST-depletion gate'),
    ('figure_07_before_after_composition', '7. Starting, removed and retained composition'),
    ('figure_08_condition_biology', '8. Cell-line and culture-condition biology'),
    ('figure_09_pareto_frontier', '9. Recovery versus SST-contamination tradeoffs'),
    ('figure_10_experimental_model', '10. Observations and prospective sorting hypothesis'),
    ('supplementary_landscape_overlays', 'Supplement. Marker overlays on the landscape'),
]
INTERACTIVE = [
    ('interactive_pv_sst_landscape.html', 'Explore the PV × SST landscape'),
    ('interactive_sample_landscape.html', 'Compare cell lines and CV/MW cultures'),
    ('interactive_gate_explorer.html', 'Move gate thresholds and compare recovery costs'),
    ('interactive_umap.html', 'Explore the preserved Loupe UMAP'),
]


def pct(value):
    return 'undefined' if value is None or not np.isfinite(float(value)) else f'{100*float(value):.1f}%'


def num(value, digits=3):
    return 'undefined' if value is None or not np.isfinite(float(value)) else f'{float(value):.{digits}f}'


def read_json(path, default=None):
    return json.loads(path.read_text()) if path.exists() else default


def fraction_table(gate):
    rows = ['| Fraction | Cells | PV-biased | Hybrid | SST-biased | Unresolved |',
            '|---|---:|---:|---:|---:|---:|']
    for fraction in ['removed', 'retained']:
        rows.append('| '+fraction.title()+' | '+str(int(gate[fraction+'_n']))+' | '+
                    ' | '.join(pct(gate[fraction+'_'+s+'_fraction']) for s in SHORT)+' |')
    return rows


def gate_comparison(summary, names):
    rows = ['| Gate | RNA capture rule | SST-biased removal | PV recovery | Hybrid recovery | Retained SST contamination | Retained N | Score |',
            '|---|---|---:|---:|---:|---:|---:|---:|']
    for name in names:
        if name not in summary:
            continue
        g = summary[name]
        label = 'No depletion: retain all' if name == 'baseline' else f'{g["action"]}: `{g["gate_label"]}`'
        rows.append(f'| {name} | {label} | {pct(g["sst_removal_recovery"])} | '
                    f'{pct(g["pv_recovery"])} | {pct(g["hybrid_recovery"])} | {pct(g["sst_contamination"])} | '
                    f'{int(g["retained_n"]):,} | {num(g["practical_score"])} |')
    for g in summary.get('phase1_erbb4_cxcr4_benchmarks', []):
        rows.append(f'| {g["source"]} (inherited) | {g["action"]}: `{g["gate_label"]}` | {pct(g["sst_removal_recovery"])} | '
                    f'{pct(g["pv_recovery"])} | {pct(g["hybrid_recovery"])} | {pct(g["sst_contamination"])} | '
                    f'{int(g["retained_n"]):,} | {num(g["practical_score"])} |')
    return rows


def state_table(cells):
    rows = ['| Operational state | Cells | All entry cells |', '|---|---:|---:|']
    for state in STATES:
        n = int(cells.state.eq(state).sum())
        rows.append(f'| {state} | {n:,} | {pct(n/len(cells))} |')
    return rows


def condition_lines(out):
    path = out / 'tables/condition_baseline.tsv'
    if not path.exists():
        raise FileNotFoundError('Run phase-two conditions before generating its report')
    d = pd.read_csv(path, sep='\t')
    d = d[d.grouping.eq('culture_operator|cell_line|sample')].sort_values('sample')
    rows = ['| Line | Culture/sample | N | PV-biased | Hybrid | SST-biased | Median PV score | Median SST score | SST RNA detected |',
            '|---|---|---:|---:|---:|---:|---:|---:|---:|']
    for r in d.itertuples():
        rows.append(f'| {r.cell_line} | {r.culture_operator} / {r.sample} | {r.n:,} | {pct(r.pv_fraction)} | '
                    f'{pct(r.hybrid_fraction)} | {pct(r.sst_fraction)} | {num(r.pv_score_median)} | '
                    f'{num(r.sst_score_median)} | {pct(r.sst_rna_detected_fraction)} |')
    rows += ['', 'CV denotes higher glucose (~2x) and MW lower glucose (~1x), from the user’s experimental metadata. '
             'The MW substring in every technical sample ID is preserved and does not encode this operator assignment. '
             'CV/MW, glucose and other operator/culture differences are confounded. Each line has one sample per condition; '
             'there is no replicated within-line causal test of glucose.', '',
             'H9 and 79B have higher median PV and lower median SST scores in CV than MW; 2E shows the opposite direction '
             'with only 37 CV entry cells. The line-specific responses therefore do not support one uniform culture effect.', '',
             'Raw marker-expression differences also reflect changing state mixtures. '
             '[The condition report](CONDITION_REPORT.md) provides equal-four-state standardization, the exact decomposition '
             'of raw CV-minus-MW differences into composition and within-state-expression components, and marker-state '
             'association differences within each condition and cell line. These remain descriptive comparisons.']
    return rows


def generate_index(out, cells, hybrid, experimental):
    links = ''.join(f'<a class="button" href="{html.escape(path)}">{html.escape(label)}</a>' for path,label in INTERACTIVE)
    cards = []
    for stem, title in FIGURES:
        formats = ' · '.join(f'<a href="figures/{stem}.{ext}">{ext.upper()}</a>' for ext in ['pdf','svg','png'])
        cards.append(f'<article><h3>{html.escape(title)}</h3><p>{formats}</p>'
                     f'<a href="figures/{stem}.svg"><img loading="lazy" src="figures/{stem}.png" alt="{html.escape(title)}"></a></article>')
    page = f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>DIV90 independent developmental programs</title><style>
body{{font-family:system-ui,sans-serif;color:#243346;background:#f4f6f8;margin:0;line-height:1.55}}
main{{max-width:1200px;margin:auto;padding:35px}}h1{{font-size:32px;line-height:1.2}}h2{{margin-top:34px}}
.lede{{font-size:19px;max-width:950px}}.note{{border-left:4px solid #39876f;background:white;padding:18px}}
.button{{display:inline-block;background:#245b81;color:white;padding:11px 16px;text-decoration:none;border-radius:7px;margin:5px}}
a{{color:#245b81}}.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(340px,1fr));gap:20px}}
article{{background:white;padding:18px;border-radius:8px;box-shadow:0 2px 8px #1d304012}}img{{width:100%;height:auto}}
.metrics{{display:flex;flex-wrap:wrap;gap:20px}}.metrics div{{background:white;padding:18px;min-width:180px;border-radius:8px}}
.metrics b{{display:block;font-size:27px}}code{{font-size:13px}}footer{{margin:35px 0;color:#526172}}
</style><main><h1>DIV90 PV × SST developmental landscape</h1>
<p class="lede">Independent continuous programs in the preserved cortical LHX6+/ERBB4+ population. Four operational groups support exploration; they are not established cell fates.</p>
<div class="metrics"><div><b>{len(cells):,}</b>frozen entry cells</div><div><b>{int(hybrid['hybrid_n']):,}</b>operational dual-high cells ({pct(hybrid['hybrid_fraction'])})</div>
<div><b>{int(hybrid['hybrid_multiple_genes_both_n']):,}</b>dual-high cells detecting ≥2 genes per program</div></div>
<p class="note">The overall dual-high population is not enriched over the sample/depth-stratified independent-gene null. Weak gene coordination and RNA complexity limit claims of a distinct hybrid population. Preserve these cells as an experimental population and evaluate maturation.</p>
<h2>Interactive exploration</h2>{links}
<h2>Scientific interpretation and tables</h2><p><a href="REPORT.md">Full report: all 12 biological questions and proposed sort</a> ·
<a href="CONDITION_REPORT.md">CV/MW condition analysis</a> · <a href="provenance/DIV90_HYPERGATE_PHASE2_EVIDENCE.md">Gene, metadata and reagent evidence</a> ·
<a href="figures/FIGURE_CAPTIONS.md">Figure captions</a></p>
<p>The prospective depletion hypothesis is <code>{html.escape(experimental['gate_label'])}</code> inside ERBB4+ cells, with the captured RNA region removed. RNA thresholds are not fluorescence thresholds, and marker staining must be validated before sorting.</p>
<p><a href="cells.tsv.gz">Every cell and its scores</a> · <a href="gate_cell_assignments.tsv.gz">Gate assignments</a> ·
<a href="tables/condition_gate_metrics.tsv">Condition and line gate metrics</a> · <a href="tables/gate_shortlist.tsv">Gate shortlist</a> ·
<a href="tables/gate_pareto.tsv">Pareto frontier</a> · <a href="gate_summary.json">Exact gate rules and metrics</a></p>
<h2>Figures</h2><div class="grid">{''.join(cards)}</div>
<footer>Phase two • Turbo runtime only • Phase-one files and original Loupe coordinates preserved • All analysis is descriptive discovery; no prospective protein or held-out sample validation is implied.</footer></main></html>'''
    (out / 'index.html').write_text(page)


def run(cells: pd.DataFrame, out: Path | str, summary: dict):
    out = Path(out).resolve()
    if out != OUT.resolve() or REPO_ROOT.resolve() in out.parents:
        raise ValueError(f'Report output must be {OUT}')
    hybrid = read_json(out / 'hybrid_validation.json')
    params = read_json(out / 'state_parameters.json')
    experimental = summary.get('experimental_depletion', summary['selected'])
    reference = summary['phase1_reference']
    for gate in [experimental, reference]:
        removed = apply_rule(cells, summary['rules'][gate['gate_id']])
        if int(removed.sum()) != int(gate['removed_n']):
            raise AssertionError('Report gate counts do not match the frozen cells')
    retained = cells.loc[~apply_rule(cells, summary['rules'][reference['gate_id']])]
    pvgenes = params['modules']['pv']['genes']
    sstgenes = params['modules']['sst']['genes']
    experimental_genes = {r['gene'] for r in summary['rules'][experimental['gate_id']]['rules']}
    reagent_note = ''
    if {'FGFR2','PTPRS'}.issubset(experimental_genes):
        reagent_note = ('A concrete first staining panel can evaluate **ERBB4 182818 + FGFR2 98739 (MAB6843) + PTPRS AF3430**. '
            'FGFR2 98739 has manufacturer human Kato III flow validation and broad FGFR2 isoform coverage. '
            'PTPRS AF3430 is a goat polyclonal against extracellular Glu30–Gly1260; its datasheet validates WB/ICC, '
            'while Bunin et al. provide human PBMC/pDC nonpermeabilized surface-flow evidence linked to this product. '
            'These support a protein pilot in DIV90 neurons. Use distinguishable direct conjugates or validated detection: '
            'ERBB4 and FGFR2 are both mouse antibodies, so a common anti-mouse secondary cannot separate them. '
            'Measure viability and maturation because receptor binding/crosslinking may affect function. '
            'No fluorescence boundary or neuron-sorting performance is established. '
            '[FGFR2 manufacturer](https://www.rndsystems.com/products/human-fgfr2-antibody-98739_mab6843); '
            '[PTPRS extracellular-immunogen datasheet](https://resources.rndsystems.com/pdfs/datasheets/af3430.pdf); '
            '[Primary human PTPRS surface flow](https://pmc.ncbi.nlm.nih.gov/articles/PMC4547994/).')
    if 'SLCO3A1' in experimental_genes:
        reagent_note = ('For **SLCO3A1/OATP3A1**, the inspected human antibodies are WB-only '
            '(Thermo PA5-42457; Sigma SAB1304633) or IHC-only (Thermo PA5-84720). No verified human '
            'nonpermeabilized antibody-flow clone was found. Thus the first concrete experiment is a protein-validation '
            'pilot before attempting the proposed depletion; a fluorescence threshold cannot yet be specified. '
            '[Thermo PA5-42457](https://www.thermofisher.com/antibody/product/SLCO3A1-Antibody-Polyclonal/PA5-42457); '
            '[Sigma SAB1304633](https://www.sigmaaldrich.com/US/en/product/sigma/sab1304633).')
    lines = ['# DIV90 Hypergate phase 2: independent programs and prospective depletion', '',
        'This run preserves the completed phase-one outputs, all 4,768 entry-cell identities and original Loupe coordinates. '
        'It evaluates PV-associated and SST-associated programs as separate continuous axes. PVALB detection is not required.', '',
        'The principal result is a heterogeneous developmental landscape with operational program overlap. '
        'The data do not establish four discrete fates or a distinct coordinated hybrid cell type. '
        'SST-negative cells are not equated with PV identity, and SST RNA is not treated as irreversible commitment.', '',
        f'The practical depletion candidate changes SST-biased contamination from '
        f'{pct(cells.state.eq("SST-biased").mean())} at baseline to {pct(experimental["sst_contamination"])} retained, '
        f'while retaining {pct(experimental["pv_recovery"])} of PV-biased and {pct(experimental["hybrid_recovery"])} '
        'of operational hybrid cells. The biological value must be judged against those recovery costs and protein feasibility.', '',
        f'This is only a {100*(cells.state.eq("SST-biased").mean()-experimental["sst_contamination"]):.2f} percentage-point '
        'reduction in SST-biased contamination. The candidate preserves much of the desired population but does not deliver clean SST depletion.', '',
        '[Open the visual results index](index.html) · [Interactive developmental map](interactive_pv_sst_landscape.html) · '
        '[Condition analysis](CONDITION_REPORT.md) · [Gene and reagent evidence](provenance/DIV90_HYPERGATE_PHASE2_EVIDENCE.md)', '',
        '## 1. Does DIV90 contain four distinct developmental states?', '',
        'Four operational groups are assigned using independent pooled median thresholds, with every cell retained. '
        'This partition is useful for discovery but does not itself prove biological clusters. Mixture fits and KDE '
        'diagnostics are saved; sparse RNA detection and score floors complicate their interpretation.', '',
        f'PV program: `{", ".join(pvgenes)}`. SST program: `{", ".join(sstgenes)}`. Scores average gene-wise '
        'population z-scores of the frozen log1p(CP10K) expression, clipping each gene to [−3,3]. '
        'No opposing score is subtracted. Generic GAD genes and the ERBB4 entry feature carry no primary score weight. '
        'All score-defining genes are excluded from surface-gate discovery.', '',
        'MEF2C has developmental PV evidence but shared regulatory functions; KCNC genes capture a fast-spiking maturation program. '
        'The gene-by-gene rationale and context limitations, including SATB1 and NR2F2, are documented in the evidence note. '
        '[Mayer et al., 2018](https://www.nature.com/articles/nature25999); '
        '[Allaway et al., 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC9316417/).', '',
        f'Primary cutoffs: PV ≥ {params["pv_threshold"]:.9g}; SST ≥ {params["sst_threshold"]:.9g}. '
        'Low/low means lower measured programs under these cutoffs; immaturity is a hypothesis, not an independently proven identity.', '']
    lines += state_table(cells)
    variance_path = out / 'tables/module_variance_contributions.tsv'
    if variance_path.exists():
        variance = pd.read_csv(variance_path, sep='\t')
        sst_share = variance.loc[variance.gene.eq('SST'), 'score_variance_share'].iloc[0]
        lines += ['', f'SST contributes {pct(sst_share)} of the SST-score variance under the covariance decomposition '
            '`Cov(gene contribution, module score) / Var(module score)`. Thus SST alone does not dominate this score’s '
            'observed variance. This diagnostic does not establish gene coordination. '
            '[All gene contributions](tables/module_variance_contributions.tsv).']
    lines += ['', '## 2. How large is the operational hybrid population?', '',
        f'**{hybrid["hybrid_n"]:,} cells ({pct(hybrid["hybrid_fraction"])})** are dual-high at the primary thresholds. '
        'This number is threshold dependent. All 35/50/65% threshold combinations and leave-one-gene-out alternatives are retained. '
        'When an alternative quantile lands at an all-undetected floor, the entire tied floor is assigned low and the achieved '
        'high fractions are reported; cells are never split by arbitrary barcode ordering.', '',
        '## 3. Is dual-high status supported by multiple genes and technical controls?', '',
        f'{hybrid["hybrid_multiple_genes_both_n"]:,} dual-high cells ({pct(hybrid["hybrid_multiple_genes_both_fraction"])}) '
        'detect at least two genes from each program. '
        f'{hybrid["hybrid_multiple_genes_both_without_SST_n"]:,} '
        f'({pct(hybrid["hybrid_multiple_genes_both_without_SST_n"]/hybrid["hybrid_n"])}) additionally detect at least two SST-module genes other than SST. '
        'These are co-detection-supported subsets, not validated future hybrid fates.', '',
        f'The independent-gene null preserving sample and UMI quintile gives a mean of {hybrid["null_dual_high_mean"]:.1f} '
        f'dual-high cells (95% empirical interval {hybrid["null_dual_high_95_interval"][0]:.1f}–{hybrid["null_dual_high_95_interval"][1]:.1f}), '
        f'versus {hybrid["hybrid_n"]:,} observed. Overall dual-high abundance is therefore not enriched over this null. '
        f'The multiple-gene dual-high count is {hybrid["hybrid_multiple_genes_both_n"]:,}, compared with null mean '
        f'{hybrid["null_multigene_mean"]:.1f} (interval {hybrid["null_multigene_95_interval"][0]:.1f}–{hybrid["null_multigene_95_interval"][1]:.1f}), '
        'a modest excess that merits follow-up rather than proof of a new population.', '',
        f'PV/SST score Spearman correlation is {num(hybrid["rho_pv_sst"])} and falls to '
        f'{num(hybrid["rho_depth_diagnostic"])} after the within-sample depth diagnostic, which preserves sample means. '
        f'Median within-module gene correlations are {num(hybrid["median_within_pv_gene_rho"])} for PV and '
        f'{num(hybrid["median_within_sst_gene_rho"])} for SST. Gene coordination is weak. '
        'Depth matching is approximate and does not eliminate every complexity effect.', '',
        'QC tables compare library counts, detected genes, mitochondrial fraction, stress/activation, cell cycle and exploratory doublets. '
        +hybrid.get('doublet_assessment','No doublet assessment available.')+' '
        f'Exploratory predicted doublets: {hybrid.get("exploratory_predicted_doublets_total", "unavailable")} among entry cells, '
        f'{hybrid.get("exploratory_predicted_doublets_hybrid", "unavailable")} dual-high. '
        'A low predicted count cannot rule out homotypic doublets; no primary state is removed for being inconvenient.', '',
        'The `phase1_genes_independent_axes` sensitivity uses the earlier gene set with phase-two score processing. '
        'It is not an exact reproduction of phase-one unbounded scores or extreme labels, which remain preserved separately.', '',
        '## 4. Which cell lines and culture conditions favor each operational state?', '']
    lines += condition_lines(out)
    lines += ['', '## 5. What does FAT3-low/PTPRM-undetected actually enrich?', '',
        f'The exact frozen RNA rule captures **{len(retained):,} cells**: `{reference["gate_label"]}`. '
        'The table below includes every captured or excluded cell, including ambiguous cells from the older binary analysis.', '']
    lines += fraction_table(reference)
    lines += ['', f'Compared with the starting population, retained PV score median changes from {num(cells.pv_score.median())} '
        f'to {num(retained.pv_score.median())}, and SST score median from {num(cells.sst_score.median())} '
        f'to {num(retained.sst_score.median())}. PV-biased recovery is {pct(reference["pv_recovery"])}, '
        f'hybrid recovery {pct(reference["hybrid_recovery"])}, and unresolved recovery {pct(reference["unresolved_recovery"])}. '
        'These four-state composition and recovery values, rather than the old ambiguous-excluded purity, define what the gate captures.', '',
        f'**The old gate loses {int(reference["hybrid_lost_n"]):,} hybrid cells ({pct(1-reference["hybrid_recovery"])}).** '
        'It enriches PV-biased cells while reducing the operational hybrid fraction; it should not be interpreted as preserving '
        'the full PV-associated developmental continuum. This hybrid reduction occurs in both pooled CV and MW conditions.', '',
        '## 6. Which single marker best identifies SST-biased cells for depletion?', '',
        f'The unconstrained best single-marker rule is `{summary["best_single"]["gate_label"]}`. '
        'Because preservation-heavy scores can favor almost no depletion, the practical comparison also requires an '
        'analyst-chosen minimum 20% removal of SST-biased cells. This threshold is a decision constraint, not a natural '
        'biological boundary or a requirement supplied by the user.', '',
        f'Under that constraint, the best single is `{summary.get("practical_single",summary["best_single"])["gate_label"]}`. '
        'The 10%, 20%, 30% and 50% minimum-removal sensitivity is saved. A marker can correlate with an operational state '
        'without being specific enough for useful live sorting.', '',
        '## 7. Which gate retains PV-biased plus hybrid cells while removing SST-biased cells?', '']
    lines += gate_comparison(summary, ['baseline','phase1_reference','selected','practical_single','practical_pair','experimental_depletion','erbb4_cxcr4'])
    lines += ['', 'The inherited phase-one ERBB4/CXCR4 benchmark is a positive **retention** rule, re-evaluated against the new four states. '
        'The `erbb4_cxcr4` row is the new positive **removal** sweep benchmark; these have different directions and purposes.', '',
        f'The experimental-depletion hypothesis is **remove `{experimental["gate_label"]}`**, retain its complement '
        'inside the existing entry population. The unconstrained `selected` optimum remains separate and visible; it may '
        'achieve better numerical score by making a much smaller intervention.', '',
        'The practical score is `sqrt(PV recovery × hybrid recovery) × (1 − SST contamination among retained resolved states)`. '
        'Unresolved cells remain present in the actual retained fraction. A separate score includes them in the purity denominator. '
        'A Pareto frontier shows combined PV/hybrid recovery versus SST contamination with total retained cells visible.', '',
        '## 8–9. How many PV-biased and hybrid cells would the experimental gate lose?', '',
        f'Predicted loss: **{int(experimental["pv_lost_n"]):,} PV-biased cells** '
        f'({pct(1-experimental["pv_recovery"])}) and **{int(experimental["hybrid_lost_n"]):,} hybrid cells** '
        f'({pct(1-experimental["hybrid_recovery"])}) from the frozen starting population. '
        'These are in-sample RNA predictions, not measured post-sort survival.', '',
        '## 10. How much SST-biased contamination remains?', '',
        f'{pct(experimental["sst_contamination"])} of all retained cells are operationally SST-biased. '
        f'The stricter high-SST/low-PV extreme definition comprises {pct(experimental.get("strong_sst_contamination"))} '
        f'of retained cells, with {pct(experimental.get("strong_sst_removal_recovery"))} of those strong SST-biased cells removed. '
        'That strict definition is an additional sensitivity, not the target used to silently re-label all SST-biased cells.', '',
        f'If unresolved cells are excluded from evaluation, remaining SST-biased contamination becomes '
        f'{pct(experimental["sst_contamination_unresolved_excluded"])}. This is a denominator/label sensitivity; '
        'no additional surface gate that selectively removes unresolved cells has been demonstrated.', '',
        '## 11. Is one depletion marker sufficient?', '',
        ('The constrained two-marker rule satisfies the analyst-defined practical improvement criterion.' if summary.get('practical_pair_is_meaningful')
         else 'The constrained two-marker rule does not satisfy the analyst-defined practical improvement criterion, so the simpler one-marker rule is preferred for the first hypothesis.'),
        'The criterion requires at least 0.02 absolute score improvement with neither PV nor hybrid recovery decreasing by '
        'more than 0.03. It is a transparent experimental preference, not a statistical discovery test. '
        'Bounded three-marker exploration and its stricter improvement criterion are documented in `gate_summary.json`; '
        'the search is not exhaustive across all surfaceome triplets.', '',
        '## 12. Proposed DIV90 sort — hypothesis to test', '',
        '**Gate 1:** viable singlets → ERBB4+ cortical interneuron-enriched cells. '
        'The RNA analysis was further restricted by recovered cortical identity and LHX6 expression; these intracellular/transcript '
        'criteria are not automatically replicated by a live ERBB4-only entry gate. The predictions below apply to that frozen '
        'cortical LHX6+/ERBB4+ population and need validation in the actual prospective entry population.', '',
        f'**Gate 2:** test the surface-protein analogue of `{experimental["gate_label"]}`; remove the positive captured region and retain its complement.', '',
        '**REMOVE — predicted composition:**', '',
        f'{pct(experimental["removed_sst_fraction"])} SST-biased; {pct(experimental["removed_hybrid_fraction"])} hybrid; '
        f'{pct(experimental["removed_pv_fraction"])} PV-biased; {pct(experimental["removed_unresolved_fraction"])} unresolved '
        f'({int(experimental["removed_n"]):,} cells).', '',
        '**RETAIN — predicted composition:**', '',
        f'{pct(experimental["retained_pv_fraction"])} PV-biased; {pct(experimental["retained_hybrid_fraction"])} hybrid; '
        f'{pct(experimental["retained_sst_fraction"])} SST-biased; {pct(experimental["retained_unresolved_fraction"])} unresolved '
        f'({int(experimental["retained_n"]):,} cells).', '',
        '**Recovery:** '+pct(experimental['pv_recovery'])+' of all PV-biased cells and '+pct(experimental['hybrid_recovery'])+
        ' of all operational hybrid cells retained; '+pct(experimental['total_recovery'])+' total entry-cell recovery.', '',
        'This is an exact **RNA-based sorting hypothesis**, not a ready-calibrated fluorescence protocol. No log1p(CP10K) '
        'cutoff translates directly to protein abundance or a fluorescence threshold. Begin with a nonpermeabilized staining '
        'pilot, extracellular epitope/clone validation, viability and dissociation controls, and linked RNA/protein validation '
        'of the proposed high/low fractions. ERBB4 clone 182818 has manufacturer human-flow evidence; the inherited FAT3/PTPRM '
        'pair lacks demonstrated human DIV90 live-sorting validation. See the reagent evidence note for precise catalog/application limits.', '',
        reagent_note, '',
        'Compare retained and removed fractions with an unsorted ERBB4+ control under matched subsequent culture conditions, '
        'preserving all cell-line and CV/MW labels. Measure immediate multi-gene composition, survival/absolute recovery, then '
        'later PVALB/SST protein and electrophysiological maturation. Fast-spiking properties and repeated multi-gene assessment '
        'are more informative than relying on later SST negativity alone. This directly tests whether the retained program '
        'mixture contains cells able to acquire the desired mature phenotype.', '',
        '## Search scope, validation, and files', '',
        f'{summary["total_surface_features"]} expressed curated surface features were evaluated after target-gene exclusion. '
        f'Actual R Hypergate fits: {summary["actual_R_hypergate_fits"]:,}; exact one-marker threshold candidates: '
        f'{summary["exact_single_candidates"]:,}; broad positive pair candidates: {summary["broad_pair_candidates"]:,}. '
        'Signed state screens cover PV-biased, SST-biased and hybrid enrichment, PV-versus-SST fitting that leaves hybrid '
        'cells unclassified, and PV-plus-hybrid retention. Python depletion sweeps are separately identified from actual R Hypergate.', '',
        summary['threshold_search'], '',
        'All feature selection and primary thresholds use the pooled population. Per-sample and CV/MW evaluations apply '
        'fixed gates and are not held-out validation. The full cell-level files, threshold grids, confusion counts and '
        'provenance enable direct re-evaluation.', '',
        '- `cells.tsv.gz`: every frozen cell, independent scores, all component expression, metadata and QC.',
        '- `gate_cell_assignments.tsv.gz`, `gate_summary.json`: exact rules and per-cell capture flags.',
        '- `tables/depletion_all_candidate_metrics.h5`: complete saved Python sweep metrics.',
        '- `tables/depletion_supplemental_candidate_metrics.h5`: constrained pair refinements and bounded three-marker exploration.',
        '- `tables/gate_shortlist.tsv`, `tables/gate_pareto.tsv`: practical candidates and frontier.',
        '- `tables/module_threshold_sensitivity.tsv`, `tables/cell_sensitivity.tsv.gz`: complete label sensitivities.',
        '- `tables/module_gene_by_state.tsv`, `hybrid_validation.json`: gene-level and null-model evidence.',
        '- `CONDITION_REPORT.md`, `tables/condition_*.tsv*`: sample, line, culture and marker-association comparisons.',
        '- `provenance/phase1_preservation.json`: SHA256 audit of unchanged phase-one files.',
        '- `figure_manifest.json`, `figures/FIGURE_CAPTIONS.md`: PDF/SVG/600-dpi PNG figures and captions.',
        '- `index.html` and the four `interactive_*.html` files: offline visual exploration.']
    validation_path = out / 'provenance/final_validation.json'
    if validation_path.exists():
        validation = read_json(validation_path)
        lines += ['', f'Independent validation checked count conservation for '
            f'{validation["total_python_candidate_rows_count_conservation_checked"]:,} saved candidate-metric rows, '
            f'recomputed {validation["random_candidate_masks_independently_reproduced"]:,} randomly selected candidate masks, '
            'and reproduced key gate counts and per-cell assignments. The audit also verifies preserved entry-cell order, '
            'original Loupe coordinates, original phase-one scores/labels and unchanged phase-one files. '
            '[Validation audit](provenance/final_validation.json).']
    (out / 'REPORT.md').write_text('\n'.join(lines)+'\n')
    generate_index(out, cells, hybrid, experimental)
    return {'report': str(out / 'REPORT.md'), 'index': str(out / 'index.html')}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, default=OUT)
    args = parser.parse_args()
    run(pd.read_csv(args.out / 'cells.tsv.gz', sep='\t'), args.out,
        json.loads((args.out / 'gate_summary.json').read_text()))
