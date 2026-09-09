#!/usr/bin/env python3
"""Assemble the biological identity report from completed, fixed-target diagnostics."""
from pathlib import Path
import json
import shutil
import html
import pandas as pd
import numpy as np
import mistune
ROOT=Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results')
OUT=ROOT/'div30_antecedent_identity_v1'
REPO=Path(__file__).resolve().parents[2]

def table(df):
    def fmt(x):
        if pd.isna(x):return '—'
        if isinstance(x,(float,np.floating)):return f'{x:.3f}'
        return str(x).replace('|',' / ').replace('\n',' ')
    rows=['| '+' | '.join(df.columns)+' |','| '+' | '.join(['---']*len(df.columns))+' |']
    rows += ['| '+' | '.join(map(fmt,row))+' |' for row in df.itertuples(index=False,name=None)]
    return '\n'.join(rows)

def read(name):return pd.read_csv(OUT/'within/tables'/name,sep='\t',float_precision='round_trip')

def main():
    g=read('group_summary.tsv').set_index('group')
    genes=read('group_gene_summary.tsv')
    models=read('model_oof_metrics.tsv').set_index(['outcome','model'])
    cfg=json.loads((OUT/'within/provenance/frozen_gate.json').read_text())
    ext_summary=(OUT/'external/ROOT_INTERPRETATION.md').read_text()
    ext_data=json.loads((OUT/'external/summary.json').read_text())
    matched_total=sum(r['n_external_matched'] for r in ext_data['matching'])
    def mm(outcome,kind,metric='roc_auc'):return float(models.loc[(outcome,kind),metric])
    def pct(x):return f'{x:.1%}'
    def a2b(outcome,base,aug,metric='roc_auc'):return f'{mm(outcome,base,metric):.3f} → {mm(outcome,aug,metric):.3f}'
    def gp(group,gene):return float(genes[(genes.group==group)&(genes.gene==gene)].positive_fraction.iloc[0])
    sections=[]
    def add(text):sections.append(text.strip())
    add('''# WHAT IS THE DIV30 ENRICHABLE ANTECEDENT STATE?

The best-supported identity is a **SLC6A1-rich, predominantly postmitotic GABAergic developmental state with enriched PV-associated transcription and substantial SST-associated expression**. The gate selects a biologically related population, but the evidence does not establish a prospective PV-producing population or a distinct PV-specific precursor. Forward and reciprocal mapping preferentially connect the retained population to DIV90 **dual-high** territory, not selectively to strict PV-biased cells.

The two hypotheses are not mutually exclusive: a broader developmental state could contain cells that later produce PV. Dual-high is also an operational transcriptional label, not proof of a stable mixed fate. The present evidence addresses RNA-state selectivity and resemblance; it cannot measure later PV production by the retained cells.

This extension preserves the original top20 consensus, every original score and the exact GRIA2/OPCML/NOTCH1 rule. It performs no large gate search. The retained 15,781 cells contain 9,210 frozen targets: 20.0% → 58.36% composition, 50.81% target recovery, 17.41% total yield. The previous held-sample **gate reselection** result remains 53.49% composition and 39.73% recovery; the new biological prediction models below are a different analysis.

“Target” means the pre-existing operational top20 score, not future PV fate. “Postmitotic” is an RNA-based operational proxy, not direct cell-cycle exit measurement. Neighbor fractions are empirical resemblance measures, not calibrated fate probabilities.''')
    answers=[
      ('1. More than generic postmitotic maturation?',f'Yes at the level of measured association, with a large generic developmental component. Flexible generic covariates predict fixed-gate membership at AUC {mm("fixed_gate_retained","A_flexible_generic"):.3f}; adding the PV module gives {mm("fixed_gate_retained","B_flexible_plus_PV"):.3f}. This is additional information, not demonstrated PV specificity.'),
      ('2. Association with an independent PV program?',f'The literature-motivated four-gene PV module is enriched. Exact requested target-model AUC is {a2b("antecedent_target","A_requested","B_requested_plus_PV")}, but that module already contributes to the target. It is not a statistically independent endpoint.'),
      ('3. Also SST-associated?',f'Yes. SST RNA is detected in {gp("retained","SST"):.1%} retained cells and {gp("A","SST"):.1%} captured targets. The composite rises modestly pooled but falls within five of six samples; this is not uniform SST-module enrichment.'),
      ('4. SLC6A1 enriched?',f'Yes: {g.loc["all","SLC6A1_positive_fraction"]:.1%} all cells, {g.loc["retained","SLC6A1_positive_fraction"]:.1%} retained, {g.loc["A","SLC6A1_positive_fraction"]:.1%} captured targets. It is a partially separable GABA-handling axis, but the old mapping/sparse consensus already uses SLC6A1.'),
      ('5. Contributions of each marker?','GRIA2 marks a neuronal transition; OPCML adds the largest conditional target enrichment; NOTCH1-low removes a cycling/progenitor-rich remainder. The third clause adds only 2.92 percentage points of composition.'),
      ('6. What are the retained non-target cells?','Mostly related inhibitory/maturing neurons: 90.47% postmitotic, 94.81% GAD1/2+, 87.84% NKX2-1+. 41.26% lie in the 70th–80th score percentiles; 69.41% in the 60th–80th; 14.73% are below the median. They are not all near-cutoff targets.'),
      ('7. Captured versus lost targets?','Captured targets are more postmitotic, less progenitor/cycling, higher in SLC6A1 and SST expression, and modestly higher in the PV module. The gate selects a subset with more mature expression and higher OPCML RNA; adhesive function was not measured.'),
      ('8. Preferential forward mapping to later PV states?','To dual-high territory, yes; to strict PV-biased territory, no. Strict-PV neighbor fraction is 11.55% retained versus 16.95% removed; dual-high is 52.95% versus 34.92%.'),
      ('9. Reciprocal strict-PV convergence?','Weak and nonspecific. Strict-PV queries have 20.37% retained DIV30 neighbors, versus 23.15% for SST and 31.60% for dual-high. Sample-balanced strict-PV enrichment is only 1.02×; essentially all reverse queries exceed the internal DIV30 held-sample 95th-percentile distance threshold.'),
      ('10. Why more separable at DIV30?','The early rule separates broad developmental states across an all-cell population. The same three genes barely separate DIV90 PV-biased from SST-biased cells. Some DIV30 marker differences persist within postmitotic and matched entry populations; age, target definition and population boundary remain entangled.'),
      ('11. Unusually abundant versus other protocols?','Not established as a protocol-level property. See the local external mapping analysis below: query resemblance, sample selection, stage and out-of-reference coverage must accompany every frequency.'),
      ('12. Distinguishing matched-state programs?',('Not estimable: no external cells passed the declared resemblance/support and matching criteria. No whole-study comparison was substituted, and no MGEO-specific pathway or causal explanation for PV protein was established.' if matched_total==0 else 'Exploratory, comparison-specific candidates are reported below with matching balance and coverage. These data do not establish MGEO-specific pathways or a causal explanation for PVALB protein.')),
      ('13. Which hypothesis is more credible?','The broader postmitotic GABAergic-state interpretation is better supported. A PV-producing component remains possible but unverified. Neither the arbitrary top20 definition nor shared-space mapping establishes lineage.')]
    add('## Answers to the 13 biological questions\n\n'+table(pd.DataFrame(answers,columns=['Question','Evidence-based answer'])))

    display=[]
    for group,label in [('A','Target + retained'),('B','Target + lost'),('C','Non-target + retained'),('D','Non-target + removed')]:
        r=g.loc[group];display.append({'Group':group+' · '+label,'Cells':int(r.n),'PV score':r.direct_pv_score_mean,'SST score':r.direct_sst_score_mean,'MGE score':r.mge_score_mean,'Postmitotic':pct(r.postmitotic_mean),'Cycling':pct(r.cycling_mean),'Progenitor':pct(r.progenitor_mean),'SLC6A1+':pct(r.SLC6A1_positive_fraction),'Detected genes':round(r.n_genes_mean)})
    add('## 1. Four-group single-cell identity\n\n'+table(pd.DataFrame(display))+'''

A and C have almost identical postmitotic rates despite their different target labels. C is therefore not primarily an immature/progenitor contamination compartment. B contains more cycling/progenitor evidence than A. All groups retain heterogeneous gene expression; none is assigned a new exclusive lineage.

The frozen MGE composite is **not enriched** in the target or gate. Its component genes do not move together: NKX2-1/SOX6 and progenitor-associated components can fall while LHX6/ERBB4 and MAF/MAFB rise. Do not replace this with the blanket claim “stronger MGE identity.” The retained cells nevertheless show substantial NKX2-1 and GABA-synthesis evidence.

The table below reports mean log1p(counts per 10,000) RNA followed by percent detected. Direct PV/SST scores are means of clipped gene z-scores standardized to the frozen DIV90 reference; MGE, maturation and supplementary context scores use their documented DIV30 gene standardization. The antecedent consensus is the mean of four percentile ranks. Scores are not on one interchangeable biological scale. All explicitly requested genes are assayed. TUBB3, present as a zero placeholder in the old context panel, is not assayed in the source H5AD; its historical contribution remains unchanged and is flagged in the availability audit.''')
    requested=['NKX2-1','LHX6','SOX6','DLX1','DLX2','DLX5','DLX6','MAF','MAFB','MEF2C','ERBB4','KCNC1','KCNC2','SLC6A1','SST','SATB1','TAC1','PVALB']
    rows=[]
    for gene in requested:
        row={'Gene':gene}
        for group in 'ABCD':
            r=genes[(genes.gene==gene)&(genes.group==group)].iloc[0];row[group]=f'{r["mean"]:.3f}; {r.positive_fraction:.2%}'
        rows.append(row)
    add(table(pd.DataFrame(rows))+'''

PVALB RNA is rare and is not enriched: 9/18,127 targets versus 224/72,504 non-targets. This argues against calling the target an already PVALB-expressing population; its rarity at DIV30 does not rule out future PV production.

Full distributions of maturation, progenitor, progression, GABA, cell cycle, immediate-early/stress context, detected genes and mitochondrial fraction are in figures 01–02 and the gene/program summary tables. Immediate-early/stress RNA is not proof of pathological stress, and RNA complexity differences are covariates rather than fate measurements.''')

    modelrows=[]
    for outcome,label in [('antecedent_target','Frozen target'),('fixed_gate_retained','Fixed gate')]:
        for model in ['A_requested','B_requested_plus_PV','A_core','B_core_plus_PV','A_extended_generic','B_extended_plus_PV','A_flexible_generic','B_flexible_plus_PV','B_extended_nonmodule_context']:
            r=models.loc[(outcome,model)];modelrows.append({'Outcome':label,'Model':model,'ROC AUC':r.roc_auc,'Average precision':r.average_precision,'Log loss':r.log_loss,'Loss reduction vs null':pct(r.predictive_log_loss_reduction_fraction)})
    add('## 2. How much is generic maturation?\n\n'+f'''The exact requested Model A uses progenitor score, neuronal maturation score and operational postmitotic status. Model B adds the direct MEF2C/KCNC1/KCNC2/TAC1 PV score. The DCX-augmented core is a sensitivity; extended models also use cycle scores/phase, GABA context excluding SLC6A1, RNA counts, detected genes and mitochondrial fraction. Flexible models use additive cubic splines. Scalers, spline knots and model coefficients are fit on five samples and evaluated on the sixth, without retuning the gate or target.

For fixed-gate membership, the exact A→B AUC is {a2b('fixed_gate_retained','A_requested','B_requested_plus_PV')}; extended {a2b('fixed_gate_retained','A_extended_generic','B_extended_plus_PV')}; flexible {a2b('fixed_gate_retained','A_flexible_generic','B_flexible_plus_PV')}. The extended log loss changes {a2b('fixed_gate_retained','A_extended_generic','B_extended_plus_PV','log_loss')}. This supports a large generic developmental contribution and a smaller residual PV-program association. The residual is not evidence of selective future PV fate.

For the antecedent label, a much larger PV-module gain is expected by construction: direct and ranked PV scores are already two of its four components. Prediction cannot serve as independent biological validation of that label. “Loss reduction versus null” quantifies predictive information in this model family, not a percentage of causal biology or variance explained.

'''+table(pd.DataFrame(modelrows))+'''

The MAF/MAFB/ERBB4/LHX6/SOX6 sensitivity is excluded from the explicit PV module and original mapping features, but is broad developmental context, not an independent PV-specific assay. Its held-sample behavior is reported even when adding it worsens prediction. Sample and cell-line stratified metrics, coefficients, calibration and log-loss increments are supplied; cell numbers are not biological replication.''')

    slc=read('slc6a1_group_summary.tsv');slc=slc[['group','n','positive_fraction','mean','median','q75','q90']].copy();slc.positive_fraction=slc.positive_fraction.map(pct)
    add('## 3. SLC6A1 as a separately measured axis\n\n'+table(slc)+'''

SLC6A1 is enriched but heterogeneous: roughly one in eight retained cells has no detected transcript. Its Spearman association with the consensus is 0.569, with the direct PV score 0.416, and with rank PV 0.369. After measured maturation/GABA/RNA adjustment, residual associations are approximately 0.301, 0.200 and 0.207, respectively. These support partially shared and partially separable axes rather than identity between SLC6A1 and the PV program.

**Independence correction:** SLC6A1 is absent from the explicit four-gene PV module but present in the frozen 1,000-feature mapping and sparse model, with coefficient +0.010314168. Nothing has been added to or removed from the target. A separately labeled 999-feature omit-SLC6A1 mapping sensitivity preserves the main forward/reverse interpretation. SLC6A1 RNA is relevant to GABA handling, but does not measure transporter protein or uptake activity. [Primary functional study](https://academic.oup.com/brain/article/144/8/2237/6324634).''')

    ab=read('gate_ablations.tsv');disp=ab[['markers','retained_n','target_retained_n','post_gate_target_fraction','target_recovery','fold_enrichment']].copy()
    for col in ['post_gate_target_fraction','target_recovery']:disp[col]=disp[col].map(pct)
    add('## 4. What each marker contributes\n\n'+table(disp)+'''

These are **fixed-threshold ablations** of the original triplet. They are not the previous separately optimized single/pair results. OPCML adds 13.05 composition points to GRIA2/NOTCH1; GRIA2 adds 4.96 points to OPCML/NOTCH1. NOTCH1 adds 2.92 points to GRIA2/OPCML, removing 1,347 non-targets and 643 targets; target recovery falls 54.36% → 50.81%.

GRIA2 is consistent with a neuronal receptor-development axis, OPCML with developmental cell adhesion, and low NOTCH1 with departure from a progenitor program. These are interpretations of the joint RNA evidence. None is a PV-specific fate marker. In particular, mouse MGE interneuron synapses can be dominated by GluA2-lacking AMPA receptors; GRIA2-high RNA is not a universal signature of mature PV neurons. [Matta et al.](https://www.nature.com/articles/nn.3459), [IgLON developmental expression](https://www.nature.com/articles/s41598-021-97768-5), [Notch1 neurogenesis experiments](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1000565).

The smaller NOTCH1 increment is still biologically coherent: rejected NOTCH1-high cells among GRIA2/OPCML passers are 43.52% operationally postmitotic and 34.92% cycling, versus 90.56% and 8.09% in the retained branch. Full pairwise and conditional relations with MGE, PV, SST, SLC6A1, ERBB4/LHX6 and cycle are supplied in `axis_correlations.tsv`, `surface_marker_pairwise.tsv` and `marker_conditional.tsv`.''')

    exact=' AND '.join(f'{r["gene"]} {r["op"]} {r["threshold"]!r}' for r in cfg['rules'])
    branches=read('gate_decision_branches.tsv')[['branch','n','antecedent_fraction','pv_score_mean','sst_score_mean','SLC6A1_positive_fraction','mge_score_mean','postmitotic_fraction']].copy()
    for col in ['antecedent_fraction','SLC6A1_positive_fraction','postmitotic_fraction']:branches[col]=branches[col].map(pct)
    add('## 5. Exact rule and terminal-branch biology\n\n`'+exact+'`\n\n'+'''The rule keeps GRIA2-high / OPCML-high / NOTCH1-low RNA cells. “High” and “low” refer only to these log1p(counts per 10,000) RNA thresholds. The sequential tree is for explanation; the AND rule is order-independent. Figure 05 shows every terminal branch. The full eight Boolean marker combinations are also tabulated.

'''+table(branches))
    add('''## 6. Does the target form a distinct three-dimensional region?

Figures 06 contain all three marker planes, each colored independently by antecedent, PV, SST, SLC6A1, neuronal maturation and MGE scores. All 90,631 cells and complete expression ranges enter the hexagonal-bin summaries; companion density panels include the RNA-zero mass. The thresholds accompany the distributions rather than replacing them.

The rule exploits overlapping developmental gradients and conditional combinations. GRIA2–NOTCH1 Spearman correlation is −0.619 across all cells but −0.155 within the target; GRIA2–OPCML is +0.400 overall and +0.013 within the target. Thus some global structure tracks neuronal transition while OPCML adds a comparatively separate within-target axis. A conjunction with incremental predictive value is not evidence of a discrete biological island. The modest third-marker increment and the substantial retained non-target fraction argue against describing a uniquely separated PV precursor compartment.''')
    add('''## 7. Retained non-targets: related cells, not all cutoff neighbors

Of 6,571 retained non-targets, 2,711 (41.26%) lie in the 70th–80th antecedent percentiles, and 4,561 (69.41%) in the 60th–80th. However, 968 (14.73%) are below the median. Thus categorical 58.36% composition understates broader continuous enrichment, but the entire remainder cannot be rebranded near-threshold targets.

These cells are 90.47% postmitotic, 94.81% GAD1/2-positive, 87.84% NKX2-1-positive and 69.30% LHX6-or-ERBB4-positive. SST is detected in 40.74%. Only 2.15% meet the operational progenitor proxy. Glutamate-transporter RNA is detected in 2.65%; two-gene astrocyte and oligodendroglial context flags occur in 0.02% and 1.02%, respectively; CHAT-or-SLC18A3 is detected in 5.24%. These flags overlap and are not exclusive lineage annotations: RNA nondetection cannot establish absence of a lineage. Their profile supports a broader related inhibitory developmental population with lower PV-module expression than A.''')
    add('''## 8. Lost targets: a different developmental subset

The gate loses 8,917/18,127 targets. Compared with captured targets, they have lower postmitotic fraction (81.93% versus 90.63%), more cycling evidence (13.87% versus 8.50%), more progenitor evidence (9.12% versus 0.72%), lower SLC6A1 detection (81.36% versus 92.70%) and lower SST detection (40.26% versus 60.97%). Their mean direct PV score is 0.398 versus 0.466 and antecedent score 0.804 versus 0.836.

This is not a simple loss of unrelated cells: all B cells meet the original target, and most are already operationally postmitotic. The rule captures a particular maturation and OPCML-expression phenotype within that target; adhesive function was not measured. Sample, line and CV/MW distributions accompany these comparisons; no stable fate subtype is inferred from the differences.''')

    con=read('continuous_score_summary.tsv');con=con[con.group.isin(['all','retained','removed','A','B','C','D'])][['group','n','mean','median','q75','q90','q95']]
    top=read('score_top_fraction_recovery.tsv')[['nominal_top_fraction','target_n','target_retained_n','post_gate_target_fraction','target_recovery','fold_enrichment']].copy()
    for col in ['nominal_top_fraction','post_gate_target_fraction','target_recovery']:top[col]=top[col].map(pct)
    add('## 9. Continuous enrichment beyond the top20 cutoff\n\n'+table(con)+'\n\n'+table(top)+'''

The gate captures 67.70% of top5 cells and 59.90% of top10 cells, despite retaining only 17.41% of all cells. Top-decile composition is 34.40%, or 3.44× baseline. Thus enrichment is visible across the score continuum, not solely at one binary threshold. These comparisons preserve the original target and use descriptive score quantiles with ties retained.''')
    prog=read('program_contrasts.tsv');prog=prog[prog.comparison.isin(['target_vs_non_target','retained_vs_removed'])&prog.feature.isin(['direct_pv_score','direct_sst_score','mge_score','progression_score','maturation_score','supplementary_gaba_score','supplementary_nonmodule_pv_context_score'])][['comparison','feature','mean_a','mean_b','standardized_mean_difference','rank_biserial']]
    add('## 10. PV, SST, MGE and broader interneuron programs\n\n'+table(prog)+'''

PV-associated RNA is enriched, but so are GABA development and selected SST genes. The SST composite is only modestly shifted and the MGE composite is lower, demonstrating why individual genes must accompany module scores. MAF/MAFB, LHX6 and ERBB4 support interneuron developmental context without selectively validating PV fate. Program enrichment is not independent evidence when the same genes define the operational target.

Pooled and within-sample effects are not interchangeable. The target SST composite rises pooled but is lower within five of six samples; mean SST RNA expression itself rises in four of six. SLC6A1 and maturation rise in all six, while the MGE composite falls in all six. Pooled LHX6/ERBB4 elevations are less uniform. These patterns describe sample composition and developmental heterogeneity; they do not establish a stable novel subtype or a condition-consistency pass/fail criterion.''')

    forward=pd.read_csv(OUT/'mapping/forward_group_distributions.tsv',sep='\t');f=forward[(forward.method=='archived_unbalanced')&forward.group.isin(['retained','removed'])].pivot(index='state',columns='group',values='mean').reindex(['strict_PV','dual_high','SST','unresolved']).reset_index()[['state','retained','removed']]
    for col in ['retained','removed']:f[col]=f[col].map(pct)
    f['state']=f.state.map({'strict_PV':'Strict PV-biased','dual_high':'Dual-high PV/SST','SST':'SST-biased','unresolved':'Unresolved'})
    f=f.rename(columns={'state':'DIV90 neighbor state','retained':'Retained DIV30','removed':'Removed DIV30'})
    add('## 11. Forward mapping: dual-high association, not strict-PV preference\n\n'+table(f)+'''

The strict-PV decrease is −5.40 percentage points (Hedges g −0.81); the dual-high increase is +18.03 points (g +1.54). The strict-PV decrease appears in all six samples and all three lines and persists within postmitotic cells, within-reference cells, equal-sample neighbors and the omit-SLC6A1 diagnostic space. Captured targets themselves have 10.75% strict-PV / 56.81% dual-high neighbors, versus 13.17% / 50.82% for lost targets.

Retained cells are closer to the DIV90 reference overall and less frequently above its frozen held-sample 95th-percentile distance threshold (5.14% versus 28.48%). This indicates stronger resemblance to represented later neuronal states; it does not establish preferential access to strict PV fate. Reference cells are preselected LHX6+/ERBB4+ cortical cells, and the original consensus already includes mapping/sparse PV prediction. Figures show full empirical neighbor-fraction distributions. The frozen reconstruction reproduced 4,096 checked queries' archived neighbor order and distances exactly.''')

    reverse=pd.read_csv(OUT/'mapping/reverse_query_group_distributions.tsv',sep='\t');r=reverse[(reverse.method.isin(['unbalanced','equal_sample_neighbors']))&reverse.query_group.isin(['all_DIV90_reference','strict_PV','dual_high','SST'])][['method','query_group','n_cells','mean','DIV30_gate_yield_baseline','fold_vs_DIV30_gate_yield']].copy()
    for col in ['mean','DIV30_gate_yield_baseline']:r[col]=r[col].map(pct)
    add('## 12. Reciprocal mapping: limited, nonspecific convergence\n\n'+table(r)+'''

Relative to the 17.41% all-DIV30 yield, raw strict-PV queries show modest 1.17× capture. However, capture is lower than for SST or all DIV90 reference cells and much lower than dual-high. After taking five neighbors per DIV30 source sample, strict-PV capture is 16.62% versus a matched 16.23% baseline (1.02×). The strongest-quarter strict-PV subset gives higher raw capture, but the strongest-decile subset does not increase further; convergence is nonmonotonic.

Two severe limits prevent a lineage conclusion: essentially 100% of strict-PV queries exceed the internal DIV30 held-sample 95th-percentile distance threshold, and the top 1% of DIV30 cells receive about 71.4% of strict-PV neighbor edges. This is an internal distance diagnostic, not calibrated cross-time probability validation. Inverse-degree weighting weakens apparent enrichment. Postmitotic-only recipient and unique-counterpart sensitivities are supplied. These distant, density-sensitive neighbors neither establish nor rule out later PV production. They do not independently validate the original consensus.''')

    temporal=pd.read_csv(OUT/'mapping/temporal_marker_program_effects.tsv',sep='\t');tem=temporal[temporal.context.isin(['all','LHX6_ERBB4_positive','generic_context_residual'])&temporal.marker.isin(['GRIA2','OPCML','NOTCH1'])][['stage','context','marker','hedges_g','rank_biserial']]
    add('## 13. Why DIV30 separates more readily than DIV90\n\n'+table(tem)+'''

GRIA2, OPCML and NOTCH1 show strong DIV30 target/non-target differences, but almost no separation between DIV90 strict PV-biased and SST-biased cells. GRIA2 is detected in roughly 78–81% of both later groups, OPCML in 16–21%, and NOTCH1 in only 2–4%. The early receptor/adhesion/progenitor partition is therefore largely shared or depleted across those later states.

Some DIV30 differences persist after restricting to the same LHX6+/ERBB4+ entry boundary and after generic-state adjustment, so a real early transcriptional phenotype remains. Nevertheless, comparing an all-cell DIV30 top20 developmental score against later subtype-biased labels is an easier and different classification problem. Cross-sectional data cannot distinguish causal loss of separability during maturation from population composition, technical depth and target-definition effects. DIV90's best observed strict-PV composition at ≥50% recovery remains 31.3% within the completed search, not a theoretical limit on every possible marker combination. No new DIV90 search was performed.

FGFR2/PTPRS remains the weak DIV90 purification benchmark: approximately 89% recovery, with SST-biased composition changing only 22.5% → 21.3%. Independent statistical significance is not established; this extension does not restore it as the preferred gate.''')

    perf=read('line_condition_performance.tsv');p=perf[perf.grouping.isin(['cell_line','condition','cell_line+condition'])][['grouping','group','n_cells','starting_target_fraction','post_gate_target_fraction','fold_enrichment','target_recovery','retained_n','target_retained_n']].copy()
    for col in ['starting_target_fraction','post_gate_target_fraction','target_recovery']:p[col]=p[col].map(pct)
    add('## 14. Line and CV/MW structure\n\n'+table(p)+'''

The pooled result is dominated in absolute counts by 2E, especially 2E-MW: that sample supplies 64.96% of retained targets. Starting target frequency ranges from 4.66% in 2E-CV to 42.58% in 2E-MW. Excluding 2E-MW, the **same fixed rule** still enriches 12.53% → 44.91%, with 37.81% recovery (3.58×). This is a descriptive omission analysis, not new held-sample gate selection.

These are differences in frequency and state composition, not a consistency pass/fail test. There is one cached sample per line/condition, and glucose/operator/batch distinctions are confounded. Only six source samples are present; metadata-listed samples 7–9 remain absent. The documented 9583/9853 sample alias is preserved. Original held-sample gate reselection remains unchanged.''')

    add('## 15–16. External existence, frequency and matched-state biology\n\n'+ext_summary)
    add('''## 17. Surface-protein translation after the biological identity assessment

The state is biologically coherent enough to justify a bounded feasibility review as a broader developmental enrichment. That does not turn it into a validated prospective PV sort. The following review distinguishes extracellular biology, manufacturer application claims and actual viable human-cell validation; no antibody purchase is recommended.

'''+(OUT/'mapping/surface_translation_review.md').read_text().split('\n',1)[1].replace('\n## ','\n### '))
    add('''## Evidence for and against the two interpretations

| Interpretation | Supporting evidence | Evidence against / unresolved |
| --- | --- | --- |
| Prospective PV-producing antecedent | PV-associated module enriched; continuous high-score cells preferentially retained; some marker association survives measured maturation adjustment; substantial GABAergic/interneuron context; later PV production is biologically possible. | PV module partly defines the target; SST expression is common; PVALB RNA is rare and not enriched; direct and reciprocal mapping prefer dual-high rather than strict PV; reverse mapping is distant and dominated by hubs; no tracked fate, protein endpoint or independent prospective validation. |
| Broader postmitotic neuronal state resembling later PV-associated cells | ~91% operational postmitotic; strong SLC6A1/GABA evidence; GRIA2/OPCML/NOTCH1 partitions neuronal transition, adhesion and progenitor/cycle context; A and C share postmitotic identity; substantial SST expression; dual-high mapping preference; generic covariates explain much gate membership. | Generic maturation does not explain all measured PV-module association; retained non-targets have heterogeneous scores; marker separation persists in restricted developmental comparisons. The state is therefore more structured than an arbitrary sample of every postmitotic neuron. |

The defensible working label is **GRIA2-high / OPCML-high / NOTCH1-low, SLC6A1-rich postmitotic GABAergic developmental state with PV-associated transcription**. “PV-producing antecedent” remains a hypothesis about a possible component, not the established identity of the retained population. Distinguishing the hypotheses requires prospective retained-versus-removed maturation with matched input developmental states and measured later outcomes; these cross-sectional computations do not supply that experiment.

## Methods, audit and review files

The five prior DIV30/DIV90 analysis directories, including the completed DIV90 PI figure package, were treated as read-only. Full-content SHA256 preservation covers their 2,245 files totaling 7.47 GB. New target/gate source-order validation independently reproduces all four groups and all 362,524 archived per-cell state fractions. Large searches and original target definitions remain unchanged.

Within-study tables include full gene/module distributions, effect sizes, per-sample directions, fixed clause subsets, truth patterns, conditional marker relationships, continuous quantile recovery, six held-sample model metrics, coefficients and calibration. Cell-level rank-test/FDR columns are explicitly descriptive; they are not independent biological-replicate significance. Effect sizes from different scoring panels should be interpreted alongside their gene membership and variance.

Open `index.html` for this report with figures. `DIV30_Biological_Identity_Figures.pdf` contains the within-study figure plates; mapping and external figures are in their respective `figures/` directories. The review ZIP contains compact reports, tables, figures and code; large expression arrays, per-cell archives and caches remain in this output directory. `provenance/USER_REQUEST.txt` preserves the authoritative instruction verbatim. See `BIOLOGICAL_CONTEXT.md` for primary-source interpretation and `provenance/` for independent validation and preservation records.''')
    report='\n\n'.join(sections)+'\n'
    (OUT/'REPORT.md').write_text(report)
    shutil.copy2(REPO/'python_notebooks/DIV30_IDENTITY_BIOLOGICAL_CONTEXT.md',OUT/'BIOLOGICAL_CONTEXT.md')
    shutil.copy2(__file__,OUT/'provenance'/Path(__file__).name)
    md=mistune.create_markdown(plugins=['table','strikethrough','url'])
    body=md(report)
    gallery=[]
    for base in [OUT/'figures',OUT/'mapping/figures',OUT/'external/figures']:
        capfile=base/'figure_captions.json'
        captions={r['figure']:r['caption'] for r in json.loads(capfile.read_text())} if capfile.exists() else {}
        for p in sorted(base.rglob('*.png')):
            rel=p.relative_to(OUT)
            pdf=(base/'pdf'/f'{p.stem}.pdf').relative_to(OUT) if p.parent.name=='png' else rel.with_suffix('.pdf')
            svg=(base/'svg'/f'{p.stem}.svg').relative_to(OUT) if p.parent.name=='png' else rel.with_suffix('.svg')
            caption=captions.get(p.stem,p.stem.replace('_',' '))
            gallery.append(f'<figure><a href="{html.escape(str(pdf))}"><img loading="lazy" src="{html.escape(str(rel))}" alt="{html.escape(p.stem.replace("_"," "))}"></a><figcaption>{html.escape(caption)} · <a href="{html.escape(str(pdf))}">PDF</a> · <a href="{html.escape(str(svg))}">SVG</a></figcaption></figure>')
    style='''body{font:16px/1.6 system-ui,sans-serif;color:#243647;background:#fbfcfd;max-width:1240px;margin:40px auto;padding:0 24px}h1{font-size:32px;line-height:1.2}h2{margin-top:2.4em;font-size:24px}a{color:#087e8b}table{border-collapse:collapse;width:100%;display:block;overflow-x:auto;font-size:13px;background:white}td,th{padding:8px 10px;border:1px solid #d8e0e5;text-align:left;vertical-align:top}th{background:#e9eff3}code{overflow-wrap:anywhere;background:#eef2f5}img{width:100%;height:auto}figure{margin:35px 0;padding:12px;background:white;border:1px solid #dae3e8}figcaption{font-size:14px;color:#536675}.toolbar{padding:16px;background:#e2f1ef;border-radius:8px} @media print{body{max-width:none;margin:0;font-size:10pt}h2{break-after:avoid}figure{break-inside:avoid}.toolbar{display:none}}'''
    catalogs=[]
    for base,label in [(OUT/'within/tables','Within-study tables'),(OUT/'mapping','Cross-time mapping tables'),(OUT/'external/tables','External-comparison tables')]:
        links=''.join(f'<li><a href="{html.escape(str(p.relative_to(OUT)))}">{html.escape(p.name)}</a></li>' for p in sorted(base.glob('*.tsv')))
        catalogs.append(f'<details><summary>{label}</summary><ul>{links}</ul></details>')
    page='<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>WHAT IS THE DIV30 ENRICHABLE ANTECEDENT STATE?</title><style>'+style+'</style><body><p class="toolbar"><a href="REPORT.md">Markdown report</a> · <a href="DIV30_Biological_Identity_Figures.pdf">Within-study figure book</a> · <a href="DIV30_Biological_Identity_Review.zip">Review archive</a> · <a href="#figures">All figure plates</a></p>'+body+'<h2>Supporting tables</h2>'+''.join(catalogs)+'<h2 id="figures">Figure plates</h2>'+''.join(gallery)+'</body></html>'
    (OUT/'index.html').write_text(page)
    print('Report written:',OUT/'REPORT.md',len(report.split()),'words',len(gallery),'figure plates')

if __name__=='__main__':main()
