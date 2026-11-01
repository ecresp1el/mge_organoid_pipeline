"""Integrate the complete exploratory screen, including null and selection-control results."""
import csv
import gzip
import json
import itertools
import hashlib
from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
OUT=RUNTIME_ROOT/'results/mosaic_deep_dive_20261031'
def read(p):
    with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:return list(csv.DictReader(f,delimiter='\t'))
def write(name,rows):
    if not rows:return
    with (OUT/name).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def mean(x):return float(np.mean(x))
def fnum(v):
    try:return float(v)
    except (ValueError,TypeError):return np.nan
plan=json.loads((OUT/'EXPLORATORY_PLAN.json').read_text())
seal=json.loads((OUT/'PREPARATION_COMPLETE.json').read_text())
summary=[];tops=[];programs=[]
for st in plan['strata']:
    d=OUT/('DE_'+st)
    assert (d/'COMPLETE.txt').exists(),st+' is incomplete'
    summary.extend(read(d/'summary.tsv'))
    if (d/'top30_all_definitions.tsv').exists():tops.extend(read(d/'top30_all_definitions.tsv'))
    if (d/'program_effects.tsv').exists():programs.extend(read(d/'program_effects.tsv'))
assert not any(r['status']=='FAILED' for r in summary)
write('all_definition_DE_summary.tsv',summary);write('top30_all_definitions.tsv',tops);write('all_program_effects.tsv',programs)
comp=read(OUT/'sample_gate_stratum_composition.tsv')
scores=read(OUT/'sample_gate_program_scores.tsv')
metadata=read(OUT/'pseudobulk_metadata.tsv')
females=['15662-JZ-'+str(i) for i in range(4,10)]
cs={(r['sample'],r['gate'],r['stratum']):r for r in comp}
ss={(r['sample'],r['gate'],r['program']):r for r in scores}
qc=[]
for s in females:
    for gate in plan['gates']:
        row=cs.get((s,gate,'all_cells'))
        if not row:continue
        qc.append(dict(sample=s,design_group='WT_F' if s in females[:3] else 'HET_F',gate=gate,cells=int(row['cells']),
          mean_UMI=float(ss[s,gate,'total_counts']['mean']),mean_genes=float(ss[s,gate,'n_genes']['mean']),
          vascular_percent=100*float(cs.get((s,gate,'vascular_marked'),{'fraction':0})['fraction']),
          mixed_percent=100*float(cs.get((s,gate,'mixed_profile'),{'fraction':0})['fraction']),
          neural_vascular_low_percent=100*float(cs.get((s,gate,'neural_vascular_low'),{'fraction':0})['fraction'])))
write('female_gate_identity_and_QC.tsv',qc)

# Continuous marker programs are summarized at the sample level, never using cells as replicates.
score_contrasts=[]
for gate in plan['gates']:
    for p in plan['programs']:
        if not all((s,gate,p) in ss for s in females[3:]):continue
        a=[float(ss[s,'all',p]['mean']) for s in females[:3]]
        b=[float(ss[s,gate,p]['mean']) for s in females[3:]]
        diff=mean(b)-mean(a);vals=a+b
        null=[mean([vals[i] for i in ix])-mean([vals[i] for i in range(6) if i not in ix]) for ix in itertools.combinations(range(6),3)]
        score_contrasts.append(dict(gate=gate,program=p,WT_mean=mean(a),HET_mean=mean(b),difference=diff,
              exact_p=sum(abs(v)>=abs(diff)-1e-12 for v in null)/20,complete_sample_separation=max(a)<min(b) or max(b)<min(a)))
write('sample_program_score_contrasts.tsv',score_contrasts)

# Summarize repeated gene effects across related definitions; recurrence is sensitivity evidence, not replication.
recurrence=defaultdict(list)
for st in ['neural_vascular_low','C0001','C0006','C0001_neural_vascular_low','neural_depth_trimmed']:
    for contrast in ['HET_Aneg_BC1_vs_WTall','HET_Aneg_BC2_vs_WTall','HET_Aneg_BC3_vs_WTall','HET_KOtail10_vs_WTall','HET_KOtail20_vs_WTall',
                     'matched_A_positive','matched_Aneg_BC1','within_HET_Aneg_BC1_vs_Apositive','within_HET_Aneg_BC2_vs_Apositive','gate_interaction_BC1_vs_Apositive']:
        path=OUT/('DE_'+st)/(contrast+'.tsv.gz')
        if not path.exists():continue
        for r in read(path):
            if r['gene_symbol']=='Pcdh19':continue
            if float(r['PValue'])<.05 and abs(float(r['logFC']))>=.5 and r['direction_concordant']=='TRUE':
                recurrence[r['gene_id'],r['gene_symbol']].append((st,contrast,float(r['logFC']),float(r['PValue']),float(r['FDR'])))
ranked=[]
for (gene,symbol),evidence in recurrence.items():
    ranked.append(dict(gene_id=gene,gene_symbol=symbol,definitions=len(evidence),strata=len(set(e[0] for e in evidence)),
                       median_logFC=mean([e[2] for e in evidence]),same_sign=all(e[2]>0 for e in evidence) or all(e[2]<0 for e in evidence),
                       best_P=min(e[3] for e in evidence),best_FDR=min(e[4] for e in evidence),
                       contexts=';'.join(e[0]+':'+e[1] for e in evidence)))
ranked.sort(key=lambda r:(not r['same_sign'],-r['strata'],-r['definitions'],r['best_P']))
write('recurring_exploratory_gene_candidates.tsv',ranked)

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.linewidth':.7,'pdf.fonttype':42,'svg.fonttype':'none'})
plotgroups=[('WT_F','all','WT all'),('HET_F','all','HET all'),('HET_F','A_positive','HET A+'),
            ('HET_F','Aneg_BC1','HET BC≥1'),('HET_F','Aneg_BC2','HET BC≥2'),('HET_F','Aneg_BC3','HET BC≥3'),('HET_F','frozen_KO','HET strict KO')]
fig,axes=plt.subplots(1,3,figsize=(10,3.8))
for ax,field,label in zip(axes,['cells','mean_UMI','vascular_percent'],['Selected cells per sample','Mean RNA counts per cell','Vascular-marker-positive cells (%)']):
    for i,(group,gate,title) in enumerate(plotgroups):
        vals=[r[field] for r in qc if r['design_group']==group and r['gate']==gate]
        ax.scatter(i+np.array([-.1,0,.1]),vals,s=24,c='#549E89' if group=='WT_F' else '#A678AC')
        ax.plot([i-.24,i+.24],[mean(vals),mean(vals)],color='black',lw=1)
    ax.set_xticks(range(len(plotgroups)));ax.set_xticklabels([x[2] for x in plotgroups],rotation=55,ha='right',fontsize=7)
    ax.set_ylabel(label);ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
    if field=='cells':ax.set_yscale('log')
    else:ax.set_ylim(bottom=0)
fig.subplots_adjust(left=.075,right=.99,bottom=.35,top=.9,wspace=.55)
fig.text(.075,.03,'Each dot is one sample. BC gates require A absent. Marker positivity indicates expression, not a definitive cell identity.',fontsize=7)
for ext in ['pdf','png','svg']:fig.savefig(OUT/('gate_selection_diagnostics.'+ext),dpi=350,bbox_inches='tight',facecolor='white')
plt.close(fig)

comparisons=['HET_all_vs_WTall','HET_A_positive_vs_WTall','HET_Aneg_BC1_vs_WTall','HET_Aneg_BC2_vs_WTall','HET_Aneg_BC3_vs_WTall',
             'HET_frozen_KO_vs_WTall','matched_A_positive','matched_Aneg_BC1','within_HET_Aneg_BC1_vs_Apositive','gate_interaction_BC1_vs_Apositive']
matrix=np.full((len(comparisons),len(plan['strata'])),np.nan)
index={(r['stratum'],r['contrast']):r for r in summary}
for i,c in enumerate(comparisons):
    for j,st in enumerate(plan['strata']):
        r=index.get((st,c))
        if r and r['status']=='COMPLETE':matrix[i,j]=float(r['nominal01_concordant_FC05'])
fig,ax=plt.subplots(figsize=(11,6));im=ax.imshow(np.ma.masked_invalid(np.log1p(matrix)),aspect='auto',cmap='viridis')
ax.set_yticks(range(len(comparisons)));ax.set_yticklabels(comparisons,fontsize=7)
ax.set_xticks(range(len(plan['strata'])));ax.set_xticklabels(plan['strata'],rotation=55,ha='right',fontsize=7)
for i in range(len(comparisons)):
    for j in range(len(plan['strata'])):
        if np.isfinite(matrix[i,j]):ax.text(j,i,str(int(matrix[i,j])),ha='center',va='center',fontsize=7,color='white' if matrix[i,j]<25 else 'black')
ax.set_title('Exploratory signals by definition: P < 0.01, |log2FC| ≥ 0.5, consistent sample direction',fontsize=10)
fig.colorbar(im,ax=ax,label='log(1 + candidate genes)');fig.subplots_adjust(left=.33,right=.94,bottom=.28,top=.92)
for ext in ['pdf','png','svg']:fig.savefig(OUT/('definition_signal_map.'+ext),dpi=350,bbox_inches='tight',facecolor='white')
plt.close(fig)

lines=['# Exploratory PCDH19 mosaic-definition deep dive','',
       'Exploratory, hypothesis-generating results. This screen deliberately includes broad probe gates and nominal gene-level signals; absence of FDR significance does not remove a result from the tables. All completed and skipped definitions are retained. Comparisons selected after this screen are not independent confirmation.',
       '', '## What was examined','',
       '- 14 probe definitions, including A-positive WT-like cells; A-negative cells with B+C ≥1, ≥2 or ≥3; frozen control-derived calls; and 5%, 10%, 20% score tails with cutoff ties retained.',
       '- 11 cell strata: full population, neural/vascular-marker strata, five consensus clusters, vascular-low C0001, and a depth-trimmed neural stratum.',
       '- WT females versus HET subsets; same-gate WT/HET controls; paired within-HET comparisons; and genotype-by-gate interactions. Interaction contrasts ask whether the HET gate contrast differs from the average same-gate contrast in WT females.',
       '- 25 independently specified marker programs. Raw counts were aggregated by sample; three registered samples per genotype are used as replicates. Independent embryo/litter structure and preparation covariates remain unverified.',
       '', '## Gate identity and RNA-content diagnostics','',
       '| Group / gate | Mean cells/sample | Mean RNA/cell | Vascular-marker positive % | Mixed neural/vascular % |','|---|---:|---:|---:|---:|']
for group,gate,label in plotgroups:
    r=[r for r in qc if r['design_group']==group and r['gate']==gate]
    lines.append('| %s | %.1f | %.0f | %.2f | %.2f |'%(label,mean([x['cells'] for x in r]),mean([x['mean_UMI'] for x in r]),mean([x['vascular_percent'] for x in r]),mean([x['mixed_percent'] for x in r])))
lines += ['','## Expression screen','', '| Stratum | Comparison | FDR <0.05 genes | FDR <0.10 genes | Nominal consistent candidates | Minimum cells/library |', '|---|---|---:|---:|---:|---:|']
order=sorted([r for r in summary if r['status']=='COMPLETE'],key=lambda r:(-int(r['nominal01_concordant_FC05']),-int(r['FDR05'])))
for r in order[:40]:lines.append('| {stratum} | {contrast} | {FDR05} | {FDR10} | {nominal01_concordant_FC05} | {min_cells} |'.format(**r))
lines += ['','Nominal candidate definition: P <0.01, |log2FC| ≥0.5 and concordant direction across samples. For unpaired contrasts, concordance requires all three case sample values to lie on the same side of all three reference sample values. These are exploratory screens, not multiplicity-corrected discoveries. Pcdh19 is excluded from candidate counts because its probe evidence defines the groups.',
          '', '## Recurring candidates in neural or cluster-restricted analyses','', '| Gene | Definitions | Cell strata | Mean log2FC | Best FDR | Same sign |','|---|---:|---:|---:|---:|---|']
for r in ranked[:30]:lines.append('| {gene_symbol} | {definitions} | {strata} | {median_logFC:+.2f} | {best_FDR:.3g} | {same_sign} |'.format(**r))
lines += ['','Recurring-candidate screen uses nominal P <0.05, |log2FC| ≥0.5 and concordant direction across selected neural/cluster-restricted comparisons. Related gates are correlated; recurrence measures sensitivity across definitions, not independent replication.',
          '', '## Interpretation limits','',
          'A-positive and A-negative B/C-positive groups are probe-evidence groups. Broad KO-enriched gates and percentile tails do not inherit the frozen classifier precision. Marker positivity is not an anatomical identity, and mixed expression does not prove a doublet. Neural/vascular exclusion and depth trimming are sensitivity definitions, not approved changes to the primary cell population. Genotype/gate differences may reflect cell mixture, capture depth, assay evidence or genuine expression changes. No causal or allele-specific DNA inference is made.',
          '', '## Files','',
          '- gate_selection_diagnostics.pdf: coverage, RNA content and vascular-marker positivity.',
          '- definition_signal_map.pdf: where exploratory gene differences appear across definitions.',
          '- all_definition_DE_summary.tsv and top30_all_definitions.tsv: complete screen and top genes.',
          '- DE_<stratum>/*.tsv.gz: full gene results for each tested definition.',
          '- sample_program_score_contrasts.tsv and all_program_effects.tsv: marker-program context.',
          '- per_cell_marker_diagnostics.h5, raw_pseudobulk_counts.npz and per-stratum counts: reusable underlying evidence.']
(OUT/'REPORT.md').write_text('\n'.join(lines)+'\n')
checks=dict(cells=seal['cells'],genes=seal['genes'],definitions=len(summary),completed=sum(r['status']=='COMPLETE' for r in summary),skipped=sum(r['status']=='SKIPPED' for r in summary),failed=0,all_strata_complete=True)
(OUT/'REVIEW_COMPLETE.json').write_text(json.dumps(checks,indent=2))
print(json.dumps(checks));print('\n'.join(lines[:27]));print('TOP_DEFINITIONS')
for r in order[:15]:print(r)
print('TOP_RECURRING_GENES')
for r in ranked[:15]:print(r)
