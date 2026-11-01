"""Requested asymmetric comparison: unfiltered WT females versus gated HET subsets."""
import csv
from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import textwrap

ROOT=Path(__file__).resolve().parents[1]
OUT=RUNTIME_ROOT/'results/unfiltered_WT_vs_gated_HET_20261031'
OUT.mkdir(exist_ok=True,parents=True)
WT=RUNTIME_ROOT/'results/composition_review_20261031/sample_proportions.tsv'
HET=RUNTIME_ROOT/'results/probe_gated_composition_20261031/gated_sample_proportions.tsv'
NAMES=RUNTIME_ROOT/'results/composition_review_20261031/full_consensus_working_names.tsv'
def read(p):return list(csv.DictReader(p.open(),delimiter='\t'))
names={r['cluster']:r['working_name'] for r in read(NAMES)}
data=[]
for r in read(WT):
    if r['level']=='full_consensus' and r['group']=='WT_F':
        data.append(dict(group='WT_F_all',sample=r['sample'],cluster=r['category'],percent=100*float(r['fraction']),count=int(r['count']),denominator=int(r['total'])))
for r in read(HET):
    if r['sample'] not in ['15662-JZ-7','15662-JZ-8','15662-JZ-9']:continue
    if r['gate'] not in ['all','WT_like','KO_like','any_probe_positive']:continue
    data.append(dict(group='HET_F_'+r['gate'],sample=r['sample'],cluster=r['cluster'],percent=float(r['percent']),count=int(r['cluster_cells']),denominator=int(r['gate_cells'])))
def write(name,rows):
    with (OUT/name).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
write('sample_proportions.tsv',data)
summaries=[]
for c in sorted(set(r['cluster'] for r in data)):
    wt=[r['percent'] for r in data if r['group']=='WT_F_all' and r['cluster']==c]
    assert len(wt)==3
    for g in ['HET_F_all','HET_F_WT_like','HET_F_KO_like','HET_F_any_probe_positive']:
        vals=[r['percent'] for r in data if r['group']==g and r['cluster']==c]
        assert len(vals)==3
        effect=np.mean(vals)-np.mean(wt)
        leave=[np.mean(vals)-np.mean(wt[:i]+wt[i+1:]) for i in range(3)]
        leave += [np.mean(vals[:i]+vals[i+1:])-np.mean(wt) for i in range(3)]
        summaries.append(dict(cluster=c,working_name=names.get(c,'Unannotated'),comparison=g+'_minus_WT_F_all',WT_F_mean_pct=np.mean(wt),HET_subset_mean_pct=np.mean(vals),difference_pp=effect,
                              complete_sample_separation=max(wt)<min(vals) or max(vals)<min(wt),leave_one_sample_direction_stable=bool(effect and all(v*effect>0 for v in leave))))
summaries.sort(key=lambda r:(r['comparison'],-abs(r['difference_pp'])))
write('ranked_comparisons.tsv',summaries)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.linewidth':.7,'pdf.fonttype':42,'svg.fonttype':'none'})
chosen=['C0001','C0006','C0007','C0011']
groups=['WT_F_all','HET_F_all','HET_F_WT_like','HET_F_KO_like']
labels=['WT\nall cells','HET\nall cells','HET\nWT-like','HET\nKO-like']
colors=['#549E89','#999999','#337EAC','#A678AC']
fig,axes=plt.subplots(1,4,figsize=(10,3.6))
for i,(ax,c) in enumerate(zip(axes,chosen)):
    for j,g in enumerate(groups):
        vals=[r['percent'] for r in data if r['group']==g and r['cluster']==c]
        ax.scatter(j+np.array([-.12,0,.12]),vals,s=28,color=colors[j],edgecolors='white',linewidths=.5,zorder=3)
        ax.plot([j-.25,j+.25],[np.mean(vals),np.mean(vals)],color=colors[j],lw=1.7)
    ax.set_title('\n'.join(textwrap.wrap(names.get(c,c),20))+'\n'+c,fontsize=9,pad=12)
    ax.text(-.2,1.39,'abcd'[i],transform=ax.transAxes,fontweight='bold',fontsize=12)
    ax.set_xticks(range(4));ax.set_xticklabels(labels,fontsize=7)
    ax.set_xlim(-.5,3.5);ax.set_ylim(bottom=0);ax.set_ylabel('Cells (% of sample or subset)')
    for side in ['top','right']:ax.spines[side].set_visible(False)
fig.subplots_adjust(left=.065,right=.99,bottom=.32,top=.68,wspace=.65)
fig.text(.065,.045,'Female samples; WT remains unfiltered. Each dot: one sample; line: mean. HET KO-like subset: 26, 10 and 103 cells/sample.',fontsize=7)
fig.text(.065,.005,'Different selection rules change the measured population; these contrasts do not estimate unbiased genotype-specific tissue abundance.',fontsize=7)
for ext in ['pdf','svg','png']:fig.savefig(OUT/('unfiltered_WT_vs_gated_HET.'+ext),dpi=400,bbox_inches='tight',facecolor='white')
plt.close(fig)
report=['# Unfiltered WT females versus probe-gated HET females','',
        'WT-female cells are not probe-filtered. HET gates reuse the frozen classifier unchanged. The unfiltered HET baseline is retained for context. Values are equal-weight means of three sample-specific percentages; WT denominators are all retained cells, while gated HET denominators are cells within that sample and gate.',
        '', 'The HET WT-like subset contains 3,947 retained cells; the KO-like subset contains 139 (26, 10, 103). WT-like and KO-like are inferred probe-evidence states, not definitive genotypes. Control precision was 99.684% and 76.316%, respectively.',
        '', 'Asymmetric gating changes the target population. Higher probe detection in a cluster can enrich both HET subsets relative to unfiltered WT. The WT negative control also showed gate-dependent shifts. These comparisons are descriptive enrichment screens, not proof of genotype-specific abundance changes; no inferential P values are assigned.',
        '', '| Cluster | Working name | WT all % | HET all % | HET WT-like % | HET KO-like % |','|---|---|---:|---:|---:|---:|']
for c in chosen:
    vals=[np.mean([r['percent'] for r in data if r['group']==g and r['cluster']==c]) for g in groups]
    report.append('| '+c+' | '+names[c]+' | '+' | '.join('%.2f'%v for v in vals)+' |')
(OUT/'REPORT.md').write_text('\n'.join(report)+'\n')
print('\n'.join(report[-6:]))
