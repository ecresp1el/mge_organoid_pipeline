"""Show female composition candidates and separate direct PCDH19 probe evidence."""
import csv
from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

ROOT=Path(__file__).resolve().parents[1]
OUT=RUNTIME_ROOT/'results/composition_review_20261031'
rows=list(csv.DictReader((OUT/'sample_proportions.tsv').open(),delimiter='\t'))
tests=list(csv.DictReader((OUT/'all_contrasts.tsv').open(),delimiter='\t'))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.linewidth':.7,'pdf.fonttype':42,'svg.fonttype':'none'})
groups=['WT_F','HET_F'];colors=['#549E89','#A678AC']
def draw(ax,values,title,effect,letter,ylabel='Cells (% of sample)'):
    for i,g in enumerate(groups):
        y=np.array(values[g]);assert len(y)==3
        ax.scatter(i+np.array([-.1,0,.1]),y,s=30,color=colors[i],edgecolors='white',linewidths=.5,zorder=3)
        ax.plot([i-.23,i+.23],[y.mean(),y.mean()],color=colors[i],lw=1.8)
    ax.set_title(title,fontsize=9,pad=25)
    ax.text(.5,1.04,'HET − WT: %+.2f pp'%effect,transform=ax.transAxes,ha='center',fontsize=8)
    ax.text(-.22,1.37,letter,transform=ax.transAxes,fontweight='bold',fontsize=12)
    ax.set_xticks([0,1]);ax.set_xticklabels(['WT female','HET female'])
    ax.set_xlim(-.5,1.5);ax.set_ylim(bottom=0);ax.set_ylabel(ylabel)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
    for side in ['top','right']:ax.spines[side].set_visible(False)
def save(fig,name):
    for suffix in ['pdf','svg','png']:fig.savefig(OUT/(name+'.'+suffix),dpi=400,facecolor='white')
    plt.close(fig)
chosen=[('full_consensus','C0006','Immature MGE interneurons\nC0006 (full dataset)'),
        ('full_consensus','C0012','MGE neurogenic progenitors\nC0012 (full dataset)'),
        ('full_consensus','C0005','Immature GABAergic neurons\nC0005 (full dataset)'),
        ('pilot_subcluster','C0001.F0005','Immature MGE interneurons\nC0001.F0005 (pilot)'),
        ('pilot_subcluster','C0002.F0011','Notch-high RG-like progenitors\nC0002.F0011 (pilot)'),
        ('pilot_cell_type','IPC','Intermediate progenitors\nGrouped IPCs (pilot)')]
fig,axes=plt.subplots(2,3,figsize=(8.5,5.8))
plotted=[]
for i,(ax,(level,category,title)) in enumerate(zip(axes.flat,chosen)):
    selected=[r for r in rows if r['level']==level and r['category']==category and r['group'] in groups]
    values={g:[100*float(r['fraction']) for r in selected if r['group']==g] for g in groups}
    t=next(t for t in tests if t['level']==level and t['category']==category and t['reference']=='WT_F' and t['case']=='HET_F')
    draw(ax,values,title,float(t['difference_pp']),'abcdef'[i]);plotted.extend(selected)
fig.subplots_adjust(left=.08,right=.99,bottom=.13,top=.79,wspace=.7,hspace=1.25)
fig.text(.08,.025,'Each dot: one sample; line: mean. n = 3/group. Exploratory candidates; no FDR-significant differences.',fontsize=7)
save(fig,'WT_female_vs_HET_female_composition_named')
with (OUT/'female_composition_plotted_values.tsv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(plotted[0]),delimiter='\t');w.writeheader();w.writerows(plotted)
source=Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/step_08_pcdh19_developmental_state_probe_detectability/step_08_per_sample_state_probe_detectability.tsv')
probe=list(csv.DictReader(source.open(),delimiter='\t'))
# Pool disjoint developmental strata within each sample, recovering all-cell detection.
aggregated={}
for r in probe:
    if r['design_group'] not in groups or r['detection_metric'] not in ['A_detected','B_detected','C_detected']:continue
    key=(r['biological_sample_id'],r['design_group'],r['detection_metric'])
    n,d=aggregated.get(key,(0,0));aggregated[key]=(n+int(r['detected_cells']),d+int(r['state_cells']))
fig,axes=plt.subplots(1,3,figsize=(6.8,2.9))
pr=[]
for i,(ax,metric) in enumerate(zip(axes,['A_detected','B_detected','C_detected'])):
    values={g:[] for g in groups}
    for (sample,g,m),(n,d) in aggregated.items():
        if m==metric:
            values[g].append(100*n/d);pr.append(dict(sample=sample,group=g,probe=metric,detected_cells=n,total_cells=d,percent=100*n/d))
    effect=np.mean(values['HET_F'])-np.mean(values['WT_F'])
    draw(ax,values,'PCDH19 probe '+metric[0],effect,'abc'[i],'Probe-positive cells (%)')
fig.subplots_adjust(left=.09,right=.99,bottom=.25,top=.67,wspace=.65)
fig.text(.09,.045,'All cells per sample; n = 3/group. Probe detection is separate from cell type abundance.',fontsize=7)
save(fig,'WT_female_vs_HET_female_PCDH19_probes')
with (OUT/'female_probe_plotted_values.tsv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(pr[0]),delimiter='\t');w.writeheader();w.writerows(pr)
print('Saved named female composition and all-cell probe figures with source values.')
for metric in ['A_detected','B_detected','C_detected']:
    print(metric,[(g,round(np.mean([r['percent'] for r in pr if r['probe']==metric and r['group']==g]),3)) for g in groups])
