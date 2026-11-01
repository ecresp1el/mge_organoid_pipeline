"""Publication-style display of the largest consistent sex-matched shifts."""
import csv
from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

ROOT = Path(__file__).resolve().parents[1]
OUT = RUNTIME_ROOT / 'results/composition_review_20261031'
ROWS = list(csv.DictReader((OUT/'sample_proportions.tsv').open(), delimiter='\t'))
TESTS = list(csv.DictReader((OUT/'all_contrasts.tsv').open(), delimiter='\t'))
ranked = [r for r in TESTS if r['level']=='full_consensus' and r['reference']=='WT_M'
          and r['case']=='KO_M' and r['complete_sample_separation']=='True'
          and r['all_leave_out_directions_agree']=='True']
ranked.sort(key=lambda r:-abs(float(r['difference_pp'])))
selected = ranked[:4]
assert [r['category'] for r in selected] == ['C0001','C0002','C0003','C0005']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.linewidth':.7,
                     'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none',
                     'xtick.major.width':.7,'ytick.major.width':.7,
                     'xtick.major.size':3,'ytick.major.size':3})
colors={'WT_M':'#337EAC','KO_M':'#D36B45','WT_F':'#549E89','HET_F':'#A678AC'}
labels={'WT_M':'WT male','KO_M':'KO male','WT_F':'WT female','HET_F':'HET female'}
names={'C0001':'MGE neurogenic\nprogenitors', 'C0002':'Radial glia-like\nprogenitors',
       'C0003':'MGE IPCs', 'C0005':'Immature\nGABAergic neurons',
       'C0006':'Immature MGE\ninterneurons','C0012':'MGE neurogenic\nprogenitors'}
export=[]
for name,groups,shape,size in [
    ('largest_consistent_differences',['WT_M','KO_M'],(1,4),(7.2,2.65)),
    ('largest_differences_all_groups',['WT_M','WT_F','HET_F','KO_M'],(2,2),(6.8,5.6))]:
    fig,axes=plt.subplots(*shape,figsize=size,squeeze=False)
    for idx,(ax,test) in enumerate(zip(axes.flat,selected)):
        cluster=test['category']
        for x,g in enumerate(groups):
            rows=[r for r in ROWS if r['level']=='full_consensus' and r['category']==cluster and r['group']==g]
            y=np.array([100*float(r['fraction']) for r in rows])
            ax.scatter(x+np.array([-.11,0,.11]),y,s=30,c=colors[g],edgecolors='white',linewidths=.5,zorder=3)
            ax.plot([x-.23,x+.23],[y.mean(),y.mean()],color=colors[g],lw=1.8,zorder=4)
            if name=='largest_consistent_differences':
                export.extend(rows)
        ax.set_title(names[cluster]+'\n'+cluster,fontsize=9,pad=23)
        ax.text(.5,1.04,'KO \u2212 WT: %+.2f pp'%float(test['difference_pp']),
                transform=ax.transAxes,ha='center',fontsize=8,color='#444444')
        ax.text(-.22,1.37,'abcd'[idx],transform=ax.transAxes,fontweight='bold',fontsize=12)
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels([labels[g] for g in groups],rotation=25 if len(groups)>2 else 0,ha='right' if len(groups)>2 else 'center')
        ax.set_xlim(-.5,len(groups)-.5)
        ax.set_ylim(bottom=0)
        ax.set_ylabel('Cells (% of sample)')
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    fig.subplots_adjust(left=.09,right=.985,bottom=.27 if shape[0]==1 else .16,top=.65 if shape[0]==1 else .80,wspace=.65,hspace=1.25)
    fig.text(.09,.035 if shape[0]==1 else .025,'Each dot: one sample; horizontal line: mean. n = 3 per group.',fontsize=7,color='#444444')
    for suffix in ['pdf','svg','png']:
        fig.savefig(OUT/(name+'.'+suffix),dpi=400,facecolor='white',bbox_inches='tight')
    plt.close(fig)
with (OUT/'largest_differences_plotted_values.tsv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(export[0]),delimiter='\t');w.writeheader();w.writerows(export)
(OUT/'FIGURE_LEGEND.txt').write_text(
    'Cell composition differences between WT and Pcdh19 KO males. Panels a-d show the four largest absolute mean shifts among full-data consensus clusters with strict separation of all three WT and all three KO sample proportions and consistent direction under leave-one-sample and leave-one-per-group analyses. Proportions are calculated within each sample using all retained cells as the denominator. Dots represent registered samples; horizontal lines indicate unweighted group means (n=3/group). pp, percentage points. Cluster IDs refer to the 34-cluster full-data consensus and do not denote the annotated pilot clusters. The companion figure includes WT and HET females for context. Panels were selected after exploratory screening; exact two-sided permutation P=0.10 for each displayed WT-male/KO-male contrast, with no FDR-significant finding in the full screen. Replicate independence and preparation covariates remain to be confirmed.\n')
print('Saved two figures as PDF, SVG and 400-dpi PNG, plus plotted values and legend.')
