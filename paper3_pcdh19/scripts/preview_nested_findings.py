"""Standalone proposal; does not modify the approved multipanel figure."""
from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import csv,json,shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
SRC=RUNTIME_ROOT/'results/nested_proliferation_20261101'
OUT=RUNTIME_ROOT/'results/nested_figure_preview_20261101';OUT.mkdir(parents=True,exist_ok=True)
rows=list(csv.DictReader((SRC/'group_summary.tsv').open(),delimiter='\t'))
def get(mode,parent,endpoint,group):return next(r for r in rows if r['mode']==mode and r['parent']==parent and r['endpoint']==endpoint and r['group']==group)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
fig=plt.figure(figsize=(9.2,3.4));ink='#26292E'
fig.text(.025,.93,'A',fontsize=16,weight='bold');fig.text(.065,.93,'Female IPCs · cycling-marker detection',fontsize=10,weight='bold')
fig.text(.565,.93,'B',fontsize=16,weight='bold');fig.text(.605,.93,'Male progenitor composition',fontsize=10,weight='bold')
ax=fig.add_axes([.09,.27,.30,.53]);genes=['Pcna','Mcm2','Birc5'];groups=['WT_F','HET_WT','HET_KO'];v=np.array([[float(get('CE_cell_groups','IPCs',g+' detected',gr)['mean']) for gr in groups] for g in genes])
im=ax.imshow(v,cmap='Greys',vmin=0,vmax=100,aspect='auto')
for i in range(3):
 for j in range(3):ax.text(j,i,'%.1f%%'%v[i,j],ha='center',va='center',fontsize=11,color='white' if v[i,j]>52 else ink)
ax.set_xticks([0,1,2]);ax.set_xticklabels(['WT\nfemale','HET\nWT cells','HET\nKO cells']);ax.set_yticks([0,1,2]);ax.set_yticklabels(genes,fontstyle='italic');ax.tick_params(length=0,pad=6)
ax.set_xticks(np.arange(-.5,3,.1) if False else [-.5,.5,1.5,2.5],minor=True);ax.set_yticks([-.5,.5,1.5,2.5],minor=True);ax.grid(which='minor',color='white',lw=2);ax.tick_params(which='minor',length=0)
for sp in ax.spines.values():sp.set_visible(False)
delta=fig.add_axes([.405,.27,.12,.53]);delta.set_ylim(2.5,-.5);delta.set_xlim(0,1);delta.axis('off');delta.text(.5,1.06,'HET KO − HET WT',transform=delta.transAxes,ha='center',fontsize=8)
for i in range(3):delta.text(.5,i,'%+.1f pp'%(v[i,2]-v[i,1]),ha='center',va='center',weight='bold',fontsize=10,color=ink)
cax=fig.add_axes([.09,.115,.30,.025]);cb=fig.colorbar(im,cax=cax,orientation='horizontal',ticks=[0,50,100]);cb.ax.tick_params(length=0,labelsize=8,pad=2);cb.outline.set_visible(False);fig.text(.09,.025,'IPCs with detected transcript (%)',fontsize=8)
ax=fig.add_axes([.65,.27,.32,.53]);cats=['RG','IPCs','Neurogenic\nprogenitors'];endpoints=['Share: RG','Share: IPCs','Share: Neurogenic progenitors'];cols=['#437C90','#D18A43'];values=[]
for i,ep in enumerate(endpoints):
 rr=[get('genotype_sex','All progenitors / cycling neuroblasts',ep,g) for g in ['WT_M','KO_M']];means=[float(r['mean']) for r in rr];values.append(means)
 ax.plot(means,[i,i],color='#B8BCC1',lw=1.5,zorder=1)
 for j,r in enumerate(rr):
  x=float(r['mean']);ax.scatter(x,i,s=55,color=cols[j],edgecolors='white',linewidths=.6,zorder=3)
  ax.text(x,i+(.20 if j==0 else -.20),'%.1f%%'%x,ha='center',va='center',fontsize=8,color=cols[j])
ax.set_yticks(range(3));ax.set_yticklabels(cats);ax.set_ylim(2.55,-.55);ax.set_xlim(0,70);ax.set_xticks([0,20,40,60]);ax.set_xlabel('Cells within progenitor pool (%)',fontsize=9);ax.tick_params(axis='y',length=0,pad=6);ax.spines['left'].set_visible(False);ax.spines['bottom'].set_color('#969CA3');ax.tick_params(axis='x',color='#969CA3')
ax.legend(handles=[Line2D([0],[0],marker='o',ls='',color=cols[i],label=n,markersize=6) for i,n in enumerate(['WT males','KO males'])],loc='upper center',bbox_to_anchor=(.5,1.20),ncol=2,frameon=False,handletextpad=.2,columnspacing=1,fontsize=8)
for ext in ['png','pdf','svg']:fig.savefig(str(OUT/('nested_findings_preview.'+ext)),dpi=350,facecolor='white')
(OUT/'PREVIEW_METHODS.json').write_text(json.dumps({'source_summary':str(SRC/'group_summary.tsv'),'female_values_percent':v.tolist(),'male_values_percent':values,'female_denominator':'all IPCs within each sample/probe-defined group','male_denominator':'all broad progenitor/cycling-neuroblast cells within each sample, no probe gate','aggregation':'equal-weight mean across 3 samples; no sample IDs displayed','heatmap_scale':'absolute marker-detection percentage, 0-100; no row z-scoring','pp':'percentage points, HET KO minus HET WT','status':'standalone preview; main figure unchanged','inference':'selected exploratory findings, not FDR-confirmed; marker detection not division rate'},indent=2))
shutil.copy2(__file__,str(OUT/Path(__file__).name))
print('PREVIEW_COMPLETE',OUT)
