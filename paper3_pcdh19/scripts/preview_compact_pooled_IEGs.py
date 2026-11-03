"""Compact descriptive pooled IEG ratios with an explicit pooling-sensitivity range."""
from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import csv,json,shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=RUNTIME_ROOT/'results/compact_pooled_IEG_preview_20261101';OUT.mkdir(parents=True,exist_ok=True)
SRC=RUNTIME_ROOT/'results/pooled_PCDH19_IEGs_20261101/all_pooling_sensitivities.tsv'
r=list(csv.DictReader(SRC.open(),delimiter='\t'));genes=['Egr1','Npas4','Nr4a1'];values=[]
for gene in genes:
 rows=[x for x in r if x['gene']==gene];primary=next(x for x in rows if x['variant']=='WT_all_positive');variants=[float(x['cell_weighted_KO_over_WT']) for x in rows]
 values.append(dict(gene=gene,primary_percent_change=100*(float(primary['cell_weighted_KO_over_WT'])-1),sensitivity_min=100*(min(variants)-1),sensitivity_max=100*(max(variants)-1)))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
fig=plt.figure(figsize=(4.8,2.8));ax=fig.add_axes([.18,.25,.75,.53]);col='#A66F50'
ax.axvline(0,color='#91979E',lw=.9,ls=(0,(3,3)),zorder=0)
for i,row in enumerate(values):
 x=row['primary_percent_change'];ax.plot([row['sensitivity_min'],row['sensitivity_max']],[i,i],lw=2.2,color='#C3C7CC',solid_capstyle='round',zorder=1);ax.scatter(x,i,s=65,color=col,zorder=3,edgecolors='white',linewidths=.6)
 ax.text(x,i-.23,'%+.0f%%'%x,ha='center',va='center',fontsize=10,color=col,weight='bold')
ax.set_yticks(range(3));ax.set_yticklabels(genes,fontstyle='italic');ax.tick_params(axis='y',length=0,pad=8);ax.set_ylim(2.50,-.55);ax.set_xlim(-20,45);ax.set_xticks([-20,0,20,40]);ax.set_xticklabels(['−20','0','+20','+40']);ax.set_xlabel('Pooled KO vs WT expression change (%)',fontsize=9);ax.spines['left'].set_visible(False);ax.spines['bottom'].set_color('#91979E');ax.tick_params(axis='x',length=3,color='#91979E')
fig.text(.06,.92,'Pooled IEG expression',fontsize=12,weight='bold');fig.text(.06,.835,'Same immature inhibitory neuronal cluster',fontsize=9,color='#5E6670')
fig.text(.18,.035,'Dot: requested pool   Line: pooling sensitivity',fontsize=8,color='#5E6670')
for ext in ['png','pdf','svg']:fig.savefig(str(OUT/('compact_pooled_IEGs.'+ext)),dpi=350,facecolor='white')
(OUT/'METHODS.json').write_text(json.dumps({'source':str(SRC),'values':values,'dot':'primary requested cell-weighted mean linear normalized expression ratio: pooled KO males plus HET KO vs Pcdh19-positive WT male/female plus HET WT','line':'min/max cell-weighted ratio over four previously specified WT-pool definitions; not confidence interval','zero':'equal pooled KO and WT expression','same_cluster':'C0006','claim':'descriptive pooled expression increases, not a causal PCDH19-loss estimate or consistent embryo-level activation','prior_sensitivities':'mean log-normalized Egr1 essentially unchanged, Nr4a1 slightly lower; no change to prior results','main_figure_changed':False},indent=2));shutil.copy2(__file__,str(OUT/Path(__file__).name));print('COMPACT_PREVIEW_COMPLETE')
