"""Same neuronal subtype, embryo-paired IEG diagnostic; no supportive panel unless consistent."""
from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import csv,gzip,json,shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
OUT=RUNTIME_ROOT/'results/paired_HET_IEGs_20261101'
rows=list(csv.DictReader((OUT/'paired_IEG_expression.tsv').open(),delimiter='\t'))
SRC=RUNTIME_ROOT/'results/mosaic_deep_dive_20261031/DE_C0006/within_HET_Aneg_BC1_vs_Apositive.tsv.gz'
with gzip.open(str(SRC),'rt') as f:de={r['gene_symbol']:r for r in csv.DictReader(f,delimiter='\t')}
genes=['Egr1','Npas4','Nr4a1'];samples=['15662-JZ-'+str(i) for i in [7,8,9]];colors=['#557C70','#AA814B','#8F7898'];summary=[];paired=[]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
fig,axs=plt.subplots(1,3,figsize=(8.2,3.4))
for ax,gene in zip(axs,genes):
 changes=[]
 for j,s in enumerate(samples):
  rr=[next(r for r in rows if r['gene']==gene and r['sample']==s and r['group']==group) for group in ['WT cells','KO cells']];val=np.array([float(r['CPM']) for r in rr]);ld=float(rr[1]['log2CPM'])-float(rr[0]['log2CPM']);changes.append(ld)
  ax.plot([0,1],val,color=colors[j],lw=1.3,alpha=.85,zorder=1);ax.scatter([0,1],val,color=colors[j],s=32,zorder=2,edgecolors='white',linewidths=.5)
  paired.append(dict(gene=gene,embryo='HET embryo '+str(j+1),sample=s,WT_CPM=val[0],KO_CPM=val[1],paired_delta_log2CPM=ld,KO_higher=bool(ld>0),WT_raw_counts=int(rr[0]['raw_pseudobulk_counts']),KO_raw_counts=int(rr[1]['raw_pseudobulk_counts']),WT_cells=int(rr[0]['cells']),KO_cells=int(rr[1]['cells'])))
 up=int((np.array(changes)>0).sum());summary.append(dict(gene=gene,KO_higher_embryos=up,KO_lower_embryos=int((np.array(changes)<0).sum()),consistent_increase=up==3,existing_DE=de.get(gene,{'status':'not tested after expression filtering'})))
 ax.set_title(gene,style='italic',weight='bold',fontsize=12);ax.text(.5,.96,'KO higher: %d/3 embryos'%up,transform=ax.transAxes,ha='center',va='top',fontsize=8);ax.set_xticks([0,1]);ax.set_xticklabels(['WT cells','KO cells']);ax.set_xlim(-.20,1.20);ax.set_ylim(bottom=0,top=ax.get_ylim()[1]*1.20);ax.tick_params(length=3)
axs[0].set_ylabel('Pseudobulk expression (CPM)')
fig.suptitle('Female HET · same immature inhibitory neuronal subtype',weight='bold',fontsize=11,y=.98)
fig.legend(handles=[Line2D([0],[0],marker='o',color=colors[j],lw=1.2,label='HET embryo '+str(j+1)) for j in range(3)],loc='lower center',bbox_to_anchor=(.5,.015),ncol=3,frameon=False,fontsize=8)
fig.subplots_adjust(left=.085,right=.99,top=.75,bottom=.23,wspace=.35)
for ext in ['png','pdf','svg']:fig.savefig(str(OUT/('paired_HET_IEG_diagnostic.'+ext)),dpi=350,facecolor='white')
with (OUT/'paired_IEG_directions.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(paired[0]),delimiter='\t');w.writeheader();w.writerows(paired)
(OUT/'IEG_DECISION.json').write_text(json.dumps({'neuronal_subtype':'C0006, same provisional immature inhibitory neuronal subtype for all comparisons','gates':'HET WT A-positive vs HET KO A=0/B+C>=1','sample_mapping':dict(zip(['HET embryo 1','HET embryo 2','HET embryo 3'],samples)),'normalization':'six within-HET pseudobulks; same filtered library sizes and TMM factors as existing paired edgeR analysis; low-count excluded genes displayed diagnostically','summary':summary,'supportive_panel_added':False,'reason':'Requested coordinated consistent KO increases are not present; diagnostic only, not evidence for increased IEG regulation','interpretation':'IEG expression alone does not identify a direct PCDH19-sensitive mechanism'},indent=2))
for name in ['plot_paired_HET_IEGs.py','prepare_paired_HET_IEGs.R']:shutil.copy2(str(Path(__file__).parent/name),str(OUT/name))
print(json.dumps(summary,indent=2))
