"""Display supported marker effects relative to the same WT-female reference."""
import csv
from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
OUT=RUNTIME_ROOT/'results/probe_subset_DE_20261031'
cases=['HET_WT_like','HET_KO_like']
tables={c:{r['gene_symbol']:r for r in csv.DictReader((OUT/(c+'_vs_WT_F_all_genes.tsv')).open(),delimiter='\t')} for c in cases}
genes=['Epas1','Igfbp7','Cavin1','Kdr','Cdh5','Erg','Sox17','Rgs5','Abcc9','Cspg4',
       'Gad1','Gad2','Lhx6','Dcx','Ascl1','Mki67','Top2a']
values=np.array([[0]+[float(tables[c][g]['logFC']) for c in cases] for g in genes])
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'svg.fonttype':'none'})
fig,ax=plt.subplots(figsize=(5.5,7.2))
im=ax.imshow(values,cmap='RdBu_r',vmin=-3,vmax=3,aspect='auto')
ax.set_yticks(range(len(genes)));ax.set_yticklabels(genes)
ax.set_xticks([0,1,2]);ax.set_xticklabels(['WT female\nreference','WT-like\nin HET female','KO-like\nin HET female'])
ax.xaxis.tick_top();ax.tick_params(length=0,pad=8)
ax.axhline(9.5,color='white',lw=3)
for i,g in enumerate(genes):
    for j,c in enumerate(cases,1):
        if float(tables[c][g]['FDR'])<=.05:ax.text(j,i,'*',ha='center',va='center',fontsize=15,color='white')
fig.colorbar(im,ax=ax,fraction=.06,pad=.07,label='log2 fold change versus WT females')
fig.subplots_adjust(left=.19,right=.82,top=.84,bottom=.19)
fig.text(.19,.95,'Three-group biological comparison',fontsize=12,fontweight='bold')
fig.text(.19,.025,'* FDR < 0.05; n = 3 samples/group.\nWT reference is zero by definition; Pcdh19 excluded.\nUpper genes: vascular/perivascular markers.\nLower genes: neural identity and proliferation markers.',fontsize=8)
for ext in ['pdf','svg','png']:fig.savefig(OUT/('three_group_biological_summary.'+ext),dpi=400,bbox_inches='tight',facecolor='white')
with (OUT/'three_group_biological_summary_values.tsv').open('w') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['gene','comparison','log2FC','FDR'])
    for g in genes:
        for c in cases:w.writerow([g,c,tables[c][g]['logFC'],tables[c][g]['FDR']])
print('Saved three-group biological marker summary.')
