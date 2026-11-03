"""Restrict the IEG screen to postmitotic-consistent inhibitory neuronal identities."""
from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import csv,json,itertools,shutil
import h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
BASE=RUNTIME_ROOT/'results';OUT=BASE/'postmitotic_inhibitory_IEGs_20261101';OUT.mkdir(parents=True,exist_ok=True)
def st(x):return np.array([v.decode() if isinstance(v,bytes) else str(v) for v in x[:]])
with h5py.File(str(BASE/'mosaic_deep_dive_20261031/per_cell_marker_diagnostics.h5'),'r') as f:
 ids=st(f['cell_id']);sn=st(f['sample']);cl=st(f['cluster']);total=f['total_counts'][:];mg=list(st(f['marker_gene']));mr=f['marker_raw_counts'][:];gn=list(st(f['gate_names']));gates=f['gate_membership'][:]
cc=np.load(str(BASE/'nested_proliferation_20261101/raw_cycle_markers.npz'));assert np.array_equal(cc['cell_id'],ids);cg=list(cc['genes']);cycle=cc['counts'];mi=cycle[:,cg.index('Mki67')];to=cycle[:,cg.index('Top2a')]
gc=np.load(str(BASE/'IEG_neuronal_maturation_20261101/all_cell_IEG_counts.npz'));assert np.array_equal(gc['cell_id'],ids);gene_names=list(gc['genes']);genes=['Egr1','Npas4','Nr4a1'];counts=gc['counts'][:,[gene_names.index(g) for g in genes]];log=np.log1p(counts/total[:,None]*10000.)
core=(mr[:,[mg.index(g) for g in ['Snap25','Syt1','Rbfox3']]]>0).sum(1);gaba=(mr[:,mg.index('Gad1')]+mr[:,mg.index('Gad2')])>0;mge=(mr[:,mg.index('Lhx6')]+mr[:,mg.index('Sox6')])>0;neur=np.isin(cl,['C0004','C0005','C0006','C0020']);base=neur&gaba&(core>=1)
subsets={'Annotated inhibitory neurons':neur,'Inhibitory neurons with neuronal-core support':base,'Cycling-marker-negative inhibitory neurons':base&(mi==0)&(to==0),'Low-cycle inhibitory neurons':base&(mi<=1)&(to<=1),'Cycling-negative with stronger neuronal support':base&(mi==0)&(to==0)&(core>=2),'Cycling-negative MGE-supported inhibitory neurons':base&(mi==0)&(to==0)&mge,'MGE interneuron compartment, cycling-negative':(cl=='C0006')&gaba&(core>=1)&(mi==0)&(to==0),'GABAergic neuron compartment, cycling-negative':(cl=='C0005')&gaba&(core>=1)&(mi==0)&(to==0)}
het=[('HET embryo '+str(j+1),'15662-JZ-'+str(i)) for j,i in enumerate([7,8,9])];wt=gates[:,gn.index('A_positive')].astype(bool);ko=gates[:,gn.index('Aneg_BC1')].astype(bool)
plan={'scope':'annotated inhibitory neuronal clusters C0004/5/6/20 only; exclude C0001 progenitors','postmitotic_consistent':'Gad1 or Gad2 detected; >=1 Snap25/Syt1/Rbfox3 detected; neither Mki67 nor Top2a detected','sensitivities':'<=1 UMI low-cycle, >=2 neuronal core markers, Lhx6 or Sox6 detected, individual compartments','interpretation':'marker-consistent postmitotic state, not direct lineage/division proof; E14.5 dissected MGE does not establish cortical location','literature':'PMID35613587 human NGN2 excitatory cortical neurons, not inhibitory; rat hippocampal neurons','IEG_gating':'none','source_counts_cache':str(BASE/'IEG_neuronal_maturation_20261101/all_cell_IEG_counts.npz')}
(OUT/'PLAN.json').write_text(json.dumps(plan,indent=2));rows=[];summary=[];coverage=[]
for label,mask in subsets.items():
 for embryo,s in het:
  for group,gate in [('WT cells',wt),('KO cells',ko)]:
   m=mask&gate&(sn==s);n=int(m.sum());coverage.append(dict(subset=label,embryo=embryo,sample=s,group=group,cells=n))
   for j,g in enumerate(genes):rows.append(dict(subset=label,embryo=embryo,sample=s,group=group,gene=g,cells=n,raw_counts=int(counts[m,j].sum()),CPM=float(counts[m,j].sum()/total[m].sum()*1e6) if n else float('nan'),mean_log_expression=float(log[m,j].mean()) if n else float('nan')))
 for g in genes:
  rr=[r for r in rows if r['subset']==label and r['gene']==g];aa=[r for r in rr if r['group']=='WT cells'];bb=[r for r in rr if r['group']=='KO cells'];d=np.array([b['mean_log_expression']-a['mean_log_expression'] for a,b in zip(aa,bb)]);dc=np.array([b['CPM']-a['CPM'] for a,b in zip(aa,bb)])
  summary.append(dict(subset=label,gene=g,WT_cells=sum(r['cells'] for r in aa),KO_cells=sum(r['cells'] for r in bb),minimum_WT_cells=min(r['cells'] for r in aa),minimum_KO_cells=min(r['cells'] for r in bb),WT_raw_counts=sum(r['raw_counts'] for r in aa),KO_raw_counts=sum(r['raw_counts'] for r in bb),KO_up_log_embryos=int((d>0).sum()),KO_up_CPM_embryos=int((dc>0).sum()),mean_log_expression_change=float(d.mean())))
def write(name,rr):
 with (OUT/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rr[0]),delimiter='\t');w.writeheader();w.writerows(rr)
write('paired_postmitotic_IEGs.tsv',rows);write('postmitotic_subset_summary.tsv',summary);write('cell_coverage.tsv',coverage)
# Primary independently defined postmitotic-consistent gate, not a result-selected gate.
primary='Cycling-marker-negative inhibitory neurons';plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False});fig,axs=plt.subplots(1,3,figsize=(8.4,3.5));colors=['#557C70','#AA814B','#8F7898']
for ax,g in zip(axs,genes):
 for j,(embryo,s) in enumerate(het):
  rr=[next(r for r in rows if r['subset']==primary and r['sample']==s and r['gene']==g and r['group']==group) for group in ['WT cells','KO cells']];ax.plot([0,1],[r['CPM'] for r in rr],marker='o',color=colors[j],lw=1.1,ms=5)
 ax.set_title(g,fontstyle='italic',weight='bold');ax.set_xticks([0,1]);ax.set_xticklabels(['WT cells','KO cells']);ax.set_ylim(bottom=0)
axs[0].set_ylabel('Pseudobulk expression (CPM)');fig.suptitle('Postmitotic-consistent inhibitory neurons · female HET',weight='bold',fontsize=11);fig.tight_layout(rect=[0,.10,1,.88]);fig.text(.07,.03,'Same neuronal identities; Gad1/2+, neuronal-core+, Mki67−/Top2a−. Colors track three HET embryos.',fontsize=8)
for ext in ['png','pdf','svg']:fig.savefig(str(OUT/('postmitotic_paired_IEGs.'+ext)),dpi=300,facecolor='white')
shutil.copy2(__file__,str(OUT/Path(__file__).name));print(json.dumps([r for r in summary if r['subset']==primary],indent=2),flush=True)
