"""Complete exploratory neuronal/maturation sensitivity screen for prespecified IEGs."""
from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import csv,json,itertools,shutil
import h5py,numpy as np
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];BASE=RUNTIME_ROOT/'results';OUT=BASE/'IEG_neuronal_maturation_20261101';OUT.mkdir(exist_ok=True,parents=True)
DIAG=BASE/'mosaic_deep_dive_20261031/per_cell_marker_diagnostics.h5';RAW=BASE/'primary_processing/06_technical_sample_batch_diagnostics/06_technical_sample_batch_diagnostics_20260903_114110_6087eeb/objects/pcdh19_step06_unintegrated_diagnostics.h5ad'
def strings(d):return np.array([x.decode() if isinstance(x,bytes) else str(x) for x in d[:]])
with h5py.File(str(DIAG),'r') as f:
 ids=strings(f['cell_id']);sample=strings(f['sample']);cluster=strings(f['cluster']);total=f['total_counts'][:];gn=list(strings(f['gate_names']));gates=f['gate_membership'][:];mg=list(strings(f['marker_gene']));markers=f['marker_raw_counts'][:];pn=list(strings(f['program_names']));program=f['program_mean_log1p_CP10K'][:]
key=list(csv.DictReader((ROOT/'config/sample_key.csv').open()));mapping={r['technical_sample_id']:r['design_group'] for r in key};dg=np.array([mapping[s] for s in sample]);samples=[r['technical_sample_id'] for r in key];het_samples=['15662-JZ-'+str(i) for i in [7,8,9]]
positive=gates[:,gn.index('informative')].astype(bool);apos=gates[:,gn.index('A_positive')].astype(bool);khet=gates[:,gn.index('Aneg_BC1')].astype(bool)
group_masks={'WT_F':(dg=='WT_F')&apos,'WT_M_positive':(dg=='WT_M')&positive,'HET_WT':(dg=='HET_F')&apos,'HET_KO':(dg=='HET_F')&khet,'KO_M':dg=='KO_M'}
pooled_wt=(np.isin(dg,['WT_F','WT_M'])&positive)|group_masks['HET_WT'];pooled_ko=group_masks['KO_M']|group_masks['HET_KO']
primary=['Egr1','Npas4','Nr4a1'];secondary=['Fos','Junb','Arc','Egr2','Egr3','Fosb'];genes=primary+secondary
maturation=program[:,pn.index('later_neuronal_maturation')];gad=(markers[:,mg.index('Gad1')]+markers[:,mg.index('Gad2')])>0
neural_detected=(markers[:,[mg.index(g) for g in ['Snap25','Syt1','Rbfox3','Dcx','Tubb3','Elavl3','Stmn2']]]>0).sum(1)
core=(markers[:,[mg.index(g) for g in ['Snap25','Syt1','Rbfox3']]]>0).sum(1)>=2
broad=np.isin(cluster,['C0004','C0005','C0006','C0020']);gaba=gad&(neural_detected>=2)
parents={'Neurogenic progenitors':cluster=='C0001','Cycling neuroblasts':cluster=='C0011','Annotated inhibitory neurons':broad,'MGE interneuron compartment':cluster=='C0006','GABAergic neuron compartment':cluster=='C0005','Ventral GABAergic neuron compartment':cluster=='C0004','Additional GABAergic neuron compartment':cluster=='C0020','Marker-defined GABAergic-lineage cells':gaba}
subsets={};thresholds=[]
for name,parent in parents.items():
 subsets[name]=parent
 for q,label in [(50,'higher maturation'),(75,'top-quarter maturation')]:
  mask=np.zeros(len(sample),bool)
  for sn in samples:
   ix=(sample==sn)&parent
   if ix.any():
    cutoff=float(np.percentile(maturation[ix],q));mask[ix]=maturation[ix]>cutoff;thresholds.append(dict(parent=name,sample=sn,percentile=q,cutoff=cutoff))
  subsets[name+' · '+label]=mask
 subsets[name+' · neuronal-core positive']=parent&core
# More focused exploratory endpoints, still defined without IEG expression.
for name in ['Annotated inhibitory neurons','MGE interneuron compartment']:
 parent=parents[name];mask=np.zeros(len(sample),bool)
 for sn in samples:
  ix=(sample==sn)&parent
  if ix.any():mask[ix]=maturation[ix]>np.percentile(maturation[ix],90)
 subsets[name+' · top-decile maturation']=mask
subsets['Neurogenic progenitors · GABAergic-lineage higher maturation']=(cluster=='C0001')&subsets['Marker-defined GABAergic-lineage cells · higher maturation']
cyclecache=np.load(str(BASE/'nested_proliferation_20261101/raw_cycle_markers.npz'));assert np.array_equal(cyclecache['cell_id'],ids)
cc=cyclecache['counts'];cgenes=list(cyclecache['genes']);noncycle=(cc[:,cgenes.index('Mki67')]==0)&(cc[:,cgenes.index('Top2a')]==0)
for name in ['Annotated inhibitory neurons','MGE interneuron compartment','Marker-defined GABAergic-lineage cells']:subsets[name+' · cycling-marker negative']=parents[name]&noncycle
(OUT/'EXPLORATORY_PLAN.json').write_text(json.dumps({'primary_genes':primary,'secondary_genes':secondary,'subsets':list(subsets),'maturation_marker_genes':['Snap25','Syt1','Rbfox3','Map2','Stmn2','Nefl','Nefm','Syp'],'gating_excludes_primary_IEGs':True,'maturation_thresholds':'per sample, within the full parent population, shared across probe states; strict greater than percentile','primary_comparison':'within each HET embryo A-negative/B+C>=1 versus A-positive; same parent/subset','secondary_comparisons':'pooled KO males+HET KO vs known WT Pcdh19-positive+HET WT; WT male positive vs KO male','all_results_retained':True,'extension_provenance':'Neurogenic and cycling progenitor parents added after the initial 29-subset screen; initial plan/results retained with initial_ prefix; this is post hoc exploration','ranking':'prioritize >=10 cells in each HET library, >=5 raw counts summed per primary gene/state, number of primary genes higher in all 3 embryos, then embryo consistency of fixed 3-gene mean log-normalized score','interpretation':'relative maturation in E14.5 mouse; neither adult maturity nor direct PCDH19 mechanism is inferred'},indent=2))
cache=OUT/'all_cell_IEG_counts.npz'
if cache.exists():
 c=np.load(str(cache));assert np.array_equal(c['cell_id'],ids);counts=c['counts']
else:
 counts=np.zeros((len(ids),len(genes)),dtype=np.int32)
 missing=[]
 for j,g in enumerate(genes):
  if g in mg:counts[:,j]=markers[:,mg.index(g)]
  else:missing.append(j)
 with h5py.File(str(RAW),'r') as f:
  assert np.array_equal(strings(f['obs/cell_id']),ids);v=f['var/gene_symbol']
  if isinstance(v,h5py.Group):cats=strings(v['categories']);symbols=cats[v['codes'][:]]
  else:symbols=strings(v)
  gi=[list(symbols).index(genes[j]) for j in missing];x=f['X'];ptr=x['indptr'][:]
  for start in range(0,len(ids),4000):
   end=min(start+4000,len(ids));lo,hi=int(ptr[start]),int(ptr[end]);mat=sparse.csr_matrix((x['data'][lo:hi],x['indices'][lo:hi],ptr[start:end+1]-lo),shape=(end-start,len(symbols)));counts[start:end,missing]=mat[:,gi].toarray()
   if start%100000==0:print('IEG_RAW',start,flush=True)
 np.savez_compressed(str(cache),cell_id=ids,genes=np.array(genes),counts=counts)
cp=counts/total[:,None]*10000.;log=np.log1p(cp);module=log[:,:3].mean(1)
def write(name,rr):
 with (OUT/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rr[0]),delimiter='\t');w.writeheader();w.writerows(rr)
write('maturation_thresholds.tsv',thresholds)
rows=[];coverage=[];pooled=[]
for subset,mask in subsets.items():
 for group,gate in group_masks.items():
  for sn in sorted(set(sample[gate])):
   m=mask&gate&(sample==sn);n=int(m.sum());coverage.append(dict(subset=subset,group=group,sample=sn,cells=n))
   if not n:continue
   for j,gene in enumerate(genes):rows.append(dict(subset=subset,group=group,sample=sn,gene=gene,cells=n,raw_counts=int(counts[m,j].sum()),mean_log_expression=float(log[m,j].mean()),mean_CP10K=float(cp[m,j].mean()),pseudobulk_CPM=float(counts[m,j].sum()/total[m].sum()*1e6),positive_cells=int((counts[m,j]>0).sum()),detection_percent=float((counts[m,j]>0).mean()*100)))
   rows.append(dict(subset=subset,group=group,sample=sn,gene='Three-gene IEG score',cells=n,raw_counts=int(counts[m,:3].sum()),mean_log_expression=float(module[m].mean()),mean_CP10K=float(cp[m,:3].mean()),pseudobulk_CPM=float(counts[m,:3].sum()/total[m].sum()*1e6/3),positive_cells=int((counts[m,:3].sum(1)>0).sum()),detection_percent=float((counts[m,:3].sum(1)>0).mean()*100)))
 for j,gene in enumerate(genes):
  a=mask&pooled_wt;b=mask&pooled_ko
  if not (a.any() and b.any()):continue
  wr=float(cp[a,j].mean());kr=float(cp[b,j].mean());pooled.append(dict(subset=subset,gene=gene,WT_cells=int(a.sum()),KO_cells=int(b.sum()),WT_mean_CP10K=wr,KO_mean_CP10K=kr,linear_KO_over_WT=kr/wr if wr else None,log_expression_difference=float(log[b,j].mean()-log[a,j].mean())))
write('all_subset_sample_expression.tsv',rows);write('subset_cell_coverage.tsv',coverage);write('all_pooled_subset_effects.tsv',pooled)
lookup={(r['subset'],r['group'],r['sample'],r['gene']):r for r in rows};paired=[];rank=[]
for subset in subsets:
 lead=[]
 for gene in genes+['Three-gene IEG score']:
  pp=[]
  for sn in het_samples:
   aa=lookup.get((subset,'HET_WT',sn,gene));bb=lookup.get((subset,'HET_KO',sn,gene))
   if aa and bb:pp.append((aa,bb))
  if len(pp)!=3:continue
  a=np.array([x[0]['mean_log_expression'] for x in pp]);b=np.array([x[1]['mean_log_expression'] for x in pp]);d=b-a;effect=float(d.mean());null=[float((d*np.array(signs)).mean()) for signs in itertools.product([-1,1],repeat=3)]
  araw=sum(x[0]['raw_counts'] for x in pp);braw=sum(x[1]['raw_counts'] for x in pp);minden=min(x['cells'] for pair in pp for x in pair);ap=np.array([x[0]['pseudobulk_CPM'] for x in pp]);bp=np.array([x[1]['pseudobulk_CPM'] for x in pp]);fold=float(bp.mean()/ap.mean()) if ap.mean() else None
  r=dict(subset=subset,gene=gene,WT_mean_log_expression=float(a.mean()),KO_mean_log_expression=float(b.mean()),difference=effect,KO_higher_embryos=int((d>0).sum()),KO_lower_embryos=int((d<0).sum()),per_embryo_differences=';'.join(str(x) for x in d),WT_raw_counts=araw,KO_raw_counts=braw,minimum_cells_per_library=minden,adequate_coverage=bool(minden>=10 and min(araw,braw)>=5),mean_pseudobulk_KO_over_WT=fold,exact_p=sum(abs(x)>=abs(effect)-1e-12 for x in null)/len(null),KO_higher_CPM_embryos=int((bp>ap).sum()));paired.append(r)
  if gene in primary:lead.append(r)
 if len(lead)!=3:continue
 mr=next(x for x in paired if x['subset']==subset and x['gene']=='Three-gene IEG score')
 rank.append(dict(subset=subset,all_primary_genes_adequate=all(r['adequate_coverage'] for r in lead),primary_genes_up_in_all_3=sum(r['KO_higher_embryos']==3 for r in lead),primary_total_increasing_pairs=sum(r['KO_higher_embryos'] for r in lead),module_up_embryos=mr['KO_higher_embryos'],module_difference=mr['difference'],minimum_cells_per_library=mr['minimum_cells_per_library']))
# Complete-screen correction; subset overlaps are retained, not treated as replication.
order=sorted(range(len(paired)),key=lambda i:paired[i]['exact_p']);q=1.
for k in range(len(order)-1,-1,-1):
 i=order[k];q=min(q,paired[i]['exact_p']*len(order)/(k+1));paired[i]['BH_q']=q
rank.sort(key=lambda r:(not r['all_primary_genes_adequate'],-r['primary_genes_up_in_all_3'],-r['primary_total_increasing_pairs'],-r['module_up_embryos'],-r['module_difference']))
write('all_paired_HET_effects.tsv',paired);write('ranked_subsets.tsv',rank)
# Overview of the top ranked subsets, with all three target genes and pair agreement visible.
show=rank[:12];matrix=np.array([[next(r['difference'] for r in paired if r['subset']==x['subset'] and r['gene']==g) for g in primary] for x in show]);limit=max(float(np.abs(matrix).max()),.001)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42,'svg.fonttype':'none'})
fig=plt.figure(figsize=(8.5,5.8));ax=fig.add_axes([.48,.14,.40,.71]);im=ax.imshow(matrix,cmap='RdBu_r',vmin=-limit,vmax=limit,aspect='auto')
for i,x in enumerate(show):
 for j,g in enumerate(primary):
  r=next(r for r in paired if r['subset']==x['subset'] and r['gene']==g);ax.text(j,i,'%d/3 up'%r['KO_higher_embryos'],ha='center',va='center',fontsize=7,color='white' if abs(matrix[i,j])>.65*limit else '#26292E')
ax.set_yticks(range(len(show)));ax.set_yticklabels([x['subset'].replace(' · ','\n') for x in show],fontsize=7);ax.set_xticks(range(3));ax.set_xticklabels(primary,fontstyle='italic');ax.tick_params(length=0)
for sp in ax.spines.values():sp.set_visible(False)
cax=fig.add_axes([.91,.14,.018,.71]);cb=fig.colorbar(im,cax=cax);cb.set_label('HET KO − HET WT mean log-normalized expression',fontsize=8)
fig.text(.05,.94,'Neuronal maturation exploration · within HET females',weight='bold',fontsize=11);fig.text(.05,.89,'Subsets defined without target IEGs; numbers show embryo direction agreement.',fontsize=8)
for ext in ['png','pdf','svg']:fig.savefig(str(OUT/('maturation_IEG_overview.'+ext)),dpi=300,facecolor='white')
# Best ranked adequate subset, without switching to whichever individual gene happens to increase.
best=rank[0]['subset'];fig,axs=plt.subplots(1,3,figsize=(8.4,3.6));colors=['#557C70','#AA814B','#8F7898']
for ax,g in zip(axs,primary):
 for j,sn in enumerate(het_samples):
  rr=[lookup[(best,gr,sn,g)] for gr in ['HET_WT','HET_KO']];val=[r['pseudobulk_CPM'] for r in rr];ax.plot([0,1],val,color=colors[j],lw=1.2);ax.scatter([0,1],val,color=colors[j],s=28)
 r=next(r for r in paired if r['subset']==best and r['gene']==g);ax.set_title(g+' · %d/3 CPM up'%r['KO_higher_CPM_embryos'],fontstyle='italic',fontsize=10);ax.set_xticks([0,1]);ax.set_xticklabels(['WT cells','KO cells']);ax.set_ylim(bottom=0);ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
axs[0].set_ylabel('Pseudobulk expression (CPM)');fig.suptitle(best+'\nSame subset paired within each HET embryo',weight='bold',fontsize=10);fig.tight_layout(rect=[0,.07,1,.86]);fig.text(.08,.02,'Green / ochre / purple: HET embryos 1 / 2 / 3. Exploratory selected subset; no independent validation.',fontsize=8)
for ext in ['png','pdf','svg']:fig.savefig(str(OUT/('best_subset_paired_IEGs.'+ext)),dpi=300,facecolor='white')
(OUT/'SCREEN_SUMMARY.json').write_text(json.dumps({'subsets':len(subsets),'paired_gene_tests':len(paired),'best_subset':best,'top_ranked':rank[:8],'any_adequate_subset_all_3_genes_up_all_3_embryos':any(r['all_primary_genes_adequate'] and r['primary_genes_up_in_all_3']==3 for r in rank),'raw_cell_ids_aligned':True,'no_depth_matching':True},indent=2));shutil.copy2(__file__,str(OUT/Path(__file__).name));print('SCREEN_COMPLETE',json.dumps(rank[:8],indent=2),flush=True)

# Check the selected broad marker-defined subset within each existing consensus cluster.
bestmask=subsets[best];composition=[];within=[]
for group in ['HET_WT','HET_KO']:
 for sn in het_samples:
  base=bestmask&group_masks[group]&(sample==sn);den=int(base.sum())
  for cl in sorted(set(cluster[base])):
   n=int((base&(cluster==cl)).sum());composition.append(dict(group=group,sample=sn,cluster=cl,cells=n,denominator=den,percent=100*n/den))
for cl in sorted(set(cluster[bestmask])):
 for j,gene in enumerate(primary):
  diffs=[];cpm_diffs=[];nc=[];raws=[]
  for sn in het_samples:
   a=bestmask&group_masks['HET_WT']&(sample==sn)&(cluster==cl);b=bestmask&group_masks['HET_KO']&(sample==sn)&(cluster==cl);nc.extend([int(a.sum()),int(b.sum())])
   if not (a.any() and b.any()):continue
   diffs.append(float(log[b,j].mean()-log[a,j].mean()));cpm_diffs.append(float(counts[b,j].sum()/total[b].sum()-counts[a,j].sum()/total[a].sum()));raws.extend([int(counts[a,j].sum()),int(counts[b,j].sum())])
  if len(diffs)==3:within.append(dict(cluster=cl,gene=gene,minimum_cells=min(nc),KO_higher_log_embryos=sum(x>0 for x in diffs),KO_higher_CPM_embryos=sum(x>0 for x in cpm_diffs),mean_log_expression_difference=float(np.mean(diffs)),raw_counts_total=sum(raws)))
write('best_subset_cluster_composition.tsv',composition);write('best_subset_within_cluster_checks.tsv',within)
