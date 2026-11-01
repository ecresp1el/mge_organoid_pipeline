from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import runpy,csv,json,collections
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=RUNTIME_ROOT/'results/female_male_signaling_20261031'
d=runpy.run_path(str(ROOT/'scripts/prepare_HET_crossstate_networks.py'));cts=d['cts'];pairs=d['pairs'];expr=d['expr'];libs=d['libs'];routes=d['routes'];names=d['names']
genes=sorted(d['genes']);gi={g:i for i,g in enumerate(genes)};den=np.zeros((2,5,3));raw=np.zeros((2,5,len(genes),3),dtype=np.int64)
for state in range(2):
 for ci,ct in enumerate(cts):
  gate=['A_positive','Aneg_BC1'][state];rr=sorted([r for r in libs if r['stratum']==ct and r['gate']==gate],key=lambda r:r['sample']);assert len(rr)==3
  den[state,ci]=[float(r['UMI']) for r in rr]
  for gene in genes:raw[state,ci,gi[gene]]=np.rint(expr.get((state,ct,gene),np.zeros(3))*den[state,ci]/1e6).astype(int)
maxsub=max(len(p['ligand_genes'].split(';')+p['receptor_genes'].split(';')) for p in pairs);lig=[];rec=[]
for p in pairs:
 l=[gi[g] for g in p['ligand_genes'].split(';')];r=[gi[g] for g in p['receptor_genes'].split(';')];lig.append(l+[l[0]]*(maxsub-len(l)));rec.append(r+[r[0]]*(maxsub-len(r)))
lig=np.array(lig);rec=np.array(rec)
def support_counts(cpm):
 ll=np.all(cpm[:,:,lig,:]>=1,axis=3);rr=np.all(cpm[:,:,rec,:]>=1,axis=3)
 # state, type, pair, sample
 result=np.zeros((4,5,5));per_sample=np.zeros((4,5,5,3))
 for ri,(a,b) in enumerate(routes):
  for i in range(5):
   for j in range(5):
    present=ll[a,i]&rr[b,j];result[ri,i,j]=np.sum(present.sum(1)>=2);per_sample[ri,i,j]=present.sum(0)
 return result,per_sample
cpm=raw/den[:,:,None,:]*1e6;counts,per_sample=support_counts(cpm)
original={(r['group'],r['source'],r['receiver']):r['pairs'] for r in d['edges'] if r['CPM_threshold']==1}
labels=[d['groups'][a][0]+' → '+d['groups'][b][0] for a,b in routes]
assert all(counts[k,i,j]==original[(labels[k],a,b)] for k in range(4) for i,a in enumerate(cts) for j,b in enumerate(cts))
rng=np.random.RandomState(20261101);minimum=den.min(0);balanced=[];balanced_samples=[]
for iteration in range(50):
 thin=rng.binomial(raw,minimum[None,:,None,:]/den[:,:,None,:]);bc,bps=support_counts(thin/minimum[None,:,None,:]*1e6);balanced.append(bc);balanced_samples.append(bps)
balanced=np.array(balanced);balanced_samples=np.array(balanced_samples);diff=counts-counts[0];bd=balanced-balanced[:,0:1]
rows=[]
for ri in range(1,4):
 for i,a in enumerate(cts):
  for j,b in enumerate(cts):
   pd=per_sample[ri,i,j]-per_sample[0,i,j];bb=bd[:,ri,i,j];rows.append(dict(route=labels[ri],source=a,receiver=b,source_name=names[i],receiver_name=names[j],WT_reference=int(counts[0,i,j]),route_count=int(counts[ri,i,j]),difference=int(diff[ri,i,j]),sample_differences=';'.join(str(int(v)) for v in pd),consistent_samples=bool(np.all(pd>0) or np.all(pd<0)),depth_balanced_mean_difference=float(bb.mean()),depth_balanced_p025=float(np.percentile(bb,2.5)),depth_balanced_p975=float(np.percentile(bb,97.5))))
with open(str(OUT/'HET_network_difference_audit.tsv'),'w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
np.savez_compressed(str(OUT/'HET_network_difference_matrices.npz'),counts=counts,differences=diff,depth_balanced_differences=bd,per_sample_counts=per_sample,cts=np.array(cts),names=np.array(names),routes=np.array(labels))
# Continuous individual-pair effects, independent of network-count threshold crossing.
continuous=[]
for i,a in enumerate(cts):
 for j,b in enumerate(cts):
  for pi,pair in enumerate(pairs):
   l=cpm[:,i,[gi[g] for g in pair['ligand_genes'].split(';')],:];r=cpm[:,j,[gi[g] for g in pair['receptor_genes'].split(';')],:]
   if np.any(l<1) or np.any(r<1):continue
   l=np.exp(np.log(l).mean(1));r=np.exp(np.log(r).mean(1));potential=np.array([np.sqrt(l[x]*r[y]) for x,y in routes]);changes=np.log2(potential[3]/potential[0]);same=bool(np.all(changes>0) or np.all(changes<0))
   continuous.append(dict(source=a,receiver=b,source_name=names[i],receiver_name=names[j],interaction=pair['interaction'],ligand=pair['ligand'],receptor=pair['receptor'],module=pair['pathway'],median_log2_change=float(np.median(changes)),mean_log2_change=float(changes.mean()),consistent_direction=same,sample_log2_changes=';'.join(str(v) for v in changes),minimum_component_CPM=float(min(cpm[:,i,lig[pi],:].min(),cpm[:,j,rec[pi],:].min()))))
continuous.sort(key=lambda r:abs(r['median_log2_change']),reverse=True)
with open(str(OUT/'HET_network_continuous_pair_audit.tsv'),'w') as f:
 w=csv.DictWriter(f,fieldnames=list(continuous[0]),delimiter='\t');w.writeheader();w.writerows(continuous)
report={'audit':'exact independent vectorized reconstruction matches all 100 exported raw edge counts','overlap':'descriptive interaction IDs, not significant communication','depth_check':'50 binomial thinnings to matched UMI depth per cell type and HET sample across probe states; intervals are simulation ranges, not biological confidence intervals','largest_raw_difference':max(rows,key=lambda r:abs(r['difference'])),'within_state_total_counts':counts[[0,3]].sum((1,2)).tolist(),'within_state_depth_balanced_mean_total_difference':float(bd[:,3].sum((1,2)).mean()),'within_state_depth_balanced_total_difference_range':np.percentile(bd[:,3].sum((1,2)),[2.5,97.5]).tolist()}
(OUT/'HET_NETWORK_AUDIT.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2),flush=True)
print('TOP_CONTINUOUS_CONSISTENT',[(r['source_name'],r['receiver_name'],r['interaction'],round(r['median_log2_change'],3)) for r in continuous if r['consistent_direction']][:15],flush=True)
