from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import csv,gzip,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1];SRC=RUNTIME_ROOT/'results/mosaic_deep_dive_20261031';SIG=RUNTIME_ROOT/'results/paired_HET_signaling_20261031';OUT=RUNTIME_ROOT/'results/female_male_signaling_20261031'
def read(p):
 with (gzip.open(str(p),'rt') if str(p).endswith('.gz') else open(str(p))) as f:return list(csv.DictReader(f,delimiter='\t'))
cts=['C0002','C0003','C0006','C0011','C0001'];names=['RG-like progenitors','MGE IPCs','Immature inhibitory neurons','Cycling neuroblasts','Neurogenic progenitors']
groups=[('WT female','WT_F','A_positive'),('HET WT-like','HET_F','A_positive'),('HET KO-like','HET_F','Aneg_BC1')]
meta=read(SRC/'pseudobulk_metadata.tsv');pairs=read(SIG/'expanded_signaling_pairs.tsv');pairs=list({r['interaction']:r for r in pairs if 'Pcdh19' not in (r['ligand_genes']+';'+r['receptor_genes']).split(';')}.values())
genes=set(g for r in pairs for k in ['ligand_genes','receptor_genes'] for g in r[k].split(';'))
expr={};libs=[]
for ct in cts:
 selected=[sorted([r for r in meta if r['stratum']==ct and r['design_group']==design and r['gate']==gate],key=lambda r:r['sample']) for label,design,gate in groups]
 assert all(len(rr)==3 and all(int(r['cells'])>=20 for r in rr) for rr in selected)
 ll=[r for rr in selected for r in rr];libs.extend(ll);denom=np.array([float(r['UMI']) for r in ll]);assert np.all(denom>0)
 accum={}
 with gzip.open(str(SRC/(ct+'_counts.tsv.gz')),'rt') as f:
  reader=csv.DictReader(f,delimiter='\t')
  for r in reader:
   if r['gene_symbol'] not in genes:continue
   vals=np.array([float(r[x['library']]) for x in ll]);accum[r['gene_symbol']]=accum.get(r['gene_symbol'],np.zeros(9))+vals
 for gene,vals in accum.items():
  cpm=vals/denom*1e6
  for gi in range(3):expr[(gi,ct,gene)]=cpm[gi*3:gi*3+3]
edges=[];details=[]
for threshold in [.5,1.,2.]:
 for gi,(label,design,gate) in enumerate(groups):
  for source in cts:
   for target in cts:
    kept=[]
    for pair in pairs:
     components=[expr.get((gi,source,g),np.zeros(3)) for g in pair['ligand_genes'].split(';')]+[expr.get((gi,target,g),np.zeros(3)) for g in pair['receptor_genes'].split(';')]
     present=np.all(np.array(components)>=threshold,axis=0)
     if int(present.sum())<2:continue
     kept.append(pair['interaction'])
     if threshold==1.:details.append(dict(group=label,source=source,receiver=target,interaction=pair['interaction'],pathway=pair['pathway'],ligand=pair['ligand'],receptor=pair['receptor'],support_samples=int(present.sum())))
    edges.append(dict(group=label,source=source,receiver=target,source_name=names[cts.index(source)],receiver_name=names[cts.index(target)],CPM_threshold=threshold,pairs=len(kept)))
for filename,rows in [('celltype_LR_network_counts.tsv',edges),('celltype_LR_network_pairs.tsv',details),('celltype_LR_network_libraries.tsv',libs)]:
 with open(str(OUT/filename),'w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
method=dict(groups=groups,nodes=dict(zip(cts,names)),scope='five separately analyzed consensus populations; provisional neural annotations, not all atlas cells',interaction_database='existing CellChatDB.mouse expanded curated pairs',unique_database_pairs=len(pairs),filter='all ligand/receptor complex subunits >=1 raw-count CPM in the same >=2/3 biological samples; >=20 cells in every selected library; Pcdh19 excluded',edge='number of distinct curated ligand-receptor interaction IDs meeting expression support for each ordered cell-type pair',normalization='gene counts / total library UMI * 1e6; libraries independent by cell type and group',display='same node order, fixed node sizes, common linear edge-width scale across all groups; no edge filtering by DE',interpretation='expression-supported candidate networks, not significance-tested CellChat inference, physical interactions, or established functional signaling',sensitivity_CPM=[.5,1,2],female_WT='A-positive reference, matching exon1 selection for HET WT-like')
(OUT/'CELLTYPE_NETWORK_METHODS.json').write_text(json.dumps(method,indent=2))
for label,_,_ in groups:
 print(label,{str(t):sum(r['pairs'] for r in edges if r['group']==label and r['CPM_threshold']==t) for t in [.5,1,2]},flush=True)
print('NETWORK_COUNTS_COMPLETE',flush=True)
