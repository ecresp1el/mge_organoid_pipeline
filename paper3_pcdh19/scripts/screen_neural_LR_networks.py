from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import csv,json,math,collections
ROOT=Path(__file__).resolve().parents[1];SRC=RUNTIME_ROOT/'results/paired_HET_signaling_20261031';OUT=RUNTIME_ROOT/'results/female_male_signaling_20261031'
def read(p):
 with open(str(p)) as f:return list(csv.DictReader(f,delimiter='\t'))
cts=['C0001','C0002','C0003','C0006','C0011'];expr={(r['stratum'],r['gene'],r['gate']):r for r in read(SRC/'signaling_gene_expression_support.tsv')};de={(r['stratum'],r['gene_symbol']):r for r in read(SRC/'paired_DE_all_celltypes_with_signaling_annotations.tsv')}
def support(st,genes,gate):
 r=[expr.get((st,g,gate)) for g in genes]
 if any(x is None or int(x['samples_CPM_ge1'])<2 for x in r):return None
 return math.exp(sum(math.log(max(float(x['mean_CPM']),1e-8)) for x in r)/len(r))
edges=[];pathgenes=collections.defaultdict(dict)
for pair in read(SRC/'expanded_signaling_pairs.tsv'):
 lg=pair['ligand_genes'].split(';');rg=pair['receptor_genes'].split(';')
 for sender in cts:
  for receiver in cts:
   vals=[support(sender,lg,g) for g in ['A_positive','Aneg_BC1']]+[support(receiver,rg,g) for g in ['A_positive','Aneg_BC1']]
   if any(v is None for v in vals):continue
   changed=[]
   for st,gs in [(sender,lg),(receiver,rg)]:
    for gene in gs:
     d=de.get((st,gene))
     if d and d['flag']=='' and float(d['PValue'])<.05:
      changed.append((gene,st,float(d['PValue']),float(d['FDR'])))
   if not changed:continue
   wt=math.sqrt(vals[0]*vals[2]);ko=math.sqrt(vals[1]*vals[3]);change=math.log(ko/wt,2)
   edges.append(dict(pathway=pair['pathway'],interaction=pair['interaction'],sender=sender,receiver=receiver,ligand=pair['ligand'],receptor=pair['receptor'],WTlike_score=wt,KOlike_score=ko,log2_score_change=change,changed_components=';'.join(sorted(set(x[0] for x in changed))),best_P=min(x[2] for x in changed),best_FDR=min(x[3] for x in changed)))
   for gene,st,p,q in changed:pathgenes[pair['pathway']][(gene,st)]=(p,q)
with open(str(OUT/'neural_LR_edges.tsv'),'w') as f:
 w=csv.DictWriter(f,fieldnames=list(edges[0]),delimiter='\t');w.writeheader();w.writerows(edges)
summary=[]
for path,gs in pathgenes.items():
 ee=[e for e in edges if e['pathway']==path]
 summary.append(dict(pathway=path,component_strata=len(gs),genes=';'.join(sorted(set(g for g,st in gs))),best_P=min(p for p,q in gs.values()),best_FDR=min(q for p,q in gs.values()),edges=len(ee),score=sum(-math.log10(p) for p,q in gs.values()),median_absolute_change=sorted(abs(e['log2_score_change']) for e in ee)[len(ee)//2]))
summary.sort(key=lambda r:r['score'],reverse=True)
with open(str(OUT/'neural_LR_pathway_screen.tsv'),'w') as f:
 w=csv.DictWriter(f,fieldnames=list(summary[0]),delimiter='\t');w.writeheader();w.writerows(summary)
for r in summary[:20]:print(r,flush=True)
(OUT/'LR_method.json').write_text(json.dumps({'database':'curated CellChatDB.mouse interactions; existing downloaded database','scope':'five neural consensus clusters, both HET probe-defined groups','expression_support':'every ligand/receptor subunit >=1 CPM in >=2/3 HET samples in each group','score':'sqrt(geometric mean ligand CPM * geometric mean receptor CPM)','contrast':'log2 KO-like score / WT-like score','component_filter':'at least one non-gate, non-blood, non-Y gene nominal P<0.05 in paired HET DE','pathway_rank':'sum -log10 paired P over unique gene-stratum pairs; exploratory screen','not_CellChat_probability':True,'not_measured_communication':True},indent=2))
