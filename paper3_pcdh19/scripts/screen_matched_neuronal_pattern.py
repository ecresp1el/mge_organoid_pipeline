from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import csv,gzip,json,numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=RUNTIME_ROOT/'results/neuronal_matched_selection_20261031';BASE=RUNTIME_ROOT/'results/mosaic_deep_dive_20261031'
def read(p):
 with (gzip.open(str(p),'rt') if str(p).endswith('.gz') else open(str(p))) as f:return list(csv.DictReader(f,delimiter='\t'))
md=read(OUT/'nine_sample_metadata.tsv');order=[('WT_F','A_positive'),('HET_F','A_positive'),('HET_F','Aneg_BC1')];md=sorted(md,key=lambda x:(order.index((x['design_group'],x['gate'])),x['sample']));libs=[r['library'] for r in md]
paired={r['gene_symbol']:r for r in read(BASE/'DE_C0006/within_HET_Aneg_BC1_vs_Apositive.tsv.gz')};wt={r['gene_symbol']:r for r in read(OUT/'KOlike_vs_matched_WT_DE.tsv')};baseline={r['gene_symbol']:r for r in read(BASE/'DE_C0006/matched_A_positive.tsv.gz')}
programs={'adhesion_synapse','neurite_guidance','neural_core','MGE_interneuron','early_neuronal_maturation','later_neuronal_maturation','immature_inhibitory','neuroblast'}
neuronal={r['gene'] for r in read(BASE/'marker_program_gene_coverage.tsv') if r['program'] in programs}
neuronal.update('Ank3 Kcnc2 Gnao1 Ptprd Epha5 Snap47 Pcdha5 Pcdha10 Pcdhga5 Lhx8 Sptbn2 Mctp1 Dlx1 Robo2 Erbb4 Gria2 Gria3 Grin1 Grin2b Gabra2 Gabrb2 Gabrb3 Gphn Nrxn1 Nrxn2 Nlgn1 Nlgn2 Slc12a5 Cntnap2 Gabrr2 Cacng8'.split())
rows=[]
for r in read(OUT/'nine_sample_log2CPM.tsv'):
 g=r['gene_symbol']
 if g not in paired or g not in wt or g=='Pcdh19':continue
 x=np.array([float(r[l]) for l in libs]);a,b,c=[float(x[k:k+3].mean()) for k in [0,3,6]]
 p=paired[g];w=wt[g];d=x[6:9]-x[3:6];gap=abs(b-a);drop=min(a-c,b-c)
 pattern=gap<=.5 and drop>=.2 and float(p['logFC'])<0 and float(w['logFC'])<0
 passed=pattern and float(p['PValue'])<.05 and sum(d<0)>=2
 rows.append(dict(gene=g,neuronal=g in neuronal,WT_mean=a,HET_WTlike_mean=b,HET_KOlike_mean=c,WTlike_minus_WT=b-a,KOlike_minus_WTlike=c-b,KOlike_minus_WT=c-a,paired_logFC=p['logFC'],paired_P=p['PValue'],paired_FDR=p['FDR'],KO_vs_WT_logFC=w['logFC'],KO_vs_WT_P=w['PValue'],KO_vs_WT_FDR=w['FDR'],WTlike_vs_WT_P=baseline.get(g,{}).get('PValue',''),paired_samples_lower=int(sum(d<0)),pattern=pattern,exploratory_candidate=passed,ranking_P=max(float(p['PValue']),float(w['PValue']))))
rows.sort(key=lambda r:(not r['exploratory_candidate'],not r['pattern'],r['ranking_P']))
with open(str(OUT/'D_pattern_screen_all_genes.tsv'),'w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
for label,selected in [('neuronal',[r for r in rows if r['neuronal'] and r['exploratory_candidate']]),('all',[r for r in rows if r['exploratory_candidate']])]:
 print(label,len(selected),flush=True)
 for r in selected[:25]:print(r['gene'],round(r['WTlike_minus_WT'],2),round(r['KOlike_minus_WTlike'],2),round(r['KOlike_minus_WT'],2),'p',r['paired_P'],r['KO_vs_WT_P'],'q',r['paired_FDR'],r['KO_vs_WT_FDR'],flush=True)
(OUT/'pattern_screen_definition.json').write_text(json.dumps({'pattern':'observed WT and HET WT-like means within 0.5 log2CPM; HET KO-like >=0.2 log2CPM lower than both; negative modeled contrasts','exploratory_nominal_filter':'paired P<0.05, decrease in >=2/3 HET pairs','ranking':'larger nominal P of KO versus WT-like and KO versus matched WT','WT_similarity_is_not_equivalence':True,'genes_screened':len(rows),'neural_annotation':'existing frozen marker programs plus previously listed neuronal candidates','Gria1_excluded_from_figure':True},indent=2))
