"""DE-first exploratory signaling report, with expression-supported curated pairs."""
import csv
import itertools
import json
import copy
from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
from collections import defaultdict
import numpy as np
from scipy.stats import hypergeom
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=RUNTIME_ROOT/'results/paired_HET_signaling_20261031'
def read(name):return list(csv.DictReader((OUT/name).open(),delimiter='\t'))
def write(name,rows):
    if not rows:return
    with (OUT/name).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
assert (OUT/'ANALYSIS_COMPLETE.txt').exists()
data=read('paired_DE_all_celltypes_with_signaling_annotations.tsv')
roles=read('signaling_gene_annotations.tsv');pairs=read('expanded_signaling_pairs.tsv')
by=defaultdict(list)
for r in data:by[r['stratum']].append(r)
lookup={(r['stratum'],r['gene_symbol']):r for r in data}
path_genes=defaultdict(set)
for r in roles:path_genes[r['pathway']].add(r['gene'])

# Include pathways with zero hits in the multiplicity correction.
enrich=[]
for st,rows in by.items():
    universe={r['gene_symbol'] for r in rows if r['flag']==''}
    hit={r['gene_symbol'] for r in rows if r['flag']=='' and float(r['FDR'])<.05}
    here=[]
    for path,genes in path_genes.items():
        gs=genes & universe
        if len(gs)<3:continue
        found=hit & gs
        here.append(dict(stratum=st,pathway=path,tested_genes=len(gs),DE_genes=len(found),total_DE=len(hit),universe=len(universe),
                         PValue=float(hypergeom.sf(len(found)-1,len(universe),len(gs),len(hit))),FDR=1.,genes=';'.join(sorted(found))))
    here.sort(key=lambda r:r['PValue']);q=1.
    for i in range(len(here)-1,-1,-1):q=min(q,here[i]['PValue']*len(here)/(i+1));here[i]['FDR']=q
    enrich.extend(here)
write('signaling_pathway_overrepresentation.tsv',enrich)

expr={(r['stratum'],r['gene'],r['gate']):r for r in read('signaling_gene_expression_support.tsv')}
celltypes=['vascular_marked','C0001','C0002','C0003','C0006','C0011']
def support(st,genes,gate):
    rows=[expr.get((st,g,gate)) for g in genes]
    if any(r is None or int(r['samples_CPM_ge1'])<2 for r in rows):return None
    return float(np.exp(np.mean([np.log(max(float(r['mean_CPM']),1e-12)) for r in rows])))
routes=[]
for pair in pairs:
    ligands=pair['ligand_genes'].split(';');receptors=pair['receptor_genes'].split(';')
    for sender,receiver in itertools.product(celltypes,repeat=2):
        if 'vascular_marked' not in [sender,receiver] or sender==receiver:continue
        component_rows=[lookup.get((sender,g)) for g in ligands]+[lookup.get((receiver,g)) for g in receptors]
        changed=[r for r in component_rows if r and r['flag']=='' and (float(r['FDR'])<.1 or float(r['PValue'])<.01)]
        if not changed:continue
        for source_gate,target_gate in [('A_positive','Aneg_BC1'),('Aneg_BC1','A_positive')]:
            l=support(sender,ligands,source_gate);r=support(receiver,receptors,target_gate)
            if l is None or r is None:continue
            routes.append(dict(sender_celltype=sender,receiver_celltype=receiver,sender_gate=source_gate,receiver_gate=target_gate,
                pathway=pair['pathway'],interaction=pair['interaction'],ligand=pair['ligand'],receptor=pair['receptor'],
                ligand_CPM_geomean=l,receptor_CPM_geomean=r,
                changed_component_genes=';'.join(x['gene_symbol'] for x in changed),
                best_component_P=min(float(x['PValue']) for x in changed),best_component_FDR=min(float(x['FDR']) for x in changed)))
routes.sort(key=lambda r:(r['best_component_FDR'],r['best_component_P'],-np.sqrt(r['ligand_CPM_geomean']*r['receptor_CPM_geomean'])))
write('cross_celltype_signaling_expression_candidates.tsv',routes)

allrows=by['all_cells'];top=[r for r in allrows if r['flag']=='' and r['roles']!=''][:20]
write('top20_all_cell_paired_signaling_DEGs.tsv',top)
primary=['Cxcl12','Flt1','Unc5b','Esam','Cdh5','Lamb1','Cspg4','Eng','Col4a1','Vegfc','Dkk2','Sema3g']
genes=[g for g in primary if ('all_cells',g) in lookup]
strata=['all_cells','vascular_marked','C0001','C0002','C0003','C0006','C0011']
names=['All cells','Vascular+','MGE prog.','RG-like','IPC','MGE neuron','Cycling NB']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.linewidth':.7,'pdf.fonttype':42,'svg.fonttype':'none'})
fig=plt.figure(figsize=(10,5.2));gs=fig.add_gridspec(1,2,width_ratios=[1,1.5],left=.08,right=.97,top=.88,bottom=.22,wspace=.65)
ax=fig.add_subplot(gs[0,0])
x=np.array([float(r['logFC']) for r in allrows]);y=-np.log10(np.maximum([float(r['FDR']) for r in allrows],1e-300))
colors=['#A678AC' if r['roles']!='' and r['flag']=='' else '#c4c7ca' for r in allrows]
ax.scatter(x,y,s=5,c=colors,linewidths=0,rasterized=True)
ax.axhline(-np.log10(.05),c='#999999',lw=.7,ls='--')
for g,offset in [('Cxcl12',(-30,12)),('Flt1',(15,15)),('Unc5b',(-30,-13)),('Pcdh19',(25,-7))]:
    r=lookup['all_cells',g];xx=float(r['logFC']);yy=-np.log10(float(r['FDR']))
    ax.annotate(g,(xx,yy),xytext=offset,textcoords='offset points',fontsize=7,arrowprops=dict(arrowstyle='-',lw=.5,color='#666666'))
ax.set_xlabel('log2 fold change\nKO-enriched / WT-like within HET');ax.set_ylabel('−log10 FDR')
ax.set_title('Paired all-cell DEGs',fontsize=10);ax.text(-.25,1.09,'a',transform=ax.transAxes,fontweight='bold',fontsize=12)
ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
ax.set_xlim(-4,4)
ax=fig.add_subplot(gs[0,1]);matrix=np.full((len(genes),len(strata)),np.nan)
for i,g in enumerate(genes):
    for j,st in enumerate(strata):
        r=lookup.get((st,g))
        if r:matrix[i,j]=float(r['logFC'])
cmap=copy.copy(plt.get_cmap('RdBu_r'));cmap.set_bad('#dddddd')
im=ax.imshow(np.ma.masked_invalid(matrix),aspect='auto',cmap=cmap,vmin=-2,vmax=2)
ax.set_yticks(range(len(genes)));ax.set_yticklabels(genes)
ax.set_xticks(range(len(strata)));ax.set_xticklabels(names,rotation=45,ha='right',fontsize=7)
for i,g in enumerate(genes):
    for j,st in enumerate(strata):
        r=lookup.get((st,g))
        if r and float(r['FDR'])<.05:ax.text(j,i,'*',ha='center',va='center',fontsize=11,color='white' if abs(float(r['logFC']))>1 else 'black')
ax.set_title('Where the signaling differences appear',fontsize=10);ax.text(-.2,1.09,'b',transform=ax.transAxes,fontweight='bold',fontsize=12)
fig.colorbar(im,ax=ax,fraction=.05,pad=.04,label='log2 fold change',extend='both')
fig.text(.08,.045,'Paired sample-level edgeR tests, n = 3 HET samples. Purple: curated ligand/receptor-associated genes; * FDR < 0.05.',fontsize=7)
fig.text(.08,.012,'Gray: not tested/expressed. Cell-type names are provisional. Expression support does not demonstrate cell–cell communication.',fontsize=7)
for ext in ['pdf','png','svg']:fig.savefig(OUT/('paired_DE_and_celltype_signaling.'+ext),dpi=400,bbox_inches='tight',facecolor='white')
plt.close(fig)

overview=read('paired_DE_overview.tsv');nr=next(r for r in overview if r['stratum']=='all_cells')
lines=['# Paired WT-like versus KO-enriched HET cells: DE-first signaling review','',
    'Primary comparison: A-negative/B-or-C-positive KO-enriched versus A-positive WT-like cells within the same HET sample. Raw counts were summed per sample and gate and fitted with the existing paired sample-level edgeR model. There are three registered HET samples. Unfiltered WT females remain the figure reference but are not the primary DE contrast here.',
    '', '## All-cell result','',
    '**%s non-gate genes pass FDR <0.05; %s are curated signaling-associated genes.** The strongest changes point to vascular-associated guidance, extracellular matrix and cell-contact expression. Pcdh19 is displayed but treated separately because it defines selection.'%(nr['DE_FDR05'],nr['signaling_DE_FDR05']),
    '', '| Gene | Role | Pathways | log2FC (KO-enriched / WT-like) | FDR |','|---|---|---|---:|---:|']
for r in top[:15]:lines.append('| {gene_symbol} | {roles} | {pathways} | {v:+.2f} | {q:.3g} |'.format(v=float(r['logFC']),q=float(r['FDR']),**r))
lines+=['','## Cell-type localization','', '| Cell stratum | All non-gate DEGs, FDR <0.05 | Signaling DEGs | Nominal signaling candidates, P <0.01 |','|---|---:|---:|---:|']
for r in overview:lines.append('| {stratum} | {DE_FDR05} | {signaling_DE_FDR05} | {signaling_nominal01} |'.format(**r))
lines+=['','The strongest corrected signaling localization is in vascular-marker-positive and mixed-profile cells. Neural-restricted comparisons preserve smaller exploratory leads; they should not be summarized as a global neuronal signaling defect. Blood-associated and Y-linked genes are flagged in the full results and are excluded from signaling prioritization.',
    '', '## Curated pathway enrichment','', '| Stratum | Pathway | DE / tested pathway genes | Enrichment FDR |','|---|---|---:|']
for st in ['all_cells','vascular_marked']:
    rows=sorted([r for r in enrich if r['stratum']==st and r['DE_genes']>0],key=lambda r:r['PValue'])[:8]
    for r in rows:lines.append('| {stratum} | {pathway} | {DE_genes}/{tested_genes} | {FDR:.3g} |'.format(**r))
lines+=['','Enrichment uses an overrepresentation test over expressed/tested genes; all eligible pathways, including those with zero hits, enter BH correction. It annotates gene lists and does not measure pathway activity.',
    '', '## Candidate intercellular axes','',
    'CellChatDB mouse was downloaded from the official jinworks/CellChat repository; its exact bytes and SHA256 are saved. Ligand/receptor complexes are expanded into required components. Non-protein interactions are excluded. A candidate pair requires every ligand and receptor component to be expressed at CPM ≥1 in at least two of three samples of its respective group/cell type, plus a DE or nominally changed component. Optional cofactors, protein abundance, physical contact and spatial adjacency are not established.',
    '', '| Source | Target | Axis | Direction | Changed component(s) | Best component FDR |','|---|---|---|---|---|---:|']
seen=set()
for r in routes:
    if r['pathway'] not in ['CXCL','VEGF','Netrin','SEMA3']:continue
    key=(r['interaction'],r['sender_celltype'],r['receiver_celltype'])
    if key in seen:continue
    seen.add(key)
    lines.append('| {sender_celltype} | {receiver_celltype} | {interaction} | {sender_gate} → {receiver_gate} | {changed_component_genes} | {best_component_FDR:.3g} |'.format(**r))
    if len(seen)>=15:break
lines+=['','These routes are **expression-supported hypotheses**, not CellChat communication probabilities or proof of signaling. Component FDRs are gene-level statistics, not interaction-level significance. Highly overlapping cell strata are sensitivity views, not independent populations.',
    '', '## Figure addition','',
    'Pcdh19 has been added as a fourth expression panel in the UMAP figure. It combines probe-derived gene counts. Both HET groups are selected for probe positivity, whereas WT females are unfiltered; its increase relative to WT females is selection-dependent and must not be presented as an independent disease DEG.',
    '', '## Files','',
    '- paired_DE_and_celltype_signaling.pdf: all-cell DE followed by cell-type localization.',
    '- all_cell_DEGs_followed_by_celltype.tsv: every all-cell significant gene, followed through cell strata.',
    '- top20_signaling_genes_by_celltype.tsv and paired_signaling_DE_all_celltypes.tsv: ranked exploratory signaling leads.',
    '- cross_celltype_signaling_expression_candidates.tsv: curated pair expression availability across HET labels and cell types.',
    '- signaling_pathway_overrepresentation.tsv: complete enrichment results.',
    '- database_source.json: source URL and checksum.']
(OUT/'REPORT.md').write_text('\n'.join(lines)+'\n')
(OUT/'REVIEW_COMPLETE.json').write_text(json.dumps(dict(all_cell_DE_non_gate=int(nr['DE_FDR05']),all_cell_signaling_DE=int(nr['signaling_DE_FDR05']),cell_strata=len(by),curated_interactions=len(pairs),cross_celltype_candidates=len(routes),pathway_BH_includes_zero_hit_sets=True),indent=2))
print('\n'.join(lines[:29]));print('SIGNALING_REVIEW_COMPLETE')
