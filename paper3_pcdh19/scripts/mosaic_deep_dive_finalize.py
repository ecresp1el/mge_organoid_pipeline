"""Write the biological interpretation and transparent three-group candidate plots."""
import csv
import json
import math
from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
OUT=RUNTIME_ROOT/'results/mosaic_deep_dive_20261031'
def read(p):return list(csv.DictReader(p.open(),delimiter='\t'))
values=read(OUT/'selected_biological_candidates_sample_values.tsv')
stats=read(OUT/'selected_biological_candidates_statistics.tsv')
def stat(g,st,co):return next(r for r in stats if r['gene']==g and r['stratum']==st and r['contrast']==co)
paired='within_HET_Aneg_BC1_vs_Apositive';interaction='gate_interaction_BC1_vs_Apositive'
selected=[('Gabrr2','neural_vascular_low','Neural; vascular-low'),('Kcnc2','C0006','Immature MGE neurons'),
          ('Pcdhga12','C0003','MGE IPCs'),('Cacng8','C0003','MGE IPCs; low-count lead'),
          ('Gria1','neural_depth_trimmed','Neural; depth-trimmed'),('Gadd45g','C0001_neural_vascular_low','MGE progenitors; vascular-low')]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.linewidth':.7,'pdf.fonttype':42,'svg.fonttype':'none'})
for name,candidates,shape,size in [('three_group_mosaic_neuronal_leads',selected[:3],(1,3),(9,3.7)),
        ('exploratory_neuronal_mosaic_candidates',selected,(2,3),(9,6.8)),
        ('exploratory_vascular_mosaic_candidates',[(g,'vascular_marked','Vascular-marker positive') for g in ['Unc5b','Cxcl12','Htra1']],(1,3),(9,3.7))]:
    fig,axes=plt.subplots(*shape,figsize=size,squeeze=False)
    for i,(ax,(g,st,label)) in enumerate(zip(axes.flat,candidates)):
        r=[r for r in values if r['gene']==g and r['stratum']==st]
        assert len(r)==9
        for s in ['15662-JZ-7','15662-JZ-8','15662-JZ-9']:
            y=[float(next(t for t in r if t['sample']==s and t['group']==group)['TMM_log2CPM']) for group in ['A_positive','Aneg_BC1']]
            ax.plot([1,2],y,color='#c4c4c4',lw=.8,zorder=1)
        for x,(group,color) in enumerate(zip(['WT_F_all','A_positive','Aneg_BC1'],['#549E89','#337EAC','#A678AC'])):
            y=[float(t['TMM_log2CPM']) for t in r if t['group']==group]
            ax.scatter(x+np.array([-.08,0,.08]),y,s=27,color=color,edgecolors='white',linewidths=.4,zorder=3)
            ax.plot([x-.19,x+.19],[np.mean(y),np.mean(y)],color=color,lw=1.7,zorder=4)
        p=stat(g,st,paired)
        ax.set_title(g+'\n'+label,fontsize=9,pad=28)
        ax.text(.5,1.04,'HET paired P %.2g; q %.2g'%(float(p['PValue']),float(p['FDR'])),transform=ax.transAxes,ha='center',fontsize=7)
        ax.text(-.2,1.3,'abcdef'[i],transform=ax.transAxes,fontweight='bold',fontsize=12)
        ax.set_xticks([0,1,2]);ax.set_xticklabels(['WT F','HET F\nWT-like','HET F\nKO-enriched'],fontsize=7)
        ax.set_xlim(-.4,2.4);ax.set_ylabel('log2 normalized CPM')
        for side in ['top','right']:ax.spines[side].set_visible(False)
    fig.subplots_adjust(left=.08,right=.99,bottom=.16 if shape[0]>1 else .26,top=.84 if shape[0]>1 else .68,wspace=.5,hspace=1.1)
    fig.text(.08,.035,'Each dot is one sample; gray lines pair HET subsets. WT-like = A+; KO-enriched = A−, B/C+. Exploratory selection.',fontsize=7)
    if shape[0]>1:fig.text(.08,.008,'Cacng8 has 0, 0 and 1 raw counts in the WT-like IPC libraries; its large fold change needs particular caution.',fontsize=7)
    for ext in ['pdf','svg','png']:fig.savefig(OUT/(name+'.'+ext),dpi=400,bbox_inches='tight',facecolor='white')
    plt.close(fig)

lines=['# Biological findings: exploratory mosaic-definition screen','',
       '**The most informative definition is broad, paired, cell-type-restricted comparison within HET embryos, with unfiltered WT females retained as the three-group reference.** A-positive HET cells are the WT-like group; A-negative cells with any B/C evidence are the broad KO-enriched group. This retains 3,949 and 16,172 cells, respectively, rather than only 139 strict KO-like cells.',
       '', 'The screen tested 14 probe definitions and 11 cell strata. Of 160 expression comparisons, 152 completed and eight lacked enough cells in at least one library. All results, including null definitions, are retained. These are exploratory leads selected after screening, not independently replicated discoveries.',
       '', '## Main neuronal leads','',
       '| Gene / context | HET WT-like vs WT-group log2FC | HET KO-enriched vs WT-group log2FC | Paired KO-enriched vs WT-like log2FC | Paired P | Paired FDR | Gate-interaction P | Gate-interaction FDR |',
       '|---|---:|---:|---:|---:|---:|---:|---:|']
for g,st,label in selected:
    a=stat(g,st,'HET_A_positive_vs_WTall');b=stat(g,st,'HET_Aneg_BC1_vs_WTall');p=stat(g,st,paired);ix=stat(g,st,interaction)
    lines.append('| %s / %s | %+.2f | %+.2f | %+.2f | %.3g | %.3g | %.3g | %.3g |'%(g,label,float(a['logFC']),float(b['logFC']),float(p['logFC']),float(p['PValue']),float(p['FDR']),float(ix['PValue']),float(ix['FDR'])))
lines += ['',
    '- **Gabrr2:** lower in KO-enriched neural cells than both the WT-like HET group and the unfiltered WT-female group. The same-gate WT/HET contrast is also lower (nominal P=0.0046). This is a candidate inhibitory-receptor expression difference, not a demonstrated physiological deficit.',
    '- **Kcnc2 and Pcdhga12:** WT-like cells in HET tissue sit below the WT-group baseline; KO-enriched cells sit closer to the WT baseline. This pattern is relevant to the hypothesis that the WT-like population in mosaic tissue may itself be altered. The tests are exploratory and do not establish cellular interference or allele-specific causation.',
    '- **Cacng8:** the same WT-like-specific reduction appears in IPCs, but the WT-like libraries have only 0, 0 and 1 raw counts, from 23, 47 and 60 cells. Treat its large fold change as a low-count validation target, not a headline discovery.',
    '- **Gria1:** decreases within paired HET neural subsets after excluding vascular-marker-positive cells and the top 10% RNA-content cells. The paired FDR is approximately 0.002, but the WT/HET gate interaction is not significant; this does not establish a uniquely HET-specific effect.',
    '- **Gadd45g:** a small, directionally consistent progenitor signal. Its interaction is nominally strong but does not pass FDR 0.05.',
    '', 'Gene-function references: [Gabrr2, GABA receptor rho-2](https://www.uniprot.org/uniprotkb/P56476/entry); [Kcnc2, potassium channel](https://www.ncbi.nlm.nih.gov/gene/268345); [Pcdhga12, cell adhesion](https://www.ncbi.nlm.nih.gov/gene/93724); [Cacng8, AMPA-receptor regulation](https://www.ncbi.nlm.nih.gov/gene/81905). These annotations motivate hypotheses; no pathway activity or electrophysiology was measured.',
    '', '## Stronger vascular-associated mosaic-pattern lead','',
    'Within vascular-marker-positive cells, broad KO-enriched versus WT-like HET subsets show lower Unc5b, Cxcl12 and Htra1. WT-like and KO-enriched HET means fall on opposite sides of the WT-group mean for these genes. The paired design is more informative than treating the two HET subsets as unrelated samples.',
    '', '| Gene | Paired HET log2FC | Paired FDR | WT/HET gate-interaction log2FC | Interaction FDR |', '|---|---:|---:|---:|---:|']
for g in ['Unc5b','Cxcl12','Htra1']:
    p=stat(g,'vascular_marked',paired);ix=stat(g,'vascular_marked',interaction)
    lines.append('| %s | %+.2f | %.3g | %+.2f | %.3g |'%(g,float(p['logFC']),float(p['FDR']),float(ix['logFC']),float(ix['FDR'])))
lines += ['',
    'Unc5b and Cxcl12 retain FDR <0.05 in the genotype-by-gate interaction; Htra1 is an exploratory FDR <0.10 lead. This supports prioritizing a vascular-associated guidance/signaling hypothesis. The vascular-marker-positive stratum is not a purified endothelial population, so different vascular subtypes or mixed profiles remain possible explanations.',
    '', '## Findings that should not be turned into the biological story','',
    '- Whole HET versus whole WT females has no FDR <0.05 expression differences after excluding the gate gene in the main whole-population analysis. Some subset signals can therefore be missed by the whole-group average, but the lack of whole-group significance is not proof of cancellation.',
    '- Strong probe-count gates enrich high-RNA, vascular-marker-positive cells in known WT controls too. Whole WT and HET vascular-marker-positive fractions are both approximately 1.9%. The original strict-subset vascular enrichment alone is not evidence of increased vascular-cell abundance in HET tissue.',
    '- Several strong hits are Y-linked genes, hemoglobins or platelet-associated genes. They remain in the complete results but are not interpreted as neuronal mechanisms. The prior sample-sex/probe audit remains unresolved.',
    '- Recurrence across related gates uses the same embryos and overlapping cells. It is sensitivity evidence, not independent biological replication. Broad KO-enriched groups are not definitive DNA-genotype calls.',
    '', '## Recommended reading order','',
    '1. [Neuronal candidate figure](exploratory_neuronal_mosaic_candidates.pdf): all three groups with individual samples and paired HET lines.',
    '2. [Vascular candidate figure](exploratory_vascular_mosaic_candidates.pdf): the stronger corrected guidance/signaling lead.',
    '3. [Gate diagnostics](gate_selection_diagnostics.pdf): selection, RNA content and vascular-marker positivity.',
    '4. [Full definition map](definition_signal_map.pdf) and [complete screen](REPORT.md).',
    '5. selected_biological_candidates_statistics.tsv and selected_biological_candidates_sample_values.tsv: exact statistics, library cell counts, raw gene UMIs and plotted values.']
(OUT/'BIOLOGICAL_FINDINGS.md').write_text('\n'.join(lines)+'\n')
# Fix an imprecise field name in the automatically generated sensitivity table:
# the underlying computation and report heading were means, not medians.
p=OUT/'recurring_exploratory_gene_candidates.tsv'
text=p.read_text();p.write_text(text.replace('median_logFC','mean_logFC',1))
assert all(len([r for r in values if r['gene']==g and r['stratum']==st])==9 for g,st,_ in selected)
validation=json.loads((OUT/'REVIEW_COMPLETE.json').read_text());validation.update(candidate_panels=9,candidate_values_verified=True,interpretation_complete=True)
(OUT/'REVIEW_COMPLETE.json').write_text(json.dumps(validation,indent=2))
print('BIOLOGICAL_REVIEW_COMPLETE',json.dumps(validation))
