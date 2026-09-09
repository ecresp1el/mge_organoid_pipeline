#!/usr/bin/env python3
"""Render the bounded cross-study identity extension from completed tables."""
import os
os.environ['MPLBACKEND']='Agg'
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from div30_identity_external import OUT,PROJECT,OLD,STUDIES,PROGRAMS,MAP,transform
os.environ['MPLCONFIGDIR']=str(OUT/'cache/matplotlib')
LABELS={'mgeo_div30':'MGEO DIV30','walsh':'Walsh DIV75','bershteyn_2023':'Bershteyn 2023','bershteyn_2025':'Bershteyn 2025','siebert_2026':'Siebert registered object','samarasinghe_2021':'Samarasinghe Ctrl'}
COLORS={'mgeo_div30':'#232323','walsh':'#24788c','bershteyn_2023':'#ac6534','bershteyn_2025':'#8b7930','siebert_2026':'#73518a','samarasinghe_2021':'#398052'}
def table(name):return pd.read_csv(OUT/'tables'/name,sep='\t')
def md(x):
    def f(v):return f'{v:.4g}' if isinstance(v,(float,np.floating)) else str(v).replace('|','\\|')
    return '\n'.join(['| '+' | '.join(x.columns)+' |','| '+' | '.join(['---']*len(x.columns))+' |']+['| '+' | '.join(f(v) for v in row)+' |' for row in x.itertuples(index=False,name=None)])
def save(fig,name):
    for ext in ['png','pdf','svg']:
        p=OUT/'figures'/ext;p.mkdir(exist_ok=True,parents=True)
        fig.savefig(p/(name+'.'+ext),dpi=600 if ext=='png' else 300,bbox_inches='tight',facecolor='white')
    plt.close(fig)
def main():
    freq=table('external_reference_resemblance_frequency.tsv');cal=table('mgeo_held_sample_resemblance_calibration.tsv')
    rank=table('matched_program_rankings.tsv');match=table('matched_state_overlap.tsv');bal=table('matched_state_balance.tsv');genes=table('matched_genes_by_sample.tsv')
    noslc=table('matched_genes_by_sample_omit_SLC6A1.tsv') if (OUT/'tables/matched_genes_by_sample_omit_SLC6A1.tsv').exists() else pd.DataFrame()
    inv=table('local_input_inventory.tsv');coverage=table('gene_coverage.tsv')
    with __import__('h5py').File(OLD/'cache/mapping_model.h5') as f:fg=list(f['features'].asstr()[:])
    cov=coverage[coverage.gene.isin(fg)].groupby('study').present.sum().to_dict()
    details=[
      dict(study='mgeo_div30',category='MGEO-SOSR DIV30 original frozen six-sample population',capture='Whole-organoid fixation with10x Fixed RNA Sample Preparation Kit1000414; paper calls subsequent libraries v.3 but exact probe-versus3-prime chemistry linkage is not resolved here',sampling='Fixed then dissociated and barcoded/pooled;18,082-gene original catalog; do not assume assay equivalence to dissociated fresh/thawed3-prime cohorts',source='https://pmc.ncbi.nlm.nih.gov/articles/PMC12236662/',accession='Original local frozen study object'),
      dict(study='walsh',category='Organoid: ventral DIV75; dorsal DIV75 is a separate regional control',capture="10x Chromium 3-prime v3.1 whole cells",sampling='One MEL1 ventral sample and one dorsal sample in this frozen object; inherited stress-cell removal; no independent ventral replicate here',source='https://pmc.ncbi.nlm.nih.gov/articles/PMC12447774/',accession='GSE250482'),
      dict(study='bershteyn_2023',category='2D/adherent manufactured MGE-pIN differentiation, D0/D14 and DIV42 EOP',capture='10x 3-prime v3 whole cells after thaw',sampling='Fourteen source aliquots; sorted/unsorted and early-stage samples; related lots are not independent protocols',source='https://pmc.ncbi.nlm.nih.gov/articles/PMC10993865/',accession='GSE208672'),
      dict(study='bershteyn_2025',category='Manufactured MGE-pIN EOP, 14 sorted pretransplant lots',capture='10x 3-prime v3 whole cells after thaw',sampling='ERBB4-selected pretransplant products; not xenograft nuclei from the companion GSE283776; same protocol family as2023',source='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283775',accession='GSE283775'),
      dict(study='siebert_2026',category='Repository-registered organoid object; source-protocol provenance incomplete',capture='Exact object chemistry not verified; parent NeMO collection describes multiome nuclei',sampling='Old1/2 andYoung1/2, four cellLine values; numeric ages and collection-to-object assay link unresolved',source='https://assets.nemoarchive.org/collection/nemo%3Adat-htzat9t',accession='nemo:dat-htzat9t'),
    ]
    if 'samarasinghe_2021' in STUDIES:
        details.append(dict(study='samarasinghe_2021',category='Cx+GE organoid controls; D56 unfused andD70/D100docked',capture='10xChromium SingleCell3-primev3 whole cells',sampling='Official processed49,942-cell object read once onSlurm; GenotypeCtrl26,935 mapped; Rett23,007 preserved but excluded from primary control comparison. No per-cell ventral/dorsal label available.',source='https://zenodo.org/records/5732813',accession='GSE165577; Zenodo5732813'))
    pd.DataFrame(details).to_csv(OUT/'tables/study_protocol_provenance.tsv',sep='\t',index=False)
    spare=[('shi_2019','results/python_anndata/shi_2019_paper_qc.h5ad','Existing human developmental GE reference; not an in-vitro protocol-frequency comparator'),('xiang_2018','results/xiang_2018/xiang_2018_seurat.rds','Source default configuration excludes Xiang pending validated sample biology; no fresh sample reinterpretation undertaken'),('samarasinghe_2021','results/samarasinghe_2021_zenodo_processed_object/samarasinghe_2021_zenodo_seurat.rds','Bounded read-only selected-gene export added through usable R/Seurat; Genotype Ctrl primary cohort' if 'samarasinghe_2021' in STUDIES else 'Bounded R/Seurat selected-gene export pending' )]
    pd.DataFrame([dict(study=s,path=str(PROJECT/p),exists=(PROJECT/p).exists(),size_bytes=(PROJECT/p).stat().st_size if (PROJECT/p).exists() else 0,reason=r) for s,p,r in spare if s not in STUDIES]).to_csv(OUT/'tables/additional_local_datasets_not_mapped.tsv',sep='\t',index=False)
    # Independent source snapshots and a direct check that frozen caches reproduced scores.
    original=pd.read_csv(OLD/'cells.tsv.gz',sep='\t',usecols=['cell_id','mge_score','maturation_score','early_pv_score'])
    c=table('mgeo_div30_mapped_cells.tsv.gz');assert np.array_equal(original.cell_id,c.cell_id)
    err=max(abs(original.mge_score-c.mge_score).max(),abs(original.maturation_score-c.maturation_score).max());assert err<1e-5
    reference_space=np.load(MAP/'frozen_reference_space.npz');gx=json.loads((OUT/'cache/genes.json').read_text());mx=np.load(OUT/'cache/mgeo_div30_expression.npy',mmap_mode='r')
    pcaerr=float(np.max(abs(transform(mx[:500],gx,reference_space)-reference_space['div30_pcs'][:500])));assert pcaerr<1e-5
    (OUT/'provenance/frozen_mapping_reproduction.json').write_text(json.dumps(dict(n_query_cells=500,maximum_absolute_PC_error=pcaerr,maximum_frozen_context_score_error=float(err),source_target_unchanged=True),indent=2))
    before=json.loads((OUT/'provenance/input_before.json').read_text());unchanged=[]
    for x in before:
        p=Path(x['path']);unchanged.append(p.stat().st_size==x['size_bytes'] and p.stat().st_mtime_ns==x['mtime_ns'])
    assert all(unchanged)
    dup={}
    xx=[]
    for study in ['bershteyn_2023','bershteyn_2025']:
        a=np.load(OUT/'cache'/f'{study}_expression.npy',mmap_mode='r')
        xx.append(set(hashlib.blake2b(row.tobytes(),digest_size=16).digest() for row in a))
    dup['exact_shared_extracted_expression_fingerprints']=len(xx[0]&xx[1]);dup['interpretation']='Exact fingerprints across ~1100 standardized CP10K features suggest reused cells if positive; zero does not exclude related manufacturing lots or different processing. No cells removed.'
    (OUT/'provenance/bershteyn_overlap_audit.json').write_text(json.dumps(dup,indent=2))
    aggregate=freq.query("level=='all'").copy();aggregate['mapping_features_available']=aggregate.study.map(cov)
    slc=genes.query("gene=='SLC6A1'").copy();slc['detection_delta']=slc.mgeo_positive_fraction-slc.external_positive_fraction
    slcsummary=slc.groupby('study').agg(n_samples=('sample','nunique'),mgeo_detection=('mgeo_positive_fraction','mean'),external_detection=('external_positive_fraction','mean'),mean_expression_difference=('mgeo_minus_external_log1p_mean','mean'),min_sample_expression_difference=('mgeo_minus_external_log1p_mean','min'),max_sample_expression_difference=('mgeo_minus_external_log1p_mean','max')).reset_index()
    slcsummary.to_csv(OUT/'tables/matched_slc6a1_summary.tsv',sep='\t',index=False)
    if not noslc.empty:
        q=noslc.query("gene=='SLC6A1'").groupby('study').agg(no_slc_n_samples=('sample','nunique'),no_slc_mgeo_detection=('mgeo_positive_fraction','mean'),no_slc_external_detection=('external_positive_fraction','mean'),no_slc_mean_expression_difference=('mgeo_minus_external_log1p_mean','mean')).reset_index()
        slcsummary=slcsummary.merge(q,on='study',how='left')
    diagpath=OUT/'tables/omit_SLC6A1_subsample_diagnostic.tsv'
    diag=table(diagpath.name) if diagpath.exists() else pd.DataFrame()
    slc_text='The bounded no-SLC6A1 diagnostic is pending.'
    if not diag.empty:
        slc_text=f'A bounded stratified diagnostic omitted SLC6A1 from the reference features in {diag.n_sampled.sum():,} external cells, using3,000 held-sample MGEO queries to recalibrate distance support. It rescued {diag.no_slc_supported_candidate_0_4_n.sum():,} supported≥0.4 target-neighbor candidates and yielded {diag.no_slc_majority_n.sum():,} majority-target-neighbor cells. This is a subsample diagnostic, not a full-cohort frequency estimate. Full external no-SLC6A1 remapping and matched expression contrasts were not run when no candidate was rescued.'
    if not noslc.empty:slc_text='The omit-SLC6A1 diagnostic removes the transporter from reference PCA features, leaves the original target fixed, repeats held-sample calibration and matching, and reduces direct conditioning of the matched axis on SLC6A1.'
    slc_text+=' The original consensus remains mathematically related to SLC6A1 because its frozen mapping/sparse model includes the gene. Low detection is not proof of biological absence.'
    summary={'original_target_prevalence':float(c.target.mean()),'mgeo_held_sample_proxy_prevalence':float(c.majority_target_resemblance.mean()),'mgeo_held_sample_resemblance_auc_range':cal.groupby('sample').auc.first().agg(['min','max']).to_dict(),'external_frequency':aggregate.to_dict('records'),'matching':match.groupby('study').agg(n_external_matched=('matched_external_n','sum'),max_unique_mgeo_per_sample=('unique_mgeo_cells','max')).reset_index().to_dict('records'),'matched_program_rankings':rank.to_dict('records'),'matched_slc6a1':slcsummary.to_dict('records'),'omit_slc6a1_bounded_diagnostic':diag.to_dict('records'),'source_preserved':all(unchanged),'frozen_context_maximum_absolute_error':float(err),'unusual_abundance_established':False,'causal_program_unique_to_mgeo_established':False,'sources':details,'cell_fingerprint_overlap':dup}
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x)))
    # Figure1: distributions rather than translating a cohort-specific target threshold.
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(2,2,figsize=(10.6,7.6),layout='constrained')
    allc={'mgeo_div30':c}
    for study in STUDIES:allc[study]=table(f'{study}_mapped_cells.tsv.gz')
    for study,z in allc.items():
        a=np.sort(z.neighbor_target_resemblance.to_numpy());axs[0,0].plot(a,np.arange(1,len(a)+1)/len(a),color=COLORS[study],label=LABELS[study])
    axs[0,0].axvline(.5,color='grey',ls='--',lw=.8);axs[0,0].set(xlabel='Fraction of DIV30 neighbors in original target',ylabel='Cumulative cell fraction',title='A  Full developmental-resemblance distributions');axs[0,0].legend(fontsize=6,loc='upper left')
    for j,study in enumerate(['mgeo_div30']+STUDIES):
        z=allc[study];rows=[]
        for sample,g in z.groupby('sample'):
            raw=(g.neighbor_target_resemblance>=.5).mean();support=((g.neighbor_target_resemblance>=.5)&g.in_div30_distance_support).mean();rows.append((raw,support))
        offs=np.linspace(-.22,.22,len(rows))
        for off,(raw,support) in zip(offs,rows):
            axs[0,1].plot([j+off]*2,[support,raw],color=COLORS[study],lw=.6);axs[0,1].scatter(j+off,raw,s=15,facecolors='white',edgecolors=COLORS[study]);axs[0,1].scatter(j+off,support,s=9,color=COLORS[study])
    axs[0,1].set(xticks=range(len(STUDIES)+1),xticklabels=[LABELS[x] for x in ['mgeo_div30']+STUDIES],ylabel='Fraction of all assayed cells',title='B  Sample frequency: open=all; filled=distance-supported');axs[0,1].tick_params(axis='x',rotation=25,labelsize=6);axs[0,1].yaxis.set_major_formatter(PercentFormatter(1))
    for study in STUDIES:
        z=freq.query("level=='sample' and study==@study");axs[1,0].scatter(z.median_counts,z.in_div30_distance_support_fraction,label=LABELS[study],color=COLORS[study],s=22)
    axs[1,0].set(xscale='log',xlabel='Median raw library counts per sample',ylabel='Within DIV30 95th-percentile distance support',title='C  Sampling and mapping support');axs[1,0].yaxis.set_major_formatter(PercentFormatter(1));axs[1,0].legend(fontsize=6)
    for study in STUDIES:
        z=freq.query("level=='sample' and study==@study");axs[1,1].scatter(z.mean_target_neighbor_resemblance,z.div90_out_of_reference_fraction,label=LABELS[study],color=COLORS[study],s=22)
    axs[1,1].set(xlabel='Mean DIV30 target-neighbor resemblance',ylabel='Outside frozen DIV90 distance support',title='D  Later reference is a narrow cortical population');axs[1,1].yaxis.set_major_formatter(PercentFormatter(1))
    fig.suptitle('External mapping tests developmental resemblance and reference support',fontsize=12)
    save(fig,'External_F1_reference_resemblance')
    # Figure2: matching balance, per-study program panels and SLC6A1 sensitivity.
    fig,axs=plt.subplots(1,3,figsize=(13,5.2),layout='constrained',gridspec_kw={'width_ratios':[1.1,1.25,1]})
    piv=rank.pivot(index='panel',columns='study',values='sample_mean_effect').reindex(PROGRAMS).reindex(columns=STUDIES).astype(float)
    im=axs[0].imshow(piv,vmin=-.7,vmax=.7,cmap='RdBu_r',aspect='auto');axs[0].set(xticks=range(len(piv.columns)),xticklabels=[LABELS[z] for z in piv.columns],yticks=range(len(piv)),yticklabels=[z.replace('_',' ') for z in piv.index],title='A  Matched program expression');axs[0].tick_params(axis='x',rotation=45,labelsize=6)
    cb=fig.colorbar(im,ax=axs[0],label='MGEO minus external, DIV30-scaled panel mean',shrink=.6)
    metrics=['mge_score','maturation_score','postmitotic_continuous','neighbor_target_resemblance','unmatched_QC_log1p_total_counts','unmatched_QC_log1p_n_genes']
    for j,metric in enumerate(metrics):
        z=bal[bal.metric==metric]
        for study,g in z.groupby('study'):
            off=(STUDIES.index(study)-(len(STUDIES)-1)/2)*.12
            axs[1].scatter(g.after_standardized_mean_difference,j+off,s=13,color=COLORS[study],alpha=.7)
    axs[1].axvline(0,color='grey',lw=.7);axs[1].axvline(.1,color='grey',ls=':',lw=.7);axs[1].axvline(-.1,color='grey',ls=':',lw=.7)
    axs[1].set(yticks=range(6),yticklabels=['MGE','maturity','postmitotic','resemblance','counts (unmatched)','genes (unmatched)'],xlabel='Post-match standardized mean difference\n(external minus matched MGEO)',title='B  Balance and remaining assay differences');axs[1].invert_yaxis()
    for j,study in enumerate(STUDIES):
        z=slc[slc.study==study];axs[2].scatter(z.mgeo_minus_external_log1p_mean,np.full(len(z),j)-.10,s=18,color=COLORS[study])
        if not noslc.empty:
            z=noslc.query("gene=='SLC6A1' and study==@study");axs[2].scatter(z.mgeo_minus_external_log1p_mean,np.full(len(z),j)+.10,s=18,facecolors='white',edgecolors=COLORS[study])
    axs[2].axvline(0,color='grey',ls='--',lw=.8);axs[2].set(yticks=range(len(STUDIES)),yticklabels=[LABELS[z] for z in STUDIES],xlabel='Matched SLC6A1 mean log1p difference\n(MGEO minus external)',title='C  Filled=frozen mapping; open=omit SLC6A1');axs[2].invert_yaxis()
    if rank.empty:
        cb.remove()
        for ax in axs:ax.clear()
        axs[0].axis('off');axs[0].text(.03,.83,'Program contrasts not estimable',fontsize=12,weight='bold');axs[0].text(.03,.66,'No sample supplied ≥30 cells meeting\nreference-distance support, antecedent\nresemblance, exact postmitotic status,\nand every state-axis caliper.',fontsize=10,va='top')
        nums=match.groupby('study').matched_external_n.sum().reindex(STUDIES).fillna(0)
        axs[1].barh(range(len(STUDIES)),nums,color=[COLORS[s] for s in STUDIES]);axs[1].set(yticks=range(len(STUDIES)),yticklabels=[LABELS[s] for s in STUDIES],xlabel='Matched external cells across samples',title='Declared support and matching criteria')
        if nums.max()==0:
            axs[1].set(xlim=(0,1),xticks=[0,1])
            for j in range(len(STUDIES)):axs[1].text(.025,j,'0',va='center',fontsize=10)
        axs[2].axis('off');axs[2].text(.05,.8,'Unresolved comparison',fontsize=12,weight='bold');axs[2].text(.05,.63,'Low or zero mapping support cannot\nestablish biological absence.\n\nSLC6A1 and pathway differences\nrequire a validated cross-study\nrepresentation or matched experiments.',fontsize=10,va='top')
    fig.suptitle('Matched-state contrasts remain descriptive and assay-confounded' if not rank.empty else 'Cross-study matching did not establish comparable states',fontsize=12);save(fig,'External_F2_matched_programs')
    # Concise source-grounded synthesis for parent report.
    choose=['study','n_cells','mapping_features_available','candidate_fraction_0.5','in_div30_distance_support_fraction','supported_candidate_fraction_0.5','div90_out_of_reference_fraction']
    lines=['# External comparison: developmental resemblance, not prospective fate','',
      '**Unusually abundant in MGEO-SOSR is not established.** External cells were tested for resemblance to the original DIV30 target. A target-neighbor majority is only a proxy: the original target is a within-MGEO consensus rank, and its 20% frequency cannot be transferred as a literal external identity label. No raw GRIA2/OPCML/NOTCH1 rule was applied across studies.',
      f'The proxy labels {c.majority_target_resemblance.mean():.1%} of MGEO DIV30 cells when all neighbors are held out by sample, versus the frozen {c.target.mean():.1%} original target. Held-sample discrimination AUC ranges {cal.auc.min():.3f}–{cal.auc.max():.3f}; this is a same-assay calibration, not cross-study validation. Threshold sensitivities0.4/0.5/0.6 are tabulated. Fractions below are fractions of each assayed cohort, not differentiation yield.',
      '', 'DIV30 distance support was '+ '; '.join(f'{LABELS[r.study]}: {r.in_div30_distance_support_fraction:.3%}' for r in aggregate.itertuples())+'. These low fractions limit cross-study interpretation before any target-state claim.',md(aggregate[choose]),'',
      'MGEO tissue was fixed before dissociation/barcoding using the10x Fixed RNA Sample Preparation Kit; the paper later calls the libraries v.3. This does not establish assay equivalence with the external3-prime whole-cell preparations. The exact fixed/probe-versus3-prime protocol linkage is unresolved here. Normalization itself was verified: canonical inputs declare raw counts, extracted log1pCP10K exactly reproduces those counts, and available source lognorm layers agree within4.77e-7. [MGEO methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC12236662/)',
      'Walsh has one ventral and one dorsal DIV75 sample in this cache; the dorsal sample is shown separately and excluded from matched-state program contrasts. Its published assembloid system produces PVALB-positive fast-spiking cells, so it is not a PV-protein-negative control. [Walsh et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC12447774/)',
      'Bershteyn datasets represent the same manufactured MGE-pIN protocol family: early cells and sorted/unsorted DIV42 products in2023, sorted EOP products in2025. The protocol uses NOTCH/CDK inhibition to promote postmitotic differentiation and ERBB4 positive selection. Their selected-product frequencies are not directly comparable to whole MGEOs, and related lots/papers do not supply independent protocol replication. [Bershteyn et al.2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10993865/), [GSE283775](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283775).',
      'Siebert is a repository-registered organoid object, but the available NeMO parent metadata describes multiome nuclei and does not resolve the exact assay, ages, or provenance linkage of this Old/Young object. Its numerical projection is exploratory; it cannot establish an absence or a protocol ranking. [NeMO collection](https://assets.nemoarchive.org/collection/nemo%3Adat-htzat9t).',
      'The Samarasinghe official processed object was exported in a bounded read-only Slurm job after the localRDS read exceeded memory. The primary comparison retains26,935 Ctrl cells atD56(unfused),D70(docked),andD100(docked);23,007 Rett cells remain in the source export. This is a Cx+GE organoid/assembloid mixture, and no per-cell regional label is supplied; a protocol-wide MGE production frequency cannot be recovered from that denominator. Gene IDs come from the official object, never anonymous rawCSV rows. [Official processed object](https://zenodo.org/records/5732813), [GEO study](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165577).' if 'samarasinghe_2021' in STUDIES else '',
      '', '**Matched biological programs:** external candidates first required≥0.4 target-neighbor resemblance and distance within the95th percentile of held-sample MGEO mapping. Each was paired to one original MGEO target cell with replacement, requiring every MGE, maturation, postmitotic-continuous, and resemblance axis to be within0.5 target SD, with exact matching on a harmonized postmitotic proxy (maturation>progenitor and cycle_score<0.5). This common proxy does not replace the original phase-derived postmitotic flag or target. Both before/after balance and unbalanced RNA-complexity axes are reported. Repeated MGEO matches and cells within a sample are not independent replicates.', '',md(match.groupby('study').agg(n_sample_groups=('sample','size'),supported_candidates=('candidate_supported_n','sum'),matched_external=('matched_external_n','sum'),max_unique_mgeo_per_sample=('unique_mgeo_cells','max')).reset_index()),'',
      'Program panels were declared before contrasts and exclude the original PV target genes and three gate markers. These are ranked operational expression panels, not formal pathway-enrichment tests. Effects use the same available genes in each matched pair, standardized to all DIV30 means/SD and averaged per external sample. Positive effects mean higher in matched MGEO. The panels cannot distinguish biological pathway activity from library, chemistry, dissociation, cryopreservation, prior sorting, gene-catalog, and unmatched maturation differences.', '',md(rank.sort_values(['study','within_study_rank'])[['study','panel','n_samples','sample_mean_effect','sample_min_effect','sample_max_effect']]),'',
      '**SLC6A1 after state matching:**', '',md(slcsummary),'',
      slc_text,
      '', '**Why PV protein appears later:** these cross-sectional contrasts do not resolve the mechanism. Low external matching support cannot establish whether the same state is present, absent, or unusually abundant in any protocol. Any surviving MGEO expression differences are candidates for follow-up, not unique causal explanations of later PV protein. A prospective retained/removed comparison remains necessary to distinguish fate enrichment from a generic postmitotic state.',
      '', '**Preservation and scope:** external inputs and original MGEO H5AD retain their recorded size and modification time; frozen-context reproduction passed. Extra local RDS datasets and the developmental Shi reference are inventoried with selection reasons. No large pipeline, gate search, original target change, or new dataset download was performed.']
    (OUT/'findings.md').write_text('\n\n'.join(lines))
    if rank.empty:
        text=(OUT/'findings.md').read_text()
        text=text.replace('**Matched biological programs:**','**Matched biological programs: not estimable under the declared criteria.** No external sample supplied at least30 supported, status-matched, caliper-matched candidate cells. The empty program tables represent a failed comparability check, not absent pathways. No whole-study expression contrast was substituted.\n\nMatching attempted:')
        text=text.replace('Each was paired to one original MGEO target cell with replacement','Each eligible candidate would be paired to one original MGEO target cell with replacement')
        text=text.replace('Both before/after balance and unbalanced RNA-complexity axes are reported.','No qualifying matched set exists, so before/after balance and RNA-complexity differences cannot be estimated.')
        text=text.replace('These are ranked operational expression panels, not formal pathway-enrichment tests. Effects use the same available genes in each matched pair, standardized to all DIV30 means/SD and averaged per external sample. Positive effects mean higher in matched MGEO.','These are predeclared operational expression panels. No expression effect or ranking could be estimated because the comparability check failed.')
        text=text.replace(md(rank.sort_values(['study','within_study_rank'])[['study','panel','n_samples','sample_mean_effect','sample_min_effect','sample_max_effect']]),'No matched program rankings are estimable.')
        text=text.replace(md(slcsummary),'No matched SLC6A1 contrast is estimable.')
        (OUT/'findings.md').write_text(text)
    captions='''# External figures and tables

External_F1_reference_resemblance: All query cells projected into the same reconstructed frozen 30-PC DIV90 space;1000 mapping features (999 in the separate SLC6A1 sensitivity), with missing features imputed to reference means. A shows full 30-neighbor target-resemblance distributions, not probabilities of later fate. B shows sample-specific majority-target-neighbor frequency (open) and the fraction also within DIV30 held-sample distance support (filled); MGEO benchmark excludes same-sample neighbors. C–D show sequencing-depth and out-of-reference diagnostics. Studies differ in stage, purification, chemistry and regional sampling. MGEO's20% operational target and this proxy have different calibration. No raw gate transferred.

External_F2_matched_programs: Matching uses four state axes with0.5-SD calipers and exact harmonized postmitotic-proxy status,1 nearest MGEO target per external candidate with replacement, and≥30 matched pairs per sample for expression summaries. A displays equal-external-sample mean expression-panel differences, with no inferential cell-level p-values. B shows remaining axis imbalance and unadjusted RNA complexity; nonzero QC imbalance limits interpretation. C shows per-sample SLC6A1 differences with original mapping and independently reconstructed no-SLC6A1 feature-space sensitivity, holding original targets fixed. These are descriptive transcriptomic contrasts, not proof of unique biochemical pathway activity or causes of later PV protein.

local_input_inventory/gene_coverage: filepaths, dimensions, input timestamps, and per-gene measurable-catalog coverage. Missing genes are not biological zeros; the frozen TUBB3 contribution is a known zero placeholder preserved across studies.

external_reference_resemblance_frequency: All-study and sample rows; raw and distance-supported fractions for0.4/0.5/0.6 target-neighbor fractions. Stage/selection differences mean frequencies are assayed-cohort proportions, not production yields.

mgeo_held_sample_resemblance_calibration: operational target agreement when each query's entire sample is excluded from neighbors, including AUC, precision/recovery and proxy prevalence. No holdout external assay exists.

matched_state_overlap/balance: candidate counts, caliper overlap, unique matched MGEO counts, and standardized differences. No cell or source-sample exclusion is hidden; small or absent match groups do not become absence claims.

matched_program_rankings/effects_by_sample: Predeclared panels, per-sample means/ranges and descriptive ranks, no statistical significance implied. matched_genes_by_sample includes individual measured genes and detection fractions; unknown low-sensitivity measurements cannot support absence claims.
'''
    (OUT/'FIGURE_CAPTIONS.md').write_text(captions)
    if rank.empty:
        p=OUT/'FIGURE_CAPTIONS.md';t=p.read_text();start=t.index('External_F2_matched_programs:');stop=t.index('\n\nlocal_input_inventory',start)
        t=t[:start]+('External_F2_matched_programs: No external sample provided a supported candidate set meeting the predeclared resemblance, exact harmonized postmitotic-status, and0.5-SD state-axis calipers. The figure displays zero eligible matched cells and explains why program/SLC6A1 contrasts cannot be estimated. The empty analysis tables do not show absence of a state, gene, or pathway. No whole-study marker contrast substitutes for missing matched-state evidence.')+t[stop:];p.write_text(t)
    for name in ['findings.md','FIGURE_CAPTIONS.md']:
        p=OUT/name;text=p.read_text()
        replacements={'sensitivities0.4':'sensitivities 0.4','the10x':'the 10x','external3-prime':'external 3-prime','probe-versus3-prime':'probe-versus-3-prime','within4.77':'within 4.77','in2023':'in 2023','in2025':'in 2025','al.2023':'al. 2023','least30':'least 30','required≥':'required ≥','the95th':'the 95th','within0.5':'within 0.5','uses four state axes with0.5':'uses four state axes with 0.5','using3,000':'using 3,000','supported≥':'supported ≥','majority-target-neighbor':'≥50%-target-neighbor','target-neighbor majority':'target-neighbor fraction of at least 50%','retains26,935':'retains 26,935','atD56':'at D56','D56(unfused),D70(docked),andD100(docked)':'D56 (unfused), D70 (docked), and D100 (docked)',';23,007':'; 23,007','localRDS':'local RDS','rawCSV':'raw CSV','30-PC DIV90 space;1000':'30-PC DIV90 space; 1,000','(999 in':'(999 in','MGEO\'s20%':'MGEO’s 20%','calipers,1':'calipers, 1','and≥30':'and ≥30','for0.4':'for 0.4','per sample for expression summaries':'per sample before estimating expression summaries'}
        for a,b in replacements.items():text=text.replace(a,b)
        while '\n\n\n' in text:text=text.replace('\n\n\n','\n\n')
        p.write_text(text)
    print('External report complete',flush=True)
if __name__=='__main__':main()
