"""Audit an existing 38-to-17 partition comparison; never fit or change clusters.

Inputs are a byte-verified snapshot of completed Step 07 tables plus baseline
cell metadata and the saved baseline final pairwise DE audit. Outputs include
all nonzero cell flows, all 291 sibling boundaries, a lost/eroded-boundary subset,
canonical gene contrasts, phase effects, 17 destination summaries and merge-flow
figures. A flow diagram is not a inferred clustering tree: cells can split.

Run the frozen copy with --package. Review thresholds are operational descriptions,
not statistical significance, causal attribution or accepted biological labels.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch

REGION = ['MGE','MGE_interneuron','LGE','CGE','POA','basal_forebrain_alternative']
DEVELOPMENT = ['rg_stemness','apical_RG_like','basal_RG_like','ipc_neurogenic','neuroblast','immature_inhibitory','early_neuronal_maturation','later_neuronal_maturation']
NONNEURAL = ['endothelial','pericyte','erythroid']


def digest(path):
    """Stream a file hash to verify the exact evidence snapshot."""
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for x in iter(lambda:f.read(8*1024*1024),b''):h.update(x)
    return h.hexdigest()


def read(root,name,index=False):
    """Read a saved TSV, optionally using its first column as a stable ID."""
    return pd.read_csv(root/name,sep='\t',index_col=0 if index else None)


def save(root,name,data):
    """Save a source table with explicit columns and no implicit row number."""
    data.to_csv(root/name,sep='\t',index=False)


def truth(series):
    """Decode saved Boolean fields without treating the string False as true."""
    return series.astype(str).str.lower().eq('true')


def effect_summary(frame,names,threshold):
    """Return largest all-cell effect and retained phase evidence for an axis."""
    f=frame[frame.program.isin(names)]
    allf=f[f.phase.eq('all')]
    top=allf.loc[allf.standardized_effect.abs().idxmax()] if len(allf) else None
    keep=f[~f.phase.eq('all') & truth(f.eligible) & truth(f.same_direction_as_unstratified)
           & truth(f.unstratified_distinguishing) & f.standardized_effect.abs().ge(threshold)]
    return (abs(float(top.standardized_effect)) if top is not None else 0,
            str(top.program) if top is not None else '',
            int(keep.phase.nunique()),';'.join(sorted(keep.program.unique())))


def flow_plot(counts,path):
    """Draw all exact flows by parent; node height is cells, edges retain splits."""
    with PdfPages(path/'merge_tree_38_to_17.pdf') as pdf:
        for parent in sorted({x.split('.')[0] for x in counts.index}):
            c=counts.loc[counts.index.str.startswith(parent),counts.columns.str.startswith(parent)]
            fig,ax=plt.subplots(figsize=(15,max(7,len(c)*.43)))
            total=c.values.sum(); gap=total*.015
            def nodes(sizes,x):
                pos={};cursor=0
                for name,n in sizes.items():
                    pos[name]=[cursor,cursor+n];ax.plot([x,x],[cursor,cursor+n],lw=7,c='#333333')
                    ax.text(x+(-.035 if x==0 else .035),cursor+n/2,'%s  n=%d'%(name,n),ha='right' if x==0 else 'left',va='center',fontsize=9)
                    cursor+=n+gap
                return pos,cursor
            left,height=nodes(c.sum(axis=1),0);right,h2=nodes(c.sum(axis=0),1)
            lc={k:v[0] for k,v in left.items()};rc={k:v[0] for k,v in right.items()}
            for i,a in enumerate(c.index):
                for b in c.columns:
                    n=int(c.loc[a,b])
                    if not n:continue
                    y,z=lc[a],rc[b];lc[a]+=n;rc[b]+=n
                    verts=[(0,y),(.45,y),(.55,z),(1,z),(1,z+n),(.55,z+n),(.45,y+n),(0,y+n),(0,y)]
                    codes=[MPath.MOVETO,MPath.CURVE4,MPath.CURVE4,MPath.CURVE4,MPath.LINETO,MPath.CURVE4,MPath.CURVE4,MPath.CURVE4,MPath.CLOSEPOLY]
                    ax.add_patch(PathPatch(MPath(verts,codes),facecolor=plt.cm.tab20(i%20),alpha=.55,lw=0))
            ax.set(xlim=(-.37,1.37),ylim=(max(height,h2),-gap*3));ax.axis('off')
            ax.set_title(parent+' | Original fine clusters → existing comparison clusters\nAll cell flows; branching shows splits, not a literal agglomerative merge history',fontsize=13)
            fig.tight_layout();pdf.savefig(fig);fig.savefig(path/(parent+'_merge_tree.png'),dpi=150);plt.close(fig)


def run(package):
    """Verify source bytes, join evidence, classify boundary associations and export."""
    start=time.time();root=package/'inputs';out=package/'outputs';out.mkdir()
    cfg=json.loads((package/'review_config.json').read_text())
    manifest=read(package,'input_manifest.tsv')
    for row in manifest.itertuples():
        if digest(root/row.file)!=row.sha256:raise ValueError('Changed snapshot: '+row.file)
    print('INPUT_HASHES_PASS',flush=True)
    counts=read(root,'counts.tsv',True).astype(int)
    obs=read(root,'baseline_cells.tsv',True)
    dest=read(root,'assignments.tsv.gz',True).iloc[:,0]
    assert obs.index.is_unique and dest.index.is_unique and set(obs.index)==set(dest.index)
    dest=dest.reindex(obs.index);labels=obs.hicat_fine_baseline
    actual=pd.crosstab(labels,dest).reindex(index=counts.index,columns=counts.columns,fill_value=0)
    assert np.array_equal(actual.values,counts.values) and counts.shape==(38,17) and counts.values.sum()==12000
    assert all(a.split('.')[0]==b.split('.')[0] for a,b in zip(labels,dest))
    sizes=counts.sum(axis=1);dsize=counts.sum(axis=0);fractions=counts.div(sizes,axis=0)
    flow=[]
    for a,b in itertools.product(counts.index,counts.columns):
        n=int(counts.loc[a,b])
        if n:flow.append(dict(original_cluster=a,comparison_cluster=b,n_cells=n,original_n=int(sizes[a]),comparison_n=int(dsize[b]),fraction_of_original=n/sizes[a],fraction_of_comparison=n/dsize[b]))
    flows=pd.DataFrame(flow);save(out,'all_cell_contributions.tsv',flows)
    save(out,'cell_memberships.tsv.gz',pd.DataFrame({'cell_id':obs.index,'original':labels.values,'comparison':dest.values,'sample':obs.technical_sample_id.values}))
    phase=read(root,'phase_pairs.tsv').set_index(['cluster_a','cluster_b'])
    effects=read(root,'phase_effects.tsv');effectgroups=dict(tuple(effects.groupby(['cluster_a','cluster_b'])))
    seed=read(root,'seed.tsv',True);sample=read(root,'sample_fractions.tsv',True)
    assert np.allclose(sample.sum(axis=1),1), 'Sample fractions must sum to one within each cluster'
    review=read(root,'cluster_review.tsv').set_index('fine_cluster')
    program_means=read(root,'program_means.tsv',True)
    de=read(root,'baseline_de.tsv');de=de.rename(columns={de.columns[0]:'a',de.columns[1]:'b'}).set_index(['a','b'])
    coverage=read(root,'gene_coverage.tsv');pconf=json.loads((root/'program_config.json').read_text())
    sets={x['name']:set(coverage.loc[coverage.program.eq(x['name']),'measured_symbol'].dropna()) for x in pconf['programs']}
    genes=read(root,'canonical_expression.tsv.gz',True).reindex(obs.index)
    assert genes.notna().all().all()
    means=genes.groupby(labels).mean();detect=genes.gt(0).groupby(labels).mean()
    save(out,'canonical_cluster_means.tsv',means.reset_index());save(out,'canonical_cluster_detection.tsv',detect.reset_index())
    symbolmap=dict(zip(coverage.feature_id,coverage.measured_symbol))
    cycle=set.union(sets['S_phase'],sets['G2M'])
    regiongenes=set.union(*(sets[n] for n in REGION)) - {'Gad1','Gad2','Arx','Dlx1','Dlx2','Dlx5','Dlx6','Sp9','Prox1'}
    devgenes=set.union(*(sets[n] for n in DEVELOPMENT))-cycle
    nongenes=set.union(*(sets[n] for n in NONNEURAL))
    rows=[];contrasts=[]
    for a,b in itertools.combinations(counts.index,2):
        if a.split('.')[0]!=b.split('.')[0]:continue
        p=phase.loc[(a,b)];ef=effectgroups[(a,b)]
        shared=counts.columns[(counts.loc[a]>0)&(counts.loc[b]>0)]
        substantial=[d for d in shared if min(fractions.loc[a,d],fractions.loc[b,d])>cfg['partial_source_fraction']]
        major=[d for d in shared if min(fractions.loc[a,d],fractions.loc[b,d])>=cfg['majority_source_fraction']]
        near=[d for d in shared if min(fractions.loc[a,d],fractions.loc[b,d])>=cfg['near_complete_source_fraction']]
        extent='near_complete' if near else 'majority' if major else 'partial' if substantial else 'minor_leakage' if len(shared) else 'retained'
        pair=a+'__'+b
        reg,regtop,regph,regkeep=effect_summary(ef,REGION,cfg['effect_sd'])
        dev,devtop,devph,devkeep=effect_summary(ef,DEVELOPMENT,cfg['effect_sd'])
        non,nontop,nonph,nonkeep=effect_summary(ef,NONNEURAL,cfg['effect_sd'])
        delta=means.loc[a]-means.loc[b];dd=detect.loc[a]-detect.loc[b]
        coherent=delta.index[(delta.abs()>=cfg['marker_mean_ln1pcpm_difference'])&(dd.abs()>=cfg['marker_detection_difference'])&(np.sign(delta)==np.sign(dd))]
        for gene in delta.index:
            contrasts.append(dict(boundary=pair,cluster_a=a,cluster_b=b,gene=gene,mean_a=means.loc[a,gene],mean_b=means.loc[b,gene],detection_a=detect.loc[a,gene],detection_b=detect.loc[b,gene],mean_difference_a_minus_b=delta[gene],detection_difference_a_minus_b=dd[gene],directionally_concordant_marker=gene in coherent))
        regmarkers=sorted(set(coherent)&regiongenes);devmarkers=sorted(set(coherent)&devgenes);nonmarkers=sorted(set(coherent)&nongenes)
        derow=de.loc[(a,b)] if (a,b) in de.index else de.loc[(b,a)]
        deids=set(json.loads(derow.up_genes))|set(json.loads(derow.down_genes));desymbols={symbolmap.get(g,g) for g in deids}
        cyclede=len(desymbols&cycle)/max(1,len(deids));sampletv=float((sample.loc[a]-sample.loc[b]).abs().sum()/2)
        cycletv=float(p.cell_cycle_total_variation)
        # A pair with no initial distinguishing program can have adequately
        # sampled phases even though the persistence summary counts zero
        # eligible *distinguishing* phases. Preserve that distinction.
        eligible=int(ef.loc[~ef.phase.eq('all') & truth(ef.eligible),'phase'].nunique())
        regional=reg>=cfg['effect_sd'] and regph>=1 and len(regmarkers)>=2
        developmental=dev>=cfg['effect_sd'] and devph>=1 and len(devmarkers)>=2
        # Require a coherent single non-neural program, not a mixture of an
        # erythroid gene, an endothelial gene and a pericyte gene. Two strongly
        # endothelial groups are not a neural-versus-vascular boundary.
        nonneural_programs=[]
        for program in NONNEURAL:
            hi,lo=(a,b) if program_means.loc[a,program]>program_means.loc[b,program] else (b,a)
            markers=[g for g in set(coherent)&sets[program] if means.loc[hi,g]>means.loc[lo,g] and detect.loc[hi,g]>=.3]
            if program_means.loc[hi,program]>=1 and program_means.loc[lo,program]<.5 and len(markers)>=3:
                nonneural_programs.append(program)
        nonneural=bool(nonneural_programs)
        opposing_region=((program_means.loc[a,'LGE']>.2 and program_means.loc[b,'LGE']<0 and program_means.loc[b,'MGE']>.3)
                         or (program_means.loc[b,'LGE']>.2 and program_means.loc[a,'LGE']<0 and program_means.loc[a,'MGE']>.3))
        # Differences in MGE-lineage program strength alone can reflect
        # maturation, not different regional origins. Preserve the contrast
        # but do not classify it as a compelling regional boundary by default.
        regional_identity=regional and opposing_region
        ipc_transition=(max(program_means.loc[a,'ipc_neurogenic'],program_means.loc[b,'ipc_neurogenic'])>=1
                        and min(program_means.loc[a,'ipc_neurogenic'],program_means.loc[b,'ipc_neurogenic'])<=0)
        compelling=nonneural or (regional_identity and regph>=2) or (developmental and devph>=2 and ipc_transition)
        # Axes are non-exclusive. Never turn insufficient phase coverage into
        # evidence that biology disappears, or sample enrichment into batch.
        if nonneural:axis='regional_lineage_identity';sub='coherent_'+','.join(nonneural_programs)+'_versus_other_compartment'
        elif regional_identity:axis='regional_lineage_identity';sub=regtop
        elif developmental:axis='developmental_maturation';sub=devtop
        elif regional:axis='weak_unresolved_transcriptional_structure';sub='lineage_program_strength_difference_origin_vs_maturation_unresolved'
        elif cycletv>=cfg['cycle_tv'] and eligible>0 and int(p.n_retained_phases)==0:axis='cell_cycle';sub='cycle_associated_without_retained_identity_at_declared_threshold'
        elif sampletv>=cfg['sample_tv'] and (bool(review.loc[a,'sample_concern']) or bool(review.loc[b,'sample_concern'])):axis='sample_composition';sub='descriptive_sample_association_not_batch'
        else:axis='weak_unresolved_transcriptional_structure';sub='limited_phase_coverage' if eligible<2 else 'no_dominant_independent_axis'
        low=seed.loc[a,'stability']=='LOW_STABILITY' and seed.loc[b,'stability']=='LOW_STABILITY'
        favorable=(eligible>=1 and int(p.n_unstratified_distinguishing_programs)==0 and not nonneural
                   and (low or min(seed.loc[a,'jaccard'],seed.loc[b,'jaccard'])<.5))
        decision='compelling_biological_boundary_lost' if compelling else 'biological_boundary_at_risk' if (regional_identity or developmental or nonneural) else 'merge_resolves_weak_subdivision' if favorable else 'unresolved'
        rows.append(dict(boundary=pair,cluster_a=a,cluster_b=b,parent=a.split('.')[0],extent=extent,is_lost_or_eroded=extent in ['near_complete','majority','partial'],shared_destinations=';'.join(shared),substantial_shared_destinations=';'.join(substantial),shared_contributions=json.dumps({d:{'a_cells':int(counts.loc[a,d]),'a_fraction':float(fractions.loc[a,d]),'b_cells':int(counts.loc[b,d]),'b_fraction':float(fractions.loc[b,d])} for d in shared}),cross_pair_coassignment_probability=float((fractions.loc[a]*fractions.loc[b]).sum()),predominant_association=axis,association_detail=sub,interpretation=decision,regional_effect_sd=reg,regional_program=regtop,regional_retained_phases=regph,regional_retained_programs=regkeep,developmental_effect_sd=dev,developmental_program=devtop,developmental_retained_phases=devph,developmental_retained_programs=devkeep,nonneural_effect_sd=non,nonneural_retained_phases=nonph,cycle_composition_tv=cycletv,sample_composition_tv=sampletv,n_eligible_phases=eligible,phase_persistence=p.identity_difference_persists_after_cell_cycle_stratification,phase_interpretation=p.interpretation_reason,original_de_score=derow.score,original_de_gene_count=derow['num'],original_de_separation=derow.meets_declared_separation,original_de_standard_cycle_gene_fraction=cyclede,regional_contrast_markers=';'.join(regmarkers),developmental_contrast_markers=';'.join(devmarkers),nonneural_contrast_markers=';'.join(nonmarkers),seed_a=seed.loc[a,'stability'],seed_b=seed.loc[b,'stability'],seed_jaccard_a=seed.loc[a,'jaccard'],seed_jaccard_b=seed.loc[b,'jaccard'],both_seed_unstable=low,sample_a=review.loc[a,'top_sample'],sample_b=review.loc[b,'top_sample'],sample_fraction_a=review.loc[a,'top_sample_fraction'],sample_fraction_b=review.loc[b,'top_sample_fraction'],hypothesis_a=review.loc[a,'provisional_identity'],hypothesis_b=review.loc[b,'provisional_identity']))
    boundaries=pd.DataFrame(rows);lost=boundaries[boundaries.is_lost_or_eroded].copy()
    assert len(boundaries)==291 and boundaries.boundary.is_unique
    save(out,'all_original_sibling_boundaries.tsv',boundaries);save(out,'lost_boundaries.tsv',lost)
    save(out,'canonical_boundary_gene_contrasts.tsv',pd.DataFrame(contrasts))
    save(out,'matched_phase_program_evidence.tsv',effects)
    blocking=lost[lost.interpretation.isin(['compelling_biological_boundary_lost','biological_boundary_at_risk'])]
    save(out,'decision_blocking_boundaries.tsv',blocking)
    save(out,'merges_resolving_weak_subdivisions.tsv',lost[lost.interpretation.eq('merge_resolves_weak_subdivision')])
    destinations=[]
    for b in counts.columns:
        f=flows[flows.comparison_cluster.eq(b)].sort_values('n_cells',ascending=False)
        related=lost[lost.substantial_shared_destinations.str.split(';').apply(lambda x:b in x)]
        destinations.append(dict(comparison_cluster=b,n_cells=int(dsize[b]),n_original_contributors=len(f),contributors='; '.join('%s: %d (%.2f%% destination; %.2f%% source)'%(x.original_cluster,x.n_cells,100*x.fraction_of_comparison,100*x.fraction_of_original) for x in f.itertuples()),compelling_boundaries=';'.join(related.loc[related.interpretation.eq('compelling_biological_boundary_lost'),'boundary']),at_risk_boundaries=';'.join(related.loc[related.interpretation.eq('biological_boundary_at_risk'),'boundary']),favorable_boundaries=';'.join(related.loc[related.interpretation.eq('merge_resolves_weak_subdivision'),'boundary'])))
    ds=pd.DataFrame(destinations);save(out,'comparison_cluster_summary.tsv',ds)
    flow_plot(counts,out)
    fig,ax=plt.subplots(figsize=(13,12));im=ax.imshow(fractions,aspect='auto',vmin=0,vmax=1,cmap='viridis');ax.set_xticks(range(17));ax.set_xticklabels(counts.columns,rotation=90);ax.set_yticks(range(38));ax.set_yticklabels(counts.index);fig.colorbar(im,ax=ax,label='Fraction of original cluster');ax.set_title('Exact membership redistribution: 38 → 17');fig.tight_layout();fig.savefig(out/'membership_heatmap.png',dpi=150);fig.savefig(out/'membership_heatmap.pdf');plt.close(fig)
    summary=dict(status='IN_REVIEW',comparison='still_unresolved',original_cells=12000,original_clusters=38,comparison_clusters=17,nonzero_flows=len(flows),original_clusters_with_multiple_destinations=int((counts.gt(0).sum(axis=1)>1).sum()),boundaries=len(boundaries),lost_or_eroded_boundaries=len(lost),extent_counts=boundaries.extent.value_counts().to_dict(),lost_associations=lost.predominant_association.value_counts().to_dict(),lost_interpretations=lost.interpretation.value_counts().to_dict(),decision_blocking_boundaries=len(blocking),new_clustering=False,parameters_changed=False,seconds=time.time()-start)
    (out/'summary.json').write_text(json.dumps(summary,indent=2))
    lines=['# Existing Step 07: 38 → 17 boundary audit','', '**Disposition: STILL UNRESOLVED; neither partition is final.**', '', 'Dissected E14.5 mouse MGE; same 12,000 pilot cells. No new fit, regression, normalization, DE test, cell removal or annotation assignment. Existing baseline DE evidence is reused. This is a membership comparison, not the chronological merge history inside HiCAT.', '', '## Read first','', '- [Merge-flow tree PDF](merge_tree_38_to_17.pdf): four parent pages, every nonzero flow, width = cell count. Branches expose splits.', '- [17-group contribution table](comparison_cluster_summary.tsv) and [every cell contribution](all_cell_contributions.tsv). Both percentage denominators are explicit.', '- [One row per lost/eroded boundary](lost_boundaries.tsv), with all underlying [291 sibling pairs](all_original_sibling_boundaries.tsv).', '- [Exact decision-blocking boundaries](decision_blocking_boundaries.tsv) and [potentially useful merges](merges_resolving_weak_subdivisions.tsv).','', '## What counts as a lost boundary?','', 'A pair is near-completely coalesced when one destination receives >=90% of each original group; majority-coalesced at >=50% each; partially eroded when a common destination receives >10% of each. Smaller shared flows are minor leakage, retained in the exhaustive table. These overlap cutoffs do not tune clusters. All nonzero contributions are exported, including tiny populations. Multiple lost pairwise boundaries can describe one merger; their number is not 38 minus 17.','', '## Evidence and interpretation rules','', 'Predominant association is a conservative descriptive synthesis, not a causal result. Regional/developmental evidence requires >=0.5 pilot SD program contrast, retained direction in an eligible phase (>=20 cells per group), and >=2 canonical genes differing by >=0.5 ln1pCPM and >=10 percentage points in detection, in concordant directions. Compelling regional contrasts additionally require opposing MGE/LGE program support and >=2 retained phases. Compelling developmental contrasts require an IPC transition (one IPC mean >=1, the other <=0) and >=2 retained phases. A coherent non-neural contrast requires one named program mean >=1 versus <0.5 and >=3 concordant markers from that same program detected in >=30% of cells on the higher side; phase limitations remain explicit. These are operational review flags, not independently validated cell types. Single eligible phases remain limited evidence. The full program effects, canonical gene differences and alternative axes are retained.', '', 'Cell-cycle predominance requires phase-composition TV >=0.4 with evaluable phases but no retained independent identity difference. Sample association uses sample TV >=0.4 plus an existing cluster flag; this is not batch causality. A potentially useful merge requires at least one original group to be seed-low, an adequately sampled phase, and no original >=0.5 SD identity-program contrast. This supports simplifying subdivisions; it does not negate their original DE evidence. Otherwise evidence is left unresolved. A negative screen does not prove absent biology. Standard-cycle overlap in the DE list is a lower bound from the saved canonical reference, not all possible cycle genes.', '', 'Existing original pairwise DE scores establish evidence under the original criteria only. The comparison is a refit with altered DE criteria, not a direct application of the new thresholds to each original pair. This audit cannot prove which specific internal DE decision caused a lost boundary. No new DE tests are run.','', '## Counts','', '```json',json.dumps(summary,indent=2),'```','', '## Every comparison cluster','']
    for row in ds.itertuples():
        lines += ['### '+row.comparison_cluster+' (n='+str(row.n_cells)+')','',row.contributors,'','Compelling boundary flags: '+(row.compelling_boundaries or 'none')+'.','At-risk boundary flags: '+(row.at_risk_boundaries or 'none')+'.','Potentially useful merge flags: '+(row.favorable_boundaries or 'none')+'.','']
    lines+=['## Exact boundaries preventing a decision','','The following flags require biological review; weak seed agreement does not cancel coherent canonical identity evidence. Partial erosion means only subsets combine, not both entire clusters.','']
    for row in blocking.itertuples():lines+=['- **'+row.boundary+'** ['+row.extent+'] in '+row.substantial_shared_destinations+': '+row.predominant_association+'; '+row.interpretation+'; retained programs '+row.regional_retained_programs+' / '+row.developmental_retained_programs+'; original DE genes='+str(row.original_de_gene_count)+'.']
    lines+=['','The 17-group partition can simplify some weak subdivisions, but the listed biological contrasts and split flows prevent accepting it as globally preferable. Coherent population/state boundaries can be lost alongside weak boundaries. The original 38-group partition is also seed-sensitive. **Conclusion: still unresolved; do not declare 17 final.**']
    (out/'README.md').write_text('\n'.join(lines)+'\n')
    checks=dict(input_hashes=True,independent_crosstab_exact=True,cells=12000,unique_cells=True,fixed_parent_nesting=True,all_291_pairs=True,all_17_destinations=True,all_38_sources=True,flow_cell_sum=int(flows.n_cells.sum()),new_fits=0)
    (out/'checks.json').write_text(json.dumps(checks,indent=2))
    save(out,'output_manifest.tsv',pd.DataFrame([dict(file=str(p.relative_to(out)),bytes=p.stat().st_size,sha256=digest(p)) for p in sorted(out.rglob('*')) if p.is_file()]))
    (package/'COMPUTATION_SUCCESS.txt').write_text('IN_REVIEW\n')
    print(json.dumps(summary),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--package',type=Path,required=True)
    run(parser.parse_args().package)
