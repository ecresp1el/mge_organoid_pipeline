"""Concise decision table for saved 38-to-17 cell flows; no HiCAT execution.

Read frozen inputs, evaluate destination-specific contributing cells, retain
original DE/phase evidence, and distinguish boundary reproducibility from whole
cluster stability. Save all flow-conditioned score and marker summaries. Labels,
expression values, score coefficients and clustering parameters never change.
"""
import argparse,hashlib,html,json
from pathlib import Path
import numpy as np
import pandas as pd

REGION=['MGE','MGE_interneuron','LGE','CGE','POA','basal_forebrain_alternative']
DEV=['rg_stemness','apical_RG_like','basal_RG_like','ipc_neurogenic','neuroblast','immature_inhibitory','early_neuronal_maturation','later_neuronal_maturation']


def sha(path):
    """Stream one file's byte identity."""
    h=hashlib.sha256()
    with path.open('rb') as f:
        for x in iter(lambda:f.read(8*1024*1024),b''):h.update(x)
    return h.hexdigest()


def read(root,name,index=False):
    """Read a source table with explicit cell/cluster IDs when requested."""
    return pd.read_csv(root/name,sep='\t',index_col=0 if index else None)


def boundary_seed(a,b):
    """Describe repeat-seed boundary retention, accounting for fragmentation.

    Strong: different dominant destinations each retain >=80% and cross-pair
    coassignment <=0.1. Moderate: different dominant destinations each retain
    >=50% and cross coassignment <=0.2. Otherwise low/unresolved. Low cross
    coassignment alone is not stability: fragmentation can also lower it.
    """
    cross=float((a*b).sum());different=a.idxmax()!=b.idxmax()
    status='strong' if different and min(a.max(),b.max())>=.8 and cross<=.1 else 'moderate' if different and min(a.max(),b.max())>=.5 and cross<=.2 else 'low/unresolved'
    return status,cross


def run(package):
    """Join existing evidence and export the requested 11-column decision table."""
    root=package/'inputs';out=package/'outputs';out.mkdir()
    manifest=read(package,'input_manifest.tsv')
    assert all(sha(root/x.file)==x.sha256 for x in manifest.itertuples())
    d=read(root,'lost.tsv');obs=read(root,'cells.tsv',True);dest=read(root,'assignments.tsv.gz',True).iloc[:,0].reindex(obs.index)
    genes=read(root,'genes.tsv.gz',True).reindex(obs.index);scores=read(root,'scores.tsv.gz',True).reindex(obs.index);identity=read(root,'identity.tsv.gz',True).reindex(obs.index)
    phases=read(root,'cycle.tsv.gz',True).validation_phase.reindex(obs.index)
    phase_names=sorted(phases.dropna().unique())
    assert set(phase_names)=={'G1-like','S','G2/M'}, 'Unexpected stored phase labels'
    original_effects=read(root,'original_phase_effects.tsv').set_index(['cluster_a','cluster_b','phase','program'])
    exact_phase_checks=0
    seed=read(root,'seed_counts.tsv',True);seed=seed.div(seed.sum(axis=1),axis=0)
    assert all(x.index.equals(obs.index) for x in [genes,scores,identity,phases,dest]) and len(obs)==12000
    assert obs.index.is_unique and dest.notna().all() and genes.notna().all().all()
    scale=identity.std(ddof=1);labels=obs.hicat_fine_baseline
    groups={};programrows=[];markerrows=[]
    for (a,b),ix in obs.groupby([labels,dest]).groups.items():
        ix=pd.Index(ix);groups[(a,b)]=ix
        for prog in scores:
            programrows.append(dict(original=a,comparison=b,n=len(ix),program=prog,mean=scores.loc[ix,prog].mean(),median=scores.loc[ix,prog].median()))
        for gene in genes:
            markerrows.append(dict(original=a,comparison=b,n=len(ix),gene=gene,mean_ln1pcpm=genes.loc[ix,gene].mean(),fraction_detected=genes.loc[ix,gene].gt(0).mean()))
    pd.DataFrame(programrows).to_csv(out/'contributing_subset_programs.tsv',sep='\t',index=False)
    pd.DataFrame(markerrows).to_csv(out/'contributing_subset_markers.tsv',sep='\t',index=False)
    rows=[];detail=[];floweffects=[]
    for r in d.itertuples():
        a,b=r.cluster_a,r.cluster_b;state,cross=boundary_seed(seed.loc[a],seed.loc[b]);shared=json.loads(r.shared_contributions)
        targetdata=[];vascular=False;regional=False;regional_phases=0;dev_phases=0;max_reg=0;max_dev=0;adequate=0;dev_top='';marker_evidence={}
        for target in r.substantial_shared_destinations.split(';'):
            ia,ib=groups[(a,target)],groups[(b,target)]
            delta=(identity.loc[ia].mean()-identity.loc[ib].mean())/scale
            max_reg=max(max_reg,float(delta[REGION].abs().max()))
            if float(delta[DEV].abs().max())>max_dev:
                max_dev=float(delta[DEV].abs().max());dev_top=delta[DEV].abs().idxmax()
            for gene in ['Nkx2-1','Lhx6','Lhx8','Sox6','Erbb4','Meis2','Isl1','Ebf1','Gsx2','Pecam1','Cldn5','Kdr']:
                gd=genes.loc[ia,gene].mean()-genes.loc[ib,gene].mean()
                dd=genes.loc[ia,gene].gt(0).mean()-genes.loc[ib,gene].gt(0).mean()
                if abs(gd)>=.5 and abs(dd)>=.1 and np.sign(gd)==np.sign(dd):
                    if gene not in marker_evidence or abs(gd)>abs(marker_evidence[gene]):marker_evidence[gene]=float(gd)
            sa,sb=scores.loc[ia].mean(),scores.loc[ib].mean()
            reg_here=((sa.LGE>.2 and sb.LGE<0 and sb.MGE>.3) or (sb.LGE>.2 and sa.LGE<0 and sa.MGE>.3))
            # Coherent vascular contrast in the actual merged subsets. Three
            # independent markers must be commonly detected on the high side.
            for hi,lo,hs,ls in [(ia,ib,sa,sb),(ib,ia,sb,sa)]:
                anchors=['Pecam1','Cdh5','Cldn5','Kdr']
                high=genes.loc[hi,anchors].gt(0).mean();low=genes.loc[lo,anchors].gt(0).mean()
                if hs.endothelial>=1 and ls.endothelial<.5 and ((high>=.5)&((high-low)>=.3)).sum()>=3:vascular=True
            rph=[];dph=[];ok=[]
            for phase in phase_names:
                ja=ia[phases.loc[ia].eq(phase)];jb=ib[phases.loc[ib].eq(phase)]
                eligible=min(len(ja),len(jb))>=20
                ed=(identity.loc[ja].mean()-identity.loc[jb].mean())/scale if eligible else delta*np.nan
                if len(ia)==int((labels==a).sum()) and len(ib)==int((labels==b).sum()):
                    # When entire original groups share a destination, these
                    # summaries must reproduce saved Step 07 phase effects.
                    for prog in REGION+DEV:
                        expected=original_effects.loc[(a,b,phase,prog),'standardized_effect']
                        assert np.isclose(ed[prog],expected,atol=1e-9,equal_nan=True), (a,b,phase,prog)
                        exact_phase_checks+=1
                if eligible:ok.append(phase)
                keep=ed.abs().ge(.5)&delta.abs().ge(.5)&np.sign(ed).eq(np.sign(delta))
                if keep.reindex(REGION).any():rph.append(phase)
                if keep.reindex(DEV).any():dph.append(phase)
                for prog in REGION+DEV:
                    floweffects.append(dict(boundary=r.boundary,comparison=target,phase=phase,program=prog,n_a=len(ja),n_b=len(jb),eligible=eligible,all_subset_effect_sd=delta[prog],matched_effect_sd=ed[prog],same_direction_retained=bool(keep[prog])))
            regional=regional or reg_here;regional_phases=max(regional_phases,len(rph) if reg_here else 0);dev_phases=max(dev_phases,len(dph));adequate=max(adequate,len(ok))
            targetdata.append(dict(comparison=target,n_a=len(ia),n_b=len(ib),region_opposition=bool(reg_here),region_phases=rph,development_phases=dph,adequate_phases=ok,means_a={k:float(sa[k]) for k in ['MGE','LGE','rg_stemness','ipc_neurogenic','neuroblast','endothelial']},means_b={k:float(sb[k]) for k in ['MGE','LGE','rg_stemness','ipc_neurogenic','neuroblast','endothelial']}))
        # Recommendations are for biological distinctions, not automatic
        # restoration of original membership or adoption of comparison labels.
        if vascular:
            recommendation='keep distinction';importance='Coherent vascular versus neural contrast survives in transferred cells; small-n phase limits remain.'
        elif regional and regional_phases>=2:
            recommendation='keep distinction' if state in ['strong','moderate'] else 'unresolved'
            importance='Opposing regional programs persist within matched phases; exact boundary seed '+state+'.'
        elif r.interpretation=='merge_resolves_weak_subdivision' and max(max_reg,max_dev)<.5 and adequate>=1 and r.original_de_gene_count<=50:
            recommendation='merge';importance='No strong canonical identity contrast in contributing cells; simplification plausible for discrete identities, not proof of identical states.'
        else:
            recommendation='unresolved'
            if regional:importance='Regional-program contrast present, but matched-phase or seed evidence is limited.'
            elif max_dev>=.5:importance='Developmental/state contrast; continuous maturation versus discrete identity is not resolved.'
            elif r.original_de_gene_count>50:importance='Canonical identity evidence is weak but original DE remains substantial; inspect unmeasured programs/QC before merging.'
            else:importance='No decisive independent identity axis; evidence insufficient to lock this boundary.'
        fate='; '.join('%s: %.1f%%/%.1f%%'%(t,100*shared[t]['a_fraction'],100*shared[t]['b_fraction']) for t in r.substantial_shared_destinations.split(';'))
        cycle=('strong' if r.cycle_composition_tv>=.4 else 'limited')+' composition difference (TV %.2f)'%r.cycle_composition_tv
        maturation=(dev_top if max_dev>=.5 else 'weak program contrast')+'; %.2f SD; %d/%d phases retain contrast'%(max_dev,dev_phases,adequate)
        lineage=('vascular–neural' if vascular else 'opposing regional programs' if regional else 'no decisive regional-origin contrast')+'; max %.2f SD'%max_reg
        lineage+='; '+', '.join(g+(' higher A' if val>0 else ' higher B') for g,val in sorted(marker_evidence.items(),key=lambda x:-abs(x[1]))[:4])
        sample='TV %.2f; %s %.0f%% / %s %.0f%%'%(r.sample_composition_tv,str(r.sample_a).replace('15662-',''),100*r.sample_fraction_a,str(r.sample_b).replace('15662-',''),100*r.sample_fraction_b)
        rows.append({'original sibling pair / lost boundary':a+' ↔ '+b,'parent coarse cluster':r.parent,'seed stability of boundary':state+'; cross-coassignment %.2f'%cross,'DE support under original settings':('%d genes; score %.1f; '%(r.original_de_gene_count,r.original_de_score))+('PASS' if r.original_de_separation else 'FAIL'),'cell-cycle association':cycle,'maturation association':maturation,'regional/lineage marker difference':lineage,'sample imbalance':sample,'fate under stricter criteria':r.extent+'; '+fate,'biological importance of lost distinction':importance,'recommendation':recommendation})
        detail.append(dict(boundary=r.boundary,seed_boundary=state,seed_cross_coassignment=cross,seed_dominant_a=seed.loc[a].idxmax(),seed_dominant_b=seed.loc[b].idxmax(),seed_retention_a=seed.loc[a].max(),seed_retention_b=seed.loc[b].max(),whole_cluster_seed_a=r.seed_a,whole_cluster_seed_b=r.seed_b,vascular_contrast=vascular,regional_opposition=regional,regional_matched_phases=regional_phases,developmental_matched_phases=dev_phases,max_regional_effect_sd=max_reg,max_developmental_effect_sd=max_dev,adequate_subset_phases=adequate,recommendation=recommendation,subset_evidence=json.dumps(targetdata)))
    table=pd.DataFrame(rows);details=pd.DataFrame(detail)
    order={'keep distinction':0,'unresolved':1,'merge':2};table=table.sort_values('recommendation',key=lambda s:s.map(order),kind='stable')
    table.to_csv(out/'decision_table.tsv',sep='\t',index=False);details.to_csv(out/'decision_evidence.tsv',sep='\t',index=False)
    pd.DataFrame(floweffects).to_csv(out/'contributing_subset_phase_effects.tsv',sep='\t',index=False)
    assert len(table)==84 and table.iloc[:,0].is_unique and set(table.recommendation)<=set(order)
    counts=table.recommendation.value_counts().to_dict()
    intro='<h1>38 → 17: boundary decision table</h1><p><b>Still unresolved. Do not lock 17.</b> Dissected E14.5 mouse MGE; 12,000 existing pilot cells. No new HiCAT fit or parameters.</p><p>All 84 substantially lost/eroded sibling boundaries. Percentages in the fate column are fractions of original A/B groups, not destination fractions. The 17 destination contribution tables and merge-flow tree are linked below.</p>'
    methods='<details><summary>How to read the evidence and recommendations</summary><p>Seed stability here describes separation in the existing repeat fit: different dominant repeat groups retaining ≥80% each and cross-coassignment ≤0.1 are strong; ≥50% each and cross-coassignment ≤0.2 are moderate. Otherwise low/unresolved. Low cross-coassignment alone is not proof because fragmentation lowers it. Whole-cluster seed results remain in decision_evidence.tsv.</p><p>Cell-cycle/sample TV measures composition difference (0 identical, 1 disjoint). Neither proves causality or batch. Program differences use the existing fixed pilot SD and existing scores. Destination-specific matched phases require ≥20 cells from each contributing subset and ≥0.5 SD same-direction differences. Small subsets remain unresolved. These are descriptive summaries of saved scores, not regression or new DE.</p><p>Coherent vascular contrast requires endothelial mean ≥1 versus &lt;0.5 and ≥3 of Pecam1/Cdh5/Cldn5/Kdr detected in ≥50% on the high side with ≥30 percentage-point detection difference. Opposing regional support requires LGE mean &gt;0.2 versus &lt;0 and MGE &gt;0.3 on the opposing side; this is a marker-program hypothesis, not proven anatomical origin. Both cycle and maturation can coexist.</p><p>Keep distinction means preserve that biological separation for review, not restore exact original memberships. Merge is a conservative candidate for discrete identities: weak canonical contrasts in both the original and contributing subsets, adequate phases, seed uncertainty, and ≤50 original DE genes. This gene-count limit is a review screen, not a HiCAT parameter. Substantial original DE is not dismissed just because selected canonical programs are similar. Maturation contrasts do not establish continuous trajectories or discrete types.</p></details>'
    links='<p><a href="../../outputs/merge_tree_38_to_17.pdf">38→17 merge-flow tree</a> · <a href="../../outputs/comparison_cluster_summary.tsv">All 17 destination contributions</a> · <a href="../../outputs/all_cell_contributions.tsv">Every nonzero flow</a> · <a href="decision_table.tsv">Decision TSV</a> · <a href="decision_evidence.tsv">Detailed evidence</a></p>'
    page='<!doctype html><meta charset="utf-8"><title>38 to 17 boundary decisions</title><style>body{font:15px system-ui;margin:24px;color:#18232d}table{border-collapse:collapse;width:100%;font-size:12px}th,td{padding:9px;border:1px solid #ccd3db;vertical-align:top;min-width:90px}th{position:sticky;top:0;background:#e7eef6}tr:nth-child(even){background:#f5f7fa}input{padding:10px;width:50%;margin:12px 0}summary{cursor:pointer}td:last-child{font-weight:bold}</style>'+intro+links+methods+'<p>Recommendations: '+html.escape(json.dumps(counts))+'</p><input id="filter" placeholder="Filter by cluster, recommendation, or evidence" oninput="document.querySelectorAll(\'tbody tr\').forEach(r=>r.hidden=!r.textContent.toLowerCase().includes(this.value.toLowerCase()))">'+table.to_html(index=False,escape=True)
    (out/'decision_table.html').write_text(page)
    summary=dict(status='IN_REVIEW',partition_assessment='still_unresolved',rows=len(table),recommendations=counts,new_hicat_fits=0,clustering_parameters_changed=False,source_cells=12000,source_flows=len(groups))
    assert exact_phase_checks>0
    (out/'summary.json').write_text(json.dumps(summary,indent=2));(out/'checks.json').write_text(json.dumps(dict(all_source_hashes=True,unique_12000_cells=True,all_84_boundaries=True,contribution_groups=len(groups),original_data_unchanged=True,phase_labels=phase_names,original_phase_effects_reproduced=exact_phase_checks),indent=2))
    pd.DataFrame([dict(file=p.name,sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()]).to_csv(out/'output_manifest.tsv',sep='\t',index=False)
    (package/'COMPUTATION_SUCCESS.txt').write_text('IN_REVIEW\n');print(json.dumps(summary),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--package',type=Path,required=True);run(p.parse_args().package)
