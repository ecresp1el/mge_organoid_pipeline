#!/usr/bin/env python3
"""Independent developmental axes, frozen DIV90 identities, Turbo-only runtime.

The phase-one directory is read-only input. See the phase-two evidence note for
gene-set rationale. Scores and state assignments are operational, not fates.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats, signal
from sklearn.mixture import GaussianMixture

from div90_hypergate_paths import PROJECT_ROOT, REPO_ROOT

OUT = PROJECT_ROOT / 'results/div90_hypergate_sst_pv_phase2'
OLD = PROJECT_ROOT / 'results/div90_hypergate_sst_pv'
STATES = ['PV-biased', 'PV/SST hybrid', 'SST-biased', 'unresolved/immature']
PV = ['MEF2C', 'KCNC1', 'KCNC2', 'TAC1']
SST = ['SST', 'SATB1', 'NR2F2', 'CDK14', 'CACNG3']
EXTRA = ['MAF','MAFB','GRIK1','TMEM91','ST18','ELMO1','IGFBP4','SOX6','CRABP1','PVALB','GAD1','GAD2']
STRESS = ['FOS','JUN','JUNB','ATF3','DDIT3','HSPA1A','HSPA1B','HSPB1']
RATIONALE = {
 'MEF2C': 'PV-associated developmental regulator; shared and context dependent, not fate specific.',
 'KCNC1': 'Kv3 fast-spiking maturation-associated effector; functional program rather than fate.',
 'KCNC2': 'Kv3 fast-spiking maturation-associated effector; functional program rather than fate.',
 'TAC1': 'PV-associated cortical transcript in developmental reference; also expressed outside PV.',
 'SST': 'SST-associated peptide transcript; equal 1/5 module weight, neither necessary nor sufficient for commitment.',
 'SATB1': 'Activity-dependent SST maturation regulator; also participates in shared interneuron maturation.',
 'NR2F2': 'Early MGE SST-associated regulator; broader CGE/subpallial roles and sparse detection require sensitivity.',
 'CDK14': 'Early SST-associated transcript in Mi developmental MGE reference.',
 'CACNG3': 'SST-associated transcript in Allaway developing cortical reference; excluded from gate features.',
 'MAF': 'Shared, stage-dependent PV/SST regulator; expanded PV sensitivity only.',
 'MAFB': 'Shared, stage-dependent PV/SST regulator; expanded PV sensitivity only.',
 'GRIK1': 'Broad interneuron expression and sparse DIV90 detection; phase-one-compatible sensitivity only.',
 'TMEM91': 'Alternative SST-associated reference transcript; alternative module sensitivity only.',
 'ST18': 'Noncortical MGE/globus pallidus specificity complicates cortical PV interpretation; excluded.',
 'ELMO1': 'Immature PV evidence strongest in mouse hippocampus; exploratory expanded sensitivity only.',
 'IGFBP4': 'Immature PV evidence strongest in mouse hippocampus; exploratory expanded sensitivity only.',
 'SOX6': 'Shared MGE interneuron regulator; excluded from discriminatory modules.',
 'CRABP1': 'Existing precursor annotation aid; not used to define independent scores.',
 'PVALB': 'Only three detected cells; not required for PV-associated programming.',
 'GAD1': 'Generic GABAergic identity; excluded from primary discrimination.',
 'GAD2': 'Generic GABAergic identity; excluded from primary discrimination.',
}


def clean_json(x):
    if isinstance(x, dict): return {str(k): clean_json(v) for k,v in x.items()}
    if isinstance(x, (list,tuple)): return [clean_json(v) for v in x]
    if isinstance(x, np.ndarray): return clean_json(x.tolist())
    if isinstance(x, np.generic): return clean_json(x.item())
    if isinstance(x, float) and not np.isfinite(x): return None
    return x


def write_json(path, value):
    path.write_text(json.dumps(clean_json(value), indent=2, allow_nan=False)+'\n')


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''): h.update(block)
    return h.hexdigest()


def initialize():
    turbo = Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder').resolve()
    if OUT.resolve() != turbo / 'results/div90_hypergate_sst_pv_phase2':
        raise ValueError('This analysis is authorized for the specified Turbo phase-two runtime only.')
    for p in ['tables','provenance','figures/png','figures/pdf','figures/svg','interactive','cache','logs']:
        (OUT/p).mkdir(parents=True, exist_ok=True)
    manifest = OUT/'provenance/phase1_sha256_before.json'
    if not manifest.exists():
        write_json(manifest, {str(p.relative_to(OLD)):sha(p) for p in OLD.rglob('*') if p.is_file()})
    for name in ['run_div90_hypergate_phase2.py','div90_hypergate_phase2_gates.py','div90_hypergate_phase2_visuals.py']:
        src = REPO_ROOT/'python_notebooks/scripts'/name
        if src.exists(): shutil.copy2(src, OUT/'provenance'/name)
    src = REPO_ROOT/'python_notebooks/DIV90_HYPERGATE_PHASE2_EVIDENCE.md'
    if src.exists(): shutil.copy2(src, OUT/'provenance'/src.name)


def state_labels(pv, sst, pcut, scut):
    ph, sh = np.asarray(pv)>=pcut, np.asarray(sst)>=scut
    return np.select([ph & ~sh, ph & sh, ~ph & sh], STATES[:3], default=STATES[3])


def bounded_z(frame, genes):
    x = frame[genes].astype(float)
    mean, sd = x.mean(), x.std(ddof=0)
    if (sd<=0).any() or x.isna().any().any(): raise ValueError('Nonvarying/missing module gene')
    z = ((x-mean)/sd).clip(-3,3)
    return z, dict(genes=genes, means=mean.to_dict(), std_population=sd.to_dict(), clipping=[-3,3])


def examine_density(cells):
    fits, modes = [], []
    for label, cols in [('pv',['pv_score']),('sst',['sst_score']),('joint',['sst_score','pv_score'])]:
        a = cells[cols].to_numpy()
        for k in [1,2,3,4]:
            model=GaussianMixture(k, random_state=90, n_init=5, reg_covar=1e-5).fit(a)
            fits.append(dict(axis=label,k=k,bic=model.bic(a),aic=model.aic(a),weights=json.dumps(model.weights_.tolist()),means=json.dumps(model.means_.tolist())))
        if len(cols)==1:
            v=a[:,0];grid=np.linspace(v.min(),v.max(),600)
            for mult in [.7,1.,1.4]:
                kde=stats.gaussian_kde(v,bw_method=lambda x: x.scotts_factor()*mult);d=kde(grid)
                peaks,_=signal.find_peaks(d,prominence=.08*d.max())
                for j in peaks: modes.append(dict(axis=label,bandwidth_multiplier=mult,mode=grid[j],density=d[j]))
    pd.DataFrame(fits).to_csv(OUT/'tables/mixture_diagnostics.tsv',sep='\t',index=False)
    pd.DataFrame(modes).to_csv(OUT/'tables/density_modes.tsv',sep='\t',index=False)
    return fits,modes


def prepare():
    if (OUT/'cells.tsv.gz').exists():
        print('Reusing phase-two cells checkpoint',flush=True)
        return pd.read_csv(OUT/'cells.tsv.gz',sep='\t')
    cells=pd.read_csv(OLD/'hypergate_input.tsv.gz',sep='\t')
    baseline=pd.read_csv(OLD/'target_labels.tsv.gz',sep='\t')
    assert len(cells)==4768 and cells.cell_id.is_unique and cells.cell_id.equals(baseline.cell_id)
    hpath=PROJECT_ROOT/'results/python_anndata/varela_div90.h5ad'
    a=ad.read_h5ad(hpath,backed='r')
    idx=a.obs_names.get_indexer(cells.cell_id)
    assert (idx>=0).all()
    obs=a.obs.iloc[idx].reset_index(drop=True)
    assert np.array_equal(obs['orig.ident'].astype(str),cells['sample'])
    genes=list(dict.fromkeys(PV+SST+EXTRA+STRESS))
    present=[g for g in genes if g in a.var_names]
    x=a.X[idx,:][:,a.var_names.get_indexer(present)].toarray()
    for j,g in enumerate(present): cells[g]=x[:,j]
    for g in set(genes)-set(present): cells[g]=np.nan
    for target,source in [('total_counts','nCount_RNA'),('n_genes','nFeature_RNA'),('S_score','S.Score'),('G2M_score','G2M.Score')]:
        cells[target]=pd.to_numeric(obs[source],errors='raise')
    cells['mito_fraction']=pd.to_numeric(obs['percent.mt'],errors='raise')/100
    cells['cell_cycle_phase']=obs['Phase'].astype(str)
    cells['cell_cycle_score']=cells[['S_score','G2M_score']].max(axis=1)
    cells['doublet_score']=np.nan
    cells['doublet_status']='not available in original H5AD'
    mapping=pd.read_csv(REPO_ROOT/'metadata/div30_div90_sample_id_to_biolabel_map.tsv',sep='\t').query("DIV == 'DIV90'").set_index('run_sample_id')
    cells['biological_label']=cells['sample'].map(mapping.biological_label)
    assert cells.biological_label.notna().all()
    cells['cell_line']=cells.biological_label.str.replace(r'_rep\d+$','',regex=True)
    # User-provided experimental metadata, checked against the six recorded IDs.
    # The MW substring in every technical run ID is not the culture operator.
    expected={f'10496-MW-{i}' for i in range(1,7)}
    assert set(cells['sample'])==expected and set(mapping.index)==expected
    sample_number=cells['sample'].str.extract(r'-(\d+)$')[0].astype(int)
    cells['culture_operator']=np.where(sample_number%2==1,'CV','MW')
    cells['relative_glucose']=np.where(sample_number%2==1,2.,1.)
    cells['condition']=np.where(sample_number%2==1,'CV / higher glucose (~2x)','MW / lower glucose (~1x)')
    cells['annotation']=cells.loupe_label
    mapfile=OUT/'sample_conditions.tsv'
    if mapfile.exists():
        cm=pd.read_csv(mapfile,sep='\t').set_index('sample')
        cells['condition']=cells['sample'].map(cm['condition']).fillna('not recorded')
    cells[['sample','biological_label','cell_line','condition','culture_operator','relative_glucose']].drop_duplicates().to_csv(OUT/'tables/sample_metadata.tsv',sep='\t',index=False)
    write_json(OUT/'provenance/condition_metadata.json',dict(source='User experimental metadata update in this conversation, 2026-09-09',div90_assignment='CV odd samples1,3,5; MW even samples2,4,6',recorded_id_verification='All six IDs and paired H9/79B/2E line assignments confirmed against repository mapping and H5AD orig.ident; no operator/condition field in H5AD to independently verify wet-lab glucose.',existing_biological_labels_preserved=True,conflict='No explicit operator/glucose field contradicts user note. MW within all run identifiers is a technical sample naming convention; it is not used as culture assignment.',interpretation='CV versus MW culture condition, including an approximately 2-fold glucose difference; glucose confounded with operator/other culture differences. One sample per line x condition; cells are not biological replicates.'))
    modulepars={}
    for name,gs in [('pv',PV),('sst',SST)]:
        z,pars=bounded_z(cells,gs);cells[name+'_score']=z.mean(axis=1)
        cells[name+'_n_detected']=(cells[gs]>0).sum(axis=1)
        for g in gs: cells['z_'+g]=z[g]
        modulepars[name]=pars
    stressgenes=[g for g in STRESS if g in present and cells[g].std(ddof=0)>0]
    z,_=bounded_z(cells,stressgenes);cells['stress_score']=z.mean(axis=1)
    # Independent axes: neither score is subtracted from the other.
    pcut=float(cells.pv_score.median());scut=float(cells.sst_score.median())
    cells['state']=state_labels(cells.pv_score,cells.sst_score,pcut,scut)
    cells['strong_sst_biased']=(cells.state==STATES[2]) & (cells.sst_score>=cells.sst_score.quantile(.65)) & (cells.pv_score<=cells.pv_score.quantile(.35))
    cells['sst_n_detected_without_SST']=(cells[[g for g in SST if g!='SST']]>0).sum(axis=1)
    cells['multi_gene_both']=(cells.pv_n_detected>=2)&(cells.sst_n_detected>=2)
    cells['multi_gene_both_without_SST']=(cells.pv_n_detected>=2)&(cells.sst_n_detected_without_SST>=2)
    oldrule=json.loads((OLD/'best_gate_rules.json').read_text())['best_shared_pair']['rules']
    cap=np.ones(len(cells),bool)
    for r in oldrule: cap &= cells[r['gene']].to_numpy()<=r['threshold'] if r['op']=='<=' else cells[r['gene']].to_numpy()>=r['threshold']
    cells['phase1_gate_retained']=cap;cells['fat3_ptprm_retained']=cap
    assert cap.sum()==3168 and (cells.SST>0).sum()==2919 and (cells.PVALB>0).sum()==3
    assert np.allclose(cells[['loupe_x','loupe_y']],baseline[['loupe_x','loupe_y']],rtol=0,atol=0)
    rows=[]
    for g in genes:
        rows.append(dict(gene=g,program='pv' if g in PV else 'sst' if g in SST else 'stress' if g in STRESS else 'sensitivity_or_audit',included=g in PV+SST,detected=int(cells[g].gt(0).sum()),detection_fraction=cells[g].gt(0).mean(),weight=1/len(PV) if g in PV else 1/len(SST) if g in SST else 0,rationale=RATIONALE.get(g,'Stress/activation sensitivity only; not lineage-specific.')))
    pd.DataFrame(rows).to_csv(OUT/'module_genes.tsv',sep='\t',index=False)
    fits,modes=examine_density(cells)
    write_json(OUT/'state_parameters.json',dict(pv_threshold=pcut,sst_threshold=scut,thresholds={'pv':pcut,'sst':scut},modules=modulepars,scoring='Equal-weight mean of gene-wise population z-scores of existing log1p(CP10K), bounded at [-3,3] per gene.',state_rule='Independent pooled median thresholds; high >= median. Every cell assigned. Operational partitions, not natural clusters or fates.',threshold_rationale='Mixture fits and KDE diagnostics are saved. Non-Gaussian sparse-expression mixtures do not establish biological natural breaks. Transparent median partitions are primary; 35/50/65% full-grid sensitivity retains every cell.',fate_inference=False,stress_genes=stressgenes,seed=90))
    cells.to_csv(OUT/'cells.tsv.gz',sep='\t',index=False)
    write_json(OUT/'provenance/input_validation.json',dict(n_cells=len(cells),n_sst_detected=int((cells.SST>0).sum()),n_pvalb_detected=int((cells.PVALB>0).sum()),phase1_gate_capture=int(cap.sum()),cell_ids_order_preserved=True,loupe_coordinates_unchanged=True,source_h5ad=str(hpath),source_h5ad_size=hpath.stat().st_size,source_h5ad_mtime_ns=hpath.stat().st_mtime_ns,condition_source='sample_conditions.tsv when supplied; otherwise unknown',doublet_original_available=False))
    a.file.close()
    print(cells.state.value_counts().to_string(),flush=True)
    return cells


def sensitivity_and_qc(cells):
    doubletfile=OUT/'tables/doublet_scores_all_div90.tsv.gz'
    if doubletfile.exists():
        doublets=pd.read_csv(doubletfile,sep='\t').set_index('cell_id')
        assert cells.cell_id.isin(doublets.index).all()
        cells['doublet_score']=cells.cell_id.map(doublets.doublet_score)
        cells['predicted_doublet']=cells.cell_id.map(doublets.predicted_doublet).astype(bool)
        cells['doublet_status']=np.where(cells.predicted_doublet,'exploratory predicted doublet','below exploratory doublet threshold')
        cells['doublet_within_sample_percentile']=cells.groupby('sample').doublet_score.rank(pct=True)
    params=json.loads((OUT/'state_parameters.json').read_text())
    pcut,scut=params['pv_threshold'],params['sst_threshold']
    scores=cells[['cell_id']].copy(); rows=[]
    def record(name,pv,sst,pc=None,sc=None):
        pc=float(np.median(pv)) if pc is None else pc;sc=float(np.median(sst)) if sc is None else sc
        # Quantiles at the all-undetected score floor cannot split tied cells.
        # Assign that complete floor to low, rather than calling every cell high
        # or breaking ties by barcode. Primary medians are above both floors.
        adjusted=[]
        if pc<=float(np.min(pv))+1e-12:
            pc=float(np.min(pv))+1e-12;adjusted.append('pv')
        if sc<=float(np.min(sst))+1e-12:
            sc=float(np.min(sst))+1e-12;adjusted.append('sst')
        labels=state_labels(pv,sst,pc,sc);hybrid=labels==STATES[1];base=cells.state.eq(STATES[1]).to_numpy()
        scores[name+'_pv_score']=pv;scores[name+'_sst_score']=sst;scores[name+'_state']=labels
        row=dict(variant=name,pv_threshold=pc,sst_threshold=sc,floor_tie_assigned_low=';'.join(adjusted),pv_high_fraction=float((np.asarray(pv)>=pc).mean()),sst_high_fraction=float((np.asarray(sst)>=sc).mean()),hybrid_n=hybrid.sum(),hybrid_fraction=hybrid.mean(),baseline_hybrid_retained=(hybrid&base).sum()/base.sum(),hybrid_jaccard=(hybrid&base).sum()/(hybrid|base).sum(),all_state_agreement=(labels==cells.state).mean())
        for state in STATES:row[state]=int((labels==state).sum())
        rows.append(row)
    for qp in [.35,.5,.65]:
        for qs in [.35,.5,.65]: record(f'q{int(qp*100)}_q{int(qs*100)}',cells.pv_score,cells.sst_score,float(cells.pv_score.quantile(qp)),float(cells.sst_score.quantile(qs)))
    for program,genes in [('pv',PV),('sst',SST)]:
        for gene in genes:
            pv=cells.pv_score if program!='pv' else cells[['z_'+g for g in genes if g!=gene]].mean(axis=1)
            sst=cells.sst_score if program!='sst' else cells[['z_'+g for g in genes if g!=gene]].mean(axis=1)
            record('without_'+gene,pv,sst)
    variants={'expanded_shared_PV':(PV+['MAF','MAFB'],SST),'minimal_PV':(['MEF2C','KCNC1','KCNC2'],SST),'hippocampal_PV_exploratory':(PV+['ELMO1','IGFBP4'],SST),'alternative_SST_TMEM91':(PV,['SST','SATB1','TMEM91','CDK14','CACNG3']),'phase1_genes_independent_axes':(['MEF2C','MAF','MAFB','KCNC1','KCNC2','GAD1','GAD2'],['SST','NR2F2','GRIK1'])}
    for name,(pg,sg) in variants.items():record(name,bounded_z(cells,pg)[0].mean(axis=1),bounded_z(cells,sg)[0].mean(axis=1))
    rawpv=(cells[PV]-cells[PV].mean())/cells[PV].std(ddof=0);rawsst=(cells[SST]-cells[SST].mean())/cells[SST].std(ddof=0)
    record('unbounded_z',rawpv.mean(axis=1),rawsst.mean(axis=1))
    # Diagnostic only: center depth covariates within sample, preserving each sample's score mean.
    depth=np.log1p(cells[['total_counts','n_genes']]).copy()
    depth-=depth.groupby(cells['sample']).transform('mean')
    for program in ['pv','sst']:
        v=cells[program+'_score'].to_numpy();center=v-cells.groupby('sample')[program+'_score'].transform('mean')
        coef=np.linalg.lstsq(depth,center,rcond=None)[0]
        cells[program+'_depth_diagnostic']=v-depth.to_numpy()@coef
    record('within_sample_depth_diagnostic',cells.pv_depth_diagnostic,cells.sst_depth_diagnostic)
    pd.DataFrame(rows).to_csv(OUT/'tables/module_threshold_sensitivity.tsv',sep='\t',index=False)
    scores.to_csv(OUT/'tables/cell_sensitivity.tsv.gz',sep='\t',index=False)
    qc=['total_counts','n_genes','mito_fraction','stress_score','cell_cycle_score','doublet_score']
    cells['high_complexity_tail']=cells.total_counts>=cells.groupby('sample').total_counts.transform(lambda x:x.quantile(.95))
    cells['high_mito_tail']=cells.mito_fraction>cells.groupby('sample').mito_fraction.transform(lambda x:x.quantile(.95))
    cells['high_stress_activation_tail']=cells.stress_score>cells.groupby('sample').stress_score.transform(lambda x:x.quantile(.95))
    summaries=[]
    for state,d in cells.groupby('state',sort=False):
        row=dict(state=state,n=len(d),fraction=len(d)/len(cells),pv_median=d.pv_score.median(),sst_median=d.sst_score.median(),multi_gene_both_fraction=d.multi_gene_both.mean(),multi_gene_both_without_SST_fraction=d.multi_gene_both_without_SST.mean(),pv_genes_detected_median=d.pv_n_detected.median(),sst_genes_detected_median=d.sst_n_detected.median())
        for c in qc:
            for stat,val in [('median',d[c].median()),('q25',d[c].quantile(.25)),('q75',d[c].quantile(.75))]:row[c+'_'+stat]=val
        for c in ['high_complexity_tail','high_mito_tail','high_stress_activation_tail','predicted_doublet']:
            if c in d:row[c+'_fraction']=float(d[c].mean())
        summaries.append(row)
    pd.DataFrame(summaries).to_csv(OUT/'tables/state_qc_summary.tsv',sep='\t',index=False)
    corr=[]
    for x in ['pv_score','sst_score','pv_depth_diagnostic','sst_depth_diagnostic']:
        for y in ['pv_score','sst_score',*qc]:
            d=cells[[x,y]].dropna()
            rho=stats.spearmanr(d.iloc[:,0],d.iloc[:,1]).statistic if len(d)>2 and d.iloc[:,1].nunique()>1 else np.nan
            corr.append(dict(x=x,y=y,spearman_rho=rho,n=len(d)))
    pd.DataFrame(corr).to_csv(OUT/'tables/score_qc_correlations.tsv',sep='\t',index=False)
    gene=[]
    for state,d in cells.groupby('state'):
        for g in PV+SST+EXTRA:
            gene.append(dict(state=state,gene=g,program='pv' if g in PV else 'sst' if g in SST else 'audit',mean_expression=d[g].mean(),fraction_expressing=d[g].gt(0).mean(),mean_z=cells.loc[d.index,'z_'+g].mean() if 'z_'+g in cells else np.nan))
    pd.DataFrame(gene).to_csv(OUT/'tables/module_gene_by_state.tsv',sep='\t',index=False)
    for group in ['sample','cell_line','condition','loupe_label']:
        rows=[]
        for label,d in cells.groupby(group):
            row={group:label,'n':len(d),'median_pv_score':d.pv_score.median(),'median_sst_score':d.sst_score.median()}
            for s in STATES:row[s+'_n']=int(d.state.eq(s).sum());row[s+'_fraction']=d.state.eq(s).mean()
            rows.append(row)
        pd.DataFrame(rows).to_csv(OUT/f'tables/composition_by_{group}.tsv',sep='\t',index=False)
    # Permute individual genes within sample and depth quintile: null preserves
    # sample composition and approximate complexity while breaking coordination.
    rng=np.random.default_rng(90);bins=cells.groupby('sample').total_counts.transform(lambda v:pd.qcut(v.rank(method='first'),5,labels=False))
    strata=list(cells.groupby([cells['sample'],bins]).indices.values())
    z=cells[['z_'+g for g in PV+SST]].to_numpy();det=cells[PV+SST].gt(0).to_numpy()
    null=[]
    for b in range(200):
        znull=z.copy();dnull=det.copy()
        for ix in strata:
            for j in range(z.shape[1]):
                perm=rng.permutation(ix);znull[ix,j]=z[perm,j];dnull[ix,j]=det[perm,j]
        ph=znull[:,:len(PV)].mean(1)>=pcut;sh=znull[:,len(PV):].mean(1)>=scut
        multi=(dnull[:,:len(PV)].sum(1)>=2)&(dnull[:,len(PV):].sum(1)>=2)
        null.append(dict(iteration=b,dual_high_n=int((ph&sh).sum()),multi_gene_dual_high_n=int((ph&sh&multi).sum())))
    pd.DataFrame(null).to_csv(OUT/'tables/hybrid_depth_stratified_null.tsv',sep='\t',index=False)
    hybrid=cells.state.eq(STATES[1]);observed=int(hybrid.sum());obs_multi=int((hybrid&cells.multi_gene_both).sum())
    pv_c=cells[PV].corr(method='spearman');sst_c=cells[SST].corr(method='spearman')
    cells[PV+SST].corr(method='spearman').to_csv(OUT/'tables/module_gene_spearman.tsv',sep='\t')
    validation=dict(hybrid_n=observed,hybrid_fraction=hybrid.mean(),hybrid_multiple_genes_both_n=obs_multi,hybrid_multiple_genes_both_fraction=obs_multi/observed,hybrid_multiple_genes_both_without_SST_n=int((hybrid&cells.multi_gene_both_without_SST).sum()),permutation='200 gene-wise permutations within sample x UMI quintile; fixed original scores calibration and cutoffs; descriptive conditional null',null_dual_high_mean=np.mean([v['dual_high_n'] for v in null]),null_dual_high_95_interval=np.quantile([v['dual_high_n'] for v in null],[.025,.975]),null_multigene_mean=np.mean([v['multi_gene_dual_high_n'] for v in null]),null_multigene_95_interval=np.quantile([v['multi_gene_dual_high_n'] for v in null],[.025,.975]),rho_pv_sst=stats.spearmanr(cells.pv_score,cells.sst_score).statistic,rho_depth_diagnostic=stats.spearmanr(cells.pv_depth_diagnostic,cells.sst_depth_diagnostic).statistic,median_within_pv_gene_rho=np.median(pv_c.to_numpy()[np.triu_indices(len(PV),1)]),median_within_sst_gene_rho=np.median(sst_c.to_numpy()[np.triu_indices(len(SST),1)]),doublet_assessment='Original doublet scores absent; exploratory simulation score may be appended separately. High complexity alone does not establish doublets.')
    write_json(OUT/'hybrid_validation.json',validation)
    if 'predicted_doublet' in cells:
        validation.update(exploratory_predicted_doublets_total=int(cells.predicted_doublet.sum()),exploratory_predicted_doublets_hybrid=int((cells.predicted_doublet&hybrid).sum()),hybrid_doublet_score_median=cells.loc[hybrid,'doublet_score'].median(),doublet_assessment='Per-sample Scrublet on all 22338 available DIV90 cells; expected rate assumed 5%, automated thresholds, exploratory only. No primary exclusion. Homotypic doublets and prefiltering remain limitations.')
        write_json(OUT/'hybrid_validation.json',validation)
    qrows=[]
    for label,keep in [('all',np.ones(len(cells),bool)),('without_high_complexity_tail',~cells.high_complexity_tail),('without_high_mito_tail',~cells.high_mito_tail),('without_high_stress_activation_tail',~cells.high_stress_activation_tail),('without_exploratory_doublets',~cells.predicted_doublet if 'predicted_doublet' in cells else np.ones(len(cells),bool)),('without_top10pct_doublet_score',cells.doublet_within_sample_percentile<=.9 if 'doublet_within_sample_percentile' in cells else np.ones(len(cells),bool))]:
        d=cells.loc[keep];r={'sensitivity':label,'n':len(d),'rho_pv_sst':stats.spearmanr(d.pv_score,d.sst_score).statistic}
        for s in STATES:r[s+'_n']=int(d.state.eq(s).sum());r[s+'_fraction']=float(d.state.eq(s).mean())
        qrows.append(r)
    pd.DataFrame(qrows).to_csv(OUT/'tables/qc_exclusion_sensitivity.tsv',sep='\t',index=False)
    cells.to_csv(OUT/'cells.tsv.gz',sep='\t',index=False)
    print(json.dumps(clean_json(validation),indent=2),flush=True)
    return cells


def verify_phase1():
    before=json.loads((OUT/'provenance/phase1_sha256_before.json').read_text())
    after={str(p.relative_to(OLD)):sha(p) for p in OLD.rglob('*') if p.is_file()}
    assert before==after, 'Phase-one files changed during phase two'
    write_json(OUT/'provenance/phase1_preservation.json',dict(all_phase1_files_unchanged=True,n_files=len(before),algorithm='SHA256'))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--stage',choices=['prepare','qc','gates','figures','verify','all'],default='all');args=ap.parse_args()
    initialize()
    if args.stage=='verify':verify_phase1();return
    cells=prepare()
    if args.stage in ['qc','all']:cells=sensitivity_and_qc(cells)
    if args.stage in ['gates','all']:
        from div90_hypergate_phase2_gates import run
        run(cells,OUT)
    if args.stage in ['figures','all']:
        from div90_hypergate_phase2_visuals import run
        run(cells,OUT,json.loads((OUT/'gate_summary.json').read_text()))
    verify_phase1()


if __name__=='__main__':main()
