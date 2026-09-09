#!/usr/bin/env python3
"""Read-only local cross-study extension of the frozen DIV30 identity analysis.

No gate search, cross-study raw gate, pooled integration, or target redefinition.
Run stages prepare, map, compare. All generated files belong to new external/.
"""
from __future__ import annotations
import os
for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[name]='2'
from pathlib import Path
import argparse, json, hashlib, gzip
import h5py
import numpy as np
import pandas as pd
from scipy import sparse,stats
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import roc_auc_score

PROJECT=Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder')
OLD=PROJECT/'results/div30_pv_antecedent_v1'
OUT=PROJECT/'results/div30_antecedent_identity_v1/external'
MAP=OUT.parent/'mapping/cache'
CAN=PROJECT/'paper2_mge_organoid_atlas/inputs/canonical'
STUDIES=['walsh','bershteyn_2023','bershteyn_2025','siebert_2026']
if (OUT/'tables/samarasinghe_2021_cells_prepared.tsv.gz').exists():STUDIES.append('samarasinghe_2021')
PROGRAMS={
 'GABA_handling':['SLC6A1','SLC6A11','ABAT','ALDH5A1','GABRA1','GABRB2','GABRG2'],
 'synaptic_maturation':['SYP','SYN1','SYN2','STX1A','VAMP2','CPLX1','DLG4','RIMS1'],
 'ion_channel_maturation':['SCN1A','SCN2A','SCN8A','KCNB1','KCND2','KCNQ2','KCNQ3','HCN1','CACNA1A'],
 'cell_adhesion':['NCAM1','NRXN1','NRXN2','NLGN2','CNTNAP2','CDH2','PCDH9','PCDH19'],
 'ECM':['BCAN','NCAN','TNR','HAPLN1','VCAN','AGRN','LAMA2','COL4A1'],
 'growth_factor_signaling':['FGFR1','FGFR2','IGF1R','EGFR','IGFBP2','IGFBP3','IGFBP5','SPRY2','DUSP6'],
 'metabolic':['NDUFA4','NDUFB8','SDHB','UQCRC1','COX5A','ATP5F1A','PDHA1','IDH3A','CS','LDHA','PKM','HK1'],
 'activity_related':['FOS','JUN','EGR1','ARC','NPAS4','NR4A1','DUSP1','BTG2'],
 'axon_dendrite':['GAP43','DPYSL2','DPYSL3','NEFL','NEFM','MAP1B','STMN1','RND3'],
 'trophic_signaling':['BDNF','NTRK2','NTRK3','NGFR','GDNF','RET','CNTF','CNTFR','NRG1'],
}
PANELS=json.loads((OLD/'provenance/scoring_config.json').read_text())['panels']
PV=json.loads((OLD/'provenance/scoring_config.json').read_text())['pv']
EXTRA=sorted(set(sum(PROGRAMS.values(),[])+sum(PANELS.values(),[])+PV+['ERBB4','SST','PVALB']))

def dump(name,x):
    (OUT/name).write_text(json.dumps(x,indent=2,default=lambda z:z.item() if isinstance(z,np.generic) else str(z))+'\n')
def tab(name,x):
    comp={'method':'gzip','compresslevel':1} if str(name).endswith('.gz') else None
    x.to_csv(OUT/'tables'/name,sep='\t',index=False,compression=comp)
def read(x):
    if isinstance(x,h5py.Group):
        if 'values' in x:
            v=read(x['values']).astype(object);v[x['mask'][:]]='';return v
        c=read(x['categories']);return np.array([c[i] if i>=0 else '' for i in x['codes'][:]])
    return x.asstr()[:] if h5py.check_string_dtype(x.dtype) is not None else x[:]
def metadata(p):
    with h5py.File(p) as f:
        obs=pd.DataFrame({k:read(v) for k,v in f['obs'].items()})
        v=f['var'];genes=list(read(v['source_feature_id'] if 'source_feature_id' in v else v[v.attrs['_index']]))
    return obs,genes
def stream(p,layer='X',block=1500):
    with h5py.File(p) as f:
        x=f[layer];ip=x['indptr'][:];n,g=x.attrs['shape']
        for s in range(0,n,block):
            t=min(s+block,n);a,b=ip[s],ip[t]
            z=sparse.csr_matrix((x['data'][a:b].astype('float32'),x['indices'][a:b],ip[s:t+1]-a),shape=(t-s,g))
            yield s,t,z
def getgenes():
    with h5py.File(OLD/'cache/mapping_model.h5') as f:features=list(f['features'].asstr()[:])
    return sorted(set(features+EXTRA))
def statscores(x,genes,cells):
    pars=json.loads((OLD/'provenance/context_score_parameters.json').read_text())
    for panel,q in pars.items():
        ix=pd.Index(genes).get_indexer(q['genes'])
        mu=np.array([q['means'][g] for g in q['genes']]);sd=np.array([q['std_population'][g] for g in q['genes']])
        v=np.array(x[:,ix],copy=True)
        # Frozen TUBB3 is absent in MGEO's catalog: old analysis inserted 0.
        # Preserve the same denominator and constant contribution across studies.
        constant=(mu==0)&(sd==1);v[:,constant]=0
        cells[panel+'_score']=np.mean(np.clip((np.nan_to_num(v,nan=0)-mu)/sd,-3,3),axis=1)
    cells['postmitotic_continuous']=cells.maturation_score-cells.progenitor_score
    cells['postmitotic_proxy']=(cells.postmitotic_continuous>0)&(cells.cycle_score<.5)
    return pd.concat([cells,pd.DataFrame({g:x[:,genes.index(g)] for g in EXTRA},index=cells.index)],axis=1)

def prepare():
    for p in ['cache','tables','figures','provenance']:(OUT/p).mkdir(parents=True,exist_ok=True)
    genes=getgenes();dump('cache/genes.json',genes)
    dump('provenance/predeclared_programs.json',{'panels':PROGRAMS,'interpretation':'Operational unweighted expression panels, not formal pathway enrichment. Gene lists fixed before study contrasts. Primary target-defining genes, GRIA2/OPCML/NOTCH1, and context-panel genes excluded from program panels except STMN1 generic neurites, which is not in matching scores. SLC6A1 independent of direct PV score but present in original mapping; no-SLC6A1 sensitivity required.'})
    manifest=[];coverage=[]
    for study in ['mgeo_div30','walsh','bershteyn_2023','bershteyn_2025','siebert_2026']:
        path=PROJECT/'results/python_anndata/varela_div30.h5ad' if study=='mgeo_div30' else CAN/study/(study+'_minimal.h5ad')
        obs,sourcegenes=metadata(path);ix=pd.Index(sourcegenes).get_indexer(genes)
        manifest.append(dict(study=study,path=str(path),n_cells=len(obs),n_genes=len(sourcegenes),size_bytes=path.stat().st_size,mtime_ns=path.stat().st_mtime_ns))
        for g,j in zip(genes,ix):coverage.append(dict(study=study,gene=g,present=j>=0,source_feature_index=int(j)))
        exprpath=OUT/'cache'/f'{study}_expression.npy'
        if exprpath.exists() and (OUT/'tables'/f'{study}_cells_prepared.tsv.gz').exists():continue
        x=np.lib.format.open_memmap(exprpath,mode='w+',dtype='float32',shape=(len(obs),len(genes)));x[:]=np.nan
        counts=np.empty(len(obs));ngenes=np.empty(len(obs)); detected=np.zeros(len(genes),int)
        for s,t,z in stream(path):
            totals=np.asarray(z.sum(axis=1)).ravel();counts[s:t]=totals;ngenes[s:t]=z.getnnz(axis=1)
            if study!='mgeo_div30':
                z=sparse.diags(10000/np.maximum(totals,1)).dot(z).tocsr();z.data=np.log1p(z.data)
            x[s:t,np.where(ix>=0)[0]]=z[:,ix[ix>=0]].toarray()
            if s%15000==0:print('prepared',study,s,'/',len(obs),flush=True)
        x.flush()
        if study=='mgeo_div30':
            cells=pd.read_csv(OLD/'cells.tsv.gz',sep='\t');assert np.array_equal(cells.cell_id,obs.cell_id)
            cells['target']=cells.early_pv_score>=cells.early_pv_score.quantile(.8)
            cells=cells[['cell_id','sample','cell_line','condition','target','early_pv_score','postmitotic','cycling','direct_pv_score','total_counts','n_genes']].copy()
        else:
            cells=obs.rename(columns={'canonical_cell_id':'cell_id','canonical_sample_id':'sample'}).copy()
            # Source sample is retained explicitly when a canonical sample field shadows it.
            if cells.columns.duplicated().any():cells=cells.loc[:,~cells.columns.duplicated()].copy()
            cells['raw_total_counts']=counts;cells['raw_n_genes']=ngenes
        cells['study']=study;cells=statscores(x,genes,cells)
        tab(f'{study}_cells_prepared.tsv.gz',cells)
    tab('local_input_inventory.tsv',pd.DataFrame(manifest));tab('gene_coverage.tsv',pd.DataFrame(coverage))
    dump('provenance/input_before.json',manifest)

def transform(x,genes,space):
    fi=pd.Index(genes).get_indexer(space['feature_names'])
    pc=np.empty((len(x),len(space['pca_components'])),dtype='float32')
    for start in range(0,len(x),3000):
        stop=min(start+3000,len(x));z=np.array(x[start:stop][:,fi],dtype='float32');miss=~np.isfinite(z)
        z=np.where(miss,space['feature_means'][None,:],z).astype('float32')
        z-=space['feature_means'];z/=space['feature_scales'];np.clip(z,-5,5,out=z)
        pc[start:stop]=(z-space['pca_mean'])@space['pca_components'].T
    return pc
def rescore():
    genes=json.loads((OUT/'cache/genes.json').read_text())
    for study in ['mgeo_div30']+STUDIES:
        path=OUT/'tables'/f'{study}_cells_prepared.tsv.gz';c=pd.read_csv(path,sep='\t')
        x=np.load(OUT/'cache'/f'{study}_expression.npy',mmap_mode='r')
        c=statscores(x,genes,c.drop(columns=[g for g in EXTRA if g in c]));tab(path.name,c)
        print('rescored frozen context',study,flush=True)
def query(knn,x,k=30):
    dd=[];nn=[]
    for s in range(0,len(x),5000):
        d,n=knn.kneighbors(x[s:s+5000],n_neighbors=k);dd.append(d);nn.append(n)
    return np.vstack(dd),np.vstack(nn)

def mapping(variant='primary'):
    suffix='' if variant=='primary' else '_omit_SLC6A1'
    genes=json.loads((OUT/'cache/genes.json').read_text());cells=pd.read_csv(OUT/'tables/mgeo_div30_cells_prepared.tsv.gz',sep='\t')
    if 'total_counts' not in cells:
        original=pd.read_csv(OLD/'cells.tsv.gz',sep='\t',usecols=['cell_id','total_counts','n_genes'])
        assert np.array_equal(original.cell_id,cells.cell_id)
        cells['total_counts']=original.total_counts;cells['n_genes']=original.n_genes
    essential=['cell_id','sample','cell_line','condition','study','target','early_pv_score','postmitotic','cycling','direct_pv_score','total_counts','n_genes','orig.ident','samples','process','celltype','cellLine','donor_batch','raw_total_counts','raw_n_genes','mge_score','progenitor_score','maturation_score','cycle_score','progression_score','postmitotic_continuous','postmitotic_proxy','SLC6A1','ERBB4','LHX6','SST','PVALB']
    cells=cells[[g for g in essential if g in cells]].copy()
    space=np.load(MAP/('frozen_reference_space.npz' if variant=='primary' else 'omit_SLC6A1_reference_space.npz'));x30=space['div30_pcs'];target=cells.target.to_numpy()
    # Every MGEO benchmark neighbor comes from a different sample.
    prior=OUT/'tables'/('mgeo_div30_mapped_cells'+suffix+'.tsv.gz')
    if prior.exists():
        cells=pd.read_csv(prior,sep='\t');held=cells.neighbor_target_resemblance.to_numpy();helddist=cells.distance_to_div30.to_numpy()
    else:
        held=np.zeros(len(cells));helddist=np.zeros(len(cells));cal=[]
        for sample,ii in cells.groupby('sample').indices.items():
            tr=np.flatnonzero(cells['sample'].to_numpy()!=sample);knn=NearestNeighbors(n_neighbors=30,n_jobs=2).fit(x30[tr])
            d,n=query(knn,x30[ii]);pred=target[tr[n]].mean(axis=1);held[ii]=pred;helddist[ii]=d.mean(axis=1)
            for threshold in [.4,.5,.6]:
                selected=pred>=threshold
                cal.append(dict(sample=sample,threshold=threshold,n=len(ii),original_target_fraction=target[ii].mean(),resemblance_fraction=selected.mean(),precision=target[ii][selected].mean() if selected.any() else np.nan,recovery=selected[target[ii]].mean(),auc=roc_auc_score(target[ii],pred)))
        tab('mgeo_held_sample_resemblance_calibration'+suffix+'.tsv',pd.DataFrame(cal));cells['neighbor_target_resemblance']=held;cells['distance_to_div30']=helddist
    limit=np.quantile(helddist,.95);cells['in_div30_distance_support']=helddist<=limit
    cells['majority_target_resemblance']=held>=.5
    tab('mgeo_div30_mapped_cells'+suffix+'.tsv.gz',cells)
    knn30=NearestNeighbors(n_neighbors=30,n_jobs=2).fit(x30)
    knn90=NearestNeighbors(n_neighbors=30,n_jobs=2).fit(space['div90_pcs'])
    ref=pd.read_csv(OLD/'tables/reference_cell_index.tsv',sep='\t')
    out90=pd.read_csv(OLD/'tables/reference_oof_predictions.tsv.gz',sep='\t').query("grouping=='sample'").neighbor_distance.quantile(.95)
    cov=pd.read_csv(OUT/'tables/gene_coverage.tsv',sep='\t');summ=[]
    for study in STUDIES:
        path=OUT/'tables'/f'{study}_mapped_cells{suffix}.tsv.gz'
        if path.exists():c=pd.read_csv(path,sep='\t')
        else:
            c=pd.read_csv(OUT/'tables'/f'{study}_cells_prepared.tsv.gz',sep='\t',usecols=lambda g:g in essential);x=np.load(OUT/'cache'/f'{study}_expression.npy',mmap_mode='r')
            pc=transform(x,genes,space);np.save(OUT/'cache'/f'{study}_pcs{suffix}.npy',pc)
            d,n=query(knn30,pc);c['neighbor_target_resemblance']=target[n].mean(axis=1);c['distance_to_div30']=d.mean(axis=1)
            c['in_div30_distance_support']=c.distance_to_div30<=limit;c['majority_target_resemblance']=c.neighbor_target_resemblance>=.5
            c['neighbor_antecedent_score']=cells.early_pv_score.to_numpy()[n].mean(axis=1)
            d90,n90=query(knn90,pc);c['distance_to_div90']=d90.mean(axis=1);c['out_of_div90_reference']=c.distance_to_div90>out90
            c['neighbor_pv_score']=ref.pv_score.to_numpy()[n90].mean(axis=1)
            for state in sorted(ref.state.unique()):c['neighbor_state_'+state]=(ref.state.to_numpy()[n90]==state).mean(axis=1)
        for level,col in [('all',None),('sample','sample')]:
            groups=[('all',np.arange(len(c)))] if col is None else c.groupby(col).indices.items()
            for sample,ii in groups:
                a=c.iloc[ii];rec=dict(study=study,level=level,sample=sample,n_cells=len(a),mge_score_mean=a.mge_score.mean(),maturation_score_mean=a.maturation_score.mean(),slc6a1_positive_fraction=(a.SLC6A1>0).mean(),median_counts=a.raw_total_counts.median(),median_genes=a.raw_n_genes.median(),div90_out_of_reference_fraction=a.out_of_div90_reference.mean(),in_div30_distance_support_fraction=a.in_div30_distance_support.mean(),mean_target_neighbor_resemblance=a.neighbor_target_resemblance.mean())
                for thresh in [.4,.5,.6]:
                    sel=a.neighbor_target_resemblance>=thresh;rec[f'candidate_fraction_{thresh}']=sel.mean();rec[f'supported_candidate_fraction_{thresh}']=(sel&a.in_div30_distance_support).mean()
                summ.append(rec)
        tab(f'{study}_mapped_cells{suffix}.tsv.gz',c)
        print('mapped',study,len(c),'candidate',(c.neighbor_target_resemblance>=.5).mean(),'supported',((c.neighbor_target_resemblance>=.5)&c.in_div30_distance_support).mean(),flush=True)
    tab('external_reference_resemblance_frequency'+suffix+'.tsv',pd.DataFrame(summ))
    dump('provenance/mapping_definitions'+suffix+'.json',dict(variant=variant,neighbors=30,thresholds=[.4,.5,.6],div30_held_sample_distance_95=limit,div90_held_sample_distance_95=out90,feature_missingness='Missing source features imputed to frozen reference mean (zero standardized expression); reported per gene; no absent-gene detection interpreted.',normalization='External raw counts library-normalized across all source features to CP10K then log1p. MGEO uses existing frozen log1p normalized expression. Frozen TUBB3 zero contribution preserved in maturity score for every study.',interpretation='Neighbor fractions are developmental resemblance, never calibrated fate probabilities. Majority-target-neighbor candidate frequency is not the original target frequency. Cross-study shifts can remain even within distance support. No raw surface thresholds applied across studies.'))

def compare(variant='primary'):
    suffix='' if variant=='primary' else '_omit_SLC6A1'
    genes=json.loads((OUT/'cache/genes.json').read_text());c30=pd.read_csv(OUT/'tables'/('mgeo_div30_mapped_cells'+suffix+'.tsv.gz'),sep='\t');x30=np.load(OUT/'cache/mgeo_div30_expression.npy',mmap_mode='r')
    targets=np.flatnonzero(c30.target);metrics=['mge_score','maturation_score','postmitotic_continuous','neighbor_target_resemblance']
    a=c30[metrics].to_numpy();mu=a[targets].mean(axis=0);sd=a[targets].std(axis=0);a=(a-mu)/sd
    zmu=np.full(x30.shape[1],np.nan);zsd=np.full(x30.shape[1],np.nan);assayed=np.isfinite(x30[0])
    for start in range(0,x30.shape[1],40):
        stop=min(start+40,x30.shape[1]);ix=np.arange(start,stop)[assayed[start:stop]]
        if len(ix):zmu[ix]=x30[:,ix].mean(axis=0);zsd[ix]=x30[:,ix].std(axis=0)
    zsd[zsd<.05]=.05
    effects=[];balance=[];matches=[];individual=[];source_rows=[]
    for study in STUDIES:
        c=pd.read_csv(OUT/'tables'/f'{study}_mapped_cells{suffix}.tsv.gz',sep='\t');x=np.load(OUT/'cache'/f'{study}_expression.npy',mmap_mode='r')
        b=(c[metrics].to_numpy()-mu)/sd
        for sample,ii in c.groupby('sample').indices.items():
            # Matching restricted to positively resembling query cells with local DIV30 support.
            ii=ii[(c.neighbor_target_resemblance.to_numpy()[ii]>=.4)&c.in_div30_distance_support.to_numpy()[ii]]
            if study=='walsh':ii=ii[c['orig.ident'].to_numpy()[ii]=='vFB']
            if not len(ii):
                matches.append(dict(study=study,sample=sample,candidate_supported_n=0,matched_external_n=0,matched_fraction=0.,unique_mgeo_cells=0,max_axis_caliper=.5));continue
            matched=np.full(len(ii),-1,dtype=int)
            for flag in [False,True]:
                refs=targets[c30.postmitotic_proxy.to_numpy()[targets]==flag]
                qq=np.flatnonzero(c.postmitotic_proxy.to_numpy()[ii]==flag)
                if not len(refs) or not len(qq):continue
                model=NearestNeighbors(n_neighbors=1,n_jobs=2).fit(a[refs]);d,n=model.kneighbors(b[ii[qq]]);matched[qq]=refs[n[:,0]]
            delta=np.abs(b[ii]-a[np.maximum(matched,0)]);ok=(matched>=0)&(delta<=.5).all(axis=1)
            ei=ii[ok];mi=matched[ok]
            matches.append(dict(study=study,sample=sample,candidate_supported_n=len(ii),matched_external_n=len(ei),matched_fraction=float(ok.mean()),unique_mgeo_cells=len(np.unique(mi)),max_axis_caliper=.5))
            if len(ei)<30:continue
            assert np.array_equal(c.postmitotic_proxy.to_numpy()[ei],c30.postmitotic_proxy.to_numpy()[mi])
            for mgeo_sample in c30['sample'].unique():
                mask=c30['sample'].to_numpy()[mi]==mgeo_sample
                source_rows.append(dict(study=study,external_sample=sample,mgeo_sample=mgeo_sample,n_pairs=int(mask.sum()),fraction_pairs=mask.mean(),unique_mgeo_cells=len(np.unique(mi[mask]))))
            for j,metric in enumerate(metrics):balance.append(dict(study=study,sample=sample,metric=metric,n_pairs=len(ei),before_standardized_mean_difference=b[ii,j].mean()-a[targets,j].mean(),after_standardized_mean_difference=(b[ei,j]-a[mi,j]).mean(),mean_absolute_pair_difference=np.abs(b[ei,j]-a[mi,j]).mean()))
            balance.append(dict(study=study,sample=sample,metric='exact_harmonized_postmitotic_status',n_pairs=len(ei),before_standardized_mean_difference=c.postmitotic_proxy.iloc[ii].mean()-c30.postmitotic_proxy.iloc[targets].mean(),after_standardized_mean_difference=0.,mean_absolute_pair_difference=0.))
            for internal,external in [('total_counts','raw_total_counts'),('n_genes','raw_n_genes')]:
                aa=np.log1p(c30[internal].to_numpy());bb=np.log1p(c[external].to_numpy());scale=aa[targets].std()
                balance.append(dict(study=study,sample=sample,metric='unmatched_QC_log1p_'+internal,n_pairs=len(ei),before_standardized_mean_difference=(bb[ii].mean()-aa[targets].mean())/scale,after_standardized_mean_difference=(bb[ei]-aa[mi]).mean()/scale,mean_absolute_pair_difference=np.abs(bb[ei]-aa[mi]).mean()/scale))
            for panel,pgenes in PROGRAMS.items():
                ix=[genes.index(g) for g in pgenes if np.isfinite(x[0,genes.index(g)]) and np.isfinite(x30[0,genes.index(g)])]
                if len(ix)<3:continue
                mg=np.clip((x30[np.ix_(mi,ix)]-zmu[ix])/zsd[ix],-3,3).mean(axis=1)
                ex=np.clip((x[np.ix_(ei,ix)]-zmu[ix])/zsd[ix],-3,3).mean(axis=1)
                diff=mg-ex
                effects.append(dict(study=study,sample=sample,panel=panel,n_pairs=len(ei),n_genes=len(ix),mgeo_minus_external_mean_score=diff.mean(),mgeo_score_mean=mg.mean(),external_score_mean=ex.mean(),interpretation='Descriptive paired expression contrast; protocol and assay confounded; matched cells are not independent replicates'))
            for gene in sum(PROGRAMS.values(),[]):
                j=genes.index(gene);mg=x30[mi,j];ex=x[ei,j]
                if not np.isfinite(ex).any():continue
                individual.append(dict(study=study,sample=sample,gene=gene,n_pairs=len(ei),mgeo_mean=mg.mean(),external_mean=ex.mean(),mgeo_positive_fraction=(mg>0).mean(),external_positive_fraction=(ex>0).mean(),mgeo_minus_external_log1p_mean=(mg-ex).mean()))
    tab('matched_state_overlap'+suffix+'.tsv',pd.DataFrame(matches));tab('matched_state_balance'+suffix+'.tsv',pd.DataFrame(balance,columns=['study','sample','metric','n_pairs','before_standardized_mean_difference','after_standardized_mean_difference','mean_absolute_pair_difference']));tab('matched_program_effects_by_sample'+suffix+'.tsv',pd.DataFrame(effects));tab('matched_genes_by_sample'+suffix+'.tsv',pd.DataFrame(individual,columns=['study','sample','gene','n_pairs','mgeo_mean','external_mean','mgeo_positive_fraction','external_positive_fraction','mgeo_minus_external_log1p_mean']))
    tab('matched_mgeo_sample_composition'+suffix+'.tsv',pd.DataFrame(source_rows,columns=['study','external_sample','mgeo_sample','n_pairs','fraction_pairs','unique_mgeo_cells']))
    e=pd.DataFrame(effects)
    if not e.empty:
        summary=e.groupby(['study','panel']).agg(n_samples=('sample','nunique'),sample_mean_effect=('mgeo_minus_external_mean_score','mean'),sample_min_effect=('mgeo_minus_external_mean_score','min'),sample_max_effect=('mgeo_minus_external_mean_score','max'),fraction_samples_mgeo_higher=('mgeo_minus_external_mean_score',lambda z:np.mean(z>0))).reset_index()
        summary['within_study_rank']=summary.groupby('study').sample_mean_effect.rank(ascending=False,method='min');tab('matched_program_rankings'+suffix+'.tsv',summary)
    else:tab('matched_program_rankings'+suffix+'.tsv',pd.DataFrame(columns=['study','panel','n_samples','sample_mean_effect','sample_min_effect','sample_max_effect','fraction_samples_mgeo_higher','within_study_rank']))
    dump('provenance/matching_definitions'+suffix+'.json',dict(variant=variant,axes=metrics,exact_status_match='Harmonized postmitotic proxy=(frozen maturation>frozen progenitor) AND frozen cycle_score<0.5, computed identically in MGEO and external. This is not the frozen original phase-derived postmitotic flag and never changes original targets.',reference='Original DIV30 target cells only; external candidates require >=0.4 target-neighbor resemblance and DIV30-distance support',caliper='Every matched axis within 0.5 DIV30-target standard deviations; 1 nearest target with replacement per external cell and exact postmitotic-proxy status',expression='Each gene standardized to all DIV30 means/SD, SD floor0.05, clipped[-3,3], panel mean; same available genes both studies',inference='Sample-level descriptive ranges only; n cell pairs does not create independent replicates. Poor support and unmatched cells excluded, not imputed. Program scores overlap biology of matching and remaining chemistry effects; not causal or proof of unique MGEO programs.'))

def slc_diagnostic():
    """Bounded stratified diagnostic, never an external cohort frequency estimate."""
    genes=json.loads((OUT/'cache/genes.json').read_text());space=np.load(MAP/'omit_SLC6A1_reference_space.npz')
    c30=pd.read_csv(OUT/'tables/mgeo_div30_mapped_cells.tsv.gz',sep='\t');p30=space['div30_pcs'];target=c30.target.to_numpy();rng=np.random.default_rng(3090)
    distances=[];cal=[]
    for sample,ix in c30.groupby('sample').indices.items():
        ix=rng.choice(ix,min(500,len(ix)),replace=False);tr=np.flatnonzero(c30['sample'].to_numpy()!=sample)
        model=NearestNeighbors(n_neighbors=30,n_jobs=2).fit(p30[tr]);d,n=query(model,p30[ix]);distances.extend(d.mean(axis=1))
        cal.append(dict(sample=sample,n_queries=len(ix),mean_distance=d.mean(),majority_target_fraction=(target[tr[n]].mean(axis=1)>=.5).mean()))
    limit=float(np.quantile(distances,.95));tab('omit_SLC6A1_subsample_mgeo_calibration.tsv',pd.DataFrame(cal))
    model=NearestNeighbors(n_neighbors=30,n_jobs=2).fit(p30);datasets={study:pd.read_csv(OUT/'tables'/f'{study}_mapped_cells.tsv.gz',sep='\t') for study in STUDIES};ngroups=sum(x['sample'].nunique() for x in datasets.values());per_sample=min(250,10000//ngroups)
    rows=[];cells=[]
    for study,c in datasets.items():
        source=np.load(OUT/'cache'/f'{study}_expression.npy',mmap_mode='r')
        for sample,ix in c.groupby('sample').indices.items():
            ix=np.sort(rng.choice(ix,min(per_sample,len(ix)),replace=False));pc=transform(source[ix],genes,space);d,n=query(model,pc);score=target[n].mean(axis=1);dist=d.mean(axis=1)
            q=pd.DataFrame(dict(study=study,sample=sample,cell_id=c.cell_id.iloc[ix].to_numpy(),primary_target_resemblance=c.neighbor_target_resemblance.iloc[ix].to_numpy(),no_slc_target_resemblance=score,no_slc_distance=dist,no_slc_supported=dist<=limit));cells.append(q)
            rows.append(dict(study=study,sample=sample,n_sampled=len(ix),primary_majority_n=int((q.primary_target_resemblance>=.5).sum()),no_slc_majority_n=int((score>=.5).sum()),no_slc_candidate_0_4_n=int((score>=.4).sum()),no_slc_supported_n=int((dist<=limit).sum()),no_slc_supported_candidate_0_4_n=int(((score>=.4)&(dist<=limit)).sum()),max_no_slc_resemblance=score.max(),mean_resemblance_change=(score-q.primary_target_resemblance).mean()))
    tab('omit_SLC6A1_subsample_diagnostic.tsv',pd.DataFrame(rows));tab('omit_SLC6A1_subsample_cells.tsv.gz',pd.concat(cells,ignore_index=True))
    dump('provenance/omit_SLC6A1_subsample_diagnostic.json',dict(seed=3090,per_external_sample=per_sample,total_external_queries=sum(x['n_sampled'] for x in rows),div30_held_sample_calibration_queries=len(distances),recalibrated_div30_distance_95=limit,interpretation='Stratified capped subsample to test rescue after SLC6A1 feature omission; not a full-cohort frequency estimate. Original target fixed. Full external remapping and program matching unnecessary if no supported candidate rescued.'))

def import_samarasinghe():
    folder=OUT.parent/'external_extra/samarasinghe';provenance=json.loads((folder/'provenance.json').read_text())
    genes=json.loads((OUT/'cache/genes.json').read_text());sourcegenes=json.loads((folder/'genes.json').read_text());src=np.load(folder/'expression.npy',mmap_mode='r');c=pd.read_csv(folder/'cells.tsv.gz',sep='\t')
    assert src.shape==(len(c),len(sourcegenes));assert 'Genotype' in c and 'cell_id' in c
    selected=np.flatnonzero(c.Genotype.eq('Ctrl'));assert len(selected)==26935, 'Unexpected control cohort: reconcile source rather than silently change selection.'
    order=pd.Index(sourcegenes).get_indexer(genes);assert (order>=0).all()
    path=OUT/'cache/samarasinghe_2021_expression.npy';x=np.lib.format.open_memmap(path,mode='w+',dtype='float32',shape=(len(selected),len(genes)))
    for start in range(0,len(selected),3000):
        stop=min(start+3000,len(selected));x[start:stop]=src[np.ix_(selected[start:stop],order)]
    x.flush();q=c.iloc[selected].copy().reset_index(drop=True);q['study']='samarasinghe_2021';q=statscores(x,genes,q);tab('samarasinghe_2021_cells_prepared.tsv.gz',q)
    cov=pd.read_csv(OUT/'tables/gene_coverage.tsv',sep='\t');cov=cov[cov.study!='samarasinghe_2021'];rr=pd.DataFrame(dict(study='samarasinghe_2021',gene=genes,present=np.isfinite(x[0]),source_feature_index=order));tab('gene_coverage.tsv',pd.concat([cov,rr],ignore_index=True))
    inventory=pd.read_csv(OUT/'tables/local_input_inventory.tsv',sep='\t');inventory=inventory[inventory.study!='samarasinghe_2021']
    source=PROJECT/'results/samarasinghe_2021_zenodo_processed_object/samarasinghe_2021_zenodo_seurat.rds'
    rec=dict(study='samarasinghe_2021',path=str(source),n_cells=len(q),n_genes=int(provenance['counts_shape'][0]),size_bytes=source.stat().st_size,mtime_ns=source.stat().st_mtime_ns)
    tab('local_input_inventory.tsv',pd.concat([inventory,pd.DataFrame([rec])],ignore_index=True))
    manifest=json.loads((OUT/'provenance/input_before.json').read_text());manifest=[r for r in manifest if r['study']!='samarasinghe_2021']+[rec];dump('provenance/input_before.json',manifest)
    dump('provenance/samarasinghe_import.json',dict(source_export=provenance,source_cells=len(c),selected_control_cells=len(q),selection='Genotype==Ctrl per pre-existing project convention; retainsD56unfused,D70docked,D100docked separately; no inferred ventral/dorsal cells',missing_gene_rule='NaN in source subset remainsNaN; no anonymous rawCSVgene identities used',mapping_features=int(rr.set_index('gene').loc[np.load(MAP/'frozen_reference_space.npz')['feature_names'],'present'].sum())))
    print('Imported Samarasinghe controls',len(q),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=['prepare','rescore','map','compare','slc-diagnostic','import-samarasinghe'],required=True);p.add_argument('--variant',choices=['primary','omit_SLC6A1'],default='primary');a=p.parse_args()
    if a.stage in ['map','compare']:{'map':mapping,'compare':compare}[a.stage](a.variant)
    else:{'prepare':prepare,'rescore':rescore,'slc-diagnostic':slc_diagnostic,'import-samarasinghe':import_samarasinghe}[a.stage]()
