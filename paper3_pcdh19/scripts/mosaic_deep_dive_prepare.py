"""Exploratory gate/identity sensitivity analysis with auditable raw pseudobulks.

Scores and marker-defined strata are exploratory descriptions, not new genotypes
or a locked cell taxonomy. No-probe cells never receive an inferred allele label.
"""
import csv
import gzip
import hashlib
import itertools
import json
from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
from collections import Counter
import numpy as np
import h5py
from scipy import sparse

ROOT=Path(__file__).resolve().parents[1]
OUT=RUNTIME_ROOT/'results/mosaic_deep_dive_20261031'
BASE=Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results')
H5=BASE/'primary_processing/06_technical_sample_batch_diagnostics/06_technical_sample_batch_diagnostics_20260903_114110_6087eeb/objects/pcdh19_step06_unintegrated_diagnostics.h5ad'
ASSIGN=BASE/'hicat/03_full_data_consensus_benchmark/hicat_consensus_review_20260909_142348_efd9e13d/outputs/tables/all_cell_assignments.tsv.gz'
MODEL=BASE/'step_05_pcdh19_logistic_regression_baseline/count_informed_vs_binary_validation/step_05_count_informed_full_fit_coefficients.tsv'
GATES=['all','informative','A_positive','Aneg_BC1','Aneg_BC2','Aneg_BC3','frozen_WT','frozen_KO',
       'WTtail05','WTtail10','WTtail20','KOtail05','KOtail10','KOtail20']
STRATA=['all_cells','neural_vascular_low','vascular_marked','mixed_profile','C0001','C0006',
        'C0002','C0003','C0011','C0001_neural_vascular_low','neural_depth_trimmed']
SAMPLES=['15662-JZ-'+str(i) for i in range(1,13)]
FEMALES=SAMPLES[3:9]
def read(p):
    with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:yield from csv.DictReader(f,delimiter='\t')
def strings(n):
    if isinstance(n,h5py.Group):
        cats=strings(n['categories']);return np.array([cats[i] if i>=0 else '' for i in n['codes'][:]])
    return np.array([v.decode() if isinstance(v,bytes) else str(v) for v in n[:]])
def write(name,rows):
    with (OUT/name).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)

def main():
    OUT.mkdir(exist_ok=True,parents=True)
    assert not (OUT/'PREPARATION_COMPLETE.json').exists(),'Preparation already complete; reuse outputs.'
    cfg=json.loads((ROOT/'config/hicat_validation_programs.json').read_text())
    programs={r['name']:r['genes'] for r in cfg['programs'] if r['name'] not in ['mitochondrial','ribosomal','S_phase','G2M']}
    programs.update(neural_core=['Sox2','Fabp7','Ascl1','Dlx1','Dlx2','Gad1','Gad2','Lhx6','Dcx','Tubb3','Sox6','Erbb4'],
      adhesion_synapse=['Cdh2','Cdh8','Cdh11','Pcdh7','Pcdh8','Pcdh9','Pcdh10','Pcdh17','Nrxn1','Nrxn2','Nrxn3','Nlgn1','Nlgn2','Nlgn3','Cntn1','Cntn2','Cntnap2','Dscam','L1cam','Ncam1','Itgb1'],
      neurite_guidance=['Reln','Dab1','Slit1','Slit2','Slit3','Robo1','Robo2','Sema3a','Sema3c','Sema3f','Nrp1','Nrp2','Plxna1','Plxna2','Epha4','Epha5','Ephb2','Efna1','Efnb2'],
      immune=['Aif1','P2ry12','Cx3cr1','Tyrobp','C1qa','C1qb','Csf1r'],
      oligodendroglial=['Olig1','Olig2','Sox10','Pdgfra','Plp1','Mbp'])
    # Freeze the exploratory definitions before looking at their results.
    plan=dict(gates=GATES,strata=STRATA,programs=programs,
      vascular_rule='>=3 detected endothelial markers and >=10 combined marker UMIs, OR >=3 pericyte markers and >=10 combined marker UMIs',
      neural_rule='>=3 detected neural_core genes',depth_rule='exclude upper 10% raw total-UMI cells within each original sample',
      percentile_rule='within each sample among informative non-000 cells, retaining all cutoff ties',
      expression_scope='female samples only; sample pseudobulk; broad exploratory expression filtering; all tested results retained')
    (OUT/'EXPLORATORY_PLAN.json').write_text(json.dumps(plan,indent=2))
    coef={r['term']:float(r['coefficient']) for r in read(MODEL)}
    assert hashlib.sha256(MODEL.read_bytes()).hexdigest()=='2596c0c2d7a6fca539b7af3b2ea19792dc884f7c5a874f3c72322d2e6ea3215a'
    assignments={r['cell_id']:r['diagnostic_cluster'] for r in read(ASSIGN)}
    probes={};cutoff_rows=[]
    sources=[MODEL,ASSIGN]
    for s in SAMPLES:
        path=BASE/'pcdh19_probe_audit/per_sample'/s/'pcdh19_probe_patterns.tsv';sources.append(path)
        vals={}
        for r in read(path):
            a,b,c=[int(r[k]) for k in ['A_UMI','B_UMI','C_UMI']]
            z=coef['intercept']+coef['A_detected']*(a>0)+coef['B_UMI']*b+coef['C_UMI']*c
            p=round(float(1/(1+np.exp(-z))),12)
            vals[s+'_'+r['barcode']]=(a,b,c,p)
        informative=np.array([v[3] for v in vals.values() if sum(v[:3])>0]);informative.sort()
        cutoffs={}
        for pct in [5,10,20]:
            k=int(np.ceil(len(informative)*pct/100));lo=float(informative[k-1]);hi=float(informative[-k]);cutoffs[pct]=(lo,hi)
            cutoff_rows.append(dict(sample=s,percentile=pct,informative_cells=len(informative),WT_max_PKO=lo,KO_min_PKO=hi,
                                    WT_cells_including_ties=int(sum(informative<=lo)),KO_cells_including_ties=int(sum(informative>=hi))))
        for cell,(a,b,c,p) in vals.items():
            info=a+b+c>0
            gates=[True,info,a>0,a==0 and b+c>=1,a==0 and b+c>=2,a==0 and b+c>=3,
                   info and p<=.301037619832,info and p>=.911554020713]
            gates += [info and p<=cutoffs[pct][0] for pct in [5,10,20]]
            gates += [info and p>=cutoffs[pct][1] for pct in [5,10,20]]
            probes[cell]=(a,b,c,p,np.array(gates,dtype=bool))
        print('PROBES_READY',s,len(vals),flush=True)
    write('percentile_cutoffs.tsv',cutoff_rows)
    libs=list(itertools.product(FEMALES,GATES,STRATA));libidx={key:i for i,key in enumerate(libs)}
    cell_counts=Counter();features=[];marker_rows=[]
    score_sums=Counter();score_counts=Counter();composition=Counter()
    with h5py.File(H5,'r') as f:
        ids=strings(f['obs/cell_id']);genes=strings(f['var/gene_id']);symbols=strings(f['var/gene_symbol'])
        total=f['obs/total_counts'][:];n_genes=f['obs/n_genes_by_counts'][:];mt=f['obs/pct_counts_mt'][:]
        sample=np.array([c.split('_')[0] for c in ids]);cluster=np.array([assignments[c] for c in ids])
        assert len(ids)==446349 and len(genes)==19071
        symbol_index={s.lower():i for i,s in enumerate(symbols)}
        panel_indices=sorted(set(symbol_index[g.lower()] for gs in programs.values() for g in gs if g.lower() in symbol_index))
        panel_lookup={g:i for i,g in enumerate(panel_indices)}
        program_indices={p:[panel_lookup[symbol_index[g.lower()]] for g in gs if g.lower() in symbol_index] for p,gs in programs.items()}
        gene_mapping=[dict(program=p,gene=g,present=g.lower() in symbol_index) for p,gs in programs.items() for g in gs]
        write('marker_program_gene_coverage.tsv',gene_mapping)
        trim={s:float(np.quantile(total[sample==s],.9)) for s in SAMPLES}
        counts=np.zeros((len(genes),len(libs)),dtype=np.int64)
        x=f['X'];assert x.attrs['encoding-type']=='csr_matrix' and x['data'].dtype.kind in 'iu'
        ptr=x['indptr'][:]
        object_path=OUT/'per_cell_marker_diagnostics.h5'
        with h5py.File(object_path,'w') as side:
            string=h5py.string_dtype('utf-8')
            for name,values in [('cell_id',ids),('sample',sample),('cluster',cluster),('marker_gene',symbols[panel_indices]),('gene_id',genes)]:
                side.create_dataset(name,data=np.asarray(values,dtype=object),dtype=string)
            side.create_dataset('total_counts',data=total);side.create_dataset('n_genes',data=n_genes);side.create_dataset('pct_counts_mt',data=mt)
            rawpanel=side.create_dataset('marker_raw_counts',shape=(len(ids),len(panel_indices)),dtype='int32',compression='gzip',chunks=(1000,len(panel_indices)))
            scores=side.create_dataset('program_mean_log1p_CP10K',shape=(len(ids),len(programs)),dtype='float32',compression='gzip')
            side.create_dataset('program_names',data=np.asarray(list(programs),dtype=object),dtype=string)
            prob=side.create_dataset('P_KO',shape=(len(ids),),dtype='float64')
            gate_store=side.create_dataset('gate_membership',shape=(len(ids),len(GATES)),dtype='bool',compression='gzip')
            stratum_store=side.create_dataset('stratum_membership',shape=(len(ids),len(STRATA)),dtype='bool',compression='gzip')
            side.create_dataset('gate_names',data=np.asarray(GATES,dtype=object),dtype=string)
            side.create_dataset('stratum_names',data=np.asarray(STRATA,dtype=object),dtype=string)
            for start in range(0,len(ids),2000):
                end=min(start+2000,len(ids));lo,hi=int(ptr[start]),int(ptr[end]);data=x['data'][lo:hi]
                assert np.all(data>=0)
                mat=sparse.csr_matrix((data,x['indices'][lo:hi],ptr[start:end+1]-lo),shape=(end-start,len(genes)))
                assert np.array_equal(np.asarray(mat.sum(axis=1)).ravel(),total[start:end])
                panel=mat[:,panel_indices].toarray();rawpanel[start:end]=panel
                norm=np.log1p(panel*10000/total[start:end,None])
                sc=np.column_stack([norm[:,ix].mean(axis=1) if ix else np.zeros(end-start) for ix in program_indices.values()]);scores[start:end]=sc
                ep=panel[:,program_indices['endothelial']];pc=panel[:,program_indices['pericyte']];nc=panel[:,program_indices['neural_core']]
                vascular=((ep>0).sum(axis=1)>=3)&(ep.sum(axis=1)>=10) | (((pc>0).sum(axis=1)>=3)&(pc.sum(axis=1)>=10))
                neural=(nc>0).sum(axis=1)>=3
                gm=np.array([probes[c][4] for c in ids[start:end]])
                sm=np.column_stack([np.ones(end-start,dtype=bool),neural&~vascular,vascular,neural&vascular,
                    *[cluster[start:end]==c for c in ['C0001','C0006','C0002','C0003','C0011']],
                    (cluster[start:end]=='C0001')&neural&~vascular,
                    neural&~vascular&(total[start:end]<=np.array([trim[s] for s in sample[start:end]]))])
                gate_store[start:end]=gm;stratum_store[start:end]=sm;prob[start:end]=[probes[c][3] for c in ids[start:end]]
                row_ids=[];col_ids=[];lib_ids=[];local={}
                for s in set(sample[start:end]):
                    samplemask=sample[start:end]==s
                    for gi,g in enumerate(GATES):
                        gate_mask=samplemask&gm[:,gi]
                        if not np.any(gate_mask):continue
                        key=(s,g);n=int(sum(gate_mask));score_counts[key]+=n
                        for pi,p in enumerate(programs):score_sums[s,g,p]+=float(sc[gate_mask,pi].sum())
                        for name,values in [('total_counts',total[start:end]),('n_genes',n_genes[start:end]),('pct_counts_mt',mt[start:end])]:score_sums[s,g,name]+=float(values[gate_mask].sum())
                        for si,stratum in enumerate(STRATA):
                            chosen=gate_mask&sm[:,si];n=int(sum(chosen));composition[s,g,stratum]+=n
                            if s not in FEMALES or not n:continue
                            lib=libidx[s,g,stratum];cell_counts[lib]+=n
                            if lib not in local:local[lib]=len(local)
                            ix=np.flatnonzero(chosen);row_ids.extend([local[lib]]*len(ix));col_ids.extend(ix)
                if local:
                    membership=sparse.csr_matrix((np.ones(len(row_ids),dtype=np.int64),(row_ids,col_ids)),shape=(len(local),end-start))
                    aggregate=(membership@mat).toarray().T
                    order=sorted(local,key=local.get);counts[:,order]+=aggregate
                if start%20000==0:print('RAW_ROWS_COMPLETE',end,flush=True)
    metadata=[]
    for i,(s,g,st) in enumerate(libs):
        metadata.append(dict(library='__'.join([s,g,st]),sample=s,design_group='WT_F' if s in FEMALES[:3] else 'HET_F',gate=g,stratum=st,cells=cell_counts[i],UMI=int(counts[:,i].sum())))
    write('pseudobulk_metadata.tsv',metadata)
    np.savez_compressed(OUT/'raw_pseudobulk_counts.npz',counts=counts,genes=genes,symbols=symbols)
    # One compact counts file per stratum supports independent SLURM array fitting.
    for st in STRATA:
        ix=[i for i,lib in enumerate(libs) if lib[2]==st]
        with gzip.open(OUT/(st+'_counts.tsv.gz'),'wt') as f:
            w=csv.writer(f,delimiter='\t');w.writerow(['gene_id','gene_symbol']+[metadata[i]['library'] for i in ix])
            for j in range(len(genes)):w.writerow([genes[j],symbols[j]]+list(counts[j,ix]))
    comp=[];score=[]
    for (s,g,st),n in sorted(composition.items()):
        den=score_counts[s,g];comp.append(dict(sample=s,gate=g,stratum=st,cells=n,gate_total=den,fraction=n/den))
    for (s,g,p),summed in sorted(score_sums.items()):score.append(dict(sample=s,gate=g,program=p,mean=summed/score_counts[s,g],cells=score_counts[s,g]))
    write('sample_gate_stratum_composition.tsv',comp);write('sample_gate_program_scores.tsv',score)
    assert [cell_counts[libidx[s,'frozen_KO','all_cells']] for s in FEMALES[3:]]==[26,10,103]
    assert [cell_counts[libidx[s,'frozen_WT','all_cells']] for s in FEMALES[3:]]==[467,1069,2411]
    (OUT/'provenance.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},indent=2))
    (OUT/'PREPARATION_COMPLETE.json').write_text(json.dumps(dict(cells=len(ids),genes=len(genes),libraries=len(libs),marker_genes=len(panel_indices),strata=STRATA,gates=GATES,all_raw_row_totals_verified=True,frozen_HET_counts_reproduced=True),indent=2))
    print('PREPARATION_COMPLETE',flush=True)
if __name__=='__main__':main()
