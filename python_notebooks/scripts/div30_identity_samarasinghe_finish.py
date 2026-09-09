#!/usr/bin/env python3
"""Convert the bounded official-object CSC export into a row-aligned NPY."""
from pathlib import Path
import json
import shutil
import numpy as np
import pandas as pd
from scipy import sparse

ROOT=Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_antecedent_identity_v1')
OUT=ROOT/'external_extra/samarasinghe'
genes=json.loads((ROOT/'external/cache/genes.json').read_text())
provenance=json.loads((OUT/'provenance.json').read_text())
availability=pd.read_csv(OUT/'gene_availability.tsv',sep='\t')
present=json.loads((OUT/'present_genes.json').read_text())
assert availability.gene.tolist()==genes
assert availability.loc[availability.assayed,'gene'].tolist()==present
shape=tuple(provenance['selected_shape'])
values=np.fromfile(OUT/'expression_values.f32',dtype='<f4')
indices=np.fromfile(OUT/'expression_indices.i32',dtype='<i4')
pointers=np.fromfile(OUT/'expression_indptr.i32',dtype='<i4')
assert len(values)==len(indices)==provenance['selected_nnz'] and len(pointers)==shape[1]+1 and pointers[-1]==len(values)
assert np.all(np.diff(pointers)>=0) and np.isfinite(values).all() and (values>=0).all()
matrix=sparse.csc_matrix((values,indices,pointers),shape=shape,copy=False)
cells=pd.read_csv(OUT/'cells.tsv.gz',sep='\t',float_precision='round_trip')
assert len(cells)==shape[1] and cells.cell_id.is_unique
assert cells.cell_id.equals(cells.source_cell_id)
assert (cells.raw_total_counts>0).all() and (cells.raw_n_genes>0).all()
expr=np.lib.format.open_memmap(OUT/'expression.npy',mode='w+',dtype='float32',shape=(len(cells),len(genes)))
expr[:]=np.nan
cols=np.flatnonzero(availability.assayed)
for start in range(0,len(cells),2000):
    stop=min(start+2000,len(cells));expr[start:stop,cols]=matrix[:,start:stop].T.toarray()
expr.flush()
(OUT/'genes.json').write_text(json.dumps(genes,indent=2)+'\n')
summary=[]
for gene in ['SLC6A1','NKX2-1','LHX6','SOX6','ERBB4','MEF2C','KCNC1','KCNC2','TAC1','PVALB','SST','GRIA2','OPCML','NOTCH1']:
    if gene in genes and gene in present:
        x=expr[:,genes.index(gene)]
        for label,mask in [('all',np.ones(len(cells),bool)),('Ctrl',cells.Genotype.eq('Ctrl').to_numpy())]:
            summary.append(dict(group=label,gene=gene,n=int(mask.sum()),detected_fraction=float((x[mask]>0).mean()),mean_log1p_cp10k=float(x[mask].mean())))
pd.DataFrame(summary).to_csv(OUT/'diagnostic_gene_summary.tsv',sep='\t',index=False)
group_cols=[c for c in ['orig.ident','Time','Genotype'] if c in cells]
cells.groupby(group_cols,dropna=False).size().rename('n_cells').reset_index().to_csv(OUT/'source_sample_group_time_counts.tsv',sep='\t',index=False)
cache_path=ROOT.parents[0]/'cross_study_marker_expression/cross_study_marker_expression_v12_pv_precursors_final_candidate_plus_vipr2_slc6a1/tables/per_study/samarasinghe_2021_marker_expression.tsv.gz'
cached=pd.read_csv(cache_path,sep='\t')
assert cached.cell_id.equals(cells.cell_id) and cached['sample'].equals(cells['sample'])
assert cells.raw_total_counts.equals(cells.nCount_RNA) or np.array_equal(cells.raw_total_counts,cells.nCount_RNA)
assert np.array_equal(cells.raw_n_genes,cells.nFeature_RNA)
checks=[]
for gene in sorted(set(genes)&set(cached.columns)):
    a=np.asarray(expr[:,genes.index(gene)],np.float64);b=cached[gene].to_numpy(float)
    checks.append(dict(gene=gene,n_cells=len(a),max_abs_difference=float(np.max(np.abs(a-b))),detected_mismatch_n=int(np.sum((a>0)!=(b>0)))))
check=pd.DataFrame(checks);assert check.max_abs_difference.max()<1e-6 and check.detected_mismatch_n.sum()==0
check.to_csv(OUT/'validation_against_official_existing_gene_cache.tsv',sep='\t',index=False)
provenance.update(status='validated_row_by_gene_npy_ready',npy_shape=[len(cells),len(genes)],npy_dtype='float32',npy_axis0='same ordered cells.tsv.gz cell_id',npy_axis1='genes.json exact requested feature order',missing_features='NaN, not biological zero',cell_selection='All official processed cells exported; use Genotype==Ctrl for primary healthy-organoid comparison',validation=dict(csc_shapes_and_nnz=True,unique_original_cell_ids=True,finite_nonnegative_assayed_log_expression=True,raw_library_and_gene_counts_positive=True,source_size_mtime_unchanged=True))
provenance['validation'].update(original_official_cache_cell_order_and_samples=True,raw_library_totals_and_detected_genes_equal_original_metadata=True,existing_gene_cache_compared_genes=len(check),existing_gene_cache_max_abs_difference=float(check.max_abs_difference.max()),existing_gene_cache_detection_mismatches=int(check.detected_mismatch_n.sum()))
(OUT/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
for script in ['div30_identity_samarasinghe_export.R','div30_identity_samarasinghe_finish.py']:
    shutil.copy2(Path(__file__).parent/script,OUT/script)
print(json.dumps(dict(status=provenance['status'],shape=provenance['npy_shape'],assayed_genes=len(present),controls=int(cells.Genotype.eq('Ctrl').sum()),metadata=group_cols,source_samples=cells.groupby(group_cols,dropna=False).size().rename('n_cells').reset_index().to_dict('records')),default=str),flush=True)
