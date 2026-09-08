"""Execute the real full-population Allen consensus resource benchmark.

Stages are separate SLURM processes: prepare -> fit80 -> fullfit -> consensus.
Only ONE independently sampled 80% fit is authorized here. A second fit on
all cells is Allen's large-graph initialization, not a second parameter choice.
No function submits production iterations. All outputs are diagnostic and end
IN_REVIEW. Read FULL_DATA_CONSENSUS_BENCHMARK.md for fields and omissions.
"""
import argparse
import gc
from itertools import combinations
import json
from pathlib import Path
import resource
import time
import h5py
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
import transcriptomic_clustering as tc
from .consensus_core import Membership, map_heldout
from .fullscale_engine import FullScaleEngine, allocation_adapters
from .provenance import Progress, sha256, write_json, manifest


def decode(node):
    """Read required AnnData metadata without interpreting unrelated uns slots."""
    if isinstance(node,h5py.Group):
        categories=decode(node['categories']);codes=node['codes'][:]
        if np.any(codes<0):
            return np.array([categories[k] if k>=0 else '' for k in codes])
        return categories[codes]
    return node.asstr()[:] if node.dtype.kind in 'OSU' else node[:]


class Benchmark:
    """Own one immutable benchmark package; run exactly one requested stage."""

    def __init__(self, run_dir, stage):
        """Load frozen controls and create a separate single-writer progress log."""
        self.root=Path(run_dir)
        self.cfg=json.loads((self.root/'config/benchmark.json').read_text())
        self.stage=stage
        self.progress=Progress(self.root/'progress'/stage)
        self.start=time.monotonic()

    def verify_frozen(self):
        """Check executable/config/reference hashes before every stage."""
        checks=json.loads((self.root/'frozen_manifest.json').read_text())
        for relative,expected in checks.items():
            if sha256(self.root/relative)!=expected:
                raise ValueError('Frozen run content changed: '+relative)

    def prepare(self):
        """Validate approved counts and save one shared sparse normalized H5AD.

        X is float64 ln(1+CPM), never raw counts. Counts remain in the linked
        approved Step02 asset. HDF5 blocks keep preparation RAM independent of
        total nnz. No gene, sample or cell filtering is performed.
        """
        cfg=self.cfg;source=Path(cfg['source_path'])
        with self.progress.track('input.verify',source=str(source)) as event:
            status=pd.read_csv(source.parent.parent/'STEP_STATUS.tsv',sep='\t')
            if len(status)!=1 or status.iloc[0]['status']!='APPROVED' or status.iloc[0]['run_id']!=cfg['step02_run_id']:
                raise ValueError('Expected approved Step02 checkpoint')
            observed=sha256(source)
            if observed!=cfg['step02_sha256']:
                raise ValueError('Approved source checksum mismatch')
            event.update(sha256=observed)
        output=self.root/'inputs/full_log1p_cpm.h5ad'
        temporary=output.with_suffix('.partial.h5ad')
        with self.progress.track('input.normalize',formula='ln(1 + counts/cell_total * 1e6)',
                                 block_cells=cfg['block_cells'],batch_correction=False) as event:
            with h5py.File(source,'r') as src:
                x=src['X'];shape=tuple(map(int,x.attrs['shape']))
                if shape!=(cfg['expected_cells'],cfg['expected_genes']):
                    raise ValueError('Full source shape mismatch')
                obs_node=src['obs'];var_node=src['var']
                ids=decode(obs_node[obs_node.attrs['_index']])
                genes=decode(var_node[var_node.attrs['_index']])
                obs=pd.DataFrame(index=pd.Index(ids,name='cell_id'))
                # Metadata are saved for review; the fitting objects strip them.
                for key in ['technical_sample_id','submitted_sample_name','genotype','sex','design_group',
                            'total_counts','n_genes_by_counts','pct_counts_mt']:
                    if key in obs_node:obs[key]=decode(obs_node[key])
                if obs.technical_sample_id.nunique()!=cfg['expected_samples'] or not obs.index.is_unique:
                    raise ValueError('Sample count or unique-cell contract failed')
                var=pd.DataFrame({'gene_symbol':decode(var_node['gene_symbol'])},index=genes)
                obj=ad.AnnData(sparse.csr_matrix(shape,dtype=np.float64),obs=obs,var=var,dtype=np.float64)
                obj.uns['consensus_benchmark']=dict(status='DIAGNOSTIC_NOT_FINAL_K',
                    biological_context='dissected E14.5 mouse MGE',X_representation='natural log(1+CPM)',
                    count_source=str(source),count_source_sha256=observed,raw_counts_saved_here=False,
                    normalization_target=1000000,batch_correction=False,cell_cycle_regression=False,
                    annotations_used_for_fitting=False,full_pair_matrix_saved=False)
                obj.write_h5ad(temporary)
                write_json(self.root/'inputs/anndata_uns_inventory.json',dict(obj.uns))
                write_json(self.root/'inputs/anndata_slots_inventory.json',dict(shape=list(shape),X='float64 CSR ln(1+CPM)',
                    layers=[],raw=None,obsm=[],obsp=[],varm=[],varp=[],obs=list(obs.columns),var=list(var.columns),
                    uns=['consensus_benchmark']))
                pd.DataFrame({'cell_id':ids,'source_row':np.arange(shape[0]),
                              'technical_sample_id':obs.technical_sample_id.to_numpy()}).to_csv(
                    self.root/'inputs/full_cell_index.tsv.gz',sep='\t',index=False)
                obs.technical_sample_id.value_counts().sort_index().to_csv(self.root/'inputs/sample_counts.tsv',sep='\t')
                nnz=x['data'].shape[0];indptr=x['indptr'][:]
                if indptr[-1]!=nnz:raise ValueError('CSR metadata mismatch')
                with h5py.File(temporary,'r+') as dst:
                    del dst['X'];target=dst.create_group('X')
                    target.attrs.update({'encoding-type':'csr_matrix','encoding-version':'0.1.0','shape':shape})
                    target.create_dataset('indptr',data=indptr,compression='lzf')
                    out_data=target.create_dataset('data',shape=(nnz,),dtype='float64',chunks=True,compression='lzf')
                    out_indices=target.create_dataset('indices',shape=(nnz,),dtype='int32',chunks=True,compression='lzf')
                    for start in range(0,shape[0],cfg['block_cells']):
                        stop=min(start+cfg['block_cells'],shape[0]);left,right=int(indptr[start]),int(indptr[stop])
                        data=x['data'][left:right];indices=x['indices'][left:right]
                        if np.any(data<0):raise ValueError('Negative raw counts')
                        block=sparse.csr_matrix((data,indices,indptr[start:stop+1]-left),shape=(stop-start,shape[1]))
                        totals=np.asarray(block.sum(axis=1)).ravel()
                        if np.any(totals<=0):raise ValueError('Zero-count source cell')
                        normalized=block.astype(np.float64)
                        normalized.data *= np.repeat(1e6/totals,np.diff(normalized.indptr))
                        normalized.data=np.log1p(normalized.data)
                        if not np.isfinite(normalized.data).all():raise ValueError('Invalid normalized values')
                        out_data[left:right]=normalized.data;out_indices[left:right]=indices
                        if start%(cfg['block_cells']*10)==0:
                            self.progress.note('input.normalize.block','COMPLETE',dict(start=start,stop=stop),dict(nnz_written=right))
                event.update(shape=list(shape),nnz=nnz,samples=obs.technical_sample_id.nunique())
        temporary.rename(output)
        with self.progress.track('input.save_identity',output=str(output)) as event:
            identity=dict(source=str(source),source_sha256=observed,normalized_sha256=sha256(output),
                          shape=list(shape),normalization='ln(1+counts/cell_total*1e6)',
                          raw_count_copy=False,source_order_preserved=True)
            write_json(self.root/'inputs/input_identity.json',identity)
            event.update(normalized_sha256=identity['normalized_sha256'])

    def load_expression(self):
        """Load the shared derived asset; strip all sample/biological metadata."""
        path=self.root/'inputs/full_log1p_cpm.h5ad'
        expected=json.loads((self.root/'inputs/input_identity.json').read_text())['normalized_sha256']
        with self.progress.track('input.load_normalized',source=str(path)) as event:
            if sha256(path)!=expected:raise ValueError('Normalized input bytes changed')
            obj=ad.read_h5ad(path)
            obj.obs=pd.DataFrame(index=obj.obs_names)
            obj.uns.clear()
            event.update(shape=list(obj.shape),dtype=str(obj.X.dtype),nnz=int(obj.X.nnz))
        return obj

    def fit(self, fraction):
        """Fit actual cells, then save labels and held-out Pearson assignments."""
        cfg=self.cfg
        if fraction==1:
            mapped=np.load(self.root/'fits/fit80/all_cell_labels.npy')
            counts=np.unique(mapped,return_counts=True)[1]
            graph_size=sum(int(x)**2 for x in counts)
            write_json(self.root/'fits/wrapper_branch.json',dict(first_iteration_graph_proxy=graph_size,
                branch='large_full_fit' if graph_size>=10**9 else 'iter_consensus_clust'))
            if graph_size<10**9:
                raise RuntimeError('Measured smaller-graph branch requires separate qualification before benchmark continuation')
        obj=self.load_expression()
        directory=self.root/'fits'/self.stage;directory.mkdir()
        n=obj.n_obs
        if fraction==1:
            selected=np.arange(n)
        else:
            selected=np.sort(np.random.default_rng(cfg['sampling_seed']).choice(n,round(n*fraction),replace=False))
        np.save(directory/'sampled_source_rows.npy',selected)
        pd.DataFrame({'cell_id':obj.obs_names[selected],'source_row':selected}).to_csv(
            directory/'sampled_cells.tsv.gz',sep='\t',index=False)
        candidate=dict(name=self.stage,seed=cfg['fit_seed']+(fraction==1),score_thresh=cfg['thresholds']['score_thresh'])
        engine=FullScaleEngine(cfg,candidate,directory/'engine',self.progress)
        with self.progress.track('fit.complete_iterative',cells=len(selected),fraction=fraction,
                                 sampling_seed=cfg['sampling_seed'],fit_seed=candidate['seed']) as event:
            working=obj if fraction==1 else obj[selected,:].copy()
            labels,summary,evidence,means=engine.run(working)
            keys,numeric=np.unique(labels,return_inverse=True)
            numeric=numeric.astype(np.int32)+1
            np.save(directory/'sampled_labels.npy',numeric)
            pd.DataFrame({'cell_id':working.obs_names,'cluster':numeric,'display_label':labels}).to_csv(
                directory/'sampled_assignments.tsv.gz',sep='\t',index=False)
            event.update(summary=summary)
            del working,means,evidence
            gc.collect()
        markers=sorted(engine.markers)
        columns=obj.var_names.get_indexer(markers)
        if np.any(columns<0):raise ValueError('Unknown fit marker')
        np.save(directory/'marker_columns.npy',columns)
        pd.DataFrame({'gene_id':markers,'gene_symbol':obj.var.gene_symbol.iloc[columns].to_numpy()}).to_csv(
            directory/'mapping_markers.tsv',sep='\t',index=False)
        with self.progress.track('fit.heldout_assignment',heldout_cells=n-len(selected),markers=len(markers)) as event:
            if fraction==1:
                all_labels=numeric;scores=np.full(n,np.nan);prototype=np.empty((0,0));prototype_labels=np.array([])
            else:
                all_labels,scores,prototype,prototype_labels=map_heldout(obj.X,selected,numeric,columns,cfg['block_cells'])
            np.save(directory/'all_cell_labels.npy',all_labels)
            np.savez_compressed(directory/'mapping_model.npz',prototype_means=prototype,cluster_ids=prototype_labels,
                                marker_columns=columns,source_shape=np.asarray(obj.shape))
            sampled=np.zeros(n,dtype=bool);sampled[selected]=True
            pd.DataFrame(dict(cell_id=obj.obs_names,cluster=all_labels,sampled_for_fit=sampled,
                              heldout_correlation=scores)).to_csv(directory/'all_cell_assignments.tsv.gz',sep='\t',index=False)
            sizes=pd.Series(all_labels).value_counts().sort_index()
            sizes.to_csv(directory/'all_cell_cluster_sizes.tsv',sep='\t',header=['cells'],index_label='cluster')
            event.update(clusters=len(sizes),first_iteration_graph_proxy=sum(int(v)**2 for v in sizes))
        write_json(directory/'fit_identity.json',dict(fraction=fraction,cells_fit=len(selected),population=n,
            seed=candidate['seed'],sampling_seed=cfg['sampling_seed'],summary=summary,
            all_labels_sha256=sha256(directory/'all_cell_labels.npy'),marker_policy='first20_up_down_each_persisted_fit_DE_call',
            annotation_fields_used=False,status='BENCHMARK_NOT_FINAL_TAXONOMY'))

    def consensus(self):
        """Time the measured Allen large-graph branch and final Python DE merge.

        R=1 here is a computation benchmark only. The separate 100-column
        repeated-label exercise measures storage/one affinity pass; it is
        explicitly synthetic and cannot support stability or a final K.
        """
        cfg=self.cfg
        original=np.load(self.root/'fits/fit80/all_cell_labels.npy')
        initial=np.load(self.root/'fits/fullfit/all_cell_labels.npy')
        output=self.root/'outputs';output.mkdir()
        with self.progress.track('consensus.compile_membership',real_iterations=1,population=len(original)) as event:
            membership=Membership([original],cfg['block_cells'])
            sparse.save_npz(output/'real_R1_membership_B.npz',membership.b)
            write_json(output/'real_R1_membership_index.json',membership.block_index)
            event.update(shape=list(membership.b.shape),nnz=int(membership.b.nnz))
        with self.progress.track('consensus.merge_by_co',diff_th=.25) as event:
            merged,actions=membership.merge_by_co(initial)
            write_json(output/'consensus_merge_actions.json',actions)
            np.save(output/'labels_before_refinement.npy',merged)
            event.update(clusters=len(np.unique(merged)),eligible_pairs=len(actions))
        with self.progress.track('consensus.refine',tol_th=.01,confusion_th=.6,min_cells=20) as event:
            refined,trace,cells,medians=membership.refine(merged,min_cells=cfg['thresholds']['cluster_size_thresh'])
            np.save(output/'labels_after_refinement.npy',refined)
            write_json(output/'refinement_trace.json',trace)
            cells.to_csv(output/'R1_cell_consensus_diagnostics.tsv.gz',sep='\t',index_label='source_row')
            medians.to_csv(output/'R1_cluster_consensus_diagnostics.tsv',sep='\t',index_label='cluster')
            event.update(clusters=len(np.unique(refined)),trace_events=len(trace))
        # A capacity experiment with the correct matrix dimensions/nnz, not
        # 100 successful biological fits. Correlation structure is artificial.
        with self.progress.track('capacity.synthetic_R100_compile',synthetic=True,real_fits=1) as event:
            capacity=Membership([original]*100,cfg['block_cells'])
            event.update(shape=list(capacity.b.shape),nnz=int(capacity.b.nnz),
                         array_bytes=capacity.b.data.nbytes+capacity.b.indices.nbytes+capacity.b.indptr.nbytes)
        with self.progress.track('capacity.synthetic_R100_affinity',synthetic=True,passes=1) as event:
            a,unique=capacity.affinity(initial)
            expected,_=membership.affinity(initial)
            if not np.allclose(a,expected,atol=1e-12,rtol=1e-12):
                raise ValueError('R100 repeated-label capacity fixture changed exact affinity')
            event.update(query_cells=a.shape[0],target_clusters=a.shape[1],output_bytes=a.nbytes)
            del a,expected,capacity
            gc.collect()
        obj=self.load_expression()
        engine=FullScaleEngine(cfg,dict(name='post_consensus_DE',seed=cfg['fit_seed'],score_thresh=150),
                               output/'final_DE_engine',self.progress)
        columns=np.load(self.root/'fits/fullfit/marker_columns.npy')
        groups={int(k):np.flatnonzero(refined==k).tolist() for k in np.unique(refined)}
        with self.progress.track('consensus.final_DE_merge',representation='fit-marker expression, no PCA',
                                 thresholds=cfg['thresholds'],input_clusters=len(groups)) as event:
            # The R wrapper uses marker expression directly to prioritize pairs.
            # Keep the existing Python eBayes/two-neighbor merge implementation.
            projected=obj[:,columns].copy()
            with allocation_adapters():
                groups=engine._merge(obj,projected,groups,refined,engine.directory)
                labels=engine.labels(list(groups.values()),obj.n_obs)
                groups={key:np.flatnonzero(labels==key) for key in sorted(set(labels))}
                means,present,variances=tc.get_cluster_means(obj,groups,labels,low_th=cfg['low_thresh'])
            del projected
            event.update(output_clusters=len(groups))
        thresholds=dict(cfg['thresholds'])
        for key in ['score_thresh','low_thresh','min_genes']:thresholds.pop(key)
        with self.progress.track('consensus.final_DE_audit',pairs=len(groups)*(len(groups)-1)//2) as event:
            pairs=list(combinations(groups,2))
            evidence=tc.de_pairs_ebayes(pairs,means,variances,present,{k:len(v) for k,v in groups.items()},thresholds) if pairs else pd.DataFrame()
            if len(evidence):
                evidence['meets_separation']=(evidence.score>=150)&(evidence.num>=cfg['thresholds']['min_genes'])
                for col in ['up_genes','down_genes']:evidence[col]=evidence[col].map(lambda x:json.dumps(list(x)))
            unresolved=int((~evidence.meets_separation).sum()) if len(evidence) else 0
            evidence.to_csv(output/'R1_final_pairwise_DE.tsv',sep='\t')
            event.update(unresolved_pairs=unresolved,accepted_final_K=False)
        index=pd.read_csv(self.root/'inputs/full_cell_index.tsv.gz',sep='\t')
        index['benchmark_R1_cluster']=labels
        index.to_csv(output/'benchmark_R1_all_cell_assignments.tsv.gz',sep='\t',index=False)
        means.to_csv(output/'R1_cluster_mean_log1p_cpm.tsv.gz',sep='\t')
        present.to_csv(output/'R1_cluster_detection.tsv.gz',sep='\t')
        write_json(output/'benchmark_summary.json',dict(population=len(labels),real_80percent_fits=1,
            real_consensus_iterations=1,full_population_fits=1,diagnostic_cluster_count=len(groups),
            residual_DE_pairs=unresolved,final_K_accepted=False,annotations_locked=False,
            synthetic_R100_capacity_exercise=True,production_iterations_launched=0,status='IN_REVIEW'))
        self.report()

    def report(self):
        """Save measured stage costs and cautious 100-iteration projections."""
        output=self.root/'outputs'
        events=[]
        for path in (self.root/'progress').glob('*/hicat_progress_events.jsonl'):
            for line in path.read_text().splitlines():
                row=json.loads(line)
                if row['status']=='COMPLETE':events.append(row)
        frame=pd.DataFrame([dict(stage=e['stage'],seconds=e['operation_seconds'],
                     peak_rss_gib=e['peak_rss_gib'],inputs=json.dumps(e['inputs']),outputs=json.dumps(e['outputs'])) for e in events])
        frame.to_csv(output/'measured_operation_resources.tsv',sep='\t',index=False)
        self.make_plots(frame)
        write_json(output/'review_state.json',dict(status='IN_REVIEW',production_launch_authorized=False,
            approval_required='User requested stop after measured full-data benchmark report',
            report='measured_operation_resources.tsv + scheduler accounting; CPU projections computed after all stages finish',
            full_pairwise_object_created=False,normalization_asset='inputs/full_log1p_cpm.h5ad',
            raw_counts_in_normalized_asset=False,consensus_h5ad_saved=False,all_assignments_saved_as_tsv=True))
        manifest(output).to_csv(output/'output_manifest.tsv',sep='\t',index=False)

    def resource_report(self):
        """Finalize measured resource accounting after analysis processes exit.

        Observed CPU seconds are user+system process time. Allocated core-hours
        are wall hours times requested CPUs and can exceed observed CPU time.
        The R100 synthetic pass is a capacity timing, not a convergence model.
        """
        output=self.root/'outputs'
        rows=[json.loads((self.root/'progress'/stage/'stage_resources.json').read_text())
              for stage in ['prepare','fit80','fullfit','consensus']]
        if any(row['status']!='COMPLETED' for row in rows):
            raise RuntimeError('Cannot report an incomplete benchmark as completed')
        frame=pd.DataFrame(rows).set_index('stage')
        cpus=dict(prepare=1,fit80=8,fullfit=8,consensus=8)
        frame['requested_cpus']=[cpus[k] for k in frame.index]
        frame['allocated_core_hours']=frame.wall_seconds*frame.requested_cpus/3600
        frame['observed_CPU_hours']=frame.cpu_seconds/3600
        frame.to_csv(output/'measured_stage_resources.tsv',sep='\t')
        disk={}
        for key in ['inputs','fits/fit80','fits/fullfit','outputs','code','reference']:
            disk[key]=sum(p.stat().st_size for p in (self.root/key).rglob('*') if p.is_file())
        projected_fit_cpu=100*float(frame.loc['fit80','cpu_seconds'])/3600
        projected_fit_core=100*float(frame.loc['fit80','allocated_core_hours'])
        projected_disk=disk['inputs']+100*disk['fits/fit80']+disk['fits/fullfit']+disk['outputs']
        observed_peak=float(frame.peak_rss_gib.max())
        recommended_mem=int(np.ceil(max(32,observed_peak*1.5)/8)*8)
        recommended_hours=int(np.ceil(max(frame.loc['fit80','wall_seconds'],frame.loc['fullfit','wall_seconds'])*2/3600))+1
        estimate=dict(status='IN_REVIEW_STOP_FOR_USER_APPROVAL',real_subsample_size=357079,
            full_population_size=446349,real_iterations=1,production_iterations=100,
            production_jobs_submitted=0,measured_stage_resources=rows,measured_disk_bytes=disk,
            projected_100_fit_observed_CPU_hours=projected_fit_cpu,
            projected_100_fit_allocated_core_hours=projected_fit_core,
            measured_fullfit_observed_CPU_hours=float(frame.loc['fullfit','observed_CPU_hours']),
            measured_consensus_R1_observed_CPU_hours=float(frame.loc['consensus','observed_CPU_hours']),
            GPU_hours=0,projected_primary_disk_bytes=projected_disk,
            projected_disk_with_50percent_margin_bytes=projected_disk*1.5,
            recommended_cpus_per_fit=8,recommended_mem_GiB=recommended_mem,
            recommended_wall_hours_per_fit=recommended_hours,recommended_initial_fit_concurrency=2,
            serial_100_fit_wall_hours=100*float(frame.loc['fit80','wall_seconds'])/3600,
            ideal_2_concurrent_fit_wall_hours=50*float(frame.loc['fit80','wall_seconds'])/3600,
            uncertainty='One seed only; tree depth, K, graph fill-in and consensus refinement count can vary. '
                        'CPU/core-hour values for 100 fits are projections. R1 consensus time is not a measured R100 convergence cost. '
                        'Synthetic R100 capacity pass has repeated partitions and optimistic compression/structure. '
                        'Queue time, filesystem contention, safety margin and consensus convergence are additional.',
            consensus_storage='Sparse all-cell B, exact N-by-K affinities; no full N-by-N matrix',
            method='audited Allen large-graph control flow with disclosed existing Python fitting/DE backend')
        write_json(output/'RESOURCE_REPORT.json',estimate)
        text=('# Full-data HiCAT benchmark resource report\n\n'
              '**Completed benchmark only; stop for user approval before production.**\n\n'
              'Population: 446,349 dissected E14.5 mouse MGE cells. One independent 357,079-cell fit, '
              'held-out assignment and a full-population fit followed by the measured Allen large-graph branch. '
              'q1=.4, qdiff=.7, DEscore=150. No batch correction or cell-cycle regression.\n\n'
              'Measured stage resources:\n\n'+frame.to_string()+'\n\n'
              'Projected 100-fit CPU hours: %.2f; allocated core-hours: %.2f; GPU hours: 0. '
              'Measured full-data initialization CPU hours: %.2f.\n\n'%
              (projected_fit_cpu,projected_fit_core,estimate['measured_fullfit_observed_CPU_hours'])+
              'Primary disk projection: %.2f GiB; with 50%% margin: %.2f GiB.\n\n'%
              (projected_disk/2**30,projected_disk*1.5/2**30)+
              'Recommended starting request: %d CPUs, %d GiB, %d hours per fit, at most 2 concurrent fits.\n\n'%
              (8,recommended_mem,recommended_hours)+estimate['uncertainty']+'\n\n'
              'See `RESOURCE_REPORT.json`, `measured_operation_resources.tsv`, `benchmark_review.pdf`, '
              'the frozen audit and source manifests. R=1 labels/plots are diagnostic, not final K. '
              'The 100-column repeated-label capacity exercise is synthetic and is not 100 biological iterations.\n')
        (output/'RESOURCE_REPORT.md').write_text(text)
        manifest(output).to_csv(output/'output_manifest.tsv',sep='\t',index=False)
        (self.root/'BENCHMARK_COMPUTATION_SUCCESS.txt').write_text('IN_REVIEW; no production jobs submitted\n')

    def make_plots(self, resources):
        """Make resource bars and diagnostic cluster/sample heatmap and dendrogram."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        from scipy.cluster.hierarchy import linkage, dendrogram
        output=self.root/'outputs'
        with PdfPages(output/'benchmark_review.pdf') as pdf:
            selected=resources[resources.seconds>1].sort_values('seconds',ascending=False).head(20)
            fig,axes=plt.subplots(1,2,figsize=(13,8))
            axes[0].barh(selected.stage,selected.seconds/60);axes[0].set_xlabel('Measured operation minutes (nested timings overlap)')
            axes[1].barh(selected.stage,selected.peak_rss_gib);axes[1].set_xlabel('Lifetime process peak RSS, GiB')
            fig.suptitle('Real full-data benchmark — one 80% fit; no final K');fig.tight_layout();pdf.savefig(fig);plt.close(fig)
            assigned=pd.read_csv(output/'benchmark_R1_all_cell_assignments.tsv.gz',sep='\t')
            counts=pd.crosstab(assigned.benchmark_R1_cluster,assigned.technical_sample_id)
            fig,axes=plt.subplots(1,2,figsize=(13,max(5,len(counts)*.25)))
            axes[0].barh(counts.index,counts.sum(axis=1));axes[0].set_xlabel('Cells in R=1 diagnostic endpoint')
            axes[1].imshow(counts.div(counts.sum(axis=1),axis=0),aspect='auto',vmin=0,vmax=1)
            axes[1].set_xticks(np.arange(counts.shape[1]));axes[1].set_xticklabels(counts.columns,rotation=90)
            axes[1].set_yticks(np.arange(len(counts)));axes[1].set_yticklabels(counts.index)
            axes[1].set_title('Within-cluster sample fractions; all cells')
            fig.tight_layout();pdf.savefig(fig);plt.close(fig)
            means=pd.read_csv(output/'R1_cluster_mean_log1p_cpm.tsv.gz',sep='\t',index_col=0)
            if len(means)>1:
                # Display-only dendrogram: not the iterative fitting tree.
                z=linkage(means.to_numpy(),method='average',metric='correlation')
                np.save(output/'R1_expression_display_linkage.npy',z)
                fig,ax=plt.subplots(figsize=(12,5));dendrogram(z,labels=means.index.to_list(),ax=ax)
                ax.set_title('Expression-mean display dendrogram — R=1 diagnostic, not final consensus')
                fig.tight_layout();pdf.savefig(fig);plt.close(fig)

    def run(self):
        """Record success/failure and per-process measured CPU, wall time and RSS."""
        self.verify_frozen()
        status='FAILED'
        try:
            if self.stage=='prepare':self.prepare()
            elif self.stage=='fit80':self.fit(.8)
            elif self.stage=='fullfit':self.fit(1.)
            elif self.stage=='consensus':self.consensus()
            elif self.stage=='report':self.resource_report()
            else:raise ValueError(self.stage)
            status='COMPLETED'
        finally:
            usage=resource.getrusage(resource.RUSAGE_SELF)
            write_json(self.root/'progress'/self.stage/'stage_resources.json',dict(stage=self.stage,status=status,
                wall_seconds=time.monotonic()-self.start,cpu_seconds=usage.ru_utime+usage.ru_stime,
                peak_rss_gib=usage.ru_maxrss/1024**2,production=False))


def main():
    """Run a frozen benchmark stage; no production/array submission interface."""
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-dir',required=True,type=Path)
    p.add_argument('--stage',choices=['prepare','fit80','fullfit','consensus','report'],required=True)
    args=p.parse_args();Benchmark(args.run_dir,args.stage).run()


if __name__=='__main__':main()
