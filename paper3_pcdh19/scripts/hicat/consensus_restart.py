"""Restart-safe independent HiCAT iterations and real-output aggregation test.

The benchmark permits iteration IDs 0 and 1 only. Each fit and mapping has an
independent immutable checkpoint; final aggregation and final DE have their
own checkpoints. A completion seal is written only after reopening/hash checks.
No function in this module launches 100 production jobs or chooses final K.
"""
import argparse
import gc
from itertools import combinations
import json
import os
from pathlib import Path
import shutil
import resource
import time
import subprocess
import sys
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
import transcriptomic_clustering as tc
from .checkpoints import checkpoint,complete_attempt
from .consensus_core import Membership,map_heldout
from .fullscale_engine import FullScaleEngine,allocation_adapters
from .provenance import Progress,write_json,sha256


def link_copy(source,destination):
    """Preserve a completed input using hard links, falling back to byte copies."""
    source=Path(source);destination=Path(destination)
    if source.is_dir():
        shutil.copytree(source,destination,copy_function=link_copy)
    else:
        try:os.link(source,destination)
        except OSError:shutil.copy2(source,destination)
    return str(destination)


def validate_fit(directory,cell_ids,rows,genes):
    """Check saved partitions and each saved node's cell/model/graph alignment."""
    directory=Path(directory);labels=np.load(directory/'sampled_labels.npy')
    selected=np.load(directory/'sampled_source_rows.npy')
    if not np.array_equal(rows,selected) or len(labels)!=len(rows) or np.any(labels<1):
        raise ValueError('Fit labels do not match frozen sampled rows')
    saved=pd.read_csv(directory/'sampled_assignments.tsv.gz',sep='\t')
    if not np.array_equal(saved.cell_id,cell_ids[rows]) or not np.array_equal(saved.cluster,labels):
        raise ValueError('Fit assignment IDs/labels changed on disk')
    nodes=0
    for membership in (directory/'engine').rglob('cell_membership.tsv'):
        node=membership.parent;table=pd.read_csv(membership,sep='\t')
        local=table.pilot_row.to_numpy(dtype=int)
        if np.any(local<0) or np.any(local>=len(rows)) or not np.array_equal(table.local_row,np.arange(len(local))):
            raise ValueError('Invalid node source positions')
        if not np.array_equal(table.cell_id,cell_ids[rows[local]]):raise ValueError('Node IDs do not match sampled universe')
        if (node/'projection.npz').exists():
            with np.load(node/'projection.npz',allow_pickle=False) as archive:
                if not np.array_equal(archive['cell_ids'],table.cell_id):raise ValueError('Projection ID mismatch')
                if archive['coordinates'].shape[0]!=len(table) or not np.isfinite(archive['coordinates']).all():
                    raise ValueError('Projection shape/numeric mismatch')
                dimensions=archive['coordinates'].shape[1]
            model=pd.read_csv(node/'pca_components_retained.tsv',sep='\t',index_col=0)
            if model.shape[1]!=dimensions or not model.index.isin(genes).all():raise ValueError('Saved PCA model mismatch')
        if (node/'graph.npz').exists():
            graph=sparse.load_npz(node/'graph.npz');graph.check_format(full_check=True)
            if graph.shape!=(len(table),len(table)):raise ValueError('Graph shape mismatch')
            assignment=pd.read_csv(node/'graph_membership.tsv',sep='\t')
            if not np.array_equal(assignment.cell_id,table.cell_id):raise ValueError('Graph assignment IDs mismatch')
        nodes+=1
    return dict(sampled_cells=len(rows),population=len(cell_ids),genes=len(genes),nodes_checked=nodes,
                saved_partition_exact=True,saved_ID_alignment=True,all_artifacts_reopened_by_completion_seal=True)


class RestartRun:
    """Execute one independently resumable stage under a frozen global contract."""

    def __init__(self,root,production=False):
        """Verify frozen code/controls; load the cell universe without expression."""
        self.root=Path(root)
        self.cfg=json.loads((self.root/'config/run.json').read_text())
        self.production=bool(production)
        if self.production and not self.production_authorized():
            raise ValueError('Production approval has not been recorded')
        for relative,digest in json.loads((self.root/'FROZEN.json').read_text()).items():
            if sha256(self.root/relative)!=digest:raise ValueError('Frozen source changed: '+relative)
        installed=subprocess.check_output([sys.executable,'-m','pip','freeze'],text=True)
        if installed!=(self.root/'config/environment.txt').read_text():
            raise ValueError('Installed package versions differ from frozen environment')
        self.universe=pd.read_csv(self.root/'config/cell_universe.tsv.gz',sep='\t')
        self.ids=self.universe.cell_id.to_numpy(dtype=str)
        self.genes=np.load(self.root/'config/gene_ids.npy',allow_pickle=False)
        self.legacy=Path(self.cfg['legacy_benchmark'])
        self.identity=json.loads((self.legacy/'inputs/input_identity.json').read_text())

    def production_authorized(self):
        """Read the separate operational permission; never alter frozen science."""
        path=self.root/'AUTHORIZATION.json'
        return path.exists() and json.loads(path.read_text()).get('production_allowed') is True

    @property
    def iteration_ids(self):
        """Separate a two-real-iteration I/O test from eventual all-100 aggregation."""
        return list(range(self.cfg['production_iterations'])) if self.production else self.cfg['benchmark_iteration_ids']

    @property
    def final_root(self):
        """Keep benchmark aggregation and eventual production results independent."""
        return self.root/'final'/('production' if self.production else 'benchmark')

    def contract(self,stage,index=None,parents=None):
        """Bind output reuse to exact seeds, code, input universe and parent seals."""
        return dict(stage=stage,iteration=index,parents=parents or {},
            run_config_sha256=sha256(self.root/'config/run.json'),
            adapter_commit=self.cfg['adapter_commit'],universe_sha256=self.cfg['universe_sha256'],
            normalized_sha256=self.identity['normalized_sha256'],
            iteration_config_sha256=sha256(self.root/'iterations'/('%03d'%index)/'config.json') if index is not None else None)

    def directory(self,index,stage):
        """Map one iteration/stage to its independent attempt/checkpoint directory."""
        if index is None or index<0 or index>=self.cfg['production_iterations']:
            raise ValueError('Invalid iteration ID')
        if index not in self.cfg['benchmark_iteration_ids'] and not self.production_authorized():
            raise ValueError('Production approval required for additional iteration IDs')
        return self.root/'iterations'/('%03d'%index)/stage

    def load(self):
        """Read the shared normalized object; fitting sees IDs and gene metadata only."""
        path=self.legacy/'inputs/full_log1p_cpm.h5ad'
        if sha256(path)!=self.identity['normalized_sha256']:raise ValueError('Normalized input changed')
        obj=ad.read_h5ad(path)
        if not np.array_equal(obj.obs_names,self.ids) or not np.array_equal(obj.var_names,self.genes):
            raise ValueError('Shared normalized object has changed row/gene ordering')
        obj.obs=pd.DataFrame(index=obj.obs_names);obj.uns.clear()
        return obj

    def fit(self,index):
        """Fit/reuse one sampled partition independently of held-out assignment."""
        iteration=self.root/'iterations'/('%03d'%index)
        cfg=json.loads((iteration/'config.json').read_text())
        rows=np.load(iteration/'sampled_source_rows.npy')
        def compute(attempt):
            """Fit only once, or adopt the already-running first benchmark fit."""
            if index==0:
                old=self.legacy/'fits/fit80'
                status=json.loads((self.legacy/'progress/fit80/stage_resources.json').read_text())
                if status['status']!='COMPLETED':raise ValueError('First benchmark fit not complete')
                if not np.array_equal(np.load(old/'sampled_source_rows.npy'),rows):raise ValueError('Legacy sample differs from frozen seed')
                for name in ['engine','sampled_source_rows.npy','sampled_labels.npy','sampled_assignments.tsv.gz','marker_columns.npy']:
                    link_copy(old/name,attempt/name)
                write_json(attempt/'adopted_resource_record.json',status)
                write_json(attempt/'adoption_provenance.json',dict(source=str(old),source_config_sha256=sha256(self.legacy/'config/benchmark.json'),
                    fit_code_sha256=sha256(self.legacy/'code/hicat/fullscale_engine.py'),adopted_without_refitting=True))
                return
            obj=self.load();selected=obj[rows,:].copy();del obj;gc.collect()
            progress=Progress(attempt/'progress')
            engine=FullScaleEngine(self.cfg['fitting'],dict(name='iteration_%03d'%index,seed=cfg['fit_seed'],score_thresh=150),
                                   attempt/'engine',progress)
            labels,summary,evidence,means=engine.run(selected)
            _,numeric=np.unique(labels,return_inverse=True);numeric=numeric.astype(np.int32)+1
            np.save(attempt/'sampled_source_rows.npy',rows);np.save(attempt/'sampled_labels.npy',numeric)
            np.save(attempt/'marker_columns.npy',selected.var_names.get_indexer(sorted(engine.markers)))
            pd.DataFrame(dict(cell_id=selected.obs_names,cluster=numeric)).to_csv(attempt/'sampled_assignments.tsv.gz',sep='\t',index=False)
            write_json(attempt/'fit_summary.json',summary)
        return checkpoint(self.directory(index,'fit'),self.contract('fit',index),compute,
                          lambda attempt:validate_fit(attempt,self.ids,rows,self.genes))

    def mapping(self,index,inject_failure=False):
        """Resume held-out mapping/B generation from a sealed fit without refitting."""
        fitted,seal=complete_attempt(self.directory(index,'fit'),self.contract('fit',index))
        labels=np.load(fitted/'sampled_labels.npy');rows=np.load(fitted/'sampled_source_rows.npy')
        columns=np.load(fitted/'marker_columns.npy')
        contract=self.contract('mapping',index,dict(fit=seal['manifest_sha256']))
        def compute(attempt):
            """Publish an exact iteration block plus fitted/inferred label evidence."""
            progress=Progress(attempt/'progress')
            with progress.track('mapping_and_membership',population=len(self.ids),sampled=len(rows),markers=len(columns)):
                if index==0:
                    old=self.legacy/'fits/fit80';all_labels=np.load(old/'all_cell_labels.npy')
                    table=pd.read_csv(old/'all_cell_assignments.tsv.gz',sep='\t')
                    if not np.array_equal(table.cell_id,self.ids):raise ValueError('Legacy mapping IDs differ')
                    scores=table.heldout_correlation.to_numpy()
                    link_copy(old/'mapping_model.npz',attempt/'mapping_model.npz')
                    write_json(attempt/'adoption_provenance.json',dict(source=str(old),heldout_assignment_reused=True))
                else:
                    obj=self.load()
                    all_labels,scores,prototypes,keys=map_heldout(obj.X,rows,labels,columns,self.cfg['fitting']['block_cells'])
                    np.savez_compressed(attempt/'mapping_model.npz',prototype_means=prototypes,cluster_ids=keys,marker_columns=columns)
                    del obj;gc.collect()
                sampled=np.zeros(len(self.ids),dtype=bool);sampled[rows]=True
                np.save(attempt/'all_cell_labels.npy',all_labels)
                np.save(attempt/'sampled_mask.npy',sampled)
                np.save(attempt/'heldout_correlation.npy',scores)
                pd.DataFrame(dict(cell_id=self.ids,cluster=all_labels,sampled_for_fit=sampled,heldout_correlation=scores)).to_csv(
                    attempt/'all_cell_assignments.tsv.gz',sep='\t',index=False)
                membership=Membership([all_labels])
                sparse.save_npz(attempt/'membership_B.npz',membership.b)
                write_json(attempt/'membership_index.json',dict(universe_sha256=self.cfg['universe_sha256'],
                    population=len(self.ids),iteration=index,blocks=membership.block_index,denominator_contribution=1,
                    heldout_cells_included=True,source_fit_manifest=seal['manifest_sha256']))
                if inject_failure:
                    raise RuntimeError('INTENTIONAL_RESTART_TEST: membership saved; completion not published')
        def validate(attempt):
            """Reopen assignments/B and compare every cell against the frozen fit."""
            all_labels=np.load(attempt/'all_cell_labels.npy');sampled=np.load(attempt/'sampled_mask.npy')
            if all_labels.shape!=(len(self.ids),) or np.any(all_labels<1):raise ValueError('Incomplete all-cell labels')
            if not np.array_equal(np.flatnonzero(sampled),rows) or not np.array_equal(all_labels[rows],labels):
                raise ValueError('Held-out mapping changed sampled labels')
            b=sparse.load_npz(attempt/'membership_B.npz');b.check_format(full_check=True)
            Membership.from_blocks([b])
            unique=np.unique(all_labels)
            if not np.array_equal(unique[b.indices],all_labels):raise ValueError('Membership columns disagree with all-cell labels')
            table=pd.read_csv(attempt/'all_cell_assignments.tsv.gz',sep='\t')
            if not np.array_equal(table.cell_id,self.ids) or not np.array_equal(table.cluster,all_labels):raise ValueError('Saved mapping ID mismatch')
            return dict(population=len(self.ids),sampled=len(rows),heldout=len(self.ids)-len(rows),B_shape=list(b.shape),B_nnz=int(b.nnz),
                        all_cells_assigned_once=True,sampled_labels_preserved=True,full_pair_object_created=False)
        return checkpoint(self.directory(index,'mapping'),contract,compute,validate)

    def aggregate(self):
        """Consume two REAL completed blocks, without forming a full pair matrix."""
        inputs=[];parents={}
        for index in self.iteration_ids:
            _,fit_seal=complete_attempt(self.directory(index,'fit'),self.contract('fit',index))
            expected=self.contract('mapping',index,dict(fit=fit_seal['manifest_sha256']))
            path,seal=complete_attempt(self.directory(index,'mapping'),expected)
            info=json.loads((path/'membership_index.json').read_text())
            if info['universe_sha256']!=self.cfg['universe_sha256']:raise ValueError('Iteration cell universes differ')
            inputs.append(path);parents[str(index)]=seal['manifest_sha256']
        old=self.legacy/'fits/fullfit'
        status=json.loads((self.legacy/'progress/fullfit/stage_resources.json').read_text())
        if status['status']!='COMPLETED':raise ValueError('Allen full-data initializer not complete')
        parents['fullfit_labels']=sha256(old/'all_cell_labels.npy')
        def compute(attempt):
            """Perform the audited large-graph merge/refine branch using saved B."""
            progress=Progress(attempt/'progress')
            initial=np.load(old/'all_cell_labels.npy')
            first=np.load(inputs[0]/'all_cell_labels.npy');counts=np.unique(first,return_counts=True)[1]
            proxy=sum(int(v)**2 for v in counts)
            if proxy<10**9:raise RuntimeError('Smaller Allen branch needs qualification; no silent substitute')
            with progress.track('aggregate.real_saved_memberships',real_iterations=len(inputs)) as event:
                blocks=[sparse.load_npz(path/'membership_B.npz') for path in inputs]
                membership=Membership.from_blocks(blocks,self.cfg['fitting']['block_cells'])
                sparse.save_npz(attempt/'membership_B.npz',membership.b)
                write_json(attempt/'membership_index.json',dict(real_iteration_inputs=[str(p) for p in inputs],
                    parent_seals=parents,blocks=membership.block_index,universe_sha256=self.cfg['universe_sha256'],denominator=len(inputs)))
                event.update(shape=list(membership.b.shape),nnz=int(membership.b.nnz),full_pair_matrix=False)
            with progress.track('aggregate.merge_by_co'):
                labels,actions=membership.merge_by_co(initial);write_json(attempt/'merge_actions.json',actions)
                np.save(attempt/'labels_before_refinement.npy',labels)
            with progress.track('aggregate.refine'):
                labels,trace,cells,medians=membership.refine(labels,min_cells=20)
                np.save(attempt/'labels_after_refinement.npy',labels)
                write_json(attempt/'refinement_trace.json',trace)
                cells.to_csv(attempt/'cell_consensus_diagnostics.tsv.gz',sep='\t',index_label='source_row')
                medians.to_csv(attempt/'cluster_consensus_diagnostics.tsv',sep='\t',index_label='cluster')
                pd.DataFrame(dict(cell_id=self.ids,diagnostic_cluster=labels)).to_csv(attempt/'all_cell_assignments.tsv.gz',sep='\t',index=False)
            write_json(attempt/'scope.json',dict(real_iterations=len(inputs),population=len(self.ids),
                       graph_proxy=proxy,interpret_K=False,production=False,full_pair_matrix_created=False))
        def validate(attempt):
            """Reopen actual combined B and compare to independently loaded blocks."""
            combined=sparse.load_npz(attempt/'membership_B.npz');combined.check_format(full_check=True)
            expected=sparse.vstack([sparse.load_npz(path/'membership_B.npz') for path in inputs],format='csc')
            if (combined!=expected).nnz or not np.all(np.asarray(combined.sum(axis=0)).ravel()==len(inputs)):
                raise ValueError('Combined membership differs from completed iteration artifacts')
            table=pd.read_csv(attempt/'all_cell_assignments.tsv.gz',sep='\t')
            if not np.array_equal(table.cell_id,self.ids) or table.diagnostic_cluster.isna().any():raise ValueError('Incomplete aggregation labels')
            return dict(real_iterations=len(inputs),population=len(self.ids),membership_nnz=int(combined.nnz),
                        saved_blocks_combined_exactly=True,full_pair_matrix_created=False,K_interpreted=False)
        return checkpoint(self.final_root/'aggregation',self.contract('aggregation',parents=parents),compute,validate)

    def merge(self):
        """Resume final DE only from sealed aggregation; never refit iterations."""
        aggregate,seal=complete_attempt(self.final_root/'aggregation')
        columns=np.load(self.legacy/'fits/fullfit/marker_columns.npy')
        contract=self.contract('final_DE',parents=dict(aggregation=seal['manifest_sha256'],
                               marker_columns=sha256(self.legacy/'fits/fullfit/marker_columns.npy')))
        def compute(attempt):
            """Apply documented marker-expression neighbor prioritization and DE."""
            obj=self.load();labels=np.load(aggregate/'labels_after_refinement.npy')
            progress=Progress(attempt/'progress')
            engine=FullScaleEngine(self.cfg['fitting'],dict(name='final_DE_R2',seed=20260908,score_thresh=150),attempt/'engine',progress)
            groups={int(k):np.flatnonzero(labels==k).tolist() for k in np.unique(labels)}
            with progress.track('final_DE.merge',input_groups=len(groups)):
                projected=obj[:,columns].copy()
                with allocation_adapters():
                    groups=engine._merge(obj,projected,groups,labels,engine.directory)
                    labels=engine.labels(list(groups.values()),obj.n_obs)
                    groups={key:np.flatnonzero(labels==key) for key in sorted(set(labels))}
                    means,present,variances=tc.get_cluster_means(obj,groups,labels,low_th=1)
                del projected;gc.collect()
            thresholds=dict(self.cfg['fitting']['thresholds'])
            for key in ['low_thresh','min_genes','score_thresh']:thresholds.pop(key)
            with progress.track('final_DE.audit',pairs=len(groups)*(len(groups)-1)//2):
                pairs=list(combinations(groups,2))
                evidence=tc.de_pairs_ebayes(pairs,means,variances,present,{key:len(v) for key,v in groups.items()},thresholds) if pairs else pd.DataFrame()
                if len(evidence):
                    evidence['meets_separation']=(evidence.score>=150)&(evidence.num>=5)
                    for col in ['up_genes','down_genes']:evidence[col]=evidence[col].map(lambda x:json.dumps(list(x)))
                evidence.to_csv(attempt/'pairwise_DE.tsv',sep='\t')
            means.to_csv(attempt/'cluster_mean_log1p_cpm.tsv.gz',sep='\t');present.to_csv(attempt/'cluster_detection.tsv.gz',sep='\t')
            pd.DataFrame(dict(cell_id=self.ids,diagnostic_cluster=labels)).to_csv(attempt/'all_cell_assignments.tsv.gz',sep='\t',index=False)
            write_json(attempt/'scope.json',dict(real_iterations=len(self.iteration_ids),final_K_accepted=False,K_interpreted=False,annotations_locked=False,
                         residual_DE_pairs=int((~evidence.meets_separation).sum()) if len(evidence) else 0,status='IN_REVIEW'))
        def validate(attempt):
            """Confirm the endpoint retains each source cell once without final-K claims."""
            table=pd.read_csv(attempt/'all_cell_assignments.tsv.gz',sep='\t')
            if not np.array_equal(table.cell_id,self.ids) or table.diagnostic_cluster.isna().any():raise ValueError('Incomplete final-DE cell partition')
            return dict(population=len(self.ids),source_ID_order_exact=True,K_interpreted=False)
        return checkpoint(self.final_root/'DE',contract,compute,validate)


def main():
    """Run one benchmark stage or inspect completed outputs; no production launch."""
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir',type=Path,required=True)
    parser.add_argument('--stage',choices=['iteration','fit','mapping','aggregate','merge'],required=True)
    parser.add_argument('--iteration',type=int)
    parser.add_argument('--inject-failure',action='store_true')
    parser.add_argument('--production',action='store_true')
    args=parser.parse_args();started=time.monotonic();status='FAILED'
    try:
        run=RestartRun(args.run_dir,args.production)
        if args.stage=='iteration':
            run.fit(args.iteration);run.mapping(args.iteration,args.inject_failure)
        elif args.stage=='fit':run.fit(args.iteration)
        elif args.stage=='mapping':run.mapping(args.iteration,args.inject_failure)
        elif args.stage=='aggregate':run.aggregate()
        elif args.stage=='merge':run.merge()
        status='COMPLETED'
    finally:
        folder=args.run_dir/'job_metrics';folder.mkdir(exist_ok=True)
        usage=resource.getrusage(resource.RUSAGE_SELF)
        job=os.environ.get('SLURM_JOB_ID','local_'+str(os.getpid()))
        write_json(folder/(job+'.json'),dict(stage=args.stage,iteration=args.iteration,status=status,
            injected_failure=args.inject_failure,production=args.production,slurm_job_id=job,
            wall_seconds=time.monotonic()-started,cpu_seconds=usage.ru_utime+usage.ru_stime,
            peak_rss_gib=usage.ru_maxrss/1024**2))


if __name__=='__main__':main()
