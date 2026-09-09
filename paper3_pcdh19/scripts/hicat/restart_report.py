"""Report the measured two-real-iteration benchmark and independent reuse proof.

This report performs no fitting and submits no production jobs. All completed
iteration/final checkpoints are revalidated before resources or reuse claims
are reported. K is deliberately not interpreted in this I/O benchmark.
"""
import argparse
import json
from pathlib import Path
import subprocess
import numpy as np
import pandas as pd
from .checkpoints import complete_attempt
from .consensus_restart import RestartRun
from .provenance import write_json,sha256


def main():
    """Revalidate two real iterations, replay completed stages, and write review."""
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--run-dir',type=Path,required=True)
    args=parser.parse_args();root=args.run_dir;run=RestartRun(root)
    report=root/'report';report.mkdir(exist_ok=True)
    rows=[];reused=[]
    for index in run.iteration_ids:
        for stage in ['fit','mapping']:
            path,complete=complete_attempt(run.directory(index,stage))
            before=sha256(path/'ARTIFACTS.json')
            replay=run.fit(index) if stage=='fit' else run.mapping(index)
            after=sha256(path/'ARTIFACTS.json')
            if replay!=path or before!=after:raise ValueError('Reuse changed completed output')
            resource=json.loads((path/'resources.json').read_text())
            actual_bytes=sum(p.stat().st_size for p in path.rglob('*') if p.is_file())
            rows.append(dict(iteration=index,stage=stage,**resource,disk_bytes=actual_bytes))
            reused.append(dict(iteration=index,stage=stage,attempt=str(path),manifest_before=before,manifest_after=after,
                               reused_without_recompute=True))
    # The final-stage replay also consumes sealed inputs and returns unchanged
    # checkpoints. A later failure can therefore rerun only final stages.
    for stage,method in [('aggregation',run.aggregate),('DE',run.merge)]:
        path,_=complete_attempt(run.final_root/stage);before=sha256(path/'ARTIFACTS.json')
        replay=method()
        if replay!=path or sha256(path/'ARTIFACTS.json')!=before:raise ValueError('Final-stage reuse changed output')
        reused.append(dict(stage=stage,attempt=str(path),reused_without_recompute=True,manifest_sha256=before))
    frame=pd.DataFrame(rows);frame.to_csv(report/'checkpoint_resources.tsv',sep='\t',index=False)
    jobs=[json.loads(p.read_text()) for p in (root/'job_metrics').glob('*.json')]
    jobframe=pd.DataFrame(jobs);jobframe.to_csv(report/'job_resources.tsv',sep='\t',index=False)
    native=jobframe[(jobframe.iteration==1)&jobframe.stage.isin(['fit','mapping'])&(jobframe.status=='COMPLETED')]
    if set(native.stage)!=set(['fit','mapping']):raise ValueError('Complete native production-like iteration is missing')
    # One row per successful phase; extra validation-only replays must not inflate
    # the per-iteration estimate. Normal retries are retained in the full table.
    native=native.sort_values('wall_seconds').groupby('stage',as_index=False).tail(1)
    wall=float(native.wall_seconds.sum());cpu=float(native.cpu_seconds.sum());peak=float(native.peak_rss_gib.max())
    failed=jobframe[(jobframe.iteration==1)&(jobframe.stage=='mapping')&(jobframe.status=='FAILED')&jobframe.injected_failure]
    if len(failed)<1:raise ValueError('Real mapping failure/resubmission test did not execute')
    failed_artifacts=list((run.directory(1,'mapping')/'attempts').glob('*/FAILURE.json'))
    if not failed_artifacts:raise ValueError('Failed mapping attempt was not preserved')
    normalized=Path(run.legacy/'inputs/full_log1p_cpm.h5ad').stat().st_size
    disks=frame.groupby('iteration').disk_bytes.sum()
    typical_disk=int(disks.max());fit_peak=max(float(jobframe.peak_rss_gib.max()),peak)
    memory=int(np.ceil(max(32,fit_peak*1.5)/8)*8)
    hours=int(np.ceil(wall*2/3600))+1
    projection=normalized+typical_disk*100+sum(p.stat().st_size for p in (run.legacy/'fits/fullfit').rglob('*') if p.is_file())
    total_final_bytes=sum(p.stat().st_size for p in run.final_root.rglob('*') if p.is_file())
    metadata_bytes=sum(p.stat().st_size for p in (root/'iterations').rglob('sampled*') if p.is_file() and 'attempts' not in str(p))
    summary=dict(status='IN_REVIEW_STOP_FOR_APPROVAL',real_iterations=2,production_iterations_executed=0,
        cells=446349,sampled_per_iteration=357079,normalization='ln(1+CPM)',q1=.4,qdiff=.7,DEscore=150,
        complete_native_iteration_wall_seconds=wall,complete_native_iteration_CPU_seconds=cpu,
        complete_native_iteration_peak_RSS_GiB=peak,measured_disk_per_iteration_bytes={str(k):int(v) for k,v in disks.items()},
        projected_100_iteration_CPU_hours=cpu*100/3600,projected_100_iteration_allocated_core_hours=wall*8*100/3600,
        GPU_hours=0,recommended_CPUs_per_iteration=8,recommended_RAM_GiB=memory,recommended_time_limit_hours=hours,
        maximum_recommended_initial_fit_concurrency=2,concurrency_is_conservative_not_cluster_wide_guarantee=True,
        projected_fit_wall_hours_at_concurrency_2=wall*50/3600,
        projected_storage_base_bytes=projection,projected_storage_with_50percent_headroom_bytes=(projection+metadata_bytes)*1.5,
        measured_two_iteration_final_stage_disk_bytes=total_final_bytes,shared_normalized_object_bytes=normalized,
        full_pair_matrix_created=False,consensus_representation='iteration CSC one-hot blocks; concatenated B; exact N-by-K affinities',
        all_completed_jobs_independently_reusable=True,real_mapping_failure_and_individual_retry_passed=True,
        final_stage_replay_reused_all_completed_iterations=True,
        adapter_commit=run.cfg['adapter_commit'],Allen_R_commit=run.cfg['fitting']['consensus_reference_commit'],
        Allen_Python_commit=run.cfg['fitting']['upstream_commit'],seed_schedule_sha256=sha256(root/'config/seeds.tsv'),
        all_100_seeds_and_sample_ID_sets_frozen=True,K_interpreted=False,
        caveats='Two seeds only. Independent 100-fit CPU/wall/disk are projections, not measurements. '
                'Add full-data initialization, final consensus/DE runtime, queue time and I/O contention. '
                'Two-iteration refinement cannot establish 100-iteration convergence or final K. '
                'Final-stage storage can grow with observed K and iterations. Initial concurrency recommendation is 2, '
                'subject to account/node limits and measured filesystem behavior.')
    write_json(report/'RESOURCE_REPORT.json',summary);write_json(report/'RESTART_VALIDATION.json',dict(reused=reused,
        failed_mapping_jobs=failed.slurm_job_id.to_list(),failure_attempts=[str(p) for p in failed_artifacts],no_completed_fit_recomputed=True))
    text=('# Completed real full-data benchmark and restart review\n\n'
          '**Stop here for approval. No 100-run production submission.**\n\n'
          'Two real 357,079-cell iterations from the full 446,349-cell E14.5 mouse MGE population were saved, reopened, '
          'validated and combined using the audited Allen large-graph consensus branch. K is not interpreted.\n\n'
          'Native complete iteration: %.2f hours wall, %.2f observed CPU-hours, %.2f GiB peak RSS. '
          'These include independent fit and mapping processes plus their checkpoint validation.\n\n'%
          (wall/3600,cpu/3600,peak)+
          'Disk per completed iteration (fit + mapping):\n\n'+disks.to_string()+'\n\n'
          'Projection for 100 iterations: %.2f observed CPU-hours; %.2f allocated core-hours; '
          '%.2f hours at concurrency 2 before initialization/consensus/queue overhead; GPU-hours 0.\n\n'%
          (cpu*100/3600,wall*8*100/3600,wall*50/3600)+
          'Recommended initial per-iteration request: 8 CPUs, %d GiB RAM, %d-hour time limit; at most two simultaneous fits.\n\n'%
          (memory,hours)+
          'Projected primary storage: %.2f GiB; with 50%% headroom: %.2f GiB.\n\n'%
          (projection/2**30,(projection+metadata_bytes)*1.5/2**30)+
          'A real mapping job was deliberately interrupted after saving membership and individually retried. '
          'Its fit checkpoint was reused unchanged. Replaying all successful fit/mapping/final stages returned '
          'the same artifact seals. Failed attempts and logs are retained. All 100 seed/sample-ID sets, numerical '
          'settings, environment and adapter/source commits are frozen.\n\n'+summary['caveats']+'\n')
    (report/'RESOURCE_REPORT.md').write_text(text)
    # Scientific displays are explicitly limited to the two-iteration I/O test.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from scipy.cluster.hierarchy import linkage,dendrogram
    final,_=complete_attempt(run.final_root/'DE')
    with PdfPages(report/'benchmark_review.pdf') as pdf:
        fig,axes=plt.subplots(1,2,figsize=(12,5));successful=jobframe[jobframe.status=='COMPLETED']
        names=successful.stage+' / '+successful.iteration.astype(str)
        axes[0].barh(names,successful.wall_seconds/60);axes[0].set_xlabel('Measured job minutes')
        axes[1].barh(names,successful.peak_rss_gib);axes[1].set_xlabel('Measured process peak GiB')
        fig.suptitle('Two real iterations: resource benchmark, no K interpretation');fig.tight_layout();pdf.savefig(fig);plt.close(fig)
        assigned=pd.read_csv(final/'all_cell_assignments.tsv.gz',sep='\t')
        counts=pd.crosstab(assigned.diagnostic_cluster,run.universe.technical_sample_id)
        counts.to_csv(report/'diagnostic_sample_composition.tsv',sep='\t')
        fig,ax=plt.subplots(figsize=(12,max(5,len(counts)*.2)))
        display=ax.imshow(counts.div(counts.sum(axis=1),axis=0),aspect='auto',vmin=0,vmax=1)
        ax.set_xticks(np.arange(counts.shape[1]));ax.set_xticklabels(counts.columns,rotation=90)
        ax.set_yticks(np.arange(len(counts)));ax.set_yticklabels(counts.index)
        fig.colorbar(display,ax=ax,label='Within-cluster sample fraction')
        ax.set_title('R=2 diagnostic endpoint; all cells; not a final taxonomy');fig.tight_layout();pdf.savefig(fig);plt.close(fig)
        means=pd.read_csv(final/'cluster_mean_log1p_cpm.tsv.gz',sep='\t',index_col=0)
        if len(means)>1:
            z=linkage(means.to_numpy(),method='average',metric='correlation')
            np.save(report/'diagnostic_expression_linkage.npy',z)
            fig,ax=plt.subplots(figsize=(12,5));dendrogram(z,labels=means.index.to_list(),ax=ax)
            ax.set_title('R=2 expression display dendrogram; not final consensus hierarchy')
            fig.tight_layout();pdf.savefig(fig);plt.close(fig)
    write_json(root/'RUN_STATE.json',dict(status='BENCHMARK_COMPLETE_IN_REVIEW',production_jobs=0,
        report='report/RESOURCE_REPORT.md',next_action='Await explicit user approval'))
    (root/'BENCHMARK_SUCCESS.txt').write_text('Two real iterations validated; stop for approval\n')


if __name__=='__main__':main()
