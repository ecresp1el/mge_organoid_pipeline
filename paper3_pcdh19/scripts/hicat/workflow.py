"""Run, validate and publish the explicitly authorized HiCAT technical pilot.

Raw counts and normalized expression are both persisted. Each candidate gets
its own node models/graphs/evidence. Publication is a single output-directory
rename after structural validation. The full dataset is never clustered here.
"""
from pathlib import Path
from itertools import combinations
import ast
import json
import logging
import os
import sys
import importlib.metadata
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from scipy.cluster.hierarchy import linkage
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
import transcriptomic_clustering as tc
from .engine import PilotEngine
from .report import PilotReport
from .provenance import Progress, sha256, write_json, manifest


class PilotWorkflow:
    """Keep normalization, candidate fitting, validation and publication explicit."""

    def __init__(self,run_dir):
        """Resolve frozen config and create new staging/progress state."""
        self.run_dir=Path(run_dir)
        self.cfg=json.loads((self.run_dir/'config/hicat_pilot.json').read_text())
        if self.cfg['normalization_target'] != 1000000 or self.cfg['normalization_log_base'] != 'natural':
            raise ValueError('Pinned Allen normalizer requires natural log(1+CPM); these are contract fields')
        if not self.cfg['final_cross_branch_merge'] or not self.cfg['pairwise_final_validation']:
            raise ValueError('This review contract requires final merge and pairwise audit')
        if [x['name'] for x in self.cfg['runs']] != ['baseline','seed_repeat','stricter_merge']:
            raise ValueError('This pilot report requires the three declared candidate roles')
        self.stage=self.run_dir/'staging'
        self.stage.mkdir()
        self.progress=Progress(self.run_dir/'provenance')

    def qualify(self):
        """Record source edge-case probes before any PCDH19 fitting.

        Expected defects are recorded as observed, not silently patched in the
        upstream checkout. The adapter deliberately avoids those API paths.
        """
        rows=[]
        raw=np.array([[1.,2.,0.],[2.,1.,1.]],dtype=np.float64)
        fixture=ad.AnnData(raw.copy(),dtype=np.float64)
        normalized=tc.normalize(fixture)
        error=float(np.max(np.abs(normalized.X-np.log1p(raw/raw.sum(axis=1)[:,None]*1e6))))
        if error>1e-10:raise ValueError('Upstream normalization differs from declared ln(1+CPM)')
        rows.append(dict(probe='normalization_formula',observation=error,disposition='PASS'))
        result=tc.merge_clusters(None,None,{0:[0,1]},np.array([0,0]))
        rows.append(dict(probe='upstream_single_cluster_return',observation=type(result).__name__,
                         disposition='ADAPTER_GUARD: preserve sole group before calling merge'))
        pcs=pd.DataFrame(np.eye(3))
        filtered=tc.dimension_reduction.filter_components(pcs,np.array([3.,2.,1.]),np.array([.5,.3,.2]),method=None)
        rows.append(dict(probe='upstream_pc_filter_none',observation=str(filtered.shape),
                         disposition='AVOIDED: use explicit elbow method, fail on zero PCs'))
        test_groups=[np.array([2,3]),np.array([0,1])]
        if list(PilotEngine.labels(test_groups,4))!=['C0001','C0001','C0002','C0002']:
            raise ValueError('Stable cluster indexing check failed')
        try:
            PilotEngine.labels([np.array([0,1]),np.array([1,2])],4)
        except ValueError:
            rows.append(dict(probe='duplicate_missing_cell_guard',observation='rejected',disposition='PASS'))
        else:raise ValueError('Invalid partition was accepted')
        pd.DataFrame(rows).to_csv(self.stage/'qualification_checks.tsv',sep='\t',index=False)

    def run(self):
        """Execute three bounded candidate runs and stop at IN_REVIEW.

        Returns
        -------
        pathlib.Path
            Published output directory. A failure leaves staging and events
            for diagnosis and does not publish a successful checkpoint.

        Notes
        -----
        The public source's algorithms are called from PilotEngine. Model
        comparisons are descriptive pilot diagnostics, not bootstrap consensus,
        biological replication, or proof of a single correct number of clusters.
        """
        logging.getLogger().setLevel(logging.INFO)
        handler=logging.FileHandler(self.run_dir/'logs/upstream.log')
        handler.setFormatter(logging.Formatter('%(asctime)s %(name)s %(levelname)s %(message)s'))
        logging.getLogger().addHandler(handler)
        try:
            with self.progress.track('qualification',upstream_commit=self.cfg['upstream_commit']):self.qualify()
            source=self.run_dir/'inputs/pilot_raw_counts.h5ad'
            identity=json.loads((self.run_dir/'inputs/input_identity.json').read_text())
            if sha256(source)!=identity['input_pilot_sha256']:raise ValueError('Pilot input identity changed')
            with self.progress.track('load_and_normalize',input=str(source),target_sum=1e6,log_base='natural') as event:
                raw=ad.read_h5ad(source)
                if raw.n_obs>self.cfg['pilot_cells_per_sample']*self.cfg['expected_samples']:
                    raise ValueError('Refusing input larger than authorized technical pilot')
                if raw.n_vars!=self.cfg['expected_genes'] or not sparse.issparse(raw.X) or raw.X.dtype.kind not in 'iu':
                    raise ValueError('Pilot input must preserve all genes and sparse integer raw counts')
                if np.any(np.asarray(raw.X.sum(axis=1)).ravel()<=0):raise ValueError('Zero-count pilot cell')
                # Dense memory is deliberately bounded to the pilot; do not reuse for full data.
                working=ad.AnnData(raw.X.toarray().astype(np.float64),
                                   obs=pd.DataFrame(index=raw.obs_names),var=raw.var.copy(),dtype=np.float64)
                norm=tc.normalize(working)
                if not np.isfinite(norm.X).all():raise ValueError('Nonfinite normalized expression')
                raw.layers['log1p_cpm']=sparse.csr_matrix(norm.X)
                event.update(shape=list(norm.shape),sparse_raw_dtype=str(raw.X.dtype),normalized_dtype=str(norm.X.dtype))
            tc.memory.set_memory_limit(GB=self.cfg['pilot_scratch_memory_gb'])
            tc.memory.allow_chunking=False
            summaries=[];engines={};evidence={};means={}
            for candidate in self.cfg['runs']:
                with self.progress.track('candidate',**candidate) as event:
                    engine=PilotEngine(self.cfg,candidate,self.stage/candidate['name'],self.progress)
                    labels,summary,pairs,centers=engine.run(norm)
                    raw.obs['hicat_'+candidate['name']]=pd.Categorical(labels)
                    raw.obs['hicat_'+candidate['name']+'_recursive']=pd.Categorical(engine.premerge_labels)
                    engines[candidate['name']]=engine;evidence[candidate['name']]=pairs;means[candidate['name']]=centers
                    summaries.append(summary);event.update(summary)
            summary_frame=pd.DataFrame(summaries)
            assignments=raw.obs.copy()
            assignments.index.name='cell_id'
            assignments.to_csv(self.stage/'cell_assignments.tsv.gz',sep='\t')
            summary_frame.to_csv(self.stage/'cluster_count_summary.tsv',sep='\t',index=False)
            baseline=raw.obs['hicat_baseline'].astype(str)
            agreement=[];stability=[]
            for candidate in self.cfg['runs'][1:]:
                labels=raw.obs['hicat_'+candidate['name']].astype(str)
                agreement.append(dict(candidate=candidate['name'],adjusted_rand=adjusted_rand_score(baseline,labels),
                                      normalized_mutual_information=normalized_mutual_info_score(baseline,labels)))
                for group in sorted(baseline.unique()):
                    a=baseline.to_numpy()==group
                    overlaps={target:float(np.sum(a & (labels.to_numpy()==target))/np.sum(a | (labels.to_numpy()==target))) for target in labels.unique()}
                    target=max(overlaps,key=overlaps.get)
                    stability.append(dict(candidate=candidate['name'],baseline_cluster=group,best_match=target,best_jaccard=overlaps[target]))
            stability=pd.DataFrame(stability)
            stability.to_csv(self.stage/'cluster_stability.tsv',sep='\t',index=False)
            pd.DataFrame(agreement).to_csv(self.stage/'partition_agreement.tsv',sep='\t',index=False)
            cross=pd.crosstab(raw.obs['technical_sample_id'],raw.obs['hicat_baseline'])
            cross.to_csv(self.stage/'sample_cluster_counts.tsv',sep='\t')
            cross.div(cross.sum(axis=1),axis=0).to_csv(self.stage/'sample_cluster_within_sample_fractions.tsv',sep='\t')
            cross.div(cross.sum(axis=0),axis=1).to_csv(self.stage/'sample_cluster_within_cluster_fractions.tsv',sep='\t')
            contingency=pd.crosstab(raw.obs['hicat_baseline'],raw.obs['hicat_seed_repeat'])
            contingency.to_csv(self.stage/'baseline_seed_overlap.tsv',sep='\t')
            baseline_engine=engines['baseline']
            raw.obsm['X_pca_hicat_root']=baseline_engine.root_projected
            raw.obsp['hicat_root_connectivities']=baseline_engine.root_graph
            # Clear schema: inherited display coordinates are labeled explicitly, and
            # all current-stage metadata is under hicat rather than a stale latest_step.
            raw.uns['hicat']={
                'schema_version':'pcdh19_hicat_pilot_v1','stage':'00_method_and_technical_pilot','status':'IN_REVIEW',
                'run_id':self.run_dir.name,'scope':'technical_pilot_not_full_dataset_cluster_count',
                'full_dataset_count_estimated':False,'annotation_performed':False,
                'upstream_commit':self.cfg['upstream_commit'],'source_step02_run_id':self.cfg['step02_run_id'],
                'source_step02_sha256':self.cfg['step02_sha256'],'source_full_cells':self.cfg['expected_cells'],
                'pilot_cells':raw.n_obs,'genes':raw.n_vars,'cells_removed_within_pilot':0,
                'sampling':'100 cells per registered sample; exact IDs in inputs/pilot_cell_selection.tsv',
                'X_representation':'raw integer counts','normalized_layer':'log1p_cpm',
                'normalization':'natural log(1 + counts / all-gene cell total * 1000000)',
                'display_embedding':'X_umap_step06_display is inherited, not fitted or used for clustering',
                'fitting_graph':'hicat_root_connectivities is baseline root graph, not a universal recursive graph',
                'resolved_config_json':json.dumps(self.cfg,sort_keys=True),
                'candidate_summaries_json':json.dumps(summaries,sort_keys=True),
                'adapter_contract':'Explicit Allen primitive orchestration; see frozen protocol and engine.py',
                'model_assets':'candidate/node/pca_components_retained.tsv + pca_center.tsv + feature_selection.tsv',
                'stability_scope':'one seed repeat and one merge-threshold comparison; no bootstrap consensus'}
            write_json(self.stage/'anndata_uns_inventory.json',raw.uns)
            if len(means['baseline'])>1:
                np.save(self.stage/'baseline_cluster_mean_linkage.npy',linkage(means['baseline'],method='average',metric='euclidean'))
                pd.Series(means['baseline'].index,name='cluster_id').to_csv(self.stage/'baseline_linkage_labels.tsv',sep='\t',index=False)
            with self.progress.track('report',cells=raw.n_obs,panels='A-L',annotation=False) as event:
                figures=PilotReport(self.stage/'figures',self.cfg).publish(raw,summary_frame,evidence['baseline'],
                    means['baseline'],pd.DataFrame(baseline_engine.nodes),stability,contingency)
                event.update(assets=figures)
            output=self.stage/'pcdh19_hicat_pilot.h5ad'
            with self.progress.track('serialize_and_validate',path=str(output),shape=list(raw.shape)) as event:
                raw.write_h5ad(output,compression='lzf')
                reopened=ad.read_h5ad(output)
                inventory=json.loads((self.stage/'anndata_uns_inventory.json').read_text())
                checks={
                    'shape_preserved':reopened.shape==raw.shape,
                    'raw_integer':reopened.X.dtype.kind in 'iu',
                    'raw_counts_exact':(reopened.X!=raw.X).nnz==0,
                    'normalization_exact':(reopened.layers['log1p_cpm']!=raw.layers['log1p_cpm']).nnz==0,
                    'cell_ids_exact':reopened.obs_names.equals(raw.obs_names),
                    'gene_ids_exact':reopened.var_names.equals(raw.var_names),
                    'all_samples_present':reopened.obs['technical_sample_id'].nunique()==self.cfg['expected_samples'],
                    'no_annotations_created':'step06_provisional_state' not in reopened.obs,
                    'scope_not_full_count':not reopened.uns['hicat']['full_dataset_count_estimated'],
                    'finite_root_pca':bool(np.isfinite(reopened.obsm['X_pca_hicat_root']).all()),
                    'root_graph_shape':reopened.obsp['hicat_root_connectivities'].shape==(raw.n_obs,raw.n_obs),
                    'uns_inventory_structured':isinstance(inventory,dict),
                    'uns_inventory_matches_saved':inventory==dict(reopened.uns),
                }
                for candidate in self.cfg['runs']:
                    key='hicat_'+candidate['name']
                    checks[key+'_complete']=not reopened.obs[key].isna().any()
                    checks[key+'_labels_exact']=np.array_equal(reopened.obs[key].astype(str),raw.obs[key].astype(str))
                frame=pd.DataFrame([dict(check=k,status='PASS' if v else 'FAIL') for k,v in checks.items()])
                frame.to_csv(self.stage/'validation_checks.tsv',sep='\t',index=False)
                if not all(checks.values()):raise ValueError('Published-object validation failed')
                event.update(passed=len(checks),failed=0,bytes=output.stat().st_size)
            versions={d.metadata['Name']:d.version for d in importlib.metadata.distributions()}
            write_json(self.stage/'software_versions.json',dict(python=sys.version,packages=versions))
            summary_text=summary_frame.to_csv(sep='\t',index=False)
            report='# HiCAT technical pilot — IN_REVIEW\n\nThis is a 1,200-cell software/method pilot, not a full-data cluster-count estimate.\n\n'+\
                '```text\n'+summary_text+'```\n\nMain object: `pcdh19_hicat_pilot.h5ad`, with raw X, saved log1p_cpm layer, all candidate cluster columns, and explicit hicat metadata.\n\n'+\
                'Review plots: `figures/hicat_pilot_review.pdf`. No cell-type annotation, integration, or cell removal was performed.\n\n'+\
                'Keep node models, graphs, memberships, DE evaluations, frozen code/config, and the input-selection ledger with the object. Root PCA is not a single model for every recursive split.\n\n'+\
                'Check unresolved_final_pairs and seed/threshold variation before designing full-data runs. The final all-pair audit reports weak separations without silently changing clusters.\n'
            (self.stage/'PILOT_REPORT.md').write_text(report)
            write_json(self.stage/'STEP_STATUS.json',dict(status='IN_REVIEW',stage='00_method_and_technical_pilot',
                       run_id=self.run_dir.name,full_data_clustering_started=False,annotation_performed=False))
            manifest(self.stage).to_csv(self.stage/'output_manifest.tsv',sep='\t',index=False)
            destination=self.run_dir/'outputs'
            if destination.exists():raise FileExistsError(str(destination))
            os.replace(self.stage,destination)
            self.progress.note('workflow','COMPLETE',{},dict(outputs=str(destination),status='IN_REVIEW'))
            return destination
        except Exception as error:
            self.progress.note('workflow','FAILED',{},dict(error_type=type(error).__name__,error=str(error)))
            raise
