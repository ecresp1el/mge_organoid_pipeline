"""Execute the authorized 12,000-cell coarse/fine review checkpoint.

Input is the frozen raw-count subset prepared from approved Step 02. Normalize
once, fit coarse/fine hierarchies for two seeds, report all-gene markers and
plots for the baseline, validate saved objects, and atomically publish IN_REVIEW.
No older checkpoint is modified, no annotation is generated, and no batch
correction or full-data cluster-count claim is introduced.
"""
from pathlib import Path
import argparse
import gc
import json
import logging
import os
import sys
import importlib.metadata
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
import transcriptomic_clustering as tc
from .hierarchy import HierarchyEngine
from .markers import MarkerReporter
from .hierarchy_report import HierarchyReport
from .provenance import Progress, write_json, sha256, manifest


class HierarchyWorkflow:
    """Own a versioned larger pilot from verified raw input through review assets."""

    def __init__(self, run_dir):
        """Load the exact frozen settings and create a new staging directory."""
        self.run_dir=Path(run_dir)
        self.cfg=json.loads((self.run_dir/'config/hicat_pilot.json').read_text())
        self.stage=self.run_dir/'staging'
        self.stage.mkdir()
        self.progress=Progress(self.run_dir/'provenance')

    def run(self):
        """Fit hierarchies, export complete marker/plot assets, and validate I/O.

        The working expression is dense float64 for bounded Allen compatibility;
        a pre-allocation budget check and exact cell cap prevent accidental full
        dataset execution. The persisted normalized layer remains sparse.
        """
        cfg=self.cfg
        logging.getLogger().setLevel(logging.INFO)
        logging.getLogger().addHandler(logging.FileHandler(self.run_dir/'logs/upstream.log'))
        try:
            source=self.run_dir/'inputs/pilot_raw_counts.h5ad'
            identity=json.loads((self.run_dir/'inputs/input_identity.json').read_text())
            if sha256(source)!=identity['input_pilot_sha256']:raise ValueError('Input identity changed')
            if cfg['normalization_target']!=1000000 or cfg['normalization_log_base']!='natural':
                raise ValueError('Allen normalization contract changed')
            with self.progress.track('normalize',input=str(source),formula='ln(1 + raw/cell_total * 1e6)',batch_correction=False) as event:
                raw=ad.read_h5ad(source)
                if raw.n_obs!=cfg['pilot_cells_per_sample']*cfg['expected_samples'] or raw.n_vars!=cfg['expected_genes']:
                    raise ValueError('Exact authorized pilot dimensions required')
                if not raw.obs.technical_sample_id.value_counts().eq(cfg['pilot_cells_per_sample']).all():
                    raise ValueError('Unequal or incomplete per-sample selection')
                if not sparse.issparse(raw.X) or raw.X.dtype.kind not in 'iu' or (raw.X.data<0).any():
                    raise ValueError('Raw nonnegative sparse integer input required')
                totals=np.asarray(raw.X.sum(1)).ravel()
                if (totals<=0).any():raise ValueError('Zero-count input cell')
                # Reserve for several working copies, feature selection, and parent fits.
                budget_estimate=raw.n_obs*raw.n_vars*8*10/1024**3
                if budget_estimate>cfg['pilot_scratch_memory_gb']:
                    raise ValueError('Dense pilot memory estimate exceeds configured budget')
                norm=ad.AnnData(raw.X.toarray().astype(np.float64),obs=pd.DataFrame(index=raw.obs_names),var=raw.var.copy(),dtype=np.float64)
                tc.normalize(norm,inplace=True)
                if not np.isfinite(norm.X).all():raise ValueError('Nonfinite expression')
                raw.layers['log1p_cpm']=sparse.csr_matrix(norm.X)
                # Independently compare every stored nonzero with the documented formula.
                expected=raw.X.astype(np.float64).multiply((1e6/totals)[:,None]).tocsr()
                expected.data=np.log1p(expected.data)
                if not np.allclose(expected.data,raw.layers['log1p_cpm'].data,rtol=1e-12,atol=1e-12):
                    raise ValueError('Normalization formula mismatch')
                del expected
                event.update(cells=raw.n_obs,genes=raw.n_vars,dense_budget_estimate_gib=budget_estimate)
            tc.memory.set_memory_limit(GB=cfg['pilot_scratch_memory_gb'])
            tc.memory.allow_chunking=False
            summaries=[];baseline_parents=None
            for candidate in cfg['runs']:
                with self.progress.track('hierarchy.candidate',**candidate) as event:
                    engine=HierarchyEngine(cfg,candidate,self.stage/candidate['name'],self.progress)
                    coarse,fine,parents,summary=engine.run(norm)
                    for level,labels in [('coarse',coarse),('fine',fine)]:
                        raw.obs['hicat_'+level+'_'+candidate['name']]=pd.Categorical(labels)
                    if candidate['name']=='baseline':
                        baseline_parents=parents
                        raw.obsm['X_pca_hicat_coarse']=engine.coarse_engine.root_projected
                        raw.obsp['hicat_coarse_root_connectivities']=engine.coarse_engine.root_graph
                    summaries.append(summary);event.update(summary)
                    del engine
                    gc.collect()
            del norm
            gc.collect()
            summary=pd.DataFrame(summaries)
            summary.to_csv(self.stage/'cluster_count_summary.tsv',sep='\t',index=False)
            agreement=[]
            for level in ['coarse','fine']:
                first=raw.obs['hicat_'+level+'_baseline'].astype(str)
                repeat=raw.obs['hicat_'+level+'_seed_repeat'].astype(str)
                agreement.append(dict(level=level,adjusted_rand=adjusted_rand_score(first,repeat),
                                      normalized_mutual_information=normalized_mutual_info_score(first,repeat)))
                pd.crosstab(first,repeat).to_csv(self.stage/(level+'_seed_overlap.tsv'),sep='\t')
            agreement=pd.DataFrame(agreement)
            agreement.to_csv(self.stage/'partition_agreement.tsv',sep='\t',index=False)
            raw.obs.to_csv(self.stage/'cell_assignments.tsv.gz',sep='\t',index_label='cell_id')
            baseline_parents.to_csv(self.stage/'parent_child_mapping.tsv',sep='\t',index=False)
            report=HierarchyReport(self.stage/'figures',cfg['report_dpi'])
            try:
                report.overview(raw,baseline_parents,summary,agreement)
                all_top=[]
                for level in ['coarse','fine']:
                    context=level+'_global'
                    directory=self.stage/'markers'/context
                    reporter=MarkerReporter(cfg['marker_reporting'],directory,self.progress)
                    means,fractions,top=reporter.run(raw.layers['log1p_cpm'],raw.obs['hicat_'+level+'_baseline'].astype(str),raw.var,context)
                    all_top.append(top)
                    report.expression(context,directory,means,fractions,top)
                    report.composition(raw,level,directory)
                for parent in sorted(baseline_parents.coarse_cluster.unique()):
                    selected=np.flatnonzero(raw.obs.hicat_coarse_baseline.astype(str).to_numpy()==parent)
                    context='fine_within_'+parent
                    directory=self.stage/'markers'/context
                    reporter=MarkerReporter(cfg['marker_reporting'],directory,self.progress)
                    means,fractions,top=reporter.run(raw.layers['log1p_cpm'][selected],
                        raw.obs.hicat_fine_baseline.iloc[selected].astype(str),raw.var,context)
                    all_top.append(top)
                    report.expression(context,directory,means,fractions,top)
                pd.concat(all_top,ignore_index=True).to_csv(self.stage/'top20_genes_all_contexts.tsv',sep='\t',index=False)
            finally:
                report.close()
            raw.uns['hicat']={
                'schema_version':cfg['schema_version'],'stage':'01_coarse_fine_expanded_pilot',
                'status':'IN_REVIEW','run_id':self.run_dir.name,'scope':cfg['scope'],
                'source_step02_run_id':cfg['step02_run_id'],'source_step02_sha256':cfg['step02_sha256'],
                'upstream_commit':cfg['upstream_commit'],'pilot_cells':raw.n_obs,'genes':raw.n_vars,
                'full_dataset_count_estimated':False,'annotation_performed':False,'batch_correction_performed':False,
                'X_representation':'raw integer counts','normalized_layer':'log1p_cpm',
                'normalization':'natural log(1 + counts / all-gene cell total * 1000000)',
                'hierarchy':'fine clusters start from and remain inside their coarse parent; weak global pairs are flagged',
                'display_embedding':'X_umap_step06_display and X_pca_step06_display are inherited display only',
                'fitted_representation':'baseline coarse root PCA/graph in object; per-node models and fine graphs external',
                'marker_reporting':'all-gene exploratory Welch/BH; coarse vs rest, fine vs rest, fine vs siblings; distinct from HiCAT eBayes',
                'resolved_config_json':json.dumps(cfg,sort_keys=True),
                'candidate_summaries_json':json.dumps(summaries,sort_keys=True),
                'external_assets':'markers/, figures/, baseline/, seed_repeat/, parent_child_mapping.tsv',
                'stability_scope':'two seeds on same 12000 cells; no bootstrap consensus'}
            write_json(self.stage/'anndata_uns_inventory.json',raw.uns)
            with self.progress.track('serialize_and_validate',cells=raw.n_obs,genes=raw.n_vars) as event:
                output=self.stage/'pcdh19_hicat_coarse_fine.h5ad'
                raw.write_h5ad(output,compression='lzf')
                reopened=ad.read_h5ad(output)
                checks=dict(shape_preserved=reopened.shape==raw.shape,
                    raw_integer=reopened.X.dtype.kind in 'iu',raw_counts_exact=(raw.X!=reopened.X).nnz==0,
                    normalized_expression_exact=(raw.layers['log1p_cpm']!=reopened.layers['log1p_cpm']).nnz==0,
                    cell_ids_exact=raw.obs_names.equals(reopened.obs_names),gene_ids_exact=raw.var_names.equals(reopened.var_names),
                    inventory_matches=json.loads((self.stage/'anndata_uns_inventory.json').read_text())==dict(reopened.uns),
                    no_batch_correction=not reopened.uns['hicat']['batch_correction_performed'],
                    no_full_data_K=not reopened.uns['hicat']['full_dataset_count_estimated'],
                    no_annotations=not reopened.uns['hicat']['annotation_performed'],
                    per_sample_counts=reopened.obs.technical_sample_id.value_counts().eq(cfg['pilot_cells_per_sample']).all(),
                    display_coordinates_exact=np.array_equal(raw.obsm['X_umap_step06_display'],reopened.obsm['X_umap_step06_display']))
                for candidate in cfg['runs']:
                    name=candidate['name']
                    HierarchyEngine.validate_nesting(reopened.obs['hicat_coarse_'+name],reopened.obs['hicat_fine_'+name])
                    checks[name+'_nested']=True
                    for level in ['coarse','fine']:
                        key='hicat_'+level+'_'+name
                        checks[key+'_exact']=np.array_equal(raw.obs[key].astype(str),reopened.obs[key].astype(str))
                pd.DataFrame([dict(check=k,status='PASS' if v else 'FAIL') for k,v in checks.items()]).to_csv(self.stage/'validation_checks.tsv',sep='\t',index=False)
                if not all(checks.values()):raise ValueError('Saved-object validation failed')
                event.update(passed=len(checks),bytes=output.stat().st_size)
            write_json(self.stage/'software_versions.json',dict(python=sys.version,
                packages={d.metadata['Name']:d.version for d in importlib.metadata.distributions()}))
            write_json(self.stage/'STEP_STATUS.json',dict(status='IN_REVIEW',run_id=self.run_dir.name,
                stage='01_coarse_fine_expanded_pilot',full_data_clustering_started=False))
            (self.stage/'REVIEW_README.md').write_text('# Coarse/fine expanded pilot — IN_REVIEW\n\n'+
                '```text\n'+summary.to_csv(sep='\t',index=False)+'```'+'\n\nStart with figures/coarse_fine_review.pdf and top20_genes_all_contexts.tsv.\n'+
                'See markers/ for all-gene statistics, exact comparison definitions, heatmap matrices and dendrogram linkages.\n'+
                'See baseline/ and seed_repeat/ for the actual fitted hierarchy, models and HiCAT DE merge evidence.\n'+
                'Raw X and normalized log1p_cpm are both in pcdh19_hicat_coarse_fine.h5ad. No batch correction, annotation or full-data K.\n'+
                'Coarse/fine labels remain nested; unresolved global pairs are flagged, not silently merged across parents.\n')
            manifest(self.stage).to_csv(self.stage/'output_manifest.tsv',sep='\t',index=False)
            destination=self.run_dir/'outputs'
            if destination.exists():raise FileExistsError(str(destination))
            os.replace(self.stage,destination)
            self.progress.note('workflow','COMPLETE',{},dict(outputs=str(destination),status='IN_REVIEW'))
            return destination
        except Exception as error:
            self.progress.note('workflow','FAILED',{},dict(error_type=type(error).__name__,error=str(error)))
            raise


def main():
    """Run only the specified frozen expanded-pilot package."""
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir',type=Path,required=True)
    args=parser.parse_args()
    print(HierarchyWorkflow(args.run_dir).run(),flush=True)


if __name__=='__main__':
    main()
