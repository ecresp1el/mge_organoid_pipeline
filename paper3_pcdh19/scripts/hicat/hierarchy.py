"""Coarse-to-fine orchestration that preserves the parent of every fine cluster.

A stringent one-step coarse partition is fitted once per seed. Independent
recursive fine fits start from its exact cell subsets with relaxed DE criteria.
Fine merges occur within parents. A global final pair audit flags cross-parent
weaknesses without silently destroying the saved hierarchy. No annotation,
batch correction, or target number of clusters enters fitting.
"""
from itertools import combinations
from pathlib import Path
import copy
import json
import numpy as np
import pandas as pd
import transcriptomic_clustering as tc
from .engine import PilotEngine
from .provenance import write_json


class CoarseEngine(PilotEngine):
    """Reuse the instrumented Allen pass but stop before fitting child branches."""

    def _visit(self, full, rows, path='root', parent='', depth=0):
        """Fit the root normally; record its children as coarse terminal groups."""
        if depth == 0:
            return super()._visit(full, rows, path, parent, depth)
        node=self.directory/path
        node.mkdir()
        self._table(node/'cell_membership.tsv',pd.DataFrame({
            'local_row':np.arange(len(rows)), 'pilot_row':rows,
            'cell_id':full.obs_names[rows]}),index=False)
        record=dict(node=path,parent=parent,depth=depth,cells=len(rows),
                    graph_clusters=0,post_merge_clusters=1,
                    terminal_reason='coarse_one_step_boundary')
        self.nodes.append(record)
        self.leaves.append(rows)
        write_json(node/'node_summary.json',record)


class HierarchyEngine:
    """Fit and save a paired coarse/fine hierarchy for one random seed.

    Parameters
    ----------
    config : dict
        Common controls plus independent coarse/fine threshold dictionaries.
    candidate : dict
        Name and seed; the seed affects graph/PCA fitting, not cell selection.
    directory : pathlib.Path
        New candidate directory. Existing outputs are never replaced.
    progress : Progress
        Shared durable operation logger.
    """

    def __init__(self, config, candidate, directory, progress):
        """Resolve candidate state without fitting or modifying input data."""
        self.cfg=config
        self.candidate=candidate
        self.directory=Path(directory)
        self.directory.mkdir()
        self.progress=progress

    def _config(self, level):
        """Return an independent level config so relaxed settings cannot leak."""
        cfg=copy.deepcopy(self.cfg)
        cfg['thresholds']=dict(cfg[level+'_thresholds'])
        return cfg

    def run(self, norm):
        """Return coarse/fine labels, parent table, and global final DE audit.

        Every fine group is nested in exactly one coarse parent. A one-cluster
        parent remains a valid unsplit result. Within-parent final merging is
        completed before the global audit, which reports rather than merges
        weak pairs across parents. Neither count is a full-data estimate.
        """
        coarse_cfg=self._config('coarse')
        candidate=dict(self.candidate,score_thresh=coarse_cfg['thresholds']['score_thresh'])
        self.coarse_engine=CoarseEngine(coarse_cfg,candidate,self.directory/'coarse',self.progress)
        coarse,coarse_summary,_,_=self.coarse_engine.run(norm)
        fine=np.empty(norm.n_obs,dtype=object)
        rows=[]
        for parent in sorted(set(coarse)):
            selected=np.flatnonzero(coarse==parent)
            fine_cfg=self._config('fine')
            fine_candidate=dict(self.candidate,name=self.candidate['name']+'/'+parent,
                                score_thresh=fine_cfg['thresholds']['score_thresh'])
            with self.progress.track('hierarchy.fine_parent',candidate=self.candidate['name'],
                                     parent=parent,cells=len(selected),thresholds=fine_cfg['thresholds']) as event:
                engine=PilotEngine(fine_cfg,fine_candidate,self.directory/('fine_'+parent),self.progress)
                local,summary,_,_=engine.run(norm[selected,:].copy())
                for child in sorted(set(local)):
                    label=parent+'.F'+child[1:]
                    fine[selected[local==child]]=label
                    rows.append(dict(coarse_cluster=parent,fine_cluster=label,
                                     cells=int((local==child).sum()),
                                     parent_cells=len(selected),fine_count_in_parent=len(set(local))))
                event.update(summary)
        parent_table=pd.DataFrame(rows)
        self.validate_nesting(coarse,fine)
        parent_table.to_csv(self.directory/'parent_child_mapping.tsv',sep='\t',index=False)
        pd.DataFrame({'cell_id':norm.obs_names,'coarse_cluster':coarse,'fine_cluster':fine}).to_csv(
            self.directory/'cell_hierarchy.tsv.gz',sep='\t',index=False)
        groups={label:np.flatnonzero(fine==label) for label in sorted(set(fine))}
        with self.progress.track('hierarchy.global_fine_pair_audit',candidate=self.candidate['name'],
                                 clusters=len(groups),preserve_coarse_parents=True):
            means,present,variances=tc.get_cluster_means(norm,groups,fine,low_th=self.cfg['low_thresh'])
            thresholds=dict(self.cfg['fine_thresholds'])
            for key in ('score_thresh','low_thresh','min_genes'):thresholds.pop(key)
            pairs=list(combinations(groups,2))
            if pairs:
                evidence=tc.de_pairs_ebayes(pairs,means,variances,present,
                                          {key:len(value) for key,value in groups.items()},thresholds)
                evidence['same_coarse_parent']=[a.split('.')[0]==b.split('.')[0] for a,b in evidence.index]
                evidence['meets_declared_separation']=(evidence['score']>=self.cfg['fine_thresholds']['score_thresh']) & (evidence['num']>=self.cfg['fine_thresholds']['min_genes'])
            else:
                evidence=pd.DataFrame(columns=['score','num','up_genes','down_genes',
                                              'same_coarse_parent','meets_declared_separation'])
            serial=evidence.copy()
            for column in ('up_genes','down_genes'):
                serial[column]=serial[column].map(lambda genes:json.dumps(list(genes)))
            serial.to_csv(self.directory/'global_fine_pairwise_evidence.tsv',sep='\t')
        summary=dict(candidate=self.candidate['name'],seed=self.candidate['seed'],cells=norm.n_obs,
                     coarse_clusters=len(set(coarse)),fine_clusters=len(groups),
                     coarse_unresolved_pairs=coarse_summary['unresolved_final_pairs'],
                     fine_unresolved_pairs=int((~evidence['meets_declared_separation'].astype(bool)).sum()),
                     full_dataset_count_estimated=False)
        write_json(self.directory/'hierarchy_summary.json',summary)
        return coarse,fine,parent_table,summary

    @staticmethod
    def validate_nesting(coarse, fine):
        """Reject missing assignments or any fine cluster spanning two parents."""
        frame=pd.DataFrame({'coarse':coarse,'fine':fine})
        if frame.isna().any().any() or not frame.groupby('fine').coarse.nunique().eq(1).all():
            raise ValueError('Every fine cluster must have exactly one coarse parent')
