"""Instrument Allen Python HiCAT primitives for a bounded technical pilot.

This is an explicit orchestration adapter, not a reimplementation of DE/HVG/
graph algorithms. It records every recursive node, uses only expression for
fitting, handles a single community before the upstream return-type defect,
and performs a final cross-branch merge. See the protocol for differences
from the upstream convenience functions and R tutorial.
"""
from pathlib import Path
from itertools import combinations
import json
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
import transcriptomic_clustering as tc
from .provenance import write_json


class PilotEngine:
    """Own one candidate run's recursive clustering and saved node artifacts."""

    def __init__(self, config, candidate, directory, progress):
        """Store frozen controls and initialize trace collections; no fit yet."""
        self.cfg = config
        self.candidate = candidate
        self.directory = Path(directory)
        self.directory.mkdir(parents=True)
        self.progress = progress
        self.nodes = []
        self.pairs = []
        self.markers = set()
        self.leaves = []
        self.merge_serial = 0
        self.thresholds = dict(config['thresholds'])
        self.thresholds['score_thresh'] = candidate['score_thresh']
        self.root_projected = None
        self.root_graph = None
        self.premerge_labels = None

    def _table(self, path, frame, index=True):
        """Save a numerical TSV without silently dropping its identifying index."""
        frame.to_csv(path, sep='\t', index=index)

    def _capture_de(self, node, function, *args, **kwargs):
        """Call upstream DE unchanged and persist each evaluated pair's evidence.

        The saved score is the upstream uncapped sum of -log10 adjusted
        p-values; genes are algorithmic separation evidence, not annotations.
        The pair labels are local to this exact call/node, not final IDs.
        """
        table = function(*args, **kwargs)
        self.merge_serial += 1
        serialized = table.copy()
        for column in serialized.columns:
            if column.endswith('_genes'):
                serialized[column] = serialized[column].map(lambda x: json.dumps(list(x)))
        self._table(node/('de_evaluation_%04d.tsv'%self.merge_serial), serialized)
        assignments = getattr(self, '_active_assignments', {})
        write_json(node/('de_evaluation_%04d_membership.json'%self.merge_serial),
                   {str(k): list(map(int,v)) for k,v in assignments.items()})
        for pair,row in table.iterrows():
            self.pairs.append(dict(node=node.name, evaluation=self.merge_serial,
                                   cluster_a=str(pair[0]),cluster_b=str(pair[1]),
                                   score=float(row['score']), de_genes=int(row['num'])))
            for column in ('up_genes','down_genes'):
                self.markers.update(list(row[column])[:20])
        return table

    def _representation(self, obj, node, marker_mask=None):
        """Fit/persist local HVGs, PCA model and projected cells using Allen code.

        Parameters
        ----------
        obj : anndata.AnnData
            Dense ln(1+CPM) pilot cells by all genes, no fitted annotations.
        node : pathlib.Path
            Unique node directory for reusable model/feature artifacts.
        marker_mask : pandas.Series or None
            Final cross-branch marker union; None requests local HVG discovery.

        Returns
        -------
        anndata.AnnData
            Projected coordinates after the configured elbow filter.

        Notes
        -----
        All pilot cells fit PCA. This does not qualify full-data sparse/backed
        execution. No QC/genotype/sex modes are removed, and no old Step 06
        PCs or marker labels enter fitting. Zero retained PCs fail explicitly.
        """
        cfg=self.cfg
        with self.progress.track('hicat.features',candidate=self.candidate['name'],node=node.name,
                                 shape=list(obj.shape),low_thresh=cfg['low_thresh'],min_cells=cfg['min_gene_cells']) as event:
            if marker_mask is None:
                means,variances,eligible=tc.get_means_vars_genes(obj,low_thresh=cfg['low_thresh'],min_cells=cfg['min_gene_cells'])
                eligible=np.asarray(eligible,dtype=bool)
                if eligible.sum()<3 or not np.isfinite(means).all() or not np.isfinite(variances).all():
                    raise ValueError('Insufficient/invalid eligible gene statistics; not interpreted as a terminal cluster')
                hvgs=tc.highly_variable_genes(obj,means,variances,eligible,max_genes=cfg['max_hvgs'])
                stats=pd.DataFrame({'gene_id':obj.var_names[eligible],'linear_cpm_mean':means,'linear_cpm_variance':variances})
                self._table(node/'eligible_gene_statistics.tsv',stats,index=False)
            else:
                hvgs=marker_mask
            self._table(node/'feature_selection.tsv',pd.DataFrame({'gene_id':obj.var_names,'selected':np.asarray(hvgs,dtype=bool)}),index=False)
            event.update(selected_genes=int(np.asarray(hvgs).sum()))
        with self.progress.track('hicat.pca_and_projection',candidate=self.candidate['name'],node=node.name,
                                 n_comps=cfg['pca_components'],fit_cells=obj.n_obs,solver='randomized',
                                 pc_filter=cfg['pc_filter'],max_pcs=cfg['max_retained_pcs'],seed=self.candidate['seed']) as event:
            # Explicit indices avoid upstream cell_select=None -> len(slice) failure.
            components,ratios,variances,means=tc.pca(obj,gene_mask=hvgs,cell_select=np.arange(obj.n_obs),
                                                    n_comps=cfg['pca_components'],svd_solver='randomized',
                                                    random_state=self.candidate['seed'])
            self._table(node/'pca_components_all.tsv',components)
            self._table(node/'pca_center.tsv',means)
            self._table(node/'pca_variance.tsv',pd.DataFrame({'variance_ratio':ratios,'variance':variances}),index=False)
            selected=tc.dimension_reduction.filter_components(components,variances.copy(),ratios.copy(),
                                                              method=cfg['pc_filter'],max_pcs=cfg['max_retained_pcs'])
            if selected.shape[1]==0:
                raise ValueError('PC filter retained zero components; no one-cluster fallback')
            self._table(node/'pca_components_retained.tsv',selected)
            projected=tc.project(obj,selected,means)
            np.savez_compressed(node/'projection.npz',coordinates=projected.X,cell_ids=obj.obs_names.to_numpy(dtype=str))
            event.update(shape=list(projected.shape),saved_model='pca_components_retained.tsv + pca_center.tsv')
        return projected

    def _merge(self,obj,projected,groups,labels,node):
        """Merge using unchanged Allen criteria, handling its single-group API bug.

        Returns
        -------
        dict
            Local label -> local row positions after small-cluster and DE merge.

        Notes
        -----
        Only the instrumentation temporarily wraps tc.de_pairs_ebayes; the
        original function receives unchanged arguments and is restored on exit.
        n_markers=None skips redundant upstream marker selection; markers for
        final projection are collected from the persisted DE evaluations.
        """
        # Allen merge_two_clusters mutates membership using list.extend.
        # Accept array-valued orchestration inputs but pass native lists upstream.
        groups={key:list(map(int,values)) for key,values in groups.items()}
        if len(groups)==1:
            write_json(node/'single_cluster_terminal.json',{'reason':'single_graph_community','upstream_single_return_guard':True})
            return groups
        original=tc.de_pairs_ebayes
        active=tc.merging.merge_two_clusters
        active_de_merge=tc.merging.merge_clusters_by_de
        self._active_assignments=groups
        def traced_de(*args,**kwargs):
            """Persist the exact upstream DE result while preserving its return value."""
            return self._capture_de(node,original,*args,**kwargs)
        def traced_merge(assignments,*args,**kwargs):
            """Track current memberships so local pair labels can be resolved to cells."""
            self._active_assignments=assignments
            before={str(k):len(v) for k,v in assignments.items()}
            result=active(assignments,*args,**kwargs)
            with (node/'merge_actions.jsonl').open('a') as handle:
                handle.write(json.dumps(dict(source=str(args[0]),destination=str(args[1]),
                              sizes_before=before,sizes_after={str(k):len(v) for k,v in assignments.items()}))+'\n')
            return result
        def traced_de_merge(assignments,*args,**kwargs):
            """Capture the post-size-merge membership before evaluating DE pairs."""
            self._active_assignments=assignments
            return active_de_merge(assignments,*args,**kwargs)
        tc.merging.merge_clusters_by_de=traced_de_merge
        tc.de_pairs_ebayes=traced_de
        tc.merging.merge_two_clusters=traced_merge
        try:
            merged,_=tc.merge_clusters(obj,projected,groups,labels,thresholds=self.thresholds,
                                       k=self.cfg['merge_neighbors'],de_method=self.cfg['de_method'],
                                       n_markers=None,n_jobs=self.cfg['workers'])
            return merged
        finally:
            tc.de_pairs_ebayes=original
            tc.merging.merge_two_clusters=active
            tc.merging.merge_clusters_by_de=active_de_merge

    def _visit(self,full,rows,path='root',parent='',depth=0):
        """Visit one recursive branch and record why it splits or terminates.

        A depth guard is a failure, not successful convergence. Branches with
        fewer than min_recursive_cells stop by a stated size rule. Otherwise
        HVGs/PCA/graph/DE are recomputed using only that branch's cells.
        """
        node=self.directory/path
        node.mkdir()
        self._table(node/'cell_membership.tsv',pd.DataFrame({'local_row':np.arange(len(rows)),
                     'pilot_row':rows,'cell_id':full.obs_names[rows]}),index=False)
        row=dict(node=path,parent=parent,depth=depth,cells=len(rows),graph_clusters=0,post_merge_clusters=1,terminal_reason='')
        self.nodes.append(row)
        if len(rows)<self.cfg['min_recursive_cells']:
            row['terminal_reason']='below_min_recursive_cells'
            self.leaves.append(rows)
            write_json(node/'node_summary.json',row)
            return
        if depth>=self.cfg['max_depth_guard']:
            raise RuntimeError('Depth guard reached before convergence')
        obj=full[rows,:].copy()
        projected=self._representation(obj,node)
        # Preserve the pinned one-step API's unusual PC-dimension cap; report it.
        effective_k=min(self.cfg['graph_k'],projected.n_vars,obj.n_obs-1)
        with self.progress.track('hicat.graph',candidate=self.candidate['name'],node=path,cells=obj.n_obs,
                                 effective_k=effective_k,requested_k=self.cfg['graph_k'],metric=self.cfg['graph_metric'],
                                 backend=self.cfg['graph_backend'],weighting=self.cfg['graph_weighting'],
                                 seed=self.candidate['seed'],resolution=self.cfg['graph_resolution']) as event:
            labels,groups,graph,quality=tc.cluster_louvain(projected,k=effective_k,
                 nn_measure=self.cfg['graph_metric'],knn_method='annoy',louvain_method=self.cfg['graph_backend'],
                 weighting_method=self.cfg['graph_weighting'],annoy_trees=self.cfg['annoy_trees'],
                 n_jobs=self.cfg['workers'],resolution=self.cfg['graph_resolution'],
                 random_seed=self.candidate['seed'],annoy_index_filename=str(node/'annoy.index'))
            sparse.save_npz(node/'graph.npz',sparse.csr_matrix(graph))
            self._table(node/'graph_membership.tsv',pd.DataFrame({'cell_id':obj.obs_names,'graph_cluster':labels}),index=False)
            event.update(graph_clusters=len(groups),modularity=float(quality),graph_nnz=graph.nnz)
        if path=='root':
            self.root_projected=np.asarray(projected.X).copy()
            self.root_graph=sparse.csr_matrix(graph)
        row.update(graph_clusters=len(groups),effective_k=effective_k,retained_pcs=projected.n_vars)
        with self.progress.track('hicat.merge',candidate=self.candidate['name'],node=path,thresholds=self.thresholds,
                                 graph_clusters=len(groups)) as event:
            merged=self._merge(obj,projected,groups,labels,node)
            event.update(post_merge_clusters=len(merged),sizes=[len(x) for x in merged.values()])
        row['post_merge_clusters']=len(merged)
        write_json(node/'merged_membership.json',{str(k):list(map(int,v)) for k,v in merged.items()})
        if len(merged)==1:
            row['terminal_reason']='one_cluster_after_graph_or_de_merge'
            self.leaves.append(rows)
        else:
            row['terminal_reason']='split'
            ordered=sorted(merged.values(),key=lambda x:int(np.min(rows[np.asarray(x,dtype=int)])))
            for i,child in enumerate(ordered):
                self._visit(full,rows[np.asarray(child,dtype=int)],path+'_%03d'%i,path,depth+1)
        write_json(node/'node_summary.json',row)

    @staticmethod
    def labels(groups,n):
        """Return deterministic C0001-style labels after checking complete partition.

        Every row must occur exactly once. Final numbering sorts groups by their
        first pilot row and carries no biological or cross-run identity meaning.
        """
        groups=sorted(groups,key=lambda x:int(np.min(x)))
        rows=np.concatenate([np.asarray(x,dtype=int) for x in groups])
        if not np.array_equal(np.sort(rows),np.arange(n)):
            raise ValueError('Clustering is not an exact partition of the input rows')
        result=np.empty(n,dtype=object)
        for i,group in enumerate(groups):
            result[np.asarray(group,dtype=int)]='C%04d'%(i+1)
        return result

    def run(self,obj):
        """Run recursion, cross-branch merge and all-pair endpoint evidence.

        Returns
        -------
        tuple
            Final labels, summary dict, pairwise evidence, cluster means.

        Notes
        -----
        Final merge uses all pilot cells and a marker-union PCA. It does not
        call upstream final_merge's sampled-cell convenience wrapper. The
        final all-pair audit flags unresolved separations; it never silently
        changes the candidate partition or claims a final full-data count.
        """
        self._visit(obj,np.arange(obj.n_obs))
        self.premerge_labels=self.labels(self.leaves,obj.n_obs)
        final_groups=self.leaves
        node=self.directory/'final_merge'
        node.mkdir()
        self._table(node/'cell_membership.tsv',pd.DataFrame({'local_row':np.arange(obj.n_obs),
                     'pilot_row':np.arange(obj.n_obs),'cell_id':obj.obs_names}),index=False)
        if len(final_groups)>1:
            marker_mask=pd.Series(obj.var_names.isin(sorted(self.markers)),index=obj.var_names)
            if marker_mask.sum()<3:
                raise ValueError('Insufficient separation markers for final cross-branch PCA')
            projected=self._representation(obj,node,marker_mask)
            groups={i:np.asarray(g) for i,g in enumerate(final_groups)}
            numeric=np.zeros(obj.n_obs,dtype=int)
            for i,g in groups.items():numeric[g]=i
            with self.progress.track('hicat.final_cross_branch_merge',candidate=self.candidate['name'],
                                     recursive_leaf_count=len(groups),thresholds=self.thresholds) as event:
                final_groups=list(self._merge(obj,projected,groups,numeric,node).values())
                event.update(final_clusters=len(final_groups))
        labels=self.labels(final_groups,obj.n_obs)
        groups={label:np.flatnonzero(labels==label) for label in sorted(set(labels))}
        means,present,variances=tc.get_cluster_means(obj,groups,labels,low_th=self.cfg['low_thresh'])
        self._table(self.directory/'cluster_mean_log1p_cpm.tsv.gz',means)
        self._table(self.directory/'cluster_detection_fraction.tsv.gz',present)
        thresholds=dict(self.thresholds)
        for key in ['score_thresh','low_thresh','min_genes']:thresholds.pop(key)
        pairs=list(combinations(groups,2))
        if pairs:
            with self.progress.track('hicat.final_pairwise_audit',candidate=self.candidate['name'],pairs=len(pairs)):
                evidence=tc.de_pairs_ebayes(pairs,means,variances,present,{k:len(v) for k,v in groups.items()},thresholds)
            evidence['meets_declared_separation']=(evidence['score']>=self.thresholds['score_thresh']) & (evidence['num']>=self.thresholds['min_genes'])
        else:
            evidence=pd.DataFrame(columns=['score','num','up_genes','down_genes','meets_declared_separation'])
        serial=evidence.copy()
        for col in ['up_genes','down_genes']:
            serial[col]=serial[col].map(lambda v:json.dumps(list(v)))
        self._table(self.directory/'final_pairwise_evidence.tsv',serial)
        self._table(self.directory/'recursion_tree.tsv',pd.DataFrame(self.nodes),index=False)
        self._table(self.directory/'merge_evaluation_summary.tsv',pd.DataFrame(self.pairs),index=False)
        self._table(self.directory/'marker_union.tsv',pd.DataFrame({'gene_id':sorted(self.markers)}),index=False)
        summary=dict(candidate=self.candidate['name'],seed=self.candidate['seed'],score_threshold=self.candidate['score_thresh'],
                     cells=obj.n_obs,root_graph_clusters=self.nodes[0]['graph_clusters'],
                     root_post_merge_clusters=self.nodes[0]['post_merge_clusters'],
                     recursive_leaves=len(self.leaves),final_clusters=len(groups),
                     visited_nodes=len(self.nodes),max_depth=max(x['depth'] for x in self.nodes),
                     unresolved_final_pairs=int((~evidence['meets_declared_separation'].astype(bool)).sum()),
                     eligible_for_full_dataset_count=False)
        write_json(self.directory/'candidate_summary.json',summary)
        return labels,summary,evidence,means
