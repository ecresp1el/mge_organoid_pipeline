"""Memory-bounded adapter for the already validated Allen Python fitter.

Executed full-data parameters, definitions and interpretation:
``paper3_pcdh19/HICAT_PARAMETERS_AND_CONSENSUS.md``. Frozen run copies
remain authoritative; documentation edits here do not change those copies.

Scientific settings and Allen's HVG/graph/eBayes implementations are preserved.
Changes are allocation-only: sparse sample statistics, forced all-node-cell
randomized PCA, blocked projection, and release of parent expression before
descending the recursive tree. No new clustering or preprocessing parameters
are selected here. The benchmark records this adapter's exact source hash.
"""
from contextlib import contextmanager
import gc
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.utils.sparsefuncs import mean_variance_axis
import transcriptomic_clustering as tc
from .engine import PilotEngine
from .provenance import write_json


def sparse_cluster_means(adata, cluster_assignments, cluster_by_obs=None,
                         chunk_size=None, low_th=1):
    """Compute original means, detection fractions and ddof=1 variances.

    Uses sklearn's stable sparse population variance and applies n/(n-1).
    Dense reduced coordinates retain NumPy's original calculation. A singleton
    has zero variance, as in the pinned Allen implementation. All assigned
    cells contribute; chunk_size never means a biological downsample.
    """
    means, detection, variances = [], [], []
    for rows in cluster_assignments.values():
        x = adata.X[np.asarray(rows, dtype=int), :]
        if sparse.issparse(x):
            mean, variance = mean_variance_axis(x, axis=0)
            variance = variance*len(rows)/(len(rows)-1) if len(rows)>1 else np.zeros(x.shape[1])
        else:
            mean = np.mean(x, axis=0)
            variance = np.var(x, axis=0, ddof=1) if len(rows)>1 else np.zeros(x.shape[1])
        means.append(np.asarray(mean).ravel())
        detection.append(np.asarray((x>low_th).mean(axis=0)).ravel())
        variances.append(np.asarray(variance).ravel())
    def frame(values):
        """Keep upstream cluster and gene ordering exactly."""
        return pd.DataFrame(np.vstack(values), index=list(cluster_assignments), columns=adata.var_names)
    return frame(means), frame(detection), frame(variances)


@contextmanager
def allocation_adapters():
    """Temporarily replace allocation paths; restore upstream APIs on exit.

    Must run in a single fitting process, never across concurrent Python
    threads. Independent SLURM processes have isolated interpreters.
    """
    original_stats, original_pca, original_project = tc.get_cluster_means, tc.pca, tc.project
    def pca_all_cells(adata, *args, **kwargs):
        """Prevent memory estimation from silently selecting IncrementalPCA."""
        kwargs['chunk_size'] = adata.n_obs
        return original_pca(adata, *args, **kwargs)
    def project_blocks(adata, *args, **kwargs):
        """Force the upstream sparse-safe projection branch, preserving PCs."""
        kwargs['chunk_size'] = max(1, min(4096, adata.n_obs-1))
        return original_project(adata, *args, **kwargs)
    tc.get_cluster_means, tc.pca, tc.project = sparse_cluster_means, pca_all_cells, project_blocks
    try:
        yield
    finally:
        tc.get_cluster_means, tc.pca, tc.project = original_stats, original_pca, original_project


class FullScaleEngine(PilotEngine):
    """Run the existing recursive fitter while releasing parent allocations."""

    def run(self, obj):
        """Execute the validated fitter with qualified allocation adapters."""
        with allocation_adapters():
            return super().run(obj)

    def _visit(self, full, rows, path='root', parent='', depth=0):
        """Fit one node, save its evidence, then release it before child fits.

        This mirrors PilotEngine._visit's rules and deterministic child order.
        ``pilot_row`` is retained as the legacy table field name; in this
        adapter it indexes the current full-scale fit, not the 12k pilot.
        """
        node = self.directory/path
        node.mkdir()
        self._table(node/'cell_membership.tsv', pd.DataFrame(dict(local_row=np.arange(len(rows)),
                    pilot_row=rows, cell_id=full.obs_names[rows])), index=False)
        row = dict(node=path,parent=parent,depth=depth,cells=len(rows),graph_clusters=0,
                   post_merge_clusters=1,terminal_reason='')
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
        effective_k=min(self.cfg['graph_k'],projected.n_vars,obj.n_obs-1)
        with self.progress.track('hicat.graph',candidate=self.candidate['name'],node=path,
                                 cells=obj.n_obs,effective_k=effective_k,seed=self.candidate['seed']) as event:
            labels,groups,graph,quality=tc.cluster_louvain(projected,k=effective_k,
                 nn_measure=self.cfg['graph_metric'],knn_method='annoy',louvain_method=self.cfg['graph_backend'],
                 weighting_method=self.cfg['graph_weighting'],annoy_trees=self.cfg['annoy_trees'],
                 n_jobs=self.cfg['workers'],resolution=self.cfg['graph_resolution'],
                 random_seed=self.candidate['seed'],annoy_index_filename=str(node/'annoy.index'))
            sparse.save_npz(node/'graph.npz',sparse.csr_matrix(graph))
            self._table(node/'graph_membership.tsv',pd.DataFrame(dict(cell_id=obj.obs_names,graph_cluster=labels)),index=False)
            event.update(graph_clusters=len(groups),modularity=float(quality),graph_nnz=int(graph.nnz))
        if path=='root':
            self.root_projected=np.asarray(projected.X).copy()
            self.root_graph=sparse.csr_matrix(graph)
        row.update(graph_clusters=len(groups),effective_k=effective_k,retained_pcs=projected.n_vars)
        with self.progress.track('hicat.merge',candidate=self.candidate['name'],node=path,
                                 thresholds=self.thresholds,graph_clusters=len(groups)) as event:
            merged=self._merge(obj,projected,groups,labels,node)
            event.update(post_merge_clusters=len(merged),sizes=[len(x) for x in merged.values()])
        row['post_merge_clusters']=len(merged)
        write_json(node/'merged_membership.json',{str(k):list(map(int,v)) for k,v in merged.items()})
        children=[rows[np.asarray(child,dtype=int)] for child in
                  sorted(merged.values(),key=lambda x:int(np.min(rows[np.asarray(x,dtype=int)])))]
        del obj, projected, graph, groups, merged, labels
        gc.collect()
        if len(children)==1:
            row['terminal_reason']='one_cluster_after_graph_or_de_merge'
            self.leaves.append(rows)
        else:
            row['terminal_reason']='split'
            for i, child in enumerate(children):
                self._visit(full,child,path+'_%03d'%i,path,depth+1)
        write_json(node/'node_summary.json',row)
