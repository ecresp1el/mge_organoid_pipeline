"""Allen consensus membership operations without a full cell-by-cell matrix.

Reference: scrattch.hicat 9af2f04, R/consensusCluster.R and R/annotate.R.
Rows always follow an explicit input-cell index. Every successful iteration
labels every cell, including held-out cells classified on fit-derived markers.
This module contains no genotype, annotation, cycle, or integration logic.

The public operations return numerical arrays and trace tables; the workflow
owns file publication. ``Membership.affinity`` is the exact factorization of
``get_cell.cl.co.ratio``. It is not a graph approximation or a label vote.
"""
import numpy as np
import pandas as pd
from scipy import sparse


def integer_labels(labels):
    """Validate a complete nonnegative integer partition; preserve numeric IDs."""
    values = np.asarray(labels)
    if values.ndim != 1 or not len(values):
        raise ValueError('Labels must be a nonempty vector')
    if not np.issubdtype(values.dtype, np.integer) or np.any(values < 0):
        raise ValueError('Missing/noninteger/negative cluster labels')
    return values.astype(np.int64, copy=False)


class Membership:
    """Own Allen's sparse B: iteration-cluster rows by full-population cells.

    Parameters
    ----------
    iterations : sequence of integer vectors
        Each vector has N complete labels, in the same cell order. Cluster
        numbers need not agree across iterations and are never matched.
    block_cells : int
        Computational row-block size, not a scientific sampling parameter.
    """

    def __init__(self, iterations, block_cells=4096):
        """Construct disjoint one-hot blocks with float64-safe accumulation."""
        if not len(iterations):
            raise ValueError('No successful consensus iterations')
        self.n = len(iterations[0])
        self.r = len(iterations)
        self.block_cells = int(block_cells)
        if self.block_cells < 1:
            raise ValueError('block_cells must be positive')
        indices = np.empty((self.n, self.r), dtype=np.int32)
        offset = 0
        self.block_index = []
        for iteration, labels in enumerate(iterations):
            labels = integer_labels(labels)
            if len(labels) != self.n:
                raise ValueError('Iteration does not label the complete population')
            unique, codes = np.unique(labels, return_inverse=True)
            indices[:, iteration] = codes + offset
            self.block_index.append(dict(iteration=iteration, row_start=offset,
                                         row_stop=offset+len(unique), labels=unique.tolist()))
            offset += len(unique)
        self.b = sparse.csc_matrix((np.ones(self.n*self.r, dtype=np.float64),
                                    indices.ravel(),
                                    np.arange(self.n+1, dtype=np.int64)*self.r),
                                   shape=(offset, self.n))

    def affinity(self, labels, target_rows=None, query_rows=None):
        """Return mean consensus to each target cluster, and sorted cluster IDs.

        ``labels`` describe ``target_rows`` (all cells when omitted).
        ``query_rows`` can request full-population predictions from a
        representative-cell partition, exactly as Allen's smaller branch.
        Self-pairs are included. No N-by-N array or sparse product is formed.
        Output memory is query_count x target_K float64, plus B_rows x K.
        """
        labels = integer_labels(labels)
        target_rows = np.arange(self.n) if target_rows is None else np.asarray(target_rows)
        query_rows = np.arange(self.n) if query_rows is None else np.asarray(query_rows)
        if len(target_rows) != len(labels) or len(np.unique(target_rows)) != len(target_rows):
            raise ValueError('Target rows must uniquely identify every supplied label')
        unique, codes, sizes = np.unique(labels, return_inverse=True, return_counts=True)
        h = sparse.csr_matrix((np.ones(len(labels)), (np.arange(len(labels)), codes)),
                              shape=(len(labels), len(unique)))
        totals = (self.b[:, target_rows] @ h).toarray()
        totals /= sizes[None, :]*self.r
        out = np.empty((len(query_rows), len(unique)), dtype=np.float64)
        for start in range(0, len(query_rows), self.block_cells):
            rows = query_rows[start:start+self.block_cells]
            out[start:start+len(rows)] = self.b[:, rows].T @ totals
        return out, unique

    def co_stats(self, labels):
        """Return Allen cohesion, separability, confusion and cluster medians.

        A single-cluster endpoint has between=0 by an explicit terminal guard;
        the R implementation can fail on an empty competitor matrix there.
        """
        labels = integer_labels(labels)
        a, unique = self.affinity(labels)
        positions = np.searchsorted(unique, labels)
        own = a[np.arange(len(labels)), positions].copy()
        a[np.arange(len(labels)), positions] = -np.inf
        between = np.max(a, axis=1) if len(unique) > 1 else np.zeros(len(labels))
        cells = pd.DataFrame(dict(cohesion=own, separability=own-between,
                                  confusion=np.divide(between, own, out=np.zeros_like(own), where=own!=0)))
        medians = cells.groupby(labels, sort=True).median()
        return cells, medians

    def merge_by_co(self, labels, diff_th=.25):
        """Port merge_cl_by_co, including R character-ID ordering and one pass.

        Eligible pairs require between > .1 and max(within)-between < diff_th.
        The source creates a column-major table and sorts its first character
        label descending before relabeling; recoding that as a numeric sort
        would change the result for IDs such as 2 and 10.
        """
        labels = integer_labels(labels).copy()
        a, unique = self.affinity(labels)
        means = np.vstack([a[labels==k].mean(axis=0) for k in unique])
        diagonal = np.diag(means)
        pairs = []
        for j, right in enumerate(unique):
            for i, left in enumerate(unique):
                gap = max(diagonal[i], diagonal[j])-means[i, j]
                if str(left) < str(right) and means[i, j] > .1 and gap < diff_th:
                    pairs.append(dict(destination=int(left), source=int(right),
                                      between=float(means[i,j]), within_destination=float(diagonal[i]),
                                      within_source=float(diagonal[j]), difference=float(gap)))
        pairs.sort(key=lambda row: str(row['destination']), reverse=True)
        for row in pairs:
            labels[labels==row['source']] = row['destination']
        return labels, pairs

    def refine(self, labels, confusion_th=.6, min_cells=20, niter=50, tol_th=.01):
        """Port refine_cl's update and stopping order, preserving every cell.

        Tolerance is tested BEFORE accepting proposed changes, as in R.
        Confused/small groups are removed as groups and their cells reassigned.
        Each returned trace row records a reassignment, stopping condition, or
        redistribution. No biological acceptance decision is made here.
        """
        labels = integer_labels(labels).copy()
        if len(labels) != self.n:
            raise ValueError('Refinement requires every population cell')
        trace = []
        outer = 0
        while len(np.unique(labels)) > 1:
            previous_correct = 0
            for iteration in range(niter):
                a, unique = self.affinity(labels)
                predicted = unique[a.argmax(axis=1)]
                correct = int(np.sum(predicted==labels))
                reason = 'update'
                if correct <= previous_correct:
                    reason = 'no_increase_in_correct_count'
                elif 1-correct/len(labels) < tol_th:
                    reason = 'below_tolerance_before_update'
                trace.append(dict(outer=outer, iteration=iteration, clusters=len(unique),
                                  unchanged_cells=correct, action=reason))
                if reason != 'update':
                    break
                previous_correct = correct
                labels = predicted.copy()
            cells, medians = self.co_stats(labels)
            sizes = pd.Series(labels).value_counts()
            remove = set(medians.index[medians.confusion > confusion_th])
            remove.update(sizes.index[sizes < min_cells])
            if not remove:
                break
            if len(remove) == len(sizes):
                labels[:] = labels.min()
                trace.append(dict(outer=outer, action='all_groups_redistributed_to_one'))
                break
            # R setdiff(unique(cl), rm.cl) preserves first-occurrence ordering.
            retained = np.array([k for k in pd.unique(labels) if k not in remove])
            bad = np.isin(labels, list(remove))
            columns = np.searchsorted(unique, retained)
            labels[bad] = retained[a[bad][:, columns].argmax(axis=1)]
            trace.append(dict(outer=outer, action='redistribute_groups', removed=sorted(map(int, remove)),
                              reassigned_cells=int(bad.sum())))
            outer += 1
            if outer > self.n:
                raise RuntimeError('Refinement did not converge; no accepted result')
        cells, medians = self.co_stats(labels)
        return labels, trace, cells, medians


def map_heldout(expression, sampled_rows, sampled_labels, marker_columns, block_cells=4096):
    """Classify held-out cells using Allen's mean/Pearson marker prototypes.

    Parameters
    ----------
    expression : sparse or dense cell-by-gene matrix
        Natural log(1+CPM), under the documented existing Python convention.
    sampled_rows : integer vector
        Exact rows used to fit sampled_labels; their order must agree.
    sampled_labels : integer vector
        Sampled-cell fit result, never overwritten by classifier predictions.
    marker_columns : integer vector
        Fit-derived separation-marker columns. No annotation marker panels.

    Returns
    -------
    all_labels, max_correlation, prototype_means, cluster_ids
        Sampled-cell scores are NaN, explicitly meaning 'not classified'.
        Undefined held-out correlations become zero, matching R. Sorted numeric
        cluster order resolves exact ties, matching numeric R cluster factors.
    """
    sampled_rows = np.asarray(sampled_rows, dtype=np.int64)
    labels = integer_labels(sampled_labels)
    marker_columns = np.asarray(marker_columns, dtype=np.int64)
    if len(labels)!=len(sampled_rows) or len(np.unique(sampled_rows))!=len(sampled_rows):
        raise ValueError('Sampled IDs and labels are inconsistent')
    unique = np.unique(labels)
    all_labels = np.full(expression.shape[0], -1, dtype=np.int32)
    all_labels[sampled_rows] = labels
    heldout = np.flatnonzero(all_labels < 0)
    scores = np.full(expression.shape[0], np.nan)
    if len(marker_columns) < 2 and len(unique)>1:
        raise ValueError('Insufficient fit-derived markers for Pearson assignment')
    prototypes = np.zeros((len(unique), len(marker_columns)), dtype=np.float64)
    for index, key in enumerate(unique):
        rows = sampled_rows[labels==key]
        for start in range(0, len(rows), block_cells):
            part = expression[rows[start:start+block_cells]][:, marker_columns]
            prototypes[index] += np.asarray(part.sum(axis=0)).ravel()
        prototypes[index] /= len(rows)
    if len(unique)==1:
        all_labels[heldout] = unique[0]
        # Unique possible label: no fictitious perfect correlation score.
        return all_labels, scores, prototypes, unique
    centered = prototypes-prototypes.mean(axis=1, keepdims=True)
    norm = np.sqrt((centered*centered).sum(axis=1))
    for start in range(0, len(heldout), block_cells):
        rows = heldout[start:start+block_cells]
        x = expression[rows][:, marker_columns]
        x = x.toarray() if sparse.issparse(x) else np.asarray(x).copy()
        x -= x.mean(axis=1, keepdims=True)
        denominator = np.sqrt((x*x).sum(axis=1))[:, None]*norm[None, :]
        correlations = np.divide(x @ centered.T, denominator,
                                 out=np.zeros((len(rows),len(unique))), where=denominator!=0)
        correlations[~np.isfinite(correlations)] = 0
        choice = correlations.argmax(axis=1)
        all_labels[rows] = unique[choice]
        scores[rows] = correlations[np.arange(len(rows)), choice]
    return all_labels, scores, prototypes, unique
