"""Save all-gene cluster-versus-rest statistics and qualified top-20 genes.

These descriptive cell-level Welch tests explain discovered partitions; they
are separate from Allen's pairwise eBayes merge evidence. Testing cells used to
form clusters is exploratory, not independent confirmation or sample-replicated
genotype DE. All genes are tested before BH correction within each contrast.
"""
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import ttest_ind_from_stats
from .provenance import write_json


class MarkerReporter:
    """Rank positive separation genes while retaining all tested gene statistics."""

    def __init__(self, config, directory, progress):
        """Store explicit marker filters, separate from clustering thresholds."""
        self.cfg=config
        self.directory=Path(directory)
        self.directory.mkdir(parents=True)
        self.progress=progress

    @staticmethod
    def bh(pvalues):
        """Benjamini-Hochberg adjusted p-values in original gene order."""
        p=np.asarray(pvalues,dtype=float)
        order=np.argsort(p,kind='stable')
        result=np.empty(len(p))
        result[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
        return result

    @staticmethod
    def moments(matrix):
        """Return n, sums, squared sums and nonzero counts without dense copies."""
        matrix=sparse.csr_matrix(matrix)
        return (matrix.shape[0],np.asarray(matrix.sum(0)).ravel(),
                np.asarray(matrix.multiply(matrix).sum(0)).ravel(),
                np.asarray((matrix>0).sum(0)).ravel())

    @staticmethod
    def stats(moment):
        """Convert sufficient statistics to means, unbiased variance, detection."""
        n,total,squares,detected=moment
        mean=total/n
        variance=np.maximum((squares-total*mean)/max(n-1,1),0)
        return n,mean,variance,detected/n

    def run(self, matrix, labels, var, context):
        """Save every contrast and return means, detection fractions and top genes.

        Parameters
        ----------
        matrix : scipy.sparse matrix
            Natural log1p(CPM), restricted to this comparison universe.
        labels : array-like
            Coarse IDs globally, fine IDs globally, or sibling IDs in a parent.
        var : pandas.DataFrame
            Gene IDs as index and gene_symbol column; IDs disambiguate symbols.
        context : str
            Explicit comparison universe, persisted on every statistic row.

        Notes
        -----
        Positive genes require BH FDR <= configured maximum, mean natural-log
        difference >= minimum, and foreground detection >= minimum. Qualified
        genes rank by descending Welch t, then mean difference, then gene ID.
        Fewer than 20 qualifying genes are reported honestly, never padded.
        A sole cluster has no rest population and receives no invented DE test.
        """
        labels=np.asarray(labels)
        matrix=sparse.csr_matrix(matrix)
        total=self.moments(matrix)
        means={};fractions={};top_tables=[];coverage=[]
        with self.progress.track('markers',context=context,cells=matrix.shape[0],
                                 clusters=len(set(labels)),comparison='cluster_vs_rest',filters=self.cfg):
            for label in sorted(set(labels)):
                foreground=self.moments(matrix[labels==label])
                background=tuple(a-b for a,b in zip(total,foreground))
                n1,mean1,var1,frac1=self.stats(foreground)
                means[label]=mean1;fractions[label]=frac1
                if background[0]<2 or n1<2:
                    coverage.append(dict(cluster=label,context=context,eligible_genes=0,reported_genes=0,
                                         status='not_testable_no_rest_or_insufficient_cells'))
                    continue
                n2,mean2,var2,frac2=self.stats(background)
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore',RuntimeWarning)
                    score,pvalue=ttest_ind_from_stats(mean1,np.sqrt(var1),n1,mean2,np.sqrt(var2),n2,equal_var=False)
                # Both groups constant and identical carry zero evidence.
                same_constant=(var1==0)&(var2==0)&(mean1==mean2)
                score[same_constant]=0;pvalue[same_constant]=1
                if np.isnan(score).any() or not np.isfinite(pvalue).all():
                    raise ValueError('Undefined marker test statistics in '+context+'/'+str(label))
                frame=pd.DataFrame(dict(gene_id=var.index.astype(str),gene_symbol=var.gene_symbol.astype(str).to_numpy(),
                    context=context,cluster=label,n_foreground=n1,n_background=n2,
                    mean_log1p_cpm_foreground=mean1,mean_log1p_cpm_background=mean2,
                    delta_mean_log1p_cpm=mean1-mean2,fraction_detected_foreground=frac1,
                    fraction_detected_background=frac2,welch_t=score,p_value=pvalue,p_adjusted_bh=self.bh(pvalue)))
                frame['qualifies_positive_marker']=(frame.p_adjusted_bh<=self.cfg['max_fdr']) & (
                    frame.delta_mean_log1p_cpm>=self.cfg['min_delta_mean_log1p_cpm']) & (
                    frame.fraction_detected_foreground>=self.cfg['min_fraction_detected'])
                frame=frame.sort_values(['welch_t','delta_mean_log1p_cpm','gene_id'],ascending=[False,False,True])
                frame['rank_all_genes']=np.arange(1,len(frame)+1)
                frame.to_csv(self.directory/(str(label)+'_all_genes.tsv.gz'),sep='\t',index=False)
                top=frame[frame.qualifies_positive_marker].head(self.cfg['top_n']).copy()
                top['rank_positive_marker']=np.arange(1,len(top)+1)
                top_tables.append(top)
                coverage.append(dict(cluster=label,context=context,eligible_genes=int(frame.qualifies_positive_marker.sum()),
                                     reported_genes=len(top),status='tested'))
                self.progress.note('markers.top_genes','COMPLETE',dict(context=context,cluster=label),
                                   dict(genes=top.gene_symbol.tolist(),gene_ids=top.gene_id.tolist(),reported=len(top)))
        top=pd.concat(top_tables,ignore_index=True) if top_tables else pd.DataFrame(columns=[
            'cluster','gene_id','gene_symbol','context','rank_positive_marker'])
        top.to_csv(self.directory/'top20_genes.tsv',sep='\t',index=False)
        pd.DataFrame(coverage).to_csv(self.directory/'marker_coverage.tsv',sep='\t',index=False)
        means=pd.DataFrame(means,index=var.index)
        fractions=pd.DataFrame(fractions,index=var.index)
        means.to_csv(self.directory/'cluster_mean_log1p_cpm.tsv.gz',sep='\t')
        fractions.to_csv(self.directory/'cluster_detection_fraction.tsv.gz',sep='\t')
        write_json(self.directory/'comparison.json',dict(context=context,method='two-sided Welch t-test; BH across all genes per contrast',
            expression='natural log1p(CPM)',detection='normalized expression > 0',
            interpretation='exploratory cell-level cluster markers, not biological replicate DE or HiCAT merge tests',filters=self.cfg))
        return means,fractions,top
