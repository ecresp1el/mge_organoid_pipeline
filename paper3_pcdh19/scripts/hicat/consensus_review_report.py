"""Numbered consensus figures using the existing HiCAT report publisher."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import linkage,dendrogram
from .hierarchy_report import HierarchyReport


class ConsensusReviewReport(HierarchyReport):
    """Render saved R=98 labels with exact source tables and numbered pages."""

    def __init__(self, output, cfg):
        """Use the hierarchy report's save/close convention and PNG numbering."""
        self.output=Path(output)
        self.directory=self.output/'figures'
        self.directory.mkdir()
        self.pdf=PdfPages(self.directory/'hicat_consensus_review.pdf')
        self.dpi=cfg['dpi'];self.pages=[];self.cfg=cfg

    def save_source(self, fig, name, description, sources):
        """Publish the numbered page and record its exact plotted inputs."""
        self.save(fig,name,description)
        self.pages[-1]['source_tables']=';'.join(sources)

    def render(self, data):
        """Generate cluster/sample, marker and expression-dendrogram pages."""
        tables=self.output/'tables'
        counts=pd.crosstab(data.assignments.diagnostic_cluster,data.universe.technical_sample_id)
        fractions=counts.div(counts.sum(axis=1),axis=0)
        counts.to_csv(tables/'sample_cluster_counts.tsv',sep='\t')
        fractions.to_csv(tables/'sample_cluster_within_cluster_fraction.tsv',sep='\t')
        counts.div(counts.sum(axis=0),axis=1).to_csv(tables/'sample_cluster_within_sample_fraction.tsv',sep='\t')
        sizes=counts.sum(axis=1);sizes.rename('cells').to_csv(tables/'cluster_sizes.tsv',sep='\t')
        fig,axes=plt.subplots(1,2,figsize=(16,max(7,len(counts)*.25)))
        y=np.arange(len(counts));axes[0].barh(y,sizes.to_numpy())
        axes[0].set_yticks(y);axes[0].set_yticklabels(counts.index);axes[0].invert_yaxis()
        axes[0].set_title('A | Consensus cluster sizes');axes[0].set_xlabel('Cells')
        im=axes[1].imshow(fractions,aspect='auto',vmin=0,vmax=1)
        axes[1].set_xticks(range(len(counts.columns)));axes[1].set_xticklabels(counts.columns,rotation=90)
        axes[1].set_yticks(y);axes[1].set_yticklabels(counts.index)
        axes[1].set_title('B | Within-cluster sample contributions')
        fig.colorbar(im,ax=axes[1],label='Fraction of cluster cells')
        fig.suptitle('HiCAT 03 | %d-iteration consensus | %s cells, %d clusters | IN_REVIEW'%
                     (self.cfg['iterations'],format(len(data.assignments),','),len(counts)))
        self.save_source(fig,'consensus_composition','Saved final-DE partition; labels are not biological annotations',
                         ['tables/cluster_sizes.tsv','tables/sample_cluster_within_cluster_fraction.tsv'])
        markers=[];gene_ids=[]
        for symbol in self.cfg['marker_panel']:
            matches=[g for g in data.means.columns if data.symbols.get(g)==symbol]
            if len(matches)==1:markers.append(symbol);gene_ids.append(matches[0])
        matrix=data.means.loc[:,gene_ids].copy();matrix.columns=markers
        z=(matrix-matrix.mean())/matrix.std(ddof=0).replace(0,1)
        matrix.to_csv(tables/'marker_mean_log1p_cpm.tsv',sep='\t')
        z.to_csv(tables/'marker_heatmap_zscores.tsv',sep='\t')
        pd.DataFrame(dict(gene_id=gene_ids,gene_symbol=markers)).to_csv(tables/'marker_panel_mapping.tsv',sep='\t',index=False)
        fig,ax=plt.subplots(figsize=(14,max(7,len(matrix)*.26)))
        im=ax.imshow(z,aspect='auto',cmap='RdBu_r',vmin=-2,vmax=2)
        ax.set_xticks(range(len(markers)));ax.set_xticklabels(markers,rotation=90)
        ax.set_yticks(range(len(matrix)));ax.set_yticklabels(matrix.index)
        ax.set_title('C | Marker expression across the saved consensus clusters')
        fig.suptitle('HiCAT 03 | R=%d | display-only marker panel; no assigned cell types'%self.cfg['iterations'])
        fig.colorbar(im,ax=ax,label='Gene-wise z score of cluster mean ln(1+CPM)')
        self.save_source(fig,'consensus_marker_expression','Gene-wise z score across cluster means; not a ranked DEG panel',
                         ['tables/marker_heatmap_zscores.tsv','tables/marker_mean_log1p_cpm.tsv','tables/marker_panel_mapping.tsv'])
        if len(data.means)>1:
            tree=linkage(data.means.to_numpy(),method='average',metric='correlation')
            np.save(tables/'expression_linkage.npy',tree)
            pd.Series(data.means.index,name='cluster').to_csv(tables/'expression_linkage_labels.tsv',sep='\t',index=False)
            fig,ax=plt.subplots(figsize=(14,6));dendrogram(tree,labels=data.means.index.to_list(),ax=ax,leaf_rotation=90)
            ax.set_title('D | Cluster-mean expression dendrogram (all %d genes)'%data.means.shape[1])
            ax.set_ylabel('Average-linkage correlation distance')
            fig.suptitle('HiCAT 03 | expression display summary; not the recursive fitting tree')
            self.save_source(fig,'consensus_expression_dendrogram','Average linkage on correlation distance of all-gene cluster means',
                             ['tables/cluster_mean_log1p_cpm.tsv.gz','tables/expression_linkage.npy','tables/expression_linkage_labels.tsv'])
        self.close()
