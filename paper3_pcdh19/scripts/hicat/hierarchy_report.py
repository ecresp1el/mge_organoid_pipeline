"""Paginated hierarchy, dendrogram and per-cluster top-gene review figures.

Every plotted matrix, linkage and label order is persisted. Display UMAP comes
from Step 06; expression dendrograms summarize cluster means and are distinct
from the recorded algorithmic split tree. All pages are saved as PDF and PNG.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import linkage, dendrogram
from .provenance import write_json


class HierarchyReport:
    """Build a readable review book without silently truncating cluster coverage."""

    def __init__(self, directory, dpi=160):
        """Open the report and prepare a page manifest for every exported figure."""
        self.directory=Path(directory)
        self.directory.mkdir()
        self.pdf=PdfPages(self.directory/'coarse_fine_review.pdf')
        self.dpi=dpi
        self.pages=[]

    def save(self, figure, name, description):
        """Save one PDF page and a matching PNG, recording its interpretation."""
        figure.tight_layout(rect=(0,0,1,0.95))
        self.pdf.savefig(figure,bbox_inches='tight')
        filename='%03d_%s.png'%(len(self.pages)+1,name)
        figure.savefig(self.directory/filename,dpi=self.dpi,bbox_inches='tight')
        self.pages.append(dict(page=len(self.pages)+1,png=filename,description=description))
        plt.close(figure)

    def overview(self, raw, parents, summary, agreement):
        """Show both levels on fixed coordinates and their exact parent-child counts."""
        fig,axes=plt.subplots(2,2,figsize=(16,12))
        fig.suptitle('Coarse → fine | %s cells across 12 samples | no full-data K'%format(raw.n_obs,','),fontsize=17)
        xy=raw.obsm['X_umap_step06_display']
        for ax,key,title in zip(axes[0],['hicat_coarse_baseline','hicat_fine_baseline'],['Coarse groups','Fine groups']):
            labels=raw.obs[key].astype(str).to_numpy()
            groups=sorted(set(labels));palette=plt.get_cmap('tab20')
            for i,group in enumerate(groups):
                selected=labels==group
                ax.scatter(xy[selected,0],xy[selected,1],s=2,color=palette(i%20),alpha=.6,rasterized=True)
                position=np.median(xy[selected],axis=0)
                ax.text(*position,group,fontsize=7,bbox=dict(facecolor='white',alpha=.7,edgecolor='none'))
            ax.set_title(title+' on existing Step 06 UMAP (display only)')
            ax.set_xticks([]);ax.set_yticks([])
        counts=parents.groupby('coarse_cluster').agg(cells=('cells','sum'),fine_clusters=('fine_cluster','size'))
        counts.fine_clusters.plot.bar(ax=axes[1,0],color='steelblue')
        axes[1,0].set_title('Number of fine clusters within each coarse parent')
        axes[1,0].set_ylabel('Fine clusters')
        for i,value in enumerate(counts.fine_clusters):axes[1,0].text(i,value,str(value),ha='center',va='bottom')
        axes[1,1].axis('off')
        text='Candidate counts\n'+summary[['candidate','coarse_clusters','fine_clusters','fine_unresolved_pairs']].to_string(index=False)
        text+='\n\nRepeat-seed agreement\n'+agreement.to_string(index=False)
        text+='\n\nNo batch correction or biological annotation.\nCounts are conditional on this subset and thresholds.'
        axes[1,1].text(0,.98,text,va='top',family='monospace',fontsize=9)
        self.save(fig,'overview','Coarse and fine assignments, nested counts and seed sensitivity')
        for start in range(0,len(counts),6):
            section=counts.index[start:start+6]
            fig,ax=plt.subplots(figsize=(14,max(5,len(parents[parents.coarse_cluster.isin(section)])*.32)))
            cursor=0
            for parent in section:
                children=parents[parents.coarse_cluster==parent]
                positions=np.arange(cursor,cursor+len(children))
                center=positions.mean()
                ax.text(0,center,parent+' (n=%s)'%format(int(counts.loc[parent,'cells']),','),ha='right',va='center')
                for y,child in zip(positions,children.itertuples()):
                    ax.plot([.03,.95],[center,y],color='gray',lw=1)
                    ax.text(1,y,child.fine_cluster+' (n=%s)'%format(child.cells,','),va='center')
                cursor+=len(children)+1
            ax.set_xlim(-.65,1.8);ax.set_ylim(-1,cursor);ax.invert_yaxis();ax.axis('off')
            fig.suptitle('Saved coarse → fine membership hierarchy (not an expression dendrogram)')
            self.save(fig,'parent_children_%02d'%start,'Every fine cluster linked to its fitted coarse parent')

    def expression(self, context, directory, means, fractions, top):
        """Save a mean-expression dendrogram and one top-20 page per cluster.

        Dendrograms use Euclidean distance on unscaled mean ln(1+CPM) across
        the union of reported marker genes. Heatmaps show row z-scores; adjacent
        detection heatmaps show percentages. All exact plot matrices are saved.
        A single cluster or no qualified genes produces an explicit notice.
        """
        directory=Path(directory)
        genes=list(dict.fromkeys(top.gene_id.astype(str)))
        clusters=list(means.columns)
        order=clusters
        if len(clusters)>1 and genes:
            matrix=means.loc[genes].T
            tree=linkage(matrix.to_numpy(),method='average',metric='euclidean')
            np.save(directory/'expression_linkage.npy',tree)
            matrix.to_csv(directory/'dendrogram_input_cluster_means.tsv',sep='\t')
            pd.Series(clusters,name='cluster').to_csv(directory/'linkage_input_labels.tsv',sep='\t',index=False)
            fig,ax=plt.subplots(figsize=(max(10,len(clusters)*.35),6))
            drawing=dendrogram(tree,labels=clusters,ax=ax,leaf_rotation=90)
            order=drawing['ivl']
            ax.set_ylabel('Euclidean distance; average linkage')
            ax.set_title('Mean expression across union of qualified top-20 genes')
            fig.suptitle(context+' | expression dendrogram')
            self.save(fig,context+'_dendrogram','Cluster-mean similarity; not a confidence-supported taxonomy')
        else:
            write_json(directory/'dendrogram_absent.json',dict(reason='one cluster or no qualifying marker genes'))
        pd.Series(order,name='cluster').to_csv(directory/'heatmap_cluster_order.tsv',sep='\t',index=False)
        for cluster in clusters:
            selected=top[top.cluster==cluster]
            genes=selected.gene_id.astype(str).tolist()
            if not genes:
                fig,ax=plt.subplots(figsize=(11,4));ax.axis('off')
                ax.text(.03,.6,'No qualifying positive marker genes, or no comparison group.\nSee marker_coverage.tsv for this cluster.',fontsize=13)
                fig.suptitle(context+' | '+cluster)
                self.save(fig,context+'_'+cluster+'_no_markers','Explicit absent-marker result')
                continue
            values=means.loc[genes,order]
            sd=values.std(axis=1,ddof=0).replace(0,1)
            z=values.sub(values.mean(axis=1),axis=0).div(sd,axis=0)
            detection=fractions.loc[genes,order]*100
            z.to_csv(directory/(cluster+'_heatmap_zscore.tsv'),sep='\t')
            detection.to_csv(directory/(cluster+'_heatmap_percent_detected.tsv'),sep='\t')
            labels=[str(i+1)+'. '+str(row.gene_symbol)+' | '+str(row.gene_id) for i,row in enumerate(selected.itertuples())]
            fig,axes=plt.subplots(1,2,figsize=(max(16,len(order)*.45+9),max(6,len(genes)*.32+2)),gridspec_kw={'width_ratios':[1,1]})
            im=axes[0].imshow(z.to_numpy(),aspect='auto',cmap='coolwarm',vmin=-2,vmax=2)
            axes[0].set_yticks(range(len(genes)));axes[0].set_yticklabels(labels,fontsize=8)
            axes[0].set_title('Mean ln(1+CPM), row z-score')
            fig.colorbar(im,ax=axes[0],shrink=.6)
            im=axes[1].imshow(detection.to_numpy(),aspect='auto',cmap='viridis',vmin=0,vmax=100)
            axes[1].set_yticks(range(len(genes)));axes[1].set_yticklabels(selected.gene_symbol,fontsize=8)
            axes[1].set_title('Cells with expression >0 (%)')
            fig.colorbar(im,ax=axes[1],shrink=.6)
            for ax in axes:
                ax.set_xticks(range(len(order)));ax.set_xticklabels(order,rotation=90,fontsize=8)
            fig.suptitle(context+' | '+cluster+' | top %d qualified genes (ranked by Welch t)'%len(genes))
            self.save(fig,context+'_'+cluster+'_top20','Ranked marker mean-expression and detection heatmaps')

    def composition(self, raw, level, directory):
        """Save counts/fractions and show sample contributions at either level."""
        counts=pd.crosstab(raw.obs['technical_sample_id'],raw.obs['hicat_'+level+'_baseline'])
        counts.to_csv(Path(directory)/'sample_cluster_counts.tsv',sep='\t')
        fractions=counts.div(counts.sum(axis=0),axis=1)
        fractions.to_csv(Path(directory)/'sample_within_cluster_fractions.tsv',sep='\t')
        fig,axes=plt.subplots(1,2,figsize=(max(14,len(counts.columns)*.45),6))
        for ax,values,title in zip(axes,[counts,fractions*100],['Cell counts','Within-cluster sample contributions (%)']):
            im=ax.imshow(values,aspect='auto',cmap='viridis')
            ax.set_xticks(range(len(counts.columns)));ax.set_xticklabels(counts.columns,rotation=90,fontsize=8)
            ax.set_yticks(range(len(counts)));ax.set_yticklabels(counts.index,fontsize=8)
            ax.set_title(title);fig.colorbar(im,ax=ax,shrink=.7)
        fig.suptitle(level+' | sample composition (no batch correction)')
        self.save(fig,level+'_sample_composition','Sample contribution counts and percentages')

    def close(self):
        """Finalize the PDF and write a complete PNG/page description index."""
        self.pdf.close()
        pd.DataFrame(self.pages).to_csv(self.directory/'figure_index.tsv',sep='\t',index=False)
