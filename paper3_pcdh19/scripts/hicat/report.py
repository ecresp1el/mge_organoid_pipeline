"""Three-page pilot report showing cluster counts, separation and stability.

All panels use observed pilot outputs. The Step 06 UMAP is a fixed display
only; neither it nor its old annotations determines a HiCAT split.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import linkage, dendrogram


class PilotReport:
    """Render a reviewable report with companion numerical tables."""

    def __init__(self,directory,config):
        """Create the figure directory and store display-only configuration."""
        self.directory=Path(directory)
        self.directory.mkdir()
        self.cfg=config

    @staticmethod
    def scatter(ax,coordinates,labels,title):
        """Draw categorical labels on fixed coordinates; no embedding is fitted."""
        names=sorted(set(labels))
        palette=plt.get_cmap('tab20')
        for i,name in enumerate(names):
            mask=np.asarray(labels)==name
            ax.scatter(coordinates[mask,0],coordinates[mask,1],s=6,alpha=.65,color=palette(i%20),label=name,rasterized=True)
        ax.set(title=title,xticks=[],yticks=[])
        if len(names)<=24:ax.legend(fontsize=6,markerscale=2,ncol=3,loc='best')

    @staticmethod
    def heatmap(ax,frame,title,cmap='viridis'):
        """Display a labeled numerical matrix; full values remain in saved TSVs."""
        im=ax.imshow(frame.to_numpy(dtype=float),aspect='auto',cmap=cmap,interpolation='nearest')
        ax.set_xticks(range(len(frame.columns)),frame.columns,rotation=90,fontsize=6)
        ax.set_yticks(range(len(frame.index)),frame.index,fontsize=6)
        ax.set_title(title,fontsize=10)
        return im

    def publish(self,obj,summaries,evidence,means,tree,stability,contingency):
        """Write three PDF pages and matching PNGs from all pilot cells.

        Returns
        -------
        list[str]
            Figure filenames. No full-data extrapolation or cell-type names
            appear in the report; cluster IDs are arbitrary within each run.
        """
        labels=obj.obs['hicat_baseline'].astype(str).to_numpy()
        coordinates=obj.obsm['X_umap_step06_display']
        pages=[]
        fig,axes=plt.subplots(2,2,figsize=(16,12),layout='constrained')
        stages=['root_graph_clusters','root_post_merge_clusters','recursive_leaves','final_clusters']
        for row in summaries.to_dict('records'):
            axes[0,0].plot(range(4),[row[x] for x in stages],'o-',label=row['candidate'])
        axes[0,0].set_xticks(range(4),['Root graph','Root merged','Recursive leaves','Final merged'])
        axes[0,0].set(title='A  Cluster count through the workflow',ylabel='Number of clusters')
        axes[0,0].legend(fontsize=8)
        axes[0,1].bar(summaries['candidate'],summaries['final_clusters'],color=['#4477aa','#66ccee','#cc6677'])
        axes[0,1].set(title='B  Pilot candidate counts',ylabel='Final clusters')
        for i,value in enumerate(summaries['final_clusters']):axes[0,1].text(i,value,str(value),ha='center',va='bottom')
        self.scatter(axes[1,0],coordinates,labels,'C  HiCAT cluster IDs on fixed Step 06 UMAP')
        self.scatter(axes[1,1],coordinates,obj.obs['technical_sample_id'].astype(str).to_numpy(),'D  Same cells and coordinates, colored by sample')
        fig.suptitle('HiCAT TECHNICAL PILOT — %s cells; not a full-dataset cluster count'%format(obj.n_obs,','),fontsize=16)
        pages.append(('01_cluster_counts_and_umap',fig))

        fig,axes=plt.subplots(2,2,figsize=(17,13),layout='constrained')
        sizes=pd.Series(labels).value_counts().sort_index()
        axes[0,0].bar(range(len(sizes)),sizes.values,color='#4477aa')
        axes[0,0].set_xticks(range(len(sizes)),sizes.index,rotation=90,fontsize=7)
        axes[0,0].set(title='E  Cluster sizes (all pilot cells)',ylabel='Cells')
        cross=pd.crosstab(obj.obs['technical_sample_id'],obj.obs['hicat_baseline'])
        im=self.heatmap(axes[0,1],cross,'F  Sample × cluster counts')
        fig.colorbar(im,ax=axes[0,1],label='Cells',shrink=.75)
        markers=[]
        for _,row in evidence.iterrows():
            for column in ['up_genes','down_genes']:
                markers.extend(list(row[column])[:3])
        markers=list(dict.fromkeys(markers))[:self.cfg['preview_max_marker_genes']]
        if markers:
            matrix=means[markers].T
            matrix=matrix.sub(matrix.mean(axis=1),axis=0).div(matrix.std(axis=1).replace(0,1),axis=0)
            symbols=obj.var['gene_symbol'].to_dict()
            matrix.index=[symbols.get(g,g) for g in matrix.index]
            im=self.heatmap(axes[1,0],matrix,'G  Data-derived separation genes (row z-score)','coolwarm')
            fig.colorbar(im,ax=axes[1,0],label='Z-score of cluster mean ln(1+CPM)',shrink=.75)
        else:
            axes[1,0].text(.5,.5,'No between-cluster genes to display',ha='center');axes[1,0].set_axis_off()
        if len(sizes)>1:
            matrix=pd.DataFrame(np.nan,index=sizes.index,columns=sizes.index)
            for pair,row in evidence.iterrows():matrix.loc[pair[0],pair[1]]=matrix.loc[pair[1],pair[0]]=row['score']
            # Cap the COLOR SCALE only; uncapped scores (including inf) remain in tables.
            capped=matrix.clip(upper=2*self.cfg['thresholds']['score_thresh'])
            im=self.heatmap(axes[1,1],capped,'H  Final pairwise DE scores (display capped at 300)')
            fig.colorbar(im,ax=axes[1,1],label='DE score; baseline threshold 150',shrink=.75)
        else:
            axes[1,1].text(.5,.5,'One cluster: no between-cluster comparison',ha='center');axes[1,1].set_axis_off()
        fig.suptitle('Separation evidence — gene labels describe evidence, not cell-type annotation',fontsize=16)
        pages.append(('02_sizes_composition_and_separation',fig))

        fig,axes=plt.subplots(2,2,figsize=(17,13),layout='constrained')
        positions={row.node:(row.depth,-i) for i,row in enumerate(tree.itertuples())}
        for row in tree.itertuples():
            x,y=positions[row.node]
            if row.parent in positions:
                xp,yp=positions[row.parent];axes[0,0].plot([xp,x],[yp,y],color='#999999',lw=.8)
            axes[0,0].scatter(x,y,s=18,color='#cc6677' if row.terminal_reason=='split' else '#4477aa')
            axes[0,0].text(x+.05,y,'%s: n=%s, K=%s'%(row.node,int(row.cells),int(row.post_merge_clusters)),fontsize=6)
        axes[0,0].set(title='I  Recorded recursive split tree (not a cell-type taxonomy)',xlabel='Recursion depth',yticks=[])
        im=self.heatmap(axes[0,1],contingency,'J  Baseline × repeat-seed membership overlap')
        fig.colorbar(im,ax=axes[0,1],label='Shared cells',shrink=.75)
        repeat=stability[stability['candidate']=='seed_repeat']
        axes[1,0].bar(range(len(repeat)),repeat['best_jaccard'],color='#66ccee')
        axes[1,0].set_xticks(range(len(repeat)),repeat['baseline_cluster'],rotation=90,fontsize=7)
        axes[1,0].set(title='K  Best-match cluster overlap across seeds',ylabel='Jaccard overlap',ylim=(0,1.05))
        if len(means)>1:
            linkage_matrix=linkage(means.to_numpy(),method='average',metric='euclidean')
            dendrogram(linkage_matrix,labels=means.index.to_list(),ax=axes[1,1],leaf_rotation=90,leaf_font_size=7)
            axes[1,1].set(title='L  Cluster-mean expression similarity',ylabel='Euclidean distance; average linkage')
        else:
            axes[1,1].text(.5,.5,'One cluster: no dendrogram',ha='center');axes[1,1].set_axis_off()
        fig.suptitle('Pilot stability and hierarchy — neither establishes a full-data K',fontsize=16)
        pages.append(('03_recursion_and_stability',fig))
        names=['hicat_pilot_review.pdf']
        with PdfPages(self.directory/names[0]) as pdf:
            for name,fig in pages:
                pdf.savefig(fig,bbox_inches='tight')
                fig.savefig(self.directory/(name+'.png'),dpi=self.cfg['report_dpi'],bbox_inches='tight')
                names.append(name+'.png');plt.close(fig)
        return names
