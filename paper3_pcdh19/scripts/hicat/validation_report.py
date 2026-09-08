"""Review books and auditable figures for the E14.5 dissected mouse MGE pilot.

The main PDF contains panels A–H across readable pages. The detailed PDF adds
coarse diagnostics, every fine cluster's seven-panel review, full-data gene and
program views, and paginated membership overlays. Figure indices list numeric
source tables. Full-data coordinates are inherited from Step 06; neither labels
nor fitted partitions are projected to unselected cells.
"""
import json
import textwrap
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist, squareform
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from .validation_review import review_clusters, DEVELOPMENT, REGION

PALETTE = ['#4878a8', '#e5a44b', '#64a578', '#aa6aab']
PHASE_COLORS = ['#b9bec8', '#559dc4', '#b25b88']
TITLE = 'Biological and stability validation of the exploratory HiCAT hierarchy'
CONTEXT = 'Dissected E14.5 mouse MGE | 12,000 pilot cells | hypotheses only'


def _finite_limits(values, percentile=(2, 98), symmetric=False):
    """Return explicit robust limits; missing values remain visually distinct."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return 0., 1.
    lo, hi = np.percentile(values, percentile)
    if symmetric:
        hi = max(abs(lo), abs(hi), 1e-6)
        return -hi, hi
    if hi <= lo:
        hi = lo + 1.
    return float(lo), float(hi)


def _heatmap(ax, frame, title, cmap='viridis', limits=None, colorbar=True, fontsize=7):
    """Draw a labeled matrix with unavailable assay/program values in gray."""
    cm = plt.get_cmap(cmap).copy()
    cm.set_bad('#d6d6d6')
    kwargs = {} if limits is None else dict(vmin=limits[0], vmax=limits[1])
    im = ax.imshow(np.ma.masked_invalid(frame.to_numpy(dtype=float)), aspect='auto', cmap=cm, **kwargs)
    ax.set_xticks(np.arange(len(frame.columns)))
    ax.set_xticklabels(frame.columns, rotation=90, fontsize=fontsize)
    ax.set_yticks(np.arange(len(frame.index)))
    ax.set_yticklabels(frame.index, fontsize=fontsize)
    ax.set_title(title, fontsize=11)
    if colorbar:
        ax.figure.colorbar(im, ax=ax, fraction=.025, pad=.015)
    return im


def _stacked(ax, frame, title, colors=None, legend=True):
    """Plot a complete composition matrix without dropping zero categories."""
    bottom = np.zeros(len(frame))
    colors = colors or [plt.get_cmap('tab20')(i % 20) for i in range(len(frame.columns))]
    for i, column in enumerate(frame):
        values = frame[column].to_numpy(dtype=float)
        ax.bar(np.arange(len(frame)), values, bottom=bottom, color=colors[i % len(colors)], label=column, width=.86)
        bottom += values
    ax.set_xticks(np.arange(len(frame)))
    ax.set_xticklabels(frame.index, rotation=90, fontsize=7)
    ax.set_ylim(0, 1)
    ax.set_ylabel('Fraction of cells')
    ax.set_title(title, fontsize=11)
    if legend:
        ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1), fontsize=7, frameon=False)


def _plain_text(ax, text, size=9, width=105):
    """Wrap prose at intentional paragraph boundaries for compact report pages."""
    ax.axis('off')
    wrapped = '\n'.join(textwrap.fill(line, width=width, subsequent_indent='  ') if line else '' for line in text.split('\n'))
    ax.text(0, 1, wrapped, va='top', ha='left', fontsize=size, transform=ax.transAxes)


class ValidationReport:
    """Own figure lifecycle, separate main/detail books, and plot provenance."""

    def __init__(self, run_root, dpi=130):
        """Create output folders and open two deterministic review documents."""
        self.root = Path(run_root)
        self.out = self.root / 'outputs'
        self.figures = self.out / 'figures'
        self.sources = self.out / 'tables' / 'figure_sources'
        self.figures.mkdir(parents=True, exist_ok=True)
        self.sources.mkdir(parents=True, exist_ok=True)
        self.main = PdfPages(self.figures / 'hicat_hierarchy_validation_main.pdf')
        self.detail = PdfPages(self.figures / 'hicat_hierarchy_validation_detailed.pdf')
        self.dpi = dpi
        self.index = []
        self.main_pages = 0
        self.detail_pages = 0

    def table(self, name, frame, index=True):
        """Save exact numeric plotted values and return a package-relative path."""
        path = self.sources / (name + '.tsv')
        frame.to_csv(path, sep='\t', index=index)
        return str(path.relative_to(self.out))

    def save(self, fig, name, description, sources, main=False):
        """Save each panel/page as PNG, PDF and detailed-book page; index sources."""
        fig.subplots_adjust(top=.92, bottom=.16, hspace=.65, wspace=.40)
        fig.savefig(self.figures / (name + '.png'), dpi=self.dpi, bbox_inches='tight')
        fig.savefig(self.figures / (name + '.pdf'), bbox_inches='tight')
        self.detail.savefig(fig, bbox_inches='tight')
        self.detail_pages += 1
        main_page = ''
        if main:
            self.main.savefig(fig, bbox_inches='tight')
            self.main_pages += 1
            main_page = self.main_pages
        self.index.append(dict(name=name, description=description, main_page=main_page,
            detailed_page=self.detail_pages, png='figures/' + name + '.png', pdf='figures/' + name + '.pdf',
            source_tables=';'.join(str(s) for s in sources)))
        plt.close(fig)
        print('REPORT_PAGE_COMPLETE %s main=%s detail=%s' % (name, main_page, self.detail_pages), flush=True)

    def close(self):
        """Close document streams and persist an exact figure/source inventory."""
        self.main.close()
        self.detail.close()
        pd.DataFrame(self.index).to_csv(self.figures/'figure_index.tsv', sep='\t', index=False)


def _dotplot(ax, mean, fraction, genes, title, group_boundaries=None):
    """Canonical markers: color is mean ln1pCPM, area is percent detected."""
    available = [g for g in genes if g in mean.columns]
    values = mean[available]
    detection = fraction[available]
    x, y = np.meshgrid(np.arange(len(available)), np.arange(len(values)))
    limits = _finite_limits(values.to_numpy(), (0, 98))
    dots = ax.scatter(x.ravel(), y.ravel(), s=detection.to_numpy().ravel()*80,
        c=values.to_numpy().ravel(), cmap='viridis', vmin=limits[0], vmax=limits[1], edgecolors='none')
    ax.set_xticks(np.arange(len(available)))
    ax.set_xticklabels(available, rotation=90, fontsize=7)
    ax.set_yticks(np.arange(len(values)))
    ax.set_yticklabels(values.index, fontsize=7)
    ax.set_ylim(len(values)-.5, -.5)
    ax.set_xlim(-.5, len(available)-.5)
    ax.set_title(title, fontsize=11)
    for boundary in ([] if group_boundaries is None else group_boundaries):
        if boundary < len(available):
            ax.axvline(boundary-.5, color='#aaaaaa', lw=.5)
    colorbar = ax.figure.colorbar(dots, ax=ax, fraction=.020, pad=.008)
    colorbar.set_label('Mean ln(1+CPM)', fontsize=8)
    for fraction_value in [.25, .5, 1.]:
        ax.scatter([], [], s=fraction_value*80, color='#888888', label='%d%%' % (100*fraction_value))
    ax.legend(loc='upper center', bbox_to_anchor=(.5, -.27), ncol=3, fontsize=8, frameon=False, title='Detection (dot area)')
    return values, detection


def _tree(ax, fine, parent_order=None):
    """Plot the actual saved parent membership tree, without inferred topology."""
    order = []
    parent_order = parent_order or sorted(fine.coarse_parent.unique())
    cursor = 0
    for pindex, parent in enumerate(parent_order):
        children = fine[fine.coarse_parent == parent].sort_index()
        positions = np.arange(cursor, cursor+len(children))
        center = positions.mean()
        color = PALETTE[pindex % len(PALETTE)]
        ax.text(-.04, center, '%s\nn=%s' % (parent, format(int(children.n_cells.sum()), ',')), ha='right', va='center', fontsize=9)
        for y, (cluster, row) in zip(positions, children.iterrows()):
            ax.plot([0, .9], [center, y], color=color, lw=1)
            ax.text(.95, y, '%s  n=%s' % (cluster, row.n_cells), va='center', fontsize=7)
            order.append(cluster)
        cursor += len(children)
    ax.set_xlim(-.6, 2.1)
    ax.set_ylim(cursor-.5, -.5)
    ax.axis('off')
    ax.set_title('A. Actual coarse → fine membership\nNot an expression dendrogram', fontsize=11)
    return order


def _overlap_order(frame):
    """Maximize correspondence within retained candidate-parent blocks only."""
    rows = list(frame.index)
    def column_key(column):
        """Keep actual candidate parents together, then rank by baseline overlap."""
        parent = str(column).split('.')[0]
        best = int(np.argmax(frame[column].to_numpy())) if frame[column].sum() else len(rows)
        return parent, best, str(column)
    return frame.loc[rows, sorted(frame.columns, key=column_key)]


def _relationships(book, review, config, canonical_genes):
    """Compare canonical expression similarity and imposed parents; flag issues."""
    fine = review['fine']
    genes = [g for g in canonical_genes if g in review['gene_means']['fine']]
    matrix = review['gene_means']['fine'][genes].loc[fine.index]
    tree = linkage(matrix.to_numpy(), method='average', metric='euclidean')
    distance = pd.DataFrame(squareform(pdist(matrix)), index=matrix.index, columns=matrix.index)
    source = book.table('canonical_expression_dendrogram_input', matrix)
    tree_source = book.table('canonical_expression_linkage', pd.DataFrame(tree, columns=['left','right','distance','leaves']))
    distance_source = book.table('canonical_expression_distances', distance)
    # The thresholds describe tails of this distance distribution, not p-values.
    upper = distance.to_numpy()[np.triu_indices(len(distance), 1)]
    close_cutoff, far_cutoff = np.percentile(upper, [10, 90])
    flags = []
    for i, a in enumerate(distance.index):
        for b in distance.index[i+1:]:
            same = a.split('.')[0] == b.split('.')[0]
            value = distance.loc[a,b]
            if (not same and value <= close_cutoff) or (same and value >= far_cutoff):
                flags.append(dict(candidate='baseline',cluster_a=a,cluster_b=b,
                    flag='similar_across_parents' if not same else 'distant_siblings',
                    canonical_expression_distance=value, threshold=close_cutoff if not same else far_cutoff,
                    DE_score=np.nan, interpretation='Descriptive canonical-expression distance tail; not an independent taxonomy test'))
    for candidate in ['baseline','seed_repeat']:
        path = Path(config['baseline_run']) / candidate / 'global_fine_pairwise_evidence.tsv'
        if path.exists():
            audit = pd.read_csv(path, sep='\t', index_col=[0,1])
            saved = audit[['score','num','same_coarse_parent','meets_declared_separation']].copy()
            saved.index.names = ['cluster_a','cluster_b']
            saved.to_csv(book.out/'tables'/('%s_original_final_pair_audit.tsv' % candidate),sep='\t')
            failed = saved[~saved.meets_declared_separation.astype(bool)]
            for (a,b), row in failed.iterrows():
                flags.append(dict(candidate=candidate,cluster_a=a,cluster_b=b,
                    flag='cross_parent_DE_failure' if not row.same_coarse_parent else 'within_parent_DE_failure',
                    canonical_expression_distance=np.nan, threshold=150, DE_score=row.score,
                    interpretation='Original final audit failed declared score or gene-count requirement; baseline labels preserved'))
    flag_table = pd.DataFrame(flags, columns=['candidate','cluster_a','cluster_b','flag','canonical_expression_distance','threshold','DE_score','interpretation'])
    flag_table.to_csv(book.out/'annotation_review'/'hierarchy_relationship_flags.tsv',sep='\t',index=False)
    fig = plt.figure(figsize=(18,15))
    gs = fig.add_gridspec(1,3,width_ratios=[1.05,1,.72])
    ax = fig.add_subplot(gs[0]); _tree(ax,fine)
    ax = fig.add_subplot(gs[1]); drawing = dendrogram(tree, labels=matrix.index.tolist(), orientation='left', ax=ax, leaf_font_size=7)
    ax.set_title('Expression-similarity dendrogram\nCanonical genes (including cycle/QC), mean ln1pCPM',fontsize=11)
    ax.set_xlabel('Euclidean distance; average linkage')
    program_names = ['rg_stemness','ipc_neurogenic','G2M','neuroblast','MGE','LGE','later_neuronal_maturation','endothelial','erythroid']
    values = review['program_means']['fine'].loc[drawing['ivl'],program_names]
    # Matplotlib dendrogram leaves bottom-to-top; reverse the matrix rows so
    # the top heatmap row corresponds to the highest dendrogram leaf.
    values = values.iloc[::-1]
    z = values.sub(review['program_means']['fine'][program_names].mean()).div(review['program_means']['fine'][program_names].std(ddof=0).replace(0,np.nan))
    ax = fig.add_subplot(gs[2]); _heatmap(ax,z,'Canonical program sidebar\nColumn z-score',cmap='coolwarm',limits=(-2,2),fontsize=6)
    fig.suptitle(CONTEXT + '\nActual nesting versus canonical expression relationships',fontsize=15)
    sidebar_source=book.table('hierarchy_dendrogram_sidebar_zscore',z)
    book.save(fig,'hierarchy_vs_expression_dendrogram','Actual memberships compared with independent canonical gene means and program sidebars',[source,tree_source,distance_source,sidebar_source,'annotation_review/hierarchy_relationship_flags.tsv'])
    return flag_table


def _summary_pages(book, review, program_cfg, metadata, metrics):
    """Generate main panels A–G and complete coarse and fine matrix summaries."""
    fine = review['fine']
    all_genes = list(dict.fromkeys(g for group in program_cfg['dotplot_groups'].values() for g in group))
    order = [p['name'] for p in program_cfg['programs']]
    mean = review['program_means']['fine'].reindex(columns=order)
    sd = mean.std(axis=0,ddof=0).replace(0,np.nan)
    z = mean.sub(mean.mean()).div(sd)
    fig=plt.figure(figsize=(21,15));gs=fig.add_gridspec(1,2,width_ratios=[.95,1.75])
    _tree(fig.add_subplot(gs[0]),fine)
    _heatmap(fig.add_subplot(gs[1]),z,'B. Canonical programs by current fine cluster\nColumn z-score of cluster means; gray = unavailable',cmap='coolwarm',limits=(-2,2))
    fig.suptitle(TITLE+'\n'+CONTEXT,fontsize=16)
    hierarchy_source=book.table('baseline_parent_child_counts',fine[['coarse_parent','n_cells','fraction_of_pilot']])
    z_source=book.table('fine_program_heatmap_zscore',z)
    book.save(fig,'main_AB_hierarchy_programs','A actual nesting; B independent canonical program scores',[hierarchy_source,z_source,'marker_programs/fine_cluster_program_means.tsv'],main=True)
    for level in ['coarse','fine']:
        means=review['gene_means'][level]; fractions=review['gene_fractions'][level]
        fig,ax=plt.subplots(figsize=(23, max(6,len(means)*.27+4)))
        plotted_mean,plotted_fraction=_dotplot(ax,means,fractions,all_genes,('C. ' if level=='fine' else '')+'Canonical developmental / regional / non-neural marker expression', group_boundaries=np.cumsum([len(g) for g in program_cfg['dotplot_groups'].values()])[:-1])
        offset=0
        for group_name,group_genes in program_cfg['dotplot_groups'].items():
            ax.text(offset+len(group_genes)/2-.5,-1.3,group_name,fontsize=7,ha='center',rotation=30,clip_on=False)
            offset+=len(group_genes)
        fig.suptitle(CONTEXT+'\n'+level.capitalize()+' canonical marker dotplot (independent of top-DEG rankings)',fontsize=15)
        sources=[book.table(level+'_canonical_marker_means',plotted_mean),book.table(level+'_canonical_marker_detection_fraction',plotted_fraction)]
        book.save(fig,level+'_canonical_marker_dotplot','Color mean ln1pCPM; dot area percentage of cells with expression >0',sources,main=(level=='fine'))
        for statistic,values in [('mean',review['program_means'][level])]:
            fig,ax=plt.subplots(figsize=(17,max(5,len(values)*.28+3)))
            _heatmap(ax,values,'Canonical program absolute mean score (signal minus controls)',cmap='coolwarm',limits=_finite_limits(values,(2,98),True))
            fig.suptitle(CONTEXT+'\n'+level.capitalize()+' program scores; gray values are unavailable',fontsize=14)
            book.save(fig,level+'_program_absolute_'+statistic,'Absolute scores retain their original units and signs',['marker_programs/%s_cluster_program_means.tsv'%level])
        if 'positive' in metadata:
            # Join original assignments supplied by scorer caller, preserving cells.
            assignment=book.pilot.obs['hicat_'+level+'_baseline'].astype(str)
            positive=metadata['positive'].groupby(assignment,observed=True).mean()
            positive_source=book.table(level+'_program_positive_fraction',positive)
            positive.to_csv(book.out/'marker_programs'/(level+'_cluster_program_positive_fractions.tsv'),sep='\t',index_label='cluster')
            fig,ax=plt.subplots(figsize=(17,max(5,len(positive)*.28+3)))
            _heatmap(ax,positive,'Fraction positive: score threshold plus ≥2 detected signature genes',limits=(0,1))
            fig.suptitle(CONTEXT+'\n'+level.capitalize()+' canonical program-positive fractions',fontsize=14)
            book.save(fig,level+'_program_positive_fraction','Program-positive fraction; threshold and coverage rules in score model',[positive_source])
    fig,axes=plt.subplots(2,1,figsize=(19,14))
    _stacked(axes[0],review['phase_fractions']['fine'],'D. Fine-cluster cell-cycle composition',PHASE_COLORS)
    sample=metrics['sample_fine']['enrichment']
    if set(fine.index).issubset(sample.columns):sample=sample.T
    _heatmap(axes[1],sample,'E. Sample contribution / balanced expectation (1/12)',cmap='coolwarm',limits=(0,3))
    fig.suptitle(TITLE+'\n'+CONTEXT,fontsize=16)
    source=book.table('fine_sample_enrichment',sample)
    book.save(fig,'main_DE_cycle_sample','Cycle composition and observed/expected sample representation',['cell_cycle/fine_phase_fractions.tsv',source],main=True)
    fig,axes=plt.subplots(1,2,figsize=(22,15))
    sources=[]
    for ax,key,title in zip(axes,['seed_fine','allen_fine'],['F. Second-seed fine overlap','G. Allen-reference fine overlap']):
        frame=_overlap_order(metrics[key]['row_fractions'])
        _heatmap(ax,frame,title+'\nRow = baseline; color = retained fraction',limits=(0,1),fontsize=6)
        ax.set_xlabel('Actual candidate IDs, parent blocks retained')
        sources.append(book.table(key+'_row_fraction_plot_order',frame))
    fig.suptitle(TITLE+'\nSame 12,000 cells; no matching-ID renaming; fixed parents only for Allen condition',fontsize=15)
    book.save(fig,'main_FG_seed_Allen_overlap','Baseline-to-repeat and baseline-to-Allen overlap fractions',sources,main=True)
    return all_genes


def _extended_diagnostics(book,pilot,scores,review,metrics):
    """Expose distributions, reciprocal sample views, and all overlap counts."""
    for level in ['coarse','fine']:
        labels=pilot.obs['hicat_'+level+'_baseline'].astype(str)
        groups=sorted(labels.unique())
        fig,axes=plt.subplots(3,1,figsize=(18,14))
        _stacked(axes[0],review['phase_fractions'][level],level.capitalize()+' phase composition',PHASE_COLORS)
        for ax,score in zip(axes[1:],['S_phase','G2M']):
            values=[scores.loc[labels==g,score].dropna().to_numpy() for g in groups]
            ax.boxplot(values,showfliers=False,widths=.65,patch_artist=True,boxprops=dict(facecolor='#7ea9bc'))
            ax.set_xticks(np.arange(1,len(groups)+1));ax.set_xticklabels(groups,rotation=90,fontsize=7)
            ax.axhline(0,color='gray',ls='--',lw=.7);ax.set_ylabel('Score');ax.set_title(score+' score distributions (median and IQR)')
        fig.suptitle(CONTEXT+'\n'+level.capitalize()+' cell-cycle composition and score distributions',fontsize=15)
        percell=scores[['S_phase','G2M']].copy();percell.insert(0,'cluster',labels)
        cycle_source=book.table(level+'_cycle_per_cell',percell)
        book.save(fig,level+'_cycle_composition_distributions','Phase calls plus per-cell score distributions, no regression',[cycle_source,'cell_cycle/'+level+'_phase_fractions.tsv'])
        sample=metrics['sample_'+level]
        cluster_fraction=sample['cluster_fractions']
        if set(groups).issubset(cluster_fraction.columns):cluster_fraction=cluster_fraction.T
        reciprocal=sample['sample_fractions']
        if set(groups).issubset(reciprocal.index):reciprocal=reciprocal.T
        fig,axes=plt.subplots(2,1,figsize=(18,13))
        _stacked(axes[0],cluster_fraction,'Sample contribution within each '+level+' cluster')
        _stacked(axes[1],reciprocal,level.capitalize()+' cluster composition within each sample')
        fig.suptitle(CONTEXT+'\nAll 12 samples, two reciprocal views',fontsize=15)
        sources=[book.table(level+'_sample_within_cluster',cluster_fraction),book.table(level+'_cluster_within_sample',reciprocal)]
        book.save(fig,level+'_sample_reciprocal_composition','Both denominators are explicit; no genotype aggregation',sources)
    for key in ['seed_coarse','seed_fine','allen_fine']:
        fig,axes=plt.subplots(1,2,figsize=(21, max(8,len(metrics[key]['counts'])*.28+4)))
        sources=[]
        for ax,name in zip(axes,['counts','row_fractions']):
            frame=_overlap_order(metrics[key][name])
            _heatmap(ax,frame,key+' '+name.replace('_',' '),limits=(0,1) if name=='row_fractions' else None,fontsize=6)
            sources.append(book.table(key+'_'+name+'_matrix',frame))
        fig.suptitle(CONTEXT+'\n'+key.replace('_',' ').title(),fontsize=15)
        book.save(fig,key+'_count_fraction_overlap','Raw overlap counts and row-normalized fractions with parent blocks retained',sources)
    # Per-cell maturation versus phase shows continuous axes without fitting a trajectory.
    maturation=scores['later_neuronal_maturation']-scores['rg_stemness']
    labels=pilot.obs.hicat_fine_baseline.astype(str)
    fig,axes=plt.subplots(1,2,figsize=(20,8))
    for ax,program in zip(axes,['S_phase','G2M']):
        for i,cluster in enumerate(sorted(labels.unique())):
            select=labels==cluster
            ax.scatter(maturation[select],scores.loc[select,program],s=3,alpha=.35,color=plt.get_cmap('turbo')(i/38),rasterized=True)
            ax.text(maturation[select].median(),scores.loc[select,program].median(),cluster,fontsize=5)
        ax.set_xlabel('Later neuronal maturation − RG program score')
        ax.set_ylabel(program+' score');ax.set_title('Fine IDs labeled at cluster medians')
    fig.suptitle(CONTEXT+'\nDevelopmental maturation versus cell cycle (no pseudotime inference)',fontsize=15)
    frame=pd.DataFrame(dict(cluster=labels,maturation_minus_RG=maturation,S_score=scores.S_phase,G2M_score=scores.G2M))
    book.save(fig,'maturation_vs_cycle_scatter','Individual cells with fine memberships and cluster-median text',[book.table('maturation_cycle_scatter_values',frame)])
    pairs=metrics.get('phase_pairs',{}).get('pairs',pd.DataFrame())
    if not pairs.empty:
        category='identity_difference_persists_after_cell_cycle_stratification'
        summary=pd.crosstab(pairs.parent,pairs[category])
        fig,axes=plt.subplots(1,2,figsize=(18,7))
        summary.plot.bar(ax=axes[0],stacked=True);axes[0].set_ylabel('Number of sibling pairs');axes[0].set_title('Identity differences after phase stratification')
        tv=next((k for k in ['cycle_composition_tv','cell_cycle_total_variation'] if k in pairs),None)
        if tv and 'max_abs_unstratified_identity_effect' in pairs:
            for category_value in pairs[category].unique():
                selected=pairs[pairs[category]==category_value]
                axes[1].scatter(selected[tv],selected.max_abs_unstratified_identity_effect,label=category_value,s=15,alpha=.7)
            axes[1].set_xlabel('Cell-cycle composition total variation')
            axes[1].set_ylabel('Maximum absolute unstratified identity effect')
            axes[1].legend(fontsize=8);axes[1].set_title('Every within-parent pair; eligibility limits apply')
        fig.suptitle(CONTEXT+'\nWithin-phase canonical identity reassessment',fontsize=15)
        book.save(fig,'phase_stratified_boundary_review','Identity differences must remain within adequately populated matched phases',[book.table('phase_stratification_summary_counts',summary),book.table('phase_pair_plot_values',pairs,index=False)])
    ordered=review['fine'].copy()
    order_key={name:i for i,name in enumerate(DEVELOPMENT)}
    ordered['development_order']=ordered.dominant_developmental_program.map(order_key).fillna(99)
    ordered=ordered.sort_values(['development_order','dominant_regional_program','maturation_score'])
    programs=['rg_stemness','apical_RG_like','basal_RG_like','ipc_neurogenic','S_phase','G2M','neuroblast','immature_inhibitory','MGE','LGE','later_neuronal_maturation','endothelial','erythroid']
    values=review['program_means']['fine'].loc[ordered.index,programs]
    z=values.sub(values.mean()).div(values.std(ddof=0).replace(0,np.nan))
    fig,ax=plt.subplots(figsize=(15,15));_heatmap(ax,z,'Program-based review ordering, with regional and non-neural branches\nNot pseudotime; not anatomical localization',cmap='coolwarm',limits=(-2,2))
    fig.suptitle(CONTEXT+'\nDevelopmental coherence of the proposed hierarchy',fontsize=15)
    book.save(fig,'developmental_coherence_ordering','Order by independent leading developmental program, then region and maturation',[book.table('developmental_order_cluster_metadata',ordered[['dominant_developmental_program','dominant_regional_program','maturation_score']]),book.table('developmental_order_program_zscore',z)])


def _fine_review_pages(book,pilot,review,metrics,program_cfg,config):
    """Write every current fine cluster's seven evidence panels plus review text."""
    top_path=Path(config['baseline_run'])/'top20_genes_all_contexts.tsv'
    top=pd.read_csv(top_path,sep='\t')
    genes=list(dict.fromkeys(g for group in program_cfg['dotplot_groups'].values() for g in group))
    program_names=[p['name'] for p in program_cfg['programs'] if p['name'] not in ['ribosomal','mitochondrial','translation_proxy']]
    fine=review['fine']
    canonical_mean=review['gene_means']['fine']
    canonical_frac=review['gene_fractions']['fine']
    page_sources=[]
    for cluster,row in fine.iterrows():
        parent=row.coarse_parent
        # Marker selection is all predeclared canonical genes; never only DEGs.
        available=[g for g in genes if g in canonical_mean]
        fig=plt.figure(figsize=(20,15))
        gs=fig.add_gridspec(4,3,height_ratios=[1.1,1.2,1.35,1.8])
        ax=fig.add_subplot(gs[0,:])
        values=canonical_mean.loc[[cluster],available]
        fraction=canonical_frac.loc[[cluster],available]
        _dotplot(ax,values,fraction,available,'A. Independent canonical markers: color mean ln1pCPM; dot area detection')
        if ax.get_legend():ax.get_legend().remove()
        ax=fig.add_subplot(gs[1,0]);programs=review['program_means']['fine'].loc[cluster,program_names]
        ax.barh(np.arange(len(programs)),programs,color=['#549b81' if x>=0 else '#b4767d' for x in programs.fillna(0)])
        ax.set_yticks(np.arange(len(programs)));ax.set_yticklabels(programs.index,fontsize=6);ax.invert_yaxis()
        ax.axvline(0,color='gray',lw=.5);ax.set_title('B. Canonical program mean scores');ax.set_xlabel('Mean signal − controls',fontsize=8)
        ax=fig.add_subplot(gs[1,1]);phase=review['phase_fractions']['fine'].loc[cluster]
        ax.bar(phase.index,phase,color=PHASE_COLORS);ax.set_ylim(0,1);ax.set_title('C. Cell-cycle composition');ax.set_ylabel('Fraction of cluster')
        ax=fig.add_subplot(gs[1,2]);sample=metrics['sample_fine']['cluster_fractions']
        sample=sample.loc[cluster] if cluster in sample.index else sample[cluster]
        ax.bar(np.arange(len(sample)),sample,color='#668caa');ax.axhline(1/len(sample),color='black',ls='--',lw=.7)
        ax.set_xticks(np.arange(len(sample)));ax.set_xticklabels(sample.index,rotation=90,fontsize=6)
        ax.set_ylim(0,max(.3,float(sample.max())*1.15));ax.set_title('D. Sample composition; dashed = 1/12')
        for ax,key,title in [(fig.add_subplot(gs[2,0]),'seed_fine','E. Second-seed overlap'),(fig.add_subplot(gs[2,1]),'allen_fine','F. Allen-reference overlap')]:
            overlap=metrics[key]['row_fractions'].loc[cluster]
            overlap=overlap[overlap>0].sort_values(ascending=False)
            ax.bar(np.arange(len(overlap)),overlap,color='#648bb2' if key=='seed_fine' else '#ab7a5a')
            ax.set_xticks(np.arange(len(overlap)));ax.set_xticklabels(overlap.index,rotation=90,fontsize=6)
            ax.set_ylim(0,1);ax.set_ylabel('Fraction of baseline cluster')
            ax.set_title(title+'\nJaccard %.3f | %s' % ((row.seed_jaccard if key=='seed_fine' else row.Allen_jaccard),(row.seed_stability if key=='seed_fine' else row.Allen_boundary_survival)),fontsize=10)
        ax=fig.add_subplot(gs[2,2])
        context='fine_within_'+parent if (top.context=='fine_within_'+parent).any() else 'fine_global'
        markers=top[(top.cluster==cluster)&(top.context==context)].sort_values('rank_positive_marker').head(20)
        if markers.empty:
            _plain_text(ax,'G. Top positive inspection genes\n\nNo genes passed the existing positive-marker filters in '+context+'.\n\nAll-gene statistics remain available; this is not missing expression data.',size=10,width=40)
        else:
            ax.barh(np.arange(len(markers)),markers.delta_mean_log1p_cpm,color='#817497')
            ax.set_yticks(np.arange(len(markers)));ax.set_yticklabels(markers.gene_symbol,fontsize=6)
            ax.invert_yaxis();ax.set_xlabel('Mean ln1pCPM difference',fontsize=8)
            ax.set_title('G. Top %d positive inspection genes\n%s; exploratory Welch ranking'%(len(markers),context),fontsize=9)
        left='Hypothesis: '+str(row.provisional_identity)+'\n'
        left+='Independent support: '+str(row.hypothesis_support)+'; confidence: '+str(row.identity_confidence)+'\n'
        left+='Supporting canonical genes: '+(row.canonical_supporting_markers or 'None at the detection rule')+'\n'
        left+='Supporting programs: '+(row.supporting_programs or 'None above score threshold')+'\n'
        left+='Conflicting evidence: '+str(row.conflicting_evidence)+'\n'
        left+='Alternative interpretation: '+str(row.alternative_interpretation)
        right='Cell-cycle influence: '+str(row.phase_stratified_pair_summary)+'\n'
        right+='Sample influence: largest %.1f%% (%s); effective samples %.2f; concern=%s\n'%(100*row.top_sample_fraction,row.top_sample,row.effective_samples,row.sample_concern)
        right+='Seed: '+row.seed_stability+'; Allen boundary: '+row.Allen_boundary_survival+'\n'
        right+='Stress/QC: '+str(row.technical_flags)+'\n'
        right+='Additional markers: '+str(row.additional_markers)+'\n'
        right+='Review category: '+row.overall_validation_status+'\nRecommended action: '+row.recommended_action+' (review only; no merge or annotation applied)'
        _plain_text(fig.add_subplot(gs[3,:2]),left,size=9,width=140)
        _plain_text(fig.add_subplot(gs[3,2]),right,size=8,width=53)
        fig.suptitle('%s | parent %s | n=%s | %.2f%% of pilot\n%s'%(cluster,parent,format(row.n_cells,','),100*row.fraction_of_pilot,CONTEXT),fontsize=16)
        marker_source=book.table(cluster+'_positive_inspection_genes',markers,index=False)
        review_source=book.table(cluster+'_review_values',pd.DataFrame([row]),index=False)
        sources=['tables/figure_sources/fine_canonical_marker_means.tsv','tables/figure_sources/fine_canonical_marker_detection_fraction.tsv',
                 'marker_programs/fine_cluster_program_means.tsv','cell_cycle/fine_phase_fractions.tsv',
                 'tables/figure_sources/fine_sample_within_cluster.tsv','tables/figure_sources/seed_fine_row_fractions_matrix.tsv',
                 'tables/figure_sources/allen_fine_row_fractions_matrix.tsv',marker_source,review_source]
        book.save(fig,'fine_review_'+cluster,'Seven evidence panels and provisional-hypothesis review for '+cluster,sources)
        page_sources.append(dict(cluster=cluster,parent=parent,n_cells=row.n_cells,figure='figures/fine_review_'+cluster+'.pdf',marker_context=context,qualified_top_genes=len(markers)))
    pd.DataFrame(page_sources).to_csv(book.out/'annotation_review'/'fine_review_page_index.tsv',sep='\t',index=False)


def _scatter_full(ax,xy,values,title,limits=None,cmap='viridis'):
    """Rasterize all measured cells on inherited coordinates with robust limits."""
    values=np.asarray(values,dtype=float)
    finite=np.isfinite(values)
    if not finite.any():
        ax.scatter(xy[:,0],xy[:,1],c='#cccccc',s=.15,rasterized=True)
        ax.text(.5,.5,'Unavailable in assay',ha='center',transform=ax.transAxes)
        ax.set_title(title);ax.set_xticks([]);ax.set_yticks([])
        return (np.nan,np.nan)
    lo,hi=limits or _finite_limits(values,(2,98))
    # Draw low expression first so sparse high signal is not hidden by zeros.
    order=np.argsort(np.nan_to_num(values,nan=-np.inf),kind='stable')
    drawn=ax.scatter(xy[order,0],xy[order,1],c=values[order],s=.15,cmap=cmap,vmin=lo,vmax=hi,rasterized=True)
    ax.set_xticks([]);ax.set_yticks([]);ax.set_title(title,fontsize=10)
    ax.figure.colorbar(drawn,ax=ax,fraction=.032,pad=.02)
    return lo,hi


def _full_data_pages(book,pilot,full_context,program_cfg):
    """Display full-cell canonical expression/scores plus exact pilot overlays."""
    import anndata as ad
    path=full_context.get('object_path') if full_context else None
    if not path:
        raise ValueError('Full-data canonical object is required for Step07 review figures')
    full=ad.read_h5ad(path)
    xy=np.asarray(full.obsm['X_umap_step06_display'])
    program_names=list(full.uns.get('program_names',full.uns.get('validation_program_names')))
    scores=pd.DataFrame(full.obsm['X_validation_program_scores'],index=full.obs_names,columns=program_names)
    symbols=full.var['gene_symbol'].astype(str).to_numpy() if 'gene_symbol' in full.var else full.var_names.astype(str).to_numpy()
    genes=pd.DataFrame(full.X.toarray() if sparse.issparse(full.X) else np.asarray(full.X),index=full.obs_names,columns=symbols)
    # Save complete cell-level plot source; selected-gene sidecar already contains
    # every value. TSVs are source tables for export, not full-gene expression.
    coordinates=pd.DataFrame(xy,index=full.obs_names,columns=['UMAP1','UMAP2'])
    coordinates['is_pilot']=full.obs['is_pilot'].to_numpy()
    if 'technical_sample_id' in full.obs:coordinates['technical_sample_id']=full.obs.technical_sample_id.astype(str).to_numpy()
    coord_path=book.out/'full_data_projection'/'plot_coordinates.tsv.gz'
    coordinates.to_csv(coord_path,sep='\t',index_label='cell_id')
    gene_path=book.out/'full_data_projection'/'plot_gene_expression.tsv.gz'
    genes.to_csv(gene_path,sep='\t',index_label='cell_id',float_format='%.9g')
    score_path=book.out/'full_data_projection'/'plot_program_scores.tsv.gz'
    scores.to_csv(score_path,sep='\t',index_label='cell_id',float_format='%.9g')
    sources=['full_data_projection/plot_coordinates.tsv.gz','full_data_projection/plot_gene_expression.tsv.gz','full_data_projection/plot_program_scores.tsv.gz']
    sources.append(str(Path(path).relative_to(book.out)))
    pilot_indices=full.obs_names.get_indexer(pilot.obs_names)
    if np.any(pilot_indices<0):raise ValueError('Pilot cells missing from full-data display object')
    limits_rows=[]
    fig,axes=plt.subplots(2,2,figsize=(15,13))
    for ax,program in zip(axes.ravel()[:3],['rg_stemness','neuroblast','MGE']):
        lo,hi=_scatter_full(ax,xy,scores[program],program)
        limits_rows.append(dict(kind='program',name=program,vmin=lo,vmax=hi,percentiles='2,98'))
    ax=axes.ravel()[3];ax.scatter(xy[:,0],xy[:,1],s=.15,c='#d2d2d2',rasterized=True)
    ax.scatter(xy[pilot_indices,0],xy[pilot_indices,1],s=.65,c='#bf3f47',alpha=.6,rasterized=True)
    ax.set_title('Exact 12,000 pilot cells over all %s cells'%format(full.n_obs,','));ax.set_xticks([]);ax.set_yticks([])
    fig.suptitle(TITLE+'\nH. Fixed full-data Step 06 UMAP; canonical scores and pilot locations',fontsize=15)
    book.save(fig,'main_H_full_context','All full-data cells displayed; pilot locations highlighted without label transfer',sources,main=True)
    required_genes=['Sox2','Fabp7','Ascl1','Insm1','Mki67','Top2a','Dcx','Tubb3','Gad1','Gad2','Nkx2-1','Lhx6','Lhx8','Gsx2','Meis2','Sox11','Snap25']
    # Include every predeclared canonical gene, satisfying all minimum panels.
    gene_order=list(dict.fromkeys(required_genes+list(genes.columns)))
    for start in range(0,len(gene_order),9):
        selected=gene_order[start:start+9]
        fig,axes=plt.subplots(3,3,figsize=(17,15))
        for ax,gene in zip(axes.ravel(),selected):
            if gene in genes:
                lo,hi=_scatter_full(ax,xy,genes[gene],gene+' | ln1pCPM')
                limits_rows.append(dict(kind='gene',name=gene,vmin=lo,vmax=hi,percentiles='2,98'))
            else:_plain_text(ax,gene+' unavailable in selected assay',size=11)
        for ax in axes.ravel()[len(selected):]:ax.axis('off')
        fig.suptitle('Dissected E14.5 mouse MGE | all %s cells\nCanonical gene expression on fixed Step 06 coordinates'%format(full.n_obs,','),fontsize=15)
        book.save(fig,'full_data_genes_%02d'%(start//9+1),'Canonical gene expression; robust per-gene color limits, no new embedding',sources)
    for start in range(0,len(program_names),9):
        selected=program_names[start:start+9]
        fig,axes=plt.subplots(3,3,figsize=(17,15))
        for ax,program in zip(axes.ravel(),selected):
            lo,hi=_scatter_full(ax,xy,scores[program],program,cmap='coolwarm')
            limits_rows.append(dict(kind='program',name=program,vmin=lo,vmax=hi,percentiles='2,98'))
        for ax in axes.ravel()[len(selected):]:ax.axis('off')
        fig.suptitle('Dissected E14.5 mouse MGE | all %s cells\nFixed pilot-derived program controls; no full-data refitting'%format(full.n_obs,','),fontsize=15)
        book.save(fig,'full_data_programs_%02d'%(start//9+1),'Full-cell canonical scores using the fixed pilot-derived control model',sources)
    assignments=pilot.obs[['hicat_coarse_baseline','hicat_fine_baseline']].astype(str).copy()
    assignments.to_csv(book.out/'full_data_projection'/'pilot_overlay_memberships.tsv',sep='\t',index_label='cell_id')
    overlay_sources=['full_data_projection/plot_coordinates.tsv.gz','full_data_projection/pilot_overlay_memberships.tsv']
    for level,batch in [('coarse',4),('fine',9)]:
        labels=pilot.obs['hicat_'+level+'_baseline'].astype(str).to_numpy()
        clusters=sorted(set(labels))
        for start in range(0,len(clusters),batch):
            selected=clusters[start:start+batch]
            nrow=2 if level=='coarse' else 3;ncol=2 if level=='coarse' else 3
            fig,axes=plt.subplots(nrow,ncol,figsize=(16,14))
            for ax,cluster in zip(axes.ravel(),selected):
                ax.scatter(xy[:,0],xy[:,1],s=.12,c='#d2d2d2',rasterized=True)
                pick=pilot_indices[labels==cluster]
                ax.scatter(xy[pick,0],xy[pick,1],s=3,c='#ba394c',alpha=.7,rasterized=True)
                ax.set_title(cluster+' | n='+str(len(pick)),fontsize=11);ax.set_xticks([]);ax.set_yticks([])
            for ax in axes.ravel()[len(selected):]:ax.axis('off')
            fig.suptitle('Exact pilot '+level+' memberships on the full Step 06 manifold\nGray cells have no transferred HiCAT label; location is not identity proof',fontsize=15)
            book.save(fig,'full_data_'+level+'_overlays_%02d'%(start//batch+1),'Current pilot cluster locations only; background retains unlabeled full-data cells',overlay_sources)
    pd.DataFrame(limits_rows).drop_duplicates(['kind','name']).to_csv(book.out/'full_data_projection'/'plot_color_limits.tsv',sep='\t',index=False)
    return dict(full_cells=full.n_obs,full_gene_panels=len(gene_order),full_program_panels=len(program_names),pilot_cells_present=len(pilot_indices))


def _names(frame,column,value):
    """Compact human-readable cluster lists with explicit empty result."""
    selected=frame.loc[frame[column]==value,'cluster'].astype(str).tolist()
    return ', '.join(selected) if selected else 'None at the declared review thresholds'


def _write_readme(book,review,metrics,full_context,full_summary,relationship_flags,config):
    """Answer the user's thirteen review questions and end with all categories."""
    fine,coarse,parent=review['fine'],review['coarse'],review['parent']
    summary=review['summary']
    pairs=metrics.get('phase_pairs',{}).get('pairs',pd.DataFrame())
    issue_counts=relationship_flags.flag.value_counts().to_dict() if len(relationship_flags) else {}
    weak= fine[fine.Allen_boundary_survival=='no']
    cycle=fine[fine.cell_cycle_concern]
    sample=fine[fine.sample_concern]
    coverage=(full_context or {}).get('coverage',{})
    readiness=summary['readiness']
    text=['# Step 07: HiCAT hierarchy validation and parameter sensitivity','',
      '**'+readiness+' — current hierarchy remains IN_REVIEW.** '+('; '.join(summary['readiness_reasons']) or 'The declared minimum evidence rules pass, subject to the caveats below.')+'.','',
      'Biological context: **dissected E14.5 mouse medial ganglionic eminence (MGE)**. These samples are not organoids. Provisional names are user-supplied hypotheses and were never used to fit or tune clusters. No annotations are locked.','',
      '## Start with these assets','',
      '- [Main validation PDF](outputs/figures/hicat_hierarchy_validation_main.pdf): panels A–H across readable pages.',
      '- [Detailed validation PDF](outputs/figures/hicat_hierarchy_validation_detailed.pdf): distributions, exact overlap counts, full-data context, and all 38 seven-panel fine-cluster reviews.',
      '- [Fine decision table](outputs/annotation_review/cluster_validation_summary.tsv), [coarse decision table](outputs/annotation_review/coarse_cluster_validation_summary.tsv), [parent decision table](outputs/annotation_review/parent_validation_summary.tsv).',
      '- [Figure/source index](outputs/figures/figure_index.tsv), [hypothesis review](outputs/annotation_review/provisional_hypothesis_review.tsv), [hierarchy conflicts](outputs/annotation_review/hierarchy_relationship_flags.tsv).','',
      '## Inputs, processing and saved objects','',
      'The primary input is the frozen 12,000-cell expanded pilot: 1,000 cells from each of 12 samples, 19,071 genes, original raw-count X and saved ln(1+CPM) layer. Baseline and second-seed labels are reused without changing cell membership. The one new sensitivity fit uses the same cells and fixed baseline coarse parents, changing only fine q1=0.3→0.4 and qdiff=0.5→0.7; DE score remains 150 and the baseline seed is reused.','',
      'Canonical scores are means of predefined signature genes minus expression-bin-matched controls fitted without cluster labels. The same saved controls are applied to full-data raw counts after library-size normalization and natural log1p transformation. No integration, batch removal, regression, new full-data embedding, full-data clustering or nearest-neighbor label transfer occurs. Standard cycle phases are descriptive G1-like/S/G2M assignments; G1-like does not prove quiescence.','',
      'The saved [pilot validation AnnData](outputs/pcdh19_hicat_hierarchy_validation.h5ad) retains original all-gene raw counts and normalized expression, inherited coordinates and memberships, plus validation diagnostics and the separate Allen candidate. The [slot inventory](outputs/anndata_slot_inventory.json) lists X/layers/obs/obsm/obsp/uns; the [complete uns inventory](outputs/anndata_uns_inventory.json) exposes original provenance and the new hicat_validation metadata. The full-data [canonical sidecar](outputs/full_data_projection/canonical_full_data.h5ad) contains selected canonical genes, program scores, source metadata and inherited UMAP for all '+format(full_summary['full_cells'],',')+' cells. **It is not an all-gene normalized full-data HiCAT input.** Its selected-gene X is ln(1+CPM), not raw counts. Full source counts remain in the referenced Step 02/06 assets.','',
      '**Not saved or created:** a full-data HiCAT partition, inferred labels for nonpilot cells, locked biological annotations, batch-corrected/regressed/integrated expression, a refitted full-data UMAP, a pseudotime trajectory, or an adopted merge/removal of baseline clusters. Ribosomal Rpl/Rps features are unmeasured in this assay; ribosomal fractions/scores are unavailable, not zero. The separately named translation proxy is not a ribosomal fraction.','',
      '## Direct answers','',
      '1. **Are the coarse groups biologically interpretable?** '+ '; '.join('%s: %s; leading program %s'%(r.parent,r.recommended_status,r.dominant_developmental_state) for r in parent.itertuples())+'. Broad interpretability does not establish a stable taxonomy.','',
      '2. **Which coarse groups are primarily regional?** '+ '; '.join('%s: leading positive regional program %s'%(r.cluster,r.dominant_regional_program) for r in coarse.itertuples())+'. A leading score is not sufficient evidence that regional identity caused that coarse boundary; inspect the regional dotplot and mixed child programs.','',
      '3. **Which groups reflect progenitor-to-neuron maturation?** '+ '; '.join('%s: %s, later maturation score %.3f'%(r.cluster,r.dominant_developmental_program,r.maturation_score) for r in coarse.itertuples())+'. The developmental-order heatmap and maturation-versus-cycle scatter show these continuous relationships without inferring a trajectory. See the sibling maturation-boundary table where available.','',
      '4. **Which fine boundaries are strongly cycle-associated?** Clusters with at least one large phase-composition difference and no retained identity difference in eligible matched-phase sibling comparisons: '+(', '.join(cycle.cluster) if len(cycle) else 'none at the declared rules')+'. Pairwise results retain the exact partner IDs and minimum-cell eligibility.','',
      '5. **Do identity differences persist after stratification?** '+(str(pairs.identity_difference_persists_after_cell_cycle_stratification.value_counts().to_dict()) if len(pairs) else 'No eligible pair table was produced')+'. “Insufficient cells” is not equivalent to no identity difference. Scores exclude cycle genes for this reassessment; no expression regression was performed.','',
      '6. **Are any fine clusters dominated by a few samples?** '+(', '.join('%s (largest %.1f%%; effective samples %.2f)'%(r.cluster,100*r.top_sample_fraction,r.effective_samples) for r in sample.itertuples()) if len(sample) else 'No fine cluster exceeded the declared representation flags')+'. Sample association is descriptive and does not automatically imply a batch artifact.','',
      '7. **Which reproduce under the second seed?** HIGH_STABILITY: '+_names(fine,'seed_stability','HIGH_STABILITY')+'. MODERATE_STABILITY: '+_names(fine,'seed_stability','MODERATE_STABILITY')+'. The remaining clusters are LOW_STABILITY. These are two same-cell fits, not bootstrap or resampling stability.','',
      '8. **Which survive the Allen-reference condition?** Strong: '+_names(fine,'Allen_boundary_survival','strong')+'. Partial: '+_names(fine,'Allen_boundary_survival','partial')+'.','',
      '9. **Which disappear or merge?** No surviving boundary at the declared rule: '+(', '.join(weak.cluster) if len(weak) else 'none')+'. The [fine table](outputs/annotation_review/cluster_validation_summary.tsv) records actual best Allen matches, Jaccards, retained fractions and contributing merges/splits. “No” is an overlap-based boundary review; baseline IDs are preserved and no merge is adopted.','',
      '10. **Is the proposed developmental organization coherent?** Compare the independent program ordering with the fixed hierarchy. RG/apical/basal-like, IPC, cycling neuroblast and neuronal programs can overlap along development and across regions. Endothelial/erythroid candidates form non-neural alternatives. Anatomical vRG or oRG identity cannot be established from these RNA scores, particularly for an E14.5 mouse MGE dissection; the basal/oRG-like panel remains exploratory.','',
      '11. **Which provisional hypotheses are supported or contradicted?** '+str(summary['hypothesis_support_counts'])+'. Every coarse and fine row has independent supporting canonical genes/programs, conflicting evidence and alternatives in the [hypothesis table](outputs/annotation_review/provisional_hypothesis_review.tsv). These automated evidence rules are review aids, not final biological adjudication. C0001.F0011 has complete all-gene statistics but zero genes passing the prior positive top-marker filters; it is not a missing-data cluster.','',
      '12. **Does the pilot cover the full manifold?** All '+format(full_summary['full_cells'],',')+' cells are shown on inherited Step 06 coordinates; all '+format(full_summary['pilot_cells_present'],',')+' exact pilot cells were found. Full-data score/sample and UMAP-bin coverage diagnostics are saved in full_data_projection. Coverage summary: '+json.dumps(coverage,default=str)+'. Balanced sampling changes representation relative to the full dataset. Visual or binned coverage cannot certify that every biological population was captured.','',
      '13. **Reasons not to scale now?** '+('; '.join(summary['readiness_reasons']) if summary['readiness_reasons'] else 'No declared hard readiness criterion failed; remaining cautions still require review')+'. Hierarchy relationship flags: '+str(issue_counts)+'. Review unstable, cycle/maturation-associated, sample-skewed and non-neural boundaries before choosing a full-data strategy. This checkpoint does not start that run.','',
      '## Rules and practical limits','',
      'HIGH seed stability requires the configured high Jaccard plus precision/recall thresholds; MODERATE uses the declared intermediate Jaccard cutoff. Allen survival uses its separately declared strong/partial overlap thresholds. Exact thresholds are in the frozen config and [review rules](outputs/annotation_review/review_thresholds.json). Sample fractions/entropy, phase differences and QC contrasts are descriptive; none is interpreted from a p-value alone. Positive canonical scores and gene detection support a hypothesis but are not probabilities of identity. Top-20 Welch markers are exploratory tests on the cells used to discover clusters, not independent biological-replicate DE.','',
      'Canonical heatmap z-scores are for display only; absolute mean, median, positive fraction, cell count and per-cell source tables are saved. Dot color is mean ln1pCPM and area is detection fraction. Expression dendrograms use canonical gene means and describe similarity; the actual HiCAT nesting is shown separately. Robust full-data color limits are saved for every panel.','',
      '## Current cluster review categories','',
      '| Level | Cluster | Cells | Review category | Recommended action |','|---|---|---:|---|---|']
    for frame in [coarse,fine]:
        for row in frame.itertuples():text.append('| %s | %s | %d | %s | %s |'%(row.level,row.cluster,row.n_cells,row.overall_validation_status,row.recommended_action))
    # Package-level README intentionally ends with all 42 category rows.
    (book.root/'REVIEW_README.md').write_text('\n'.join(text)+'\n')
    (book.out/'REVIEW_README.md').write_text('\n'.join(text).replace('(outputs/','(')+'\n')


def create_report(run_root,pilot,score_frame,program_metadata,metrics,config,full_context=None):
    """Generate the complete Step07 review package from independently fit results.

    ``program_metadata`` is the complete scorer result mapping. ``config`` must
    contain ``programs_path``, ``hypotheses_path`` and ``baseline_run``. The input
    pilot must contain validation phases and the separate Allen fine candidate.
    Reports are computed after all clustering so user hypotheses cannot feed
    fitting. Returned summary is machine-readable; no approval is implied.
    """
    root=Path(run_root)
    for sub in ['tables','marker_programs','stability','sample_composition','cell_cycle','full_data_projection','parameter_sensitivity','annotation_review','validation']:
        (root/'outputs'/sub).mkdir(parents=True,exist_ok=True)
    review=review_clusters(pilot,score_frame,metrics,config['hypotheses_path'],root,config,program_metadata)
    program_cfg=json.loads(Path(config['programs_path']).read_text())
    book=ValidationReport(root,dpi=config.get('report_dpi',130));book.pilot=pilot
    try:
        genes=_summary_pages(book,review,program_cfg,program_metadata,metrics)
        _extended_diagnostics(book,pilot,score_frame,review,metrics)
        relationship_flags=_relationships(book,review,config,genes)
        # Incorporate actual unresolved cross-parent evidence in parent records.
        for index,row in review['parent'].iterrows():
            involved=relationship_flags[(relationship_flags.candidate=='baseline') & (relationship_flags.flag=='cross_parent_DE_failure') &
                ((relationship_flags.cluster_a.str.split('.').str[0]==row.parent)|(relationship_flags.cluster_b.str.split('.').str[0]==row.parent))]
            review['parent'].loc[index,'cross_parent_DE_issue']=bool(len(involved))
            review['parent'].loc[index,'cross_parent_DE_failure_count']=len(involved)
        review['parent'].to_csv(root/'outputs'/'annotation_review'/'parent_validation_summary.tsv',sep='\t',index=False)
        _fine_review_pages(book,pilot,review,metrics,program_cfg,config)
        full_summary=_full_data_pages(book,pilot,full_context,program_cfg)
        review['summary']['readiness_assessment_axes'] = dict(seed_stability=metrics['seed_fine']['summary'], Allen_boundary_survival=metrics['allen_fine']['summary'], sample_flagged_clusters=int(review['fine'].sample_concern.sum()), cycle_associated_clusters=int(review['fine'].cell_cycle_concern.sum()), maturation_associated_clusters=int(review['fine'].maturation_concern.sum()), technical_flagged_clusters=int(review['fine'].technical_concern.sum()), full_context_coverage=(full_context or {}).get('coverage', {}), hierarchy_flags=relationship_flags.flag.value_counts().to_dict(), interpretation='Conservative exploratory screen; passing overlap rules alone never produces automatic approval or a full-data run.')
        _write_readme(book,review,metrics,full_context,full_summary,relationship_flags,config)
    finally:
        book.close()
    summary=dict(review['summary'])
    summary.update(full_summary)
    summary.update(main_pdf='figures/hicat_hierarchy_validation_main.pdf',detailed_pdf='figures/hicat_hierarchy_validation_detailed.pdf',
        main_pages=book.main_pages,detailed_pages=book.detail_pages,individual_figures=len(book.index),fine_review_pages=len(review['fine']),
        hierarchy_relationship_flag_counts=relationship_flags.flag.value_counts().to_dict())
    (root/'outputs'/'annotation_review'/'review_summary.json').write_text(json.dumps(summary,indent=2,default=str)+'\n')
    return summary
