#!/usr/bin/env python3
"""One supplement: per-sample QC retention/distributions and cell-cycle UMAP comparison."""
from pathlib import Path
import argparse,hashlib,json,shutil,sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT=Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19')
QC=ROOT/'results/primary_processing/02_qc_filtering/02_qc_filtering_20260830_124611_97e1bb5/tables'
RUN=ROOT/'results/rebuild_20261103'


def render(output, min_genes=0):
    source=RUN/('cell_cycle_comparison_min%d'%min_genes if min_genes else 'cell_cycle_comparison')
    if not (source/'SUCCESS.json').exists():raise RuntimeError('Cell-cycle comparison has not completed')
    output.mkdir(parents=True,exist_ok=True)
    counts=pd.read_csv(QC/'step02_before_after_by_sample.tsv',sep='\t')
    qc=pd.read_csv(QC/'step02_per_cell_disposition.tsv.gz',sep='\t')
    counts['number']=counts['group'].str.rsplit('-',n=1).str[-1].astype(int)
    counts=counts.sort_values('number');samples=counts['group'].tolist();names=['JZ-'+str(i) for i in counts.number]
    removed=qc['step02_excluded'].astype(str).str.lower().eq('true').to_numpy()
    assert len(qc)==450788 and removed.sum()==4439
    original_removed=removed.copy()
    if min_genes:
        removed|=qc['n_genes_by_counts'].to_numpy()<min_genes
        counts['cells_after']=[int(((qc.technical_sample_id==sn).to_numpy()&~removed).sum()) for sn in samples]
        counts['cells_removed']=counts.cells_before-counts.cells_after
        counts['pct_removed']=100*counts.cells_removed/counts.cells_before
        counts['minimum_genes_additional_n']=[int(((qc.technical_sample_id==sn).to_numpy()&removed&~original_removed).sum()) for sn in samples]
        qc['below_minimum_detected_genes']=qc.n_genes_by_counts<min_genes
        qc['final_QC_excluded']=removed
        qc.to_csv(output/'per_cell_QC_disposition.tsv.gz',sep='\t',index=False,compression='gzip')
    nretained=int((~removed).sum());nremoved=int(removed.sum());additional=int((removed&~original_removed).sum())
    coords=np.load(source/'before_after_embeddings.npz',allow_pickle=False)
    obs=pd.read_csv(source/'cells_and_cycle_scores.tsv.gz',sep='\t',index_col=0)
    assert np.array_equal(obs.index.astype(str),coords['cell_id'])
    retained_ids=qc.loc[~removed,'cell_id'].to_numpy(dtype=str)
    retained_index=pd.Index(retained_ids); embedding_index=pd.Index(coords['cell_id'])
    assert retained_index.is_unique and embedding_index.is_unique
    assert len(retained_index)==len(embedding_index) and retained_index.difference(embedding_index).empty and embedding_index.difference(retained_index).empty
    assert int(counts.cells_before.sum())==450788 and int(counts.cells_after.sum())==nretained
    sample=obs.technical_sample_id.astype(str).to_numpy()
    palette=plt.get_cmap('tab20');colors={sn:palette(i) for i,sn in enumerate(samples)}
    rng=np.random.RandomState(20260903)
    display=np.concatenate([rng.choice(np.flatnonzero(sample==sn),min(6000,int((sample==sn).sum())),replace=False) for sn in samples]);rng.shuffle(display)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':9,'axes.labelsize':8,
                         'xtick.labelsize':7,'ytick.labelsize':7,'pdf.fonttype':42,'svg.fonttype':'none',
                         'axes.spines.top':False,'axes.spines.right':False})
    fig=plt.figure(figsize=(8.4,11))
    fig.text(.05,.978,'Supplement · sample QC and cell-cycle regression',fontsize=16,weight='bold')
    fig.text(.05,.953,'E14.5 mouse MGE · 12 samples',fontsize=10)
    fig.text(.05,.929,'450,788 input cells  →  {:,} excluded ({:.3f}%)  →  {:,} retained ({:.3f}%)'.format(nremoved,100*nremoved/len(qc),nretained,100*nretained/len(qc)),fontsize=11,weight='bold')
    fig.text(.05,.910,'19,071 genes retained · cell-cycle regression removes 0 additional cells',fontsize=9)
    fig.text(.065,.898,'JZ-1–3: WT Male · JZ-4–6: WT Female · JZ-7–9: HET Female · JZ-10–12: KO Male',fontsize=7)
    ax=fig.add_axes([.065,.754,.88,.118])
    x=np.arange(12)
    ax.bar(x,counts.cells_after/1000,color=[colors[s] for s in samples],width=.72,label='Retained')
    ax.bar(x,counts.cells_removed/1000,bottom=counts.cells_after/1000,color='#B74F52',width=.72,label='Excluded')
    for i,row in enumerate(counts.itertuples()):
        ax.text(i,row.cells_before/1000+1.1,'{:,} → {:,}'.format(row.cells_before,row.cells_after),rotation=35,ha='left',va='bottom',fontsize=6)
    ax.set_ylim(0,82);ax.set_xlim(-.6,11.85);ax.set_ylabel('Cells (thousands)');ax.set_xticks(x);ax.set_xticklabels(names)
    ax.set_title('A  Cells entering and retained after QC, by sample',loc='left',weight='bold',pad=7)
    fig.text(.63,.898,'Bars: retained · red caps: excluded',fontsize=6)
    for j,(metric,label,transform) in enumerate([
        ('total_counts','UMIs per cell (log₁₀)',np.log10),
        ('n_genes_by_counts','Detected genes (log₁₀)',np.log10),
        ('pct_counts_mt','Mitochondrial UMIs (%)',lambda v:v)]):
        ax=fig.add_axes([.065+j*.31,.574,.26,.125])
        for i,sn in enumerate(samples):
            mask=qc.technical_sample_id.eq(sn).to_numpy()
            for shift,selection,color in [(-.16,mask,'#CDD1D5'),(.16,mask&~removed,colors[sn])]:
                vals=transform(qc.loc[selection,metric].to_numpy(dtype=float))
                violin=ax.violinplot([vals],positions=[i+shift],widths=.28,showextrema=False,showmedians=True,points=60)
                for body in violin['bodies']:body.set_facecolor(color);body.set_edgecolor('none');body.set_alpha(.9)
                violin['cmedians'].set_color('#333333');violin['cmedians'].set_linewidth(.6)
        ax.set_xlim(-.6,11.6);ax.set_xticks(range(12));ax.set_xticklabels(names,rotation=65,ha='right',fontsize=5.8)
        if metric=='pct_counts_mt':
            ax.set_yscale('symlog',linthresh=1);ax.set_ylim(0,100)
            ax.set_yticks([0,1,10,100]);ax.set_yticklabels(['0','1','10','100'])
        ax.set_ylabel(label);ax.set_title(chr(ord('B')+j)+'  '+label.split(' (')[0],loc='left',weight='bold')
    fig.text(.065,.526,'QC distributions: left/gray = before; right/sample color = retained. All cells contribute to each density.',fontsize=7)
    qc_rule=('Per-sample 5-MAD QC removes 4,439 cells; ≥%d detected genes removes another %s cells.'%(min_genes,format(additional,',')) if min_genes else 'Per-sample 5-MAD QC: 0 low UMI; 67 low genes; 4,372 high mitochondrial exclusions.')
    fig.text(.065,.510,qc_rule,fontsize=7)
    for j,(key,title) in enumerate([('X_umap_before','E  Before cell-cycle regression · scaled'),('X_umap_after','F  After S/G2M regression · rescaled')]):
        ax=fig.add_axes([.065+j*.465,.285,.415,.18])
        xy=coords[key]
        ax.scatter(xy[display,0],xy[display,1],s=.28,c=[colors[s] for s in sample[display]],linewidths=0,rasterized=True)
        ax.set_title(title,loc='left',weight='bold');ax.set_xticks([]);ax.set_yticks([])
        ax.set_xlabel('UMAP 1');ax.set_ylabel('UMAP 2')
    handles=[Line2D([0],[0],marker='o',ls='',color=colors[sn],label=name,markersize=4) for sn,name in zip(samples,names)]
    fig.legend(handles=handles,loc='center',bbox_to_anchor=(.5,.254),ncol=12,frameon=False,fontsize=6,handletextpad=.2,columnspacing=.8)
    phase=obs.phase.astype(str).to_numpy();phase_colors={'G1':'#C4C8CD','S':'#DA9F4B','G2M':'#8755A3'}
    for j,(key,title) in enumerate([('X_umap_before','G  Before · original cell-cycle phase'),('X_umap_after','H  After · same cell-cycle phase labels')]):
        ax=fig.add_axes([.065+j*.465,.070,.415,.14])
        xy=coords[key]
        ax.scatter(xy[display,0],xy[display,1],s=.28,c=[phase_colors[s] for s in phase[display]],linewidths=0,rasterized=True)
        ax.set_title(title,loc='left',weight='bold');ax.set_xticks([]);ax.set_yticks([])
        ax.set_xlabel('UMAP 1');ax.set_ylabel('UMAP 2')
    fig.legend(handles=[Line2D([0],[0],marker='o',ls='',color=c,label=k,markersize=4) for k,c in phase_colors.items()],loc='center',bbox_to_anchor=(.5,.048),ncol=3,frameon=False,fontsize=7)
    correlations=pd.read_csv(source/'PC_cycle_correlations.tsv',sep='\t')
    correlations['absolute_r']=correlations.pearson_r.abs()
    maximum=correlations.groupby(['arm','score']).absolute_r.max()
    fig.text(.065,.030,'Max |PC–cycle score r|: S {:.3f} → {:.3f}; G2/M {:.3f} → {:.3f} (before → after)'.format(
        maximum['before','S_score'],maximum['after','S_score'],maximum['before','G2M_score'],maximum['after','G2M_score']),fontsize=7)
    fig.text(.065,.018,'Both UMAPs: same {:,} post-QC cells; display ≤6,000 cells/sample.'.format(nretained),fontsize=7)
    fig.text(.065,.006,'Original counts retained; only S/G2M scores regressed. No additional cell exclusions.',fontsize=7)
    for ext in ['png','pdf','svg']:fig.savefig(output/('PCDH19_QC_cell_cycle_supplement.'+ext),dpi=300,facecolor='white')
    counts.drop(columns='number').to_csv(output/'sample_cell_retention.tsv',sep='\t',index=False)
    methods='''Supplement: sample QC and cell-cycle correction in E14.5 dissected mouse MGE.
Input cells are the Cell Ranger filtered-barcode delivery, not all raw droplets. A: Cell totals before/after QC for each registered sample (JZ-1 through JZ-12); labels report input → retained; red segments show excluded cells.
B–D: Paired per-sample densities using all cells: gray/left before QC; sample color/right retained. UMI and gene counts displayed on log10 scales; mitochondrial percentages retain percent units on a symmetric-log axis (linear from 0–1%, logarithmic above 1%) so the retained distribution and excluded tail remain visible. Median bars are descriptive. Per-sample 5 scaled-MAD union rule excludes low counts OR low detected genes OR high mitochondrial percentage. Zero genes removed. Doublet filtering and ambient RNA correction were not applied; previously rejected doublet calls were not used.
E/F: Both UMAPs contain the identical 446349 approved post-QC cells. Before is the saved original Step06 UMAP after normalization/scaling; after refits PCA, graph and UMAP following S_score/G2M_score OLS regression on log-normalized HVGs and scaling. It is not an unscaled-versus-scaled comparison. Same raw-count-derived 3000 HVGs and baseline graph/layout parameters; seeded layouts can rotate independently. Colors identify technical samples; display is capped at 6000 cells per sample, while analysis fits all retained cells. No sample, genotype, sex, mitochondrial or library-size covariates are regressed.
The footer reports the maximum absolute Pearson correlation of either cycle score with any of the first 50 PCs, before and after, as a numerical diagnostic, not a guarantee that all biological cycle variation is removed. G/H: Original pre-regression cell-cycle phase labels shown on the two post-QC embeddings; phase labels are not recalculated to make correction look effective. Phase separation is descriptive; cell-cycle regression can also remove real developmental biology. Before-QC UMAP was not fitted or implied. Raw counts remain unchanged for later HiCAT/DE/Allen mapping.
'''
    if min_genes:
        methods=methods.replace('446349',str(nretained)).replace('Before is the saved original Step06 UMAP after normalization/scaling','Before is refitted after the additional minimum-gene filter, following normalization/scaling')
        methods+='\nAdditional QC: retain cells with at least %d detected genes, applied alongside the original per-sample MAD filters. Additional exclusions=%d; total retained=%d.\n'%(min_genes,additional,nretained)
    (output/'FIGURE_LEGEND.txt').write_text(methods)
    shutil.copy2(__file__,output/Path(__file__).name)
    files=[QC/'step02_per_cell_disposition.tsv.gz',QC/'step02_before_after_by_sample.tsv',source/'before_after_embeddings.npz',source/'cells_and_cycle_scores.tsv.gz',source/'parameters.json',source/'PC_cycle_correlations.tsv']
    def sha(p):
        h=hashlib.sha256()
        with p.open('rb') as f:
            for block in iter(lambda:f.read(1048576),b''):h.update(block)
        return h.hexdigest()
    (output/'PROVENANCE.json').write_text(json.dumps({'python':sys.version,'matplotlib':matplotlib.__version__,'numpy':np.__version__,'pandas':pd.__version__,'comparison_parameters':json.loads((source/'parameters.json').read_text()),'input_files':{str(p):sha(p) for p in files},'input_cells':450788,'excluded_cells':nremoved,'retained_cells':nretained,'minimum_detected_genes':min_genes,'additional_minimum_gene_exclusions':additional,'genes':19071,'UMAP_cell_identity_exact':True,'display_cells':len(display),'regression_covariates':['S_score','G2M_score']},indent=2))
    (output/'VALIDATION.json').write_text(json.dumps({'status':'complete','input_cells':len(qc),'excluded_cells':nremoved,'retained_cells':nretained,'minimum_detected_genes':min_genes,'regression_additional_exclusions':0,'same_post_QC_cells_in_both_embeddings':True,'exports':{ext:sha(output/('PCDH19_QC_cell_cycle_supplement.'+ext)) for ext in ['png','pdf','svg']}},indent=2))
    print('QC_CYCLE_SUPPLEMENT_COMPLETE',str(output),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=RUN/'supplementary_QC_cell_cycle')
    parser.add_argument('--min-genes',type=int,default=0)
    args=parser.parse_args();render(args.output,args.min_genes)
