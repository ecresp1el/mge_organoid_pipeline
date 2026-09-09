#!/usr/bin/env python3
"""Figures from frozen-target biological identity diagnostics; performs no search."""
from pathlib import Path
import os
import json
import shutil
import itertools
import gc
ROOT = Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results')
OUT = ROOT / 'div30_antecedent_identity_v1'
os.environ['MPLCONFIGDIR'] = str(OUT / 'cache/matplotlib')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import Normalize

GROUPS = {'A':'A · target retained', 'B':'B · target lost',
          'C':'C · non-target retained', 'D':'D · non-target removed'}
COLORS = dict(A='#087e8b', B='#2662a3', C='#df8f2d', D='#9a9fa8',
              retained='#087e8b', removed='#9a9fa8', all='#26364a')
GENES = ['NKX2-1','LHX6','SOX6','DLX1','DLX2','DLX5','DLX6','MAF','MAFB',
         'MEF2C','ERBB4','KCNC1','KCNC2','SLC6A1','SST','SATB1','TAC1','PVALB',
         'GRIA2','OPCML','NOTCH1']
MARKERS = ['GRIA2','OPCML','NOTCH1']
LABELS = dict(early_pv_score='Antecedent consensus',direct_pv_score='PV module',
              direct_sst_score='SST module',mge_score='MGE context',
              maturation_score='Neuronal maturation',progenitor_score='Progenitor',
              progression_score='Interneuron progression',cycle_score='Cell cycle',
              supplementary_gaba_score='GABA context',supplementary_stress_score='Stress context',
              n_genes='Detected genes',log_counts='log1p RNA counts',percent_mt='Mitochondrial %',
              SLC6A1='SLC6A1 expression')
UNITS = dict(early_pv_score='Mean percentile rank',direct_pv_score='DIV90 gene-z module',
             direct_sst_score='DIV90 gene-z module',SLC6A1='log1p(CP10k)',
             n_genes='Gene count',percent_mt='Percent mitochondrial RNA')

def read(name):
    return pd.read_csv(OUT/'within/tables'/name,sep='\t',float_precision='round_trip')

def ecdf(ax, x, label, color):
    x=np.sort(np.asarray(x,float));ax.plot(x,np.arange(1,len(x)+1)/len(x),label=label,color=color,lw=1.8)

def main():
    (OUT/'figures').mkdir(parents=True,exist_ok=True)
    needed=['retained','target','group','early_pv_score','direct_pv_score','direct_sst_score',
            'mge_score','maturation_score','progenitor_score','progression_score',
            'supplementary_gaba_score','cycle_score','supplementary_stress_score',
            'n_genes','percent_mt','GRIA2','OPCML','NOTCH1','SLC6A1']
    cells=pd.read_csv(OUT/'within/cells.tsv.gz',sep='\t',float_precision='round_trip',usecols=needed)
    assert len(cells)==90631 and cells.retained.sum()==15781 and (cells.retained&cells.target).sum()==9210
    config=json.loads((OUT/'within/provenance/frozen_gate.json').read_text())
    thresholds={r['gene']:r['threshold'] for r in config['rules']}
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':11,
        'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none',
        'savefig.facecolor':'white','figure.facecolor':'white'})
    inventory=[]
    with PdfPages(OUT/'DIV30_Biological_Identity_Figures.pdf') as book:
        def save(fig,name,caption):
            for ext in ['png','pdf','svg']:
                fig.savefig(OUT/'figures'/f'{name}.{ext}',dpi=600 if ext=='png' else 200,bbox_inches='tight')
            book.savefig(fig,bbox_inches='tight');plt.close(fig)
            inventory.append(dict(figure=name,caption=caption))
            gc.collect()
            print(name,flush=True)

        genes=read('group_gene_summary.tsv');genes=genes[genes.group.isin(GROUPS)&genes.gene.isin(GENES)]
        # Mean is compared only within each gene; dot size exposes RNA nondetection.
        mu=genes.pivot(index='gene',columns='group',values='mean').reindex(GENES)
        mx=mu.max(axis=1).replace(0,1)
        fig,ax=plt.subplots(figsize=(10,5))
        for j,g in enumerate(GENES):
            for i,group in enumerate(GROUPS):
                r=genes[(genes.gene==g)&(genes.group==group)].iloc[0]
                ax.scatter(j,i,s=220*r.positive_fraction+2,c=[r['mean']/mx[g]],cmap='viridis',vmin=0,vmax=1,edgecolors='none')
        ax.set(xticks=range(len(GENES)),xticklabels=GENES,yticks=range(4),yticklabels=list(GROUPS.values()),ylim=(-.6,3.6))
        ax.tick_params(axis='x',rotation=65);ax.invert_yaxis()
        ax.set_title('The four groups: individual RNA evidence',loc='left',fontweight='bold')
        fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,1),cmap='viridis'),ax=ax,label='Mean / maximum group mean for the same gene',pad=.02)
        for frac in [.1,.5,1]:ax.scatter([],[],s=220*frac+2,color='#60727c',label=f'{frac:.0%} detected')
        ax.legend(loc='upper left',bbox_to_anchor=(0,-.33),ncol=3,frameon=False)
        save(fig,'01_four_group_gene_evidence','All requested genes and three gate markers. Dot area is RNA detection fraction; color is within-gene relative mean, not cross-gene expression. A/B/C/D are frozen target and exact fixed-gate partitions.')

        axes_to_show=['early_pv_score','direct_pv_score','direct_sst_score','mge_score','maturation_score','progenitor_score','progression_score','supplementary_gaba_score','cycle_score','supplementary_stress_score','n_genes','percent_mt']
        fig,axs=plt.subplots(3,4,figsize=(14,10))
        for ax,col in zip(axs.flat,axes_to_show):
            for g,label in GROUPS.items():ecdf(ax,cells.loc[cells.group==g,col],label,COLORS[g])
            ax.set(title=LABELS.get(col,col),ylabel='Cumulative fraction',xlabel=UNITS.get(col,'DIV30 gene-z module'))
        axs[0,0].legend(fontsize=8,frameon=False);fig.suptitle('Biological programs and RNA complexity: full group distributions',fontweight='bold')
        fig.tight_layout();save(fig,'02_group_program_distributions','Empirical cumulative distributions use every cell. PV/SST are mean clipped DIV90-standardized gene z scores; MGE/maturation and other context panels use DIV30-standardized z scores. Antecedent is the mean of four percentile ranks. Supplementary GABA excludes SLC6A1. Stress is immediate-early/stress context, not a measured damage endpoint.')

        fig,axs=plt.subplots(1,3,figsize=(14,4.4))
        for g in ['all','retained','removed']:
            mask=np.ones(len(cells),bool) if g=='all' else cells.retained if g=='retained' else ~cells.retained
            ecdf(axs[0],cells.loc[mask,'early_pv_score'],g,COLORS[g])
        for g,label in GROUPS.items():
            axs[1].hist(cells.loc[cells.group==g,'early_pv_score'],bins=np.linspace(cells.early_pv_score.min(),cells.early_pv_score.max(),70),density=True,histtype='step',label=label,color=COLORS[g],lw=1.5)
        for ax in axs[:2]:ax.axvline(config['target_threshold'],color='#26364a',ls='--',lw=1);ax.set_xlabel('Frozen continuous antecedent score');ax.legend(fontsize=7,frameon=False)
        axs[0].set_ylabel('Cumulative fraction');axs[1].set_ylabel('Density')
        top=read('score_top_fraction_recovery.tsv');axs[2].plot(top.nominal_top_fraction*100,top.target_recovery*100,'o-',color=COLORS['A'])
        axs[2].axhline(cells.retained.mean()*100,color='#9a9fa8',ls=':',label='Random retention baseline')
        axs[2].set(xlabel='Top fraction of continuous score (%)',ylabel='Recovery of that fraction (%)',ylim=(0,100));axs[2].legend(fontsize=7,frameon=False)
        fig.suptitle('Enrichment across the score continuum',fontweight='bold');fig.tight_layout()
        save(fig,'03_continuous_antecedent','No target was redefined. Top 5/10/20/30/40% thresholds are descriptive sensitivity checks, with ties retained. Dashed vertical lines mark the original top20 threshold.')

        ab=read('gate_ablations.tsv').query('n_markers>0')
        fig,axs=plt.subplots(1,2,figsize=(12,5))
        y=np.arange(len(ab));axs[0].barh(y,100*ab.post_gate_target_fraction,color='#087e8b');axs[1].barh(y,100*ab.target_recovery,color='#2662a3')
        for ax in axs:ax.set_yticks(y,ab.markers);ax.invert_yaxis();ax.set_xlim(0,105)
        axs[0].set_xlabel('Target composition (%)');axs[1].set_xlabel('Target recovery (%)')
        for ax,col in zip(axs,['post_gate_target_fraction','target_recovery']):
            for j,v in enumerate(ab[col]):ax.text(v*100+1,j,f'{v:.1%}',va='center',fontsize=8)
        fig.suptitle('Which clauses contribute? Fixed-threshold ablations',fontweight='bold');fig.tight_layout()
        save(fig,'04_fixed_clause_contributions','Seven subsets of the original three exact clauses. No threshold optimization and no additional gate search. Fixed-pair results differ from the previously optimized pair.')

        branches=read('gate_decision_branches.tsv')
        fig,ax=plt.subplots(figsize=(15,8));ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
        def box(x,y,txt,fc='#e9eef3',fs=11):ax.text(x,y,txt,ha='center',va='center',fontsize=fs,bbox=dict(boxstyle='round,pad=.65',fc=fc,ec='#c5ccd4'))
        def arrow(x1,y1,x2,y2,label):
            ax.annotate('',xy=(x2,y2),xytext=(x1,y1),arrowprops=dict(arrowstyle='->',color='#586675'));ax.text((x1+x2)/2,(y1+y2)/2+.012,label,ha='center',fontsize=9)
        box(.5,.94,'90,631 DIV30 cells',fs=13)
        box(.5,.8,f'GRIA2 > {thresholds["GRIA2"]:.4f}')
        arrow(.5,.91,.5,.835,'')
        box(.64,.61,f'OPCML > {thresholds["OPCML"]:.4f}')
        arrow(.5,.765,.64,.645,'yes');arrow(.46,.765,.12,.33,'no')
        box(.78,.43,f'NOTCH1 ≤ {thresholds["NOTCH1"]:.4f}')
        arrow(.64,.575,.78,.465,'yes');arrow(.60,.575,.38,.33,'no')
        arrow(.745,.395,.64,.33,'no');arrow(.815,.395,.88,.33,'yes')
        for x,(_,r) in zip([.12,.38,.64,.88],branches.iterrows()):
            text=(f'{"RETAINED" if x==.88 else "Removed"}\nn = {int(r.n):,}\nAntecedent {r.antecedent_fraction:.1%}\nPV score {r.pv_score_mean:.2f}\nSST score {r.sst_score_mean:.2f}\nSLC6A1+ {r.SLC6A1_positive_fraction:.1%}\nMGE score {r.mge_score_mean:.2f}\nPostmitotic {r.postmitotic_fraction:.1%}')
            box(x,.17,text,'#d9efed' if x==.88 else '#f0f1f4',10)
        ax.set_title('The biology partitioned by the exact RNA rule',loc='left',fontweight='bold')
        save(fig,'05_fixed_gate_decision_tree','Sequential display order is GRIA2, OPCML, NOTCH1; the conjunction itself is order-independent. All four terminal branches sum to 90,631. Means use frozen program scores. Display thresholds rounded; full precision is in the report and frozen_gate.json.')

        color_axes=['early_pv_score','direct_pv_score','direct_sst_score','SLC6A1','maturation_score','mge_score']
        for a,b in itertools.combinations(MARKERS,2):
            fig,axs=plt.subplots(2,3,figsize=(14,9))
            for ax,col in zip(axs.flat,color_axes):
                h=ax.hexbin(cells[a],cells[b],C=cells[col],reduce_C_function=np.mean,gridsize=55,mincnt=1,cmap='viridis',rasterized=True)
                ax.axvline(thresholds[a],c='white',ls='--',lw=.8);ax.axhline(thresholds[b],c='white',ls='--',lw=.8)
                ax.set(xlabel=f'{a} · log1p CP10k',ylabel=f'{b} · log1p CP10k',title=LABELS[col])
                fig.colorbar(h,ax=ax,label='Bin mean · '+UNITS.get(col,'DIV30 gene-z module'),shrink=.8)
            fig.suptitle(f'{a} × {b}: every cell, full expression range',fontweight='bold');fig.tight_layout()
            save(fig,f'06_marker_plane_{a}_{b}','Hexagonal bins include all 90,631 cells, with no percentile clipping or subsampling; color is the mean of the labeled biological axis. White dashed lines show fixed RNA thresholds. Two-dimensional overlap cannot establish a discrete three-dimensional population.')
        fig,axs=plt.subplots(1,3,figsize=(14,4.6))
        for ax,(a,b) in zip(axs,itertools.combinations(MARKERS,2)):
            h=ax.hexbin(cells[a],cells[b],gridsize=60,bins='log',mincnt=1,cmap='magma',rasterized=True)
            ax.axvline(thresholds[a],c='white',ls='--',lw=.8);ax.axhline(thresholds[b],c='white',ls='--',lw=.8)
            ax.set(xlabel=a+' · log1p(CP10k)',ylabel=b+' · log1p(CP10k)');fig.colorbar(h,ax=ax,label='Cell count (log scale)')
        fig.suptitle('Marker-space density, including RNA-zero mass',fontweight='bold');fig.tight_layout()
        save(fig,'07_marker_plane_density','Companion density plots prevent sparsely populated colored regions being mistaken for most cells. Full ranges and every cell are included.')

        slc=read('slc6a1_group_summary.tsv').set_index('group');order=['all','target','non_target','retained','removed','A','B','C','D']
        fig,axs=plt.subplots(1,3,figsize=(15,4.5))
        axs[0].bar(range(len(order)),slc.loc[order,'positive_fraction']*100,color=[COLORS.get(g,'#526a85') for g in order])
        axs[0].set(xticks=range(len(order)),xticklabels=order,ylabel='SLC6A1 detected (%)',ylim=(0,100));axs[0].tick_params(axis='x',rotation=55)
        for g,label in GROUPS.items():ecdf(axs[1],cells.loc[cells.group==g,'SLC6A1'],label,COLORS[g])
        axs[1].set(xlabel='SLC6A1 · log1p CP10k',ylabel='Cumulative fraction');axs[1].legend(fontsize=7,frameon=False)
        h=axs[2].hexbin(cells.early_pv_score,cells.SLC6A1,gridsize=55,bins='log',mincnt=1,cmap='magma',rasterized=True)
        axs[2].set(xlabel='Frozen antecedent score',ylabel='SLC6A1 · log1p CP10k');fig.colorbar(h,ax=axs[2],label='Cell count')
        fig.suptitle('SLC6A1 is a separate readout, with partial dependence in the old consensus',fontweight='bold');fig.tight_layout()
        save(fig,'08_slc6a1_axis','SLC6A1 is excluded from the explicit four-gene PV module and from new context panels. It was already used in frozen mapping and sparse prediction, so consensus coupling is not wholly independent. No target was changed.')

        models=read('model_oof_metrics.tsv')
        wanted=['A_requested','B_requested_plus_PV','A_extended_generic','B_extended_plus_PV','A_flexible_generic','B_flexible_plus_PV']
        model_labels=['Requested maturation','Requested + PV module','Extended maturation','Extended + PV module','Flexible maturation','Flexible + PV module']
        fig,axs=plt.subplots(1,2,figsize=(13,5))
        for ax,outcome in zip(axs,['antecedent_target','fixed_gate_retained']):
            tab=models[models.outcome==outcome].set_index('model').reindex(wanted)
            ax.barh(range(len(wanted)),tab.roc_auc,color=['#9a9fa8','#087e8b']*3)
            ax.set(yticks=range(len(wanted)),yticklabels=model_labels,xlim=(.45,1),xlabel='Pooled held-sample ROC AUC',title='Antecedent target' if outcome=='antecedent_target' else 'Fixed gate membership');ax.invert_yaxis()
            for j,v in enumerate(tab.roc_auc):ax.text(v+.004,j,f'{v:.3f}',va='center',fontsize=8)
        fig.suptitle('What does the PV module add beyond measured maturation?',fontweight='bold');fig.tight_layout()
        save(fig,'09_generic_vs_pv_models','Prediction uses six held-sample folds. Requested A uses exactly progenitor, maturation and postmitotic status; B adds the direct PV module. Extended adds DCX, cycle, GABA and RNA-complexity terms; flexible uses additive splines. Grey models use generic developmental covariates; green adds PV. The PV module is part of the target definition, and the fixed gate was already selected on that target: incremental prediction is not independent fate validation. Full AP, log-loss, Brier and fold metrics are supplied.')

        perf=read('line_condition_performance.tsv').query("grouping=='cell_line+condition'")
        display_groups=[f'{r.cell_line} · {r.condition.split()[0]}' for r in perf.itertuples()]
        fig,axs=plt.subplots(1,3,figsize=(14,4.5));y=np.arange(len(perf))
        axs[0].barh(y-.17,100*perf.starting_target_fraction,height=.32,color='#9a9fa8',label='Starting');axs[0].barh(y+.17,100*perf.post_gate_target_fraction,height=.32,color='#087e8b',label='Retained')
        axs[1].barh(y,100*perf.target_recovery,color='#2662a3');axs[2].barh(y,perf.retained_n,color='#df8f2d')
        for ax in axs:ax.set_yticks(y,display_groups);ax.invert_yaxis()
        axs[0].set_xlabel('Target composition (%)');axs[0].legend(fontsize=8,frameon=False);axs[1].set_xlabel('Target recovery (%)');axs[2].set_xlabel('Retained cells')
        fig.suptitle('Cell-line and CV/MW frequency structure',fontweight='bold');fig.tight_layout()
        save(fig,'10_line_condition_structure','Each line-condition combination is represented by one cached source sample. Differences are frequency structure, not a consistency pass/fail test or independent condition effects. Existing gate-reselection held-sample results remain unchanged.')
    (OUT/'figures/figure_captions.json').write_text(json.dumps(inventory,indent=2)+'\n')
    shutil.copy2(__file__,OUT/'provenance'/Path(__file__).name)
    print('Rendered',len(inventory),'figure plates',flush=True)

if __name__=='__main__':main()
