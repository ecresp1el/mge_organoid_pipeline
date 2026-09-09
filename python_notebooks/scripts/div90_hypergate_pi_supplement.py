#!/usr/bin/env python3
"""PI-package technical supplement, rendered only from completed phase-two outputs."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import textwrap

ROOT = Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_hypergate_sst_pv_phase2')
OUT = ROOT / 'pi_figure_package_v1'
os.environ['MPLCONFIGDIR'] = str(OUT / 'cache/matplotlib_supplement')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import FancyBboxPatch, Patch
import numpy as np
import pandas as pd

STATES = ['PV-biased', 'PV/SST hybrid', 'SST-biased', 'unresolved/immature']
LABELS = ['PV-biased', 'Dual-high', 'SST-biased', 'Unresolved']
SHORT = ['pv', 'hybrid', 'sst', 'unresolved']
COLORS = ['#326DAB', '#39876F', '#CB733F', '#A6ACB6']
LINES = ['H9', '79B', '2E']
INK = '#243346'
GRAY = '#566473'
CONDITION_NOTE = 'CV = higher glucose (~2x); MW = lower glucose (~1x). Condition is confounded with operator/culture differences; no glucose causality is claimed.'


def read(name):
    return pd.read_csv(ROOT / 'tables' / name, sep='\t')


def configure():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 12, 'axes.titlesize': 15,
        'axes.labelsize': 12, 'axes.spines.top': False, 'axes.spines.right': False,
        'axes.titlepad': 12, 'axes.labelcolor': INK, 'text.color': INK,
        'xtick.labelsize': 11, 'ytick.labelsize': 11, 'pdf.fonttype': 42,
        'ps.fonttype': 42, 'svg.fonttype': 'none', 'savefig.facecolor': 'white'})


def canvas(number, title, subtitle, footer='', shape=(1, 1)):
    fig = plt.figure(figsize=(15, 9), facecolor='white')
    fig.text(.06, .956, f'S{number}  {title}', fontsize=23, weight='bold', va='top')
    fig.text(.06, .895, subtitle, fontsize=13, color=GRAY, va='top')
    fig.text(.06, .023, '\n'.join(textwrap.wrap(footer, 165)), fontsize=10.5, color=GRAY, va='bottom')
    grid = fig.add_gridspec(*shape, left=.085, right=.95, bottom=.15, top=.82, wspace=.42, hspace=.57)
    return fig, grid


def state_legend(fig, y=.075):
    fig.legend([Patch(facecolor=c) for c in COLORS], LABELS, ncol=4, loc='lower center',
               bbox_to_anchor=(.5, y), frameon=False, fontsize=12, columnspacing=2)


def condition_legend(ax):
    ax.legend([Patch(facecolor='white', edgecolor=GRAY), Patch(facecolor='white', edgecolor=GRAY, hatch='///')],
              ['CV', 'MW'], frameon=False, loc='upper left', fontsize=11, ncol=2)


def box(ax, x, y, width, height, title, body, color='#F1F4F7', fontsize=13):
    p = FancyBboxPatch((x, y), width, height, boxstyle='round,pad=0.014,rounding_size=.02',
                      transform=ax.transAxes, facecolor=color, edgecolor='#D4DAE0', linewidth=1)
    ax.add_patch(p)
    ax.text(x+.025, y+height-.035, title, transform=ax.transAxes, fontsize=15, weight='bold', va='top')
    ax.annotate(body, (x+.025, y+height-.035), xycoords='axes fraction', xytext=(0,-30),
                textcoords='offset points', fontsize=fontsize, va='top', linespacing=1.4)


def stack(ax, rows, labels, n=None, title=None):
    vals = np.array(rows)*100
    base = np.zeros(len(vals))
    for j, (color, name) in enumerate(zip(COLORS, LABELS)):
        bars = ax.bar(np.arange(len(vals)), vals[:, j], bottom=base, color=color, width=.58,
                      edgecolor='white', linewidth=1)
        for i, b in enumerate(bars):
            if vals[i, j] >= 5:
                ax.text(b.get_x()+b.get_width()/2, base[i]+vals[i,j]/2, f'{vals[i,j]:.1f}%',
                        ha='center', va='center', fontsize=12, color='white' if j != 3 else INK)
        base += vals[:, j]
    ax.set_xticks(np.arange(len(vals)), labels)
    ax.set_ylim(0, 112 if n is not None else 100)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_ylabel('Cells in each operational state (%)')
    if n is not None:
        for i, v in enumerate(n):
            ax.text(i, 103, f'n = {int(v):,}', ha='center', fontsize=12,
                    weight='bold' if int(v) == 37 else 'normal')
    if title:
        ax.set_title(title, loc='left')


def s1(cells, params):
    fig, grid = canvas(1, 'Separate scores capture overlapping developmental signals',
        'The score definitions are explicit; quadrant membership does not establish a distinct cell fate.',
        'Histograms and correlations display the completed analysis. No score fitting or biological re-analysis is performed for this package.', (2, 2))
    for j, (program, color) in enumerate(zip(['pv','sst'], [COLORS[0], COLORS[2]])):
        ax = fig.add_subplot(grid[0,j])
        ax.hist(cells[program+'_score'], bins=42, color=color, alpha=.85)
        ax.axvline(params[program+'_threshold'], color=INK, linestyle='--', linewidth=1.6)
        ax.set(xlabel=f'{program.upper()}-associated score', ylabel='Cells', title=f'{"AB"[j]}  {program.upper()} score distribution')
        ax.text(.98, .95, 'Dashed line: primary cutoff', transform=ax.transAxes, ha='right', va='top', fontsize=11)
    variance = read('module_variance_contributions.tsv')
    ax = fig.add_subplot(grid[1,0])
    ax.bar(variance.gene, variance.score_variance_share*100,
           color=[COLORS[0] if p == 'pv' else COLORS[2] for p in variance.program])
    ax.tick_params(axis='x', rotation=45)
    ax.set(ylabel='Share of score variance (%)', ylim=(0, 36), title='C  SST alone contributes 25.1% of SST-score variance')
    corr = read('module_gene_spearman.tsv').set_index('Unnamed: 0')
    ax = fig.add_subplot(grid[1,1])
    im = ax.imshow(corr, vmin=-.25, vmax=1, cmap='RdBu_r', aspect='auto')
    ax.set_xticks(range(len(corr)), corr.columns, rotation=45, ha='right', fontsize=9)
    ax.set_yticks(range(len(corr)), corr.index, fontsize=9)
    ax.set_title('D  Gene coordination is weak', loc='left')
    cb = fig.colorbar(im, ax=ax, fraction=.045, pad=.03)
    cb.set_label('Spearman correlation')
    return fig, 'Equal-weight averages of gene-wise population z-scores, clipped at −3 to +3, of frozen log1p(CP10K) RNA. PV genes: MEF2C, KCNC1, KCNC2, TAC1. SST genes: SST, SATB1, NR2F2, CDK14, CACNG3. Independent pooled medians: PV ≥ '+str(params['pv_threshold'])+'; SST ≥ '+str(params['sst_threshold'])+'. Variance share = Cov(gene contribution, score)/Var(score). Sources: state_parameters.json; module_variance_contributions.tsv; module_gene_spearman.tsv. Modules reflect developmental associations, not lineage-specific classifiers. See ../evidence/DIV90_HYPERGATE_PHASE2_EVIDENCE.md.'


def s2(cells, params, hybrid):
    fig, grid = canvas(2, 'Only a subset of dual-high cells has multiple detected genes per program',
        '1,309 operational dual-high cells; 356 detect ≥2 genes per program; 157 meet the same support without SST RNA.',
        'Co-detection supports a subset but does not establish a coordinated hybrid fate. QC exclusions below are sensitivity analyses; all 4,768 cells remain in the primary analysis.', (2,2))
    genes = params['modules']['pv']['genes']+params['modules']['sst']['genes']
    data = read('module_gene_by_state.tsv').pivot(index='state', columns='gene', values='fraction_expressing').loc[STATES, genes]*100
    ax = fig.add_subplot(grid[0,0]); im = ax.imshow(data, vmin=0, vmax=100, cmap='Blues', aspect='auto')
    ax.set_xticks(range(len(genes)), genes, rotation=45, ha='right', fontsize=10)
    ax.set_yticks(range(4), LABELS, fontsize=11); ax.set_title('A  RNA detection within each state', loc='left')
    fig.colorbar(im, ax=ax, fraction=.035, pad=.025, label='Detected (%)')
    ax = fig.add_subplot(grid[0,1])
    vals = [hybrid['hybrid_n'], hybrid['hybrid_multiple_genes_both_n'], hybrid['hybrid_multiple_genes_both_without_SST_n']]
    ax.barh([2,1,0], vals, color=[COLORS[1], '#71AB8E', '#ADD0BB'])
    ax.set_yticks([2,1,0], ['All dual-high','≥2 genes / program','≥2 SST-support genes\nwithout SST + ≥2 PV'])
    ax.set_xlim(0, 1540); ax.set_xlabel('Dual-high cells'); ax.set_title('B  Nested co-detection support', loc='left')
    for y, v in zip([2,1,0], vals):ax.text(v+25, y, f'{v:,}', va='center', weight='bold')
    ax = fig.add_subplot(grid[1,0])
    ax.boxplot([cells.loc[cells.state.eq(s), 'total_counts'] for s in STATES], tick_labels=LABELS,
               showfliers=False, patch_artist=True, boxprops={'facecolor':'#E4E9EE'}, medianprops={'color':INK})
    ax.set_yscale('log'); ax.set_ylabel('RNA counts per cell'); ax.tick_params(axis='x', rotation=20)
    ax.set_title('C  Complexity differs between operational groups', loc='left')
    qc = read('qc_exclusion_sensitivity.tsv')
    ax = fig.add_subplot(grid[1,1]); yy=np.arange(len(qc))
    names=['All cells','Exclude high-count tail','Exclude high-mito tail','Exclude high-stress tail','Exclude predicted doublets','Exclude top 10% doublet scores']
    ax.barh(yy, qc[STATES[1]+'_fraction']*100, color=COLORS[1], alpha=.8)
    ax.set_yticks(yy,names,fontsize=10); ax.invert_yaxis();ax.set_xlim(0,34)
    ax.set_xlabel('Dual-high fraction (%)');ax.set_title('D  Recorded QC sensitivities',loc='left')
    return fig, 'The support groups are nested subsets of the 1,309 primary dual-high cells. Without-SST support requires at least two other SST-module genes and at least two PV-module genes. RNA box plots show medians/interquartile ranges and whiskers, with outliers omitted only from display. Source QC sensitivities include per-sample 95th-percentile count/mitochondrial/stress tails and exploratory Scrublet estimates; primary cells are not excluded. Per-sample Scrublet used all 22,338 available DIV90 cells, assumed 5% doublets, and called five entry cells, three dual-high. Homotypic doublets are not ruled out. Sources: hybrid_validation.json; module_gene_by_state.tsv; qc_exclusion_sensitivity.tsv; cells.tsv.gz.'


def s3(hybrid):
    null=read('hybrid_depth_stratified_null.tsv')
    fig, grid=canvas(3,'Overall dual-high abundance is not enriched over the conditional null',
        'Independent gene permutations preserve sample and approximate RNA complexity; calibration and score cutoffs stay fixed.',
        '200 saved permutations; this null is conditional on sample and RNA-count quintile and does not remove every technical or developmental confound.',(2,2))
    for j,(column,observed,title) in enumerate([
        ('dual_high_n',hybrid['hybrid_n'],'A  All dual-high cells'),
        ('multi_gene_dual_high_n',hybrid['hybrid_multiple_genes_both_n'],'B  Dual-high cells with ≥2 genes per program')]):
        ax=fig.add_subplot(grid[0,j]);ax.hist(null[column],bins=20,color='#BFCBD5',edgecolor='white')
        ax.axvline(observed,color=COLORS[1],lw=2.5)
        ax.text(.03,.93,f'Observed: {observed:,}\nNull mean: {null[column].mean():,.1f}',transform=ax.transAxes,va='top',fontsize=13,
                bbox={'facecolor':'white','edgecolor':'none','alpha':.9,'pad':3})
        ax.set(xlabel='Cells',ylabel='Permutations',title=title)
    ax=fig.add_subplot(grid[1,0]);ax.axis('off')
    box(ax,0,.03,.97,.88,'What the null supports',
        'Overall abundance: no excess over the null.\n\nMultiple-gene subset: a modest excess.\n\nNeither result establishes a new cell type.',fontsize=14)
    ax=fig.add_subplot(grid[1,1]);ax.axis('off')
    box(ax,0,.03,.97,.88,'Depth adjustment reduces score correlation',
        f'PV/SST score correlation: {hybrid["rho_pv_sst"]:.3f}\nAfter within-sample depth diagnostic: {hybrid["rho_depth_diagnostic"]:.3f}\n\nMedian within-module gene correlation:\nPV {hybrid["median_within_pv_gene_rho"]:.3f}; SST {hybrid["median_within_sst_gene_rho"]:.3f}',fontsize=13)
    return fig, 'Saved null results: all dual-high observed 1,309 versus mean 1,347.97, empirical 95% interval 1,311.975–1,391.025; multi-gene dual-high observed 356 versus mean 322.57, interval 301.975–351.2. Genes were independently permuted within sample × RNA-count quintile. The depth diagnostic centers RNA-count and detected-gene covariates within sample, preserving sample score means; it is not the primary scoring model and does not regress out culture condition. No new permutations or models were run for this figure. Sources: hybrid_depth_stratified_null.tsv; hybrid_validation.json.'


def s4():
    baseline=read('condition_baseline.tsv')
    d=baseline[baseline.grouping.eq('culture_operator|cell_line')]
    fig,grid=canvas(4,'Culture associations differ across the three cell lines',
        'Each cell line is shown separately; the primary display does not pool CV and MW across lines.',CONDITION_NOTE,(1,3))
    for j,line in enumerate(LINES):
        p=d[d.cell_line.eq(line)].set_index('culture_operator').loc[['CV','MW']]
        ax=fig.add_subplot(grid[0,j]); stack(ax,p[[x+'_fraction' for x in SHORT]],['CV','MW'],p.n,title=line)
        if j:ax.set_ylabel('')
        if line=='2E':ax.text(.5,.95,'Only 37 CV entry cells',transform=ax.transAxes,ha='center',weight='bold',color='#934A2C',fontsize=12)
    state_legend(fig)
    return fig, 'User-supplied condition metadata: CV corresponds to odd DIV90 sample IDs 1/3/5 and higher glucose (~2x); MW corresponds to even IDs 2/4/6 and lower glucose (~1x). All technical IDs contain MW, which is not the operator field. H9: CV n=505, MW n=224; 79B: CV n=1,039, MW n=1,755; 2E: CV n=37, MW n=1,208. One sample per condition per line prevents a replicated within-line glucose effect estimate. Operational state thresholds are fixed globally. Source: condition_baseline.tsv; original sample and line identities were verified in condition_metadata_verified.tsv.'


def gate_condition(number, gate_name, title, subtitle):
    allgate=read('condition_gate_metrics.tsv')
    d=allgate[allgate.grouping.eq('culture_operator|cell_line') & allgate.gate_name.eq(gate_name)]
    fig,grid=canvas(number,title,subtitle,CONDITION_NOTE+'  The 2E/CV entry group has 37 cells, including only two PV-biased cells.',(2,3))
    for j,line in enumerate(LINES):
        p=d[d.cell_line.eq(line)].set_index('culture_operator').loc[['CV','MW']]
        ax=fig.add_subplot(grid[0,j]);x=np.arange(2)
        for k,condition in enumerate(['CV','MW']):
            y=np.array([p.loc[condition,'pv_recovery'],p.loc[condition,'hybrid_recovery']])*100
            ax.bar(x+(k-.5)*.32,y,width=.30,color=COLORS[:2],hatch='///' if k else '',edgecolor='white',linewidth=.8)
            for a,b in zip(x+(k-.5)*.32,y):ax.text(a,b+2,f'{b:.1f}',ha='center',fontsize=10)
        ax.set_xticks(x,['PV-biased','Dual-high']);ax.set_ylim(0,123);ax.set_yticks([0,25,50,75,100])
        ax.set_title(f'{line}: state recovery',loc='left');condition_legend(ax)
        if j==0:ax.set_ylabel('Starting state retained (%)')
        ax=fig.add_subplot(grid[1,j]);stack(ax,p[['retained_'+x+'_fraction' for x in SHORT]],['CV','MW'],p.retained_n,title='Retained composition')
        if j:ax.set_ylabel('')
    state_legend(fig)
    return fig


def s5(summary):
    fig=gate_condition(5,'phase1_reference','The old gate loses dual-high cells in every sample',
        'FAT3-low / PTPRM-undetected is too restrictive for a strategy that preserves overlapping developmental programs.')
    return fig, 'Frozen phase-one retention rule: FAT3 ≤ 0.764341 AND PTPRM ≤ 0, in log1p(CP10K) RNA. Overall PV-biased recovery is 79.1%; dual-high recovery is 47.9%, losing 682 cells. Each panel applies the same global rule; there is no per-condition reoptimization. Pooled dual-high recovery is 69.4% CV and 36.6% MW; individual line panels preserve heterogeneity. FAT3/PTPRM RNA mean shifts persist within operational states: pooled CV−MW FAT3 raw −0.398, within-state component −0.396; PTPRM raw −0.242, within-state −0.242. Associations with the dual-high state remain positive in both conditions and stronger in MW. These are descriptive, not causal glucose effects. Sources: condition_gate_metrics.tsv; condition_marker_shift_decomposition.tsv; condition_marker_state_effects.tsv.gz.'


def s6(summary):
    fig=gate_condition(6,'experimental_depletion','The current candidate preserves most dual-high cells in both conditions',
        'Remove FGFR2-high OR PTPRS-high; retain the complement. Performance varies by cell line and culture.')
    return fig, 'Experimental candidate RNA removal rule: FGFR2 > 1.2491154 OR PTPRS > 3.4076235 in log1p(CP10K). The analyst-defined comparison requires at least 20% SST-biased removal; this is a decision constraint, not a biological threshold or user-mandated purity target. Overall retains 4,024 cells: PV recovery 88.7%, dual-high 89.3%, total 84.4%. SST-biased contamination changes from 22.55% to 21.27% (exact difference 1.27 percentage points before display rounding), a modest benefit. Pooled dual-high recovery CV 93.3%, MW 87.2%. 2E/CV contains only two baseline PV-biased cells; its 50% PV recovery means one retained cell. These are in-sample RNA predictions, not measured FACS performance. Source: condition_gate_metrics.tsv; gate_summary.json.'


def s7(summary):
    frontier=read('gate_pareto.tsv').sort_values('target_recovery')
    fig,grid=canvas(7,'Preservation and SST depletion compete across candidate gates',
        'The practical FGFR2/PTPRS candidate preserves cells; it produces only modest SST-biased purification.',
        'The unconstrained preservation optimum changes very few cells. All candidates were selected and evaluated on the same frozen population.',(1,2))
    ax=fig.add_subplot(grid[0,0]);ax.plot(frontier.sst_contamination*100,frontier.target_recovery*100,color='#A7B2BD',lw=1)
    ax.scatter(frontier.sst_contamination*100,frontier.target_recovery*100,c='#6E92B4',s=35)
    keys=['baseline','phase1_reference','selected','practical_single','experimental_depletion']
    names=['No depletion','Old FAT3/PTPRM','Unconstrained optimum','SLCO3A1 single marker','FGFR2/PTPRS candidate']
    colors=[GRAY,'#8B627E','#937728','#A2784F',COLORS[1]]
    for key,name,color in zip(keys,names,colors):
        g=summary[key];ax.scatter(g['sst_contamination']*100,g['target_recovery']*100,s=95,c=color,edgecolors='white',
                                 label=name+f' (n={int(g["retained_n"]):,})',zorder=5)
    ax.set(xlabel='SST-biased among retained cells (%)',ylabel='PV-biased + dual-high recovery (%)',ylim=(0,105))
    ax.set_title('A  Saved Pareto frontier',loc='left');ax.legend(frameon=False,fontsize=10,loc='lower right')
    ax=fig.add_subplot(grid[0,1]);ys=np.arange(len(keys));width=.32
    for j,(metric,color,label) in enumerate([('pv_recovery',COLORS[0],'PV-biased'),('hybrid_recovery',COLORS[1],'Dual-high')]):
        ax.barh(ys+(j-.5)*width,[summary[k][metric]*100 for k in keys],height=width,color=color,label=label)
    shortnames=['No depletion','Old FAT3/PTPRM','Unconstrained\noptimum','SLCO3A1\nsingle marker','FGFR2/PTPRS\ncandidate']
    ax.set_yticks(ys,shortnames,fontsize=11);ax.invert_yaxis();ax.set_xlim(0,105);ax.set_xlabel('Recovery of starting state (%)')
    ax.set_title('B  Cost in the populations to preserve',loc='left')
    fig.legend(*ax.get_legend_handles_labels(),frameon=False,loc='lower center',bbox_to_anchor=(.5,.065),ncol=2)
    return fig, 'Pareto frontier uses retained SST-biased contamination and combined PV-biased/dual-high recovery. The practical score is sqrt(PV recovery × dual-high recovery) × (1 − SST contamination among retained resolved states); unresolved cells remain in actual retained fractions. Analyst-defined practical depletion requires ≥20% SST-biased removal; sensitivity to 10/20/30/50% is saved. The unconstrained SLC6A6 > 2.9582818 rule removes only 22 cells. Practical single: SLCO3A1 > 0.85467964. Practical pair: FGFR2 > 1.2491154 OR PTPRS > 3.4076235. All thresholds are RNA log1p(CP10K), not protein/fluorescence. A separate inherited ERBB4/CXCR4 positive-retention benchmark retains 1,715 cells and recovers 39.8% PV-biased and 24.9% dual-high; it differs from a new positive-removal sweep. Source: gate_pareto.tsv; gate_summary.json.'


def s8():
    d=read('module_threshold_sensitivity.tsv')
    fig,grid=canvas(8,'Dual-high membership depends on thresholds and gene choice',
        'The primary median-based partition is transparent, with all alternative definitions retained for review.',
        'If a quantile falls at the all-undetected score floor, the entire tied floor is assigned low. No cells are split by arbitrary barcode ordering.',(1,2))
    q=d[d.variant.str.match(r'q\d+_q\d+$')].copy()
    q['pv_q']=q.variant.str.extract(r'q(\d+)_')[0].astype(int);q['sst_q']=q.variant.str.extract(r'_q(\d+)')[0].astype(int)
    matrix=q.pivot(index='pv_q',columns='sst_q',values='hybrid_fraction').loc[[35,50,65],[35,50,65]]*100
    ax=fig.add_subplot(grid[0,0]);im=ax.imshow(matrix,cmap='Greens',vmin=0,vmax=50,aspect='equal')
    for i in range(3):
        for j in range(3):ax.text(j,i,f'{matrix.iloc[i,j]:.1f}%',ha='center',va='center',fontsize=18,color='white' if matrix.iloc[i,j]>32 else INK)
    ax.set_xticks(range(3),['35th','50th','65th']);ax.set_yticks(range(3),['35th','50th','65th'])
    ax.set(xlabel='SST-score quantile cutoff',ylabel='PV-score quantile cutoff',title='A  Fraction of all cells labeled dual-high')
    fig.colorbar(im,ax=ax,fraction=.045,pad=.03,label='Dual-high (%)')
    loo=d[d.variant.str.startswith('without_')]
    ax=fig.add_subplot(grid[0,1]);ax.barh(np.arange(len(loo)),loo.hybrid_fraction*100,color=COLORS[1],alpha=.85)
    ax.set_yticks(np.arange(len(loo)),loo.variant.str.replace('without_','Omit ',regex=False));ax.invert_yaxis()
    ax.axvline(1309/4768*100,color=INK,ls='--',lw=1.4,label='Primary dual-high fraction')
    ax.set(xlabel='Dual-high fraction of all cells (%)',xlim=(0,50),title='B  Leave-one-gene-out sensitivity')
    ax.legend(frameon=False,fontsize=10,loc='lower right')
    return fig, 'All values come from the completed module_threshold_sensitivity.tsv. The grid combines PV and SST quantile cutoffs at 35%, 50% and 65%; gene omission uses the documented re-calculated median thresholds without dropping cells. Primary cutoffs are above both all-undetected floors. Floor-tied alternatives report achieved high fractions, which need not equal nominal quantile masses. Expanded/shared PV, alternative SST, unbounded-score and within-sample depth diagnostics are also saved. The phase1_genes_independent_axes variant uses phase-one genes with phase-two processing, not an exact reproduction of historical phase-one scores. No thresholds were newly optimized for this display.'


def s9():
    fig,grid=canvas(9,'Human flow evidence supports a protein-staining pilot',
        'Developmental evidence supports the score modules; human surface-flow evidence supports testing the candidate panel.',
        'Gene and reagent evidence: ../evidence/DIV90_HYPERGATE_PHASE2_EVIDENCE.md. Fluorescence boundaries require direct experimental calibration.',(2,2))
    items=[('PV-associated program','MEF2C · KCNC1 · KCNC2 · TAC1\n\nDevelopmental and fast-spiking associations.\nPVALB detection is not required.\nMAF / MAFB are shared-regulator sensitivity genes.', '#EEF3F9'),
        ('SST-associated program','SST · SATB1 · NR2F2 · CDK14 · CACNG3\n\nSATB1 also reflects maturation/activity.\nNR2F2 is context dependent, not SST specific.\nNo generic GAD or gate-marker score weight.','#FBF2EA'),
        ('Concrete surface-staining pilot','ERBB4: clone 182818 / MAB11311\nFGFR2: clone 98739 / MAB6843\nPTPRS: goat polyclonal AF3430\n\nHuman flow evidence exists for these candidates.', '#EDF5F0'),
        ('What remains experimental','Extracellular staining in DIV90 neurons.\nDissociation tolerance and cell survival.\nRNA/protein correspondence and thresholds.\nReceptor-binding effects on maturation.\n\nCollect both fractions and an unsorted control.', '#F1F4F7')]
    for i,(title,body,color) in enumerate(items):
        ax=fig.add_subplot(grid[i//2,i%2]);ax.axis('off');box(ax,0,0,.98,.98,title,body,color,13)
    return fig, 'Evidence was audited in the completed phase-two analysis, with no new literature search for this package. Primary developmental sources: Mayer et al. Nature 2018, https://www.nature.com/articles/nature25999; Allaway et al. 2021, https://pmc.ncbi.nlm.nih.gov/articles/PMC9316417/; Mi et al. 2018, https://pmc.ncbi.nlm.nih.gov/articles/PMC6195193/. Manufacturer human flow: ERBB4 clone 182818, extracellular immunogen Glu26–Arg649, https://www.rndsystems.com/products/human-erbb4-her4-antibody-182818_mab11311; FGFR2 clone 98739, human Kato III flow and broad receptor isoform coverage, https://www.rndsystems.com/products/human-fgfr2-antibody-98739_mab6843. PTPRS AF3430 extracellular immunogen Glu30–Gly1260: https://resources.rndsystems.com/pdfs/datasheets/af3430.pdf; human nonpermeabilized pDC/PBMC surface-flow evidence in Bunin et al. Immunity 2015, https://pmc.ncbi.nlm.nih.gov/articles/PMC4547994/. AF3430 supplier links this primary study. ERBB4/FGFR2 candidates are both mouse antibodies, so a common anti-mouse secondary cannot distinguish them; use compatible direct conjugates or validated selective detection. PTPRS crosslinking can alter cell function. FAT3/PTPRM and SLCO3A1 lack equivalent verified neuronal live-sort evidence. Neither candidate panel nor RNA thresholds constitute a validated protein FACS protocol.'


def s10(summary):
    validation=json.loads((ROOT/'provenance/final_validation.json').read_text())
    fig,grid=canvas(10,'The completed RNA analysis defines a sort-and-mature experiment',
        'Completed analysis and new biological experiment are separate stages; only figure rendering is performed in this package.',
        'The predicted gate applies to the frozen cortical LHX6+/ERBB4+ subset. A live ERBB4-only entry gate does not automatically reproduce cortical/LHX6 transcript selection.',(1,1))
    ax=fig.add_subplot(grid[0,0]);ax.axis('off')
    stages=[(.00,.62,.29,.32,'1  Frozen entry population','4,768 recovered cortical cells\nLHX6+ / ERBB4+ RNA\nOriginal identities and UMAP preserved'),
        (.355,.62,.29,.32,'2  Independent programs','PV and SST scores, separately\nFour operational quadrants\nGene/QC/threshold sensitivities'),
        (.71,.62,.27,.32,'3  Surface RNA discovery','314 expressed allowed features\nSimple depletion rules\nPV and dual-high preservation'),
        (.00,.12,.29,.32,'4  Candidate RNA split','FGFR2-high OR PTPRS-high\nRemove captured cells\nRetain the complement'),
        (.355,.12,.29,.32,'5  Protein pilot','Viable ERBB4+ entry\nMeasure FGFR2 / PTPRS staining\nCollect BOTH fractions + control'),
        (.71,.12,.27,.32,'6  Identical maturation','PV / SST protein\nMature transcriptomic state\nFast-spiking electrophysiology')]
    for x,y,w,h,title,body in stages:box(ax,x,y,w,h,title,body,'#EDF5F0' if title.startswith(('5','6')) else '#F1F4F7',12.5)
    for a,b in [((.30,.79),(.345,.79)),((.655,.79),(.70,.79)),((.30,.29),(.345,.29)),((.655,.29),(.70,.29))]:
        ax.annotate('',xy=b,xytext=a,xycoords='axes fraction',arrowprops={'arrowstyle':'-|>','color':GRAY,'lw':1.7})
    ax.text(.49,.51,'Completed RNA results  →  prospective protein and fate test',ha='center',fontsize=15,weight='bold',transform=ax.transAxes)
    return fig, f'Frozen eligibility: recovered cortical_only membership, LHX6 > 0 and ERBB4 > 0 RNA, excluding original stressed clusters 6/7. Source labels, coordinates and phase-one results are preserved. {summary["actual_R_hypergate_fits"]:,} actual R Hypergate fits were audited separately from Python threshold sweeps; count conservation was checked for {validation["total_python_candidate_rows_count_conservation_checked"]:,} saved candidate rows, with {validation["random_candidate_masks_independently_reproduced"]:,} independently recomputed random masks. Score genes are excluded from surface-gate discovery. The analyst-defined practical comparison requires at least 20% SST-biased removal. Selected RNA rule: remove FGFR2 > 1.2491154 OR PTPRS > 3.4076235; complement retained. All feature selection and evaluations use the same cells, without held-out sample validation. Prospective experiment: viable ERBB4+ entry, FGFR2/PTPRS protein pilot, retain and remove fractions plus unsorted ERBB4+ control, identical subsequent culture, later protein/transcriptome/physiology and absolute-survival comparisons. Sources: gate_summary.json; state_parameters.json; provenance/input_validation.json; provenance/final_validation.json; evidence note.'


def run(cells=None, out=OUT, summary=None, params=None):
    out=Path(out).resolve()
    if out != OUT.resolve():raise ValueError(f'Supplement output must be {OUT}')
    configure();directory=out/'supplement';directory.mkdir(parents=True,exist_ok=True)
    (out/'qc/supplement_previews').mkdir(parents=True,exist_ok=True)
    if cells is None:cells=pd.read_csv(ROOT/'cells.tsv.gz',sep='\t')
    if summary is None:summary=json.loads((ROOT/'gate_summary.json').read_text())
    if params is None:params=json.loads((ROOT/'state_parameters.json').read_text())
    hybrid=json.loads((ROOT/'hybrid_validation.json').read_text())
    assert len(cells)==4768 and cells.cell_id.is_unique
    assert cells.state.value_counts().reindex(STATES).tolist()==[1076,1309,1075,1308]
    builders=[('S1_score_diagnostics',lambda:s1(cells,params)),('S2_dual_high_gene_support',lambda:s2(cells,params,hybrid)),
        ('S3_conditional_null',lambda:s3(hybrid)),('S4_CV_MW_by_cell_line',s4),
        ('S5_old_gate_condition_dependence',lambda:s5(summary)),('S6_candidate_condition_performance',lambda:s6(summary)),
        ('S7_gate_pareto',lambda:s7(summary)),('S8_threshold_robustness',s8),('S9_marker_reagent_evidence',s9),
        ('S10_computational_workflow',lambda:s10(summary))]
    records=[]
    with PdfPages(directory/'supplemental_figures.pdf') as pdf:
        pdf.infodict().update(Title='DIV90 PI figure package: technical supplement',Author='DIV90 analysis',Subject='Completed phase-two results, no biological refitting')
        for stem, builder in builders:
            fig,caption=builder()
            # Keep all plot marks and text vector in PDF/SVG.
            for artist in fig.findobj():
                if hasattr(artist,'get_rasterized') and artist.get_rasterized():artist.set_rasterized(False)
            pdf.savefig(fig)
            fig.savefig(directory/(stem+'.svg'))
            fig.savefig(directory/(stem+'.png'),dpi=600)
            fig.savefig(out/'qc/supplement_previews'/(stem+'.png'),dpi=110)
            plt.close(fig)
            records.append({'figure':stem,'caption':caption,'svg':'supplement/'+stem+'.svg',
                            'png':'supplement/'+stem+'.png','png_dpi':600,'width_inches':15,'height_inches':9})
            print('Rendered '+stem,flush=True)
    (directory/'SUPPLEMENTAL_CAPTIONS.md').write_text('# DIV90 technical supplement\n\nRead the four main figures first. These pages expose uncertainty, culture dependence and computational provenance. Exact thresholds are RNA values, not fluorescence gates.\n\n'+'\n\n'.join('## '+r['figure']+'\n\n'+r['caption'] for r in records)+'\n')
    metadata={'source_root':str(ROOT),'output_root':str(out),'source_cells_n':len(cells),
        'new_biological_fits':False,'pages':10,'figures':records,'consolidated_pdf':'supplement/supplemental_figures.pdf',
        'palette':dict(zip(LABELS,COLORS)),'condition_display':'H9, 79B, 2E shown separately; no pooling as primary display',
        'source_gate':summary['experimental_depletion']['gate_label']}
    (directory/'supplement_manifest.json').write_text(json.dumps(metadata,indent=2)+'\n')
    return metadata


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,default=OUT);args=parser.parse_args()
    run(out=args.out)
