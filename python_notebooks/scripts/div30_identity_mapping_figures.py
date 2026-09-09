#!/usr/bin/env python3
"""Render frozen-neighborhood identity diagnostics; reads completed tables only."""
import os
from pathlib import Path
ROOT=Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_antecedent_identity_v1/mapping')
os.environ['MPLCONFIGDIR']=str(ROOT/'cache/matplotlib')
for _name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[_name]='2'
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':11,'axes.labelsize':10,
                     'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
COLORS={'strict_PV':'#27679F','dual_high':'#8857A3','SST':'#DA8047','unresolved':'#9BA1A6'}
NAMES={'strict_PV':'Strict PV-biased','dual_high':'Dual-high (PV + SST)','SST':'SST-biased','unresolved':'Unresolved'}
GROUPS=['A_target_retained','B_target_lost','C_non_target_retained','D_non_target_removed']
GROUP_LABELS=['A Target + retained','B Target + lost','C Non-target + retained','D Non-target + removed']
GROUP_COLORS=['#08519c','#6baed6','#cb681f','#e2bc96']
CAPTIONS={}

def read(name):return pd.read_csv(ROOT/name,sep='\t')
def ecdf(ax,v,**kwargs):
    a=np.sort(np.asarray(v));ax.plot(a,np.arange(1,len(a)+1)/len(a),**kwargs)
def save(fig,stem,caption):
    for ext in ['png','pdf','svg']:
        p=ROOT/'figures'/ext;p.mkdir(parents=True,exist_ok=True)
        fig.savefig(p/f'{stem}.{ext}',dpi=600 if ext=='png' else None,bbox_inches='tight',facecolor='white')
    plt.close(fig);CAPTIONS[stem]=caption

def forward():
    d=read('forward_cell_neighbors.tsv.gz');s=read('forward_group_distributions.tsv');e=read('forward_retained_removed_effects.tsv')
    fig,axs=plt.subplots(2,3,figsize=(16,9.6),gridspec_kw={'height_ratios':[1,1.05]})
    fig.suptitle('Retained DIV30 cells resemble DIV90 dual-high territory, with fewer strict-PV neighbors',fontsize=15,y=.995)
    ax=axs[0,0];order=['retained','removed']+GROUPS;bottom=np.zeros(len(order));yy=np.arange(len(order))
    for state in COLORS:
        v=s.query("method=='archived_unbalanced' and state==@state").set_index('group').loc[order,'mean'].to_numpy()
        ax.barh(yy,100*v,left=100*bottom,color=COLORS[state],label=NAMES[state]);bottom+=v
    ax.set_yticks(yy,['All retained','All removed']+GROUP_LABELS);ax.invert_yaxis();ax.set_xlim(0,100)
    ax.set_xlabel('Mean empirical DIV90 neighbor fraction (%)');ax.set_title('A  Full frozen 30-neighbor composition',loc='left')
    ax.legend(loc='lower left',bbox_to_anchor=(-.02,-.46),ncol=2,frameon=False,fontsize=8)
    for ax,state,panel in [(axs[0,1],'strict_PV','B'),(axs[0,2],'dual_high','C'),(axs[1,0],'SST','D'),(axs[1,1],'unresolved','E')]:
        for group,color,label in zip(GROUPS,GROUP_COLORS,GROUP_LABELS):
            ecdf(ax,100*d.loc[d.four_group.eq(group),'archived_unbalanced__'+state],color=color,label=label,lw=1.8)
        ecdf(ax,100*d.loc[d.gate_retained,'archived_unbalanced__'+state],color='black',label='All retained',lw=1.3)
        ecdf(ax,100*d.loc[~d.gate_retained,'archived_unbalanced__'+state],color='black',label='All removed',lw=1.3,ls='--')
        ax.set_xlabel('Empirical DIV90 neighbor fraction (%)');ax.set_ylabel('Cumulative fraction of DIV30 cells')
        ax.set_title(f'{panel}  {NAMES[state]}: entire distributions',loc='left');ax.set_ylim(0,1)
    axs[1,0].legend(frameon=False,fontsize=8,loc='lower right')
    ax=axs[1,2]
    for keep,color,label in [(True,'#08519c','Retained'),(False,'#cb681f','Removed')]:
        ecdf(ax,d.loc[d.gate_retained.eq(keep),'archived_unbalanced__distance'],color=color,label=label,lw=2)
    raw=pd.read_csv(ROOT.parent.parent/'div30_pv_antecedent_v1/tables/reference_oof_predictions.tsv.gz',sep='\t')
    cut=raw.query("grouping=='sample'").neighbor_distance.quantile(.95)
    ax.axvline(cut,color='gray',ls='--',label='Frozen DIV90 held-sample 95%')
    ax.set_xlabel('Mean distance to 30 DIV90 neighbors');ax.set_ylabel('Cumulative fraction of DIV30 cells')
    ax.set_title('F  Retained cells are closer to the reference',loc='left');ax.legend(frameon=False,fontsize=8)
    fig.tight_layout(rect=(0,.02,1,.965),w_pad=2.8,h_pad=3)
    save(fig,'01_forward_neighbor_identity',
      'Frozen archived DIV30-to-DIV90 neighbors, k=30. A shows all four mutually exclusive DIV90 states, including dual-high PV/SST state; B–E show every cell through empirical cumulative distributions, for retained/removed and four target-by-gate groups. F uses the original mean-neighbor distance and held-sample DIV90 95th percentile threshold. Strict-PV fractions decrease in retained DIV30 cells, whereas dual-high fractions increase; the effect direction also holds in every sample (table forward_retained_removed_effects.tsv). These are descriptive effects, not independent-replicate significance. Empirical neighbor fractions are not calibrated probabilities or lineage assignments. The DIV90 reference is restricted to 4,768 LHX6+/ERBB4+ cortical cells.')

def reverse():
    d=read('reverse_cell_neighbors.tsv.gz');s=read('reverse_query_group_distributions.tsv');h=read('reverse_neighbor_hubness.tsv');ef=read('reverse_strictPV_vs_SST_effects.tsv')
    states=list(COLORS);fig,axs=plt.subplots(2,3,figsize=(16,9.6))
    fig.suptitle('Reciprocal mapping favors dual-high cells; strict-PV convergence weakens after sample balancing',fontsize=15,y=.995)
    ax=axs[0,0];xx=np.arange(4)
    for method,offset,label in [('unbalanced',-.17,'30 nearest overall'),('equal_sample_neighbors',.17,'5 nearest per sample')]:
        q=s.query('method==@method').set_index('query_group');v=q.loc[states,'mean']
        ax.bar(xx+offset,100*v,width=.32,color=[COLORS[x] for x in states],alpha=1 if method=='unbalanced' else .5,label=label)
        ax.axhline(q.iloc[0].DIV30_gate_yield_baseline*100,color='gray',ls='--' if method=='unbalanced' else ':',lw=1)
    ax.set_xticks(xx,[NAMES[x].replace(' (PV + SST)','') for x in states],rotation=20,ha='right');ax.set_ylabel('Mean retained DIV30 neighbor fraction (%)')
    ax.set_ylim(0,40)
    ax.set_title('A  Return to gate-retained DIV30 cells',loc='left');ax.legend(frameon=False,fontsize=8)
    ax.text(.02,.98,'Dashed: all-cell gate yield (17.4%)\nDotted: equal-sample gate yield (16.2%)',transform=ax.transAxes,va='top',fontsize=8)
    ax=axs[0,1]
    state_map={'strict_PV':'PV-biased','dual_high':'PV/SST hybrid','SST':'SST-biased','unresolved':'unresolved/immature'}
    for st in states:ecdf(ax,100*d.loc[d.state.eq(state_map[st]),'unbalanced__retained_fraction'],color=COLORS[st],label=NAMES[st],lw=1.8)
    ax.set_xlabel('Empirical retained DIV30 neighbor fraction (%)');ax.set_ylabel('Cumulative fraction of DIV90 cells')
    ax.set_title('B  Full reciprocal neighbor distributions',loc='left');ax.legend(frameon=False,fontsize=8)
    ax=axs[0,2];names=['strict_PV','strict_PV_within_PV_score_top50','strict_PV_within_PV_score_top25','strict_PV_within_PV_score_top10']
    for method,label,color in [('unbalanced','30 nearest overall','#27679F'),('equal_sample_neighbors','5 nearest per sample','#8857A3')]:
        v=s.query('method==@method').set_index('query_group').loc[names]
        ax.plot(np.arange(4),v['mean']*100,'o-',color=color,label=label)
    ax.set_xticks(np.arange(4),['All strict PV\nn=1,076','Top 50%\nn=538','Top 25%\nn=269','Top 10%\nn=108']);ax.set_ylabel('Retained DIV30 neighbor fraction (%)')
    ax.set_title('C  Stricter continuous-PV subsets: nonmonotonic',loc='left');ax.legend(frameon=False,fontsize=8)
    ax=axs[1,0]
    for method,label,marker in [('unbalanced','Raw degree correction','o'),('equal_sample_neighbors','Sample balanced + degree correction','s')]:
        q=s.query('method==@method').set_index('query_group').loc[states]
        ax.scatter(100*q['mean'],100*q.inverse_all_query_degree_weighted_retained_fraction,color=[COLORS[x] for x in states],marker=marker,s=60,label=label)
    ax.plot([0,40],[0,40],color='gray',ls='--');ax.set_xlim(0,40);ax.set_ylim(0,40)
    ax.set_xlabel('Unweighted edge retention (%)');ax.set_ylabel('Inverse-global-degree weighted retention (%)')
    ax.set_title('D  Hubness changes apparent convergence',loc='left');ax.legend(frameon=False,fontsize=8)
    ax=axs[1,1];q=h.query("method=='unbalanced'").set_index('query_group').loc[states]
    ax.bar(np.arange(4),q.edge_share_top_1pct_all_DIV30*100,color=[COLORS[x] for x in states])
    ax.set_xticks(np.arange(4),[NAMES[x].replace(' (PV + SST)','') for x in states],rotation=20,ha='right')
    ax.set_ylabel('Edges received by most-used 1% of DIV30 cells (%)');ax.set_ylim(0,100)
    ax.set_title('E  A small set of counterparts receives many edges',loc='left')
    ax=axs[1,2];cal=read('reverse_DIV30_held_sample_distance_calibration.tsv').mean_30nn_distance
    ecdf(ax,cal,color='#7A7A7A',label='DIV30 → other DIV30 samples',lw=2)
    ecdf(ax,d.loc[d.state.eq('PV-biased'),'unbalanced__distance'],color=COLORS['strict_PV'],label='DIV90 strict PV → DIV30',lw=2)
    ecdf(ax,d.loc[d.state.eq('PV/SST hybrid'),'unbalanced__distance'],color=COLORS['dual_high'],label='DIV90 dual-high → DIV30',lw=2)
    ax.axvline(cal.quantile(.95),color='gray',ls='--');ax.set_xlabel('Mean distance to 30 DIV30 neighbors');ax.set_ylabel('Cumulative fraction of query cells')
    ax.set_title('F  Cross-time distances exceed\n    the internal DIV30 threshold',loc='left');ax.legend(frameon=False,fontsize=8)
    fig.tight_layout(rect=(0,.02,1,.965),w_pad=3,h_pad=3)
    save(fig,'02_reverse_neighbor_convergence',
      'DIV90 queries mapped to 90,631 DIV30 cells in the same verified frozen DIV90 PCA space. A compares unrestricted 30-neighbor mapping with five neighbors from each of six DIV30 samples; horizontal lines are relevant recipient gate yields. B includes all four reference states. C preserves the frozen strict-PV definition and adds within-strict-PV continuous-score top50/25/10 sensitivity subsets. D downweights each edge by its DIV30 recipient in-degree across all DIV90 queries, diagnosing hubs without using target labels. E reports how concentrated edges are in the most-used 1% of all DIV30 cells. F compares reciprocal distances with a stratified 3,000-cell DIV30 held-sample distance calibration. All strict-PV queries exceed that internal 95th percentile; nearest counterparts therefore do not imply confident equivalence or ancestry. Forward and reverse mapping share representation and training/reference information, so reciprocal convergence is not independent validation. Equal sample representation addresses sample density but not all within-sample state density.')

def temporal():
    e=read('temporal_marker_program_effects.tsv');fig,axs=plt.subplots(2,2,figsize=(15,10.5),gridspec_kw={'width_ratios':[1,1.2]})
    fig.suptitle('Surface separation depends on age, comparison population and target definition',fontsize=15,y=.995)
    contexts=['all','postmitotic','LHX6_ERBB4_positive','postmitotic_LHX6_ERBB4_positive','generic_context_residual']
    clabels=['Whole analysis population','Operational postmitotic only','LHX6+ / ERBB4+ only','Postmitotic + LHX6+ / ERBB4+','Adjusted for generic context']
    for ax,g,panel in zip([axs[0,0],axs[0,1],axs[1,0]],['GRIA2','OPCML','NOTCH1'],['A','B','C']):
        for stage,off,color,marker in [('DIV30',-.1,'#27679F','o'),('DIV90',.1,'#DA8047','s')]:
            q=e.query('stage==@stage and marker==@g').set_index('context').loc[contexts]
            ax.plot(q.hedges_g,np.arange(len(contexts))+off,marker,ms=6,color=color,label=stage)
        ax.axvline(0,color='gray',lw=.8);ax.set_yticks(np.arange(len(contexts)),clabels);ax.invert_yaxis()
        ax.set_xlabel('Hedges g: target versus comparison');ax.set_title(f'{panel}  {g}: effect size across comparable contexts',loc='left');ax.legend(frameon=False)
    ax=axs[1,1];programs=['pv_target_component','sst_target_component','mge','PV_context_not_primary_target','gaba_progression','maturation','progenitor','cycle']
    columns=[('DIV30','all'),('DIV30','LHX6_ERBB4_positive'),('DIV30','postmitotic_LHX6_ERBB4_positive'),('DIV90','all'),('DIV90','postmitotic')]
    lookup=e.set_index(['stage','context','marker']).hedges_g
    vals=np.array([[lookup.loc[(st,ct,'program_'+program)] for st,ct in columns] for program in programs])
    im=ax.imshow(vals,cmap='RdBu_r',vmin=-2,vmax=2,aspect='auto')
    ax.set_yticks(np.arange(len(programs)),['PV target components*','SST target components*','MGE','PV contextual genes†','GABA / progression','Maturation','Progenitor','Cell cycle'])
    ax.set_xticks(np.arange(len(columns)),['DIV30\nall','DIV30\nLHX6+/ERBB4+','DIV30\npostmitotic+entry','DIV90\nall','DIV90\npostmitotic'],rotation=20,ha='right')
    for i in range(vals.shape[0]):
        for j in range(vals.shape[1]):ax.text(j,i,f'{vals[i,j]:.2f}',ha='center',va='center',fontsize=8,color='white' if abs(vals[i,j])>1.1 else 'black')
    ax.set_title('D  Programs distinguish different comparisons',loc='left');fig.colorbar(im,ax=ax,label='Hedges g',shrink=.75)
    fig.text(.03,.017,'DIV30: antecedent target vs non-target. DIV90: strict PV-biased vs SST-biased within the frozen LHX6+/ERBB4+ population.\n* Target-related components are descriptive, not independent validation. † MAF/MAFB/ELMO1/IGFBP4 is contextual, not a validated PV-specific program.',fontsize=9)
    fig.tight_layout(rect=(0,.065,1,.965),w_pad=3.5,h_pad=3)
    save(fig,'03_temporal_marker_program_effects',
      'Individual surface-marker Hedges g (A–C) and within-age gene-z program effects (D). Positive values indicate higher expression in the DIV30 antecedent target than non-target or DIV90 strict PV than SST. DIV30 starts with all 90,631 cells; DIV90 is restricted to 4,768 LHX6+/ERBB4+ cortical cells before state comparison. Restricting DIV30 to the same entry markers and restricting both ages to operationally postmitotic cells makes population context closer but does not equalize target definitions. This operational postmitotic rule uses inherited RNA cell-cycle phase and stage-relative maturation/progenitor scores; it is not a measurement of division. Only 194 strict-PV and 224 SST cells meet it at DIV90. Generic-context residuals remove linear associations with stage-specific progenitor, neuronal maturation, GABA progression, cycle, phase, log1p RNA complexity and sample. Scores using target components remain tautologically related to state labels. TUBB3 is absent from both H5AD gene indexes; its legacy zero contribution is preserved for generic-panel comparability and no TUBB3 detection claim is made. RNA normalization/depth and cross-sectional samples limit attribution to maturation itself; these effect sizes are descriptive and do not establish independent significance.')

if __name__=='__main__':
    forward();reverse();temporal()
    (ROOT/'figure_captions.md').write_text('\n\n'.join('## '+k+'\n\n'+v for k,v in CAPTIONS.items())+'\n')
    (ROOT/'figures/figure_captions.json').write_text(json.dumps([{'figure':k,'caption':v} for k,v in CAPTIONS.items()],indent=2)+'\n')
    (ROOT/'figure_manifest.json').write_text(json.dumps({k:{'caption':v,'files':{e:f'figures/{e}/{k}.{e}' for e in ['png','pdf','svg']}} for k,v in CAPTIONS.items()},indent=2)+'\n')
    print('Three mapping identity plates exported as PDF/SVG/600dpiPNG.')
