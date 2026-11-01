"""Apply identical frozen gates to WT females as a composition negative control."""
import csv
import gzip
import hashlib
import json
import math
from collections import Counter
from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import textwrap

ROOT=Path(__file__).resolve().parents[1]
BASE=Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results')
OUT=RUNTIME_ROOT/'results/wt_female_probe_gate_control_20261031'
OUT.mkdir(parents=True,exist_ok=True)
MODEL=BASE/'step_05_pcdh19_logistic_regression_baseline/count_informed_vs_binary_validation/step_05_count_informed_full_fit_coefficients.tsv'
HELD=MODEL.parent/'step_05_binary_vs_count_informed_held_out_predictions.tsv'
ASSIGN=BASE/'hicat/03_full_data_consensus_benchmark/hicat_consensus_review_20260909_142348_efd9e13d/outputs/tables/all_cell_assignments.tsv.gz'
CALLS=BASE/'step_07_pcdh19_het_female_wt_ko_like_classification/step_07_het_female_cell_classifications.tsv'
def read(p):
    with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:yield from csv.DictReader(f,delimiter='\t')
def write(name,rows):
    with (OUT/name).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
assert hashlib.sha256(MODEL.read_bytes()).hexdigest()=='2596c0c2d7a6fca539b7af3b2ea19792dc884f7c5a874f3c72322d2e6ea3215a'
coef={r['term']:float(r['coefficient']) for r in read(MODEL)}
def probability(a,b,c):
    z=coef['intercept']+coef['A_detected']*(a>0)+coef['B_UMI']*b+coef['C_UMI']*c
    return round(1/(1+math.exp(-z)),12)
def classify(p,a,b,c):
    if a+b+c==0:return 'Uncalled_000'
    if p<=.301037619832:return 'WT_like'
    if p>=.911554020713:return 'KO_like'
    return 'Uncertain'
assignment={r['cell_id']:r['diagnostic_cluster'] for r in read(ASSIGN)}
clusters=sorted(set(assignment.values()));counts=Counter();totals=Counter();excluded=Counter()
sources=[MODEL,HELD,ASSIGN,CALLS]
def add(method,s,barcode,label):
    cell=s+'_'+barcode
    if cell not in assignment:excluded[method,s]+=1;return
    c=assignment[cell]
    for g in [label,'all']:
        counts[method,s,g,c]+=1;totals[method,s,g]+=1
for i in [4,5,6]:
    s='15662-JZ-'+str(i);p=BASE/'pcdh19_probe_audit/per_sample'/s/'pcdh19_probe_patterns.tsv';sources.append(p)
    for r in read(p):
        a,b,c=[int(r[x]) for x in ['A_UMI','B_UMI','C_UMI']]
        add('WT_F_same_full_model',s,r['barcode'],classify(probability(a,b,c),a,b,c))
validated=0
for r in read(CALLS):
    a,b,c=[int(r[x]) for x in ['A_UMI','B_UMI','C_UMI']]
    label=classify(probability(a,b,c),a,b,c)
    assert label==r['final_classification'];validated+=1
    add('HET_F_same_full_model',r['biological_sample_id'],r['cell_barcode'],label)
assert validated==101102
for r in read(HELD):
    if r['design_group']!='WT_F':continue
    a,b,c=[int(r[x]) for x in ['A_UMI','B_UMI','C_UMI']]
    add('WT_F_sample_held_out',r['biological_sample_id'],r['cell_barcode'],classify(float(r['count_predicted_ko_probability']),a,b,c))
methods=['WT_F_same_full_model','HET_F_same_full_model','WT_F_sample_held_out']
samples={m:['15662-JZ-'+str(i) for i in ([7,8,9] if m.startswith('HET') else [4,5,6])] for m in methods}
coverage=[];proportions=[]
for m in methods:
    for s in samples[m]:
        assert totals[m,s,'all']==sum(totals[m,s,g] for g in ['WT_like','KO_like','Uncertain','Uncalled_000'])
        coverage.append(dict(method=m,sample=s,retained_cells=totals[m,s,'all'],WT_like=totals[m,s,'WT_like'],KO_like=totals[m,s,'KO_like'],
                             uncertain=totals[m,s,'Uncertain'],uncalled=totals[m,s,'Uncalled_000'],
                             KO_like_percent_of_all=100*totals[m,s,'KO_like']/totals[m,s,'all']))
        for g in ['WT_like','KO_like','Uncertain','Uncalled_000','all']:
            for c in clusters:
                n=counts[m,s,g,c];den=totals[m,s,g]
                proportions.append(dict(method=m,sample=s,gate=g,cluster=c,count=n,denominator=den,percent=100*n/den if den else ''))
write('gate_coverage.tsv',coverage);write('sample_gated_proportions.tsv',proportions)
names={r['cluster']:r['working_name'] for r in read(RUNTIME_ROOT/'results/composition_review_20261031/full_consensus_working_names.tsv')}
effects=[]
for c in clusters:
    r=dict(cluster=c,working_name=names.get(c,'Unannotated'))
    for m in methods:
        ds=[]
        for s in samples[m]:
            if totals[m,s,'WT_like'] and totals[m,s,'KO_like']:
                ds.append(100*counts[m,s,'KO_like',c]/totals[m,s,'KO_like']-100*counts[m,s,'WT_like',c]/totals[m,s,'WT_like'])
        r[m+'_mean_shift_pp']=sum(ds)/len(ds) if ds else ''
        r[m+'_evaluable_samples']=len(ds)
        r[m+'_all_same_direction']=bool(ds and all(x>0 for x in ds) or ds and all(x<0 for x in ds))
    effects.append(r)
write('cluster_shift_control_comparison.tsv',effects)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.linewidth':.7,'pdf.fonttype':42,'svg.fonttype':'none'})
chosen=['C0001','C0011','C0007','C0005'] # Fixed from previous HET ranking, before examining controls.
for control in ['WT_F_same_full_model','WT_F_sample_held_out']:
    fig,axes=plt.subplots(2,4,figsize=(9,6.2))
    for row,m in enumerate([control,'HET_F_same_full_model']):
        for col,c in enumerate(chosen):
            ax=axes[row,col]
            for j,s in enumerate(samples[m]):
                den0,den1=totals[m,s,'WT_like'],totals[m,s,'KO_like']
                a=100*counts[m,s,'WT_like',c]/den0 if den0 else None
                b=100*counts[m,s,'KO_like',c]/den1 if den1 else None
                if a is not None and b is not None:ax.plot([-.1+j*.1,.9+j*.1],[a,b],color='#bbbbbb',lw=.8,zorder=1)
            for x,g,color in [(0,'WT_like','#549E89'),(1,'KO_like','#A678AC')]:
                vals=[]
                for j,s in enumerate(samples[m]):
                    if not totals[m,s,g]:continue
                    v=100*counts[m,s,g,c]/totals[m,s,g];vals.append(v)
                    ax.scatter(x-.1+j*.1,v,s=25,color=color,edgecolors='white',linewidths=.4,zorder=3)
                if vals:ax.plot([x-.22,x+.22],[np.mean(vals),np.mean(vals)],color=color,lw=1.7)
                else:ax.text(x,.5,'No cells',transform=ax.get_xaxis_transform(),ha='center',rotation=90)
            ax.set_title('\n'.join(textwrap.wrap(names.get(c,c),19))+'\n'+c,fontsize=8,pad=10)
            ax.set_xticks([0,1]);ax.set_xticklabels(['WT-like','KO-like'],fontsize=7)
            ax.set_xlim(-.5,1.5);ax.set_ylim(bottom=0);ax.set_ylabel('Cells (% of gated subset)',fontsize=7)
            for side in ['top','right']:ax.spines[side].set_visible(False)
        label='WT females: same frozen model' if row==0 and control=='WT_F_same_full_model' else 'WT females: sample-held-out predictions' if row==0 else 'HET females: same frozen model'
        denominators=', '.join(str(totals[m,s,'KO_like']) for s in samples[m])
        fig.text(.08,.96 if row==0 else .49,label+' (KO-like cells/sample: '+denominators+')',fontsize=9,fontweight='bold')
    for col in range(4):
        ceiling=max(axes[0,col].get_ylim()[1],axes[1,col].get_ylim()[1])
        for row in range(2):axes[row,col].set_ylim(0,ceiling)
    fig.subplots_adjust(left=.08,right=.99,bottom=.12,top=.81,wspace=.75,hspace=1.1)
    fig.text(.08,.04,'Each line joins subsets from one sample. KO-like calls in known WT females are false positives.',fontsize=7)
    for ext in ['pdf','svg','png']:fig.savefig(OUT/(control+'_versus_HET.'+ext),dpi=400,bbox_inches='tight',facecolor='white')
    plt.close(fig)
checks=dict(model_sha256_matches_frozen_pin=True,HET_calls_reproduced=validated,all_subset_counts_partition_retained_cells=True)
(OUT/'validation.json').write_text(json.dumps(checks,indent=2))
(OUT/'provenance.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},indent=2))
report=['# WT-female probe-gate negative control','',
        'Identical frozen full-fit coefficients and thresholds were applied to WT-female probes, without refitting. All 101,102 historical HET classifications were independently reproduced exactly. A second WT analysis uses saved leave-one-sample-out probabilities to examine sensitivity to training inclusion. Thresholds themselves were previously selected on these controls; this is a diagnostic negative control, not independent validation.',
        '', 'All cells are joined to the existing consensus by exact sample/barcode; primary-QC exclusions are retained separately. WT females are genetically known WT: any KO-like gate membership is classifier false-positive evidence, not a discovered KO population. Different cluster proportions between gates are therefore possible even in WT tissue. The negative control tests whether the earlier HET gate enrichment is specific; it does not assume a null outcome.',
        '', 'Rare KO-like subsets can produce unstable percentages. Missing gates are not treated as zero composition. The four displayed clusters were fixed from the previous HET effect-size ranking before inspecting WT-control outcomes.','', '## Gate coverage','']
for r in coverage:report.append('- {method}; {sample}: WT-like {WT_like}, KO-like {KO_like}; KO-like {KO_like_percent_of_all:.4f}% of retained cells'.format(**r))
report += ['','## Mean within-sample KO-like minus WT-like proportion shift','',
           '| Cluster | Working name | WT full model pp | HET full model pp | WT held-out pp |','|---|---|---:|---:|---:|']
for c in chosen:
    r=next(r for r in effects if r['cluster']==c)
    vals=[r[m+'_mean_shift_pp'] for m in methods]
    report.append('| '+c+' | '+r['working_name']+' | '+' | '.join('%+.2f'%v if v!='' else 'not estimable' for v in vals)+' |')
(OUT/'REPORT.md').write_text('\n'.join(report)+'\n')
print('\n'.join(report[9:]))
