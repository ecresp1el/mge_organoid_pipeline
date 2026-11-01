"""Join frozen probe calls to exact consensus cell IDs and compare gated composition."""
import csv
import gzip
import itertools
import json
import hashlib
from collections import Counter
from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
BASE=Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results')
OUT=RUNTIME_ROOT/'results/probe_gated_composition_20261031'
OUT.mkdir(parents=True,exist_ok=True)
ASSIGN=BASE/'hicat/03_full_data_consensus_benchmark/hicat_consensus_review_20260909_142348_efd9e13d/outputs/tables/all_cell_assignments.tsv.gz'
CALLS=BASE/'step_07_pcdh19_het_female_wt_ko_like_classification/step_07_het_female_cell_classifications.tsv'
NAMES=RUNTIME_ROOT/'results/composition_review_20261031/full_consensus_working_names.tsv'
def read(path):
    with (gzip.open(path,'rt') if path.suffix=='.gz' else path.open()) as f:
        yield from csv.DictReader(f,delimiter='\t')
def write(name,rows):
    with (OUT/name).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def avg(x):return sum(x)/len(x)
names={r['cluster']:r['working_name'] for r in read(NAMES)}
assignment={r['cell_id']:r['diagnostic_cluster'] for r in read(ASSIGN)}
assert len(assignment)==446349
clusters=sorted(set(assignment.values()))
counts,totals,coverage=Counter(),Counter(),Counter()
sources=[ASSIGN,CALLS,NAMES]
samples=['15662-JZ-'+str(i) for i in range(4,10)]
def add(sample,barcode,gates):
    cell=sample+'_'+barcode
    coverage[sample,'input']+=1
    if cell not in assignment:
        coverage[sample,'not_in_consensus']+=1
        return
    coverage[sample,'retained']+=1
    for gate in gates:
        counts[sample,gate,assignment[cell]]+=1;totals[sample,gate]+=1
for s in samples[:3]:
    path=BASE/'pcdh19_probe_audit/per_sample'/s/'pcdh19_probe_patterns.tsv';sources.append(path)
    for r in read(path):
        gates=['all']
        if int(r['Pcdh19_total_UMI'])>0:gates.append('any_probe_positive')
        if int(r['A_UMI'])>0:gates.append('A_positive')
        add(s,r['barcode'],gates)
for r in read(CALLS):
    gates=['all',r['final_classification']]
    if r['probe_pattern']!='000':gates.append('any_probe_positive')
    if int(r['A_UMI'])>0:gates.append('A_positive')
    add(r['biological_sample_id'],r['cell_barcode'],gates)
cov=[dict(sample=s,input_cells=coverage[s,'input'],consensus_cells=coverage[s,'retained'],
          excluded_by_primary_qc=coverage[s,'not_in_consensus'],gate=g,gate_cells=n)
     for (s,g),n in sorted(totals.items())]
write('gate_coverage.tsv',cov)
props=[]
for (s,g),total in sorted(totals.items()):
    for c in clusters:
        props.append(dict(sample=s,gate=g,cluster=c,working_name=names.get(c,'Unannotated'),
                          cluster_cells=counts[s,g,c],gate_cells=total,percent=100*counts[s,g,c]/total))
write('gated_sample_proportions.tsv',props)
comparisons=[('within_HET_frozen_calls',samples[3:], 'WT_like',samples[3:],'KO_like',True),
             ('WT_F_vs_HET_F_A_positive',samples[:3],'A_positive',samples[3:],'A_positive',False),
             ('WT_F_vs_HET_F_any_probe_positive',samples[:3],'any_probe_positive',samples[3:],'any_probe_positive',False),
             ('WT_F_A_positive_vs_HET_KO_like',samples[:3],'A_positive',samples[3:],'KO_like',False)]
tests=[]
for label,sa,ga,sb,gb,paired in comparisons:
    assert all(totals[s,g]>0 for ss,g in [(sa,ga),(sb,gb)] for s in ss)
    for c in clusters:
        a=[100*counts[s,ga,c]/totals[s,ga] for s in sa]
        b=[100*counts[s,gb,c]/totals[s,gb] for s in sb]
        effect=avg(b)-avg(a)
        if paired:
            diffs=[y-x for x,y in zip(a,b)]
            null=[avg([d*z for d,z in zip(diffs,signs)]) for signs in itertools.product([-1,1],repeat=3)]
            leave=[avg(diffs[:i]+diffs[i+1:]) for i in range(3)]
            consistency=sum(d*effect>0 for d in diffs) if effect else 0
        else:
            joined=a+b
            null=[avg([joined[i] for i in ix])-avg([joined[i] for i in range(6) if i not in ix]) for ix in itertools.combinations(range(6),3)]
            leave=[avg(b)-avg(a[:i]+a[i+1:]) for i in range(3)] + [avg(b[:i]+b[i+1:])-avg(a) for i in range(3)]
            leave += [avg(b[:j]+b[j+1:])-avg(a[:i]+a[i+1:]) for i in range(3) for j in range(3)]
            consistency=sum((y-x)*effect>0 for x in a for y in b) if effect else 0
        tests.append(dict(comparison=label,cluster=c,working_name=names.get(c,'Unannotated'),reference_gate=ga,case_gate=gb,
                          reference_mean_pct=avg(a),case_mean_pct=avg(b),difference_pp=effect,
                          exact_p=sum(abs(x)>=abs(effect)-1e-10 for x in null)/len(null),bh_q=1.,
                          leave_out_direction_stable=bool(effect and all(x*effect>0 for x in leave)),
                          directional_agreement=consistency,minimum_reference_cells=min(counts[s,ga,c] for s in sa),
                          minimum_case_cells=min(counts[s,gb,c] for s in sb),paired=paired))
ordered=sorted(tests,key=lambda r:r['exact_p']);q=1.
for i in range(len(ordered)-1,-1,-1):
    q=min(q,ordered[i]['exact_p']*len(ordered)/(i+1));ordered[i]['bh_q']=q
tests.sort(key=lambda r:(r['comparison'],-abs(r['difference_pp'])))
write('all_gated_contrasts.tsv',tests)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.linewidth':.7,'pdf.fonttype':42,'svg.fonttype':'none'})
for label,sa,ga,sb,gb,paired in comparisons:
    top=sorted([r for r in tests if r['comparison']==label],key=lambda r:-abs(r['difference_pp']))[:4]
    fig,axes=plt.subplots(1,4,figsize=(8.2,3.4))
    for i,(ax,t) in enumerate(zip(axes,top)):
        c=t['cluster'];a=np.array([100*counts[s,ga,c]/totals[s,ga] for s in sa]);b=np.array([100*counts[s,gb,c]/totals[s,gb] for s in sb])
        if paired:
            for j in range(3):ax.plot([-.1+j*.1,.9+j*.1],[a[j],b[j]],color='#bbbbbb',lw=.8,zorder=1)
        for x,y,color in [(0,a,'#549E89'),(1,b,'#A678AC')]:
            ax.scatter(x+np.array([-.1,0,.1]),y,s=28,color=color,edgecolors='white',linewidths=.4,zorder=3)
            ax.plot([x-.22,x+.22],[y.mean(),y.mean()],color=color,lw=1.6)
        import textwrap
        title='\n'.join(textwrap.wrap(t['working_name'],19))+'\n'+c
        ax.set_title(title,fontsize=8,pad=24)
        ax.text(.5,1.03,'%+.2f pp'%t['difference_pp'],transform=ax.transAxes,ha='center',fontsize=8)
        ax.text(-.22,1.45,'abcd'[i],transform=ax.transAxes,fontweight='bold',fontsize=12)
        ax.set_xticks([0,1]);ax.set_xticklabels(['WT-like','KO-like'] if paired else ['WT female\n'+ga.replace('_',' '),'HET female\n'+gb.replace('_',' ')],fontsize=7)
        ax.set_xlim(-.5,1.5);ax.set_ylim(bottom=0);ax.set_ylabel('Cells (% of gated subset)')
        for side in ['top','right']:ax.spines[side].set_visible(False)
    fig.subplots_adjust(left=.08,right=.99,bottom=.28,top=.61,wspace=.75)
    note='Paired within HET females; KO-like calls have 76.3% control precision.' if paired else 'Between female samples; probe-positive selection changes the measured population.'
    fig.text(.08,.07,note,fontsize=7)
    fig.text(.08,.025,'Each dot: one sample; line: mean. n = 3 per subset. Exploratory, selected by effect size.',fontsize=7)
    for ext in ['pdf','svg','png']:fig.savefig(OUT/(label+'.'+ext),dpi=400,bbox_inches='tight',facecolor='white')
    plt.close(fig)
provenance={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
(OUT/'provenance.json').write_text(json.dumps(provenance,indent=2))
report=['# Probe-gated cell composition','',
        'Frozen Step 07 thresholds were reused unchanged; no classifier was refitted. Cells were joined by exact technical sample plus barcode to the 446,349-cell consensus. Cells excluded by primary QC are excluded here and counted separately.',
        '', 'WT-like and KO-like are probe-evidence states, not definitive DNA genotypes. WT-like held-out precision was 99.684%; KO-like precision 76.316%, with extremely low sensitivity. No-probe and uncertain HET cells are excluded from the frozen-call comparison, not relabeled. These counts cannot estimate the biological mosaic ratio.',
        '', 'Four analyses: paired WT-like versus KO-like subsets within the same HET samples; WT versus HET females selected identically for A positivity; WT versus HET females selected identically for any probe positivity; and WT-female A-positive versus HET-female KO-like subsets (different gates; descriptive only).',
        '', 'Composition is normalized within each sample and gate. This estimates composition among detected/classifiable cells, not unbiased tissue composition. Probe detectability varies by cluster; differential detectability and genotype-specific gate sensitivity may drive apparent enrichment. No adjustment for detection bias is claimed.',
        '', 'Exact paired sign-flip P values have minimum 0.25 with three HET samples. Unpaired exact P values have minimum 0.10. BH correction covers all 136 contrasts jointly. Independent biological replication and preparation covariates remain unverified.',
        '', '## Retained HET call coverage','']
for r in cov:
    if r['sample'] in samples[3:] and r['gate'] in ['WT_like','KO_like','Uncertain','Uncalled_000']:
        report.append('- {sample}: {gate} = {gate_cells}'.format(**r))
for label,*_ in comparisons:
    report.extend(['','## '+label,'','| Cluster | Working name | Reference % | Case % | Change pp | Leave-out stable |','|---|---|---:|---:|---:|---|'])
    for r in [r for r in tests if r['comparison']==label][:8]:
        report.append('| {cluster} | {working_name} | {reference_mean_pct:.2f} | {case_mean_pct:.2f} | {difference_pp:+.2f} | {leave_out_direction_stable} |'.format(**r))
(OUT/'REPORT.md').write_text('\n'.join(report)+'\n')
print('\n'.join(report[12:]))
