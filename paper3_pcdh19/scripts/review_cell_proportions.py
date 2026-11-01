"""Sample-level exploratory composition contrasts; no cell-level replication."""
import csv
import gzip
import hashlib
import itertools
import json
import math
import html
from collections import Counter
from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT

ROOT = Path(__file__).resolve().parents[1]
OUT = RUNTIME_ROOT / 'results/composition_review_20261031'
BASE = Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat')
PILOT = BASE / '02_hierarchy_validation/hicat_validation_20260908_221625_d5d06a5d/review/merge_collapse_20260908_225259/outputs/cell_memberships.tsv.gz'
FULL = BASE / '03_full_data_consensus_benchmark/hicat_consensus_review_20260909_142348_efd9e13d/outputs/tables/all_cell_assignments.tsv.gz'


def read(path):
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as f:
        return list(csv.DictReader(f, delimiter=',' if path.suffix == '.csv' else '\t'))


def write(name, rows):
    with (OUT / name).open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def mean(x):
    return sum(x) / len(x)


def analyze():
    OUT.mkdir(exist_ok=True, parents=True)
    key = read(ROOT / 'config/sample_key.csv')
    samples = [r['technical_sample_id'] for r in key]
    group = {r['technical_sample_id']: r['design_group'] for r in key}
    annotations = {r['cluster']: r for r in read(ROOT / 'config/hicat_pilot_best_bet_annotations.tsv')}
    counts, totals = Counter(), Counter()
    pilot = read(PILOT)
    assert len(pilot) == 12000 and len({r['cell_id'] for r in pilot}) == 12000
    for r in pilot:
        s, c = r['sample'], r['original']
        parent = c.split('.')[0]
        for level, label, denom in [('pilot_parent', parent, 'all'), ('pilot_subcluster', c, 'all'),
                                     ('pilot_cell_type', annotations[c]['cell_type'], 'all'),
                                     ('pilot_within_parent', c, parent)]:
            counts[level, label, denom, s] += 1
        totals['pilot', 'all', s] += 1
        totals['pilot', parent, s] += 1
    assert all(totals['pilot', 'all', s] == 1000 for s in samples)
    full = read(FULL)
    assert len(full) == 446349 and len({r['cell_id'] for r in full}) == len(full)
    for r in full:
        s = r['cell_id'].split('_')[0]
        assert s in group
        counts['full_consensus', r['diagnostic_cluster'], 'all', s] += 1
        totals['full', 'all', s] += 1
    categories = sorted({k[:3] for k in counts})
    proportions, tests = [], []
    groups = ['WT_M', 'WT_F', 'HET_F', 'KO_M']
    for level, label, denom in categories:
        dataset = 'full' if level == 'full_consensus' else 'pilot'
        values = {}
        for s in samples:
            n, total = counts[level, label, denom, s], totals[dataset, denom, s]
            values[s] = n / total if total else None
            proportions.append(dict(level=level, category=label, denominator=denom, sample=s,
                                    group=group[s], count=n, total=total, fraction=values[s]))
        for ga, gb in itertools.combinations(groups, 2):
            sa = [s for s in samples if group[s] == ga]
            sb = [s for s in samples if group[s] == gb]
            if any(values[s] is None for s in sa + sb):
                continue
            a, b = [values[s] for s in sa], [values[s] for s in sb]
            effect = mean(b) - mean(a)
            joined = a + b
            null = [mean([joined[i] for i in ix]) - mean([joined[i] for i in range(6) if i not in ix])
                    for ix in itertools.combinations(range(6), 3)]
            p = sum(abs(v) >= abs(effect) - 1e-12 for v in null) / len(null)
            # Drop each sample independently, then drop one from each group together.
            loo = [mean(b[:i] + b[i+1:]) - mean(a) for i in range(3)]
            loo += [mean(b) - mean(a[:i] + a[i+1:]) for i in range(3)]
            loo += [mean(b[:j] + b[j+1:]) - mean(a[:i] + a[i+1:]) for i in range(3) for j in range(3)]
            stable = effect != 0 and all(v * effect > 0 for v in loo)
            pairwise = sum((y-x)*effect > 0 for x in a for y in b) if effect else 0
            complete = max(a) < min(b) or max(b) < min(a)
            minimum = min(counts[level, label, denom, s] for s in sa + sb)
            tests.append(dict(level=level, category=label, denominator=denom, reference=ga, case=gb,
                              reference_mean_pct=100*mean(a), case_mean_pct=100*mean(b),
                              difference_pp=100*effect, exact_p=p, bh_q=1.,
                              all_leave_out_directions_agree=stable, concordant_sample_pairs=pairwise,
                              complete_sample_separation=complete, minimum_cells_per_sample=minimum,
                              minimum_abs_leave_out_effect_pp=100*min(abs(v) for v in loo),
                              sex_matched=(ga,gb) in [('WT_M','KO_M'),('WT_F','HET_F')]))
    # Correct jointly across every tested category and pair; overlaps are retained explicitly.
    ordered = sorted(tests, key=lambda r:r['exact_p'])
    q = 1.
    for i in range(len(ordered)-1, -1, -1):
        q = min(q, ordered[i]['exact_p']*len(ordered)/(i+1))
        ordered[i]['bh_q'] = q
    tests.sort(key=lambda r:(not r['sex_matched'], not r['complete_sample_separation'],
                             not r['all_leave_out_directions_agree'], -abs(r['difference_pp'])))
    write('sample_proportions.tsv', proportions)
    write('all_contrasts.tsv', tests)
    candidates = [r for r in tests if r['complete_sample_separation'] and r['all_leave_out_directions_agree']]
    if candidates:
        write('directionally_consistent_candidates.tsv', candidates)
    provenance = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in [PILOT,FULL,ROOT/'config/sample_key.csv', ROOT/'config/hicat_pilot_best_bet_annotations.tsv']}
    (OUT/'provenance.json').write_text(json.dumps(dict(inputs=provenance,pilot_cells=len(pilot),full_cells=len(full),tests=len(tests)),indent=2))
    def table(rows):
        lines=['| Level | Category | Case − reference | Difference (pp) | Exact P | Min cells/sample |',
               '|---|---|---|---:|---:|---:|']
        for r in rows:
            lines.append('| {level} | {category} | {case} − {reference} | {difference_pp:+.2f} | {exact_p:.2f} | {minimum_cells_per_sample} |'.format(**r))
        return '\n'.join(lines)
    matched=[r for r in candidates if r['sex_matched']]
    report = '''# Exploratory cell composition review

Twelve registered samples, three per group. Independent embryo/litter structure and preparation covariates remain unverified. Proportions use each sample's own denominator, including zero-count categories. Named types derive from working pilot marker interpretations; no pilot annotation is transferred onto full-data consensus IDs.

All six group pairs tested for pilot parents, pilot fine clusters, grouped cell types, fine clusters within their own parent, and 34 full-data consensus clusters. Parent/cell-type/fine tests overlap; BH correction covers all tests jointly. Arbitrary unions of categories were not searched. WT_F versus HET_F and WT_M versus KO_M are sex-matched; other contrasts include sex differences.

Exact two-sided permutation P values enumerate all 20 allocations of six samples. The smallest attainable P is 0.10: no contrast can reach 0.05 even before correction. Candidate means descriptive consistency, not statistical confirmation or independent replication. Candidate selection requires strict separation of all three sample proportions in each group and unchanged direction after removing any one sample or one sample from each group. Small counts are shown explicitly.

## Sex-matched directionally consistent candidates

'''+table(matched)+ '\n\n## All other directionally consistent candidates\n\n'+table([r for r in candidates if not r['sex_matched']])+ '\n'
    (OUT/'REPORT.md').write_text(report)
    # Standalone, dependency-free sample dotplots. Each point is one sample.
    chosen = []
    for level, categories_to_plot in [('full_consensus',['C0001','C0002','C0003','C0005','C0030']),
                                      ('pilot_cell_type',['RG','IPC','IPC / cycling neuroblast','Immature MGE interneuron'])]:
        for category in categories_to_plot:
            chosen.append((level, category))
    svg=['<svg xmlns="http://www.w3.org/2000/svg" width="1100" height="%d" viewBox="0 0 1100 %d">' % (100+len(chosen)*180,100+len(chosen)*180),
         '<rect width="100%" height="100%" fill="white"/>',
         '<style>text{font-family:sans-serif;font-size:13px}</style>',
         '<text x="25" y="28" font-size="19">Cell proportions: each dot is one registered sample</text>',
         '<text x="25" y="52">Consensus IDs and pilot IDs are different taxonomies. All denominators here are all cells in each sample.</text>']
    colors=['#3377aa','#55aa88','#cc8866','#aa5588']
    for index,(level,category) in enumerate(chosen):
        y=90+index*180
        selected=[r for r in proportions if r['level']==level and r['category']==category]
        ceiling=max(1.,max(100*r['fraction'] for r in selected)*1.15)
        svg.append('<text x="25" y="%d">%s: %s</text>'%(y,level,html.escape(category)))
        for gi,g in enumerate(groups):
            x=390+gi*165
            vals=[r for r in selected if r['group']==g]
            for j,r in enumerate(vals):
                py=y+115-90*100*r['fraction']/ceiling
                tip='%s: %.3f%% (%d/%d)'%(r['sample'],100*r['fraction'],r['count'],r['total'])
                svg.append('<circle cx="%d" cy="%.2f" r="6" fill="%s"><title>%s</title></circle>'%(x+(j-1)*14,py,colors[gi],tip))
            my=y+115-90*mean([100*r['fraction'] for r in vals])/ceiling
            svg.append('<path d="M%d %.2f h52" stroke="black" stroke-width="2"/>'%(x-26,my))
            svg.append('<text x="%d" y="%d" text-anchor="middle">%s</text>'%(x,y+140,g))
        for tick in [0,ceiling/2,ceiling]:
            ty=y+115-90*tick/ceiling
            svg.append('<text x="330" y="%.2f" text-anchor="end">%.2f%%</text>'%(ty+4,tick))
    svg.append('</svg>')
    (OUT/'sample_proportions.svg').write_text('\n'.join(svg))
    print(json.dumps(dict(tests=len(tests),candidates=len(candidates),sex_matched_candidates=len(matched),minimum_q=min(r['bh_q'] for r in tests))))
    print(table(matched[:25]))


if __name__ == '__main__':
    analyze()
