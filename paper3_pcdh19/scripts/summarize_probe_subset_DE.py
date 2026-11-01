"""Validate and summarize sample-level probe-subset differential expression."""
import csv
import hashlib
import json
from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
ROOT=Path(__file__).resolve().parents[1]
OUT=RUNTIME_ROOT/'results/probe_subset_DE_20261031'
def read(name):return list(csv.DictReader((OUT/name).open(),delimiter='\t'))
summary=read('DE_summary.tsv')
metadata=read('sample_metadata.tsv')
assert len(metadata)==9 and all(int(r['cells'])>0 for r in metadata)
report=['# WT-female group versus probe-defined cells in HET females','',
        'PRELIMINARY. Raw Step 06 counts were summed by registered sample and subset. WT females retain all QC-passing cells. HET subsets use the unchanged frozen Step 07 WT-like/KO-like calls. Each comparison uses three WT-female libraries and three HET-subset libraries; cells are not treated as independent replicates.',
        '', 'edgeR 4.4.2: filterByExpr, TMM normalization, robust dispersion estimation, robust quasi-likelihood GLM and a two-sided case-versus-WT test. BH FDR is calculated within each comparison. Top ten genes are ranked by P value after expression filtering. Positive log2FC means higher expression in the HET subset. Pcdh19 is retained in complete tables but omitted from the discovery top ten because its probes define the gate.',
        '', 'These are composition-unadjusted group/subset expression contrasts. Different proportions of progenitors and neurons, probe detectability, and classification error can contribute; the results do not isolate cell-intrinsic mosaic effects. Independent embryo/litter structure and preparation covariates remain unverified. KO-like calls have 76.316% control precision and are inferred expression states, not definitive DNA genotypes.',
        '', '## Libraries','', '| Group | Sample | Cells | Total UMI |','|---|---|---:|---:|']
for r in metadata:report.append('| {group} | {sample} | {cells} | {total_UMI} |'.format(**r))
checks={}
for s in summary:
    case=s['comparison']
    rows=read(case+'_vs_WT_F_top10.tsv')
    full=read(case+'_vs_WT_F_all_genes.tsv')
    assert len(rows)==10 and len({r['gene_id'] for r in rows})==10
    assert all(r['gene_symbol']!='Pcdh19' for r in rows)
    assert all(0<=float(r['PValue'])<=float(r['FDR'])<=1 for r in full)
    checks[case]=dict(genes_tested=len(full),top10_unique=True,pvalue_FDR_bounds=True,Pcdh19_excluded_from_top10=True)
    report += ['', '## '+case+' versus WT females','',
               'Genes tested: {genes_tested}; FDR ≤ 0.05: {FDR_005}; excluding gate gene: {FDR_005_non_gate_gene}.'.format(**s),
               '', '| Gene | log2FC | FDR | All sample pairs agree |','|---|---:|---:|---|']
    for r in rows:
        report.append('| %s | %+.2f | %.3g | %s |'%(r['gene_symbol'] or r['gene_id'],float(r['logFC']),float(r['FDR']),r['all_9_sample_pairs_same_direction']))
(OUT/'REPORT.md').write_text('\n'.join(report)+'\n')
(OUT/'DE_validation.json').write_text(json.dumps(checks,indent=2))
files=[OUT/'raw_pseudobulk_counts.tsv.gz',OUT/'sample_metadata.tsv',ROOT/'scripts/probe_subset_pseudobulk.py',ROOT/'scripts/probe_subset_pseudobulk.R']
(OUT/'analysis_provenance.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},indent=2))
print('\n'.join(report))
