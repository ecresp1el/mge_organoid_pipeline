#!/usr/bin/env python3
"""Create a compact, offline review bundle after scientific validation succeeds."""
from __future__ import annotations

import hashlib
import html
import json
from pathlib import Path
import shutil
import zipfile

ROOT = Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results')
OUT = ROOT / 'div30_div90_purification_comparison_v1'
D30 = ROOT / 'div30_pv_antecedent_v1'
D90 = ROOT / 'div90_purification_audit_v1'
REPO = Path(__file__).resolve().parents[2]


def digest(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda: f.read(8 << 20), b''): h.update(block)
    return h.hexdigest()


def main():
    import pandas as pd
    # Required checks represent different implementations and source preservation.
    required = [D90 / 'REPORT.md', D30 / 'REPORT.md', OUT / 'REPORT.md',
                OUT / 'provenance/source_preservation.json',
                OUT / 'provenance/independent_legacy_validation.json',
                OUT / 'provenance/independent_comparison_validation.json',
                D90 / 'provenance/report_independent_validation.json',
                D30 / 'provenance/surface_independent_validation.json',
                D30 / 'tables/end_to_end_score_sensitivity.tsv',
                D30 / 'tables/held_sample_surface_gate_summary.tsv',
                OUT / 'figures/DIV30_DIV90_developmental_schematic.pdf']
    for p in required: assert p.is_file(), p
    preservation = json.loads((OUT / 'provenance/source_preservation.json').read_text())
    assert preservation['all_existing_results_preserved'] and not preservation['added_within_frozen_directories']
    critical = pd.read_csv(D90 / 'tables/critical_best_gates.tsv', sep='\t')
    assert set(critical.criterion) == set('ABCDEFGHI')
    summary90 = json.loads((D90 / 'audit_summary.json').read_text())
    assert summary90['legacy_evaluated_rows'] == 2885531
    source_scripts = [p for p in (REPO / 'python_notebooks/scripts').glob('*.py') if any(k in p.name for k in ['div30_antecedent', 'div90_purification', 'div30_div90_purification', 'div30_div90_comparison'])]
    (OUT / 'provenance/code').mkdir(exist_ok=True)
    for p in source_scripts: shutil.copy2(p, OUT / 'provenance/code' / p.name)
    b = critical.set_index('criterion').loc['I']
    d = pd.read_csv(D30 / 'tables/surface_gate_best_at_recovery_floors.tsv', sep='\t')
    b30 = d[(d.target == 'top20') & (d.recovery_floor == .5) & (d.maximum_markers == 3)].iloc[0]
    page = f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>DIV90 purification and separate DIV30 antecedent analysis</title>
<style>body{{font:17px/1.55 system-ui,sans-serif;margin:40px auto;max-width:1080px;padding:0 24px;color:#243b49;background:#fafcfd}}h1{{font-size:32px;line-height:1.2}}h2{{font-size:23px;margin-top:32px}}a{{color:#146e98}}.cards{{display:grid;grid-template-columns:1fr 1fr;gap:22px}}article{{padding:22px;background:white;border:1px solid #dbe4e9;border-radius:10px}}img{{max-width:100%;height:auto}}.small{{font-size:14px;color:#536a78}}@media(max-width:700px){{.cards{{grid-template-columns:1fr}}}}</style>
<h1>DIV90 purification and a separate DIV30 antecedent search</h1>
<p>Completed from the frozen PI package. Original results are preserved. Each age uses its own target, feature screening and gate optimization.</p>
<div class="cards"><article><h2>DIV90: stronger purity costs recovery</h2><p>Best three-marker rule at ≥10% PV recovery: <b>{b.pv_fraction:.1%} strict PV-biased purity</b>, {b.pv_recovery:.1%} PV recovery and {b.sst_fraction:.1%} SST contamination. At ≥50% PV recovery, the best purity is {critical.set_index('criterion').loc['B','pv_fraction']:.1%}.</p><p><a href="../div90_purification_audit_v1/REPORT.md">Full DIV90 report</a> · <a href="../div90_purification_audit_v1/tables/critical_best_gates.tsv">Exact A–I gate table</a></p></article>
<article><h2>DIV30: enrich a candidate transcriptional state</h2><p>Its independent top20 consensus target rises from <b>{b30.starting_target_fraction:.1%} to {b30.post_gate_target_fraction:.1%}</b>, with {b30.target_recovery:.1%} recovery. This is developmental resemblance, not a confirmed future PV identity.</p><p><a href="../div30_pv_antecedent_v1/REPORT.md">Full DIV30 report</a> · <a href="../div30_pv_antecedent_v1/tables/surface_gate_best_at_recovery_floors.tsv">Gate tradeoffs</a></p></article></div>
<p>FGFR2/PTPRS remains a weak DIV90 benchmark: SST-biased composition changes 22.5% → 21.3%. Independent statistical significance is not established. No antibodies are recommended.</p>
<h2>Linked interpretation</h2><p><a href="REPORT.md">Comparison report</a> · <a href="tables/important_marker_temporal_classification.tsv">Temporal marker classifications</a> · <a href="tables/cross_age_gate_benchmarks.tsv">Exact cross-age diagnostic applications</a></p>
<a href="figures/DIV30_DIV90_developmental_schematic.pdf"><img src="figures/DIV30_DIV90_developmental_schematic.png" alt="Separate DIV30 and DIV90 analyses linked by dashed developmental hypotheses to user-reported later PVALB protein"></a>
<h2>DIV90 purity/recovery frontier</h2><a href="../div90_purification_audit_v1/figures/01_pv_purification_frontier.pdf"><img src="../div90_purification_audit_v1/figures/01_pv_purification_frontier.png" alt="Full observed PV purity recovery Pareto frontier and marker complexity comparison"></a>
<p class="small">All {summary90['legacy_evaluated_rows']:,} legacy rows were reconstructed; expanded DIV90 evaluation covers {summary90['total_evaluated_rows']:,} boundaries, including repeats. Gridded pairs and bounded triplets do not prove a global theoretical optimum. Complete large candidate archives remain in each analysis directory on Turbo; the review ZIP contains compact reports, tables, figures and validation evidence.</p></html>'''
    (OUT / 'index.html').write_text(page)
    (OUT / 'REVIEW_GUIDE.md').write_text('''# Review guide

Open `index.html`, then the linked comparison report and each independent age report.
The original DIV90 PI package remains in `div90_hypergate_sst_pv_phase2/pi_figure_package_v1/` and is historical context, not the updated gate recommendation.

The compact ZIP preserves three sibling directory names, so relative links work after extraction. It includes exact shortlisted rules, recovery/frontier tables, sensitivity and held-sample summaries, figures, source preservation checks and the code snapshot. Large HDF5 candidate archives, expression caches, and per-cell tables remain in the original new Turbo output directories and are explicitly excluded from the review ZIP to keep it manageable.

RNA cutoffs are saved at full precision. They are not fluorescence thresholds. DIV30 and DIV90 percentages refer to different targets and entry populations; no pooled optimization or lineage tracing was performed.
''')
    # Exclude full candidate and cell archives; keep actionable summaries and evidence.
    skip_names = {'surface_gate_candidate_registry.tsv.gz', 'surface_triplet_materiality_audit.tsv.gz',
                  'reference_oof_predictions.tsv.gz', 'reference_cell_index.tsv',
                  'cell_sensitivity.tsv.gz', 'end_to_end_cell_sensitivity.tsv.gz', 'end_to_end_cell_score_sensitivity.tsv.gz',
                  'target_cell_assignments.tsv.gz', 'gate_cell_assignments.tsv.gz'}
    files = []
    for base in [D90, D30, OUT]:
        for p in sorted(base.rglob('*')):
            if not p.is_file(): continue
            rel = p.relative_to(base)
            if rel.parts[0] in {'cache', 'logs', 'runtime_deps'}: continue
            if p.suffix in {'.h5', '.h5ad', '.npy', '.npz', '.zip'}: continue
            if p.name in skip_names or p.name in {'cells.tsv.gz', 'REVIEW_FILE_MANIFEST.json', 'REVIEW_PACKAGE_VALIDATION.json'}: continue
            if 'failed' in p.name or 'partial' in p.name: continue
            if p.stat().st_size > 25 << 20: continue
            files.append(p)
    manifest = {str(p.relative_to(ROOT)): {'sha256': digest(p), 'bytes': p.stat().st_size} for p in files}
    mf = OUT / 'REVIEW_FILE_MANIFEST.json'; mf.write_text(json.dumps(manifest, indent=2) + '\n')
    archive = OUT / 'DIV30_DIV90_Purification_Review.zip'
    with zipfile.ZipFile(archive, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for p in files + [mf]: z.write(p, str(p.relative_to(ROOT)))
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        entries = len(z.namelist())
    validation = {'required_files_present': len(required), 'critical_criteria': 'ABCDEFGHI',
                  'legacy_rows_reconstructed': summary90['legacy_evaluated_rows'],
                  'original_results_preserved': True, 'archive_entries': entries,
                  'archive_bytes': archive.stat().st_size, 'archive_sha256': digest(archive),
                  'archive_CRC_passed': True, 'large_candidate_archives_remain_on_Turbo': True}
    (OUT / 'REVIEW_PACKAGE_VALIDATION.json').write_text(json.dumps(validation, indent=2) + '\n')
    print(json.dumps(validation, indent=2), flush=True)


if __name__ == '__main__': main()
