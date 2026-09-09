#!/usr/bin/env python3
"""Verify the scientific deliverable inventory and make a portable review bundle."""
import json
from pathlib import Path
import shutil
import zipfile
import xml.etree.ElementTree as ET
from PIL import Image
from run_div90_hypergate_phase2 import OUT, REPO_ROOT, initialize, sha, write_json
from div90_hypergate_phase2_report import FIGURES, INTERACTIVE


def run():
    initialize()
    final=json.loads((OUT/'provenance/final_validation.json').read_text())
    assert final['all_phase1_files_unchanged'] and final['key_gate_counts_and_cell_flags_independently_reproduced']
    browser=json.loads((OUT/'provenance/interactive_validation.json').read_text())
    assert len(browser['assets'])==4 and browser['offline']
    inventory=[]
    for stem,_ in FIGURES:
        for ext in ['pdf','svg','png']:
            path=OUT/'figures'/f'{stem}.{ext}'
            assert path.is_file() and path.stat().st_size>1000
            if ext=='pdf':assert path.read_bytes()[:5]==b'%PDF-'
            elif ext=='svg':ET.parse(path)
            else:
                with Image.open(path) as im:
                    assert min(im.info.get('dpi',(0,0)))>=599
                    im.verify()
            inventory.append(str(path.relative_to(OUT)))
    for name,_ in INTERACTIVE:
        path=OUT/name;assert path.stat().st_size>1000
        inventory.append(name)
    for name in ['REPORT.md','CONDITION_REPORT.md','index.html','cells.tsv.gz','gate_summary.json','gate_cell_assignments.tsv.gz']:
        assert (OUT/name).is_file()
    evidence=OUT/'provenance/DIV90_HYPERGATE_PHASE2_EVIDENCE.md'
    evidence.write_text(evidence.read_text().replace('](../metadata/',']('))
    handoff=(REPO_ROOT/'python_notebooks/HANDOFF_div90_hypergate_phase2.md').read_text().replace('](scripts/','](').replace('](../metadata/','](')
    (OUT/'provenance/HANDOFF_div90_hypergate_phase2.md').write_text(handoff)
    # Preserve the exact original request alongside the structured condition update.
    request=Path('/home/elcrespo/.codex/attachments/50c3ea93-4d27-43b8-baf1-3367076fa461/pasted-text.txt')
    if request.exists():shutil.copy2(request,OUT/'provenance/original_phase2_request.txt')
    write_json(OUT/'provenance/deliverable_validation.json',dict(core_figures=10,supplementary_overlay_figures=1,additional_doublet_supplement=True,formats=['PDF','SVG','600dpi PNG'],all_figure_files_parsed=True,interactive_assets=4,browser_validation='provenance/interactive_validation.json',scientific_validation='provenance/final_validation.json',runtime_root=str(OUT),inventory=inventory))
    excluded={'depletion_all_candidate_metrics.h5','depletion_supplemental_candidate_metrics.h5','depletion_single_all.tsv.gz','gate_ranked_candidates.tsv.gz','failed_random_io_metrics.h5','gate_summary_before_constrained_refinement.json'}
    files=[]
    for path in OUT.rglob('*'):
        if not path.is_file():continue
        relative=path.relative_to(OUT)
        if relative.parts[0] in ['cache','runtime_deps','logs'] or path.suffix=='.zip' or path.name in excluded:continue
        if relative.parts[0] not in ['figures','tables','provenance'] and len(relative.parts)>1:continue
        if path.name.startswith('.') or path.name=='review_bundle_manifest.json':continue
        files.append(path)
    manifest=dict(purpose='Portable review bundle: reports, figures, offline interactives, cells, code/configuration, sensitivity and condition tables. Complete large candidate-search files remain at the Turbo runtime root and are not duplicated here.',omitted_large_search_files=sorted(excluded),files=[dict(path=str(p.relative_to(OUT)),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(files)])
    write_json(OUT/'review_bundle_manifest.json',manifest);files.append(OUT/'review_bundle_manifest.json')
    archive=OUT/'div90_hypergate_phase2_review_bundle.zip'
    if archive.exists():raise FileExistsError(f'Preserve existing bundle before creating a new version: {archive}')
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=4) as z:
        for p in sorted(files):z.write(p,p.relative_to(OUT))
    with zipfile.ZipFile(archive) as z:assert z.testzip() is None
    write_json(OUT/'provenance/review_bundle_validation.json',dict(path=str(archive),bytes=archive.stat().st_size,sha256=sha(archive),n_files=len(files),zip_crc_verified=True))
    print(f'Validated deliverables and created {archive} ({archive.stat().st_size/1e6:.1f} MB)',flush=True)


if __name__=='__main__':run()
