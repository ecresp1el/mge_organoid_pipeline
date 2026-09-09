#!/usr/bin/env python3
"""Package the completed identity extension; leaves all frozen results untouched."""
from pathlib import Path
import hashlib
import json
import shutil
import zipfile
from html.parser import HTMLParser
from urllib.parse import unquote,urlsplit

ROOT=Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results')
OUT=ROOT/'div30_antecedent_identity_v1'
REPO=Path(__file__).resolve().parents[2]

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8<<20),b''):h.update(b)
    return h.hexdigest()

class Links(HTMLParser):
    def __init__(self):super().__init__();self.links=[]
    def handle_starttag(self,tag,attrs):
        for k,v in attrs:
            if k in ['href','src'] and v:self.links.append(v)

def main():
    required=['REPORT.md','index.html','BIOLOGICAL_CONTEXT.md','DIV30_Biological_Identity_Figures.pdf',
              'within/summary.json','within/provenance/validation.json',
              'mapping/frozen_mapping_reconstruction_audit.json','mapping/surface_translation_review.md',
              'external/findings.md','external/summary.json','provenance/independent_identity_validation.json',
              'provenance/preservation_verification.json']
    for name in required:assert (OUT/name).is_file(),name
    preservation=json.loads((OUT/'provenance/preservation_verification.json').read_text())
    assert preservation['passed'] and preservation['files_checked']==2245
    assert json.loads((OUT/'provenance/independent_identity_validation.json').read_text())['passed']
    assert json.loads((OUT/'within/provenance/validation.json').read_text())['status']=='passed'
    assert json.loads((OUT/'external/summary.json').read_text())['source_preserved']
    archive=OUT/'DIV30_Biological_Identity_Review.zip'
    parser=Links();parser.feed((OUT/'index.html').read_text())
    checked=[]
    for raw in parser.links:
        q=urlsplit(raw)
        if q.scheme or not q.path:continue
        p=OUT/unquote(q.path)
        if p==archive:continue
        assert p.is_file(),raw
        checked.append(q.path)
    code=OUT/'provenance/code';code.mkdir(exist_ok=True)
    for p in (REPO/'python_notebooks/scripts').glob('*div30_identity*'):
        if p.is_file() and p.suffix in {'.py','.R','.r'}:shutil.copy2(p,code/p.name)
    shutil.copy2(REPO/'python_notebooks/HANDOFF_div30_antecedent_identity.md',
                 OUT/'HANDOFF_div30_antecedent_identity.md')
    (OUT/'REVIEW_GUIDE.md').write_text('''# Biological identity review

Open `index.html` for the full report and figure gallery. The report answers
the 13 questions and includes all 17 requested analysis areas. PDF and SVG
versions accompany the 600-dpi PNG figures. Full per-cell results and cached
reference spaces remain in the same Turbo output directory and are omitted
from this compact archive.

The original target, scores, gate thresholds, large searches and prior PI
package are unchanged. The archive contains a new biological interpretation,
not a replacement target or another gate-discovery result. Nearest neighbors
are resemblance diagnostics, not temporal tracking or calibrated fate
probabilities. Cross-study frequencies and matched expression depend on
source-cell selection, distance support, assay coverage and state matching.

The four groups are A=target retained, B=target lost, C=non-target retained,
D=non-target removed. A/B/C/D sizes are 9,210/8,917/6,571/65,933. Frozen gate
composition is 58.36% and recovery 50.81%; held-sample gate reselection remains
53.49% composition and 39.73% recovery. New held-sample prediction models
explain biological associations and are not a repeat of gate validation.
''')
    files=[]
    for p in sorted(OUT.rglob('*')):
        if not p.is_file():continue
        rel=p.relative_to(OUT)
        if any(x in {'cache','logs'} for x in rel.parts):continue
        if p.suffix in {'.npy','.npz','.h5','.h5ad','.zip','.gz','.f32','.i32'}:continue
        if p.name in {'REVIEW_FILE_MANIFEST.json','REVIEW_PACKAGE_VALIDATION.json'}:continue
        if p.stat().st_size>30<<20:continue
        files.append(p)
    manifest={str(p.relative_to(OUT)):{'sha256':sha(p),'bytes':p.stat().st_size} for p in files}
    mf=OUT/'REVIEW_FILE_MANIFEST.json';mf.write_text(json.dumps(manifest,indent=2)+'\n')
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in files+[mf]:z.write(p,str(p.relative_to(OUT)))
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        archived=set(z.namelist())
        for name in checked:
            assert name in archived, f'Linked review file omitted from archive: {name}'
        entries=len(z.namelist())
    result=dict(required_outputs_present=len(required),local_html_links_checked=len(checked),
        original_files_preserved=preservation['files_checked'],archive_entries=entries,
        archive_bytes=archive.stat().st_size,archive_sha256=sha(archive),archive_crc_passed=True,
        all_local_html_links_present_in_archive=True,
        large_gate_search_rerun=False,target_changed=False,
        excluded='Large caches and per-cell gzip tables remain in the Turbo result directory.')
    (OUT/'REVIEW_PACKAGE_VALIDATION.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
