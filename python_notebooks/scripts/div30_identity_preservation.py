#!/usr/bin/env python3
"""Read-only full-content preservation audit for the biological identity extension."""
import argparse
import hashlib
import json
from pathlib import Path
from datetime import datetime, timezone

ROOT = Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results')
OUT = ROOT / 'div30_antecedent_identity_v1' / 'provenance'
FROZEN = ['div90_hypergate_sst_pv', 'div90_hypergate_sst_pv_phase2',
          'div90_purification_audit_v1', 'div30_pv_antecedent_v1',
          'div30_div90_purification_comparison_v1']

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as fh:
        for block in iter(lambda: fh.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', choices=['snapshot', 'verify'])
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    dest = OUT / 'frozen_outputs_sha256.json'
    if args.mode == 'snapshot':
        if dest.exists():
            raise SystemExit('Refusing to replace the original preservation snapshot.')
        files = []
        for folder in FROZEN:
            for p in sorted((ROOT / folder).rglob('*')):
                if p.is_file():
                    st = p.stat()
                    files.append(dict(path=str(p.relative_to(ROOT)), size=st.st_size,
                                      mtime_ns=st.st_mtime_ns, sha256=digest(p)))
            print(folder, len(files), flush=True)
        dest.write_text(json.dumps(dict(created_utc=datetime.now(timezone.utc).isoformat(),
                                       files=files), indent=2) + '\n')
        print('Snapshotted', len(files), 'files', sum(x['size'] for x in files), 'bytes', flush=True)
    else:
        snapshot = json.loads(dest.read_text())
        changed = []
        for row in snapshot['files']:
            p = ROOT / row['path']
            if not p.exists() or p.stat().st_size != row['size'] or digest(p) != row['sha256']:
                changed.append(row['path'])
        before={r['path'] for r in snapshot['files']}
        after={str(p.relative_to(ROOT)) for folder in FROZEN for p in (ROOT/folder).rglob('*') if p.is_file()}
        added=sorted(after-before)
        result = dict(verified_utc=datetime.now(timezone.utc).isoformat(),
                      files_checked=len(snapshot['files']), changed=changed, added=added, passed=not changed and not added)
        (OUT / 'preservation_verification.json').write_text(json.dumps(result, indent=2) + '\n')
        print(json.dumps(result), flush=True)
        if changed or added:
            raise SystemExit(1)

if __name__ == '__main__':
    main()
