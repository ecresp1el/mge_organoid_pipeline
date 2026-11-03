#!/usr/bin/env python3
"""Read-only, bounded storage audit of this PCDH19 runtime, without symlink traversal."""
import argparse
import csv
import datetime
import json
import os
from pathlib import Path

DEFAULT_ROOT = Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19')


def audit(root, output, max_depth=12):
    root = root.resolve()
    if root != DEFAULT_ROOT.resolve():
        raise ValueError('This audit is restricted to the known PCDH19 study runtime')
    root_stat = root.stat()
    totals = {}
    largest = []
    skipped = []
    errors = []
    files = 0
    def walk(folder, depth):
        nonlocal files
        if depth > max_depth:
            skipped.append({'path': str(folder), 'reason': 'maximum depth'})
            return
        try:
            with os.scandir(str(folder)) as entries:
                for entry in entries:
                    path = Path(entry.path)
                    try:
                        if entry.is_symlink():
                            skipped.append({'path': str(path), 'reason': 'symlink; not followed'})
                            continue
                        stat = entry.stat(follow_symlinks=False)
                        if stat.st_dev != root_stat.st_dev:
                            skipped.append({'path': str(path), 'reason': 'different filesystem; not traversed'})
                        elif entry.is_dir(follow_symlinks=False):
                            walk(path, depth + 1)
                        elif entry.is_file(follow_symlinks=False):
                            files += 1
                            relative = path.relative_to(root)
                            key = '/'.join(relative.parts[:2]) if len(relative.parts) > 2 else str(relative.parts[0])
                            row = totals.setdefault(key, {'path': key, 'files': 0, 'bytes': 0, 'allocated_bytes': 0})
                            row['files'] += 1; row['bytes'] += stat.st_size
                            row['allocated_bytes'] += stat.st_blocks * 512
                            largest.append({'path': str(relative), 'bytes': stat.st_size})
                            if len(largest) > 200:
                                largest.sort(key=lambda r: -r['bytes']); del largest[100:]
                    except OSError as exc:
                        errors.append({'path': str(path), 'error': str(exc)})
        except OSError as exc:
            errors.append({'path': str(folder), 'error': str(exc)})
    walk(root, 0)
    output.mkdir(parents=True, exist_ok=True)
    rows = sorted(totals.values(), key=lambda r: -r['allocated_bytes'])
    with (output/'storage_by_directory.tsv').open('w') as handle:
        writer=csv.DictWriter(handle, fieldnames=['path','files','bytes','allocated_bytes'], delimiter='\t')
        writer.writeheader(); writer.writerows(rows)
    report={'root':str(root),'date_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'files':files,'bytes':sum(r['bytes'] for r in rows),
            'allocated_bytes':sum(r['allocated_bytes'] for r in rows),
            'max_depth':max_depth,'symlinks_followed':False,'skipped':skipped,'errors':errors,
            'largest_files':sorted(largest,key=lambda r:-r['bytes'])[:40],
            'note':'Read-only file-size inventory, not project quota. Hard-linked copies may be counted more than once.'}
    (output/'storage_audit.json').write_text(json.dumps(report,indent=2))
    print(json.dumps({'files':files,'GiB':round(report['bytes']/2**30,2),
                      'allocated_GiB':round(report['allocated_bytes']/2**30,2),
                      'skipped':len(skipped),'errors':len(errors),'top_directories':rows[:15]},indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();audit(DEFAULT_ROOT,args.output)
