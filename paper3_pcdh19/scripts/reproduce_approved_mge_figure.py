#!/usr/bin/env python3
"""Re-render the approved 2026-11-01 figure from its sealed NFS input bundle."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

DEFAULT_BUNDLE = Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/final_figures/approved_mge_mosaic_20261101')

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', type=Path, default=DEFAULT_BUNDLE)
    parser.add_argument('--output-root', type=Path)
    args = parser.parse_args()
    bundle = args.bundle.resolve()
    output = args.output_root or bundle.parents[1] / 'results' / ('approved_mge_reproduction_' + datetime.datetime.utcnow().strftime('%Y%m%d_%H%M%S'))
    output = output.resolve()
    checkout = Path(__file__).resolve().parents[2]
    if output == checkout or checkout in output.parents or output == bundle or bundle in output.parents:
        raise ValueError('Reproduction output must be outside checkout and sealed bundle')
    for entry in json.loads((bundle / 'frozen_manifest.json').read_text()):
        path = bundle / entry['path']
        if sha(path) != entry['sha256']:
            raise ValueError('Frozen file changed: ' + str(path))
    output.mkdir(parents=True, exist_ok=False)
    for entry in json.loads((bundle / 'input_manifest.json').read_text()):
        source = bundle / entry['path']
        target = output / source.relative_to(bundle / 'inputs')
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(str(source), str(target))
    environment = os.environ.copy()
    environment['PCDH19_RUNTIME_ROOT'] = str(output)
    renderer = bundle / 'code/paper3_pcdh19/scripts/render_HET_palette_expression.py'
    subprocess.run([sys.executable, str(renderer)], env=environment, check=True)
    rendered = output / 'results/female_male_signaling_20261031'
    checks = {}
    for name in ['HET_palette_expression_figure.png', 'HET_palette_expression_gene_values.tsv']:
        checks[name] = {'expected_sha256': sha(bundle / 'approved' / name), 'actual_sha256': sha(rendered / name)}
        checks[name]['exact_match'] = checks[name]['expected_sha256'] == checks[name]['actual_sha256']
    import numpy
    import scipy
    import matplotlib
    from matplotlib import font_manager
    report = {'bundle': str(bundle), 'output_root': str(output), 'checks': checks,
              'python': sys.version, 'numpy': numpy.__version__, 'scipy': scipy.__version__,
              'matplotlib': matplotlib.__version__, 'font': font_manager.findfont('DejaVu Sans'),
              'note': 'PNG and gene table checked byte for byte; PDF/SVG metadata and generated IDs can differ.'}
    (output / 'REPRODUCTION_VALIDATION.json').write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2), flush=True)
    if not all(check['exact_match'] for check in checks.values()):
        raise SystemExit('Reproduction differs from approved reference; inspect report.')

if __name__ == '__main__':
    main()
