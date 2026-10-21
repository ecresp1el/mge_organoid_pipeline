#!/usr/bin/env python3
"""Freeze GO/SynGO annotations and build auditable neuronal functional sets."""
import csv
import functools
import gzip
import hashlib
import io
import json
import re
from pathlib import Path
import subprocess
import sys
import zipfile
import xml.etree.ElementTree as ET

SEEDS = {
    'ion_channel': ['GO:0005216'],
    'excitability': ['GO:0042391', 'GO:0001508'],
    'synaptic': ['GO:0007268', 'GO:0050808', 'GO:0045202'],
    'maturation': ['GO:0030182', 'GO:0007411', 'GO:0043194'],
}
PANELS = {
    'canonical_immediate_early': 'FOS FOSB JUN JUNB EGR1 EGR2 EGR3 EGR4 ARC NPAS4 NR4A1 NR4A2 NR4A3 IER2 IER3 IER5 DUSP1 DUSP5 DUSP6'.split(),
    'activity_associated_regulation': 'CREB1 CREM SRF MEF2A MEF2C MEF2D BDNF NPTX2 HOMER1 RGS2'.split(),
    'gabaergic_competence': 'GAD1 GAD2 SLC6A1 SLC32A1 SLC6A11 ABAT GPHN ARHGEF9'.split(),
    'firing_and_axon_specialization': 'ANK3 SPTBN4 NFASC SCN1A SCN2A SCN3A SCN8A SCN1B SCN2B KCNC1 KCNC2 KCNC3 KCNIP1 KCNQ2 KCNQ3 ATP1A3 HCN1 HCN2 HCN4'.split(),
}

def xlsx_rows(data):
    z = zipfile.ZipFile(io.BytesIO(data))
    ns = {'s': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    strings = []
    if 'xl/sharedStrings.xml' in z.namelist():
        strings = [''.join(x.itertext()) for x in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si', ns)]
    rows = []
    for row in ET.fromstring(z.read('xl/worksheets/sheet1.xml')).findall('s:sheetData/s:row', ns):
        out = {}
        for c in row:
            col = ''.join(x for x in c.attrib['r'] if x.isalpha())
            v = c.find('s:v', ns)
            value = v.text if v is not None else ''.join(c.itertext())
            out[col] = strings[int(value)] if c.get('t') == 's' else value
        rows.append(out)
    header = rows.pop(0)
    return [{name: row.get(col, '') for col, name in header.items()} for row in rows]

def prepare(root):
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    urls = {
        'goa_human.gaf.gz': 'https://current.geneontology.org/annotations/goa_human.gaf.gz',
        'syngo1.3_complete_data.zip': 'https://www.syngoportal.org/data/syngo1.3_complete_data.zip',
    }
    manifest = {}
    for name, url in urls.items():
        print('REFERENCE_DOWNLOAD', name, flush=True)
        subprocess.run(['curl', '-fLsS', '--retry', '2', '--max-time', '180', url, '-o', str(root/name)], check=True)
        manifest[name] = {'url': url, 'sha256': hashlib.sha256((root/name).read_bytes()).hexdigest()}
    with gzip.open(root/'goa_human.gaf.gz','rt') as f:
        header = []
        for line in f:
            if not line.startswith('!'): break
            header.append(line)
    version = re.search(r'go-version: .*?/releases/(\d{4}-\d{2}-\d{2})/', ''.join(header))
    if version is None: raise ValueError('Missing source GO ontology version')
    url = 'https://release.geneontology.org/{}/ontology/go-basic.obo'.format(version.group(1))
    print('REFERENCE_DOWNLOAD ontology matching GAF',version.group(1),flush=True)
    fetched=subprocess.run(['curl','-fLsS','--retry','2','--max-time','180',url,'-o',str(root/'go-basic.obo')])
    if fetched.returncode:
        manifest['ontology_resolution_note']='Annotation-header ontology archive unavailable; frozen current ontology used; both versions recorded and unrecognized terms audited.'
        url='https://current.geneontology.org/ontology/go-basic.obo'
        subprocess.run(['curl','-fLsS','--retry','2','--max-time','180',url,'-o',str(root/'go-basic.obo')],check=True)
    manifest['go-basic.obo']={'url':url,'sha256':hashlib.sha256((root/'go-basic.obo').read_bytes()).hexdigest()}
    obo = (root/'go-basic.obo').read_text()
    manifest['go-basic.obo']['header'] = obo.split('[Term]')[0]
    terms, aliases = {}, {}
    for block in obo.split('[Term]')[1:]:
        fields = {}
        for line in block.splitlines():
            if ': ' in line:
                k, v = line.split(': ', 1)
                fields.setdefault(k, []).append(v)
        if 'id' not in fields or 'true' in fields.get('is_obsolete', []):
            continue
        term = fields['id'][0]
        parents = [x.split()[0] for x in fields.get('is_a', [])]
        parents += [x.split()[1] for x in fields.get('relationship', []) if x.startswith('part_of ')]
        terms[term] = {'name': fields['name'][0], 'parents': parents}
        for alt in fields.get('alt_id', []): aliases[alt] = term
    @functools.lru_cache(None)
    def ancestors(term):
        term = aliases.get(term, term)
        return frozenset([term]).union(*(ancestors(p) for p in terms.get(term, {}).get('parents', [])))
    for seed in sum(SEEDS.values(), []):
        if seed not in terms: raise ValueError('Missing GO seed '+seed)
    activity_terms = [t for t, v in terms.items() if any(phrase == v['name'].lower() for phrase in
        ['response to electrical stimulus', 'response to membrane depolarization',
         'response to neuronal activity', 'activity-dependent', 'cellular response to calcium ion'])]
    if activity_terms: SEEDS['activity_associated_GO'] = activity_terms
    manifest['resolved_category_seeds'] = {k: [{'id':t, 'name':terms[t]['name']} for t in v] for k,v in SEEDS.items()}
    sets, headers, unknown = {}, [], set()
    with gzip.open(root/'goa_human.gaf.gz', 'rt') as f:
        for line in f:
            if line.startswith('!'):
                headers.append(line.strip()); continue
            a = line.rstrip('\n').split('\t')
            if len(a)<15 or 'NOT' in a[3].split('|') or 'taxon:9606' not in a[12].split('|'): continue
            gene, term = a[2], aliases.get(a[4], a[4])
            if term not in terms: unknown.add(term);continue
            for t in ancestors(term):
                if any(set(seeds).intersection(ancestors(t)) for seeds in SEEDS.values()):
                    sets.setdefault(t, set()).add(gene)
    manifest['goa_human.gaf.gz']['header'] = headers
    manifest['unrecognized_annotation_terms_excluded']=sorted(unknown)
    rows = []
    for term, genes in sorted(sets.items()):
        cats = [cat for cat, seeds in SEEDS.items() if set(seeds).intersection(ancestors(term))]
        for gene in sorted(genes):
            rows.append([gene, term, terms[term]['name'], '|'.join(cats), 'GO_human', 'GO propagated via is_a/part_of; NOT excluded'])
    z = zipfile.ZipFile(root/'syngo1.3_complete_data.zip')
    for row in xlsx_rows(z.read('ontologies.xlsx')):
        name = row['name']
        cat = 'synaptic'
        if 'presynap' in name.lower(): cat += '|presynaptic'
        if 'postsynap' in name.lower(): cat += '|postsynaptic'
        for gene in row['hgnc_symbol'].split(','):
            if gene.strip(): rows.append([gene.strip(), 'SynGO:'+row['id'], name, cat, 'SynGO_1.3', 'Expert-curated; hierarchical rollup'])
    # Requested panels supplement, rather than replace, ontology-derived genes.
    # They are explicitly labelled targeted hypotheses, not independent curated discoveries.
    for panel, genes in PANELS.items():
        cat = 'activity' if panel in ('canonical_immediate_early','activity_associated_regulation') else 'maturation'
        for gene in genes:
            rows.append([gene, 'PANEL:'+panel, panel, cat, 'requested_targeted_panel_v1', 'User-requested panel; broad additions exploratory'])
    with (root/'functional_gene_sets.tsv').open('w') as f:
        w = csv.writer(f, delimiter='\t'); w.writerow(['gene','set_id','set_name','categories','source','evidence']);w.writerows(rows)
    (root/'reference_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print('REFERENCE_READY', len(rows), 'gene-set memberships', flush=True)

if __name__ == '__main__': prepare(sys.argv[1])
