#!/usr/bin/env python3
"""PI-facing display package; reuse frozen results without biological refitting."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import sys
import zipfile

PROJECT=Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder')
SOURCE=PROJECT/'results/div90_hypergate_sst_pv_phase2'
OUT=SOURCE/'pi_figure_package_v1'
REPO=Path(__file__).resolve().parents[2]
JOINED=PROJECT/'final_figures/div90_guidance_with_sst/tables/joined_with_sst.tsv.gz'
sys.dont_write_bytecode=True
os.environ['MPLCONFIGDIR']=str(OUT/'cache/matplotlib')
os.environ['NUMBA_CACHE_DIR']=str(OUT/'cache/numba')

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D

STATES=['PV-biased','PV/SST hybrid','SST-biased','unresolved/immature']
MAIN_STEMS=['Figure1_ERBB4_entry_population','Figure2_PV_SST_state_space','Figure3_old_gate_tradeoff','Figure4_prospective_sort_strategy']
PAGES=['01_PV_SST_STATE_MAP.html','02_VIRTUAL_FACS_GATE.html','03_CULTURE_CONDITION_COMPARISON.html']


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()


def dump(path,value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def setup():
    for p in ['main_figures','supplement','interactive','PI_summary','data','provenance','logs','cache']:
        (OUT/p).mkdir(parents=True,exist_ok=True)
    used=[SOURCE/n for n in ['cells.tsv.gz','gate_summary.json','state_parameters.json','gate_cell_assignments.tsv.gz','module_genes.tsv','hybrid_validation.json','hypergate_r_audit.json','REPORT.md','CONDITION_REPORT.md']]
    used.extend(p for p in (SOURCE/'tables').iterdir() if p.is_file() and p.suffix!='.h5')
    used.extend([JOINED,SOURCE/'provenance/DIV90_HYPERGATE_PHASE2_EVIDENCE.md'])
    path=OUT/'provenance/source_hashes_before.json'
    if not path.exists():dump(path,{str(p):sha(p) for p in sorted(used)})


def load():
    cells=pd.read_csv(SOURCE/'cells.tsv.gz',sep='\t',float_precision='round_trip')
    summary=json.loads((SOURCE/'gate_summary.json').read_text())
    params=json.loads((SOURCE/'state_parameters.json').read_text())
    assert len(cells)==4768 and cells.cell_id.is_unique
    assert cells.state.value_counts().to_dict()==dict(zip(STATES,[1076,1309,1075,1308]))
    assert int((cells.PVALB>0).sum())==3
    assert summary['phase1_reference']['retained_n']==3168
    assert summary['experimental_depletion']['retained_n']==4024
    return cells,summary,params


def configure():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':15,'axes.titlesize':17,'axes.labelsize':15,
        'xtick.labelsize':13,'ytick.labelsize':13,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none',
        'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})


def hide_umap(ax,d):
    ax.set_aspect('equal',adjustable='datalim')
    ax.margins(.08)
    ax.set_xticks([]);ax.set_yticks([])
    for spine in ax.spines.values():spine.set_visible(False)


def render_entry(cells,out=OUT):
    configure()
    joined=pd.read_csv(JOINED,sep='\t',float_precision='round_trip')
    cortex=joined[joined.set_id.eq('cortical_only')].copy()
    sub=joined[joined.set_id.eq('subpallial_only')].copy()
    assert cortex.cell_id.is_unique and sub.cell_id.is_unique
    saved=cortex.set_index('cell_id').loc[cells.cell_id]
    assert np.array_equal(saved[['loupe_x','loupe_y']].to_numpy(),cells[['loupe_x','loupe_y']].to_numpy())
    rows=[]
    for label,d in [('Cortical',cortex),('Subpallial',sub)]:
        valid=d.expression_available & d[['ERBB4','LHX6']].notna().all(axis=1)
        a=d.loc[valid]
        rows.append(dict(compartment=label,n_umap=len(d),n_expression_available=len(a),n_expression_missing=int((~valid).sum()),
            erbb4_detected=int(a.ERBB4.gt(0).sum()),erbb4_detected_fraction=float(a.ERBB4.gt(0).mean()),erbb4_mean=float(a.ERBB4.mean()),
            erbb4_median=float(a.ERBB4.median()),lhx6_detected=int(a.LHX6.gt(0).sum()),lhx6_detected_fraction=float(a.LHX6.gt(0).mean())))
    audit=pd.DataFrame(rows);audit.to_csv(out/'data/entry_population_audit.tsv',sep='\t',index=False)
    fig,axes=plt.subplots(2,2,figsize=(14,10))
    fig.subplots_adjust(left=.045,right=.90,bottom=.17,top=.86,wspace=.08,hspace=.30)
    fig.suptitle('ERBB4 identifies the cortical interneuron\npopulation we want to enrich',x=.045,y=.98,ha='left',fontsize=25,fontweight='bold')
    ax=axes[0,0]
    ax.scatter(cortex.loupe_x,cortex.loupe_y,s=5,c='#DCE2E8',linewidths=0)
    ax.scatter(cells.loupe_x,cells.loupe_y,s=5.5,c='#39876F',linewidths=0,alpha=.72)
    hide_umap(ax,cortex)
    ax.set_title('A  Recovered cortical recluster',loc='left',pad=10,fontweight='bold')
    ax.text(.02,.02,'4,768 entry cells highlighted',transform=ax.transAxes,fontsize=15,color='#27634F',
        bbox=dict(facecolor='white',edgecolor='none',alpha=.90,pad=3))
    cmap=LinearSegmentedColormap.from_list('PI_RNA',['#E8EDF2','#A1B5DD','#234D8C'])
    norm=Normalize(0,4)
    titles=[('B  Cortical ERBB4','ERBB4',cortex),('C  Cortical LHX6','LHX6',cortex),('D  Subpallial ERBB4','ERBB4',sub)]
    for ax,(title,gene,d) in zip([axes[0,1],axes[1,0],axes[1,1]],titles):
        available=d.expression_available & d[gene].notna()
        missing=d.loc[~available]
        ax.scatter(missing.loupe_x,missing.loupe_y,s=5,c='#CACED3',linewidths=0)
        a=d.loc[available].sort_values(gene)
        ax.scatter(a.loupe_x,a.loupe_y,c=a[gene],s=5,cmap=cmap,norm=norm,linewidths=0)
        hide_umap(ax,d)
        ax.set_title(title,loc='left',pad=10,fontweight='bold')
        if gene=='ERBB4':
            pct=100*a[gene].gt(0).mean()
            ax.text(.02,.02,f'{pct:.1f}% with detectable ERBB4 RNA',transform=ax.transAxes,fontsize=15,color='#34465C',
                bbox=dict(facecolor='white',edgecolor='none',alpha=.93,pad=3))
    cax=fig.add_axes([.925,.29,.017,.42])
    cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax,ticks=[0,4])
    if cb.solids is not None:cb.solids.set_rasterized(False)
    cb.ax.set_yticklabels(['Low','High']);cb.set_label('RNA abundance',labelpad=12)
    cb.outline.set_visible(False)
    fig.text(.045,.095,'Entry population for prospective sorting: cortical LHX6+/ERBB4+ cells',fontsize=17,fontweight='bold',color='#244D40')
    fig.text(.045,.049,'PVALB RNA is essentially absent at DIV90: detected in only 3 of 4,768 entry cells.',fontsize=15,color='#566371')
    for ext in ['pdf','svg','png']:
        fig.savefig(out/'main_figures'/f'Figure1_ERBB4_entry_population.{ext}',dpi=600,bbox_inches='tight',pad_inches=.12)
    plt.close(fig)
    caption=('Saved pre-entry cortical (10,028 cells) and subpallial (6,124 cells) local UMAPs were reused unchanged. '
        'Panel A overlays the frozen 4,768 entry cells on all cortical coordinates. RNA is available for 9,999 cortical and 6,121 subpallial cells; '
        '29 and 3 unavailable profiles respectively remain spatial context only and are excluded from detection denominators. '
        'Panels B–D use the same 0–4 log1p(CP10K) expression scale, with higher values clipped for display. '
        'The embeddings are separate saved reclusters; distances between compartments are not comparable. '
        'ERBB4 is enriched, not exclusive, in cortical cells. ERBB4-only live selection does not by itself recreate the recovered cortical/LHX6+ transcriptomic entry definition. '
        'PVALB was detected in 3/4,768 entry cells and did not define prospective PV identity. No expression, UMAP, or biological labels were recomputed.')
    dump(out/'provenance/entry_figure_validation.json',dict(n_entry=4768,n_pvalb_detected=3,entry_coordinates_exact=True,pre_entry_source=str(JOINED),
         displayed_compartments=audit.to_dict('records'),caption=caption,no_umap_recomputed=True))
    print(audit.to_string(index=False),flush=True)
    return caption


def guide():
    text='''# PI figure guide

**Central question:** Can we prospectively enrich cortical PV-associated precursors before PVALB is readily detectable?

This package uses the completed DIV90 phase-two results. PV means parvalbumin-associated; SST means somatostatin-associated. “Dual-high” denotes overlapping program scores, not a proven hybrid fate. Main figures present biological interpretation; the supplement contains computational details and culture comparisons.

## Figure 1 — ERBB4 entry population

- **One-sentence conclusion:** ERBB4 RNA is enriched in the recovered cortical compartment that supplies the 4,768-cell cortical LHX6+/ERBB4+ entry population.
- **What the PI should notice:** The saved cortical recluster has much stronger ERBB4 signal than the subpallial recluster; PVALB RNA is detected in only three entry cells.
- **One caveat:** ERBB4 is not cortex-exclusive, and ERBB4-only live selection does not automatically reproduce the recovered cortical/LHX6+ RNA selection.
- **15-second explanation:** “These are the original recovered populations. ERBB4 is much stronger in the cortical compartment. We start with 4,768 cortical LHX6+/ERBB4+ cells, because direct PVALB RNA detection identifies only three cells at this stage.”

## Figure 2 — Developmental state space

- **One-sentence conclusion:** PV-associated and SST-associated developmental information is measurable at DIV90 before PVALB becomes broadly detectable.
- **What the PI should notice:** Separate program axes retain PV-biased, dual-high, SST-biased and unresolved cells; 27.5% fall in the operational dual-high region.
- **One caveat:** The quadrants are operational score partitions, and coordinated hybrid identity has not been established for the whole dual-high population.
- **15-second explanation:** “PVALB is scarce, but developmental programs are already measurable. We keep PV and SST as separate axes so that overlapping programs remain visible. About 28% are dual-high; that is an experimental population to follow, not an established hybrid fate.”

## Figure 3 — The first gate's tradeoff

- **One-sentence conclusion:** FAT3-low/PTPRM-undetected is too restrictive for a strategy that aims to preserve overlapping developmental programs.
- **What the PI should notice:** The old gate retains 79.1% of PV-biased cells but only 47.9% of dual-high cells, losing 682 dual-high cells.
- **One caveat:** This is a different biological objective from the earlier binary analysis; the old gate's original performance and all cell labels remain preserved.
- **15-second explanation:** “The first gate enriched the PV-biased corner, but it discarded more than half of the dual-high cells. That makes it too restrictive if those overlapping programs are part of the population we want to mature.”

## Figure 4 — Prospective sort-and-mature experiment

- **One-sentence conclusion:** Removing cells high for either FGFR2 or PTPRS preserves about 89% of PV-biased and dual-high cells and creates a testable prospective experiment.
- **What the PI should notice:** Retain the complement—both markers low—and collect both fractions plus an unsorted ERBB4+ control for identical maturation and later protein, transcriptomic and electrophysiological readouts.
- **One caveat:** These are RNA-derived predictions; SST-biased contamination changes only modestly, from 22.5% to 21.3%, and neither a protein fluorescence threshold nor a validated PV sort is established.
- **15-second explanation:** “The practical candidate removes FGFR2-high or PTPRS-high cells while preserving roughly 89% of both target populations. Collect both fractions and an unsorted control, mature them identically, and ask whether the retained cells develop the desired PV phenotype.”

## Reading sequence

Show Figures 1 through 4, then the PI summary page. Keep culture dependence, gene-level support, technical diagnostics and exact thresholds in the supplement. Do not frame the small SST contamination change as dramatic purification.

## Interactive pages

- [Developmental state map](interactive/01_PV_SST_STATE_MAP.html): default view shows all 4,768 cells by operational state.
- [Virtual FGFR2/PTPRS gate](interactive/02_VIRTUAL_FACS_GATE.html): the RNA-derived removed and retained regions, with composition and recovery.
- [Culture-condition comparison](interactive/03_CULTURE_CONDITION_COMPARISON.html): defaults to H9; compare CV and MW separately within each line. The 2E/CV sample contains only 37 entry cells.

All three pages are standalone and work offline. Exact gate thresholds and computational details are in supplemental captions, not the main figures. The supplementary culture comparisons treat CV/MW as a biological condition including the reported glucose difference; operator/culture confounding prevents a glucose-specific causal claim.
'''
    (OUT/'PI_FIGURE_GUIDE.md').write_text(text)
    entry=json.loads((OUT/'provenance/entry_figure_validation.json').read_text())
    captions='# Main figure captions and display provenance\n\n## Figure 1\n\n'+entry['caption']+'\n'
    (OUT/'provenance/Figure1_caption.md').write_text(captions)
    manifest=OUT/'provenance/main_figures_2_4_manifest.json'
    if manifest.exists():
        for row in json.loads(manifest.read_text()):
            captions+='\n## '+row['figure']+'\n\n'+row['caption']+'\n'
    (OUT/'MAIN_FIGURE_CAPTIONS.md').write_text(captions)


def support_files():
    """Package the existing evidence and exact rules alongside display captions."""
    (OUT/'evidence').mkdir(exist_ok=True)
    shutil.copy2(SOURCE/'provenance/DIV90_HYPERGATE_PHASE2_EVIDENCE.md',OUT/'evidence/DIV90_HYPERGATE_PHASE2_EVIDENCE.md')
    summary=json.loads((SOURCE/'gate_summary.json').read_text())
    rules={key:summary[key] for key in ['phase1_reference','experimental_depletion']}
    dump(OUT/'data/exact_gate_rules.json',rules)
    exact=json.loads(rules['experimental_depletion']['rules'])
    marker='\n\n## Full-precision gate evaluation\n\n'
    caption=OUT/'supplement/SUPPLEMENTAL_CAPTIONS.md'
    base=caption.read_text().split(marker)[0]
    rule=' OR '.join(f"{r['gene']} {r['op']} {r['threshold']!r}" for r in exact)
    caption.write_text(base+marker+'The exact saved practical rule is **REMOVE '+rule+'**; retain cells at or below BOTH thresholds. '
        'Values are log1p(CP10K) RNA, not protein fluorescence. The old retention rule is FAT3 ≤ 0.764341 AND PTPRM ≤ 0. '
        'Use `../data/exact_gate_rules.json` to reproduce the rules; shorter numbers elsewhere are display rounding.\n')


def verify_sources():
    before=json.loads((OUT/'provenance/source_hashes_before.json').read_text())
    changed=[name for name,digest in before.items() if sha(Path(name))!=digest]
    assert not changed,changed
    dump(OUT/'provenance/source_preservation.json',dict(source_files_unchanged=True,n_files=len(before),algorithm='SHA256',no_biological_reanalysis=True))


def package():
    import xml.etree.ElementTree as ET
    from PIL import Image
    archive=OUT/'DIV90_Hypergate_PI_Figure_Package.zip'
    if archive.exists():raise FileExistsError('Keep the existing published package; use a new version for further edits.')
    setup();guide();support_files();verify_sources()
    required=[OUT/'main_figures'/f'{s}.{ext}' for s in MAIN_STEMS for ext in ['pdf','svg','png']]
    required.extend([OUT/'supplement/supplemental_figures.pdf',OUT/'PI_summary/PI_summary_one_page.pdf',OUT/'PI_summary/PI_summary_one_page.png',OUT/'PI_FIGURE_GUIDE.md'])
    required.extend(OUT/'interactive'/p for p in PAGES)
    supplement=json.loads((OUT/'supplement/supplement_manifest.json').read_text())
    assert supplement['pages']==10 and not supplement['new_biological_fits']
    for row in supplement['figures']:required.extend([OUT/row['svg'],OUT/row['png']])
    document_checks=[]
    for p in required:
        assert p.is_file() and p.stat().st_size>1000,p
        if p.suffix=='.pdf':
            raw=p.read_bytes()
            assert raw.startswith(b'%PDF-') and raw.rstrip().endswith(b'%%EOF'),p
            pages=len(re.findall(rb'/Type\s*/Page\b',raw))
            assert pages==(10 if p.name=='supplemental_figures.pdf' else 1),(p,pages)
            raster_images=len(re.findall(rb'/Subtype\s*/Image\b',raw))
            if p.parent.name=='main_figures':assert raster_images==0,p
            document_checks.append(dict(path=str(p.relative_to(OUT)),pages=pages,raster_images=raster_images))
        elif p.suffix=='.svg':
            root=ET.parse(p).getroot()
            if p.parent.name=='main_figures':assert not root.findall('.//{http://www.w3.org/2000/svg}image'),p
        elif p.suffix=='.png':
            with Image.open(p) as im:assert min(im.info.get('dpi',(0,0)))>=599;im.verify()
    assert sorted(p.name for p in (OUT/'interactive').glob('*.html'))==PAGES
    assert len(list((OUT/'main_figures').glob('*.pdf')))==4
    browser=json.loads((OUT/'provenance/interactive_validation.json').read_text())
    assert browser['offline'] and browser['all_network_blocked']
    assert [r['filename'] for r in browser['assets']]==PAGES
    for row in browser['assets']:
        assert row['hover'] and row['filtering'] and row['reset_and_zoom']
        assert not row['javascript_errors'] and not row['external_requests']
    gate=browser['assets'][1]['default_counts']
    assert (gate['n'],gate['removed'],gate['retained'])==(4768,744,4024)
    assert gate['retainedCounts']==[954,1169,856,1045]
    assert gate['removedCounts']==[122,140,219,263]
    assert browser['assets'][2]['condition_counts']['2E']['CV']['n']==37
    for src in [Path(__file__),*Path(__file__).parent.glob('div90_hypergate_pi_*.py')]:
        shutil.copy2(src,OUT/'provenance'/src.name)
    for name in ['DIV90_HYPERGATE_PHASE2_EVIDENCE.md','culture_condition_user_20260909.tsv']:
        src=SOURCE/'provenance'/name
        if src.exists():shutil.copy2(src,OUT/'provenance'/name)
    attachment=Path('/home/elcrespo/.codex/attachments/665117ba-f473-43c4-8cd1-7f71073c5deb/pasted-text.txt')
    shutil.copy2(attachment,OUT/'provenance/PI_package_request.txt')
    dump(OUT/'PACKAGE_VALIDATION.json',dict(main_figures=4,main_formats=['PDF','SVG','600dpi PNG'],interactive_pages=3,pi_summary_pages=1,supplement='supplement/supplemental_figures.pdf',
        original_results_unchanged=True,source_preservation='provenance/source_preservation.json',document_checks=document_checks,
        offline_browser_checks_passed=True,required_files=[str(p.relative_to(OUT)) for p in required]))
    files=[]
    for directory,subdirs,names in os.walk(OUT):
        subdirs[:]=[d for d in subdirs if d not in ['cache','logs','supplement_previews']]
        files.extend(Path(directory)/n for n in names if n not in ['FILE_MANIFEST.json','zip_validation.json'] and not n.endswith('.zip'))
    files.sort()
    dump(OUT/'FILE_MANIFEST.json',dict(files=[dict(path=str(p.relative_to(OUT)),bytes=p.stat().st_size,sha256=sha(p)) for p in files]))
    files.append(OUT/'FILE_MANIFEST.json')
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=5) as z:
        for p in files:z.write(p,p.relative_to(OUT))
    with zipfile.ZipFile(archive) as z:assert z.testzip() is None
    dump(OUT/'provenance/zip_validation.json',dict(file=str(archive),bytes=archive.stat().st_size,sha256=sha(archive),crc_verified=True))
    print(f'Created {archive} ({archive.stat().st_size/1e6:.1f} MB)',flush=True)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--stage',choices=['entry','guide','all','package'],default='all');args=ap.parse_args()
    setup();cells,summary,params=load()
    if args.stage in ['entry','all']:render_entry(cells)
    if args.stage=='all':
        from div90_hypergate_pi_main import run as main_run
        from div90_hypergate_pi_supplement import run as supplement_run
        from div90_hypergate_pi_interactive import run as interactive_run, verify_browser
        main_run(cells,OUT,summary,params)
        supplement_run(cells,OUT,summary,params)
        interactive_run(cells,OUT,summary,params)
        verify_browser(out=OUT,cells=cells)
    if args.stage in ['guide','all','package']:guide()
    if args.stage=='package':package()
    verify_sources()


if __name__=='__main__':main()
