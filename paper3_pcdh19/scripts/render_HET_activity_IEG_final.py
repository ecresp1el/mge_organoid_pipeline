"""Extend the display-clarity figure with an exploratory paired IEG bottom row."""
from pathlib import Path
import csv
import hashlib
import json
import shutil
from pcdh19_exploratory_paths import RUNTIME_ROOT

HERE = Path(__file__).resolve().parent
OUT = RUNTIME_ROOT / 'final_figures/mge_activity_IEG_dots_20261101'
OUT.mkdir(parents=True, exist_ok=True)
SOURCE = RUNTIME_ROOT / 'results/female_male_signaling_20261031'
IEG = RUNTIME_ROOT / 'results/postmitotic_inhibitory_IEGs_20261101'
inputs = ['atlas_sex_genotype_display.tsv.gz', 'male_total_Pcdh19_display.npz',
          'male_C0006_DE.tsv', 'matched_display_logCPM.tsv',
          'matched_female_display_metadata.tsv', 'saved_MGE_LGE_CGE_scores.npz',
          'HET_replicated_interaction_selection.tsv']
for name in inputs:
    shutil.copy2(str(SOURCE / name), str(OUT / name))
# Redirect all base-renderer exports into this new package, preserving prior figures.
base = HERE / 'render_HET_display_clarity.py'
source = base.read_text()
needle = "OUT=RUNTIME_ROOT/'results/female_male_signaling_20261031'"
assert source.count(needle) == 1
source = source.replace(needle, 'OUT=Path(' + repr(str(OUT)) + ')')
ns = {'__file__': str(base), '__name__': '__activity_base__'}
exec(compile(source, str(base), 'exec'), ns)
fig, plt, Rectangle = ns['fig'], ns['plt'], ns['Rectangle']
# Keep the established four rows; fit G/H/I across the last row.
import numpy as np
import h5py
height = 8.49
for ax in fig.axes:
    if getattr(ax, '_figure_row', None) == 3:
        p = ax.get_position(original=True)
        if p.x0 < .44:
            ax.set_position([p.x0 * .75, p.y0, p.width * .75, p.height])
        else:
            ax.set_position([.325 + (p.x0 - .435) * .60, p.y0, p.width * .60, p.height])
        for tx in ax.texts:
            tx.set_fontsize(max(5, tx.get_fontsize() * .76))
        ax.xaxis.label.set_fontsize(5)
        ax.tick_params(labelsize=5)
for tx in fig.texts:
    if getattr(tx, '_figure_row', None) == 3:
        x, y = tx.get_position()
        if x < .435:
            tx.set_position((x * .75, y))
        else:
            tx.set_position((.325 + (x - .435) * .60, y))
        if tx.get_text() == 'HET candidate network':
            tx.set_text('Ligand–receptor network'); tx.set_fontsize(7)
        elif tx.get_text().startswith('KO cells vs WT cells within'):
            tx.set_text('Ligand–receptor changes'); tx.set_fontsize(7); tx.set_x(.36)
        elif tx.get_text() not in ['G', 'H']:
            tx.set_fontsize(5.5)
for artist in fig.artists:
    if isinstance(artist, Rectangle) and getattr(artist, '_figure_row', None) == 3:
        artist.set_x(.325 + (artist.get_x() - .435) * .60)
        artist.set_width(artist.get_width() * .60)
# Compute detection fractions directly from the same cached raw IEG counts and gate.
def strings(ds):
    return np.array([v.decode() if isinstance(v, bytes) else str(v) for v in ds[:]])
base_results = RUNTIME_ROOT / 'results'
with h5py.File(str(base_results / 'mosaic_deep_dive_20261031/per_cell_marker_diagnostics.h5'), 'r') as f:
    ids = strings(f['cell_id']); clusters = strings(f['cluster']); samples = strings(f['sample'])
    total = f['total_counts'][:]; markers = list(strings(f['marker_gene']))
    raw = f['marker_raw_counts'][:]; gates = f['gate_membership'][:]
    gate_names = list(strings(f['gate_names']))
cc = np.load(str(base_results / 'nested_proliferation_20261101/raw_cycle_markers.npz'))
gc = np.load(str(base_results / 'IEG_neuronal_maturation_20261101/all_cell_IEG_counts.npz'))
assert np.array_equal(cc['cell_id'], ids) and np.array_equal(gc['cell_id'], ids)
core = (raw[:, [markers.index(g) for g in ['Snap25', 'Syt1', 'Rbfox3']]] > 0).sum(1)
gaba = raw[:, markers.index('Gad1')] + raw[:, markers.index('Gad2')] > 0
cycle_genes = list(cc['genes'])
mask = np.isin(clusters, ['C0004','C0005','C0006','C0020']) & gaba & (core >= 1)
mask &= (cc['counts'][:, cycle_genes.index('Mki67')] == 0) & (cc['counts'][:, cycle_genes.index('Top2a')] == 0)
het = np.isin(samples, ['15662-JZ-7','15662-JZ-8','15662-JZ-9'])
genes = ['Npas4','Fos','Fosb','Junb','Egr1','Egr2','Egr3','Arc','Nr4a1']; all_genes = list(gc['genes'])
groups = [
    ('HET WT cells', het & gates[:, gate_names.index('A_positive')].astype(bool)),
    ('HET KO cells', het & gates[:, gate_names.index('Aneg_BC1')].astype(bool)),
    ('WT males', np.isin(samples, ['15662-JZ-1','15662-JZ-2','15662-JZ-3'])),
    ('KO males', np.isin(samples, ['15662-JZ-10','15662-JZ-11','15662-JZ-12']))]
values = []; per_sample = []
for group, group_mask in groups:
    selected = mask & group_mask
    assert selected.any(), group
    for gene in genes:
        counts = gc['counts'][selected, all_genes.index(gene)]
        values.append(dict(group=group, gene=gene, cells=int(selected.sum()), raw_UMI=int(counts.sum()),
                           detected_pct=float((counts > 0).mean() * 100),
                           mean_log_expression=float(np.log1p(counts / total[selected] * 10000).mean())))
        for sample in sorted(set(samples[selected])):
            sel = selected & (samples == sample)
            c = gc['counts'][sel, all_genes.index(gene)]
            per_sample.append(dict(group=group,sample=sample,gene=gene,cells=int(sel.sum()),
                                   raw_UMI=int(c.sum()),detected_pct=float((c>0).mean()*100),
                                   mean_log_expression=float(np.log1p(c/total[sel]*10000).mean())))
assert [v['cells'] for v in values[:18]] == [16]*len(genes) + [151]*len(genes)
# Display each sample equally; scale each gene independently across four groups.
for v in values:
    sample_values = [r['mean_log_expression'] for r in per_sample
                     if r['group'] == v['group'] and r['gene'] == v['gene']]
    v['equal_sample_mean_log_expression'] = float(np.mean(sample_values))
ranked = []
for gene in genes:
    entries = [next(v for v in values if v['gene'] == gene and v['group'] == group)
               for group, _ in groups]
    maximum = max(v['equal_sample_mean_log_expression'] for v in entries)
    scaled = [v['equal_sample_mean_log_expression'] / maximum if maximum else 0.
              for v in entries]
    score = .5 * ((scaled[1] - scaled[0]) + (scaled[3] - scaled[2]))
    ranked.append((gene, score, maximum))
    for v, relative in zip(entries, scaled):
        v['relative_expression_0_1'] = relative
        v['ranking_score'] = score
        v['circle_area_pt2'] = 6 + 30 * relative if relative > 0 else 5
# Descending mean of female and male KO-minus-WT relative-expression changes.
genes = [g for g, score, maximum in sorted(ranked, key=lambda r: (r[2] == 0, -r[1], r[0]))]
for v in values:
    v['row_rank'] = genes.index(v['gene']) + 1
# Update only the C/E display wording; underlying selections remain unchanged.
for axes in fig.axes:
    if getattr(axes, '_figure_row', None) == 1 and axes.get_xlabel() == 'Log-normalized expression':
        axes.set_yticklabels(['WT Female','HET Female WT cells','HET Female KO cells'])
        for text in list(axes.texts):
            if '0 exon-1 signal' in text.get_text():
                text.remove()
    if getattr(axes, '_figure_row', None) == 2 and axes.get_xlabel() == 'Log-normalized expression':
        axes.set_yticklabels(['WT Male','KO Male'])
# Rebalance H to leave a wider, integrated I panel in row four.
for axes in fig.axes:
    if getattr(axes, '_figure_row', None) == 3 and axes.get_position().x0 >= .32:
        p = axes.get_position(original=True)
        axes.set_position([.31 + (p.x0-.325)*.85, p.y0, p.width*.85, p.height])
for text in fig.texts:
    if getattr(text, '_figure_row', None) == 3 and text.get_position()[0] >= .32:
        x,y=text.get_position();text.set_position((.31+(x-.325)*.85,y))
for artist in fig.artists:
    if isinstance(artist, Rectangle) and getattr(artist,'_figure_row',None)==3:
        artist.set_x(.31+(artist.get_x()-.325)*.85);artist.set_width(artist.get_width()*.85)
fig.text(.045,.07/height,'Arrows: ligand source → receptor receiver',fontsize=5.5)
fig.text(.34,1.64/height,'Female HET · KO cells vs WT cells',fontsize=6)
ns['barax'].set_xlabel('Inferred L–R change (%)',fontsize=6)
fig.text(.625, 1.82 / height, 'I', fontsize=11, weight='bold')
fig.text(.660, 1.82 / height, 'Activity-dependent gene expression', fontsize=8, weight='bold')
fig.text(.65, 1.64 / height, 'Cycling-marker-negative', fontsize=6)
fig.text(.65, 1.51 / height, 'inhibitory neurons', fontsize=6)
ax = fig.add_axes([.735, .49 / height, .14, .87 / height])
fig.text(.77,1.41/height,'HET F',ha='center',fontsize=5.5)
fig.text(.84,1.41/height,'Males',ha='center',fontsize=5.5)
vmax = 1.
for v in values:
    x = [g[0] for g in groups].index(v['group']); y = genes.index(v['gene'])
    if v['relative_expression_0_1'] == 0:
        ax.scatter([x],[y],s=5,facecolors='none',edgecolors='#BFC4CC',linewidths=.5)
    else:
        ax.scatter([x],[y],s=v['circle_area_pt2'],c=[v['relative_expression_0_1']],
                   cmap='Purples',vmin=0,vmax=vmax,edgecolors='#746284',linewidths=.4)
ax.set_xlim(-.6,3.6); ax.set_ylim(len(genes)-.45,-.55)
ax.set_xticks([0,1,2,3]); ax.set_xticklabels(['WT','KO','WT','KO'],fontsize=6)
ax.axvline(1.5,color='#E1E1E1',lw=.6)
ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes,fontstyle='italic',fontsize=6.5)
ax.tick_params(length=0,pad=3)
for spine in ax.spines.values(): spine.set_visible(False)
# Add aligned change columns without replacing the four expression columns.
delta_ax = fig.add_axes([.915,.49/height,.07,.87/height])
delta_ax.set_xlim(-.5,1.5);delta_ax.set_ylim(len(genes)-.45,-.55)
delta_ax.set_xticks([0,1]);delta_ax.set_xticklabels(['HET F','Male'],rotation=0,ha='center',fontsize=6)
delta_ax.set_yticks([]);delta_ax.tick_params(length=0,pad=3)
for spine in delta_ax.spines.values():spine.set_visible(False)
fig.text(.95,1.41/height,'KO−WT',ha='center',fontsize=6)
changes=[]
for gene in genes:
    entries=[next(v for v in values if v['gene']==gene and v['group']==group)
             for group,_ in groups]
    for x,(name,wt,ko) in enumerate([('HET female',entries[0],entries[1]),('Male',entries[2],entries[3])]):
        delta=ko['relative_expression_0_1']-wt['relative_expression_0_1']
        changes.append(dict(gene=gene,comparison=name,scaled_KO_minus_WT=delta,
                            mean_log_expression_change=ko['equal_sample_mean_log_expression']-wt['equal_sample_mean_log_expression'],
                            row_rank=genes.index(gene)+1))
        if delta==0:
            delta_ax.scatter([x],[genes.index(gene)],s=5,facecolors='none',edgecolors='#BFC4CC',linewidths=.5)
        else:
            delta_ax.scatter([x],[genes.index(gene)],s=6+30*abs(delta),c=[delta],
                             cmap='RdBu_r',vmin=-1,vmax=1,edgecolors='#888888',linewidths=.35)
# Separate keys distinguish within-gene expression from signed expression change.
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
for rect,cmap,low,high in [([.75,.09/height,.10,.045/height],'Purples',0,1),
                           ([.92,.09/height,.05,.045/height],'RdBu_r',-1,1)]:
    scale_ax=fig.add_axes(rect)
    cb=fig.colorbar(ScalarMappable(norm=Normalize(low,high),cmap=cmap),cax=scale_ax,
                    orientation='horizontal',ticks=[])
    cb.outline.set_visible(False)
fig.text(.75,.24/height,'Expression',fontsize=6)
fig.text(.915,.24/height,'Change',fontsize=6)
fig.text(.73,.08/height,'0',fontsize=6);fig.text(.86,.08/height,'1',fontsize=6)
fig.text(.895,.08/height,'−1',fontsize=6);fig.text(.975,.08/height,'+1',fontsize=6)
# A common typographic hierarchy throughout the four-row figure.
from matplotlib.text import Text
for text in fig.findobj(match=Text):
    if not text.get_text():
        continue
    if getattr(text, '_panel_heading', False):
        text.set_fontsize(11 if text.get_text() in list('ABCDEFGHI') else 8)
    elif text.get_text() in ['G', 'H', 'I'] and text in fig.texts:
        text.set_fontsize(11)
    elif text.get_text() in ['Ligand–receptor network','Ligand–receptor changes','Activity-dependent gene expression']:
        text.set_fontsize(8)
    else:
        text.set_fontsize(min(7, max(6, text.get_fontsize())))
for text in ns['labelsax'].texts:
    if text.get_text().startswith('Sema4a'):
        text.set_text('Sema4a–Nrp1/\nPlxna4')
        text.set_fontsize(5.5)
# Give the compressed network labels room around their arrows and nodes.
network = ns['main_network']
for text in network.texts:
    text.set_fontsize(5.5)
    if text.get_text().startswith('Pcdhb21'):
        text.set_position((1.45, .65))
    elif text.get_text().startswith('Ncam2'):
        text.set_x(-1.60)
# Preserve the tightly spaced atlas legend and inset hierarchy.
for axes in fig.axes:
    legend = axes.get_legend()
    if legend is not None:
        for text in legend.get_texts():
            text.set_fontsize(5.5)
for axes in fig.axes:
    if getattr(axes, '_figure_row', None) == 1 and axes.get_xlabel() == 'Log-normalized expression':
        p=axes.get_position(original=True)
        axes.set_position([p.x0+.045,p.y0,p.width-.045,p.height])
        for label in axes.get_yticklabels():
            label.set_fontsize(6)
# Final user-approved sex wording, applied to every visible text artist.
for text in fig.findobj(match=Text):
    label=text.get_text()
    for old,new in [('Females','Female'),('females','female'),('Males','Male'),('males','male')]:
        label=label.replace(old,new)
    text.set_text(label)
fig.canvas.draw()
renderer=fig.canvas.get_renderer()
scale_text=[text for text in fig.texts if text.get_text() in ['Expression','Change','0','1','−1','+1']]
label_boxes=[text.get_window_extent(renderer).padded(2) for axes in [ax,delta_ax] for text in axes.get_xticklabels()]
assert not any(box.overlaps(text.get_window_extent(renderer)) for box in label_boxes for text in scale_text), 'Panel I labels overlap scale text'
(OUT/'PANEL_I_SPACING_VALIDATION.json').write_text(json.dumps({'group_labels_vs_scale_text_no_overlap':True,'clearance_pixels':2},indent=2))
with (OUT / 'IEG_dot_values.tsv').open('w') as f:
    writer=csv.DictWriter(f,fieldnames=list(values[0]),delimiter='\t');writer.writeheader();writer.writerows(values)
with (OUT / 'IEG_dot_per_sample_values.tsv').open('w') as f:
    writer=csv.DictWriter(f,fieldnames=list(per_sample[0]),delimiter='\t');writer.writeheader();writer.writerows(per_sample)
with (OUT / 'IEG_KO_minus_WT_changes.tsv').open('w') as f:
    writer=csv.DictWriter(f,fieldnames=list(changes[0]),delimiter='\t');writer.writeheader();writer.writerows(changes)
rows = values
for ext in ['png', 'pdf', 'svg']:
    fig.savefig(str(OUT / ('HET_activity_IEG_final.' + ext)), dpi=450, facecolor='white')
legend = (OUT / 'HET_DISPLAY_CLARITY_LEGEND.txt').read_text()
legend += '\nC/E display wording updated: C WT Female / HET Female WT cells / HET Female KO cells; C zero-signal text removed. E WT Male / KO Male. All visible sex labels use singular wording. Underlying cell selections, measurements and fits remain unchanged.\n'
legend += '\nI: Nine prespecified cached neuronal IEGs (Npas4, Fos, Fosb, Junb, Egr1, Egr2, Egr3, Arc, Nr4a1), retained regardless of detection or effect direction. Size legend removed at user request; detection is no longer used for circle size. Exploratory activity-associated IEG expression, not measured intracellular pathway activation. Existing inhibitory clusters C0004/5/6/20; Gad1 or Gad2 detected; >=1 Snap25/Syt1/Rbfox3 detected; zero Mki67 and Top2a UMIs. No IEG-based gate. Dot plot shows four separate columns: probe-defined WT (A-positive) and KO (A-negative, B+C>=1) cells from three HET females, then known WT males and known KO males (three samples each). The same inhibitory identity/neuronal-core/cycle gate is applied to all columns. Male columns retain all gate-eligible neurons without a Pcdh19 expression filter; unlike the positive-only Pcdh19 violin display in E, no male expression selection is applied. Each column uses equally weighted sample means, rather than pooling cells for expression display; pooled detection and raw count summaries are retained as unplotted diagnostics. IEG_dot_per_sample_values.tsv preserves sample-level values. The exploratory footer is removed from the artwork at user request; limitations remain here. Circle size and purple color now encode relative expression, NOT detection. Compute each sample mean ln(1 + gene UMI / total UMI * 10000), including zeros, then average sample means equally within each of four groups. Divide each gene by its maximum across the four group means (0–1 scale); all-zero genes remain zero. Positive dot area = 6 + 30 * relative expression pt^2; zero uses a 5 pt^2 hollow reference circle. Row order descends by the mean of the female and male KO-minus-WT scaled expression differences, weighting both comparisons equally, with alphabetical tie-breaking; all-zero genes are placed last. Ranking is descriptive and effect-selected, not significance. Absolute expression cannot be compared across gene rows. Detection percentages remain in the TSV as unplotted source values. Displayed rows are pooled descriptive summaries, not independent replicate inference. Npas4 has 0 WT/9 KO UMIs, from 16 WT/151 KO cells; smallest WT sample has two cells. Relaxing the cycling-marker gate to <=1 UMI weakens the pattern (mean log-expression higher in KO in only 1/3 embryos). No significance or mechanistic claim. Cycling-marker absence is not direct proof of postmitotic state. E14.5 dissected MGE does not establish cortical residence.\n'
legend += '\nI right-hand change columns supplement the four expression columns. For each gene, subtract WT from KO relative expression separately in HET females and males. Red indicates higher KO, blue lower KO; area=6+30*abs(change) pt2, with hollow reference for zero. Both change columns use a shared -1 to +1 scale; this is a difference in gene-relative expression, not log fold change or a significance statistic. Gene rows align across expression/change views. Title Activity-dependent gene expression describes the gene panel, not an experimentally demonstrated stimulation response.\n'
(OUT / 'HET_ACTIVITY_IEG_LEGEND.txt').write_text(legend)
for p in [Path(__file__), base, IEG / 'paired_postmitotic_IEGs.tsv', IEG / 'postmitotic_subset_summary.tsv', IEG / 'PLAN.json', IEG / 'REPORT.md']:
    shutil.copy2(str(p), str(OUT / p.name))
manifest = {'figure_size_inches': [7.2, height], 'panel_rows': [['A','B'],['C','D'],['E','F'],['G','H','I']],
            'IEG_change_columns': ['HET female KO−WT','Male KO−WT'], 'IEG_size_encoding': 'per-gene 0–1 relative equal-sample mean log expression; positive area 6+30*value pt2', 'IEG_ranking': 'descending average of female and male KO-minus-WT scaled changes', 'IEG_groups': [g[0] for g in groups], 'male_Pcdh19_expression_filter': False, 'IEG_genes': genes, 'size_key_removed': True, 'typography_pt': {'panel_letter': 11, 'panel_title': 8, 'body': '6–7', 'atlas_legend': 5.5}, 'IEG_source': str(IEG), 'paired_rows': len(rows), 'exploratory': True,
            'files': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.iterdir() if p.is_file() and p.name != 'PROVENANCE.json'}}
(OUT / 'PROVENANCE.json').write_text(json.dumps(manifest, indent=2))
print('ACTIVITY_IEG_FINAL_COMPLETE', str(OUT), flush=True)
