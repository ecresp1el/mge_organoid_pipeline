#!/usr/bin/env python3
"""PI Figures 2–4 and summary, drawn only from frozen DIV90 phase-two outputs.

Runtime artifacts and caches belong to the new Turbo figure-package directory.
This module does not rescore cells, refit gates, or alter source artifacts.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

ROOT = Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_hypergate_sst_pv_phase2')
OUT = ROOT / 'pi_figure_package_v1'
os.environ.setdefault('MPLCONFIGDIR', str(OUT / 'cache/matplotlib'))
os.environ.setdefault('XDG_CACHE_HOME', str(OUT / 'cache'))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle
import numpy as np
import pandas as pd

STATES = ['PV-biased', 'PV/SST hybrid', 'SST-biased', 'unresolved/immature']
LABELS = ['PV-biased', 'Dual-high', 'SST-biased', 'Unresolved']
PALETTE = ['#326DAB', '#39876F', '#CB733F', '#A6ACB6']
INK = '#233442'
MUTED = '#596670'
LINE = '#CFD7DD'


def configure():
    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 13, 'axes.labelsize': 14,
        'axes.titlesize': 15, 'axes.titleweight': 'bold', 'axes.spines.top': False,
        'axes.spines.right': False, 'axes.linewidth': .8, 'xtick.labelsize': 12,
        'ytick.labelsize': 12, 'text.color': INK, 'axes.labelcolor': INK,
        'xtick.color': MUTED, 'ytick.color': MUTED, 'pdf.fonttype': 42,
        'ps.fonttype': 42, 'svg.fonttype': 'none', 'figure.facecolor': 'white',
        'savefig.facecolor': 'white', 'axes.facecolor': 'white',
    })


def save(fig, out, relative_stem, caption, extensions=('pdf', 'svg', 'png')):
    stem = out / relative_stem
    stem.parent.mkdir(parents=True, exist_ok=True)
    for artist in fig.findobj():
        if hasattr(artist, 'get_rasterized') and artist.get_rasterized():
            artist.set_rasterized(False)
    for extension in extensions:
        fig.savefig(stem.with_suffix('.' + extension), dpi=600,
                    facecolor='white', bbox_inches=None)
    # Review-sized exports keep inspection practical; these are not deliverables.
    preview = out / 'cache/previews' / (stem.name + '.png')
    preview.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(preview, dpi=110, facecolor='white')
    plt.close(fig)
    return {'figure': stem.name, 'files': [str(Path(relative_stem).with_suffix('.' + x)) for x in extensions],
            'caption': caption, 'preview': str(preview.relative_to(out))}


def apply_rule(df, rule):
    pieces = []
    for item in rule['rules']:
        values = df[item['gene']].to_numpy()
        threshold = float(item['threshold'])
        pieces.append({'<=': np.less_equal, '<': np.less, '>=': np.greater_equal,
                       '>': np.greater, '=': np.equal, '==': np.equal}[item['op']](values, threshold))
    hit = np.logical_or.reduce(pieces) if rule['logic'] == 'OR' else np.logical_and.reduce(pieces)
    return ~hit if rule['action'] == 'remove' else hit


def frozen_audit(df, summary, params):
    assert len(df) == 4768 and df.cell_id.is_unique
    counts = df.state.value_counts().reindex(STATES).to_numpy()
    assert counts.tolist() == [1076, 1309, 1075, 1308]
    high_pv = df.pv_score >= params['pv_threshold']
    high_sst = df.sst_score >= params['sst_threshold']
    expected = np.select([high_pv & ~high_sst, high_pv & high_sst, ~high_pv & high_sst],
                         STATES[:3], default=STATES[3])
    assert np.array_equal(expected, df.state.to_numpy())
    old = apply_rule(df, summary['rules'][str(summary['phase1_reference']['gate_id'])])
    candidate = apply_rule(df, summary['rules'][str(summary['experimental_depletion']['gate_id'])])
    assert old.sum() == 3168 and candidate.sum() == 4024
    old_counts = df.loc[old, 'state'].value_counts().reindex(STATES).to_numpy()
    kept_counts = df.loc[candidate, 'state'].value_counts().reindex(STATES).to_numpy()
    assert old_counts.tolist() == [851, 627, 625, 1065]
    assert kept_counts.tolist() == [954, 1169, 856, 1045]
    assert int(counts[1] - old_counts[1]) == 682
    assignments = pd.read_csv(ROOT / 'gate_cell_assignments.tsv.gz', sep='\t').set_index('cell_id').reindex(df.cell_id)
    for key, mask in [('phase1_retained', old), ('experimental_depletion_retained', candidate)]:
        saved = assignments[key].astype(str).str.lower().isin(['true', '1', '1.0']).to_numpy()
        assert np.array_equal(saved, mask), key
    return old, candidate, {
        'source': str(ROOT), 'entry_n': len(df), 'state_counts': dict(zip(LABELS, counts.tolist())),
        'old_retained_n': int(old.sum()), 'old_retained_state_counts': old_counts.tolist(),
        'old_dual_high_lost_n': 682, 'candidate_retained_n': int(candidate.sum()),
        'candidate_removed_n': int((~candidate).sum()),
        'candidate_retained_state_counts': kept_counts.tolist(),
        'exact_saved_assignments_reproduced': True, 'scores_or_gates_refit': False,
        'candidate_rule': summary['rules'][str(summary['experimental_depletion']['gate_id'])],
        'source_state_display_mapping': dict(zip(STATES, LABELS)),
    }


def panel_title(ax, letter, title):
    ax.set_title(title, loc='left', pad=13)
    ax.text(-.13, 1.055, letter, transform=ax.transAxes, fontsize=17, weight='bold')


def score_axes(ax, df, params, *, labels=True):
    ax.set_xlim(-.87, max(1.9, float(df.sst_score.max()) + .09))
    ax.set_ylim(-.82, float(df.pv_score.max()) + .14)
    ax.axvline(params['sst_threshold'], color=MUTED, lw=.9, ls=(0, (4, 4)), zorder=0)
    ax.axhline(params['pv_threshold'], color=MUTED, lw=.9, ls=(0, (4, 4)), zorder=0)
    ax.set_xticks([-.5, 0, .5, 1, 1.5])
    ax.set_yticks([-.5, 0, .5, 1, 1.5, 2])
    if labels:
        ax.set_xlabel('SST-associated program score')
        ax.set_ylabel('PV-associated program score')


def figure2(df, out, params):
    fig = plt.figure(figsize=(15, 9.0))
    fig.text(.055, .965, 'DIV90 cortical interneurons already contain distinct\nPV- and SST-associated programs',
             fontsize=23, weight='bold', va='top', linespacing=1.15)
    fig.text(.055, .85, '4,768 cortical LHX6+/ERBB4+ cells  ·  Independent continuous developmental axes',
             fontsize=14, color=MUTED)
    ax = fig.add_axes([.08, .175, .49, .60])
    for state, color in zip(STATES, PALETTE):
        sub = df.loc[df.state == state]
        ax.scatter(sub.sst_score, sub.pv_score, s=13, c=color, alpha=.57, linewidths=0)
    score_axes(ax, df, params)
    panel_title(ax, 'A', 'Developmental state space')
    for x, y, label, count, pct, color, ha, va in [
        (.025, .955, 'PV-biased', '1,076', '22.6%', PALETTE[0], 'left', 'top'),
        (.965, .955, 'Dual-high', '1,309', '27.5%', PALETTE[1], 'right', 'top'),
        (.965, .012, 'SST-biased', '1,075', '22.5%', PALETTE[2], 'right', 'bottom'),
        (.025, .012, 'Unresolved', '1,308', '27.4%', INK, 'left', 'bottom'),
    ]:
        ax.text(x, y, f'{label}\n{pct} · {count} cells', transform=ax.transAxes,
                ha=ha, va=va, fontsize=13, weight='bold', color=color, linespacing=1.5,
                bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': .92, 'pad': 4})
    cmap = LinearSegmentedColormap.from_list('pi_expression', ['#DDE3E8', '#7D85AB', '#4C3478'])
    for rect, letter, gene in [([.68, .55, .26, .225], 'B', 'SST'),
                                ([.68, .175, .26, .225], 'C', 'MEF2C')]:
        axg = fig.add_axes(rect)
        ordered = df.sort_values(gene)
        vmax = float(ordered[gene].max())
        points = axg.scatter(ordered.sst_score, ordered.pv_score, c=ordered[gene], s=5,
                             cmap=cmap, vmin=0, vmax=vmax, alpha=.75, linewidths=0)
        score_axes(axg, df, params, labels=False)
        axg.set_xlabel('SST program score', fontsize=12)
        axg.set_ylabel('PV program score', fontsize=12)
        axg.set_xticks([0, 1]); axg.set_yticks([0, 1, 2])
        axg.tick_params(labelsize=11)
        panel_title(axg, letter, f'{gene} RNA')
        cax = fig.add_axes([.68, rect[1] - .065, .26, .013])
        cb = fig.colorbar(points, cax=cax, orientation='horizontal')
        cb.set_ticks([0, vmax]); cb.set_ticklabels(['Undetected', 'Higher RNA'])
        cb.ax.tick_params(length=0, labelsize=11, pad=3)
        cb.outline.set_visible(False)
    fig.text(.055, .035, '27.5% of cells are dual-high, but coordinated hybrid identity is not yet established.',
             fontsize=14, weight='bold')
    return save(fig, out, 'main_figures/Figure2_PV_SST_state_space',
                'Frozen independent SST (horizontal) and PV (vertical) scores for all 4,768 entry cells. '
                'Dashed lines show the primary operational thresholds; the four quadrants are PV-biased '
                '(1,076; 22.6%), dual-high (1,309; 27.5%), SST-biased (1,075; 22.5%) and unresolved '
                '(1,308; 27.4%). Panels B and C show the same scores colored by existing SST and MEF2C RNA; '
                'MEF2C is a score component, so its relationship to the PV axis is not independent validation. '
                'Color endpoints span zero to the observed maximum, without a display transformation beyond '
                'the frozen normalized RNA values. Overlapping programs are operational and do not establish a hybrid fate.')


def composition_block(ax, values, label, n):
    ax.set_xlim(0, 1); ax.set_ylim(0, 100); ax.axis('off')
    # Draw top to bottom in the common PV / dual-high / SST / unresolved order.
    y = 100
    for name, value, color in zip(LABELS, values, PALETTE):
        y -= value
        ax.add_patch(Rectangle((.08, y), .84, value, facecolor=color, edgecolor='white', lw=2))
        ax.text(.5, y + value / 2, f'{name}\n{value:.1f}%', ha='center', va='center',
                color='white' if name != 'Unresolved' else INK, fontsize=15, weight='bold', linespacing=1.5)
    ax.text(.5, 1.035, f'{n:,} cells', transform=ax.transAxes, ha='center', fontsize=14, color=MUTED)


def figure3(df, out, old):
    fig = plt.figure(figsize=(14, 7.6))
    fig.text(.055, .96, 'The first FAT3/PTPRM gate discards too much\nof the overlapping developmental population',
             fontsize=23, weight='bold', va='top', linespacing=1.15)
    fig.text(.055, .835, 'FAT3-low / PTPRM-undetected retention rule', fontsize=15, color=MUTED)
    a = fig.add_axes([.07, .205, .205, .48]); b = fig.add_axes([.335, .205, .205, .48])
    composition_block(a, [22.6, 27.5, 22.5, 27.4], 'Starting population', 4768)
    composition_block(b, [26.9, 19.8, 19.7, 33.6], 'Retained population', int(old.sum()))
    for x, letter, title in [(.065, 'A', 'Starting population'), (.33, 'B', 'Retained population')]:
        fig.text(x, .765, letter, fontsize=17, weight='bold')
        fig.text(x + .024, .765, title, fontsize=15, weight='bold')
    c = fig.add_axes([.685, .43, .27, .255])
    c.barh([1, 0], [100, 100], color='#EDF0F2', height=.40)
    c.barh([1, 0], [79.1, 47.9], color=PALETTE[:2], height=.40)
    for y, value in [(1, 79.1), (0, 47.9)]:
        c.text(value - 2, y, f'{value:.1f}%', ha='right', va='center', fontsize=17, weight='bold', color='white')
    c.set_yticks([1, 0], ['PV-biased', 'Dual-high'])
    c.set_xticks([0, 50, 100]); c.set_xlim(0, 100); c.set_ylim(-.6, 1.6)
    c.set_xlabel('Recovery in retained fraction (%)', fontsize=13)
    c.spines['left'].set_visible(False); c.tick_params(axis='y', length=0)
    fig.text(.60, .765, 'C', fontsize=17, weight='bold')
    fig.text(.627, .765, 'Recovery of developmental states', fontsize=15, weight='bold')
    fig.text(.795, .30, '682', ha='center', fontsize=43, weight='bold', color=PALETTE[1])
    fig.text(.795, .245, 'dual-high cells lost', ha='center', fontsize=17, weight='bold')
    fig.text(.055, .075, 'Too restrictive for a preservation-oriented strategy.', fontsize=18, weight='bold')
    return save(fig, out, 'main_figures/Figure3_old_gate_tradeoff',
                'The frozen FAT3-low/PTPRM-undetected RNA retention rule captures 3,168 of 4,768 cells. '
                'Starting and retained four-state compositions are shown directly. The rule retains 851/1,076 '
                'PV-biased cells (79.1%) and 627/1,309 dual-high cells (47.9%), losing 682 dual-high cells. '
                'It enriches a PV-biased corner but is too restrictive for a preservation-oriented strategy. '
                'These are operational RNA states, not demonstrated cell fates.')


def box(ax, x, y, w, h, text, color='#F1F4F6', edge='none', fontsize=16, weight='normal', textcolor=INK):
    p = FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.009,rounding_size=.012',
                      facecolor=color, edgecolor=edge, linewidth=1, transform=ax.transAxes)
    ax.add_patch(p)
    ax.text(x + w / 2, y + h / 2, text, ha='center', va='center', transform=ax.transAxes,
            fontsize=fontsize, weight=weight, color=textcolor, linespacing=1.35)


def arrow(ax, start, end, color=MUTED, width=1.5, style='-|>'):
    ax.add_patch(FancyArrowPatch(start, end, transform=ax.transAxes, arrowstyle=style,
                                mutation_scale=13, color=color, linewidth=width,
                                shrinkA=1, shrinkB=1))


def figure4(df, out, candidate):
    fig = plt.figure(figsize=(15, 11.8)); ax = fig.add_axes([0, 0, 1, 1]); ax.axis('off')
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    fig.text(.045, .968, 'A prospective depletion strategy preserves ~89%\nof PV-associated and dual-high cells',
             fontsize=23, weight='bold', va='top', linespacing=1.15)
    ax.text(.39, .872, 'DIV90 viable cells', ha='center', fontsize=17, weight='bold')
    arrow(ax, (.39, .859), (.39, .836))
    ax.text(.39, .818, 'ERBB4+', ha='center', fontsize=24, weight='bold', color=PALETTE[0])
    arrow(ax, (.39, .808), (.39, .784))
    ax.text(.39, .763, 'Cortical interneuron-enriched population', ha='center', fontsize=16)
    ax.text(.39, .742, 'Reference RNA population: 4,768 cells', ha='center', fontsize=13, color=MUTED)
    arrow(ax, (.39, .729), (.39, .705))
    box(ax, .275, .656, .23, .04, 'FGFR2 / PTPRS', fontsize=19, weight='bold')
    arrow(ax, (.335, .648), (.208, .62))
    arrow(ax, (.445, .648), (.57, .62))
    box(ax, .06, .554, .30, .062, 'FGFR2-high OR PTPRS-high\nREMOVE · 744 cells',
        color='#F8EDE6', fontsize=15, weight='bold')
    box(ax, .42, .554, .30, .062, 'BOTH LOW · retain complement\nRETAIN · 4,024 cells',
        color='#E6F1ED', fontsize=15, weight='bold')
    for x, names, vals, colors in [
        (.085, ['SST-biased', 'Dual-high', 'PV-biased', 'Unresolved'],
         [29.4, 18.8, 16.4, 35.3], [PALETTE[2], PALETTE[1], PALETTE[0], PALETTE[3]]),
        (.445, LABELS, [23.7, 29.1, 21.3, 26.0], PALETTE),
    ]:
        for j, (label, val, color) in enumerate(zip(names, vals, colors)):
            y = .523 - j * .031
            ax.scatter([x], [y], c=[color], s=90, transform=ax.transAxes)
            ax.text(x + .016, y, f'{val:.1f}% {label.lower() if label == "Dual-high" else label}',
                    transform=ax.transAxes, va='center', fontsize=15)
    box(ax, .756, .425, .215, .187, '', color='#F0F5F8')
    for y, value, label in [(.565, '88.7%', 'PV-biased recovery'),
                            (.505, '89.3%', 'dual-high recovery'),
                            (.445, '84.4%', 'total cell recovery')]:
        ax.text(.778, y + .014, value, fontsize=27, weight='bold', color=PALETTE[0])
        ax.text(.778, y - .010, label, fontsize=13, color=MUTED)
    ax.plot([.045, .965], [.395, .395], color=LINE, linewidth=1, transform=ax.transAxes)
    ax.text(.045, .365, 'PROSPECTIVE VALIDATION', fontsize=17, weight='bold')
    box(ax, .08, .288, .235, .047, 'Retain fraction', color='#E6F1ED', fontsize=17, weight='bold')
    box(ax, .3825, .288, .235, .047, 'Removed fraction', color='#F8EDE6', fontsize=17, weight='bold')
    box(ax, .685, .288, .235, .047, 'Unsorted ERBB4+ control', fontsize=15, weight='bold')
    ax.text(.35, .31, 'vs', ha='center', va='center', color=MUTED, fontsize=14)
    ax.text(.653, .31, 'vs', ha='center', va='center', color=MUTED, fontsize=14)
    for x in [.1975, .5, .8025]:
        arrow(ax, (x, .276), (x, .239))
    ax.plot([.1975, .8025], [.239, .239], color=MUTED, linewidth=1.5, transform=ax.transAxes)
    arrow(ax, (.5, .239), (.5, .21))
    ax.text(.5, .187, 'Identical maturation', fontsize=21, weight='bold', ha='center')
    arrow(ax, (.5, .178), (.5, .152))
    ax.text(.5, .124, 'PV protein  ·  SST protein  ·  PVALB RNA / mature transcriptomic identity',
            ha='center', fontsize=16)
    ax.text(.5, .094, 'Fast-spiking physiology', ha='center', fontsize=16)
    ax.text(.5, .035, 'This gate is a testable enrichment hypothesis, not yet a validated PV sort.',
            ha='center', fontsize=17, weight='bold')
    return save(fig, out, 'main_figures/Figure4_prospective_sort_strategy',
                'Proposed ERBB4 entry followed by removal of FGFR2-high OR PTPRS-high cells; retain the '
                'complement (both below their frozen RNA cutoffs). Among 4,768 reference entry cells, '
                '4,024 are retained and 744 removed. Retained recoveries are 954/1,076 PV-biased (88.7%), '
                '1,169/1,309 dual-high (89.3%) and 4,024/4,768 total (84.4%). Branch percentages are '
                'within-fraction composition, rounded independently. Starting SST-biased composition is '
                '22.5% versus 21.3% retained; the principal value is preservation, not dramatic SST purification. '
                'Collect retained, removed and unsorted ERBB4+ control populations, mature identically, and compare '
                'later PV/SST protein, PVALB RNA / transcriptomic identity and fast-spiking physiology. '
                'An RNA-derived sorting hypothesis does not establish protein gate behavior; prospective '
                'protein-level validation and ERBB4 entry specificity must be tested.')


def pi_summary(out):
    fig = plt.figure(figsize=(15, 12.6)); ax = fig.add_axes([0, 0, 1, 1]); ax.axis('off')
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.text(.045, .94, 'DIV90: a prospective route toward PV-associated precursors',
            fontsize=23, weight='bold', va='top')
    statements = [
        ('ERBB4 identifies the cortical interneuron population of interest.', .855),
        ('PV-associated developmental information is present at DIV90 even though\nPVALB RNA is essentially absent.', .75),
        ('FGFR2/PTPRS depletion preserves ~89% of PV-biased and dual-high cells\nand provides a prospective sort-and-mature experiment.', .625),
    ]
    for i, (statement, y) in enumerate(statements, 1):
        ax.text(.046, y, str(i), fontsize=28, weight='bold', color=PALETTE[0], va='center')
        ax.text(.096, y, statement, fontsize=19, weight='bold', va='center', linespacing=1.55)
    ax.plot([.045, .955], [.553, .553], color=LINE, transform=ax.transAxes)
    box(ax, .39, .472, .22, .047, 'ERBB4+', color='#EAF0F7', fontsize=22, weight='bold')
    arrow(ax, (.5, .461), (.5, .428))
    box(ax, .37, .363, .26, .047, 'FGFR2 / PTPRS', fontsize=20, weight='bold')
    arrow(ax, (.445, .351), (.345, .314))
    arrow(ax, (.555, .351), (.655, .314))
    box(ax, .265, .252, .16, .045, 'remove', color='#F8EDE6', fontsize=20, weight='bold')
    box(ax, .575, .252, .16, .045, 'retain', color='#E6F1ED', fontsize=20, weight='bold')
    arrow(ax, (.655, .24), (.655, .205))
    ax.text(.655, .175, 'mature', ha='center', fontsize=21, weight='bold')
    arrow(ax, (.655, .158), (.655, .122))
    ax.text(.655, .093, 'test PV fate', ha='center', fontsize=21, weight='bold', color=PALETTE[0])
    ax.text(.5, .035, 'Computational prediction → prospective protein-level validation required.',
            ha='center', fontsize=18, weight='bold')
    return save(fig, out, 'PI_summary/PI_summary_one_page',
                'The three requested PI statements appear above the proposed ERBB4 → FGFR2/PTPRS → retained-fraction '
                'maturation schematic. The full experiment collects both fractions and an unsorted ERBB4+ '
                'control as shown in main Figure 4. Computational prediction requires prospective protein-level validation.',
                extensions=('pdf', 'png'))


def run(cells=None, out=OUT, summary=None, params=None):
    out = Path(out)
    if not out.is_relative_to(ROOT) or out == ROOT:
        raise ValueError('All figure runtime outputs must be in a new directory under the frozen phase-two Turbo root.')
    out.mkdir(parents=True, exist_ok=True)
    df = cells.copy() if isinstance(cells, pd.DataFrame) else pd.read_csv(cells or ROOT / 'cells.tsv.gz', sep='\t')
    summary = summary if summary is not None else json.loads((ROOT / 'gate_summary.json').read_text())
    params = params if params is not None else json.loads((ROOT / 'state_parameters.json').read_text())
    configure()
    old, candidate, audit = frozen_audit(df, summary, params)
    manifest = [figure2(df, out, params), figure3(df, out, old),
                figure4(df, out, candidate), pi_summary(out)]
    (out / 'provenance').mkdir(exist_ok=True)
    (out / 'provenance/main_figures_2_4_audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    (out / 'provenance/main_figures_2_4_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    return {'figures': manifest, 'audit': audit}


if __name__ == '__main__':
    result = run()
    print(json.dumps({'figures': [item['figure'] for item in result['figures']],
                      'old_retained_n': result['audit']['old_retained_n'],
                      'candidate_retained_n': result['audit']['candidate_retained_n']}, indent=2))
