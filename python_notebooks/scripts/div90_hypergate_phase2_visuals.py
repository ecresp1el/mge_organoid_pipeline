#!/usr/bin/env python3
"""Publication figures and offline exploration of the frozen DIV90 phase-two cells.

Only code belongs in the checkout. All generated figures/HTML/metadata are written
under the new Turbo phase-two output directory; phase one is never modified.
"""
from __future__ import annotations

import argparse
import html
import json
import math
import os
import sys
from pathlib import Path

OUT = Path("/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_hypergate_sst_pv_phase2")
os.environ["MPLCONFIGDIR"] = str(OUT / "cache/matplotlib")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

STATES = ["PV-biased", "PV/SST hybrid", "SST-biased", "unresolved/immature"]
SHORT = ["PV-biased", "Hybrid", "SST-biased", "Unresolved"]
COLORS = dict(zip(STATES, ["#326DAB", "#39876F", "#CB733F", "#A6ACB6"]))
SURFACE_COLORS = ["SST", "ERBB4", "FAT3", "PTPRM", "CXCR4", "ACKR3", "NRP2"]


def configure():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
        "axes.titlesize": 11, "axes.labelsize": 10, "axes.spines.top": False,
        "axes.spines.right": False, "axes.linewidth": .7, "xtick.major.width": .6,
        "ytick.major.width": .6, "pdf.fonttype": 42, "ps.fonttype": 42,
        "svg.fonttype": "none", "savefig.facecolor": "white", "figure.facecolor": "white"})


def save(fig, out, stem, caption):
    directory = out / "figures"
    directory.mkdir(parents=True, exist_ok=True)
    # Matplotlib rasterizes dense colorbar meshes by default; retain true vector gradients.
    for artist in fig.findobj():
        if hasattr(artist, "get_rasterized") and artist.get_rasterized():
            artist.set_rasterized(False)
    for extension in ("pdf", "svg", "png"):
        fig.savefig(directory / f"{stem}.{extension}", dpi=600 if extension == "png" else None,
                    bbox_inches="tight", pad_inches=.12)
    plt.close(fig)
    return {"figure": stem, "caption": caption}


def label(ax, letter, title):
    ax.set_title(title, loc="left", pad=10)
    ax.text(-.11, 1.045, letter, transform=ax.transAxes, fontsize=12, fontweight="bold")


def legend_states(fig, y=.015):
    fig.legend([Line2D([0], [0], marker="o", linestyle="", color=COLORS[s], markersize=6)
                for s in STATES], SHORT, frameon=False, loc="lower center", ncol=4,
               bbox_to_anchor=(.5, y), fontsize=9, columnspacing=1.7)


def limits(df):
    result = []
    for key in ("sst_score", "pv_score"):
        lo, hi = np.nanmin(df[key]), np.nanmax(df[key])
        result.append((lo - .04 * (hi - lo), hi + .04 * (hi - lo)))
    return result


def thresholds(df):
    parameter_path = OUT / "state_parameters.json"
    if parameter_path.exists():
        parameters = json.loads(parameter_path.read_text())
        if "sst_threshold" in parameters and "pv_threshold" in parameters:
            return [float(parameters["sst_threshold"]), float(parameters["pv_threshold"])]
    # Bounds between operational low/high groups preserve the actual assignments.
    cuts = []
    for score, high_states in [("sst_score", [STATES[1], STATES[2]]),
                               ("pv_score", [STATES[0], STATES[1]])]:
        hi = df.loc[df.state.isin(high_states), score]
        lo = df.loc[~df.state.isin(high_states), score]
        cuts.append(float((hi.min() + lo.max()) / 2) if len(hi) and len(lo) else float(df[score].median()))
    return cuts


def landscape(ax, df, lim, cuts=None, color="state", point_size=4, alpha=.62, title=None):
    if color == "state":
        for state in STATES:
            sub = df.loc[df.state == state]
            ax.scatter(sub.sst_score, sub.pv_score, s=point_size, color=COLORS[state],
                       alpha=alpha, linewidths=0)
    elif color in df and pd.api.types.is_numeric_dtype(df[color]):
        ordered = df.sort_values(color)
        values = ordered[color].to_numpy()
        vmax = max(float(np.nanquantile(values, .99)), 1e-8)
        dots = ax.scatter(ordered.sst_score, ordered.pv_score, c=values, s=point_size,
                          cmap="magma_r", vmin=float(np.nanmin(values)), vmax=vmax,
                          linewidths=0, alpha=alpha)
        plt.colorbar(dots, ax=ax, fraction=.035, pad=.02, label="log1p(CP10K)")
    elif color in df:
        for j, (name, sub) in enumerate(df.groupby(color, observed=True)):
            ax.scatter(sub.sst_score, sub.pv_score, s=point_size, c=[plt.cm.tab20(j % 20)],
                       linewidths=0, alpha=alpha, label=str(name))
        ax.legend(frameon=False, fontsize=6, ncol=2, loc="best", markerscale=2)
    if cuts is not None:
        ax.axvline(cuts[0], color="#56616D", ls=(0, (3, 3)), lw=.7, zorder=0)
        ax.axhline(cuts[1], color="#56616D", ls=(0, (3, 3)), lw=.7, zorder=0)
    ax.set(xlim=lim[0], ylim=lim[1], xlabel="SST developmental program", ylabel="PV developmental program")
    if title:
        ax.set_title(title, loc="left")


def apply_rule(df, rule):
    pieces = []
    for item in rule["rules"]:
        values = df[item["gene"]].to_numpy()
        cut = float(item["threshold"])
        op = item["op"]
        pieces.append({"<=": lambda: values <= cut, "<": lambda: values < cut,
                       ">=": lambda: values >= cut, ">": lambda: values > cut,
                       "==": lambda: values == cut, "=": lambda: values == cut}[op]())
    hit = np.logical_or.reduce(pieces) if rule.get("logic") == "OR" else np.logical_and.reduce(pieces)
    return ~hit if rule.get("action", "remove") == "remove" else hit


def masks(df, out, summary):
    result = {}
    assignment_path = out / "gate_cell_assignments.tsv.gz"
    if assignment_path.exists():
        assignment = pd.read_csv(assignment_path, sep="\t").set_index("cell_id").reindex(df.cell_id)
        for key in ["phase1", "best_single", "best_pair", "selected", "erbb4_cxcr4"]:
            col = f"{key}_retained"
            if col in assignment:
                result[key] = assignment[col].astype(str).str.lower().isin(["true", "1", "1.0"]).to_numpy()
    for key in ["best_single", "best_pair", "selected", "phase1_reference", "erbb4_cxcr4"]:
        name = "phase1" if key == "phase1_reference" else key
        if summary.get(key):
            row = summary[key]
            result[name] = apply_rule(df, summary["rules"][str(row["gate_id"])])
    if "phase1" not in result:
        result["phase1"] = (df.FAT3 <= .764341).to_numpy() & (df.PTPRM <= 0).to_numpy()
    if "selected" not in result:
        raise ValueError("Selected depletion rule/assignments are required for honest sorting figures.")
    return result


def composition(df, mask):
    counts = df.loc[mask, "state"].value_counts().reindex(STATES, fill_value=0).to_numpy()
    return counts / max(counts.sum(), 1), counts


def composition_bars(ax, df, groups, annotate=True):
    positions = np.arange(len(groups))
    base = np.zeros(len(groups))
    fractions = np.array([composition(df, group[1])[0] for group in groups])
    for j, state in enumerate(STATES):
        height = fractions[:, j] * 100
        ax.bar(positions, height, bottom=base, width=.6, color=COLORS[state], edgecolor="white", linewidth=.6)
        if annotate:
            for x, bottom, value in zip(positions, base, height):
                if value > 7:
                    ax.text(x, bottom + value / 2, f"{value:.1f}%", ha="center", va="center", fontsize=9,
                            color="white" if state != STATES[3] else "#24323E")
        base += height
    ax.set(xticks=positions, xticklabels=[f"{g[0]}\nn={int(np.sum(g[1])):,}" for g in groups],
           ylim=(0, 100), ylabel="Fraction of cells (%)")


def figure1(df, out, lim, cuts):
    fig = plt.figure(figsize=(9, 8))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 4.6], width_ratios=[4.6, 1],
                          hspace=.06, wspace=.06, left=.1, bottom=.14, right=.95, top=.91)
    main, top, side = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 1])
    xx, yy = np.meshgrid(np.linspace(*lim[0], 100), np.linspace(*lim[1], 100))
    density = gaussian_kde(np.vstack([df.sst_score, df.pv_score]))
    zz = density(np.vstack([xx.ravel(), yy.ravel()])).reshape(xx.shape)
    main.contour(xx, yy, zz, levels=np.linspace(zz.max() * .08, zz.max() * .88, 7),
                 colors="#D0D5DC", linewidths=.85, zorder=0)
    landscape(main, df, lim, cuts, point_size=5, alpha=.52)
    for state in STATES:
        sub = df.loc[df.state == state]
        top.hist(sub.sst_score, bins=50, histtype="step", color=COLORS[state], lw=1,
                 range=lim[0], density=False)
        side.hist(sub.pv_score, bins=50, histtype="step", color=COLORS[state], lw=1,
                  range=lim[1], orientation="horizontal", density=False)
    top.set(xlim=lim[0], ylabel="Cells"); top.tick_params(labelbottom=False)
    side.set(ylim=lim[1], xlabel="Cells"); side.tick_params(labelleft=False)
    names = [(STATES[0], .02, .98, "left", "top"), (STATES[1], .98, .98, "right", "top"),
             (STATES[2], .98, .02, "right", "bottom"), (STATES[3], .02, .02, "left", "bottom")]
    for state, x, y, ha, va in names:
        count = int((df.state == state).sum())
        main.text(x, y, f"{state}\n{count:,} cells ({count / len(df):.1%})", transform=main.transAxes,
                  ha=ha, va=va, color=COLORS[state], fontsize=9,
                  bbox=dict(facecolor="white", edgecolor="none", alpha=.86, pad=3))
    fig.suptitle("DIV90: independent PV and SST developmental programs", x=.1, ha="left", fontsize=15)
    fig.text(.1, .055, f"{len(df):,} frozen LHX6+/ERBB4+ cortical interneurons · dashed lines: operational pooled thresholds\n"
             "Density contours show the continuous distribution; labels describe score regions, not fates.", fontsize=9, color="#53606B")
    return save(fig, out, "figure_01_pv_sst_landscape", "Independent score axes, density contours, individual cells, marginal histograms and four operational score regions. All frozen entry cells are shown.")


def figure2(df, out, lim, cuts):
    groups = [("Cell line", key, sub) for key, sub in df.groupby("cell_line", observed=True)]
    groups += [("Condition", key, sub) for key, sub in df.groupby("condition", observed=True)]
    groups += [("Line × condition", f"{key[0]} / {key[1]}", sub) for key, sub in df.groupby(["cell_line", "condition"], observed=True)]
    columns = min(4, len(groups)); rows = math.ceil(len(groups) / columns)
    fig, axes = plt.subplots(rows, columns, figsize=(3.35 * columns, 2.9 * rows), squeeze=False)
    for ax, (kind, name, sub) in zip(axes.ravel(), groups):
        ax.scatter(df.sst_score, df.pv_score, s=1.5, color="#E8E9EC", linewidths=0)
        landscape(ax, sub, lim, cuts, point_size=5, alpha=.72)
        title = f"{kind}: {name}".replace(" / higher", "\nhigher").replace(" / lower", "\nlower")
        ax.set_title(f"{title}\nn={len(sub):,}", loc="left", fontsize=9)
    for ax in axes.ravel()[len(groups):]: ax.axis("off")
    fig.suptitle("Cell line and culture condition on the same developmental axes\nCV: higher glucose (~2×) · MW: lower glucose (~1×); coupled culture factors", fontsize=13, y=1.01)
    fig.tight_layout(rect=(0, .06, 1, .98)); legend_states(fig)
    return save(fig, out, "figure_02_sample_facets", "Same score scale across cell-line and condition facets. Light gray background shows all cells. Sample dependence is preserved; pooled thresholds are unchanged.")


def module_list(out):
    modules = pd.read_csv(out / "module_genes.tsv", sep="\t")
    if "included" in modules:
        modules = modules.loc[modules.included.astype(str).str.lower().isin(["true", "1", "1.0", "yes"])]
    return modules


def figure3(df, out):
    modules = module_list(out)
    modules = modules.loc[modules.gene.isin(df.columns)].copy()
    genes = modules.gene.tolist()
    means = df.groupby("state", observed=True)[genes].mean().reindex(STATES)
    detection = df.assign(**{g: df[g] > 0 for g in genes}).groupby("state", observed=True)[genes].mean().reindex(STATES)
    zmeans = (means - means.mean(axis=0)) / means.std(axis=0, ddof=0).replace(0, 1)
    fig = plt.figure(figsize=(12, 9))
    gs = fig.add_gridspec(3, 3, height_ratios=[1.2, 1, 1], hspace=.7, wspace=.42)
    ax = fig.add_subplot(gs[0, :])
    x, y = np.meshgrid(np.arange(len(genes)), np.arange(4))
    dots = ax.scatter(x.ravel(), y.ravel(), s=12 + detection.to_numpy().ravel() * 220,
                      c=zmeans.to_numpy().ravel(), cmap="RdBu_r", vmin=-2, vmax=2, edgecolor="#52606B", lw=.3)
    ax.set(xticks=np.arange(len(genes)), xticklabels=genes, yticks=np.arange(4), yticklabels=SHORT,
           xlim=(-.6, len(genes) - .4), ylim=(3.5, -.5))
    ax.tick_params(axis="x", rotation=30)
    for tick, program in zip(ax.get_xticklabels(), modules.program):
        tick.set_color(COLORS[STATES[0]] if "PV" in str(program).upper() else COLORS[STATES[2]])
    plt.colorbar(dots, ax=ax, fraction=.02, pad=.02, label="Gene-wise mean z score")
    size_handles = [ax.scatter([], [], s=12 + value * 220, facecolors="none", edgecolors="#52606B") for value in [.1, .5, 1]]
    ax.legend(size_handles, ["10%", "50%", "100%"], title="Detected", frameon=False, loc="upper left", bbox_to_anchor=(1.11, 1), fontsize=8)
    label(ax, "A", "Every module component: mean expression and fraction detected")
    metrics = [("n_genes", "Detected genes"), ("total_counts", "Total RNA counts"),
               ("mito_fraction", "Mitochondrial fraction"), ("doublet_score", "Scrublet doublet score"),
               ("stress_score", "Stress program"), ("cell_cycle_score", "Cell-cycle program")]
    for i, (metric, title) in enumerate(metrics):
        ax = fig.add_subplot(gs[1 + i // 3, i % 3])
        if metric in df and pd.to_numeric(df[metric], errors="coerce").notna().any():
            data = [df.loc[df.state == state, metric].dropna().to_numpy() for state in STATES]
            violins = ax.violinplot(data, showmeans=False, showmedians=True, showextrema=False)
            for body, state in zip(violins["bodies"], STATES):
                body.set_facecolor(COLORS[state]); body.set_alpha(.7); body.set_edgecolor("none")
            violins["cmedians"].set_color("#172A39")
            ax.set_xticks(range(1, 5), ["PV", "Hybrid", "SST", "Unres."])
            if metric == "total_counts": ax.set_yscale("log")
        else:
            ax.text(.5, .5, "Not available in source data\nDoublets cannot be ruled out", ha="center", va="center",
                    transform=ax.transAxes, color="#6B737B", fontsize=9)
            ax.set(xticks=[], yticks=[])
        label(ax, chr(66 + i), title)
    hybrid = df.loc[df.state == STATES[1]]
    pv_genes = modules.loc[modules.program.str.lower() == "pv", "gene"].tolist()
    sst_genes = modules.loc[modules.program.str.lower() == "sst", "gene"].tolist()
    support = ((hybrid[pv_genes] > 0).sum(axis=1) >= 2) & ((hybrid[sst_genes] > 0).sum(axis=1) >= 2)
    non_sst = [gene for gene in sst_genes if gene != "SST"]
    stringent = ((hybrid[pv_genes] > 0).sum(axis=1) >= 2) & ((hybrid[non_sst] > 0).sum(axis=1) >= 2)
    fig.suptitle("Hybrid-state checks: module components and available QC", fontsize=14, y=.995)
    fig.text(.085, .022, f"Within {len(hybrid):,} operational hybrid cells: {support.sum():,} ({support.mean():.1%}) detect ≥2 genes from each program; "
             f"{stringent.sum():,} ({stringent.mean():.1%}) meet this after excluding SST.\n"
             "Dual-high score assignment alone does not establish coordinated co-expression or validate a biological hybrid identity.",
             fontsize=8, color="#52606B")
    fig.subplots_adjust(left=.085, right=.84, top=.93, bottom=.10)
    return save(fig, out, "figure_03_hybrid_validation", "Gene-level module dot plot across all four operational states; point size is detection fraction and color is gene-wise standardized state mean. QC distributions are descriptive. Available Scrublet doublet scores are computational estimates fitted by sample with an assumed 5% expected doublet rate, not direct doublet measurements. These checks do not establish developmental fate.")


def figure4(df, out):
    ranking = pd.read_csv(out / "tables/surface_state_markers.tsv", sep="\t")
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 6), sharex=False)
    for ax, state, short in zip(axes, STATES[:3], SHORT[:3]):
        ranked = ranking.loc[ranking.state == state].sort_values("cohens_d", ascending=False).head(10).iloc[::-1]
        y = np.arange(len(ranked))
        ax.barh(y, ranked.cohens_d, color=COLORS[state], alpha=.75, height=.6)
        ax.set(yticks=y, yticklabels=ranked.gene, xlabel="Standardized mean difference (Cohen's d)")
        right = ax.twiny()
        right.scatter(ranked.detection_state * 100, y, s=35, facecolors="white", edgecolors="#172A39", lw=.9)
        right.set(xlim=(0, 105), xlabel="Detected in state (%)")
        ax.set_title(short, loc="left", pad=38, color=COLORS[state], fontsize=12)
        for pos, row in enumerate(ranked.itertuples()):
            ax.text(.98, (pos + .5) / max(len(ranked), 1), f"Δ {row.mean_difference:+.2f}",
                    transform=ax.transAxes, ha="right", va="center", fontsize=7, color="#344351")
    fig.suptitle("Surfaceome discovery by operational state", x=.07, ha="left", fontsize=15, y=1.01)
    fig.text(.07, .02, "Curated CSPA transmembrane/GPI surface annotation · top positive effects per state · Δ = mean log1p(CP10K) difference\n"
             "State-defining genes are excluded. Rankings describe RNA association; live-cell protein detection requires validation.", fontsize=9, color="#52606B")
    fig.tight_layout(rect=(.01, .1, 1, .94))
    return save(fig, out, "figure_04_surface_marker_discovery", "Top surface-associated positive effects separately for PV-biased, hybrid and SST-biased cells. Bars: Cohen's d; open circles: within-state detection; Δ: mean RNA expression difference. Full protein annotations and ranking statistics are in surface_state_markers.tsv.")


def score_distributions(ax, df, masks_, score, colors):
    grid = np.linspace(df[score].min(), df[score].max(), 240)
    for (name, mask), color in zip(masks_, colors):
        vals = df.loc[mask, score].to_numpy()
        if len(vals) > 2 and np.std(vals) > 0:
            ax.plot(grid, gaussian_kde(vals)(grid), label=name, color=color, lw=1.4)
    ax.set(xlabel="PV program" if score == "pv_score" else "SST program", ylabel="Density")
    ax.legend(frameon=False, fontsize=8)


def figure5(df, out, lim, cuts, capture):
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    ax = axes[0, 0]
    ax.scatter(df.sst_score, df.pv_score, color="#DCE0E5", s=4, linewidths=0)
    landscape(ax, df.loc[capture], lim, cuts, point_size=6, alpha=.65)
    label(ax, "A", "Captured: FAT3 ≤ 0.764341 and PTPRM = 0")
    composition_bars(axes[0, 1], df, [("Starting", np.ones(len(df), bool)), ("Captured", capture), ("Lost", ~capture)])
    label(axes[0, 1], "B", "Four-state composition")
    for ax, score, letter in zip(axes[1], ["pv_score", "sst_score"], ["C", "D"]):
        score_distributions(ax, df, [("Starting", np.ones(len(df), bool)), ("Captured", capture), ("Lost", ~capture)],
                            score, ["#647282", "#326DAB", "#CB733F"])
        label(ax, letter, "Continuous score distribution")
    fig.suptitle("Reinterpretation of the completed FAT3/PTPRM gate", fontsize=14, y=.995)
    fig.tight_layout(rect=(0, .07, 1, .96)); legend_states(fig)
    return save(fig, out, "figure_05_fat3_ptprm_reinterpretation", "Exact saved phase-one RNA thresholds applied without refitting. Captured cells overlay all starting cells; composition and continuous score densities include the full population rather than only module extremes.")


def figure6(df, out, summary, retained):
    selected = summary["selected"]
    rule = summary["rules"][str(selected["gate_id"])]; rules = rule["rules"]
    genes = list(dict.fromkeys(r["gene"] for r in rules))
    gx = genes[0]; gy = genes[1] if len(genes) > 1 else "ERBB4"
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), gridspec_kw={"width_ratios": [1.2, 1, 1]})
    rng = np.random.default_rng(903)
    x = df[gx].to_numpy(); y = df[gy].to_numpy()
    # Jitter is only a display aid for RNA-zero ties and is never used in gates.
    jx = rng.normal(0, .008 * max(x.max(), 1), len(df)); jy = rng.normal(0, .008 * max(y.max(), 1), len(df))
    for mask, color, name in [(retained, "#39876F", "Retained"), (~retained, "#CB733F", "Removed")]:
        axes[0].scatter(x[mask] + jx[mask], y[mask] + jy[mask], s=4, color=color, alpha=.45, linewidths=0, label=name)
    for item in rules:
        if item["gene"] == gx: axes[0].axvline(item["threshold"], color="#172A39", ls="--", lw=1)
        elif item["gene"] == gy: axes[0].axhline(item["threshold"], color="#172A39", ls="--", lw=1)
    axes[0].set(xlabel=f"{gx} RNA · log1p(CP10K)", ylabel=f"{gy} RNA · log1p(CP10K)")
    axes[0].legend(frameon=False, markerscale=3, fontsize=8)
    label(axes[0], "A", "Virtual RNA gate")
    composition_bars(axes[1], df, [("Removed", ~retained), ("Retained", retained)])
    label(axes[1], "B", "Both experimental fractions")
    recovery = []
    for state in STATES:
        state_mask = (df.state == state).to_numpy()
        recovery.append(100 * (retained & state_mask).sum() / max(state_mask.sum(), 1))
    axes[2].barh(np.arange(4), recovery, color=[COLORS[s] for s in STATES], height=.55)
    axes[2].set(yticks=np.arange(4), yticklabels=SHORT, xlabel="Recovery in retained fraction (%)", xlim=(0, 105))
    axes[2].invert_yaxis()
    for i, value in enumerate(recovery): axes[2].text(min(value + 1, 99), i, f"{value:.1f}", va="center", ha="left", fontsize=8)
    label(axes[2], "C", "Preservation of each state")
    constraint = "Exploratory rule with ≥20% SST-biased removal constraint" if summary.get("visual_depletion_constraint") else "Proposed depletion rule"
    fig.suptitle(f"{constraint}\n{selected['gate_label']}", fontsize=12, y=1.03)
    fig.text(.035, -.005, "Thresholds are transcriptomic, not fluorescence cutoffs. Small plotting jitter separates tied RNA measurements; classification uses exact unjittered values.", fontsize=8.5, color="#52606B")
    fig.tight_layout(rect=(0, .05, 1, .95))
    return save(fig, out, "figure_06_best_sst_depletion", "Virtual RNA FACS display of the selected practical rule and its exact thresholds, removed/retained composition, and per-state recovery. Display jitter does not affect classification. Fluorescence gates must be calibrated experimentally.")


def figure7(df, out, retained):
    fig, ax = plt.subplots(figsize=(7.4, 6))
    composition_bars(ax, df, [("Starting ERBB4+", np.ones(len(df), bool)), ("Removed", ~retained), ("Retained", retained)])
    ax.set_title("Before and after the proposed SST-depletion gate", loc="left", pad=15, fontsize=13)
    ax.spines["bottom"].set_visible(False)
    fig.tight_layout(rect=(0, .09, 1, 1)); legend_states(fig)
    return save(fig, out, "figure_07_before_after_composition", "Four-state fractions in the starting population, selected removed fraction and complementary retained fraction. All four operational states contribute to the denominators.")


def figure8(df, out, retained):
    keys = ["cell_line", "condition"]
    grouped = list(df.groupby(keys, observed=True))
    fig, axes = plt.subplots(1, 4, figsize=(14, max(4.5, .48 * len(grouped) + 1.4)), sharey=True)
    names = [" / ".join(map(str, key)) for key, _ in grouped]
    for j, (ax, state, short) in enumerate(zip(axes, STATES, SHORT)):
        for i, (_, sub) in enumerate(grouped):
            indices = df.index.get_indexer(sub.index)
            before = (sub.state == state).mean() * 100
            kept = sub.loc[retained[indices]]
            after = (kept.state == state).mean() * 100 if len(kept) else np.nan
            ax.plot([before, after], [i, i], color=COLORS[state], lw=1.4, alpha=.8)
            ax.scatter(before, i, facecolor="white", edgecolor=COLORS[state], s=38, lw=1.2, zorder=3)
            ax.scatter(after, i, color=COLORS[state], s=38, zorder=4)
        ax.set(xlim=(0, 100), xlabel="State fraction (%)", title=short, yticks=np.arange(len(names)))
        ax.grid(axis="x", color="#E9ECEF", lw=.6)
    axes[0].set_yticklabels(names); axes[0].invert_yaxis()
    fig.legend([Line2D([0], [0], marker="o", ls="", mfc="white", mec="#344351"),
                Line2D([0], [0], marker="o", ls="", color="#344351")], ["Starting", "Retained"],
               loc="lower center", ncol=2, frameon=False)
    fig.suptitle("Culture-condition biology: the same gate has different starting populations\nCV: higher glucose (~2×) · MW: lower glucose (~1×); descriptive culture contrast", fontsize=12)
    fig.tight_layout(rect=(0, .06, 1, .95))
    return save(fig, out, "figure_08_condition_biology", "Within every observed cell-line/condition combination, open points show starting state fractions and filled points show retained fractions after the same global gate. Differences are descriptive, with cell-level rather than independent-culture replication.")


def figure9(df, out, summary):
    candidates = pd.read_csv(out / "tables/gate_shortlist.tsv", sep="\t")
    pareto = pd.read_csv(out / "tables/gate_pareto.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(9, 6.8))
    size = 8 + 60 * candidates.retained_n / len(df)
    ax.scatter(candidates.target_recovery * 100, candidates.sst_contamination * 100,
               s=size, c="#C7CDD4", alpha=.45, linewidths=0, label="Shortlisted candidate gates")
    frontier = pareto.sort_values("target_recovery")
    ax.plot(frontier.target_recovery * 100, frontier.sst_contamination * 100,
            color="#243E53", lw=1.3, marker="o", markersize=3, label="Pareto frontier")
    styles = [("best_single", "Practical one-marker REMOVE", "#39876F", "o"),
              ("best_pair", "Practical two-marker REMOVE", "#326DAB", "D"),
              ("phase1_reference", "FAT3/PTPRM RETAIN reference", "#9D658C", "s"),
              ("erbb4_cxcr4", "New ERBB4/CXCR4 REMOVE", "#CB733F", "^"),
              ("unconstrained_selected", "Unconstrained REMOVE optimum", "#545A64", "P"),
              ("phase1_erbb4_cxcr4_retention", "Inherited ERBB4/CXCR4 RETAIN", "#91622E", "v")]
    positions = [(48, 24.5), (47, 17.8), (35, 22), (69, 13.5), (66, 27), (5, 25.5)]
    for (key, name, color, marker), position in zip(styles, positions):
        row = summary.get(key)
        if not row: continue
        x, y = 100 * row["target_recovery"], 100 * row["sst_contamination"]
        ax.scatter(x, y, s=110, color=color, marker=marker, edgecolor="white", linewidth=.8, zorder=5, label=name)
        ax.annotate(f"{name}\n{row['gate_label']}", (x, y), xytext=position, textcoords="data", fontsize=7.5,
                    color=color, ha="left", arrowprops=dict(arrowstyle="-", color=color, lw=.7),
                    bbox=dict(facecolor="white", edgecolor="none", alpha=.8, pad=1))
    ax.set(xlabel="PV-biased + hybrid recovery (%)", ylabel="SST-biased contamination among retained cells (%)", xlim=(-2, 106), ylim=(-1, 29))
    ax.set_title("Experimental tradeoffs: preserve both PV-associated score regions", loc="left", pad=14, fontsize=13)
    ax.grid(color="#ECEEF1", lw=.6, zorder=0)
    ax.legend(frameon=False, loc="lower right", fontsize=8)
    fig.text(.13, .015, "Lower contamination and higher target recovery are preferred. Point area scales with retained cell count.\n"
             "All gates are evaluated on the observed dataset; these are exploratory tradeoffs, not held-out predictions.", fontsize=9, color="#52606B")
    fig.tight_layout(rect=(0, .085, 1, 1))
    return save(fig, out, "figure_09_pareto_frontier", "Recovery/contamination tradeoff with point size proportional to total retained cells. Gray points display the shortlist; the frontier uses the evaluated candidate universe. Practical highlighted gates enforce the analyst-selected minimum 20% SST-biased removal constraint. The saved phase-one FAT3/PTPRM and positive ERBB4/CXCR4 RETAIN rules are explicitly distinguished from new REMOVE rules and the unconstrained preservation optimum.")


def tree_box(ax, xy, text, color="#F1F4F7", width=.27, height=.10, fontsize=10, edge="#C6CED7"):
    x, y = xy
    patch = FancyBboxPatch((x - width / 2, y - height / 2), width, height,
                           boxstyle="round,pad=.012,rounding_size=.015", facecolor=color, edgecolor=edge, lw=.8)
    ax.add_patch(patch); ax.text(x, y, text, ha="center", va="center", fontsize=fontsize, color="#203142")


def arrow(ax, start, end, text=None):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=11, lw=1, color="#52606B",
                                connectionstyle="angle3,angleA=-90,angleB=0"))
    if text: ax.text((start[0] + end[0]) / 2, (start[1] + end[1]) / 2, text, fontsize=8, ha="center", color="#52606B")


def figure10(df, out, summary, retained):
    fig, axes = plt.subplots(1, 2, figsize=(13, 8))
    for ax in axes: ax.set(xlim=(0, 1), ylim=(0, 1)); ax.axis("off")
    a, b = axes
    a.text(.02, .99, "A   Observed developmental organization", fontsize=13, weight="bold", va="top")
    tree_box(a, (.5, .88), "CV: higher glucose (~2×)\nMW: lower glucose (~1×)", width=.7, height=.085, fontsize=9)
    tree_box(a, (.5, .73), "LHX6+/ERBB4+ DIV90 cortical interneurons", width=.83, height=.065, fontsize=9)
    arrow(a, (.5, .827), (.5, .771))
    tree_box(a, (.5, .59), "OBSERVED\nIndependent PV and SST program scores", width=.8, height=.09, fontsize=9)
    arrow(a, (.5, .688), (.5, .644))
    for x, state, short in zip([.17, .5, .83], [STATES[2], STATES[1], STATES[0]], ["SST-biased", "PV/SST hybrid", "PV-biased"]):
        n = (df.state == state).sum()
        tree_box(a, (x, .41), f"{short}\n{n:,} cells", color=COLORS[state] + "26", width=.28, height=.11, fontsize=9)
        arrow(a, (.5, .535), (x, .48))
    tree_box(a, (.5, .21), f"Unresolved / immature\n{int((df.state == STATES[3]).sum()):,} cells · low on both axes", width=.7, height=.1, color="#E9ECEF")
    a.plot([.083, .005, .005], [.59, .59, .21], color="#52606B", lw=1)
    a.add_patch(FancyArrowPatch((.005, .21), (.137, .21), arrowstyle="-|>", mutation_scale=11, lw=1, color="#52606B"))
    a.text(.5, .07, "INFERRED: operational program bias\nHybrid means concurrent program scores, not lineage direction.\nOperator and glucose are coupled; causality is unresolved.",
           ha="center", va="center", fontsize=9, color="#52606B")
    b.text(.02, .99, "B   Experimental sorting hypothesis", fontsize=13, weight="bold", va="top")
    tree_box(b, (.5, .87), "Viable DIV90 → ERBB4+\nMGE cortical interneuron preparation", width=.76, height=.1)
    selected = summary["selected"]
    rule_label = selected["gate_label"].replace(" AND ", "\nAND ").replace(" OR ", "\nOR ")
    tree_box(b, (.5, .65), f"HYPOTHESIS TO TEST\nRNA · log1p(CP10K)\n{rule_label}", width=.92, height=.16, fontsize=9)
    arrow(b, (.5, .81), (.5, .74))
    fractions_remove, _ = composition(df, ~retained)
    fractions_keep, _ = composition(df, retained)
    remove_text = "REMOVE\n" + "\n".join(f"{name}: {value:.1%}" for name, value in zip(SHORT, fractions_remove))
    keep_text = "RETAIN\n" + "\n".join(f"{name}: {value:.1%}" for name, value in zip(SHORT, fractions_keep))
    tree_box(b, (.25, .32), remove_text, color="#F7E6DB", width=.43, height=.28, fontsize=9)
    tree_box(b, (.75, .32), keep_text, color="#E3F0E9", width=.43, height=.28, fontsize=9)
    arrow(b, (.4, .56), (.25, .47), "rule-positive")
    arrow(b, (.6, .56), (.75, .47), "rule-negative")
    b.text(.5, .065, "RNA-predicted composition is observed in this dataset.\nProtein separation, maturation and physiological outcomes remain hypotheses.\nTrack unresolved cells separately in both fractions.",
           ha="center", va="center", fontsize=8.5, color="#52606B")
    fig.subplots_adjust(wspace=.14, left=.025, right=.99, top=.98, bottom=.015)
    return save(fig, out, "figure_10_experimental_model", "Vector logic schematic distinguishing observed RNA organization, inferred operational bias and prospective sorting hypotheses. Hybrid is not drawn as a lineage intermediate; unresolved cells remain explicit. Branch proportions are calculated from the exact selected RNA rule.")


def supplementary_overlays(df, out, lim, cuts):
    fields = SURFACE_COLORS + ["sample", "cell_line", "condition", "loupe_label"]
    fig, axes = plt.subplots(3, 4, figsize=(16, 12))
    for ax, field in zip(axes.ravel(), fields):
        landscape(ax, df, lim, cuts, color=field, point_size=3.5, alpha=.75, title=field)
    for ax in axes.ravel()[len(fields):]: ax.axis("off")
    fig.suptitle("Independent developmental axes with expression and source-metadata overlays", fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, .96))
    return save(fig, out, "supplementary_landscape_overlays", "Same developmental map independently colored by SST, ERBB4, FAT3, PTPRM, CXCR4, ACKR3, NRP2, sample, cell line, culture condition and the preserved Loupe annotation.")


HTML_TEMPLATE = r'''<!doctype html>
<html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>__TITLE__</title><style>
*{box-sizing:border-box}[hidden]{display:none!important}body{margin:0;background:#f5f7fa;font:14px system-ui,sans-serif;color:#203142}
header{padding:22px 30px;background:white;border-bottom:1px solid #dbe1e7}h1{font-size:23px;margin:0 0 6px;font-weight:650}p{margin:5px 0;line-height:1.45}.muted{color:#657380}
main{max-width:1500px;margin:18px auto;padding:0 22px}.controls{display:flex;flex-wrap:wrap;gap:12px;background:white;padding:16px;border:1px solid #dbe1e7;border-radius:9px;margin-bottom:14px}
label{display:flex;flex-direction:column;gap:5px;font-size:12px;font-weight:600}select,input,button{font:14px system-ui;padding:7px;border:1px solid #c5d0da;border-radius:5px;background:white;max-width:235px}
button{cursor:pointer}.plotbox{background:white;border:1px solid #dbe1e7;border-radius:9px;padding:15px;display:grid;grid-template-columns:minmax(400px,1fr) 270px;gap:16px}
canvas{width:100%;height:660px;display:block;touch-action:none}#legend{font-size:12px;line-height:1.7;max-height:200px;overflow:auto;margin-bottom:14px}.swatch{display:inline-block;width:10px;height:10px;border-radius:50%;margin-right:6px}
#hover{font-size:12px;line-height:1.6;overflow-wrap:anywhere}#hover b{color:#203142}#metrics{display:flex;gap:12px;flex-wrap:wrap;margin:12px 0}.metric{background:white;border:1px solid #dbe1e7;border-radius:8px;padding:12px 16px;min-width:175px}.metric strong{font-size:23px;display:block;margin:4px 0}.metric small{color:#657380}.gate{background:#edf4f0}.hint{font-size:12px;color:#657380;padding:14px 0}.error{color:#9b362d}
#comparison{overflow:auto;margin:12px 0;background:white;border-radius:8px}table{border-collapse:collapse;width:100%;font-size:12px}th,td{text-align:left;padding:10px 12px;border-bottom:1px solid #e2e7ec}th{background:#edf2f5}caption{text-align:left;font-size:14px;font-weight:600;padding:12px}
@media(max-width:850px){.plotbox{grid-template-columns:1fr}canvas{height:480px}.controls{gap:8px}header{padding:18px}main{padding:0 10px}}
</style><header><h1>__TITLE__</h1><p>__DESCRIPTION__</p><p class="muted">Frozen DIV90 LHX6+/ERBB4+ cortical cells · independent continuous programs · operational states, not fates</p></header>
<main><div class="controls" id="controls"></div><div class="controls gate" id="gatecontrols" hidden></div><div id="metrics"></div><div id="comparison"></div>
<div class="plotbox"><canvas id="plot" aria-label="Interactive single-cell scatterplot"></canvas><aside><div id="legend"></div><div id="hover">Hover over a cell to inspect its preserved identity, expression and assignments.</div></aside></div>
<p class="hint" id="hint">Self-contained offline asset: all cell data and code are embedded. Scroll over the plot to zoom; double-click to reset. Filters preserve global score thresholds.</p></main>
<script id="payload" type="application/json">__PAYLOAD__</script><script>
'use strict';
const P=JSON.parse(document.getElementById('payload').textContent), rows=P.rows, idx=Object.fromEntries(P.columns.map((x,i)=>[x,i]));
const stateColors={'PV-biased':'#326DAB','PV/SST hybrid':'#39876F','SST-biased':'#CB733F','unresolved/immature':'#A6ACB6'};
const palette=['#326DAB','#CB733F','#39876F','#A96C9D','#B89437','#68A6AC','#834E37','#7E83B7','#5C6C76','#A8A04A','#BF6888','#699348','#A48969','#5361B5','#CE9B90','#6B9080','#934A60','#6796BD','#B89D7E','#706C61'];
let shown=rows.map((_,i)=>i),color=P.mode==='gate'?'explored_gate':'state',xfield=P.mode==='umap'?'loupe_x':'sst_score',yfield=P.mode==='umap'?'loupe_y':'pv_score',hits=[],screen=[],bounds=null;
const canvas=document.getElementById('plot'),ctx=canvas.getContext('2d'),controls=document.getElementById('controls'),gatecontrols=document.getElementById('gatecontrols');
const el=(tag,props={})=>Object.assign(document.createElement(tag),props);
const escapeHtml=v=>String(v??'Unavailable').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const value=(i,key)=>key==='explored_gate'?(hits[i]?'Removed':'Retained'):rows[i][idx[key]];
function selector(parent,name,options,selected,onchange){const l=el('label',{textContent:name}),s=el('select');s.setAttribute('aria-label',name);for(const option of options)s.append(el('option',{value:option,textContent:option}));s.value=selected;s.addEventListener('change',()=>onchange(s.value));l.append(s);parent.append(l);return s;}
const available=P.colors.filter(k=>idx[k]!==undefined);if(P.mode==='gate')available.unshift('explored_gate');
selector(controls,'Color by',available,color,v=>{color=v;draw()});
const filters={};for(const field of ['cell_line','condition','sample']){const vals=[...new Set(rows.map(r=>String(r[idx[field]])))].sort();filters[field]=selector(controls,field.replace('_',' '),['All',...vals],'All',filter);}
controls.append(el('button',{textContent:'Reset view',onclick:()=>{bounds=null;draw()}}));
controls.append(el('button',{textContent:'Download plotted cells',onclick:download}));
function filter(){shown=rows.map((_,i)=>i).filter(i=>Object.entries(filters).every(([k,s])=>s.value==='All'||String(value(i,k))===s.value));draw();}
let gene1,gene2,op1,op2,t1,t2,logic,candidate,axisMode='Developmental landscape';
function ruleHit(v,op,t){return op==='>'?v>t:op==='>='?v>=t:op==='<'?v<t:v<=t;}
function thresholdControl(name,gene){const l=el('label',{textContent:name}),wrap=el('div'),r=el('input',{type:'range',step:'any'}),n=el('input',{type:'number',step:'any'});r.style.width='160px';n.style.width='90px';function sync(src,dst){dst.value=src.value;candidate.value='Custom';recompute();}r.oninput=()=>sync(r,n);n.oninput=()=>sync(n,r);wrap.append(r,n);l.append(wrap);gatecontrols.append(l);return {range:r,number:n,set(g,t){const a=rows.map(row=>Number(row[idx[g]])).filter(Number.isFinite),max=Math.max(...a);r.min=Math.min(0,...a);r.max=max;r.value=t;n.min=r.min;n.max=max;n.value=t;},get(){return Number(n.value)}};}
function recompute(){if(P.mode!=='gate')return;if(axisMode==='Surface-marker pair'){xfield=gene1.value;yfield=gene2.value==='None'?'ERBB4':gene2.value;bounds=null;}hits=rows.map((r,i)=>{const a=ruleHit(Number(value(i,gene1.value)),op1.value,t1.get());if(gene2.value==='None')return a;const b=ruleHit(Number(value(i,gene2.value)),op2.value,t2.get());return logic.value==='AND'?a&&b:a||b;});draw();}
if(P.mode==='gate'){
 gatecontrols.hidden=false;
 selector(controls,'Plot axes',['Developmental landscape','Surface-marker pair'],axisMode,v=>{axisMode=v;xfield=v==='Surface-marker pair'?gene1.value:'sst_score';yfield=v==='Surface-marker pair'?(gene2.value==='None'?'ERBB4':gene2.value):'pv_score';bounds=null;draw()});
 candidate=selector(gatecontrols,'Saved candidate',['Custom',...P.candidates.map(g=>g.name)],P.candidates[0].name,v=>{const found=P.candidates.find(g=>g.name===v);if(found)loadRule(found.rule)});
 gene1=selector(gatecontrols,'Depletion marker 1',P.gateGenes,P.gateGenes[0],g=>{t1.set(g,0);candidate.value='Custom';recompute()});
 op1=selector(gatecontrols,'Direction 1',['>','>=','<=','<'],'>',()=>{candidate.value='Custom';recompute()});t1=thresholdControl('Threshold 1');
 logic=selector(gatecontrols,'Combine',['AND','OR'],'AND',()=>{candidate.value='Custom';recompute()});
 gene2=selector(gatecontrols,'Depletion marker 2',['None',...P.gateGenes],'None',g=>{if(g!=='None')t2.set(g,0);candidate.value='Custom';recompute()});
 op2=selector(gatecontrols,'Direction 2',['>','>=','<=','<'],'>',()=>{candidate.value='Custom';recompute()});t2=thresholdControl('Threshold 2');
 function loadRule(rule){gene1.value=rule.rules[0].gene;op1.value=rule.rules[0].op;t1.set(gene1.value,rule.rules[0].threshold);logic.value=rule.logic||'AND';const second=rule.rules[1];gene2.value=second?second.gene:'None';op2.value=second?second.op:'>';t2.set(second?second.gene:gene1.value,second?second.threshold:0);recompute();}
 loadRule(P.candidates[0].rule);
 document.getElementById('hint').textContent+=' Marker-positive cells matching the displayed rule are removed. Metrics update for the selected cell subset. RNA cutoffs are not fluorescence thresholds; this explorer does not establish protein separation.';
}
function numeric(key){return key!=='explored_gate'&&rows.some(r=>typeof r[idx[key]]==='number');}
function draw(){
 const rect=canvas.getBoundingClientRect(),dpr=window.devicePixelRatio||1,w=rect.width,h=rect.height;canvas.width=Math.round(w*dpr);canvas.height=Math.round(h*dpr);ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,w,h);
 const m={left:67,right:20,top:20,bottom:58},pw=w-m.left-m.right,ph=h-m.top-m.bottom;
 if(!bounds){const xs=rows.map((_,i)=>value(i,xfield)),ys=rows.map((_,i)=>value(i,yfield));let xmin=Math.min(...xs),xmax=Math.max(...xs),ymin=Math.min(...ys),ymax=Math.max(...ys);const px=(xmax-xmin)*.04,py=(ymax-ymin)*.04;bounds=[xmin-px,xmax+px,ymin-py,ymax+py];}
 const [xmin,xmax,ymin,ymax]=bounds,sx=x=>m.left+(x-xmin)/(xmax-xmin)*pw,sy=y=>m.top+ph-(y-ymin)/(ymax-ymin)*ph;
 ctx.font='11px system-ui';ctx.fillStyle='#5c6b78';ctx.strokeStyle='#e6eaef';ctx.lineWidth=1;
 for(let j=0;j<=5;j++){let x=xmin+(xmax-xmin)*j/5,y=ymin+(ymax-ymin)*j/5;ctx.beginPath();ctx.moveTo(sx(x),m.top);ctx.lineTo(sx(x),m.top+ph);ctx.stroke();ctx.beginPath();ctx.moveTo(m.left,sy(y));ctx.lineTo(m.left+pw,sy(y));ctx.stroke();ctx.textAlign='center';ctx.fillText(x.toFixed(2),sx(x),h-35);ctx.textAlign='right';ctx.fillText(y.toFixed(2),m.left-8,sy(y)+4);}
 if(P.mode!=='umap'){ctx.setLineDash([4,4]);ctx.strokeStyle='#73808b';ctx.beginPath();const cutx=axisMode==='Surface-marker pair'?t1.get():P.cuts[0],cuty=axisMode==='Surface-marker pair'?(gene2.value==='None'?null:t2.get()):P.cuts[1];ctx.moveTo(sx(cutx),m.top);ctx.lineTo(sx(cutx),m.top+ph);if(cuty!==null){ctx.moveTo(m.left,sy(cuty));ctx.lineTo(m.left+pw,sy(cuty));}ctx.stroke();ctx.setLineDash([]);}
 let cats=[],vmin=0,vmax=1,isnum=numeric(color);if(isnum){const vals=shown.map(i=>value(i,color)).filter(Number.isFinite).sort((a,b)=>a-b);vmin=vals[0]??0;vmax=vals[Math.floor((vals.length-1)*.99)]??1;if(vmax===vmin)vmax=vmin+1;}else cats=[...new Set(shown.map(i=>String(value(i,color))))].sort();
 function col(i){const v=value(i,color);if(isnum){if(v===null)return '#b4bcc5';const q=Math.max(0,Math.min(1,(v-vmin)/(vmax-vmin)));return `rgb(${Math.round(234-195*q)},${Math.round(239-130*q)},${Math.round(243-85*q)})`;}if(color==='state')return stateColors[v]||'#999';if(color.includes('gate')||color.includes('retained'))return ['Removed','Lost','false','0'].includes(String(v))?'#CB733F':'#39876F';return palette[cats.indexOf(String(v))%palette.length];}
 screen=[];ctx.save();ctx.beginPath();ctx.rect(m.left,m.top,pw,ph);ctx.clip();const order=isnum?[...shown].sort((a,b)=>(value(a,color)||0)-(value(b,color)||0)):shown;
 for(const i of order){const x=sx(value(i,xfield)),y=sy(value(i,yfield));ctx.globalAlpha=.73;ctx.fillStyle=col(i);ctx.beginPath();ctx.arc(x,y,2.1,0,Math.PI*2);ctx.fill();if(x>=m.left&&x<=m.left+pw&&y>=m.top&&y<=m.top+ph)screen.push([x,y,i]);}ctx.restore();ctx.globalAlpha=1;
 ctx.fillStyle='#203142';ctx.font='14px system-ui';ctx.textAlign='center';ctx.fillText(P.mode==='umap'?'Preserved Loupe UMAP 1':axisMode==='Surface-marker pair'?xfield+' RNA · log1p(CP10K)':'SST developmental program',m.left+pw/2,h-9);ctx.save();ctx.translate(17,m.top+ph/2);ctx.rotate(-Math.PI/2);ctx.fillText(P.mode==='umap'?'Preserved Loupe UMAP 2':axisMode==='Surface-marker pair'?yfield+' RNA · log1p(CP10K)':'PV developmental program',0,0);ctx.restore();
 const legend=document.getElementById('legend');legend.replaceChildren(el('b',{textContent:color}));if(isnum){legend.append(el('div',{textContent:`${vmin.toFixed(3)} to ${vmax.toFixed(3)} (99th-percentile color cap)`}));const gradient=el('div');gradient.style.cssText='height:14px;background:linear-gradient(to right,rgb(234,239,243),rgb(39,109,158));margin:8px 0';legend.append(gradient);}else for(const cat of cats){const d=el('div'),s=el('span',{className:'swatch'}),representative=shown.find(i=>String(value(i,color))===cat);s.style.background=col(representative);d.append(s,document.createTextNode(cat));legend.append(d);}
 const metrics=document.getElementById('metrics');metrics.replaceChildren();function metric(name,number,detail){const div=el('div',{className:'metric'});div.append(el('div',{textContent:name}),el('strong',{textContent:number}),el('small',{textContent:detail}));metrics.append(div);}
 metric('Cells displayed',shown.length.toLocaleString(),`${rows.length.toLocaleString()} frozen starting cells`);
 if(P.mode==='gate'){const retained=shown.filter(i=>!hits[i]),count=(arr,s)=>arr.filter(i=>value(i,'state')===s).length,basePV=count(shown,'PV-biased'),baseH=count(shown,'PV/SST hybrid'),rPV=count(retained,'PV-biased'),rH=count(retained,'PV/SST hybrid'),rS=count(retained,'SST-biased');const pct=(n,d)=>d?(100*n/d).toFixed(1)+'%':'Undefined';metric('PV-biased recovery',pct(rPV,basePV),`${basePV-rPV} PV-biased cells lost`);metric('Hybrid recovery',pct(rH,baseH),`${baseH-rH} hybrid cells lost`);metric('SST contamination',pct(rS,retained.length),`${rS} SST-biased cells retained`);metric('Total retained',retained.length.toLocaleString(),`${pct(retained.length,shown.length)} of displayed cells`);metric('Unresolved retained',count(retained,'unresolved/immature').toLocaleString(),'Kept separate from lineage-associated targets');}
 if(P.mode==='gate')comparisonTable();
}
function comparisonTable(){const container=document.getElementById('comparison'),table=el('table');table.append(el('caption',{textContent:'Same thresholds by cell line and culture condition (descriptive; operator and glucose are coupled)'}));const tr=el('tr');for(const h of ['Cell line','Culture condition','Starting n','Retained n','PV recovery','Hybrid recovery','SST contamination','Unresolved retained'])tr.append(el('th',{textContent:h}));table.append(tr);const groups={};for(const i of shown){const key=String(value(i,'cell_line'))+'\t'+String(value(i,'condition'));(groups[key]??=[]).push(i);}for(const [key,arr] of Object.entries(groups).sort()){const kept=arr.filter(i=>!hits[i]),n=(a,s)=>a.filter(i=>value(i,'state')===s).length,pct=(a,b)=>b?(100*a/b).toFixed(1)+'%':'Undefined',vals=[...key.split('\t'),arr.length,kept.length,pct(n(kept,'PV-biased'),n(arr,'PV-biased')),pct(n(kept,'PV/SST hybrid'),n(arr,'PV/SST hybrid')),pct(n(kept,'SST-biased'),kept.length),n(kept,'unresolved/immature')],row=el('tr');for(const v of vals)row.append(el('td',{textContent:v}));table.append(row);}container.replaceChildren(table);}
canvas.addEventListener('mousemove',e=>{const r=canvas.getBoundingClientRect(),x=e.clientX-r.left,y=e.clientY-r.top;let best=null,dist=100;for(const p of screen){const d=(p[0]-x)**2+(p[1]-y)**2;if(d<dist){dist=d;best=p;}}if(!best)return;const i=best[2];const fields=['cell_id','sample','cell_line','condition','loupe_label','SST','pv_score','sst_score','ERBB4','FAT3','PTPRM','CXCR4','ACKR3','NRP2','state','gate_assignment','phase1_gate_assignment'];let t='<b>Cell details</b><br>';for(const f of fields){let v=value(i,f);if(typeof v==='number')v=v.toPrecision(5);t+='<b>'+escapeHtml(f)+':</b> '+escapeHtml(v)+'<br>';}if(P.mode==='gate')t+='<b>Explored gate:</b> '+(hits[i]?'Removed':'Retained');document.getElementById('hover').innerHTML=t;});
canvas.addEventListener('wheel',e=>{e.preventDefault();const r=canvas.getBoundingClientRect(),fx=Math.max(0,Math.min(1,(e.clientX-r.left-67)/(r.width-87))),fy=1-Math.max(0,Math.min(1,(e.clientY-r.top-20)/(r.height-78))),[a,b,c,d]=bounds,cx=a+(b-a)*fx,cy=c+(d-c)*fy,k=e.deltaY>0?1.15:1/1.15;bounds=[cx+(a-cx)*k,cx+(b-cx)*k,cy+(c-cy)*k,cy+(d-cy)*k];draw();},{passive:false});
canvas.addEventListener('dblclick',()=>{bounds=null;draw()});window.addEventListener('resize',draw);
function download(){const cols=[...P.columns,...(P.mode==='gate'?['explored_gate']:[])],text=[cols.join('\t'),...shown.map(i=>cols.map(k=>String(value(i,k)??'').replace(/[\t\r\n]/g,' ')).join('\t'))].join('\n');const url=URL.createObjectURL(new Blob([text],{type:'text/tab-separated-values'})),a=el('a',{href:url,download:'div90_displayed_cells.tsv'});a.click();URL.revokeObjectURL(url);}
draw();
</script></html>'''


def interactive(df, out, summary, capture, cuts):
    cells = df.copy()
    cells["gate_assignment"] = np.where(capture["selected"], "Retained", "Removed")
    cells["phase1_gate_assignment"] = np.where(capture["phase1"], "Captured", "Lost")
    gate_genes = set(SURFACE_COLORS[1:])
    candidates = []
    for key, display in [("selected", "Selected practical gate"), ("best_single", "Best one-marker"),
                         ("best_pair", "Best two-marker"), ("phase1_reference", "Phase-one FAT3/PTPRM complement"),
                         ("erbb4_cxcr4", "New ERBB4/CXCR4 REMOVE"),
                         ("phase1_erbb4_cxcr4_retention", "Inherited ERBB4/CXCR4 RETAIN (complement shown)"),
                         ("unconstrained_selected", "Unconstrained preservation optimum")]:
        if not summary.get(key): continue
        rule = json.loads(json.dumps(summary["rules"][str(summary[key]["gate_id"])]))
        if rule.get("action") == "retain":
            rule["action"] = "remove"
            rule["logic"] = "OR" if rule.get("logic") == "AND" else "AND"
            complements = {"<=": ">", ">": "<=", ">=": "<", "<": ">="}
            for item in rule["rules"]: item["op"] = complements[item["op"]]
        if len(rule["rules"]) <= 2:
            candidates.append({"name": display, "rule": rule})
        gate_genes.update(item["gene"] for item in rule["rules"])
    ranking = pd.read_csv(out / "tables/surface_state_markers.tsv", sep="\t")
    for state in STATES[:3]:
        gate_genes.update(ranking.loc[ranking.state == state].nlargest(8, "cohens_d").gene)
    gate_genes = sorted(g for g in gate_genes if g in cells)
    fields = ["cell_id", "sample", "cell_line", "condition", "loupe_label", "loupe_x", "loupe_y",
              "pv_score", "sst_score", "state", "gate_assignment", "phase1_gate_assignment"] + SURFACE_COLORS + gate_genes
    fields = list(dict.fromkeys(f for f in fields if f in cells))
    small = cells[fields].copy()
    # Preserve full float precision: rounding would change cells on gate boundaries.
    # Object conversion permits JSON null for unavailable values rather than NaN.
    data = small.astype(object).where(pd.notna(small), None).to_dict(orient="split")
    common = {"columns": data["columns"], "rows": data["data"], "cuts": cuts,
              "colors": ["state", "pv_score", "sst_score"] + SURFACE_COLORS +
                        ["sample", "cell_line", "condition", "loupe_label", "gate_assignment", "phase1_gate_assignment"],
              "gateGenes": gate_genes, "candidates": candidates}
    assets = [("interactive_pv_sst_landscape.html", "landscape", "DIV90 developmental landscape", "Explore PV and SST programs with independent expression and metadata overlays."),
              ("interactive_umap.html", "umap", "Preserved cortical recluster UMAP", "Original recovered Loupe coordinates are displayed unchanged; geometry is not lineage evidence."),
              ("interactive_gate_explorer.html", "gate", "Interactive SST-depletion gate explorer", "Move transcript thresholds and see the biological cost in PV-biased and hybrid cells immediately."),
              ("interactive_sample_landscape.html", "sample", "Cell-line and culture-condition landscapes", "Select individual cell lines, culture conditions and samples while retaining the shared developmental axes.")]
    for filename, mode, title, description in assets:
        description += f" Displayed experimental gate: {summary['selected']['gate_label']}."
        if summary.get("visual_depletion_constraint"):
            description += " This practical comparison requires at least 20% SST-biased removal, an analyst-selected constraint."
        payload = json.dumps({**common, "mode": mode}, separators=(",", ":"), allow_nan=False).replace("</", "<\\/")
        content = HTML_TEMPLATE.replace("__TITLE__", html.escape(title)).replace("__DESCRIPTION__", html.escape(description)).replace("__PAYLOAD__", payload)
        (out / filename).write_text(content)
    return [a[0] for a in assets]


def run(cells=None, out=OUT, summary=None):
    out = Path(out).expanduser().resolve()
    if out != OUT.resolve():
        raise ValueError(f"Phase-two visuals must be written to the authorized Turbo directory: {OUT}")
    out.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(cells or out / "cells.tsv.gz", sep="\t") if not isinstance(cells, pd.DataFrame) else cells.copy()
    df = df.reset_index(drop=True)
    if len(df) != 4768 or df.cell_id.duplicated().any():
        raise ValueError("Expected exactly 4,768 unique frozen entry cells")
    if set(df.state) != set(STATES):
        raise ValueError(f"Unexpected operational state names: {set(df.state)}")
    summary = summary or json.loads((out / "gate_summary.json").read_text())
    summary = visual_summary(summary)
    configure(); lim = limits(df); cuts = thresholds(df); capture = masks(df, out, summary)
    captions = []
    calls = [lambda: figure1(df, out, lim, cuts), lambda: figure2(df, out, lim, cuts),
             lambda: figure3(df, out), lambda: figure4(df, out),
             lambda: figure5(df, out, lim, cuts, capture["phase1"]),
             lambda: figure6(df, out, summary, capture["selected"]),
             lambda: figure7(df, out, capture["selected"]), lambda: figure8(df, out, capture["selected"]),
             lambda: figure9(df, out, summary), lambda: figure10(df, out, summary, capture["selected"]),
             lambda: supplementary_overlays(df, out, lim, cuts)]
    for create in calls:
        caption = create(); captions.append(caption)
        print(f"Saved {caption['figure']} (PDF/SVG/600-dpi PNG)", flush=True)
    assets = interactive(df, out, summary, capture, cuts)
    manifest = {"n_cells": len(df), "score_boundary_sst": cuts[0], "score_boundary_pv": cuts[1],
                "png_dpi": 600, "vector_formats": ["pdf", "svg"], "figures": captions,
                "interactive_assets": assets, "offline": True, "umap_recomputed": False,
                "displayed_experimental_gate_id": summary["selected"]["gate_id"],
                "visual_depletion_constraint": summary.get("visual_depletion_constraint")}
    (out / "figure_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (out / "figures/FIGURE_CAPTIONS.md").write_text("# DIV90 phase-two figure captions\n\n" +
        "\n\n".join(f"**{item['figure']}**\n\n{item['caption']}" for item in captions) + "\n")
    print("Saved all four standalone offline HTML assets", flush=True)
    return manifest


def visual_summary(summary):
    """Use the explicitly depletion-constrained strategy for experimental panels."""
    summary = json.loads(json.dumps(summary))
    if summary.get("experimental_depletion"):
        summary["unconstrained_selected"] = summary["selected"]
        summary["selected"] = summary["experimental_depletion"]
        summary["best_single"] = summary["practical_single"]
        summary["best_pair"] = summary["practical_pair"]
        summary["visual_depletion_constraint"] = summary.get("practical_depletion_constraint", "At least 20% SST-biased removal, an analyst-selected threshold")
    inherited = summary.get("phase1_erbb4_cxcr4_benchmarks", [])
    if inherited:
        summary["phase1_erbb4_cxcr4_retention"] = next((row for row in inherited if "module" in row.get("source", "")), inherited[0])
    return summary


def verify_html_assets(out=OUT):
    """Exercise offline rendering, filtering, hover and exact live gate counts."""
    out = Path(out).resolve()
    if out != OUT.resolve():
        raise ValueError("Interactive QA artifacts must remain in the authorized Turbo directory")
    sys.path.insert(0, str(out / "runtime_deps"))
    os.environ["PLAYWRIGHT_BROWSERS_PATH"] = str(out / "cache/playwright")
    (out / "cache/browser_tmp").mkdir(parents=True, exist_ok=True)
    os.environ["TMPDIR"] = str(out / "cache/browser_tmp")
    import tempfile
    tempfile.tempdir = os.environ["TMPDIR"]
    shim = out / "provenance/playwright_autofs_shim.cjs"
    if shim.exists():
        os.environ["NODE_OPTIONS"] = f"--require={shim}"
    from playwright.sync_api import sync_playwright
    df = pd.read_csv(out / "cells.tsv.gz", sep="\t")
    summary = visual_summary(json.loads((out / "gate_summary.json").read_text()))
    selected_retained = apply_rule(df, summary["rules"][str(summary["selected"]["gate_id"])])
    report = {"offline": True, "assets": []}
    with sync_playwright() as pw:
        browser = pw.chromium.launch(headless=True, args=["--no-sandbox"])
        context = browser.new_context(offline=True, viewport={"width": 1440, "height": 1100})
        for filename in ["interactive_pv_sst_landscape.html", "interactive_umap.html",
                         "interactive_gate_explorer.html", "interactive_sample_landscape.html"]:
            page = context.new_page(); errors = []
            page.on("pageerror", lambda error: errors.append(str(error)))
            page.goto((out / filename).as_uri(), wait_until="load")
            page.wait_for_function("typeof shown !== 'undefined' && shown.length === 4768")
            assert not errors, errors
            assert page.evaluate("screen.length") > 0, filename
            page.get_by_label("Color by", exact=True).select_option("SST")
            assert "SST" in page.locator("#legend").inner_text()
            point = page.evaluate("screen[0]")
            page.locator("#plot").scroll_into_view_if_needed()
            box = page.locator("#plot").bounding_box()
            page.mouse.move(box["x"] + point[0], box["y"] + point[1])
            assert "cell_id:" in page.locator("#hover").inner_text(), filename
            condition = df.condition.iloc[0]
            page.get_by_label("condition", exact=True).select_option(condition)
            expected_n = int((df.condition == condition).sum())
            assert page.evaluate("shown.length") == expected_n
            page.get_by_label("condition", exact=True).select_option("All")
            checks = {"filename": filename, "n_cells": 4768, "coloring": True,
                      "hover": True, "condition_filter_n": expected_n, "errors": errors}
            if "gate_explorer" in filename:
                assert page.evaluate("hits.filter(x => !x).length") == int(selected_retained.sum())
                inherited = summary.get("phase1_erbb4_cxcr4_retention")
                if inherited:
                    page.get_by_label("Saved candidate", exact=True).select_option("Inherited ERBB4/CXCR4 RETAIN (complement shown)")
                    expected_inherited = int(apply_rule(df, summary["rules"][str(inherited["gate_id"])]).sum())
                    assert page.evaluate("hits.filter(x => !x).length") == expected_inherited
                    checks["inherited_erbb4_cxcr4_retained_n"] = expected_inherited
                    page.get_by_label("Saved candidate", exact=True).select_option("Selected practical gate")
                changes = page.evaluate("""(() => {
                    const marker=gene1.value, cutoff=Number(t1.range.max)*0.4;
                    gene2.value='None';op1.value='>';t1.number.value=cutoff;t1.range.value=cutoff;recompute();
                    return {marker,cutoff,retained:hits.filter(x=>!x).length};
                })()""")
                assert changes["retained"] == int((df[changes["marker"]] <= changes["cutoff"]).sum())
                page.get_by_label("condition", exact=True).select_option(condition)
                actual_retained = page.evaluate("shown.filter(i=>!hits[i]).length")
                assert actual_retained == int(((df.condition == condition) & (df[changes["marker"]] <= changes["cutoff"])).sum())
                page.get_by_label("condition", exact=True).select_option("All")
                page.get_by_label("Saved candidate", exact=True).select_option("Selected practical gate")
                page.get_by_label("Color by", exact=True).select_option("explored_gate")
                page.get_by_label("Plot axes", exact=True).select_option("Surface-marker pair")
                assert page.evaluate("xfield") == page.evaluate("gene1.value")
                page.get_by_label("Plot axes", exact=True).select_option("Developmental landscape")
                checks.update({"exact_selected_retained_n": int(selected_retained.sum()),
                               "live_threshold_counts": True, "condition_gate_counts": True,
                               "surface_axes": True, "comparison_rows": page.locator("#comparison tr").count() - 1})
                page.screenshot(path=str(out / "provenance/interactive_gate_explorer_preview.png"), full_page=True)
            assert not errors, errors
            report["assets"].append(checks)
            page.close()
        browser.close()
    (out / "provenance/interactive_validation.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", type=Path)
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--verify-html", action="store_true")
    args = parser.parse_args()
    if args.verify_html:
        print(json.dumps(verify_html_assets(args.out), indent=2))
    else:
        run(args.cells, args.out)
