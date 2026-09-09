#!/usr/bin/env python3
"""Four-state surface discovery, actual R Hypergate, and positive depletion rules.

All model selection and metrics are descriptive, in-sample calculations. This
module deliberately does not modify any phase-one rule, label, or result.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import os
from pathlib import Path
import subprocess
import time

import h5py
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from threadpoolctl import threadpool_limits

from div90_hypergate_paths import PROJECT_ROOT, REPO_ROOT

STATES = ["PV-biased", "PV/SST hybrid", "SST-biased", "unresolved/immature"]
SHORT = ["pv", "hybrid", "sst", "unresolved"]
PHASE1 = PROJECT_ROOT / "results/div90_hypergate_sst_pv"
EXPECTED_OUT = PROJECT_ROOT / "results/div90_hypergate_sst_pv_phase2"
R_BIN = "/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/Rscript"


def safe_div(a, b):
    a, b = np.broadcast_arrays(np.asarray(a, dtype=float), np.asarray(b, dtype=float))
    return np.divide(a, b, out=np.full(a.shape, np.nan), where=b != 0)


def clean_json(obj):
    if isinstance(obj, dict):
        return {str(k): clean_json(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [clean_json(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return clean_json(obj.tolist())
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (float, np.floating)):
        return float(obj) if np.isfinite(obj) else None
    if isinstance(obj, np.bool_):
        return bool(obj)
    return obj


def dump_json(path, value):
    path.write_text(json.dumps(clean_json(value), indent=2, allow_nan=False) + "\n")


def metrics(removed, totals):
    """Metrics for BOTH fractions; unresolved has two explicit treatments."""
    r = np.atleast_2d(np.asarray(removed, dtype=np.int64))
    t = np.asarray(totals, dtype=np.int64)
    k = t[None, :] - r
    rn, kn = r.sum(axis=1), k.sum(axis=1)
    d = {"removed_n": rn, "retained_n": kn}
    for j, name in enumerate(SHORT):
        d[f"removed_{name}_n"] = r[:, j]
        d[f"retained_{name}_n"] = k[:, j]
        d[f"removed_{name}_fraction"] = safe_div(r[:, j], rn)
        d[f"retained_{name}_fraction"] = safe_div(k[:, j], kn)
        d[f"{name}_recovery"] = safe_div(k[:, j], t[j])
        d[f"{name}_removal_recovery"] = safe_div(r[:, j], t[j])
    d["pv_lost_n"] = r[:, 0]
    d["hybrid_lost_n"] = r[:, 1]
    d["target_recovery"] = safe_div(k[:, 0] + k[:, 1], t[0] + t[1])
    d["total_recovery"] = safe_div(kn, t.sum())
    d["sst_contamination"] = d["retained_sst_fraction"]
    d["retained_resolved_n"] = k[:, :3].sum(axis=1)
    d["sst_contamination_unresolved_excluded"] = safe_div(k[:, 2], k[:, :3].sum(axis=1))
    d["target_fraction_unresolved_excluded"] = safe_div(k[:, 0] + k[:, 1], k[:, :3].sum(axis=1))
    d["pv_fraction_unresolved_excluded"] = safe_div(k[:, 0], k[:, :3].sum(axis=1))
    d["hybrid_fraction_unresolved_excluded"] = safe_div(k[:, 1], k[:, :3].sum(axis=1))
    d["total_recovery_unresolved_excluded"] = safe_div(k[:, :3].sum(axis=1), t.sum())
    balance = np.sqrt(d["pv_recovery"] * d["hybrid_recovery"])
    d["balanced_target_recovery"] = balance
    d["practical_score"] = balance * (1 - d["sst_contamination_unresolved_excluded"])
    d["practical_score_unresolved_retained"] = balance * (1 - d["sst_contamination"])
    return d


def apply_rule(cells, rule):
    masks = []
    for item in rule["rules"]:
        x = cells[item["gene"]].to_numpy(float)
        v = item["threshold"]
        masks.append({">": np.greater, ">=": np.greater_equal,
                      "<": np.less, "<=": np.less_equal}[item["op"]](x, v))
    if not masks:
        mask = np.zeros(len(cells), dtype=bool) if rule["action"] == "remove" else np.ones(len(cells), dtype=bool)
    elif rule["logic"] == "OR":
        mask = np.logical_or.reduce(masks)
    else:
        mask = np.logical_and.reduce(masks)
    return mask if rule["action"] == "remove" else ~mask


def rule_id(rule):
    return hashlib.sha256(json.dumps(rule, sort_keys=True).encode()).hexdigest()[:16]


def gate_label(rule):
    if not rule["rules"]:
        return "retain all entry cells"
    return (" " + rule["logic"] + " ").join(
        f"{x['gene']} {x['op']} {x['threshold']:.8g}" for x in rule["rules"])


def metric_row(rule, removed_counts, totals, source):
    result = {k: float(v[0]) for k, v in metrics(removed_counts, totals).items()}
    result.update(gate_id=rule_id(rule), source=source, action=rule["action"],
                  logic=rule["logic"], n_markers=len({r["gene"] for r in rule["rules"]}),
                  gate_label=gate_label(rule), rules=json.dumps(rule["rules"], separators=(",", ":")))
    return result


def pareto_positions(d):
    """Maximize combined recovery, minimize all-retained SST contamination."""
    x, y = d["target_recovery"], d["sst_contamination"]
    good = np.isfinite(x) & np.isfinite(y)
    ix = np.flatnonzero(good)
    order = ix[np.lexsort((y[ix], -x[ix]))]
    out, best_y = [], np.inf
    for i in order:
        if y[i] < best_y - 1e-12:
            out.append(i)
            best_y = y[i]
    return np.asarray(out, dtype=int)


def positions_to_keep(d, top=1500):
    score = np.nan_to_num(d["practical_score"], nan=-1.)
    alt = np.nan_to_num(d["practical_score_unresolved_retained"], nan=-1.)
    n = min(top, len(score))
    a = np.argpartition(score, len(score) - n)[-n:]
    b = np.argpartition(alt, len(alt) - n)[-n:]
    constrained = []
    for minimum in [.1,.2,.3,.5]:
        allowed = np.flatnonzero(d["sst_removal_recovery"]>=minimum)
        if len(allowed):
            k = min(200,len(allowed))
            constrained.extend(allowed[np.argpartition(score[allowed],len(allowed)-k)[-k:]])
    return np.unique(np.r_[a, b, pareto_positions(d),constrained]).astype(int)


def write_metric_group(h5, name, d, **definitions):
    g = h5.create_group(name)
    keys = list(d)
    g.attrs["metric_columns"] = json.dumps(keys)
    # Float64 preserves exact ratios, and all counts are far below 2**53.
    values = np.column_stack([d[k] for k in keys])
    g.create_dataset("metrics", data=values, compression="gzip", compression_opts=3, shuffle=True)
    for key, value in definitions.items():
        g.create_dataset(key, data=np.asarray(value), compression="gzip", compression_opts=3, shuffle=True)


def marker_screen(cells, genes, annotations, out):
    """Gene-level effect sizes + optimal signed one-marker state enrichments."""
    state = cells.state.to_numpy()
    effects, screens = [], []
    for gene in genes:
        x = cells[gene].to_numpy(float)
        ranks = rankdata(x)
        idx = np.argsort(x, kind="stable")
        unique, starts, freq = np.unique(x[idx], return_index=True, return_counts=True)
        for target in STATES:
            y = state == target
            n, m = int(y.sum()), int((~y).sum())
            a, b = x[y], x[~y]
            sd = np.sqrt(((n - 1) * np.var(a, ddof=1) + (m - 1) * np.var(b, ddof=1)) / (n + m - 2))
            effects.append(dict(gene=gene, state=target, n_state=n, n_other=m,
                                mean_logexpr_state=a.mean(), mean_logexpr_other=b.mean(),
                                mean_difference=a.mean()-b.mean(), detection_state=(a>0).mean(),
                                detection_other=(b>0).mean(), detection_difference=(a>0).mean()-(b>0).mean(),
                                cohens_d=(a.mean()-b.mean())/sd if sd else np.nan,
                                auc=(ranks[y].sum()-n*(n+1)/2)/(n*m)))
        queries = [(f"Q{j+1}", target, np.ones(len(cells), dtype=bool))
                   for j, target in enumerate([STATES[0], STATES[2], STATES[1]])]
        queries.append(("Q4", STATES[0], np.isin(state, [STATES[0], STATES[2]])))
        for question, target, allowed in queries:
            xx = x[allowed]
            yy = (state[allowed] == target).astype(int)
            order = np.argsort(xx, kind="stable")
            ux, st, counts = np.unique(xx[order], return_index=True, return_counts=True)
            positives = np.add.reduceat(yy[order], st)
            low_tp, low_n = np.cumsum(positives), np.cumsum(counts)
            for op, tp, selected in [("<=", low_tp, low_n), (">", yy.sum()-low_tp, len(yy)-low_n)]:
                fp, fn = selected-tp, yy.sum()-tp
                f1 = safe_div(2*tp, 2*tp+fp+fn)
                f1 = np.nan_to_num(f1, nan=0)
                best = int(np.argmax(f1))
                screens.append(dict(gene=gene, question=question, target=target, op=op,
                                    threshold=ux[best], tp=tp[best], fp=fp[best], fn=fn[best],
                                    n_evaluated=len(yy), n_hybrid_excluded=int((state==STATES[1]).sum()) if question=="Q4" else 0,
                                    purity=safe_div(tp[best], selected[best]).item(),
                                    recovery=safe_div(tp[best], yy.sum()).item(), F1=f1[best]))
    eff = pd.DataFrame(effects).merge(annotations, on="gene", how="left")
    eff.to_csv(out / "tables/surface_state_markers.tsv", sep="\t", index=False)
    screen = pd.DataFrame(screens)
    screen.to_csv(out / "tables/surface_univariate_state_screen.tsv", sep="\t", index=False)
    return eff, screen


def exact_single_screen(cells, genes, totals, state_index, h5, out):
    rows, registry = [], {}
    counts_all, gene_indices, threshold_all = [], [], []
    for gi, gene in enumerate(genes):
        x = cells[gene].to_numpy(float)
        order = np.argsort(x, kind="stable")
        thresholds, starts = np.unique(x[order], return_index=True)
        per_value = np.stack([np.add.reduceat((state_index[order]==j).astype(int), starts) for j in range(4)], axis=1)
        if thresholds[0]>0:
            thresholds = np.r_[0.,thresholds]
            per_value = np.vstack([np.zeros(4,dtype=int),per_value])
        removed = totals[None, :] - np.cumsum(per_value, axis=0)
        # All observed cut points, including RNA detection > 0 when available.
        counts_all.append(removed)
        gene_indices.append(np.full(len(thresholds), gi, dtype=np.int16))
        threshold_all.append(thresholds)
    removed = np.concatenate(counts_all)
    gidx = np.concatenate(gene_indices)
    thresh = np.concatenate(threshold_all)
    d = metrics(removed, totals)
    write_metric_group(h5, "exact_single", d, gene_index=gidx, threshold=thresh)
    # Human-readable complete exact one-marker metric table.
    table = pd.DataFrame(d)
    table.insert(0, "threshold", thresh)
    table.insert(0, "op", ">")
    table.insert(0, "gene", np.asarray(genes)[gidx])
    table.to_csv(out / "tables/depletion_single_all.tsv.gz", sep="\t", index=False, compression={"method":"gzip", "compresslevel":3})
    best_by_gene = table.groupby("gene", sort=False)["practical_score"].idxmax().to_numpy()
    chosen = np.unique(np.r_[positions_to_keep(d), best_by_gene])
    for i in chosen:
        rule = dict(action="remove", logic="AND", rules=[dict(gene=genes[gidx[i]], op=">", threshold=float(thresh[i]))])
        row = metric_row(rule, removed[i], totals, "python_exact_single")
        registry[row["gate_id"]] = rule
        rows.append(row)
    best_gene_table = table.iloc[best_by_gene].sort_values("practical_score", ascending=False)
    best_gene_table.to_csv(out / "tables/depletion_single_best_by_gene.tsv", sep="\t", index=False)
    return rows, registry, len(table)


def broad_pair_screen(cells, genes, totals, state_index, h5, out):
    definitions, binary = [], []
    for gene in genes:
        x = cells[gene].to_numpy(float)
        positive = x[x>0]
        thresholds = np.unique(np.r_[0., np.quantile(positive, [.25, .5, .75, .9])])
        for threshold in thresholds:
            definitions.append(dict(gene=gene, op=">", threshold=float(threshold)))
            binary.append(x>threshold)
    definition_table = pd.DataFrame(definitions)
    definition_table.to_csv(out / "tables/pair_threshold_grid.tsv", sep="\t", index=False)
    h5.attrs["pair_threshold_definitions"] = json.dumps(definitions)
    b = np.asarray(binary, dtype=np.float32)
    single = np.stack([b[:, state_index==j].sum(axis=1).astype(np.int32) for j in range(4)], axis=1)
    intersections = []
    with threadpool_limits(limits=4):
        for j in range(4):
            z = b[:, state_index==j]
            intersections.append(np.rint(z @ z.T).astype(np.int32))
    aa, bb = np.triu_indices(len(definitions), 1)
    gene_ids = definition_table.gene.to_numpy()
    ok = gene_ids[aa] != gene_ids[bb]
    aa, bb = aa[ok], bb[ok]
    and_counts = np.stack([z[aa, bb] for z in intersections], axis=1)
    del b, intersections
    rows, registry = [], {}
    for logic in ["AND", "OR"]:
        removed = and_counts if logic=="AND" else single[aa]+single[bb]-and_counts
        d = metrics(removed, totals)
        write_metric_group(h5, "all_pairs_" + logic.lower(), d, threshold_a_index=aa, threshold_b_index=bb)
        for i in positions_to_keep(d, top=3000):
            rule = dict(action="remove", logic=logic, rules=[definitions[aa[i]], definitions[bb[i]]])
            row = metric_row(rule, removed[i], totals, "python_full_surfaceome_pair_grid")
            registry[row["gate_id"]] = rule
            rows.append(row)
    return rows, registry, 2*len(aa)


def refine_pairs(cells, candidates, registry, totals, state_index, h5):
    selected = []
    seen = set()
    for row in sorted(candidates, key=lambda x:x["practical_score"], reverse=True):
        rule = registry[row["gate_id"]]
        key = tuple(sorted(x["gene"] for x in rule["rules"]))
        if len(key)==2 and key not in seen:
            seen.add(key)
            selected.append(key)
        if len(selected)>=12:
            break
    # ERBB4/CXCR4 is always a benchmark, independent of discovery rankings.
    if ("CXCR4", "ERBB4") not in selected:
        selected.append(("CXCR4", "ERBB4"))
    rows, rules, all_counts, definitions = [], {}, [], []
    for genes in selected:
        grids = []
        for gene in genes:
            x = cells[gene].to_numpy(float)
            positive = x[x>0]
            grids.append(np.unique(np.r_[0., np.quantile(positive, np.linspace(0, .98, 26))]))
        for t1, t2 in itertools.product(*grids):
            m1 = cells[genes[0]].to_numpy(float)>t1
            m2 = cells[genes[1]].to_numpy(float)>t2
            for logic in ["AND", "OR"]:
                removed = m1&m2 if logic=="AND" else m1|m2
                counts = np.bincount(state_index[removed], minlength=4)
                rule = dict(action="remove", logic=logic, rules=[dict(gene=genes[0],op=">",threshold=float(t1)),dict(gene=genes[1],op=">",threshold=float(t2))])
                row = metric_row(rule, counts, totals, "python_refined_pair_grid")
                rows.append(row)
                rules[row["gate_id"]] = rule
                all_counts.append(counts)
                definitions.append(rule)
    d = metrics(all_counts, totals)
    write_metric_group(h5, "refined_pairs", d)
    h5["refined_pairs"].create_dataset("rules_json", data=np.asarray([json.dumps(x) for x in definitions], dtype=h5py.string_dtype()))
    return rows, rules, len(rows)


def probe_third_marker(cells, candidates, registry, totals, state_index, h5, out):
    """Bounded complexity check; a third marker is retained only for a large gain."""
    definitions = pd.read_csv(out / "tables/pair_threshold_grid.tsv",sep="\t").to_dict("records")
    masks = np.asarray([cells[d["gene"]].to_numpy(float)>d["threshold"] for d in definitions])
    ranked = sorted(candidates,key=lambda x:x["practical_score"],reverse=True)
    parents, seen = [], set()
    for row in ranked:
        rule = registry[row["gate_id"]]
        signature = (tuple(sorted(x["gene"] for x in rule["rules"])),rule["logic"])
        if len(signature[0])==2 and signature not in seen:
            seen.add(signature)
            parents.append(rule)
        if len(parents)>=8:
            break
    allcounts, allrules = [], []
    for parent in parents:
        pairmask = apply_rule(cells,parent)
        used = {x["gene"] for x in parent["rules"]}
        ix = [i for i,x in enumerate(definitions) if x["gene"] not in used]
        combined = masks[ix]&pairmask if parent["logic"]=="AND" else masks[ix]|pairmask
        counts = np.stack([combined[:,state_index==j].sum(axis=1) for j in range(4)],axis=1)
        allcounts.append(counts)
        allrules.extend(dict(action="remove",logic=parent["logic"],rules=parent["rules"]+[definitions[i]]) for i in ix)
    removed = np.concatenate(allcounts)
    d = metrics(removed,totals)
    write_metric_group(h5,"bounded_third_marker_probe",d)
    h5["bounded_third_marker_probe"].create_dataset("rules_json",data=np.asarray([json.dumps(x) for x in allrules],dtype=h5py.string_dtype()))
    keep = positions_to_keep(d,top=500)
    rows, rules = [], {}
    for i in keep:
        row = metric_row(allrules[i],removed[i],totals,"python_bounded_third_marker_probe")
        rows.append(row)
        rules[row["gate_id"]] = allrules[i]
    return rows,rules,len(allrules)


def run_r_hypergate(cells, genes, screen, totals, state_index, out):
    """Use the existing unmodified R runner on NEW labels and NEW outputs."""
    targets = [("PV-biased", "state", "PV-biased", "Q1"),
               ("SST-biased", "state", "SST-biased", "Q2"),
               ("hybrid", "state", "PV/SST hybrid", "Q3"),
               ("PV_vs_SST", "pv_vs_sst", "PV-biased", "Q4"),
               ("retain_PV_hybrid", "retention_label", "retain", "Q5")]
    data = cells[["cell_id", "state"]+genes].copy()
    data["pv_vs_sst"] = np.where(cells.state.isin([STATES[0], STATES[2]]), cells.state, "ambiguous")
    data["retention_label"] = np.where(cells.state.isin(STATES[:2]), "retain", np.where(cells.state==STATES[2], "remove", "ambiguous"))
    data.to_csv(out / "hypergate_r_input.tsv.gz", sep="\t", index=False)
    jobs, panels = [], {}
    for name, label, target, question in targets:
        if question=="Q5":
            ranked = screen[screen.question=="Q2"].sort_values("F1", ascending=False)
        else:
            ranked = screen[screen.question==question].sort_values("F1", ascending=False)
        panel = list(dict.fromkeys(ranked.gene))[:6]
        # Reference features remain explicitly available in every bounded panel.
        panel = list(dict.fromkeys(panel+["FAT3", "PTPRM"]))
        panel = [g for g in panel if g in genes]
        panels[name] = panel
        for n in [1, 2]:
            for features in itertools.combinations(panel, n):
                for beta in [.5, 1., 2.]:
                    jobs.append(dict(job_id=f"phase2_{len(jobs):05d}",definition=name,label_column=label,target=target,beta=beta,
                                     features=";".join(features),requested_markers=n,stage="four_state_bounded_actual_R_hypergate"))
    job_table = pd.DataFrame(jobs)
    job_table.to_csv(out / "hypergate_r_jobs.tsv", sep="\t", index=False)
    fitpath = out / "hypergate_r_fits.tsv"
    # Checkpoints are reusable only if every exact requested job is present.
    if fitpath.exists():
        existing = pd.read_csv(fitpath, sep="\t")
        if existing.job_id.tolist() != job_table.job_id.tolist():
            raise RuntimeError("Incomplete or incompatible phase-two R checkpoint; preserve it and inspect before resuming.")
    else:
        env = dict(os.environ, PROJECT_ROOT=str(PROJECT_ROOT), OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2")
        with (out / "logs/hypergate_r.log").open("w") as log:
            subprocess.run([R_BIN, str(REPO_ROOT / "python_notebooks/scripts/run_div90_hypergate.R"),
                            str(out / "hypergate_r_input.tsv.gz"), str(out / "hypergate_r_jobs.tsv"),str(fitpath)],
                           stdout=log,stderr=subprocess.STDOUT,check=True,env=env)
    fits = pd.read_csv(fitpath, sep="\t")
    if len(fits)!=len(jobs):
        raise AssertionError("R checkpoint job count mismatch")
    fitted, registry, validated = [], {}, 0
    for fit in fits.itertuples():
        if isinstance(fit.error, str) and fit.error:
            continue
        conditions = json.loads(fit.rules)
        if isinstance(conditions, dict):
            conditions = [conditions]
        # R returns a captured target region. SST capture is removed, while all
        # other targets are retained; all four states are still evaluated.
        rule = dict(action="remove" if fit.target=="SST-biased" else "retain",logic="AND",rules=conditions)
        removed = apply_rule(cells, rule)
        capture = removed if rule["action"]=="remove" else ~removed
        label = data[fit.label_column].to_numpy()
        keep = label!="ambiguous"
        y = label==fit.target
        observed = (int((capture & y & keep).sum()),int((capture & ~y & keep).sum()),
                    int((~capture & y & keep).sum()),int((~capture & ~y & keep).sum()))
        expected = tuple(int(getattr(fit,k)) for k in ["tp","fp","fn","tn"])
        if observed != expected:
            raise AssertionError(f"Exact Python/R confusion mismatch {fit.job_id}: {observed} != {expected}")
        validated += 1
        row = metric_row(rule, np.bincount(state_index[removed],minlength=4),totals,"actual_R_hypergate")
        row.update(job_id=fit.job_id,hypergate_definition=fit.definition,hypergate_target=fit.target,beta=fit.beta)
        fitted.append(row)
        registry[row["gate_id"]] = rule
    pd.DataFrame(fitted).to_csv(out / "tables/hypergate_four_state_evaluations.tsv",sep="\t",index=False)
    dump_json(out / "hypergate_r_audit.json",dict(package="hypergate",panels=panels,requested_fits=len(jobs),successful_fits=validated,
               failed_fits=int(len(jobs)-validated),independently_validated_confusion_counts=validated,
               search_scope="All singles/pairs within target-specific 6-to-8-marker panels selected from the full signed univariate surface screen; beta 0.5, 1, 2. Not exhaustive all-surfaceome Hypergate pairs.",
               Q4="Only PV-biased versus SST-biased cells enter fitting. Hybrid and unresolved are excluded from fitting, retained for four-state gate evaluation.",
               Q5="PV-biased plus hybrid versus SST-biased; unresolved excluded from fitting, retained for four-state evaluation."))
    return fitted, registry, validated


def benchmark_rules(cells, totals, state_index):
    rows, registry = [], {}
    rules = [("baseline",dict(action="remove",logic="AND",rules=[]))]
    saved = json.loads((PHASE1 / "best_gate_rules.json").read_text())["best_shared_pair"]
    exact = json.loads((PHASE1 / "gate_rules.json").read_text())[saved["gate_id"]]
    rules.append(("phase1_FAT3_PTPRM_reference",dict(action="retain",logic="AND",rules=exact)))
    old = pd.read_csv(PHASE1 / "benchmark_best_by_hypothesis.tsv",sep="\t")
    er = old[(old.stage=="ERBB4_high_AND_CXCR4_high") & old.target.isin(["PV-like","non-SST/PV-candidate"])]
    if er.empty:
        er = old[old.stage.str.contains("ERBB4.*CXCR4") & old.target.isin(["PV-like","non-SST/PV-candidate"])]
    for rr in er.itertuples():
        rules.append(("phase1_ERBB4_CXCR4_"+rr.definition,dict(action="retain",logic="AND",rules=json.loads(rr.rules))))
    for source, rule in rules:
        removed = apply_rule(cells, rule)
        if source=="phase1_FAT3_PTPRM_reference" and int((~removed).sum())!=3168:
            raise AssertionError("Frozen FAT3/PTPRM capture is not exactly 3168 cells")
        row = metric_row(rule,np.bincount(state_index[removed],minlength=4),totals,source)
        rows.append(row)
        registry[row["gate_id"]] = rule
    return rows, registry


def fixed_gate_groups(cells, keyrows, registry, out):
    rows, score_rows = [], []
    for name, row in keyrows.items():
        removed = apply_rule(cells,registry[row["gate_id"]])
        for groupby in ["sample","cell_line","condition"]:
            if groupby not in cells:
                continue
            for value, indices in cells.groupby(groupby,dropna=False,sort=False).indices.items():
                state = cells.iloc[indices].state
                totals = np.asarray([(state==s).sum() for s in STATES])
                rc = np.asarray([((state.to_numpy()==s)&removed[indices]).sum() for s in STATES])
                v = {k:float(a[0]) for k,a in metrics(rc,totals).items()}
                v.update(gate_name=name,gate_id=row["gate_id"],groupby=groupby,group=str(value),baseline_n=len(indices))
                for j,s in enumerate(SHORT):
                    v[f"baseline_{s}_n"] = int(totals[j])
                    v[f"baseline_{s}_fraction"] = totals[j]/len(indices)
                rows.append(v)
        for fraction, mask in [("starting",np.ones(len(cells),bool)),("removed",removed),("retained",~removed)]:
            for score in ["pv_score","sst_score"]:
                vals = cells.loc[mask,score]
                score_rows.append(dict(gate_name=name,gate_id=row["gate_id"],fraction=fraction,score=score,n=len(vals),
                                       mean=vals.mean(),median=vals.median(),q25=vals.quantile(.25),q75=vals.quantile(.75)))
    pd.DataFrame(rows).to_csv(out / "tables/gate_per_sample_condition_descriptive.tsv",sep="\t",index=False)
    pd.DataFrame(score_rows).to_csv(out / "tables/gate_continuous_score_distributions.tsv",sep="\t",index=False)


def add_strong_sst_metrics(cells, row, registry):
    if "strong_sst_biased" not in cells:
        return row
    strong = cells.strong_sst_biased.astype(str).str.lower().isin(["true","1"]).to_numpy()
    removed = apply_rule(cells,registry[row["gate_id"]])
    rn, kn = int((strong&removed).sum()), int((strong&~removed).sum())
    row.update(strong_sst_total_n=int(strong.sum()),removed_strong_sst_n=rn,retained_strong_sst_n=kn,
               strong_sst_removal_recovery=safe_div(rn,strong.sum()).item(),
               removed_strong_sst_purity=safe_div(rn,removed.sum()).item(),
               strong_sst_contamination=safe_div(kn,(~removed).sum()).item())
    return row


def run(cells: pd.DataFrame, out: Path | str):
    out = Path(out).resolve()
    if out != EXPECTED_OUT.resolve():
        raise ValueError(f"Phase-two gate runtime must be {EXPECTED_OUT}")
    for directory in ["tables","logs"]:
        (out / directory).mkdir(parents=True,exist_ok=True)
    if len(cells)!=4768 or cells.cell_id.duplicated().any():
        raise AssertionError("Expected 4,768 distinct frozen entry cells")
    annotations = pd.read_csv(PHASE1 / "allowed_surface_markers.tsv",sep="\t")
    modules = pd.read_csv(out / "module_genes.tsv",sep="\t")
    included = modules.included.astype(str).str.lower().isin(["true","1"])
    excluded = set(modules.loc[included,"gene"])
    genes = [g for g in annotations.gene if g in cells and g not in excluded]
    if not genes or set(genes)&excluded:
        raise AssertionError("No allowed features, or developmental-label leakage")
    if not np.isfinite(cells[genes].to_numpy(float)).all():
        raise AssertionError("Nonfinite surface expression")
    state_index = pd.Categorical(cells.state,categories=STATES).codes
    if (state_index<0).any():
        raise AssertionError("Unrecognized state label")
    totals = np.bincount(state_index,minlength=4)
    if (totals==0).any():
        raise AssertionError("All four operational states must be represented")
    start = time.time()
    print(f"Phase2 surface discovery: {len(genes)} genes, state counts {dict(zip(STATES,totals))}",flush=True)
    effects, screen = marker_screen(cells,genes,annotations,out)
    print("Gene effects and exact signed state screens complete",flush=True)
    rows, registry = benchmark_rules(cells,totals,state_index)
    h5path = out / "tables/depletion_all_candidate_metrics.h5"
    if h5path.exists():
        raise FileExistsError(f"Preserve existing gate search before resuming: {h5path}")
    # The mounted Turbo filesystem can fail HDF5's read-after-write random I/O.
    # Build this <=1GB scientific table in RAM and publish one sequential image.
    # No runtime file is created in the checkout or a local scratch directory.
    with h5py.File("phase2_depletion_metrics_in_memory","w",driver="core",backing_store=False) as h5:
        h5.attrs["genes"] = json.dumps(genes)
        h5.attrs["states"] = json.dumps(STATES)
        h5.attrs["baseline_counts"] = json.dumps(totals.tolist())
        h5.attrs["description"] = "Every Python depletion candidate: exact RNA > thresholds; metric columns in each group's metric_columns attribute; pair threshold index references root pair_threshold_definitions."
        one, rr, nsingle = exact_single_screen(cells,genes,totals,state_index,h5,out)
        rows.extend(one); registry.update(rr)
        print(f"Exact positive one-marker sweep complete: {nsingle} thresholds",flush=True)
        pair, rr, npair = broad_pair_screen(cells,genes,totals,state_index,h5,out)
        rows.extend(pair); registry.update(rr)
        print(f"Broad pair sweep complete: {npair} AND/OR candidates",flush=True)
        refined, rr, nrefined = refine_pairs(cells,pair,registry,totals,state_index,h5)
        rows.extend(refined); registry.update(rr)
        print(f"Dense refinement complete: {nrefined} candidates",flush=True)
        triples, rr, ntriples = probe_third_marker(cells,pair+refined,registry,totals,state_index,h5,out)
        rows.extend(triples); registry.update(rr)
        print(f"Bounded third-marker improvement probe complete: {ntriples} candidates",flush=True)
        h5.flush()
        h5path.write_bytes(h5.id.get_file_image())
    fitted, rr, n_r = run_r_hypergate(cells,genes,screen,totals,state_index,out)
    rows.extend(fitted); registry.update(rr)
    print(f"Actual R Hypergate complete: {n_r} fits validated",flush=True)
    allrows = pd.DataFrame(rows).drop_duplicates("gate_id").sort_values("practical_score",ascending=False)
    # Discovery priority requires a positive SST-removal rule; arbitrary R
    # enrichment rectangles remain separately available for Questions 1-5.
    eligible = allrows[(allrows.action=="remove") & (allrows.n_markers>0)].copy()
    eligible = eligible[eligible.gate_id.map(lambda gid:all(x["op"] in [">",">="] for x in registry[gid]["rules"]))]
    best_single = eligible[eligible.n_markers==1].iloc[0].to_dict()
    best_pair = eligible[eligible.n_markers==2].iloc[0].to_dict()
    # A second marker must improve the balanced score by >= 0.02 absolute;
    # neither PV nor hybrid recovery may fall by more than 0.03 absolute.
    meaningful = (best_pair["practical_score"]-best_single["practical_score"]>=.02 and
                  best_pair["pv_recovery"]>=best_single["pv_recovery"]-.03 and
                  best_pair["hybrid_recovery"]>=best_single["hybrid_recovery"]-.03)
    selected = best_pair if meaningful else best_single
    best_three = eligible[eligible.n_markers==3].iloc[0].to_dict()
    meaningful_three = (best_three["practical_score"]-best_pair["practical_score"]>=.03 and
                        best_three["pv_recovery"]>=best_pair["pv_recovery"]-.03 and
                        best_three["hybrid_recovery"]>=best_pair["hybrid_recovery"]-.03)
    if meaningful_three and best_three["practical_score"]-selected["practical_score"]>=.03:
        selected = best_three
    reference = allrows[allrows.source=="phase1_FAT3_PTPRM_reference"].iloc[0].to_dict()
    baseline = allrows[allrows.source=="baseline"].iloc[0].to_dict()
    erbb = allrows[allrows.gate_id.map(lambda gid:set(x["gene"] for x in registry[gid]["rules"])=={"ERBB4","CXCR4"})]
    erbb = erbb.iloc[0].to_dict()
    old_erbb = [x for x in rows if x["source"].startswith("phase1_ERBB4_CXCR4")]
    tradeoffs = []
    for minimum in [.1,.2,.3,.5]:
        for nmarkers in [1,2,3]:
            possible = eligible[(eligible.n_markers==nmarkers)&(eligible.sst_removal_recovery>=minimum)]
            if not possible.empty:
                row = possible.iloc[0].to_dict()
                row.update(minimum_sst_removal=minimum,analyst_chosen_constraint=True)
                tradeoffs.append(add_strong_sst_metrics(cells,row,registry))
    practical_single = next(row.copy() for row in tradeoffs if row["minimum_sst_removal"]==.2 and row["n_markers"]==1)
    practical_pair = next(row.copy() for row in tradeoffs if row["minimum_sst_removal"]==.2 and row["n_markers"]==2)
    practical_pair_meaningful = (practical_pair["practical_score"]-practical_single["practical_score"]>=.02 and
                                practical_pair["pv_recovery"]>=practical_single["pv_recovery"]-.03 and
                                practical_pair["hybrid_recovery"]>=practical_single["hybrid_recovery"]-.03)
    experimental = practical_pair if practical_pair_meaningful else practical_single
    keyrows = dict(baseline=baseline,best_single=best_single,best_pair=best_pair,best_three=best_three,selected=selected,
                   phase1_reference=reference,erbb4_cxcr4=erbb,practical_single=practical_single,practical_pair=practical_pair,
                   experimental_depletion=experimental)
    keyrows = {name:add_strong_sst_metrics(cells,row,registry) for name,row in keyrows.items()}
    pd.DataFrame(tradeoffs).to_csv(out / "tables/depletion_minimum_removal_sensitivity.tsv",sep="\t",index=False)
    three_marker = dict(run=True,candidates=ntriples,meaningful_improvement=bool(meaningful_three),
                        score_gain_over_pair=best_three["practical_score"]-best_pair["practical_score"],
                        scope="Bounded improvement probe: eight best distinct gene-pair/logic combinations, plus every remaining surface gene on the original positive-threshold grid. Homogeneous AND or OR only. Not exhaustive triplets.",
                        selection_criterion="At least 0.03 absolute primary-score gain over best pair, with neither PV nor hybrid recovery dropping more than 0.03 absolute.")
    global_d = {k:allrows[k].to_numpy() for k in ["target_recovery","sst_contamination"]}
    frontier = allrows.iloc[pareto_positions(global_d)].copy()
    frontier["is_pareto"] = True
    frontier.to_csv(out / "tables/gate_pareto.tsv",sep="\t",index=False)
    shortlist = pd.concat([eligible.head(100),pd.DataFrame(list(keyrows.values())),frontier],ignore_index=True).drop_duplicates("gate_id")
    shortlist = pd.DataFrame([add_strong_sst_metrics(cells,row,registry) for row in shortlist.to_dict("records")])
    shortlist["is_pareto"] = shortlist.gate_id.isin(frontier.gate_id)
    shortlist.to_csv(out / "tables/gate_shortlist.tsv",sep="\t",index=False)
    allrows.to_csv(out / "tables/gate_ranked_candidates.tsv.gz",sep="\t",index=False,compression={"method":"gzip","compresslevel":3})
    assignment = cells[["cell_id"]].copy()
    for name,row in keyrows.items():
        if name=="baseline":
            continue
        name = "phase1" if name=="phase1_reference" else name
        assignment[name+"_retained"] = ~apply_rule(cells,registry[row["gate_id"]])
    assignment.to_csv(out / "gate_cell_assignments.tsv.gz",sep="\t",index=False)
    fixed_gate_groups(cells,keyrows,registry,out)
    # Independently verify every material rule against stored counts.
    verified = 0
    for row in shortlist.to_dict("records"):
        removed = apply_rule(cells,registry[row["gate_id"]])
        observed = np.bincount(state_index[removed],minlength=4)
        expected = np.asarray([row[f"removed_{s}_n"] for s in SHORT],dtype=int)
        if not np.array_equal(observed,expected):
            raise AssertionError(f"Candidate count audit failed: {row['gate_id']}")
        verified += 1
    summary = dict(**keyrows,rules=registry,phase1_erbb4_cxcr4_benchmarks=old_erbb,
                   total_surface_features=len(genes),excluded_module_features=sorted(excluded),
                   baseline_state_counts=dict(zip(STATES,totals)),exact_single_candidates=nsingle,
                   broad_pair_candidates=npair,refined_pair_candidates=nrefined,actual_R_hypergate_fits=n_r,
                   bounded_three_marker_candidates=ntriples,
                   shortlist_rules_independently_verified=verified,
                   score_definition="sqrt(PV-biased recovery * hybrid recovery) * (1 - retained SST-biased fraction among retained PV-biased + hybrid + SST-biased cells). This excludes unresolved from purity only, without discarding their cell identities. practical_score_unresolved_retained uses all retained cells in the purity denominator.",
                   unresolved_sensitivity="Primary retained fractions include unresolved; *_unresolved_excluded metrics describe separately excluding unresolved. All original cell IDs and scores remain available.",
                   pair_is_meaningful=bool(meaningful),pair_score_gain=best_pair["practical_score"]-best_single["practical_score"],
                   simple_gate_rule="Prefer one marker unless best pair improves primary score by >=0.02 absolute, with PV and hybrid recoveries each falling no more than 0.03 absolute.",
                   three_marker=three_marker,
                   practical_depletion_constraint="For an informative prospective depletion comparison, an analyst-chosen minimum of 20% SST-biased removal is evaluated separately from the unconstrained score optimum. Sensitivity at 10%, 20%, 30%, and 50% removal is reported. This is a decision constraint, not a fitted biological boundary or a requirement supplied by the user.",
                   practical_pair_is_meaningful=bool(practical_pair_meaningful),
                   threshold_search="Python exact observed-value upper-threshold sweep for every one marker. Every pair of allowed genes with AND and OR of positive thresholds: 0 and 25th/50th/75th/90th percentiles of positive expression. Top 12 distinct pairs and ERBB4/CXCR4 then refined on a 27-point threshold grid. Actual R Hypergate results are separate.",
                   complete_metrics_file=str(h5path),
                   limitations=["All results are in-sample descriptive discovery, not held-out validation.","Per-sample fixed global gate results describe biological/sample composition and are not cross-validation.","RNA thresholds do not establish extracellular protein expression or fluorescence thresholds.","Unresolved exclusion is a labeling sensitivity, not an additional demonstrated surface sorting gate.","The four labels are operational developmental states, not established fates."],
                   elapsed_seconds=time.time()-start)
    dump_json(out / "gate_summary.json",summary)
    (out / "tables/README_gate_metrics.md").write_text(
        "# Phase-two gate metrics\n\nEvery Python one-marker and pair candidate is in `depletion_all_candidate_metrics.h5`. "
        "Each group has a `metrics` matrix and a JSON `metric_columns` attribute. Exact singles refer to the root `genes` attribute using `gene_index` plus a `threshold` vector (remove expression > threshold). "
        "Broad pairs refer to the root JSON `pair_threshold_definitions` attribute using `threshold_a_index`/`threshold_b_index`; group names specify AND or OR. "
        "Refined pairs contain aligned `rules_json`. Counts are exact; metrics are float64.\n\n"
        "`depletion_single_all.tsv.gz` is the complete human-readable exact one-marker sweep. `gate_ranked_candidates.tsv.gz` contains top candidates under both unresolved treatments and all potentially nondominated candidates from each exhaustive search stage, plus all refined candidates, saved references, and successful R gates. "
        "`gate_pareto.tsv` is the pooled frontier. `gate_shortlist.tsv` is for figures. "
        "Gate rules in `gate_summary.json` state whether the captured region is removed or retained. All RNA thresholds use the frozen log1p(CP10K) matrix.\n\n"
        +summary["score_definition"]+"\n\n"+summary["threshold_search"]+"\n")
    return clean_json(summary)


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--out",type=Path,default=EXPECTED_OUT)
    args = p.parse_args()
    cells = pd.read_csv(args.out / "cells.tsv.gz",sep="\t")
    run(cells,args.out)
