#!/usr/bin/env python3
"""Independent DIV30 surface discovery, confound diagnostics, and figures.

The caller supplies barcode-aligned DIV30 expression and scores. No DIV90
gate is used to select DIV30 features or thresholds. Exact RNA cut points
are descriptive discovery thresholds, never fluorescence gate settings.
"""
from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder')
OLD = ROOT / 'results/div90_hypergate_sst_pv'
DIV90 = ROOT / 'results/div90_hypergate_sst_pv_phase2'
FLOORS = (.1, .2, .3, .5, .7, .8, .9)
DIAGNOSTIC = ['FGFR2', 'PTPRS', 'FAT3', 'PTPRM', 'ERBB4', 'CXCR4', 'ACKR3', 'NRP1', 'NRP2', 'PLXNA2']
HOUSEKEEPING = {'B2M', 'TFRC', 'SLC2A1', 'ATP1A1', 'ATP1B1', 'ATP1B3', 'SLC3A2',
                'SLC7A5', 'LAMP1', 'LAMP2', 'HSPA5', 'HSP90B1', 'P4HB'}


def clean(value):
    if isinstance(value, dict): return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)): return [clean(v) for v in value]
    if isinstance(value, np.ndarray): return clean(value.tolist())
    if isinstance(value, np.generic): return clean(value.item())
    if isinstance(value, float) and not np.isfinite(value): return None
    return value


def write_json(path, value):
    path.write_text(json.dumps(clean(value), indent=2, allow_nan=False) + '\n')


def enrichment_bin(x):
    return 'negligible' if x < 1.2 else 'modest' if x < 1.5 else 'potentially meaningful' if x < 2 else 'strong'


def metric(tp, retained, n_target, n_total):
    tp, retained = np.asarray(tp, float), np.asarray(retained, float)
    purity = np.divide(tp, retained, out=np.full_like(tp, np.nan), where=retained > 0)
    baseline = n_target / n_total if n_total else np.nan
    return dict(starting_target_fraction=np.full_like(tp, baseline), post_gate_target_fraction=purity,
                fold_enrichment=purity / baseline if baseline else np.full_like(tp, np.nan),
                target_recovery=tp / n_target if n_target else np.full_like(tp, np.nan),
                total_cell_yield=retained / n_total if n_total else np.full_like(tp, np.nan),
                target_retained_n=tp, retained_n=retained,
                composition_change_pp=100 * (purity - baseline))


def gate(clauses=(), logic='AND', source='independent_DIV30_search'):
    clauses = [dict(gene=g, op=op, threshold=float(t)) for g, op, t in clauses]
    obj = dict(clauses=clauses, logic=logic)
    ident = hashlib.sha256(json.dumps(obj, sort_keys=True).encode()).hexdigest()[:16]
    return dict(gate_id=ident, gate_label=(' ' + logic + ' ').join(f'{c["gene"]} {c["op"]} {c["threshold"]:.9g}' for c in clauses) or 'retain all DIV30 cells',
                n_markers=len({c['gene'] for c in clauses}), rules=json.dumps(clauses), logic=logic, source=source)


def apply(row, expression, lookup):
    clauses = json.loads(row['rules'])
    # NumPy 2 treats a Python float as a weak scalar and can round the cut to
    # float32. A strong float64 scalar preserves the archived RNA threshold.
    masks = [(expression[:, lookup[c['gene']]] > np.float64(c['threshold'])) if c['op'] == '>' else
             (expression[:, lookup[c['gene']]] <= np.float64(c['threshold'])) for c in clauses]
    if not masks: return np.ones(expression.shape[0], bool)
    return (np.logical_or.reduce(masks) if row['logic'] == 'OR' else np.logical_and.reduce(masks))


def add_row(spec, tp, retained, y, target):
    row = {**spec, **{k: float(v) for k, v in metric(tp, retained, int(y.sum()), len(y)).items()}, 'target': target}
    row['enrichment_interpretation'] = enrichment_bin(row['fold_enrichment'])
    row['evaluation'] = 'in-sample descriptive discovery'
    return row


def selected_positions(tp, retained, y):
    m = metric(tp, retained, int(y.sum()), len(y))
    recovery, purity = m['target_recovery'], m['post_gate_target_fraction']
    chosen = []
    for floor in FLOORS:
        allowed = np.flatnonzero((recovery >= floor) & np.isfinite(purity))
        if len(allowed):
            order = np.lexsort((-recovery[allowed], -purity[allowed]))
            chosen.append(int(allowed[order[0]]))
    # Full nondominated frontier, with deterministic tie handling.
    ix = np.flatnonzero(np.isfinite(purity))
    ix = ix[np.lexsort((-purity[ix], -recovery[ix]))]
    previous_best = np.r_[-np.inf, np.maximum.accumulate(purity[ix])[:-1]]
    frontier = ix[purity[ix] > previous_best + 1e-12].tolist()
    # Preserve the complete frontier in count arrays; compact searchable registry.
    if len(frontier) > 40:
        frontier = [frontier[i] for i in np.linspace(0, len(frontier)-1, 40).astype(int)]
    return np.unique(chosen + frontier)


def annotation_catalog(expression, genes, excluded_genes, out):
    cspa = pd.read_csv(OLD / 'provenance/cspa_annotation_audit.tsv', sep='\t')
    tm = pd.to_numeric(cspa['Phobius TM predicted yes/no'], errors='coerce').eq(1)
    gpi = pd.to_numeric(cspa['GPI'], errors='coerce').eq(1)
    cspa = cspa.loc[tm | gpi].dropna(subset=['gene']).drop_duplicates('gene')
    ann = cspa.set_index('gene').to_dict('index')
    curated = pd.read_csv(OLD / 'provenance/surface_feature_filter_audit.tsv', sep='\t')
    # Prior expression exclusions are discarded; only biological annotation enters.
    for row in curated.to_dict('records'):
        if row['gene'] not in ann and (row['source'] == 'user_core_panel' or bool(row.get('allowed', False))):
            ann[row['gene']] = {'CSPA category': 'existing curated surface annotation', 'UP_Protein_name': row.get('protein_name', '')}
    lookup, excluded = {g:i for i,g in enumerate(genes)}, set(excluded_genes)
    rows = []
    for g, item in sorted(ann.items()):
        reason, detection, sd = '', np.nan, np.nan
        if g not in lookup: reason = 'not matched in DIV30 input'
        else:
            x = expression[:, lookup[g]]
            detection, sd = float((x > 0).mean()), float(np.std(x))
            if g in excluded: reason = 'score_mapping_or_supervised_feature_excluded_to_prevent_circularity'
            elif g in HOUSEKEEPING or g.startswith(('RPL', 'RPS', 'MT-', 'HSP')): reason = 'housekeeping_or_intracellular_ambiguity'
            elif g in {'BCAN', 'FRAS1'}: reason = 'extracellular_matrix_or_isoform_ambiguity'
            elif detection < .02 or sd < 1e-7: reason = 'DIV30_detection_below_2pct_or_no_variation'
        rows.append(dict(gene=g, allowed=not bool(reason), exclusion_reason=reason, DIV30_detection_fraction=detection,
                         DIV30_log_expression_sd=sd, annotation_category=item.get('CSPA category', ''),
                         protein_name=item.get('UP_Protein_name', ''), source='CSPA_TM_or_GPI_all_categories_union_prior_curated',
                         surface_protein_confirmation='not measured in this scRNA-seq dataset'))
    table = pd.DataFrame(rows)
    table.to_csv(out / 'tables/surface_feature_audit.tsv', sep='\t', index=False)
    return table, table.loc[table.allowed, 'gene'].tolist()


def targets(cells, sst_stable):
    score = cells.early_pv_score.to_numpy(float)
    out = {f'top{int(q*100)}': score >= np.quantile(score, 1-q) for q in (.1, .2, .3)}
    out['top20_strong_MGE'] = out['top20'] & (cells.mge_score.to_numpy(float) >= cells.mge_score.median())
    out['top20_postmitotic'] = out['top20'] & cells.postmitotic.to_numpy(bool)
    if sst_stable:
        out['top20_lowSST'] = out['top20'] & (cells.direct_sst_score.to_numpy(float) < cells.direct_sst_score.median())
    return {k:v for k,v in out.items() if 20 <= v.sum() < len(v)}


def exact_singles(expression, lookup, feature_genes, y, target, h5=None, purity_by_target_count=None):
    rows = [add_row(gate(), int(y.sum()), len(y), y, target)]
    hg = h5.create_group('exact_single') if h5 is not None else None
    for g in feature_genes:
        x = expression[:, lookup[g]]
        order = np.argsort(x, kind='stable')
        thresholds, starts, counts = np.unique(x[order], return_index=True, return_counts=True)
        low_n = np.cumsum(counts); low_tp = np.cumsum(np.add.reduceat(y[order].astype(np.int32), starts))
        if hg is not None:
            gg = hg.create_group(g)
            for k,v in [('threshold', thresholds.astype(np.float64)), ('retained_n_le', low_n.astype(np.int32)), ('target_n_le', low_tp.astype(np.int32))]:
                gg.create_dataset(k, data=v, compression='gzip', compression_opts=3)
        for op, tp, n in [('<=', low_tp, low_n), ('>', y.sum()-low_tp, len(y)-low_n)]:
            if purity_by_target_count is not None:
                valid=n>0
                np.maximum.at(purity_by_target_count,tp[valid].astype(int),tp[valid]/n[valid])
            for i in selected_positions(tp, n, y):
                rows.append(add_row(gate([(g, op, thresholds[i])]), int(tp[i]), int(n[i]), y, target))
    return rows


def shortlist(rows, maximum=18):
    df = pd.DataFrame(rows)
    df = df[df.n_markers.eq(1)]
    genes = []
    for floor in (.8, .5, .9, .7, .3):
        d = df[df.target_recovery >= floor].sort_values(['post_gate_target_fraction', 'target_recovery'], ascending=False)
        for rules in d.rules:
            g = json.loads(rules)[0]['gene']
            if g not in genes:
                genes.append(g)
                if len(genes) % max(1, maximum//5) == 0: break
        if len(genes) >= maximum: break
    if len(genes) < maximum:
        for rules in df.sort_values('post_gate_target_fraction', ascending=False).rules:
            g = json.loads(rules)[0]['gene']
            if g not in genes: genes.append(g)
            if len(genes) >= maximum: break
    return genes[:maximum]


def grid(x):
    positive = x[x > 0]
    return np.unique(np.r_[0., np.quantile(positive, [.05, .1, .2, .35, .5, .65, .8, .9, .95])]) if len(positive) else np.array([0.])


def pair_search(expression, lookup, candidate_genes, y, target, h5=None, purity_by_target_count=None):
    rows = []; hg = h5.create_group('screened_pairs') if h5 is not None else None
    grids = {g:grid(expression[:, lookup[g]]) for g in candidate_genes}
    bins = {g:np.searchsorted(grids[g], expression[:, lookup[g]], side='left') for g in candidate_genes}
    for ga,gb in itertools.combinations(candidate_genes, 2):
        a,b = grids[ga], grids[gb]
        shape = (len(a)+1, len(b)+1)
        flat = bins[ga]*shape[1] + bins[gb]
        all_hist = np.bincount(flat, minlength=np.prod(shape)).reshape(shape)
        target_hist = np.bincount(flat, weights=y.astype(int), minlength=np.prod(shape)).reshape(shape).astype(int)
        all_cum, target_cum = all_hist.cumsum(0).cumsum(1), target_hist.cumsum(0).cumsum(1)
        gg = hg.create_group(ga+'__'+gb) if hg is not None else None
        if gg is not None:
            gg.create_dataset('threshold_a',data=a); gg.create_dataset('threshold_b',data=b)
            gg.attrs['gene_a']=ga; gg.attrs['gene_b']=gb
        for opa,opb in itertools.product(('<=','>'),repeat=2):
            vals=[]
            for c,total in [(all_cum,len(y)), (target_cum,int(y.sum()))]:
                lower=c[:-1,:-1]; margin_a=c[:-1,-1,None]; margin_b=c[-1,None,:-1]
                aa=margin_a if opa=='<=' else total-margin_a
                bb=margin_b if opb=='<=' else total-margin_b
                intersect = lower if opa=='<=' and opb=='<=' else margin_a-lower if opa=='<=' else margin_b-lower if opb=='<=' else total-margin_a-margin_b+lower
                vals.append((intersect, aa+bb-intersect))
            for li,logic in enumerate(('AND','OR')):
                n,tp = vals[0][li].ravel(),vals[1][li].ravel()
                if purity_by_target_count is not None:
                    valid=n>0
                    np.maximum.at(purity_by_target_count,tp[valid].astype(int),tp[valid]/n[valid])
                if gg is not None:
                    name=('le' if opa=='<=' else 'gt')+'_'+('le' if opb=='<=' else 'gt')+'_'+logic
                    ds=gg.create_dataset(name,data=np.stack([n,tp],axis=1).astype(np.int32),compression='gzip',compression_opts=3)
                    ds.attrs['columns']='retained_n,target_n'; ds.attrs['reshape']=(len(a),len(b))
                for i in selected_positions(tp,n,y):
                    ia,ib=np.unravel_index(i,(len(a),len(b)))
                    rows.append(add_row(gate([(ga,opa,a[ia]),(gb,opb,b[ib])],logic),int(tp[i]),int(n[i]),y,target))
    return rows


def best_at_floors(rows):
    df=pd.DataFrame(rows); selected=[]
    for floor in FLOORS:
        for k in (1,2,3):
            d=df[(df.target_recovery >= floor) & (df.n_markers <= k)]
            if len(d):
                row=d.sort_values(['post_gate_target_fraction','n_markers','target_recovery'],ascending=[False,True,False]).iloc[0].to_dict()
                row.update(recovery_floor=floor,maximum_markers=k);selected.append(row)
    return selected


def triplet_search(expression, lookup, candidate_genes, y, target, base_rows, matched_frontier=None):
    df=pd.DataFrame(base_rows); seeds=[]
    for floor in FLOORS:
        d=df[(df.n_markers==2)&(df.target_recovery>=floor)].sort_values('post_gate_target_fraction',ascending=False)
        seeds.extend(d.head(2).to_dict('records'))
    seeds=list({r['gate_id']:r for r in seeds}.values())
    evaluated=[]; qualifying=[]
    for seed in seeds:
        present={c['gene'] for c in json.loads(seed['rules'])}; initial=apply(seed,expression,lookup)
        clauses=[(c['gene'],c['op'],c['threshold']) for c in json.loads(seed['rules'])]
        for g in candidate_genes[:12]:
            if g in present: continue
            x=expression[:,lookup[g]]
            for threshold in grid(x)[::2]:
                for op in ('<=','>'):
                    one=x<=np.float64(threshold) if op=='<=' else x>np.float64(threshold)
                    keep=initial&one if seed['logic']=='AND' else initial|one
                    tp,n=int(y[keep].sum()),int(keep.sum())
                    if not n: continue
                    row=add_row(gate(clauses+[(g,op,threshold)],seed['logic'],'conditional_DIV30_triplet_extension'),tp,n,y,target)
                    matched=float(matched_frontier[tp]) if matched_frontier is not None else df[df.target_recovery>=row['target_recovery']-1e-12].post_gate_target_fraction.max()
                    row['purity_gain_vs_1_or_2_markers_at_equal_or_higher_recovery']=row['post_gate_target_fraction']-matched
                    row['material_triplet_improvement']=bool(row['purity_gain_vs_1_or_2_markers_at_equal_or_higher_recovery']>=.02)
                    evaluated.append(row)
                    if row['material_triplet_improvement']: qualifying.append(row)
    return evaluated,qualifying


def project_residual(values, covariates):
    x=np.asarray(values,float)
    c=np.column_stack([np.ones(len(x)),np.asarray(covariates,float)])
    return x-c@np.linalg.lstsq(c,x,rcond=None)[0]


def correlation(a,b):
    a,b=np.asarray(a,float),np.asarray(b,float)
    if np.std(a)<1e-12 or np.std(b)<1e-12: return np.nan
    return float(np.corrcoef(a,b)[0,1])


def classify_association(raw, partial, neighbor_partial, sample_partial, progen, mat, mge):
    if (np.isfinite(partial) and abs(partial)>=.1 and abs(partial)>=.5*abs(raw) and
        np.isfinite(neighbor_partial) and abs(neighbor_partial)>=.05 and partial*neighbor_partial>0 and
        np.isfinite(sample_partial) and abs(sample_partial)>=.05 and partial*sample_partial>0):
        return 'early PV-associated independent of measured maturation (exploratory)'
    if progen>=.25 and abs(progen)>=abs(mat): return 'progenitor marker'
    if abs(mat)>=.25: return 'generic neuronal maturation'
    if abs(mge)>=.25: return 'MGE lineage marker'
    return 'ambiguous'


def confound_marker_table(cells, expression, lookup, genes, targets_by_name, out):
    y=cells.early_pv_score.to_numpy(float); neigh=cells.neighbor_pv_score.to_numpy(float)
    covs=cells[['progenitor_score','maturation_score','mge_score','log_counts']].to_numpy(float)
    ncovs=cells[['progenitor_score','maturation_score','log_counts']].to_numpy(float)
    dummies=pd.get_dummies(cells['sample'],drop_first=True,dtype=float).to_numpy()
    samplecov=np.column_stack([covs,dummies])
    # A projection basis avoids refitting the same nuisance design per marker.
    q=np.linalg.qr(np.column_stack([np.ones(len(y)),covs]),mode='reduced')[0]
    nq=np.linalg.qr(np.column_stack([np.ones(len(y)),ncovs]),mode='reduced')[0]
    sq=np.linalg.qr(np.column_stack([np.ones(len(y)),samplecov]),mode='reduced')[0]
    ry=y-q@(q.T@y); rn=neigh-nq@(nq.T@neigh); sy=y-sq@(sq.T@y)
    top=targets_by_name['top20']; rows=[]
    for g in genes:
        x=np.asarray(expression[:,lookup[g]],float); rx=x-q@(q.T@x); nx=x-nq@(nq.T@x); sx=x-sq@(sq.T@x)
        a,b=x[top],x[~top]; pooled=np.sqrt(((len(a)-1)*a.var(ddof=1)+(len(b)-1)*b.var(ddof=1))/(len(x)-2))
        edges=np.unique(np.r_[0.,np.quantile(x[x>0],np.linspace(0,1,21))]) if (x>0).any() else np.array([0.,1.])
        if len(edges)<2: edges=np.array([0.,max(1.,float(x.max()))])
        edges=np.r_[-np.inf,np.nextafter(0.,1.),edges[edges>0],np.inf]
        hist_a=np.histogram(a,bins=edges)[0]/len(a);hist_b=np.histogram(b,bins=edges)[0]/len(b)
        partial=correlation(rx,ry); neighbor_partial=correlation(nx,rn); sample_partial=correlation(sx,sy)
        raw=correlation(x,y); mat=correlation(x,cells.maturation_score);progen=correlation(x,cells.progenitor_score);mge=correlation(x,cells.mge_score)
        classification=classify_association(raw,partial,neighbor_partial,sample_partial,progen,mat,mge)
        df=max(1,len(x)-q.shape[1]-1)
        p=2*stats.t.sf(abs(partial)*np.sqrt(df/max(1e-15,1-partial**2)),df) if np.isfinite(partial) else np.nan
        rows.append(dict(gene=g,mean_logexpr_top20=a.mean(),mean_logexpr_other=b.mean(),mean_difference=a.mean()-b.mean(),
                         cohens_d=(a.mean()-b.mean())/pooled if pooled else np.nan,distribution_overlap=float(np.minimum(hist_a,hist_b).sum()),
                         detection_top20=float((a>0).mean()),detection_other=float((b>0).mean()),
                         pearson_early_pv=raw,spearman_early_pv=stats.spearmanr(x,y).statistic,
                         partial_r_early_pv=partial,partial_r_neighbor_pv=neighbor_partial,partial_r_early_pv_plus_sample=sample_partial,
                         standardized_adjusted_beta_early_pv=float(np.dot(rx,ry)/np.dot(rx,rx)*x.std()/y.std()) if np.dot(rx,rx)>1e-10 else np.nan,
                         standardized_adjusted_beta_neighbor_pv=float(np.dot(nx,rn)/np.dot(nx,nx)*x.std()/neigh.std()) if np.dot(nx,nx)>1e-10 else np.nan,
                         pearson_progenitor=progen,pearson_maturation=mat,pearson_MGE=mge,
                         pearson_SST=correlation(x,cells.direct_sst_score),nominal_cell_level_p=p,classification=classification,
                         inference='descriptive cell-level association; sample adjustment is sensitivity, not independent biological replication'))
    table=pd.DataFrame(rows)
    valid=table.nominal_cell_level_p.notna(); ps=table.loc[valid,'nominal_cell_level_p'].to_numpy();order=np.argsort(ps)
    adj=np.minimum.accumulate((ps[order]*len(ps)/np.arange(1,len(ps)+1))[::-1])[::-1]; result=np.empty_like(adj);result[order]=np.minimum(adj,1)
    table.loc[valid,'nominal_cell_level_BH_q']=result
    table.to_csv(out/'tables/surface_marker_confound_models.tsv',sep='\t',index=False)
    return table


def gate_confound_table(cells, expression, lookup, rows, out):
    unique={r['gate_id']:r for r in rows};result=[]
    cov=cells[['progenitor_score','maturation_score','mge_score','log_counts']].to_numpy(float)
    ncov=cells[['progenitor_score','maturation_score','log_counts']].to_numpy(float)
    samples=pd.get_dummies(cells['sample'],drop_first=True,dtype=float).to_numpy()
    q=np.linalg.qr(np.column_stack([np.ones(len(cells)),cov]),mode='reduced')[0]
    nq=np.linalg.qr(np.column_stack([np.ones(len(cells)),ncov]),mode='reduced')[0]
    sq=np.linalg.qr(np.column_stack([np.ones(len(cells)),cov,samples]),mode='reduced')[0]
    y=cells.early_pv_score.to_numpy(float);n=cells.neighbor_pv_score.to_numpy(float)
    ry=y-q@(q.T@y);rn=n-nq@(nq.T@n);sy=y-sq@(sq.T@y)
    for row in unique.values():
        keep=apply(row,expression,lookup).astype(float)
        partial=correlation(keep-q@(q.T@keep),ry);neighbor_partial=correlation(keep-nq@(nq.T@keep),rn)
        sample_partial=correlation(keep-sq@(sq.T@keep),sy)
        raw=correlation(keep,y);progen=correlation(keep,cells.progenitor_score)
        mat=correlation(keep,cells.maturation_score);mge=correlation(keep,cells.mge_score)
        result.append(dict(gate_id=row['gate_id'],gate_label=row['gate_label'],
                           partial_r_early_pv=partial,partial_r_neighbor_pv=neighbor_partial,
                           partial_r_early_pv_plus_sample=sample_partial,
                           classification=classify_association(raw,partial,neighbor_partial,sample_partial,progen,mat,mge),
                           pearson_progenitor=progen,pearson_maturation=mat,pearson_MGE=mge,
                           retained_progenitor_fraction=float(cells.progenitor[keep>0].mean()),
                           retained_postmitotic_fraction=float(cells.postmitotic[keep>0].mean()),
                           interpretation='adjusted cross-sectional association; no lineage or independent significance established'))
    pd.DataFrame(result).to_csv(out/'tables/surface_gate_confound_models.tsv',sep='\t',index=False)


def diagnostic_single_metrics(expression, lookup, target_dict, annotation, out):
    rows=[];eligible=set(annotation.loc[annotation.allowed,'gene'])
    for target,y in target_dict.items():
        candidates=exact_singles(expression,lookup,[g for g in DIAGNOSTIC if g in lookup],y,target)
        for g in DIAGNOSTIC:
            one=[r for r in candidates if r['n_markers']==1 and json.loads(r['rules'])[0]['gene']==g]
            if not one: continue
            for row in best_at_floors(one):
                if row['maximum_markers']!=1: continue
                row.update(gene=g,eligible_independent_gate_feature=g in eligible,
                           threshold_status='DIV30 descriptive single-marker optimum, separate from frozen DIV90 benchmark')
                rows.append(row)
    pd.DataFrame(rows).to_csv(out/'tables/DIV90_marker_DIV30_diagnostic_gate_metrics.tsv',sep='\t',index=False)


def held_sample_validation(cells, expression, lookup, feature_genes, y, out):
    rows=[]; pooled={floor:dict(tp=0,n=0,targets=0,total=0) for floor in FLOORS}
    sample_values=cells['sample'].astype(str).to_numpy()
    for sample in sorted(set(sample_values)):
        test=sample_values==sample;train=~test
        if y[train].sum()<20 or y[test].sum()<5: continue
        # Selection and thresholds are trained anew without this sample.
        # Only surface columns are copied from the read-only expression memmap.
        columns=[lookup[g] for g in feature_genes]
        train_x=expression[np.ix_(train,columns)];test_x=expression[np.ix_(test,columns)]
        fold_lookup={g:i for i,g in enumerate(feature_genes)}
        train_y=y[train];test_y=y[test]
        one=exact_singles(train_x,fold_lookup,feature_genes,train_y,'top20')
        selected_genes=shortlist(one,12)
        two=pair_search(train_x,fold_lookup,selected_genes,train_y,'top20')
        df=pd.DataFrame(one+two)
        for floor in FLOORS:
            d=df[df.target_recovery>=floor].sort_values(['post_gate_target_fraction','n_markers','target_recovery'],ascending=[False,True,False])
            chosen=d.iloc[0].to_dict();keep=apply(chosen,test_x,fold_lookup)
            tp,n=int(test_y[keep].sum()),int(keep.sum())
            row=add_row({k:chosen[k] for k in ('gate_id','gate_label','n_markers','rules','logic','source')},tp,n,test_y,'top20')
            row.update(held_out_sample=sample,recovery_floor=floor,train_purity=chosen['post_gate_target_fraction'],
                       train_recovery=chosen['target_recovery'],evaluation='leave-one-sample-out: surface features and thresholds selected on training cells only',
                       train_candidate_genes=';'.join(selected_genes))
            rows.append(row);pooled[floor]['tp']+=tp;pooled[floor]['n']+=n;pooled[floor]['targets']+=int(test_y.sum());pooled[floor]['total']+=len(test_y)
        print(f'DIV30 held-sample gate validation: {sample}',flush=True)
    pd.DataFrame(rows).to_csv(out/'tables/held_sample_surface_gate_validation.tsv',sep='\t',index=False)
    summary=[]
    for floor,v in pooled.items():
        if v['total']:
            row={k:float(a) for k,a in metric(v['tp'],v['n'],v['targets'],v['total']).items()}
            row.update(recovery_floor=floor,evaluation='pooled held-out sample predictions; different training-selected gates across folds',
                       limitation='scores/reference mapping fixed before surface CV; no external prospective labels; correlated samples are not independent biological replicates')
            summary.append(row)
    pd.DataFrame(summary).to_csv(out/'tables/held_sample_surface_gate_summary.tsv',sep='\t',index=False)
    return summary


def benchmark(cells, expression, lookup, target_dict, out):
    path=DIV90/'gate_summary.json';data=json.loads(path.read_text())
    item=data['experimental_depletion']
    rules=json.loads(item['rules']) if isinstance(item['rules'],str) else item['rules']
    assert item['action']=='remove' and item['logic']=='OR'
    assert {c['gene'] for c in rules}=={'FGFR2','PTPRS'} and all(c['op']=='>' for c in rules)
    rows=[];spec=None
    if {'FGFR2','PTPRS'}<=set(lookup):
        # Frozen DIV90 depletion OR is equivalent to retained <= AND.
        spec=gate([(c['gene'],'<=',c['threshold']) for c in rules],'AND','frozen_DIV90_FGFR2_PTPRS_depletion_benchmark')
        keep=apply(spec,expression,lookup)
        for target,y in target_dict.items(): rows.append(add_row(spec,int(y[keep].sum()),int(keep.sum()),y,target))
        write_json(out/'provenance/frozen_DIV90_benchmark_rule.json',dict(source=str(path),source_key='experimental_depletion',source_gate_id=item['gate_id'],original_rule=rules,retained_rule=spec,
                                                                      caution='RNA threshold transfer is diagnostic only; no stage equivalence assumed'))
    pd.DataFrame(rows).to_csv(out/'tables/frozen_DIV90_gate_on_DIV30.tsv',sep='\t',index=False)
    return rows


def temporal_table(marker_table, out):
    p=DIV90/'tables/surface_state_markers.tsv'
    if not p.exists(): return pd.DataFrame()
    old=pd.read_csv(p,sep='\t');old=old[old.state=='PV-biased'].set_index('gene')
    rows=[]
    for r in marker_table.to_dict('records'):
        g=r['gene']; d90=float(old.loc[g,'cohens_d']) if g in old.index else np.nan
        d30=r['cohens_d']; measured=r['classification']; early=abs(d30)>=.3;late=np.isfinite(d90) and abs(d90)>=.3
        if measured in ('generic neuronal maturation','progenitor marker'): label='MATURATION-RELATED'
        elif early and late and d30*d90>0: label='EARLY-STABLE'
        elif early and np.isfinite(d90) and not late: label='EARLY-TRANSIENT'
        elif not early and late: label='LATE-EMERGING'
        elif early: label='SURFACE-INFORMATIVE BUT NOT SUBTYPE-SPECIFIC'
        else: label='NOT ESTABLISHED'
        rows.append(dict(gene=g,DIV30_top20_cohens_d=d30,DIV90_PV_biased_cohens_d=d90,
                         DIV30_adjusted_classification=measured,temporal_classification=label,
                         rule='descriptive |d|>=0.3, same direction for early-stable; measured maturation/progenitor association takes precedence',
                         limitation='different age-specific targets and denominators; cross-sectional heuristic, no lineage or formal age-interaction test'))
    table=pd.DataFrame(rows);table.to_csv(out/'tables/temporal_surface_marker_classification.tsv',sep='\t',index=False)
    return table


def subgroup_gate_metrics(cells, expression, lookup, target_dict, best, bench, out):
    winners={(r['target'],r['gate_id']):r for r in best+bench}; rows=[]
    for (target,_),row in winners.items():
        keep=apply(row,expression,lookup);y=target_dict[target]
        for field in ('sample','cell_line','condition'):
            if field not in cells: continue
            for label,idx in cells.groupby(field,observed=True,dropna=False).indices.items():
                idx=np.asarray(idx);tp=int((keep[idx]&y[idx]).sum());n=int(keep[idx].sum())
                r={k:float(v) for k,v in metric(tp,n,int(y[idx].sum()),len(idx)).items()}
                r.update(target=target,gate_id=row['gate_id'],gate_label=row['gate_label'],stratum_type=field,stratum=str(label),
                         stratum_n=len(idx),stratum_target_n=int(y[idx].sum()),
                         evaluation='fixed global discovery threshold; subgroup descriptive sensitivity, not held-out validation')
                rows.append(r)
    pd.DataFrame(rows).to_csv(out/'tables/fixed_gate_by_sample_line_condition.tsv',sep='\t',index=False)


def surface_report(out, table, best, validation, bench, marker_table, feature_count):
    best=pd.DataFrame(best);main=best[(best.target=='top20')&(best.maximum_markers==3)]
    top20=table[table.target=='top20']
    unconstrained=[]
    for target,frame in table.groupby('target',sort=False):
        for floor in (0.,.1,.2,.3):
            allowed=frame[frame.target_recovery>=floor]
            if len(allowed):
                r=allowed.sort_values(['post_gate_target_fraction','target_recovery','n_markers'],ascending=[False,False,True]).iloc[0].to_dict()
                r['recovery_floor']=floor;unconstrained.append(r)
    pd.DataFrame(unconstrained).to_csv(out/'tables/surface_maximum_and_low_recovery_summary.tsv',sep='\t',index=False)
    lines=['# DIV30 surface discovery: quantitative interpretation','',
           f'The independent DIV30 search considered {feature_count} annotated, expressed surface candidates. All score, mapping, and supervised-model genes were excluded from gate features. Targets are individual cells ranked by developmental resemblance; they are not confirmed future PV cells.', '',
           '## 6. Does FGFR2/PTPRS contain useful information at DIV30?', '']
    b=next((r for r in bench if r['target']=='top20'),None)
    if b:
        clauses=json.loads(b['rules'])
        exact=' AND '.join(f"{r['gene']} <= {r['threshold']:.17g}" for r in clauses)
        label='weak purification' if b['fold_enrichment']<1.2 or b['composition_change_pp']<5 else b['enrichment_interpretation']+' enrichment'
        lines.extend([f'The frozen DIV90 depletion benchmark retains `{exact}`. The original source is `gate_summary.json["experimental_depletion"]`; all threshold precision is retained.', '',
            f'At DIV30 it changes the top20 target fraction from **{b["starting_target_fraction"]:.1%} to {b["post_gate_target_fraction"]:.1%}** '
            f'({b["fold_enrichment"]:.2f}×; {b["composition_change_pp"]:+.2f} percentage points), with {b["target_recovery"]:.1%} target recovery and {b["total_cell_yield"]:.1%} total yield. '
            f'This is **{label}**. Direction agreement alone is not replication. The RNA cut points are a cross-age diagnostic, not translated fluorescence thresholds.',''])
    diagnostic=marker_table[marker_table.gene.isin(['FGFR2','PTPRS','FAT3','PTPRM','ERBB4'])]
    lines.extend(['| Diagnostic marker | Cohen d, top20 vs other | Distribution overlap | Adjusted early PV r | Classification |',
                  '|---|---:|---:|---:|---|'])
    for r in diagnostic.to_dict('records'):
        lines.append(f'| {r["gene"]} | {r["cohens_d"]:.3f} | {r["distribution_overlap"]:.1%} | {r["partial_r_early_pv"]:.3f} | {r["classification"]} |')
    lines.extend(['','Each marker’s optimized DIV30 single-threshold composition, enrichment, recovery, and yield are recorded separately in `tables/DIV90_marker_DIV30_diagnostic_gate_metrics.tsv`. Diagnostic genes excluded from independent search remain explicitly flagged.','',
                  '## 7. Is there a stronger independently discovered DIV30 surface gate?','',
                  '| Required recovery | Retained DIV30 RNA rule | Target composition | Enrichment | Actual recovery | Total yield |',
                  '|---:|---|---:|---:|---:|---:|'])
    for r in main[main.recovery_floor.isin([.5,.7,.8])].sort_values('recovery_floor').to_dict('records'):
        lines.append(f'| {r["recovery_floor"]:.0%} | {r["gate_label"]} | {r["post_gate_target_fraction"]:.1%} | {r["fold_enrichment"]:.2f}× | {r["target_recovery"]:.1%} | {r["total_cell_yield"]:.1%} |')
    lines.extend(['','These are the best **evaluated** gates at each recovery floor, selected for actual retained composition. One-marker thresholds cover every observed cut point. Two-marker search covers all signed AND/OR combinations among 18 independently shortlisted genes on zero plus nine positive-expression quantiles. Triplets extend up to 14 pair seeds using up to 12 genes and enter results only with at least a 2-percentage-point purity gain over the complete evaluated one/two-marker frontier at equal or higher recovery. This bounded search does not establish a global optimum.','',
                  '## 8. Are the best markers subtype-associated or developmental confounds?',''])
    important=set()
    for r in main[main.recovery_floor.isin([.5,.7,.8])].to_dict('records'):
        important.update(c['gene'] for c in json.loads(r['rules']))
    lines.extend(['| Selected marker | Early PV partial r | Neighbor PV partial r | Early PV r plus sample adjustment | Classification |',
                  '|---|---:|---:|---:|---|'])
    for r in marker_table[marker_table.gene.isin(important)].to_dict('records'):
        lines.append(f'| {r["gene"]} | {r["partial_r_early_pv"]:.3f} | {r["partial_r_neighbor_pv"]:.3f} | {r["partial_r_early_pv_plus_sample"]:.3f} | {r["classification"]} |')
    lines.extend(['','The direct/consensus model adjusts progenitor, maturation, MGE, and library-depth scores; the neighbor model adjusts progenitor, maturation, and depth. Sample fixed effects provide an additional sensitivity. “Independent of measured maturation” denotes residual cross-sectional association under these models, not fate specificity or independent statistical significance. Nominal cell-level p/q values are descriptive and do not turn cells into independent biological replicates. Gate-level adjusted models and fixed-gate sample/line/condition sensitivities accompany the marker models.','',
                  '## 9. What is the strongest achievable enrichment in this search?',''])
    for floor in (0.,.1,.2,.3):
        r=next(v for v in unconstrained if v['target']=='top20' and v['recovery_floor']==floor)
        desc='Unconstrained maximum' if floor==0 else f'At least {floor:.0%} target recovery'
        lines.append(f'- {desc}: {r["post_gate_target_fraction"]:.1%} target composition, {r["fold_enrichment"]:.2f}× enrichment, {r["target_recovery"]:.1%} actual recovery, {r["total_cell_yield"]:.1%} total yield ({int(r["retained_n"]):,} retained cells), `{r["gate_label"]}`.')
    lines.extend(['','A tiny retained fraction can achieve an impressive purity maximum; composition, recovery, and absolute retained counts must be read together. The low-recovery table covers all independently defined target sensitivities.','',
                  '## 10–11. Is DIV30 worth an experimental sort, and which gate is a candidate?','',
                  '| Training recovery floor | Pooled held-sample target composition | Enrichment | Target recovery | Total yield |',
                  '|---:|---:|---:|---:|---:|'])
    for r in validation:
        lines.append(f'| {r["recovery_floor"]:.0%} | {r["post_gate_target_fraction"]:.1%} | {r["fold_enrichment"]:.2f}× | {r["target_recovery"]:.1%} | {r["total_cell_yield"]:.1%} |')
    useful=[r for r in validation if r['fold_enrichment']>=1.5 and r['target_recovery']>=.5]
    candidate=main.loc[main.recovery_floor.eq(.7)].iloc[0]
    if useful:
        lines.extend(['',f'The held-sample surface-selection exercise supports **exploratory testing of a DIV30 candidate surface phenotype**: the practical in-sample example is `{candidate.gate_label}`. '
                      'This is a hypothesis for surface-protein measurement and a prospective endpoint experiment. RNA thresholds cannot define an actual live-cell sort, and no single panel has yet been validated against later PV protein.'])
    else:
        lines.extend(['','The current held-sample results do **not** establish a practical DIV30 subtype sort that combines at least 1.5× enrichment with 50% recovery. '
                      f'The strongest practical in-sample candidate, `{candidate.gate_label}`, remains exploratory. '
                      'DIV30 can still locate a candidate developmental window and guide surface-protein measurement; these data cannot establish that waiting until DIV90 or changing culture is superior.'])
    lines.extend(['','Every held sample receives a gate chosen using the other samples, with the surface shortlist and thresholds retrained. The pooled row therefore combines different fold-specific gates; it is not validation of one fixed assay. '
                  'Scores and the DIV90 reference are frozen before this surface validation. Samples can share lines and culture context, and none supplies a prospective PV-fate label. No lineage tracing, independent statistical significance, or validated experimental gate is claimed.','',
                  'The linked hypothesis is a DIV30 transcriptional antecedent → a DIV90 PV-associated state → the user-reported later PV protein endpoint in this MGEO system. Cross-sectional resemblance does not prove those transitions for individual cells.',''])
    (out/'SURFACE_REPORT.md').write_text('\n'.join(lines))


def figures(cells, expression, lookup, marker_table, best, validation, bench, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    def save(fig,name):
        for fmt in ('png','pdf','svg'):
            (out/f'figures/{fmt}').mkdir(parents=True,exist_ok=True)
            fig.savefig(out/f'figures/{fmt}/{name}.{fmt}',dpi=190,bbox_inches='tight')
        plt.close(fig)
    fields=['progenitor_score','mge_score','maturation_score','progression_score','NKX2-1','LHX6','ERBB4','MAF','MAFB','MEF2C','DCX','log_counts']
    fig,axes=plt.subplots(3,4,figsize=(15,10));binrows=[]
    for ax,col in zip(axes.flat,fields):
        x=cells[col].to_numpy(float) if col in cells else expression[:,lookup[col]]
        y=cells.early_pv_score.to_numpy(float)
        if np.unique(x).size<3: ax.text(.1,.5,'Insufficient variation');continue
        hb=ax.hexbin(x,y,gridsize=35,mincnt=1,bins='log',cmap='Blues',rasterized=True)
        edges=np.unique(np.quantile(x,np.linspace(0,1,16)));bid=np.searchsorted(edges[1:-1],x,side='right')
        bx,by=[],[]
        for j in np.unique(bid):
            take=bid==j
            bx.append(float(np.median(x[take])));by.append(float(np.mean(y[take])))
            binrows.append(dict(axis=col,quantile_bin=int(j),n=int(take.sum()),axis_median=bx[-1],early_pv_mean=by[-1],
                                early_pv_q25=float(np.quantile(y[take],.25)),early_pv_q75=float(np.quantile(y[take],.75))))
        if len(bx)>=3 and len(set(bx))==len(bx):
            from scipy.interpolate import PchipInterpolator
            smooth_x=np.linspace(min(bx),max(bx),180)
            ax.plot(smooth_x,PchipInterpolator(bx,by)(smooth_x),color='#B6472B',lw=2)
        else: ax.plot(bx,by,color='#B6472B',lw=2)
        ax.set(xlabel=col,ylabel='Early PV consensus rank')
    fig.suptitle('DIV30 candidate antecedent signal across continuous developmental axes\nCross-sectional quantile trends; timing and lineage are not established',y=1.02)
    save(fig,'DIV30_developmental_continuous_axes')
    pd.DataFrame(binrows).to_csv(out/'tables/developmental_axis_quantile_trends.tsv',sep='\t',index=False)
    d=pd.DataFrame(best);d=d[d.target=='top20']
    fig,axes=plt.subplots(1,3,figsize=(15,4.8))
    for k,c in [(1,'#668C72'),(2,'#286294'),(3,'#BD803B')]:
        q=d[d.maximum_markers==k].sort_values('recovery_floor')
        axes[0].plot(q.target_recovery,q.post_gate_target_fraction,'o-',color=c,label=f'At most {k} markers')
    base=(cells.early_pv_score>=cells.early_pv_score.quantile(.8)).mean()
    axes[0].axhline(base,color='grey',ls='--',label=f'Baseline {base:.1%}')
    if validation:
        v=pd.DataFrame(validation);axes[0].plot(v.target_recovery,v.post_gate_target_fraction,'s--',color='#982B37',label='Held-sample selected gates')
    if bench:
        b=[r for r in bench if r['target']=='top20'][0]
        axes[0].scatter([b['target_recovery']],[b['post_gate_target_fraction']],marker='X',s=100,color='black',label='Frozen DIV90 benchmark')
    axes[0].set(xlabel='Target recovery',ylabel='Retained early PV top20 fraction',title='Composition at useful recovery');axes[0].legend(fontsize=8)
    strongest=marker_table.assign(effect=marker_table.partial_r_early_pv.abs()).nlargest(12,'effect').sort_values('partial_r_early_pv')
    axes[1].barh(strongest.gene,strongest.partial_r_early_pv,color='#286294',label='Measured developmental axes')
    axes[1].scatter(strongest.partial_r_early_pv_plus_sample,strongest.gene,color='#BD803B',label='Also sample adjusted',zorder=3)
    axes[1].axvline(0,color='grey',lw=.6);axes[1].set(xlabel='Partial correlation with early PV score',title='Residual association is exploratory');axes[1].legend(fontsize=8)
    diagnostic=marker_table[marker_table.gene.isin(DIAGNOSTIC)].sort_values('cohens_d')
    axes[2].barh(diagnostic.gene,diagnostic.cohens_d,color='#668C72');axes[2].axvline(0,color='grey',lw=.6)
    axes[2].set(xlabel="Cohen's d: top20 versus other DIV30 cells",title='DIV90 markers as DIV30 diagnostics')
    fig.suptitle('Independent DIV30 surface search: purity first; no confirmed future PV labels',y=1.03);fig.tight_layout();save(fig,'DIV30_surface_search_and_confound_summary')


def run_surface(out, cells, expression, genes, excluded_genes, sst_stable=False):
    out=Path(out); (out/'tables').mkdir(parents=True,exist_ok=True);(out/'provenance').mkdir(parents=True,exist_ok=True)
    if len(cells)!=expression.shape[0] or len(genes)!=expression.shape[1] or not cells.cell_id.is_unique:
        raise ValueError('DIV30 metadata/expression dimensions or barcode uniqueness invalid')
    lookup={g:i for i,g in enumerate(genes)}
    annotation,feature_genes=annotation_catalog(expression,genes,excluded_genes,out)
    if set(feature_genes)&set(excluded_genes): raise AssertionError('Circular target/mapping gene in gate search')
    target_dict=targets(cells,sst_stable)
    write_json(out/'provenance/DIV30_surface_search_design.json',dict(n_cells=len(cells),eligible_genes=feature_genes,
        excluded_score_mapping_supervised_genes=sorted(set(excluded_genes)),surface_feature_count=len(feature_genes),
        targets={k:dict(n=int(v.sum()),fraction=float(v.mean())) for k,v in target_dict.items()},
        signed_single_thresholds='every observed expression cut point; > and <=',
        pair_search='18 independently DIV30-screened genes per target; all signed AND/OR pairs at zero plus 9 positive-expression quantiles',
        triplets='extend up to 14 best pair seeds with up to 12 genes; retain only >=0.02 purity gain at equal-or-higher recovery vs every evaluated <=2-marker gate',
        recovery_floors=FLOORS,ranking='maximize retained target composition subject to target recovery floor',
        validation='top20 leave-one-sample-out retraining of all univariate candidates, 12-gene pair shortlist, and thresholds',
        limitations=['no prospective fate labels','sample units may not be independent biological replicates','scRNA surface transcripts do not establish live-cell protein accessibility','bounded pair/triplet search, not a proven global optimum']))
    all_rows=[];all_best=[];triplets=[]
    archive=out/'tables/surface_gate_complete_counts.h5'
    with h5py.File(archive,'w') as h5:
        h5.attrs['metric_definition']='For all encoded gates: starting_target_fraction=T/N; post_gate_target_fraction=tp/n; fold_enrichment=(tp/n)/(T/N); target_recovery=tp/T; total_cell_yield=n/N. Single > is complement of stored <= counts. Pair arrays include all signed AND/OR rules.'
        for target,y in target_dict.items():
            group=h5.create_group(target);group.attrs['n_total']=len(y);group.attrs['n_target']=int(y.sum())
            purity_by_target_count=np.full(int(y.sum())+1,-np.inf)
            purity_by_target_count[-1]=float(y.mean())
            one=exact_singles(expression,lookup,feature_genes,y,target,group,purity_by_target_count)
            selected_genes=shortlist(one)
            two=pair_search(expression,lookup,selected_genes,y,target,group,purity_by_target_count)
            group.attrs['pair_genes']=json.dumps(selected_genes)
            matched_frontier=np.maximum.accumulate(purity_by_target_count[::-1])[::-1]
            group.create_dataset('best_single_or_pair_purity_at_minimum_target_count',data=matched_frontier,compression='gzip')
            evaluated,accepted=triplet_search(expression,lookup,selected_genes,y,target,one+two,matched_frontier)
            triplets.extend(evaluated);rows=one+two+accepted
            all_rows.extend(rows);all_best.extend(best_at_floors(rows))
            print(f'DIV30 {target}: {int(y.sum())} target cells, {len(feature_genes)} single-marker candidates, {len(selected_genes)} pair features, {len(accepted)} material triplets',flush=True)
    table=pd.DataFrame(all_rows).drop_duplicates(['target','gate_id'])
    table.to_csv(out/'tables/surface_gate_candidate_registry.tsv.gz',sep='\t',index=False)
    pd.DataFrame(all_best).to_csv(out/'tables/surface_gate_best_at_recovery_floors.tsv',sep='\t',index=False)
    pd.DataFrame(triplets).to_csv(out/'tables/surface_triplet_materiality_audit.tsv.gz',sep='\t',index=False)
    (out/'tables/README_surface_gate_complete_counts.md').write_text(
        '# Complete DIV30 gate archive\n\n`surface_gate_complete_counts.h5` records every evaluated single and pair gate exactly. '
        'Each target group has attributes `n_total=N` and `n_target=T`. Single marker groups store every expression `threshold`, '
        '`retained_n_le=n`, and `target_n_le=tp` for retained expression <= threshold. The > gate has counts N-n and T-tp. '
        'Pair groups name both genes, list their threshold arrays, and store each signed AND/OR combination as columns retained_n,target_n. '
        'Rows reshape to (number of thresholds for gene_a, number of thresholds for gene_b), C order.\n\n'
        'For every gate: baseline=T/N; post-gate composition=tp/n; enrichment=(tp/n)/(T/N); recovery=tp/T; yield=n/N. '
        'Empty retained populations have undefined composition. All counts are integers and thresholds are saved at float64 precision. '
        'The compressed TSV registry is a compact subset of frontiers and recovery-floor winners, not the exhaustive archive. '
        'All triplet evaluations and the matched-recovery materiality flag are in surface_triplet_materiality_audit.tsv.gz. '
        'Only triplets improving purity by at least 0.02 enter the final frontier.\n')
    diagnostics=sorted(set(feature_genes)|({g for g in DIAGNOSTIC if g in lookup}))
    marker_table=confound_marker_table(cells,expression,lookup,diagnostics,target_dict,out)
    marker_table.merge(annotation[['gene','allowed','exclusion_reason','annotation_category']],on='gene',how='left').to_csv(out/'tables/surface_marker_evidence.tsv',sep='\t',index=False)
    gate_confound_table(cells,expression,lookup,all_best,out)
    diagnostic_single_metrics(expression,lookup,target_dict,annotation,out)
    validation=held_sample_validation(cells,expression,lookup,feature_genes,target_dict['top20'],out)
    bench=benchmark(cells,expression,lookup,target_dict,out)
    subgroup_gate_metrics(cells,expression,lookup,target_dict,all_best,bench,out)
    temporal=temporal_table(marker_table,out)
    surface_report(out,table,all_best,validation,bench,marker_table,len(feature_genes))
    figures(cells,expression,lookup,marker_table,all_best,validation,bench,out)
    summary=dict(n_cells=len(cells),eligible_surface_features=len(feature_genes),n_targets=len(target_dict),
                 best_at_recovery_floors=all_best,held_sample_validation=validation,frozen_DIV90_benchmark=bench,
                 full_count_archive=str(archive),complete=True,
                 limitation='Exploratory cross-sectional antecedent discovery; no lineage tracing, protein gate validation, or independent statistical significance established.')
    write_json(out/'surface_summary.json',summary)
    return summary
