"""Full-data marker/display context without fitting clusters or transferring labels.

Read the immutable Step 06 count object in sparse row chunks, normalize using
all measured genes, and apply the fixed pilot signature weights. Only selected
canonical genes, diagnostic scores, metadata and existing UMAP coordinates are
saved in a new small display AnnData. This sidecar is NOT a replacement for the
full-gene raw input and contains no inferred full-data cluster memberships.
"""
from pathlib import Path
import json
import h5py
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import wasserstein_distance
from .provenance import write_json
from .validation_programs import normalize_log1p_cpm


class FullDataProjector:
    """Stream fixed program scores onto the existing full-data display manifold.

    Parameters
    ----------
    full_source : path
        Immutable Step 06 H5AD with raw CSR X and existing obsm['X_umap'].
    output_dir : path
        New validation run's full_data_projection directory.
    chunk_size : int
        Number of raw rows read/normalized at once; controls memory only.
    """

    def __init__(self, full_source, output_dir, chunk_size=4096):
        """Record paths; no source file is opened for writing."""
        self.full_source = Path(full_source)
        self.output_dir = Path(output_dir)
        self.chunk_size = int(chunk_size)
        if self.chunk_size < 1:
            raise ValueError('chunk_size must be positive')

    @staticmethod
    def _note(stage, **fields):
        """Print a flushed progress event for the containing frozen job's log."""
        print(json.dumps(dict(stage='full_data_projection.'+stage, **fields)), flush=True)

    def run(self, pilot, scorer):
        """Save measured full-data context and verify exact pilot agreement.

        Parameters
        ----------
        pilot : AnnData
            The unchanged 12,000-cell pilot, with its original cell IDs.
        scorer : CanonicalProgramScorer
            Already fitted exclusively to the pilot. No controls or thresholds
            are refitted on full data.

        Returns
        -------
        dict
            Saved object/source-table paths, cell/gene counts, maximum pilot
            score/expression disagreement, and numerical coverage diagnostics.

        Outputs and omissions
        ---------------------
        ``canonical_full_data.h5ad``: selected canonical normalized X only;
        obs sample/genotype/sex/QC, is_pilot and transcriptional phase proxy;
        obsm original Step06 UMAP and fixed program scores; uns score-model
        JSON and explicit content inventory. There is no .raw, counts layer,
        graph, PCA, newly fitted UMAP, full-gene normalized matrix, transferred
        cluster label, locked biological annotation or batch correction.
        Compressed TSVs contain every plotted gene/score/coordinate value;
        binned UMAP and per-sample program tables quantify pilot coverage.
        """
        if scorer.result is None:
            raise ValueError('Scorer must be fitted on the pilot first')
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self._note('start', source=str(self.full_source), chunk_size=self.chunk_size)
        full = ad.read_h5ad(self.full_source, backed='r')
        try:
            if not np.array_equal(full.var_names.astype(str), scorer.gene_ids):
                raise ValueError('Full-data feature IDs/order do not match pilot score model')
            if not np.array_equal(full.var['gene_symbol'].astype(str).to_numpy(), scorer.symbols):
                raise ValueError('Full-data gene symbols do not match pilot score model')
            if not full.obs_names.is_unique or not pilot.obs_names.is_unique:
                raise ValueError('Cell IDs must be unique')
            positions = full.obs_names.get_indexer(pilot.obs_names)
            if np.any(positions < 0):
                raise ValueError('Pilot cells are absent from full display source')
            # Matching cell IDs is necessary but insufficient: the same cell
            # must retain its registered sample identity in both checkpoints.
            full_pilot_samples = full.obs['technical_sample_id'].iloc[positions].astype(str).to_numpy()
            if not np.array_equal(full_pilot_samples, pilot.obs['technical_sample_id'].astype(str).to_numpy()):
                raise ValueError('Full-data technical sample identity differs from pilot')
            required = ['technical_sample_id','submitted_sample_name','genotype','sex','design_group',
                        'total_counts','n_genes_by_counts','pct_counts_mt']
            obs = full.obs[[c for c in required if c in full.obs]].copy()
            for column in required[:5]:
                if column in obs:
                    obs[column] = obs[column].astype(str)
            # Preserve the source coordinate values AND dtype; a float32 cast
            # would unnecessarily alter a float64 source manifold.
            coords = np.asarray(full.obsm['X_umap']).copy()
            var = full.var.iloc[scorer.canonical_indices][['gene_symbol']].copy()
            n_cells, n_genes = full.shape
        finally:
            full.file.close()
        if not np.all(np.isfinite(coords)):
            raise ValueError('Full-data UMAP coordinates contain missing/nonfinite values')
        obs['is_pilot'] = False
        obs.iloc[positions, obs.columns.get_loc('is_pilot')] = True
        selected_expression = np.empty((n_cells,len(scorer.canonical_indices)),dtype=np.float32)
        scores = np.empty((n_cells,len(scorer.names)),dtype=np.float32)
        detected = np.empty_like(scores)
        with h5py.File(self.full_source, 'r') as handle:
            group = handle['X']
            if group.attrs.get('encoding-type') != 'csr_matrix':
                raise ValueError('Full-data X must be raw CSR, not a dense or CSC matrix')
            indptr = group['indptr'][:]
            for chunk, start in enumerate(range(0,n_cells,self.chunk_size)):
                stop = min(n_cells,start+self.chunk_size)
                left,right = int(indptr[start]),int(indptr[stop])
                raw = sparse.csr_matrix((group['data'][left:right],group['indices'][left:right],indptr[start:stop+1]-left),shape=(stop-start,n_genes))
                if not np.issubdtype(raw.dtype,np.integer) or np.any(raw.data<0):
                    raise ValueError('Full-data source is not nonnegative integer raw counts')
                if np.any(np.asarray(raw.sum(axis=1)).ravel() <= 0):
                    raise ValueError('Full-data source contains an empty-count cell')
                normalized = normalize_log1p_cpm(raw,scorer.settings['normalization_target'])
                scores[start:stop] = scorer.transform(normalized)
                selected_expression[start:stop] = normalized[:,scorer.canonical_indices].toarray()
                for j,signal in enumerate(scorer.signal_indices):
                    detected[start:stop,j] = np.asarray((raw[:,signal]>0).sum(axis=1)).ravel() if signal else np.nan
                if chunk % 16 == 0 or stop == n_cells:
                    self._note('chunk_complete', rows_complete=stop, rows_total=n_cells)
        score_error = float(np.nanmax(np.abs(scores[positions].astype(float)-scorer.result['scores'].to_numpy())))
        expression_error = float(np.max(np.abs(selected_expression[positions].astype(float)-scorer.result['canonical_expression'].to_numpy())))
        coordinate_error = float(np.max(np.abs(coords[positions]-pilot.obsm['X_umap_step06_display'])))
        if score_error > 2e-5 or expression_error > 2e-5 or not np.array_equal(coords[positions],pilot.obsm['X_umap_step06_display']):
            raise ValueError('Pilot/full marker or coordinate alignment failed: %s' % dict(score=score_error,expression=expression_error,coordinates=coordinate_error))
        s_index,g_index = scorer.names.index('S_phase'),scorer.names.index('G2M')
        phase = np.full(n_cells,'S',dtype=object)
        phase[scores[:,g_index]>scores[:,s_index]] = 'G2/M'
        phase[(scores[:,g_index]<0)&(scores[:,s_index]<0)] = 'G1-like'
        obs['validation_phase'] = pd.Categorical(phase,categories=['G1-like','S','G2/M'])
        obj = ad.AnnData(X=sparse.csr_matrix(selected_expression), obs=obs, var=var)
        obj.obsm['X_umap_step06_display'] = coords
        obj.obsm['X_validation_program_scores'] = scores
        obj.obsm['X_validation_signal_detected_genes'] = detected
        obj.uns['validation_program_names'] = np.asarray(scorer.names,dtype=str)
        obj.uns['program_names'] = np.asarray(scorer.names,dtype=str)
        obj.uns['score_model_json'] = json.dumps(scorer.model)
        inventory = dict(X='natural log1p CPM of selected canonical genes only; denominator uses all 19071 measured genes',
                         obs='source sample/design/QC metadata, is_pilot, transcriptional cell-cycle proxy',
                         obsm=dict(X_umap_step06_display='unchanged original Step06 UMAP',X_validation_program_scores='fixed pilot score coefficients',X_validation_signal_detected_genes='number of signature genes with positive raw counts'),
                         not_saved=['all-gene raw counts','all-gene normalized expression','layers','raw','graphs','PCA','new UMAP','full-data cluster labels','transferred annotations'],
                         biological_context='dissected mouse MGE, E14.5; not organoids')
        obj.uns['asset_inventory_json'] = json.dumps(inventory)
        object_path = self.output_dir/'canonical_full_data.h5ad'
        self._note('write_h5ad', shape=list(obj.shape))
        obj.write_h5ad(object_path,compression='lzf')
        self._note('write_numeric_source_tables')
        coordinate_table = obs.copy()
        coordinate_table['umap_1'] = coords[:,0]
        coordinate_table['umap_2'] = coords[:,1]
        coordinate_table.rename_axis('cell_id').to_csv(self.output_dir/'full_data_cells_and_umap.tsv.gz',sep='\t',compression='gzip')
        pd.DataFrame(scores,index=obs.index,columns=scorer.names).rename_axis('cell_id').to_csv(self.output_dir/'full_data_program_scores.tsv.gz',sep='\t',compression='gzip')
        pd.DataFrame(selected_expression,index=obs.index,columns=var['gene_symbol'].astype(str)).rename_axis('cell_id').to_csv(self.output_dir/'full_data_canonical_expression.tsv.gz',sep='\t',compression='gzip')
        pd.DataFrame(detected,index=obs.index,columns=scorer.names).rename_axis('cell_id').to_csv(self.output_dir/'full_data_signature_detection.tsv.gz',sep='\t',compression='gzip')
        self._note('coverage_diagnostics')
        coverage = self._coverage_tables(obs,coords,scores,detected,scorer)
        summary = dict(status='completed_display_and_marker_scoring_only',object_path=str(object_path),n_cells=n_cells,n_selected_genes=len(var),
                       n_measured_genes_for_normalization=n_genes,n_pilot_cells=len(positions),
                       pilot_score_max_abs_difference=score_error,pilot_expression_max_abs_difference=expression_error,
                       pilot_coordinate_max_abs_difference=coordinate_error,pilot_sample_identity_exact=True,
                       coordinate_dtype=str(coords.dtype),coordinates_saved_without_cast=True,
                       normalization_refitted=False,score_controls_refitted_on_full_data=False,full_data_clustering=False,
                       label_transfer=False,umap_refitted=False,batch_correction=False,regression=False,
                       full_source=str(self.full_source),asset_inventory=inventory,coverage=coverage,
                       tables=dict(cells='full_data_cells_and_umap.tsv.gz',scores='full_data_program_scores.tsv.gz',
                                   expression='full_data_canonical_expression.tsv.gz',detection='full_data_signature_detection.tsv.gz',
                                   manifold_coverage='umap_bin_pilot_coverage.tsv',program_coverage='full_vs_pilot_program_coverage.tsv'))
        write_json(self.output_dir/'full_data_projection_summary.json',summary)
        self._note('complete',**{k:summary[k] for k in ['n_cells','n_selected_genes','pilot_score_max_abs_difference']})
        return summary

    def _coverage_tables(self, obs, coords, scores, detected, scorer):
        """Quantify representation without assigning or extrapolating cell labels.

        UMAP bins describe coverage of a fixed display only. Expected pilot
        counts account for the deliberately balanced 1,000/sample selection,
        rather than treating the pilot as a uniform sample from all cells.
        Program Wasserstein distances are divided by full-data standard
        deviation for scale comparability; they are descriptive, not tests.
        """
        samples = obs['technical_sample_id'].astype(str).to_numpy()
        selected = obs['is_pilot'].to_numpy(dtype=bool)
        full_counts = pd.Series(samples).value_counts()
        pilot_counts = pd.Series(samples[selected]).value_counts()
        probability = np.array([pilot_counts.get(sample,0)/full_counts[sample] for sample in samples])
        xedges = np.linspace(coords[:,0].min(),coords[:,0].max(),41)
        yedges = np.linspace(coords[:,1].min(),coords[:,1].max(),41)
        full_hist = np.histogram2d(coords[:,0],coords[:,1],bins=(xedges,yedges))[0]
        pilot_hist = np.histogram2d(coords[selected,0],coords[selected,1],bins=(xedges,yedges))[0]
        expected_hist = np.histogram2d(coords[:,0],coords[:,1],bins=(xedges,yedges),weights=probability)[0]
        rows = []
        for i in range(40):
            for j in range(40):
                n,e,k = full_hist[i,j],expected_hist[i,j],pilot_hist[i,j]
                rows.append(dict(x_bin=i,y_bin=j,x_low=xedges[i],x_high=xedges[i+1],y_low=yedges[j],y_high=yedges[j+1],
                                 full_cells=int(n),pilot_cells=int(k),expected_balanced_pilot_cells=e,
                                 observed_expected_ratio=k/e if e>0 else np.nan,
                                 well_populated_full_bin=bool(n>=100),adequate_expected_pilot=bool(e>=5),
                                 underrepresented=bool(n>=100 and e>=5 and k/e<0.5)))
        bins = pd.DataFrame(rows)
        bins.to_csv(self.output_dir/'umap_bin_pilot_coverage.tsv',sep='\t',index=False)
        program_rows = []
        for sample in ['ALL']+sorted(full_counts.index):
            mask = np.ones(len(obs),dtype=bool) if sample=='ALL' else samples==sample
            for j,name in enumerate(scorer.names):
                values=scores[mask,j].astype(float)
                subset=scores[mask&selected,j].astype(float)
                valid=np.isfinite(values)
                if valid.any():
                    standard_deviation=np.std(values[valid])
                    full_positive=(values>scorer.settings['positive_score_threshold'])&(detected[mask,j]>=scorer.settings['min_detected_genes_positive'])
                    pilot_positive=(subset>scorer.settings['positive_score_threshold'])&(detected[mask&selected,j]>=scorer.settings['min_detected_genes_positive'])
                    row=dict(full_mean=np.nanmean(values),pilot_mean=np.nanmean(subset),full_sd=standard_deviation,
                             standardized_wasserstein=wasserstein_distance(values[valid],subset[np.isfinite(subset)])/standard_deviation if standard_deviation>0 else 0.0,
                             full_positive_fraction=float(np.mean(full_positive)),pilot_positive_fraction=float(np.mean(pilot_positive)))
                    for q in [0.05,0.25,0.5,0.75,0.95]:
                        row['full_q%02d'%int(q*100)]=np.nanquantile(values,q)
                        row['pilot_q%02d'%int(q*100)]=np.nanquantile(subset,q)
                else:
                    row=dict(full_mean=np.nan,pilot_mean=np.nan,full_sd=np.nan,standardized_wasserstein=np.nan,
                             full_positive_fraction=np.nan,pilot_positive_fraction=np.nan)
                row.update(sample=sample,program=name,n_full=int(mask.sum()),n_pilot=int((mask&selected).sum()),
                           interpretation='balanced pilot differs in sample weighting' if sample=='ALL' else 'within-sample descriptive comparison')
                program_rows.append(row)
        pd.DataFrame(program_rows).to_csv(self.output_dir/'full_vs_pilot_program_coverage.tsv',sep='\t',index=False)
        fraction=float(bins.loc[bins['underrepresented'],'full_cells'].sum()/len(obs))
        return dict(n_bins=1600,well_populated_full_bins=int(bins['well_populated_full_bin'].sum()),
                    underrepresented_bins=int(bins['underrepresented'].sum()),fraction_full_cells_in_underrepresented_bins=fraction,
                    rule='40x40 fixed UMAP bins; full n>=100, expected balanced-pilot n>=5, observed/expected<0.5',
                    limitation='UMAP coverage is descriptive manifold coverage, not proof of exhaustive biological identities. ALL program comparison has deliberately different sample weighting.')
