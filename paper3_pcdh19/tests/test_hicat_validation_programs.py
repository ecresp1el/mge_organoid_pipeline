"""Independent reference and streamed-source checks for Step07 score assets."""
import json
import hashlib
import h5py
from pathlib import Path
import tempfile
import unittest
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from hicat.validation_programs import CanonicalProgramScorer, normalize_log1p_cpm
from hicat.validation_full import FullDataProjector, read_full_metadata


class ProgramScoringContracts(unittest.TestCase):
    """Protect reference compatibility, missing-assay NA and pilot/full alignment."""

    def fixture(self):
        """Return count data and a small signature configuration with missing genes."""
        rng=np.random.RandomState(29)
        symbols=['Gene%d'%i for i in range(300)]
        raw=sparse.csr_matrix(rng.poisson(np.linspace(.02,2,300),(80,300)).astype(np.int32))
        obs=pd.DataFrame({'technical_sample_id':['a','b']*40},index=['cell%d'%i for i in range(80)])
        var=pd.DataFrame({'gene_symbol':symbols},index=['feature%d'%i for i in range(300)])
        obj=ad.AnnData(raw,obs=obs,var=var,dtype=raw.dtype)
        obj.obsm['X_umap']=rng.normal(size=(80,2)).astype(np.float64)
        cfg=dict(scoring=dict(normalization_target=1000000,n_bins=25,control_size=12,random_seed=17,
                             min_available_genes=3,positive_score_threshold=.1,min_detected_genes_positive=2),
                 programs=[dict(name='S_phase',genes=symbols[3:8],group='cycle'),
                           dict(name='G2M',genes=symbols[12:17],group='cycle'),
                           dict(name='identity',genes=symbols[4:6]+symbols[30:35],group='region'),
                           dict(name='ribosomal',genes=['Rpl3','Rps2'],group='quality')],
                 cycle_reference={'renamed_aliases':{}},identity_cycle_exclusions=symbols[3:8]+symbols[12:17],
                 dotplot_groups={'canonical':symbols[3:8]+symbols[30:35]})
        return obj,cfg

    def test_scanpy_cell_cycle_reference_and_missing_program(self):
        """Our fixed coefficients agree with actual Scanpy, not a mirrored test formula."""
        obj,cfg=self.fixture()
        obj.layers['log1p_cpm']=normalize_log1p_cpm(obj.X)
        with tempfile.TemporaryDirectory() as directory:
            scorer=CanonicalProgramScorer(cfg)
            result=scorer.fit_transform(obj,directory)
            reference=ad.AnnData(obj.layers['log1p_cpm'].copy(),obs=obj.obs.copy(),var=obj.var.copy(),dtype=np.float64)
            sc.tl.score_genes_cell_cycle(reference,s_genes=list(obj.var_names[3:8]),g2m_genes=list(obj.var_names[12:17]),
                                        n_bins=25,random_state=17,use_raw=False)
            np.testing.assert_allclose(result['scores']['S_phase'],reference.obs['S_score'],atol=1e-12)
            np.testing.assert_allclose(result['scores']['G2M'],reference.obs['G2M_score'],atol=1e-12)
            expected_phase=reference.obs['phase'].astype(str).replace({'G1':'G1-like','G2M':'G2/M'})
            self.assertEqual(result['phase'].tolist(),expected_phase.tolist())
            self.assertTrue(result['scores']['ribosomal'].isna().all())
            self.assertTrue(result['positive']['ribosomal'].isna().all())
            self.assertTrue(result['identity_scores']['S_phase'].isna().all())
            excluded=set(scorer.model['identity_cycle_exclusion_feature_ids'])
            for model in scorer.model['identity_programs']:
                self.assertFalse(excluded.intersection(model['control_feature_ids']))
                self.assertFalse(excluded.intersection(model['signal_feature_ids']))
            restored=np.zeros_like(scorer.weights)
            weights=pd.read_csv(Path(directory)/'score_coefficients.tsv',sep='\t')
            for row in weights[weights['model']=='standard'].itertuples():
                restored[list(obj.var_names).index(row.feature_id),scorer.names.index(row.program)]=row.weight
            np.testing.assert_allclose(restored,scorer.weights,atol=1e-15)

    def test_streamed_full_source_matches_pilot_and_saves_no_labels(self):
        """Independent chunk normalization must reproduce pilot scores and expression."""
        full,cfg=self.fixture()
        with tempfile.TemporaryDirectory() as directory:
            directory=Path(directory)
            source=directory/'raw_full.h5ad'
            full.write_h5ad(source)
            pilot=full[:40].copy()
            pilot.obsm['X_umap_step06_display']=pilot.obsm['X_umap'].copy()
            pilot.layers['log1p_cpm']=normalize_log1p_cpm(pilot.X)
            scorer=CanonicalProgramScorer(cfg)
            scorer.fit_transform(pilot,directory/'programs')
            summary=FullDataProjector(source,directory/'full',chunk_size=13).run(pilot,scorer)
            self.assertLess(summary['pilot_score_max_abs_difference'],2e-5)
            self.assertEqual(summary['n_cells'],80)
            saved=ad.read_h5ad(summary['object_path'])
            self.assertEqual(saved.shape,(80,10))
            self.assertEqual(int(saved.obs['is_pilot'].sum()),40)
            self.assertIsNone(saved.raw)
            self.assertEqual(len(saved.layers),0)
            self.assertFalse(any('hicat' in key for key in saved.obs))
            self.assertEqual(len(saved.obsp),0)
            np.testing.assert_array_equal(saved.obsm['X_umap_step06_display'],full.obsm['X_umap'])
            self.assertEqual(saved.obsm['X_umap_step06_display'].dtype,np.float64)
            self.assertTrue(np.isnan(saved.obsm['X_validation_program_scores'][:,-1]).all())
            original=ad.read_h5ad(source)
            self.assertEqual(len(original.uns),0)
            self.assertEqual((original.X!=full.X).nnz,0)
            pilot.obs['technical_sample_id']=pilot.obs['technical_sample_id'].astype(str)
            pilot.obs.iloc[0,pilot.obs.columns.get_loc('technical_sample_id')]='wrong_sample'
            with self.assertRaisesRegex(ValueError,'technical sample identity'):
                FullDataProjector(source,directory/'mismatch',chunk_size=13).run(pilot,scorer)


    def test_modern_null_uns_is_ignored_without_changing_source(self):
        """Selective metadata reading survives modern uns unsupported by old Allen env."""
        full,cfg=self.fixture()
        full.obs['submitted_sample_name']=['sample α','sample β']*40
        full.obs['n_genes_by_counts']=pd.array([None]+list(range(1,80)),dtype='Int64')
        with tempfile.TemporaryDirectory() as directory:
            directory=Path(directory)
            source=directory/'modern_uns_full.h5ad'
            full.write_h5ad(source)
            with h5py.File(source,'r+') as handle:
                marker=handle['uns'].create_dataset('modern_null',data=h5py.Empty('f'))
                marker.attrs['encoding-type']='null'
                marker.attrs['encoding-version']='0.1.0'
            before=hashlib.sha256(source.read_bytes()).hexdigest()
            obs,var,coords,shape=read_full_metadata(source)
            self.assertEqual(shape,full.shape)
            self.assertEqual(obs['submitted_sample_name'].astype(str).iloc[0],'sample α')
            self.assertTrue(pd.isna(obs['n_genes_by_counts'].iloc[0]))
            self.assertEqual(obs['n_genes_by_counts'].iloc[1],1)
            np.testing.assert_array_equal(coords,full.obsm['X_umap'])
            # Reproduce the actual compatibility failure in the pinned 0.8 env.
            if tuple(int(v) for v in ad.__version__.split('.')[:2]) < (0,10):
                with self.assertRaises(Exception) as error:
                    ad.read_h5ad(source,backed='r')
                self.assertIn('null',str(error.exception))
            pilot=full[:40].copy()
            pilot.obsm['X_umap_step06_display']=pilot.obsm['X_umap'].copy()
            pilot.layers['log1p_cpm']=normalize_log1p_cpm(pilot.X)
            scorer=CanonicalProgramScorer(cfg)
            scorer.fit_transform(pilot,directory/'programs')
            summary=FullDataProjector(source,directory/'full',chunk_size=13).run(pilot,scorer)
            self.assertEqual(summary['source_metadata_reader'],'selective_h5py_without_uns')
            self.assertEqual(hashlib.sha256(source.read_bytes()).hexdigest(),before)
            self.assertLess(summary['pilot_score_max_abs_difference'],2e-5)


if __name__=='__main__':
    unittest.main()
