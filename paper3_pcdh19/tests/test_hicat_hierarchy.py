"""Numerical marker checks and end-to-end synthetic hierarchy export validation."""
from pathlib import Path
import json
import tempfile
import unittest
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from scipy.stats import ttest_ind
from hicat.hierarchy import HierarchyEngine
from hicat.hierarchy_workflow import HierarchyWorkflow
from hicat.markers import MarkerReporter
from hicat.provenance import Progress, write_json, sha256


class HierarchyContracts(unittest.TestCase):
    """Exercise invalid hierarchy, marker statistics, and real saved workflow assets."""

    def test_nesting_rejects_cross_parent_cluster(self):
        """A fine label spanning coarse parents must fail explicitly."""
        with self.assertRaises(ValueError):
            HierarchyEngine.validate_nesting(['A','B'],['F1','F1'])

    def test_marker_statistics_match_direct_welch(self):
        """Saved p-values match direct tests, absent DE is not padded to 20 genes."""
        rng=np.random.default_rng(1)
        x=rng.poisson(1,size=(40,8)).astype(float)
        x[:20,0]+=8
        labels=np.repeat(['A','B'],20)
        var=pd.DataFrame({'gene_symbol':['gene%d'%i for i in range(8)]},index=['g%d'%i for i in range(8)])
        cfg=dict(top_n=20,max_fdr=.05,min_delta_mean_log1p_cpm=.25,min_fraction_detected=.25)
        with tempfile.TemporaryDirectory(prefix='hicat_markers_',dir='/tmp') as tmp:
            root=Path(tmp)
            reporter=MarkerReporter(cfg,root/'markers',Progress(root/'progress'))
            _,_,top=reporter.run(sparse.csr_matrix(x),labels,var,'synthetic')
            table=pd.read_csv(root/'markers/A_all_genes.tsv.gz',sep='\t').set_index('gene_id').loc[var.index]
            expected=ttest_ind(x[:20],x[20:],axis=0,equal_var=False)
            np.testing.assert_allclose(table.p_value,expected.pvalue,rtol=1e-10)
            np.testing.assert_allclose(table.welch_t,expected.statistic,rtol=1e-10)
            self.assertEqual(top[top.cluster=='A'].iloc[0].gene_id,'g0')
            self.assertLess(len(top),40)
        np.testing.assert_allclose(MarkerReporter.bh([.01,.04,.03]),[.03,.04,.04])

    def test_end_to_end_saved_hierarchy(self):
        """Run a real synthetic coarse/fine fit and verify models, plots and marker scope."""
        cfg=json.loads((Path(__file__).parents[1]/'config/hicat_coarse_fine.json').read_text())
        cfg.update(expected_cells=180,expected_genes=80,expected_samples=12,pilot_cells_per_sample=15,
                   max_hvgs=60,pca_components=10,max_retained_pcs=6,graph_k=5,min_recursive_cells=30,report_dpi=50)
        for key,score in [('coarse_thresholds',80),('fine_thresholds',40)]:
            cfg[key].update(cluster_size_thresh=10,score_thresh=score)
        rng=np.random.default_rng(12)
        rates=rng.uniform(.001,.05,size=(180,80))
        for group in range(3):rates[group*60:(group+1)*60,group*15:(group+1)*15]+=15
        counts=sparse.csr_matrix(rng.poisson(rates).astype(np.int32))
        raw=ad.AnnData(counts,obs=pd.DataFrame({'technical_sample_id':np.tile(['S%d'%i for i in range(12)],15)},
                         index=['cell%d'%i for i in range(180)]),
                         var=pd.DataFrame({'gene_symbol':['Gene%d'%i for i in range(80)]},index=['g%d'%i for i in range(80)]),dtype=np.int32)
        raw.obsm['X_umap_step06_display']=rng.normal(size=(180,2))
        raw.obsm['X_pca_step06_display']=rng.normal(size=(180,5))
        with tempfile.TemporaryDirectory(prefix='hicat_hierarchy_',dir='/tmp') as tmp:
            root=Path(tmp)
            for name in ['config','inputs','logs','provenance']:(root/name).mkdir()
            write_json(root/'config/hicat_pilot.json',cfg)
            path=root/'inputs/pilot_raw_counts.h5ad'
            raw.write_h5ad(path)
            write_json(root/'inputs/input_identity.json',dict(input_pilot_sha256=sha256(path)))
            output=HierarchyWorkflow(root).run()
            saved=ad.read_h5ad(output/'pcdh19_hicat_coarse_fine.h5ad')
            mapping=pd.read_csv(output/'parent_child_mapping.tsv',sep='\t')
            self.assertEqual(mapping.cells.sum(),180)
            self.assertGreater(mapping.coarse_cluster.nunique(),1)
            for candidate in cfg['runs']:
                name=candidate['name']
                HierarchyEngine.validate_nesting(saved.obs['hicat_coarse_'+name],saved.obs['hicat_fine_'+name])
            self.assertTrue((output/'figures/coarse_fine_review.pdf').is_file())
            self.assertTrue((output/'markers/coarse_global/expression_linkage.npy').is_file())
            tree=pd.read_csv(output/'baseline/coarse/recursion_tree.tsv',sep='\t')
            self.assertTrue(tree[tree.depth>0].terminal_reason.eq('coarse_one_step_boundary').all())
            self.assertEqual((raw.X!=saved.X).nnz,0)
            self.assertEqual(json.loads((output/'anndata_uns_inventory.json').read_text()),dict(saved.uns))
            self.assertTrue(pd.read_csv(output/'validation_checks.tsv',sep='\t').status.eq('PASS').all())


if __name__=='__main__':
    unittest.main()
