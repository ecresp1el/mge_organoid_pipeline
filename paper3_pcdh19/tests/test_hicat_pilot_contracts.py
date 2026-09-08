"""Small synthetic checks for the pilot's partition and upstream-adapter contracts.

Run with frozen/source upstream on PYTHONPATH in the isolated Python 3.8 env.
These exercise numerical code but are not evidence about PCDH19 cluster count.
"""
import json
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import transcriptomic_clustering as tc
from hicat.engine import PilotEngine
from hicat.provenance import Progress, write_json


class PilotContracts(unittest.TestCase):
    """Check that adapter outputs preserve identity and make failures explicit."""

    def test_uns_inventory_preserves_mapping(self):
        """AnnData's uns wrapper must export structured metadata, not repr text."""
        obj=ad.AnnData(np.ones((2,2)),dtype=np.float64)
        obj.uns['hicat']={'pilot_cells':np.int64(2),'annotation_performed':False}
        with tempfile.TemporaryDirectory(prefix='hicat_uns_',dir='/tmp') as tmp:
            path=Path(tmp)/'inventory.json'
            write_json(path,obj.uns)
            self.assertEqual(json.loads(path.read_text()),dict(obj.uns))

    def test_partition_rejects_duplicates_and_missing(self):
        """An invalid membership must raise instead of writing incomplete labels."""
        with self.assertRaises(ValueError):
            PilotEngine.labels([np.array([0,1]),np.array([1,2])],4)

    def test_normalization_formula(self):
        """The actual upstream formula must be natural log of all-gene CPM+1."""
        x=np.array([[1.,2.,0.],[2.,1.,1.]],dtype=np.float64)
        obj=ad.AnnData(x.copy(),dtype=np.float64)
        out=tc.normalize(obj)
        np.testing.assert_allclose(out.X,np.log1p(x/x.sum(axis=1)[:,None]*1e6))
        np.testing.assert_array_equal(obj.X,x)

    def test_array_memberships_can_actually_merge(self):
        """Two identical-distribution groups supplied as arrays must merge via the adapter."""
        cfg=json.loads((Path(__file__).parents[1]/'config/hicat_pilot.json').read_text())
        rng=np.random.default_rng(7)
        counts=rng.poisson(2,size=(30,40)).astype(float)
        matrix=np.vstack([counts,counts])
        norm=tc.normalize(ad.AnnData(matrix,dtype=np.float64))
        with tempfile.TemporaryDirectory(prefix='hicat_merge_',dir='/tmp') as tmp:
            directory=Path(tmp)/'candidate'
            engine=PilotEngine(cfg,dict(name='merge_test',seed=1,score_thresh=150),directory,Progress(Path(tmp)/'provenance'))
            groups={0:np.arange(30),1:np.arange(30,60)}
            merged=engine._merge(norm,norm[:,:5].copy(),groups,np.repeat([0,1],30),directory)
            self.assertEqual(len(merged),1)
            self.assertEqual(sorted(next(iter(merged.values()))),list(range(60)))

    def test_synthetic_engine_preserves_cells_and_writes_models(self):
        """Run separable synthetic data through recursion/final audit and persistence."""
        cfg=json.loads((Path(__file__).parents[1]/'config/hicat_pilot.json').read_text())
        cfg.update(max_hvgs=60,pca_components=10,max_retained_pcs=6,graph_k=5,min_recursive_cells=30)
        cfg['thresholds'].update(cluster_size_thresh=10,score_thresh=40)
        rng=np.random.default_rng(12)
        rates=rng.uniform(.001,.05,size=(180,80))
        for group in range(3):rates[group*60:(group+1)*60,group*15:(group+1)*15]+=15
        counts=rng.poisson(rates).astype(np.float64)
        raw=ad.AnnData(counts,obs=pd.DataFrame(index=['c%d'%i for i in range(180)]),
                       var=pd.DataFrame(index=['g%d'%i for i in range(80)]),dtype=np.float64)
        norm=tc.normalize(raw)
        saved=norm.X.copy()
        with tempfile.TemporaryDirectory(prefix='hicat_contract_',dir='/tmp') as tmp:
            engine=PilotEngine(cfg,dict(name='synthetic',seed=123,score_thresh=40),Path(tmp)/'candidate',Progress(Path(tmp)/'provenance'))
            labels,summary,evidence,means=engine.run(norm)
            self.assertEqual(len(labels),180)
            self.assertGreater(len(set(labels)),1)
            self.assertTrue((Path(tmp)/'candidate/root/pca_components_retained.tsv').is_file())
            self.assertTrue((Path(tmp)/'candidate/final_pairwise_evidence.tsv').is_file())
            self.assertFalse(summary['eligible_for_full_dataset_count'])
        np.testing.assert_array_equal(norm.X,saved)


if __name__=='__main__':
    unittest.main()
