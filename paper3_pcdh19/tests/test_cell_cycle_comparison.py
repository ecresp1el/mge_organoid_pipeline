"""Numerical regression and audited gene-resolution contracts for the rebuild."""
import importlib.util
from pathlib import Path
import numpy as np

spec=importlib.util.spec_from_file_location('cycle_compare',Path(__file__).parents[1]/'scripts/compare_pcdh19_cell_cycle_umaps.py')
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)

def test_block_ols_matches_joint_least_squares():
    rng=np.random.RandomState(4);cov=rng.normal(size=(100,2));matrix=rng.normal(size=(100,17)).astype(np.float32)
    matrix+=cov[:,0,None].astype(np.float32)*2
    design=np.column_stack([np.ones(100),cov]);beta=np.linalg.lstsq(design,matrix,rcond=None)[0]
    expected=matrix-design@beta;actual=matrix.copy()
    fitted=mod.regress_blocks(actual,cov,block_genes=4)
    np.testing.assert_allclose(actual,expected,atol=1e-6)
    np.testing.assert_allclose(fitted,beta,atol=1e-10)
    np.testing.assert_allclose(design.T@actual,0,atol=1e-5)

def test_gene_resolution_rejects_ambiguous_features_and_records_aliases():
    names=['MCM'+str(i) for i in range(12)]
    cfg={'cycle_reference':{'renamed_aliases':{'OLD':'New'}},'programs':[
        {'name':phase,'genes':names+['OLD','DUP','MISSING']} for phase in ['S_phase','G2M']]}
    resolved,audit=mod.resolve_cycle_genes(names+['New','DUP','dup'],[str(i) for i in range(15)],cfg)
    assert len(resolved['S_phase'])==13
    assert any(r['requested']=='OLD' and r['status']=='matched' for r in audit)
    assert any(r['requested']=='DUP' and r['status']=='ambiguous' for r in audit)
    assert any(r['requested']=='MISSING' and r['status']=='missing' for r in audit)
