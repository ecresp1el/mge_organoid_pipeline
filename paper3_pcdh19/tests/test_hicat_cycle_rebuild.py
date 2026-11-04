"""Numerical checks for the new regression/model-mapping boundary."""
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from hicat.cycle_rebuild import fit_cycle_model,corrected_block,map_corrected
from hicat.consensus_core import map_heldout


def test_sparse_OLS_and_scaling_match_dense():
    rng=np.random.RandomState(123);x=rng.gamma(2,2,(120,30));x[x<1.5]=0
    design=np.column_stack([np.ones(len(x)),rng.normal(size=(len(x),2))])
    beta,sd=fit_cycle_model(sparse.csr_matrix(x),design)
    expected=np.linalg.lstsq(design,x,rcond=None)[0];residual=x-design@expected
    np.testing.assert_allclose(beta,expected,rtol=1e-10,atol=1e-10)
    np.testing.assert_allclose(sd,residual.std(axis=0,ddof=1),rtol=1e-10,atol=1e-10)
    columns=np.arange(0,30,2)
    actual=corrected_block(sparse.csr_matrix(x),design,beta,sd,columns)
    np.testing.assert_allclose(actual,np.minimum(residual[:,columns]/sd[columns],10),atol=1e-10)


def test_block_mapping_matches_original_on_explicit_residuals():
    rng=np.random.RandomState(321);x=rng.gamma(2,2,(80,30));x[x<1.5]=0
    design=np.column_stack([np.ones(len(x)),rng.normal(size=(len(x),2))])
    beta,sd=fit_cycle_model(sparse.csr_matrix(x),design);model=dict(design=design,beta=beta,sd=sd)
    rows=np.sort(rng.choice(80,64,replace=False));labels=np.tile(np.arange(1,5),16);columns=np.arange(25)
    obj=ad.AnnData(sparse.csr_matrix(x),obs=pd.DataFrame(index=['c%d'%i for i in range(80)]))
    corrected=corrected_block(obj.X,design,beta,sd,np.arange(30))
    expected=map_heldout(corrected,rows,labels,columns,block_cells=11)
    actual=map_corrected(obj,rows,labels,columns,model,block_cells=11)
    for a,b in zip(actual,expected):np.testing.assert_allclose(a,b,equal_nan=True,atol=1e-12)
    np.testing.assert_array_equal(actual[0][rows],labels)
