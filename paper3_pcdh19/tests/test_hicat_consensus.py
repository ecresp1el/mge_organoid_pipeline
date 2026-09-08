"""Independent numerical qualification for the full-data benchmark adapters.

These fixtures are synthetic numerical checks, never a biological pilot or a
replacement for the real 357,079-cell benchmark. The R oracle executes pinned
Allen functions. The sparse fitter check compares the existing dense adapter
against the memory-bounded one using identical expression and settings.
"""
from pathlib import Path
import json
import os
import subprocess
import tempfile
import unittest
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from scipy.stats import pearsonr
from sklearn.metrics import adjusted_rand_score
import transcriptomic_clustering as tc
from hicat.consensus_core import Membership, map_heldout
from hicat.fullscale_engine import FullScaleEngine, sparse_cluster_means
from hicat.engine import PilotEngine
from hicat.provenance import Progress, sha256, write_json
from hicat.consensus_benchmark import Benchmark


class ConsensusQualification(unittest.TestCase):
    """Check source fidelity, label invariance, mapping and fitting numerics."""

    def test_factorization_and_label_permutation(self):
        """Explicit pair frequencies equal factorized affinities; IDs arbitrary."""
        rng=np.random.default_rng(21)
        iterations=[rng.integers(0,5,77) for _ in range(9)]
        labels=rng.integers(0,4,77)
        membership=Membership(iterations,block_cells=11)
        pair=sum(x[:,None]==x[None,:] for x in iterations)/len(iterations)
        a,keys=membership.affinity(labels)
        expected=np.column_stack([pair[:,labels==k].mean(axis=1) for k in keys])
        np.testing.assert_allclose(a,expected,rtol=1e-14,atol=1e-14)
        permuted=[np.array([71,5,2,99,13])[x] for x in iterations]
        observed,_=Membership(permuted).affinity(labels)
        np.testing.assert_allclose(a,observed,rtol=1e-14,atol=1e-14)
        rows=np.arange(0,77,3)
        subset,keys=membership.affinity(labels[rows],target_rows=rows)
        np.testing.assert_allclose(subset,np.column_stack([pair[:,rows[labels[rows]==k]].mean(axis=1) for k in keys]))

    def test_missing_iteration_cells_rejected(self):
        """A missing held-out assignment cannot change the denominator silently."""
        with self.assertRaises(ValueError):Membership([np.array([1,1]),np.array([1])])
        with self.assertRaises(ValueError):Membership([np.array([1,-1])])

    def test_heldout_pearson(self):
        """Assignments and scores match direct SciPy Pearson with explicit means."""
        rng=np.random.default_rng(33);x=rng.uniform(0,9,(40,20));x[x<5]=0
        sampled=np.arange(30);labels=np.repeat([1,2,10],10);markers=np.arange(15)
        observed,scores,prototypes,keys=map_heldout(sparse.csr_matrix(x),sampled,labels,markers,7)
        expected=np.array([[pearsonr(x[row,markers],x[sampled[labels==key]][:,markers].mean(axis=0))[0]
                            for key in keys] for row in range(30,40)])
        np.testing.assert_array_equal(observed[:30],labels)
        np.testing.assert_array_equal(observed[30:],keys[expected.argmax(axis=1)])
        np.testing.assert_allclose(scores[30:],expected.max(axis=1),atol=1e-14)
        self.assertTrue(np.isnan(scores[:30]).all())
        x[30]=0
        observed,scores,_,_=map_heldout(sparse.csr_matrix(x),sampled,labels,markers)
        self.assertEqual(observed[30],1);self.assertEqual(scores[30],0)

    def test_sparse_statistics_and_DE(self):
        """Sparse ddof=1/detection and eBayes evidence agree with dense upstream."""
        rng=np.random.default_rng(56)
        rates=np.full((120,50),.03)
        rates[:40,:10]+=8;rates[40:80,10:20]+=8;rates[80:,20:30]+=8
        x=np.log1p(rng.poisson(rates).astype(float)*100)
        dense=ad.AnnData(x,dtype=np.float64);csr=ad.AnnData(sparse.csr_matrix(x),dtype=np.float64)
        groups={k:np.arange(k*40,(k+1)*40) for k in range(3)};labels=np.repeat(np.arange(3),40)
        original=tc.get_cluster_means(dense,groups,labels,low_th=1)
        observed=sparse_cluster_means(csr,groups,labels,low_th=1)
        for a,b in zip(original,observed):np.testing.assert_allclose(a,b,atol=1e-12,rtol=1e-12)
        thresholds=dict(q1_thresh=.4,q2_thresh=None,cluster_size_thresh=20,qdiff_thresh=.7,padj_thresh=.05,lfc_thresh=1)
        def evidence(values):
            """Use unchanged upstream eBayes on independently computed statistics."""
            means,detection,variance=values
            return tc.de_pairs_ebayes([(0,1),(0,2),(1,2)],means,variance,detection,{k:40 for k in groups},thresholds)
        first,second=evidence(original),evidence(observed)
        self.assertTrue((first.score.astype(float)>150).all())
        self.assertTrue((first.num.astype(int)>=5).all())
        np.testing.assert_allclose(first.score.astype(float),second.score.astype(float),atol=1e-8,rtol=1e-10)
        np.testing.assert_array_equal(first.num,second.num)
        for col in ['up_genes','down_genes']:
            self.assertEqual(first[col].map(list).to_list(),second[col].map(list).to_list())

    def test_pinned_R_oracle(self):
        """Execute Allen compile/affinity/merge/refine and mean-Pearson mapping."""
        reference=Path(__file__).resolve().parents[1]/'references/allen_consensus_audit_20260908/allen_r/R'
        rscript='/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/Rscript'
        rng=np.random.default_rng(112)
        base=np.repeat([1,2,10],32)
        iterations=[]
        for repeat in range(6):
            value=base.copy();rows=rng.choice(len(base),8,replace=False);value[rows]=rng.choice([1,2,10],len(rows))
            iterations.append(value)
        initial=base.copy();initial[[0,1,33,65]]=[2,10,10,2]
        membership=Membership(iterations)
        with tempfile.TemporaryDirectory(prefix='hicat_R_oracle_') as tmp:
            root=Path(tmp)
            np.savetxt(root/'iterations.tsv',np.column_stack(iterations),fmt='%d',delimiter='\t')
            np.savetxt(root/'initial.tsv',initial,fmt='%d')
            x=rng.normal(size=(96,15));np.savetxt(root/'expression.tsv',x,delimiter='\t')
            script=root/'oracle.R'
            script.write_text('''args=commandArgs(TRUE)
library(Matrix); library(matrixStats)
source(file.path(args[1],"util.R"));source(file.path(args[1],"consensusCluster.R"));source(file.path(args[1],"annotate.R"))
setwd(args[2]); values=as.matrix(read.table("iterations.tsv")); ids=paste0("c",seq_len(nrow(values)))
cl.list=lapply(seq_len(ncol(values)),function(j)setNames(values[,j],ids)); names(cl.list)=paste0("i",seq_along(cl.list))
b=compile_cl_mat(cl.list,ids); cl=setNames(scan("initial.tsv",quiet=TRUE),ids)
write.table(as.matrix(b),"B.tsv",sep="\\t",col.names=FALSE,row.names=FALSE)
write.table(get_cell.cl.co.ratio(cl,cl.mat=b),"affinity.tsv",sep="\\t",col.names=FALSE,row.names=FALSE)
merged=merge_cl_by_co(cl,cl.mat=b)
write.table(merged,"merged.tsv",col.names=FALSE,row.names=FALSE)
split=setNames(values[,1]*10+rep(rep(c(0,1),each=16),3),ids)
write.table(split,"split.tsv",col.names=FALSE,row.names=FALSE)
write.table(merge_cl_by_co(split,cl.mat=b),"split_merged.tsv",col.names=FALSE,row.names=FALSE,quote=FALSE)
refined=refine_cl(cl,cl.mat=b,min.cells=4,tol.th=.01,confusion.th=.6)
write.table(refined$cl[ids],"refined.tsv",col.names=FALSE,row.names=FALSE,quote=FALSE)
stats=get_cl_co_stats(cl,cl.mat=b)$cell.co.stats
write.table(stats[ids,],"stats.tsv",sep="\\t",col.names=FALSE,row.names=FALSE)
x=t(as.matrix(read.table("expression.tsv")));colnames(x)=ids
mapping=map_by_cor(x[,1:75],cl[1:75],x[,76:96],method="mean")$pred.df
write.table(mapping,"mapping.tsv",sep="\\t",col.names=TRUE,row.names=FALSE)
''')
            subprocess.run([rscript,str(script),str(reference),str(root)],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            np.testing.assert_array_equal(membership.b.toarray(),np.loadtxt(root/'B.tsv'))
            a,_=membership.affinity(initial);np.testing.assert_allclose(a,np.loadtxt(root/'affinity.tsv'),rtol=1e-12,atol=1e-12)
            merged,_=membership.merge_by_co(initial);np.testing.assert_array_equal(merged,np.loadtxt(root/'merged.tsv'))
            split=np.loadtxt(root/'split.tsv',dtype=int)
            merged,actions=membership.merge_by_co(split)
            self.assertGreater(len(actions),0)
            np.testing.assert_array_equal(merged,np.loadtxt(root/'split_merged.tsv'))
            refined,_,_,_=membership.refine(initial,min_cells=4);np.testing.assert_array_equal(refined,np.loadtxt(root/'refined.tsv'))
            cells,_=membership.co_stats(initial);np.testing.assert_allclose(cells,np.loadtxt(root/'stats.tsv'),rtol=1e-12,atol=1e-12)
            mapped,scores,_,_=map_heldout(x,np.arange(75),initial[:75],np.arange(15))
            expected=pd.read_csv(root/'mapping.tsv',sep='\t')
            np.testing.assert_array_equal(mapped[75:],expected['pred.cl'])
            np.testing.assert_allclose(scores[75:],expected['pred.score'],rtol=1e-12,atol=1e-12)

    def test_sparse_full_fitter_matches_existing(self):
        """Same synthetic expression yields the existing dense partition and tree."""
        cfg=json.loads((Path(__file__).parents[1]/'config/hicat_coarse_fine.json').read_text())
        cfg['thresholds']=dict(cfg['fine_thresholds'],q1_thresh=.4,qdiff_thresh=.7)
        cfg.update(max_hvgs=60,pca_components=10,max_retained_pcs=6,graph_k=5)
        rng=np.random.default_rng(82);rates=np.full((180,80),.2)
        for k in range(3):rates[k*60:(k+1)*60,k*20:(k+1)*20]+=8
        counts=rng.poisson(rates).astype(float)
        x=np.log1p(counts/counts.sum(axis=1)[:,None]*1e6)
        candidate=dict(name='qualification',seed=20260908,score_thresh=150)
        with tempfile.TemporaryDirectory(prefix='hicat_fit_qualification_') as tmp:
            root=Path(tmp)
            first=PilotEngine(cfg,candidate,root/'dense',Progress(root/'dense_progress'))
            second=FullScaleEngine(cfg,candidate,root/'sparse',Progress(root/'sparse_progress'))
            a=first.run(ad.AnnData(x,dtype=np.float64))[0]
            b=second.run(ad.AnnData(sparse.csr_matrix(x),dtype=np.float64))[0]
            self.assertEqual(adjusted_rand_score(a,b),1)
            self.assertEqual([(r['cells'],r['post_merge_clusters']) for r in first.nodes],
                             [(r['cells'],r['post_merge_clusters']) for r in second.nodes])

    def test_blockwise_prepare_saved_object(self):
        """Reopened block-written AnnData equals direct normalization and IDs."""
        rng=np.random.default_rng(52);counts=rng.poisson(2,(120,35)).astype(np.int32)
        obs=pd.DataFrame({'technical_sample_id':np.tile(['S%d'%i for i in range(12)],10)},index=['c%d'%i for i in range(120)])
        var=pd.DataFrame({'gene_symbol':['g%d'%i for i in range(35)]})
        with tempfile.TemporaryDirectory(prefix='hicat_prepare_qualification_') as tmp:
            root=Path(tmp)
            for sub in ['config','inputs','source/objects']:(root/sub).mkdir(parents=True)
            source=root/'source/objects/counts.h5ad'
            ad.AnnData(sparse.csr_matrix(counts),obs=obs,var=var,dtype=np.int32).write_h5ad(source)
            pd.DataFrame([dict(status='APPROVED',run_id='fixture')]).to_csv(root/'source/STEP_STATUS.tsv',sep='\t',index=False)
            cfg=dict(source_path=str(source),step02_run_id='fixture',step02_sha256=sha256(source),
                     expected_cells=120,expected_genes=35,expected_samples=12,block_cells=17)
            write_json(root/'config/benchmark.json',cfg);write_json(root/'frozen_manifest.json',{})
            Benchmark(root,'prepare').run()
            reopened=ad.read_h5ad(root/'inputs/full_log1p_cpm.h5ad')
            expected=np.log1p(counts/counts.sum(axis=1)[:,None]*1e6)
            np.testing.assert_allclose(reopened.X.toarray(),expected,atol=1e-12,rtol=1e-12)
            np.testing.assert_array_equal(reopened.obs_names,obs.index)
            self.assertEqual(dict(reopened.uns),json.loads((root/'inputs/anndata_uns_inventory.json').read_text()))


if __name__=='__main__':unittest.main()
