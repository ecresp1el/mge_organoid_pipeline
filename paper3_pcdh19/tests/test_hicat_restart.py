"""Failure/reuse and disk-consumption contracts for independent iterations."""
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from scipy import sparse
from hicat.checkpoints import checkpoint,complete_attempt
from hicat.consensus_core import Membership
from hicat.provenance import sha256


class RestartQualification(unittest.TestCase):
    """Check checkpoint isolation and exact aggregation of saved sparse blocks."""

    def test_failed_mapping_preserves_completed_fit(self):
        """Retry only a failed mapping attempt; fit call count and hashes stay fixed."""
        with tempfile.TemporaryDirectory(prefix='hicat_restart_test_') as tmp:
            root=Path(tmp);calls=dict(fit=0,mapping=0)
            def fit(attempt):
                """Minimal numerical fit fixture, completed once."""
                calls['fit']+=1;np.save(attempt/'labels.npy',np.array([1,1,2,2]))
            def verify(attempt):
                """Read the actual saved arrays, not cached writer objects."""
                self.assertTrue((attempt/'labels.npy').is_file())
                return dict(cells=len(np.load(attempt/'labels.npy')))
            fitted=checkpoint(root/'fit',{'seed':71},fit,verify)
            before=sha256(fitted/'labels.npy');seal=sha256(fitted/'ARTIFACTS.json')
            def bad_mapping(attempt):
                """Inject failure after writing output but before publishing a seal."""
                calls['mapping']+=1;complete_attempt(root/'fit')
                np.save(attempt/'labels.npy',np.array([1,1,2,2,1]))
                raise RuntimeError('injected mapping failure')
            with self.assertRaises(RuntimeError):checkpoint(root/'mapping',{'fit':seal},bad_mapping,verify)
            self.assertFalse((root/'mapping/CURRENT.json').exists())
            self.assertEqual(len(list((root/'mapping/attempts').glob('*/FAILURE.json'))),1)
            def good_mapping(attempt):
                """Recovery consumes the completed fit and writes a new attempt."""
                calls['mapping']+=1;complete_attempt(root/'fit')
                np.save(attempt/'labels.npy',np.array([1,1,2,2,1]))
            result=checkpoint(root/'mapping',{'fit':seal},good_mapping,verify)
            again=checkpoint(root/'fit',{'seed':71},fit,verify)
            same=checkpoint(root/'mapping',{'fit':seal},good_mapping,verify)
            self.assertEqual(fitted,again);self.assertEqual(result,same)
            self.assertEqual(calls,dict(fit=1,mapping=2));self.assertEqual(sha256(fitted/'labels.npy'),before)
            self.assertEqual(sha256(fitted/'ARTIFACTS.json'),seal)
            np.save(result/'labels.npy',np.array([99]))
            with self.assertRaises(ValueError):complete_attempt(root/'mapping')
            complete_attempt(root/'fit')

    def test_saved_iteration_memberships_combine_exactly(self):
        """The aggregator consumes reopened files, not an in-memory surrogate."""
        with tempfile.TemporaryDirectory(prefix='hicat_membership_test_') as tmp:
            root=Path(tmp);labels=[np.array([1,1,2,2,9,9]),np.array([8,8,5,5,5,3])]
            blocks=[]
            for index,value in enumerate(labels):
                path=root/('iteration%d.npz'%index);sparse.save_npz(path,Membership([value]).b)
                blocks.append(sparse.load_npz(path))
            combined=Membership.from_blocks(blocks)
            oracle=Membership(labels)
            self.assertEqual(combined.b.shape,oracle.b.shape)
            self.assertEqual((combined.b!=oracle.b).nnz,0)
            observed,_=combined.affinity(np.array([1,1,2,2,3,3]))
            expected,_=oracle.affinity(np.array([1,1,2,2,3,3]))
            np.testing.assert_array_equal(observed,expected)
            with self.assertRaises(ValueError):Membership.from_blocks([blocks[0][:,:5],blocks[1]])
            corrupted=blocks[0].copy();corrupted.data[0]=0
            with self.assertRaises(ValueError):Membership.from_blocks([corrupted])


if __name__=='__main__':unittest.main()
