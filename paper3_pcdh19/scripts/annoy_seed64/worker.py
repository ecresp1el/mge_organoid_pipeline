"""Run the unchanged frozen HiCAT worker with an exact uint64 Annoy seed bridge.

Frozen seeds are never truncated, wrapped, regenerated or replaced. Seeds that
fit the original Python binding use that original binding. Only larger seeds
use the native uint64 setter in the same installed Annoy binary. The bridge and
its provenance are separately frozen and included in each new attempt manifest.
Existing valid checkpoints retain their original contracts and are reused.
"""
import ctypes
import hashlib
import json
from pathlib import Path
import runpy
import sys
import annoy
import annoy.annoylib


def sha(path):
    """Hash the small runtime artifacts attested by the compatibility ledger."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def install(folder):
    """Verify the bridge/native binary, then patch only the Annoy constructor."""
    folder=Path(folder);manifest=json.loads((folder/'MANIFEST.json').read_text())
    for name,expected in manifest['files'].items():
        if sha(folder/name)!=expected:raise ValueError('Compatibility source/binary changed: '+name)
    if sha(annoy.annoylib.__file__)!=manifest['installed_annoy_binary_sha256']:
        raise ValueError('Annoy binary differs from the qualified native ABI')
    bridge=ctypes.PyDLL(str(folder/'seed64.so'))
    setter=bridge.pcdh19_annoy_set_seed64
    setter.argtypes=[ctypes.py_object,ctypes.c_uint64];setter.restype=ctypes.c_int
    original=annoy.AnnoyIndex
    class Seed64Annoy:
        """Delegate every operation to the original index; widen only seed input."""
        def __init__(self,*args,**kwargs):self.native=original(*args,**kwargs)
        def __getattr__(self,name):return getattr(self.native,name)
        def set_seed(self,seed):
            """Forward the unchanged positive seed, without modulo or sign casting."""
            seed=int(seed)
            if not 0<=seed<2**64:raise ValueError('Seed outside native uint64 domain')
            if seed<=2**31-1:return self.native.set_seed(seed)
            if setter(self.native,seed)!=0:raise RuntimeError('Native Annoy seed forwarding failed')
    import transcriptomic_clustering.clustering as clustering
    clustering.AnnoyIndex=Seed64Annoy
    return manifest,Seed64Annoy


def main():
    """Record bridge provenance inside new checkpoints, then use frozen main."""
    folder=Path(__file__).resolve().parent
    manifest,_=install(folder)
    import hicat.consensus_restart as worker
    original_checkpoint=worker.checkpoint
    def checkpoint(stage,contract,compute,validate):
        """Keep scientific contracts unchanged; seal runtime provenance as an asset."""
        def recorded_compute(attempt):
            (attempt/'annoy_seed_binding_provenance.json').write_text(json.dumps(dict(
                compatibility_manifest_sha256=sha(folder/'MANIFEST.json'),
                compatibility_manifest=str(folder/'MANIFEST.json'),
                native_annoy_binary_sha256=manifest['installed_annoy_binary_sha256'],
                exact_frozen_seed_preserved=True,seed_modulo_or_regeneration=False,
                original_scientific_worker_unchanged=True),indent=2)+'\n')
            compute(attempt)
        return original_checkpoint(stage,contract,recorded_compute,validate)
    worker.checkpoint=checkpoint
    worker.main()


if __name__=='__main__':main()
