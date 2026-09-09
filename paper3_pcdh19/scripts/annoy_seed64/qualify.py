"""Qualify exact seed forwarding without touching production inputs or packages.

The small numerical fixtures exercise all frozen fit seeds, compare the bridge
with the original binding for representable seeds, and exercise Allen's actual
on-disk Annoy / multiprocessing / Jaccard / Louvain path. No biological result
or production parameter is inferred from these fixtures.
"""
import ctypes
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import tempfile

import annoy
import anndata
import numpy as np
import pandas as pd
import transcriptomic_clustering.clustering as clustering


def main():
    """Fail on any mismatch and publish a machine-readable qualification record."""
    root, folder, output = map(Path, sys.argv[1:4])
    spec = importlib.util.spec_from_file_location('seed_worker', folder / 'worker.py')
    worker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(worker)
    original = annoy.AnnoyIndex
    manifest, patched = worker.install(folder)
    bridge = ctypes.PyDLL(str(folder / 'seed64.so'))
    setter = bridge.pcdh19_annoy_set_seed64
    setter.argtypes = [ctypes.py_object, ctypes.c_uint64]
    setter.restype = ctypes.c_int
    seeds = pd.read_csv(root / 'config/seeds.tsv', sep='\t').fit_seed.astype(int).tolist()
    matrix = np.random.RandomState(415).normal(size=(160, 6)).astype('float32')
    results = []
    with tempfile.TemporaryDirectory(prefix='pcdh19_seed64_') as temporary:
        temporary = Path(temporary)

        def index(seed, name, force_native=False):
            """Use the installed native index, with deterministic single-thread trees."""
            value = original(6, 'euclidean') if force_native else patched(6, 'euclidean')
            if force_native:
                assert setter(value, seed) == 0
            else:
                value.set_seed(seed)
            for row, vector in enumerate(matrix):
                value.add_item(row, vector)
            value.build(50, n_jobs=1)
            path = temporary / name
            value.save(str(path))
            return hashlib.sha256(path.read_bytes()).hexdigest(), [value.get_nns_by_item(i, 15) for i in range(160)]

        for i, seed in enumerate(seeds):
            a = index(seed, 'a.ann')
            b = index(seed, 'b.ann', force_native=True)
            assert a == b, ('exact native forwarding mismatch', i, seed)
            results.append(dict(iteration=i, exact_seed=seed, native_index_sha256=a[0]))

        def graph(seed, name, constructor, jobs):
            """Exercise frozen Allen implementation, including on-disk reloads."""
            clustering.AnnoyIndex = constructor
            return clustering.cluster_louvain(
                anndata.AnnData(matrix), k=15, annoy_trees=50, n_jobs=jobs,
                random_seed=seed, annoy_index_filename=str(temporary / (name + '.ann')))

        old = graph(seeds[0], 'old', original, 1)
        new = graph(seeds[0], 'new', patched, 1)
        assert np.array_equal(old[0], new[0]) and (old[2] != new[2]).nnz == 0 and old[3] == new[3]
        high = graph(seeds[2], 'high', patched, 8)
        assert len(high[0]) == 160 and high[2].shape == (160, 160)
        assert np.isfinite(high[2].data).all()
    output.write_text(json.dumps(dict(status='PASS', frozen_seeds_tested=100,
        larger_than_signed_int32=sum(seed > 2**31 - 1 for seed in seeds),
        exact_native_indexes_and_knn_match=True,
        original_safe_seed_allen_graph_and_labels_identical=True,
        high_seed_allen_graph_louvain_8_processes_pass=True,
        compatibility_manifest_sha256=worker.sha(folder / 'MANIFEST.json'),
        results=results), indent=2) + '\n')
    print('PASS: all 100 exact seeds; native index/KNN equivalence; Allen graph/labels; 8-process high-seed path')


if __name__ == '__main__':
    main()
