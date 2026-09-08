# Frozen source material for the Allen consensus audit

These files were extracted using `git show COMMIT:PATH`; source repositories and commits are in `reference_manifest.json`. Original licenses are retained with each source snapshot.

`allen_r/` is the current audited consensus reference; `allen_r_2018/` establishes what the historical sparse collector actually did; `allen_python/` records the existing fitting backend and its scaling limits. These snapshots are for inspection, not an executable new pipeline or benchmark.

`source_object_metadata.json` records a read-only HDF5 metadata inspection. Its expected count-object checksum came from the existing approved-input contract and was not recomputed in this metadata-only audit. The benchmark must recompute it before fitting.

`proposed_method_contract.json` is explicitly proposed. No scientific output or new cluster labels are saved here.
