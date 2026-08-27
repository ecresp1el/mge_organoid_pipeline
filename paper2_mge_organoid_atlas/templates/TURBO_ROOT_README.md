# Paper 2 MGE organoid atlas: Turbo workstream

This directory is the output root for the Paper 2 cross-study in-vitro MGE
organoid atlas.

Canonical code and configuration are version controlled at:

```text
/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline/paper2_mge_organoid_atlas
```

Directory roles:

- `results/`: immutable, versioned major analysis packages.
- `logs/`: convenient links to run-specific SLURM logs.
- `jobs/`: exact timestamped `.sbatch` files submitted to SLURM.
- `final_figures/`: finalized or explicitly labeled candidate figure packages.

Each major package must contain exact code/configuration, input manifests and
checksums, commands, SLURM information/logs, provenance, validation tables,
figures when applicable, a README, and a terminal `SUCCESS.txt` or
`FAILED.txt` marker.

The operational handoff is copied into this directory as `HANDOFF.md`. Its
canonical version-controlled source is:

```text
/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline/paper2_mge_organoid_atlas/HANDOFF.md
```
