# Paper 3 Ziobro PCDH19 MGE single-cell RNA-seq: Turbo workstream

This directory is the output root for the independent Paper 3 Ziobro PCDH19
single-cell workstream. It follows the same layout and reproducibility contract
as the neighboring Paper 2 workstream.

Canonical code and configuration are version controlled at:

```text
/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline/paper3_pcdh19
```

The upstream source delivery is:

```text
/nfs/turbo/umms-ziobroj/Ziobro Lab/MGE scRNA-seq
```

That source remains on the separate Ziobro Turbo allocation and must remain
read-only. Paper 3 reads it in place while writing analysis products here under
the existing `umms-parent` MGE project.

Directory roles:

- `inputs/`: selected and validated input registrations or canonical copies;
  the canonical biological sample key is version controlled in the repository
  at `paper3_pcdh19/config/sample_key.csv`.
- `results/`: versioned major analysis packages.
- `logs/`: run-specific or linked scheduler logs.
- `jobs/`: exact submitted scheduler scripts.
- `final_figures/`: finalized or explicitly labeled candidate figure packages.

Major new workflows follow the neighboring Paper 2 run contract: default
submissions create versioned packages, exact runnable files are copied into
each package before submission, and SLURM executes those copies. A workflow may
offer an explicit guarded `--replace-run RUN_ID` option to regenerate one named
inactive run inside that same step. Such replacement does not authorize
changes to frozen steps or shared read-only/cached inputs.

The primary-processing workflow has its own Steps 00–07 namespace beneath
`results/primary_processing/`. Its exact Python, shell, SLURM, configuration,
and registered metadata are frozen into every run. Computation success leaves
a step `IN_REVIEW`; only explicit user approval promotes a checkpoint for use
by the next step. Existing mapping, probe, and classification results are
downstream-only and do not enter primary preprocessing.

Current checkpoint: Step 00 run
`00_input_validation_and_canonical_anndata_20260830_113749_d8b6bf7` is
approved. Step 01 run `01_qc_metrics_20260830_115715_2b57907` completed in
Great Lakes job `59281063`, retained all 450,788 cells and 19,071 genes, and
passed 40/40 validations without thresholds or filtering. It remains
`IN_REVIEW`; Step 02 is not authorized.

The operational handoff is copied here as `HANDOFF.md`; its canonical source is
the version-controlled file in the repository.

The complete input/output and biological-interpretation contract is maintained
in the repository at:

```text
/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline/paper3_pcdh19/PIPELINE_IO_AND_BIOLOGICAL_SCOPE.md
```

The frozen probe-audit outputs contain technical sample IDs and observed raw
probe counts only. Genotype and sex are supplied by the separate registered
sample key; neither that mapping nor an A-negative probe pattern by itself
establishes an individual cell's genotype or mutant-cell identity.

The separate `results/xgfp_probe_audit/` package validates whether the custom
Flex EGFP probe sequences match the Nagy/Kalantry D4/XEGFP construct-level
reporter sequence. It establishes sequence compatibility, not reporter
expression or individual-cell identity.
