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
  initially empty until the sample key is confirmed.
- `results/`: versioned major analysis packages.
- `logs/`: run-specific or linked scheduler logs.
- `jobs/`: exact submitted scheduler scripts.
- `final_figures/`: finalized or explicitly labeled candidate figure packages.

The operational handoff is copied here as `HANDOFF.md`; its canonical source is
the version-controlled file in the repository.
