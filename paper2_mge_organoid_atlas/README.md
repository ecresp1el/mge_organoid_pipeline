# Paper 2: cross-study MGE organoid atlas

This workstream builds a cross-study single-cell atlas from six processed
in-vitro MGE organoid datasets while keeping its code in the existing
`mge_organoid_pipeline` Git repository.

Included datasets:

1. Varela DIV30
2. Varela DIV90
3. Walsh
4. Bershteyn 2025
5. Bershteyn 2023
6. Siebert 2026

Primary fetal datasets are explicitly out of scope for the initial atlas
integration. They will be introduced later as independent mapping and
benchmarking references.

The workstream does not repeat study-level QC. It starts from the registered
processed objects in [`config/input_objects.tsv`](config/input_objects.tsv),
audits their expression representations and metadata, harmonizes them, and
then freezes a pre-integration input used identically by all integration
methods.

Step 00 also regenerates all-cell cluster UMAP inventories for all six inputs,
including dedicated Varela DIV30/DIV90 panels with cluster numbers and names.
It uses the established DIV30 paper/manual mapping and the DIV90
`cluster_number_name` metadata when available. The familiar DIV90 vertical
plotting orientation is applied without excluding clusters 6/7 or any other
cells from the audit.

## Locations

- Version-controlled code and configuration: this directory
- Great Lakes/Turbo output root:
  `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper2_mge_organoid_atlas`
- Operational handoff: [`HANDOFF.md`](HANDOFF.md)

## Initialize and run the first audit

On Great Lakes:

```bash
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline/paper2_mge_organoid_atlas
./bin/initialize_turbo.sh
./bin/submit.sh 00_input_audit --dry-run
./bin/submit.sh 00_input_audit
```

The submit command prints the run directory. Monitor it with:

```bash
./bin/status.sh /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper2_mge_organoid_atlas/results/00_input_audit/<run_id>
```

No input object is modified by the audit. UMAPs are regenerated as plots from
the saved embeddings; the audit does not recompute the UMAP embedding.

## Working reruns versus frozen runs

The default command creates a new versioned run. During active figure or audit
development, an existing completed working run can instead be explicitly
replaced without creating another results folder:

```bash
./bin/submit.sh 00_input_audit \
  --replace-run 00_input_audit_20260827_130337_6ce39f6 \
  --svg false
```

Replacement is allowed only for a run ID inside the matching Paper 2 step
directory, and is refused while its prior SLURM job is pending or running.
The previous job ID is recorded in the replacement package. Use the default
versioned mode for frozen milestones that must never be overwritten.

PNG and PDF are always generated. SVG is optional (`--svg true|false`) and is
disabled by default. When SVG is enabled, rasterized UMAP point layers are
always exported at 300 dpi while text remains vector/editable. PNG and PDF are
also always exported with 300 dpi rasterized layers.

## Reproducible output contract

Every major result package uses this structure (the `svg/` directory can be
empty when SVG export is disabled):

```text
<output_package>/
├── README.md
├── code/
├── config/
├── figures/
│   ├── png/
│   ├── pdf/
│   └── svg/
├── tables/
├── logs/
├── provenance/
└── SUCCESS.txt or FAILED.txt
```

Each package must contain the exact code and configuration used, registered
inputs and checksums, command and SLURM information, software/session details,
validation tables, and a README. Final figures use the same contract under the
workstream's `final_figures/` directory.
