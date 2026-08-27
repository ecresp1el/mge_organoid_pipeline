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

The submit command prints the immutable run directory. Monitor it with:

```bash
./bin/status.sh /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper2_mge_organoid_atlas/results/00_input_audit/<run_id>
```

No input object is modified by the audit.

## Reproducible output contract

Every major result package uses this structure:

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

