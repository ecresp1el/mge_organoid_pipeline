# Paper 2 MGE organoid atlas: operational handoff

Last updated: 2026-08-27

## Scientific objective

Build a cross-study in-vitro MGE organoid single-cell atlas from exactly six
processed datasets:

- Varela DIV30
- Varela DIV90
- Walsh
- Bershteyn 2025
- Bershteyn 2023
- Siebert 2026

The initial workflow is:

```text
six processed organoid objects
  -> audit expression layers and metadata
  -> harmonize genes and metadata
  -> freeze one pre-integration master input
  -> compare integration methods on that same input
  -> evaluate batch removal and biological preservation
  -> select a supported integrated representation
  -> freeze the in-vitro MGE organoid atlas
```

Primary fetal datasets are not members of the initial integration. They will
be used later as independent references for mapping and biological evaluation.

## Non-negotiable scope decisions

1. Do not repeat the six studies' original QC by default.
2. Do not silently remove cells, genes, assays, layers, or metadata during the
   input audit.
3. Do not include fetal cells in the frozen pre-integration organoid object.
4. Do not let individual integration methods receive different starting cells,
   genes, or harmonized covariates unless a method-specific requirement is
   explicitly documented and evaluated.
5. Do not overwrite a frozen input or selected atlas. Create a new versioned
   package and document why it supersedes the prior version.

## Canonical locations

```text
Code repository:
/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline

Paper 2 code/configuration:
/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline/paper2_mge_organoid_atlas

Existing MGE project/data root:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder

Paper 2 Turbo root:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper2_mge_organoid_atlas
```

The Paper 2 Turbo root owns its own `results/`, `logs/`, `jobs/`, and
`final_figures/`. Existing processed objects remain read-only upstream inputs;
they are not copied or altered merely to establish the workstream.

## Registered input objects

The authoritative candidate-input registry is
[`config/input_objects.tsv`](config/input_objects.tsv). It currently selects
the annotated Walsh final object because that is the object used by the
existing cross-study workflows; the registry documents its relationship to
the post-QC unannotated object.

At this stage, "registered" does not mean "frozen." The six objects become a
frozen pre-integration input only after the following are reviewed:

- file identity and SHA-256 checksum;
- object class and dimensions;
- assays and available counts/data/scale-data layers;
- feature namespaces and cross-study gene overlap;
- metadata columns, types, missingness, and candidate harmonization fields;
- reductions already present, recorded for provenance but not assumed to be
  valid shared integration inputs;
- per-study cell identifiers and a collision-safe atlas cell ID plan.

## Independent Paper 2 step numbering

Step numbers are local to `paper2_mge_organoid_atlas/`; they do not continue
the historical top-level pipeline numbering.

| Step | Status | Purpose |
| --- | --- | --- |
| `00_input_audit` | Implemented, not yet submitted | Register and checksum the six processed objects; inventory assays, layers, features, metadata, and reductions. |
| `01_harmonize_genes` | Planned | Resolve gene identifiers, duplicates, shared/union feature policies, and method-compatible matrices without integration. |
| `02_harmonize_metadata` | Planned | Create a documented common schema for study, dataset, sample, replicate, age/time point, cell labels, and QC provenance. |
| `03_freeze_preintegration` | Planned | Produce and validate the immutable six-study pre-integration master object/package. |
| `10_scvi` | Planned | Run scVI from the frozen input. |
| `11_liger` | Planned | Run LIGER from the frozen input. |
| `12_scpoli` | Planned | Run scPoli from the frozen input. |
| `13_seurat_rpca` | Planned | Run Seurat RPCA from the frozen input. |
| `20_integration_benchmark` | Planned | Compare batch removal and biological preservation using common metrics and panels. |
| `30_freeze_atlas` | Planned | Select, justify, validate, and package the supported organoid atlas representation. |
| `40_fetal_reference_mapping` | Deferred | Map or benchmark the frozen organoid atlas against independent fetal GE references. |

New steps should receive a new number and descriptive name. Existing step
numbers should not be repurposed after they have generated a reviewed output.

## Major-output reproducibility contract

Every major `results/` package and every `final_figures/` package must contain:

```text
README.md                 scientific purpose, inputs, methods, outputs, status
code/                     exact submitted scripts/wrappers used for this run
config/                   submitted and resolved configuration/input registry
figures/png/              high-resolution raster exports when applicable
figures/pdf/              publication/vector exports when applicable
figures/svg/              editable SVG exports when applicable
tables/                   machine-readable results and validation summaries
logs/                     SLURM stdout/stderr and pipeline logs
provenance/               command, job ID, Git state, environment, checksums
SUCCESS.txt or FAILED.txt unambiguous terminal state
```

For expensive analyses, figure formatting must use cached, validated tables or
objects through a plot-only step. A cosmetic change must not silently rerun or
replace the integration.

## SLURM submission and monitoring

Initialize the Turbo directories idempotently:

```bash
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline/paper2_mge_organoid_atlas
./bin/initialize_turbo.sh
```

Inspect the first submission without creating or submitting a run:

```bash
./bin/submit.sh 00_input_audit --dry-run
```

Submit:

```bash
./bin/submit.sh 00_input_audit
```

Submission creates a new timestamped output package, snapshots code and
configuration, copies the exact submitted `.sbatch` into the Paper 2 `jobs/`
directory, records the current Git commit and working-tree status, and captures
the parsable SLURM job ID. The printed run directory is the durable handle for
monitoring:

```bash
./bin/status.sh <run_directory>
```

`status.sh` reports `sacct`, terminal markers, and recent stdout/stderr. A job
is complete only when SLURM reports successful completion, `SUCCESS.txt`
exists, required outputs exist and are nonempty, and the output README and
provenance are sufficient to reconstruct the run.

## Step 00 expected outputs

`00_input_audit` is read-only with respect to all six source objects. Its
versioned package will include:

- audited input registry with file size, modification time, and SHA-256;
- object and default-assay summaries;
- assay/layer inventory;
- metadata-column type, missingness, cardinality, and example-value inventory;
- reduction inventory;
- per-study default-assay feature lists;
- pairwise and six-way feature-overlap summaries;
- R/Seurat session information;
- exact input registry, R script, SLURM script, submit command, job ID, Git
  commit/status, package checksum manifest, and logs.

The audit does not decide which expression layer should be used. That decision
belongs to `01_harmonize_genes` and must cite the audit evidence.

## Resume point

Current state as of 2026-08-27:

- Paper 2 scope is defined.
- Six candidate processed objects are registered and exist on Turbo.
- The code and Turbo directory scaffold are established.
- Step `00_input_audit` is implemented but has not been submitted.
- No master dataset has been created or frozen.
- No integration method has been run for Paper 2.

The next safe action is to dry-run, submit, monitor, and review
`00_input_audit`. Do not implement integration-specific preprocessing until
that audit has established the available expression and metadata contract.
