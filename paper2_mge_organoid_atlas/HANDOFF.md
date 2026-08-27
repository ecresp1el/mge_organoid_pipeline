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
5. Do not overwrite a frozen input or selected atlas. Active working packages
   may use the explicit replacement mode described below; frozen milestones
   require a new versioned package documenting what they supersede.

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

### Step 00 UMAP visualization contract

Step 00 does not recompute UMAP embeddings. It reads each processed object's
saved UMAP coordinates and regenerates publication-quality cluster inventory
plots in mandatory PNG/PDF and optional editable-text SVG. Rasterized layers
in every enabled format are fixed at 300 dpi.

- Every registered cell with finite saved UMAP coordinates is shown.
- No study-level QC is rerun.
- DIV90 current clusters 6/7 are retained in the input-audit plots even though
  they were excluded from some prior publication-style final figures.
- DIV90 receives the established plotting-only transform
  `UMAP1_plot = UMAP1_original`, `UMAP2_plot = -UMAP2_original`.
- Original coordinates are retained unchanged in the per-cell audit tables.
- DIV30 uses the documented seven-raw-cluster to five-paper/manual-class
  mapping used by the existing cross-study final figures. Both raw and
  displayed cluster counts are reported.
- DIV90 uses `cluster_number_name` from the object when present, while also
  retaining `seurat_clusters` as the raw ID.
- Other studies use the first available documented cluster-name metadata
  field; if none exists, the plot and tables state that only raw IDs were
  available.

Required UMAP outputs include a six-study grid, a dedicated DIV30/DIV90 grid,
one figure per study, cluster-count and label-source tables, a figure manifest,
and the per-cell original-coordinate tables.

## Independent Paper 2 step numbering

Step numbers are local to `paper2_mge_organoid_atlas/`; they do not continue
the historical top-level pipeline numbering.

| Step | Status | Purpose |
| --- | --- | --- |
| `00_input_audit` | Completed, replaced in place, and visually checked; current job `58956196` | Registered and checksummed the six processed objects; inventoried assays, layers, features, metadata, reductions, and saved UMAP/cluster labels; regenerated all-cell UMAP inventory figures. Job `58956196` replaced the original job `58955368` package in the same run directory. |
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
figures/svg/              optional editable SVG exports when enabled
tables/                   machine-readable results and validation summaries
logs/                     SLURM stdout/stderr and pipeline logs
provenance/               command, job ID, Git state, environment, checksums
SUCCESS.txt or FAILED.txt unambiguous terminal state
```

For expensive analyses, figure formatting must use cached, validated tables or
objects through a plot-only step. A cosmetic change must not silently rerun or
replace the integration.

There are two explicit output modes:

- Default versioned mode creates a new timestamped run directory and is used
  for frozen milestones.
- Working replacement mode uses `--replace-run RUN_ID`, clears and regenerates
  that exact completed run directory, and does not create another results
  folder. It is refused for targets outside the step directory and while the
  prior job is pending or running. The replacement records the prior job ID.

PNG and PDF are mandatory for UMAP renders. SVG is controlled by
`--svg true|false` and defaults to false to reduce asset accumulation. When
enabled, SVG rasterized point layers are fixed at 300 dpi; SVG text remains
editable. PNG and PDF rasterized layers are also fixed at 300 dpi.

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

Replace the existing completed working audit instead of making a new folder:

```bash
./bin/submit.sh 00_input_audit \
  --replace-run 00_input_audit_20260827_130337_6ce39f6 \
  --svg false
```

Enable SVG when it is specifically needed:

```bash
./bin/submit.sh 00_input_audit \
  --replace-run 00_input_audit_20260827_130337_6ce39f6 \
  --svg true
```

Versioned submission creates a new timestamped output package; replacement
submission reuses only the explicitly named completed working directory. Both
snapshot code and configuration, copy the exact submitted `.sbatch` into the
Paper 2 `jobs/` directory, record the current Git commit and working-tree
status, and capture the parsable SLURM job ID. The printed run directory is the
durable handle for monitoring:

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
- all-cell UMAP source/cluster-label inventory and raw/display cluster counts;
- six-study, Varela-paired, and per-study cluster UMAP figures in mandatory
  PNG/PDF and optional SVG;
- per-cell saved UMAP coordinates with original IDs, collision-safe candidate
  atlas IDs, sample provenance, raw cluster IDs, displayed names, and mapping
  sources;
- per-study default-assay feature lists;
- pairwise and six-way feature-overlap summaries;
- R/Seurat session information;
- exact input registry, R script, SLURM script, submit command, job ID, Git
  commit/status, package checksum manifest, and logs.

The audit does not decide which expression layer should be used. That decision
belongs to `01_harmonize_genes` and must cite the audit evidence.

## Step 00 completed evidence and findings

The current completed audit package is:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper2_mge_organoid_atlas/results/00_input_audit/00_input_audit_20260827_130337_6ce39f6
```

The replacement run loaded all six registered objects and found 404,789 cells
with saved UMAP coordinates:

| Study | Cells | Default assay | Raw cluster IDs | Display groups in audit |
| --- | ---: | --- | ---: | ---: |
| Varela DIV30 | 90,631 | RNA | 7 | 5 paper/manual classes |
| Varela DIV90 | 22,338 | RNA | 13 | 13 `cluster_number_name` classes |
| Walsh | 4,519 | RNA | 24 | 156 observed raw-cluster/`walsh_group` combinations |
| Bershteyn 2025 | 124,583 | RNA | 9 | 9 raw IDs; no candidate name column was present |
| Bershteyn 2023 | 98,042 | RNA | 6 | 9 observed raw-cluster/`celltype` combinations |
| Siebert 2026 | 64,676 | SCT | 28 | 28 raw IDs; no candidate name column was present |

All registered cells with finite saved UMAP coordinates were retained in the
audit plots. In particular, the DIV90 audit includes current clusters 6 and 7;
their exclusion from some older published-style figures is not an input-audit
filter.

The default-assay feature union contains 54,746 genes, of which 13,965 occur
in all six objects. The two Varela objects have sample-split Seurat v5 RNA
count/data layers. Siebert defaults to SCT while the other five objects default
to RNA. These are completed audit observations, not harmonization decisions;
the expression layer, feature-identity policy, and method-compatible matrix
must be decided in Step 01.

### Varela DIV30 raw-cluster gap

The DIV30 object contains 90,631 cells and exactly seven observed raw Seurat
cluster IDs: `0, 1, 2, 3, 4, 6, 7`. Both `seurat_clusters` and
`RNA_snn_res.0.2` contain those same values and agree for all cells. There are
no missing cluster assignments. Raw cluster ID `5` is absent from the source
object and from the derived AnnData categorical levels; it was not removed by
the Paper 2 plotting code.

Do not confuse the absent raw cluster ID `5` with publication/display class 5.
Publication class `5 - MGE subpallial neurons` is present and contains 17,287
cells from raw Seurat cluster `2`. The full DIV30 mapping is:

| Raw Seurat cluster | Cells | Paper/manual display class |
| ---: | ---: | --- |
| 0 | 19,148 | 1 - Radial glia |
| 1 | 18,907 | 3 - SST+ cIN |
| 2 | 17,287 | 5 - MGE subpallial neurons |
| 3 | 15,932 | 1 - Radial glia |
| 4 | 13,750 | 4 - PV neuron precursor |
| 5 | 0 / absent | not present in the registered source object |
| 6 | 3,696 | 2 - Inhibitory progenitors |
| 7 | 1,911 | 1 - Radial glia |

The nonconsecutive raw IDs are consistent with an upstream removal, merge, or
manual relabeling before the legacy `Day30.rds` was saved, but the current
repository and Turbo legacy directory do not contain a pre-change object or
source script that establishes which event occurred or why. This unresolved
history must remain an explicit provenance limitation. Do not manufacture a
cluster 5 or infer its biology without locating an earlier object or original
analysis record.

## Resume point

Current state as of 2026-08-27:

- Paper 2 scope is defined.
- Six candidate processed objects are registered and exist on Turbo.
- The code and Turbo directory scaffold are established.
- Step `00_input_audit` was successfully replaced in place from the six source
  objects as SLURM job `58956196` (`COMPLETED`, exit `0:0`, elapsed
  `00:05:44`, peak batch RSS approximately 60.3 GiB).
- Its current completed run directory is
  `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper2_mge_organoid_atlas/results/00_input_audit/00_input_audit_20260827_130337_6ce39f6`.
- The replacement configuration records `OUTPUT_MODE=replace`,
  `REPLACED_PREVIOUS_JOB_ID=58955368`, `PNG_DPI=300`, `PDF_DPI=300`,
  `SVG_DPI=300`, and `MAKE_SVG=false`.
- The run has a `SUCCESS.txt` marker, complete audit tables, 8 PNG files,
  8 PDF files, no SVG files, and 49 verified package checksum entries.
- All eight PNG files report embedded resolution of approximately
  299.9994 dpi, and every figure-manifest row records PNG/PDF at 300 dpi.
- Exactly one Step 00 results directory exists; the replacement did not create
  another timestamped output directory.
- The Varela paired UMAP was visually checked: DIV30 contains 90,631 cells,
  seven raw clusters mapped to five paper/manual classes; DIV90 contains
  22,338 cells and all 13 current `cluster_number_name` clusters, including
  clusters 6/7, with the documented plotting-only vertical flip.
- No master dataset has been created or frozen.
- No integration method has been run for Paper 2.

The next safe action is to review the completed assay/layer, feature-overlap,
metadata, and cluster-label evidence and use it to specify Step
`01_harmonize_genes`. Do not implement integration-specific preprocessing
until the expression-layer and gene-identity decisions are documented from
this audit.
