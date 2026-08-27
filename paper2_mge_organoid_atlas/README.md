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
audits their expression representations and metadata, creates clean canonical
copies, harmonizes them, and then freezes a pre-integration input used
identically by all integration methods.

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

## Frozen canonical inputs

Step `01_canonical_inputs` is complete. The permanent input layer is:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper2_mge_organoid_atlas/inputs/canonical
```

It contains six minimal Seurat `.rds` files and six equivalent AnnData `.h5ad`
files. They preserve the same QC'd cells, genes, counts, available source RNA
normalized data, and selected provenance metadata as the registered objects.
They contain no old PCA/UMAP reductions, graphs, neighbors, integration/SCT
artifacts, scaled matrices, or command history. No cells or genes were
filtered, and no integration, normalization, gene harmonization, or cell-type
harmonization was performed during this step.

All six RDS/H5AD pairs passed exact expression, identifier, and metadata
equivalence validation. The build used R 4.4.1, Seurat 5.1.0, SeuratObject
5.0.2, Python 3.11.15, Scanpy 1.11.5, and AnnData 0.12.14. Array job `58958446`
and finalizer job `58958447` both completed successfully. Future Paper 2
workflows must read only from this canonical directory.

The canonical layer is protected by an overwrite refusal in the launcher, a
`FROZEN.txt` marker, and SHA-256 manifests. Turbo preserves project-group write
mode bits despite `chmod`, so the workflow safeguards and checksums—not POSIX
mode bits—define the freeze.

To inspect the completed jobs and markers:

```bash
./bin/canonical_status.sh \
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper2_mge_organoid_atlas/inputs/canonical
```

## Step 02 gene-identifier reports: review stop

The report-only feature/gene harmonization audit completed as job `58978281`:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper2_mge_organoid_atlas/results/02_harmonize_genes/02_harmonize_genes_20260827_171616_2483850
```

Five datasets use gene-symbol feature IDs; Siebert 2026 uses a mixture of
22,413 symbol-like IDs and 9,718 Ensembl gene IDs. Features were mapped to
versionless GENCODE 50 Ensembl gene IDs using exact GENCODE matches followed by
uniquely resolvable HGNC approved, previous, and alias symbols. Ambiguous,
unresolved, and duplicate-to-one mappings were reported without guessing.

The six-way exact raw-feature intersection is 14,152. Mapping produces 14,483
common identities present across all six; applying a strict one-source-feature
to-one-common-gene rule yields 14,112. These are reports only—no matrix was
subset, collapsed, concatenated, normalized, HVG-selected, or integrated, and
the canonical inputs were not modified.

The workflow is stopped for review. Before Step 03, choose the strict 14,112
intersection or define a validated duplicate-resolution policy for the 14,483
identity-level set, and decide whether unresolved historical lncRNA symbols
need a separate legacy-reference mapping pass. See [`HANDOFF.md`](HANDOFF.md)
for the complete counts and review questions.

A plot-only extension completed as job `58981144` using only the cached Step 02
tables. It added four 300-dpi PNG/PDF figure pairs: identity-level and strict
UpSet-style intersection diagrams, pairwise overlap/Jaccard heatmaps, and a
mapping-coverage/non-overlap dashboard. SVG remained disabled. The figures
show that 14,483 identities are shared by all six, while 10,569 are unique to
Siebert's mapped feature catalog. The Varela pair has identical mapped feature
sets, as does the Bershteyn pair; Varela–Siebert is the least-overlapping
non-identical pair by mapped Jaccard (0.513). These are feature-catalog
relationships, not expression or biological-similarity results.

The reproducible submission entrypoint for a new versioned report run is:

```bash
./bin/submit_gene_harmonization.sh --dry-run
./bin/submit_gene_harmonization.sh
```

### Step 02b: completed legacy gene-ID recovery extension

The unresolved-ID question was handled as a separately versioned,
report-only extension rather than by rewriting Step 02. Step 02b searches for
original source feature tables and compares exact symbols across GENCODE 27,
32, 35, 44, and 50. Varela DIV30/DIV90 have confirmed original Cell Ranger
feature tables; Walsh, Bershteyn 2025, Bershteyn 2023, and Siebert 2026 do not
currently have an exact original feature table/reference bundle locally.

Consequently, Step 02b distinguishes three evidence classes: confirmed source
feature-table mappings, historical-reference consensus candidates, and
unresolved/ambiguous names. Historical consensus is not treated as proof of a
study's original reference bundle. All proposed mappings and revised overlap
counts remain `REPORT_ONLY_NOT_APPLIED`; canonical inputs and Step 02 are read
only, and the Step 02 review stop remains active.

Job `58986448` completed successfully. It proposed IDs for 22,947 additional
study-feature rows. The proposed identity-level six-way intersection is 14,496
(versus 14,483 in Step 02), and the proposed strict intersection is 14,122
(versus 14,112). These small net gains occur because most recovered historical
features are not shared by all six studies, and confirmed Varela source IDs
also replace some non-original current-symbol assignments. All 34 package
checksums pass; nothing was applied to the canonical objects or Step 02.

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper2_mge_organoid_atlas/results/02b_legacy_gene_id_recovery/02b_legacy_gene_id_recovery_20260827_181028_fd488f5
```

```bash
./bin/submit_legacy_gene_id_recovery.sh --dry-run
./bin/submit_legacy_gene_id_recovery.sh
```

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
