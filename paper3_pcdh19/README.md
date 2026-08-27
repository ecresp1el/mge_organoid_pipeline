# Paper 3: Ziobro PCDH19 MGE single-cell RNA-seq

This workstream is reserved for the Ziobro PCDH19 MGE single-cell paper. It is
separate from the Paper 2 cross-study atlas and uses the Ziobro Turbo
allocation, independent step numbering, and independent outputs.

## Correct Turbo allocation

The current account is a member of `umms-ziobroj-turbo` and can read and
traverse:

```text
/nfs/turbo/umms-ziobroj
```

The source dataset discovered there is:

```text
/nfs/turbo/umms-ziobroj/Ziobro Lab/MGE scRNA-seq
```

This is service request `15662-JZ`: a mouse 10x Flex v2 Gene Expression
experiment with 12 probe-barcoded samples. It was sequenced on the Ultima
Genomics UG 100 platform and processed against the GRCm39 mouse probe set.

The earlier Miranda bulk RNA-seq directories are not Julia Ziobro data and are
not Paper 3 inputs.

## Observed data layout

The source directory contains:

- two raw paired FASTQs for the pooled gene-expression library;
- a pooled Cell Ranger multi output with raw and filtered matrices, QC reports,
  Loupe files, and alignment products;
- 12 demultiplexed sample-level filtered feature-barcode matrices;
- sample-level Cell Ranger web summaries, metrics, and Loupe files;
- a second copy of the 12 sample matrices inside the AGC Dropbox delivery;
- MD5 and demultiplexing records.

The Cell Ranger QC summary reports 450,788 called cells across the 12 samples.
Each delivered sample matrix contains 19,071 features. Sample-level called-cell
counts range from 11,085 to 60,680.

## Locations

- Version-controlled Paper 3 code/configuration: this directory
- Source registry: [`config/input_candidates.tsv`](config/input_candidates.tsv)
- Technical sample inventory:
  [`config/sample_manifest_draft.tsv`](config/sample_manifest_draft.tsv)
- Operational handoff: [`HANDOFF.md`](HANDOFF.md)
- Complete pipeline I/O and interpretation contract:
  [`PIPELINE_IO_AND_BIOLOGICAL_SCOPE.md`](PIPELINE_IO_AND_BIOLOGICAL_SCOPE.md)
- Paper 3 output root:
  `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19`

On a workstation where the same Turbo allocation is mounted under `/Volumes`,
the corresponding output path is:

```text
/Volumes/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19
```

The source directory remains read-only. Paper 3 analysis products use the same
`inputs/`, `results/`, `logs/`, `jobs/`, and `final_figures/` layout as Paper 2.

## Safe access check

```bash
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline/paper3_pcdh19
./bin/check_candidate_access.sh
```

This checks path metadata only. It does not open, checksum, copy, or modify the
large scientific files.

## Loupe files

The delivery includes a pooled 5.3 GB `raw_cloupe.cloupe` and one
`sample_cloupe.cloupe` for each of the 12 technical samples. The per-sample
files are also duplicated in the AGC Dropbox delivery directory.

Inspection of the embedded Loupe v9 metadata found a blank run description,
no custom cell-data tables, and no PCDH19 genotype or condition labels. The
files contain the expression matrices, vendor UMAP/t-SNE, graph and k-means
clusters, and technical Cell Ranger metadata, but not the missing biological
sample key.

## Current stop point

The technical dataset is identified, but the biological sample key is not yet
present in the files inspected. Before QC or clustering, obtain the mapping
from `15662-JZ-1` through `15662-JZ-12` to genotype, line, condition, time
point, differentiation batch, replicate, and any PCDH19 mosaic design fields.

## Locked Pcdh19 probe audit

The technical Pcdh19 Flex audit is implemented as a biological-metadata-free,
version-locked pipeline step. It checks the delivered panel, exact 10x v2.0.0
BED and metadata, and GRCm39-2024-A GTF by SHA-256 before extracting the three
validated Pcdh19 probes from every technical sample.

Run all 12 samples locally or in an interactive Great Lakes allocation:

```bash
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
./paper3_pcdh19/bin/run_pcdh19_probe_audit_all.sh
```

Submit the same entry point through SLURM:

```bash
sbatch paper3_pcdh19/slurm/pcdh19_probe_audit_all.sbatch
```

The runner always validates `15662-JZ-1` against the frozen prototype before
processing samples 2 through 12. Existing validated outputs are checksum-
verified and retained; incomplete or different outputs cause a hard failure.
Final products are written to:

```text
results/pcdh19_probe_audit/
  references/
    pcdh19_probe_coordinates.tsv
    reference_manifest.tsv
  per_sample/<technical_sample_id>/
    pcdh19_probe_patterns.tsv
    pcdh19_probe_summary.tsv
    validation.tsv
    checksums.sha256
  combined/
    all_samples_pcdh19_probe_summary.tsv
    all_samples_pcdh19_pattern_summary.tsv
    all_samples_validation.tsv
  software_environment.tsv
  output_manifest.tsv
```

The lock is
[`config/pcdh19_probe_audit.lock.json`](config/pcdh19_probe_audit.lock.json),
and the implementation is
[`scripts/pcdh19_probe_audit.py`](scripts/pcdh19_probe_audit.py). This step
contains technical sample IDs only and must not be joined to genotype, sex,
condition, or other biological annotations until the sample key is recovered.

For the exact role of every shell, SLURM, configuration, Python input, output
column, validation, and biological boundary, read
[`PIPELINE_IO_AND_BIOLOGICAL_SCOPE.md`](PIPELINE_IO_AND_BIOLOGICAL_SCOPE.md).
The validated production files record the SHA-256 of the frozen Python
implementation; documentation changes must not be represented as having
generated those existing results.
