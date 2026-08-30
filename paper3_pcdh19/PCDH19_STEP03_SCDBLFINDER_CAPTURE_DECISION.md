# PCDH19 Step 03 scDblFinder capture decision

## Decision status

**PROPOSED — awaiting explicit user approval. scDblFinder has not run.**

The proposed `scDblFinder(samples=...)` unit is one capture, identified as
`GEX_1`. The 12 `15662-JZ-*` records must remain cell metadata for reporting,
but must not be passed as 12 independently processed captures.

## Local experimental evidence

The original read-only Cell Ranger root is:

`/nfs/turbo/umms-ziobroj/Ziobro Lab/MGE scRNA-seq/10x_analysis_15662-JZ/Sample_15662-JZ-P01`

Its `config.csv` contains:

- one `[libraries]` row: FASTQ ID `15662-JZ-P01-CNVT01`, feature type
  `Gene Expression`;
- 12 `[samples]` rows, `15662-JZ-1` through `15662-JZ-12`; and
- exactly one Probe Barcode assignment per sample (`B-A03` through `B-D04`).

Its `qc_library_metrics.csv` contains one Gene Expression
`Physical library ID`, `GEX_1`, with 450,788 called cells. The same file then
reports the 12 sample IDs under `Grouped By = Probe barcode ID`, rather than as
12 physical libraries. All per-sample cells sum to the one `GEX_1` total.

Therefore the Cell Ranger word “sample” here denotes a Probe Barcode
demultiplexing assignment inside one multiplexed physical library. It does not
establish 12 independent GEM captures.

## scDblFinder definition

The scDblFinder manual defines `samples` as independently processed units and
says independent captures—not biological samples—should be supplied for
multiplexed data. The current package vignette is even more explicit for Flex:
`samples` indicates different captures, and multiplexed samples processed
together should use the shared batch/capture. It also notes that the large
cell count of one Flex capture can create substantial computational demand;
that engineering issue does not change the biological capture definition.

Sources:

- [scDblFinder package manual](https://bioconductor.org/packages/release/bioc/manuals/scDblFinder/man/scDblFinder.pdf)
- [scDblFinder vignette: samples and 10x Flex](https://plger.github.io/scDblFinder/)
- [Germain et al., scDblFinder method](https://doi.org/10.12688/f1000research.73600.2)

## Consequence for the future Step 03 call

After—and only after—the exact Step 02 checkpoint and this capture definition
are approved, Step 03 should:

1. load the approved Step 02 cells and sparse raw integer gene-count matrix;
2. construct one `SingleCellExperiment` with a constant capture field
   `GEX_1` while preserving `technical_sample_id` separately;
3. run scDblFinder as one capture, either by omitting `samples` or explicitly
   passing the constant `GEX_1` field;
4. use `clusters=TRUE` and `dbr.sd=1`;
5. leave `dbr` unset and leave other model parameters at package defaults
   unless a documented issue requires a reviewed change; and
6. preserve scores, calls, cluster/origin fields, metadata, software versions,
   and RNG provenance without removing any called cell.

The preferred explicit representation is to pass the constant capture field.
It makes the one-capture decision visible in executed code and output
provenance, even though omitting `samples` would also analyze all cells
together.

## Required post-run review; no deletion

The future run must stop with calls retained as metadata. Review must include:

- score distributions overall and by technical sample/design group;
- called count and fraction overall and by technical sample/design group;
- bimodality/separation, including the expected low-score majority and
  distinct high-score population when the classifier separates well;
- stability across frozen reproducibility checks;
- technical-sample composition of calls; and
- call/score location in a diagnostic expression space generated within Step
  03 and clearly distinguished from a final integrated analysis.

No scDblFinder call is an automatic exclusion.

## Important citation qualification

The requested `clusters=TRUE`, `dbr.sd=1` Flex strategy appears in the Howitt
et al. preprint and is also supported directly by current scDblFinder
documentation. However, the Howitt preprint associated with DOI
`10.1101/2024.10.03.616596` has been withdrawn: its authors report that the
original combinatorial model over-predicted multiplets and state that they do
not want the work cited as a project reference. The future workflow can use
the requested settings, but its method justification should rely on the
scDblFinder documentation and Germain et al.; the withdrawn paper should be
recorded as historical context, not affirmative scientific evidence.

Withdrawal notice:
[Howitt et al. preprint status and author statement](https://sciety.org/articles/activity/10.1101/2024.10.03.616596).
