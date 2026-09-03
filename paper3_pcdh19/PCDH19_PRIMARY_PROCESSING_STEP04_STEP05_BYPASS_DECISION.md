# PCDH19 primary-processing Steps 04–05 bypass decision

## Decision

On 2026-09-03, the user explicitly directed the primary-processing workflow to
skip both:

- Step 04, ambient RNA and contamination assessment; and
- Step 05, broad biological contaminant assessment.

Neither step will run or publish a checkpoint. No ambient correction,
corrected-count layer, broad-population exclusion, or other cell/gene removal
is authorized by this bypass.

## Authorized Step 06 lineage

Step 06 technical, sample, and batch diagnostics is authorized to consume the
exact approved Step 02 checkpoint directly:

- run ID: `02_qc_filtering_20260830_124611_97e1bb5`;
- object: `objects/pcdh19_step02_qc_filtered.h5ad`;
- dimensions: 446,349 cells x 19,071 genes;
- SHA-256:
  `fadba4a25a7b6b7320219b21c189b6325687493519ba0fe1bd27efd79606b103`.

Step 06 must not read the rejected Step 03 H5AD or tables. It must not expect
Step 04 or Step 05 outputs, because none exist.

## Step 06 boundary

Diagnostic normalization, feature selection, dimensional reduction,
neighbor-graph construction, and visualization may be used to evaluate
technical, sample, and batch structure. They do not themselves authorize
batch correction, integration, annotation, or biological removal. Any
corrective transformation must first be justified by the diagnostics and
explicitly approved.
