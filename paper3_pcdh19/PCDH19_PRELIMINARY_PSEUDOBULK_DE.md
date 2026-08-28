# PRELIMINARY Step 09: sample-level pseudobulk differential expression

## Status

**Every result in this step is preliminary.** Step 09 is a rapid first-pass
transcriptional analysis, not a finalized differential-expression result. It
does not modify or extend the frozen Step 03–07 classifier.

## Design

Raw integer gene-level Flex UMI-deduplicated probe-pair ligation counts from
Cell Ranger-filtered cells were summed by biological sample and analysis
stratum. These values are assay counts, not literal transcript numbers.

The three prespecified comparisons are:

- `WT_F vs WT_M`: sex-associated difference within WT;
- `HET_F vs WT_F`: HET-associated difference in females; and
- `KO_M vs WT_M`: KO-associated difference in males.

Each was run for:

- all filtered cells;
- transferred proliferating neural progenitors; and
- transferred postmitotic immature neurons.

Each test contains three case and three reference pseudobulk libraries. edgeR
4.4.2 used `filterByExpr`, TMM normalization, a `~0 + design_group` design,
robust dispersion estimation, a quasi-likelihood negative-binomial GLM, and
Benjamini–Hochberg FDR within each contrast/state. Cells were never treated as
replicates. WT-like/KO-like HET classifications were not opened.

## Immediate results

| Contrast | Stratum | Genes tested | FDR ≤ 0.05 | Up in case | Down in case | FDR ≤ 0.05 and abs(log2FC) ≥ 1 |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| WT-F vs WT-M | All cells | 14,835 | 1,839 | 1,058 | 781 | 11 |
| WT-F vs WT-M | Progenitor | 14,153 | 2,452 | 1,490 | 962 | 11 |
| WT-F vs WT-M | Immature neuron | 14,042 | 0 | 0 | 0 | 0 |
| HET-F vs WT-F | All cells | 14,782 | 0 | 0 | 0 | 0 |
| HET-F vs WT-F | Progenitor | 14,100 | 0 | 0 | 0 | 0 |
| HET-F vs WT-F | Immature neuron | 13,803 | 0 | 0 | 0 | 0 |
| KO-M vs WT-M | All cells | 14,892 | 0 | 0 | 0 | 0 |
| KO-M vs WT-M | Progenitor | 14,223 | 0 | 0 | 0 | 0 |
| KO-M vs WT-M | Immature neuron | 14,107 | 0 | 0 | 0 | 0 |

The WT female/male signal is dominated by small effects: among significant
genes, only 47 in all cells and 64 in progenitors have absolute log2 fold
change at least 0.5; only 11 in either analysis reach 1.0. It must not yet be
reported as a finalized biological sex program.

No HET-F/WT-F or KO-M/WT-M gene passes FDR 0.05 in any stratum. This means the
current six-library tests do not resolve a reproducible disease-group signal;
it does **not** prove that no transcriptional effect exists. The MDS plots show
substantial within-group sample variation, including prominent HET JZ-7 and
other sample-specific positions.

Gene-level `Pcdh19` total counts are not different in these tests. This is not
in conflict with the genotype-specific probe result: the gene-level feature
combines A/B/C evidence, while KO loses A but gains downstream B/C evidence.
The engineered-allele question therefore remains probe-specific.

## Marker sanity check

TMM logCPM marker summaries move in the expected broad direction:

- `Mki67` and `Top2a` are higher in progenitor pseudobulks;
- `Gad1`, `Gad2`, `Dcx`, and `Lhx6` are generally higher in immature-neuron
  pseudobulks; and
- these remain broad states, not mature PV/SST interneuron identities.

The sex-marker behavior requires dedicated QC before interpreting the WT-F
versus WT-M result. Several Y-linked features are variably detected in
female-labeled libraries. In a Flex assay this could reflect low-level ambient
signal, probe cross-reactivity, or sample/metadata issues; the current analysis
does not distinguish those possibilities.

## Required caution before finalization

- Confirm biological replicate and preparation/batch metadata.
- Audit X/Y marker probes and sample sex assignments.
- Examine sample influence and leave-one-sample-out stability.
- Decide whether a confirmed batch covariate belongs in the model.
- Treat the reference-transferred states as preliminary annotations.
- Do not interpret zero FDR hits as equivalence or absence of biological
  effect with `n=3` per group.

## Outputs

```text
results/step_09_pcdh19_preliminary_pseudobulk_differential_expression/
```

The package contains the 36 raw integer pseudobulk libraries, metadata, TMM
factors, nine complete edgeR tables, a nine-test summary, 33 preliminary plots,
validation, exact code/dependency hashes, and R session information.

Plot coverage includes three all-sample MDS plots, nine volcano plots, nine MA
plots, nine top-30 sample heatmaps, marker context, library-size QC, and the
nine-test DEG summary.

The output-manifest SHA-256 is
`4d6a09d98f921e937c56c0ab0676c42abb0ea35c29e1dced3959b9110b905287`.

