# DIV30 antecedent biological identity extension

Authoritative request: `provenance/USER_REQUEST.txt` in the output directory.
This extends the completed DIV30/DIV90 purification work without another large
gate search, changed target labels, or replacement of the original PI package.

## Outputs

All new results are under:

`/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_antecedent_identity_v1/`

- `index.html`: biological report, linked tables and complete figure gallery.
- `REPORT.md`: **WHAT IS THE DIV30 ENRICHABLE ANTECEDENT STATE?**, answering all
  13 final questions and 17 requested analysis areas.
- `DIV30_Biological_Identity_Review.zip`: portable reports, summary tables,
  PDF/SVG/600-dpi PNG figures, source code and provenance.
- `DIV30_Biological_Identity_Figures.pdf`: within-study figure book. Mapping
  and external figures have separate PDF/SVG/PNG exports.
- `within/`, `mapping/`, `external/`: component results and findings.
- `external_extra/samarasinghe/`: validated official Seurat extraction and
  source metadata; controls separated from Rett cells.
- `provenance/`: frozen-output hashes, independent validation, original request
  and code. Large expression arrays and per-cell compressed tables remain on
  Turbo rather than in the portable review archive.

## Biological conclusion

The supported working identity is a **GRIA2-high / OPCML-high / NOTCH1-low,
SLC6A1-rich, predominantly postmitotic GABAergic developmental state with
PV-associated transcription and substantial SST expression**. It is not an
established prospective PV-producing population. A broader state could include
a PV-producing component; the data do not track fate.

- Frozen target composition: 20.0% to 58.36%, recovery 50.81%, total yield
  17.41%. The old held-sample gate-reselection result remains 53.49%
  composition and 39.73% recovery.
- A/B/C/D groups contain 9,210/8,917/6,571/65,933 cells. Captured targets and
  captured non-targets are similarly postmitotic (90.63% and 90.47%).
- SLC6A1 detection rises from 48.16% of all DIV30 cells to 87.21% retained;
  SST is detected in 52.54% retained. The MGE composite is not enriched.
- Flexible generic developmental covariates predict gate membership at
  held-sample AUC 0.828; adding the direct PV module gives 0.858. The larger
  gain for predicting the original target is partly definitional because the
  target already includes direct/ranked PV scores. These gains are predictive
  associations, not causal variance explained or independent fate validation.
- OPCML gives the largest conditional target-composition increment. The
  NOTCH1 clause adds 2.92 percentage points while losing 643 targets.
- Retained non-targets are mostly related inhibitory/maturing neurons;
  41.26% lie in the score's 70th–80th percentiles, but 14.73% lie below median.
- 2E-MW supplies 64.96% of retained targets. Omitting that sample while keeping
  the fixed rule still enriches 12.53% to 44.91%, with 37.81% recovery. This
  descriptive omission is not another held-sample gate-selection analysis.
- Forward strict-PV neighbor fraction decreases (11.55% retained vs 16.95%
  removed), while dual-high increases (52.95% vs 34.92%). The strict-PV
  decrease survives within-sample, postmitotic, distance-support,
  sample-balanced and omit-SLC6A1 sensitivities.
- Reverse strict-PV capture is weak/nonspecific, about 1.02-fold after source
  sample balancing. Essentially all reverse strict-PV queries exceed the
  internal DIV30 held-sample 95th-percentile distance threshold. Hubness is
  substantial; nearest neighbors do not establish cross-time equivalents.
- External mappings do not establish equivalent states, unusual MGEO state
  abundance, or a unique pathway explaining later PV protein. Failed mapping
  support is not proof of external biological absence. Matched-state program
  contrasts are not estimable under the declared criteria; no whole-study
  comparison is substituted.
- Bounded antibody feasibility is documented after the biological analysis.
  Live human neural sorting remains unvalidated; no purchase is recommended.

## Frozen definitions and important caveats

The exact RNA rule, using log1p(counts per 10,000), is:

```text
GRIA2 > 1.1180085897445677
AND OPCML > 0.47731003165245056
AND NOTCH1 <= 0.88839703798294067
```

Compare stored float32 expression after promotion to float64 with the frozen
float64 thresholds. Implicit float32 scalar comparison changes boundary cells.
The original top20 consensus target has 18,127 of 90,631 DIV30 cells from six
source samples. The 9583/9853 alias remains reconciled; metadata samples 7–9
are absent. Do not redefine the target, postmitotic flag or old module scores.

The PV module is MEF2C/KCNC1/KCNC2/TAC1. SLC6A1 is absent from this explicit
module but was already included in the frozen 1,000 mapping features and sparse
model (coefficient +0.010314168). Its association with the old consensus is not
fully independent. The 999-feature mapping is a labeled sensitivity, with
original labels unchanged. TUBB3 is unassayed, although the old context panel
contains a preserved zero placeholder; requested identity genes are assayed.

DIV90's best observed >=50% recovery solution remains 31.3% strict-PV
composition within the completed search. FGFR2/PTPRS retains about 89% but
changes SST-biased composition only 22.5% to 21.3%; independent statistical
significance is not established. This extension does not promote it again.

Cell-level rank-test/FDR values describe sampled cells and are not independent
biological-replicate significance. Only six DIV30 samples are available, one
per line/condition. External cell selection, age, assay and sample provenance
limit comparisons. Walsh already reports later PVALB protein and fast-spiking
cells; it is not a PV-protein-negative protocol control. Xiang remains excluded
because its sample biology/metadata validation is unresolved, not because RDS
files are unusable. Samarasinghe was successfully extracted through the
established Slurm route after local memory failure; do not repeat extraction.

## Reproduction and validation

Python: `/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/python`.
Use `PYTHONDONTWRITEBYTECODE=1`, `OPENBLAS_NUM_THREADS=2`, `OMP_NUM_THREADS=2`
and `MKL_NUM_THREADS=2`. Code is in `python_notebooks/scripts/` and copied into
the review archive. The independent analysis components are
`div30_identity_within.py`, `div30_identity_mapping.py`, and
`div30_identity_external.py`; their outputs already exist and should be reused.

To rebuild only presentation after reviewing the completed outputs:

```bash
/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/python python_notebooks/scripts/build_div30_identity_report.py
/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/python python_notebooks/scripts/package_div30_identity.py
```

`validate_div30_identity.py` independently checks original column preservation,
all four groups, fixed clause subsets and all 362,524 archived forward
neighbor-state fractions. Mapping reconstruction checks 4,096 queries against
the frozen neighbor order/distances. External normalization checks original
raw counts, existing log-normalized layers and reconstruction of the MGEO
reference. Samarasinghe normalization/cell-order/QC checks are separately
recorded in its export provenance.

`div30_identity_preservation.py verify` confirms the five old result directories
against the original SHA256 snapshot: 2,245 files, 7,473,453,634 bytes, unchanged.
Do not regenerate that snapshot over previous results. Package validation
checks required artifacts, every local HTML link and its presence in the ZIP,
and ZIP CRC; its result is `REVIEW_PACKAGE_VALIDATION.json`.
