# Step 08: PCDH19 detectability across transferred developmental states

## Purpose and frozen boundary

Step 08 is the first biological analysis after the frozen Step 03–07 genotype-
classification unit. It asks whether PCDH19 Flex probe detection varies across
the broad GSE94641-transferred MGE developmental states and whether the known
WT-versus-KO probe behavior remains visible within comparable states.

Step 08 does **not** fit, select, calibrate, or evaluate a classifier. It does
not open Step 07 as an input, change either frozen threshold, or produce new
HET classifications. WT-female and HET-female samples are descriptive context;
the only inferential comparison is sex-matched WT-male versus KO-male.

## Inputs and analysis unit

- Step 00 E15.5-primary GSE94641 label-transfer table: all 450,788 cells, all
  raw primary and all-age contextual labels retained.
- Frozen Step 02a PCDH19 audit: exact raw A/B/C UMI-deduplicated probe-pair
  ligation evidence for the same 450,788 barcodes.
- Registered sample key: WT-M JZ-1–3, WT-F JZ-4–6, HET-F JZ-7–9, and KO-M
  JZ-10–12.

The joined per-cell table preserves the entire Step 00 annotation record plus
raw A/B/C counts, binary detection, any-probe detection, and pattern. Probe
UMIs are ligation evidence, not literal transcript counts.

All plotted points and tests use biological sample ID as the replicate. The
WT-M versus KO-M tests enumerate all 20 assignments of three samples to each
group and report exact two-sided label-permutation P values. Benjamini–Hochberg
values are reported separately for the three composition comparisons and 12
probe-detectability comparisons. With only three samples per group, the
smallest attainable two-sided exact P value is 0.10; results are effect-size
and consistency evidence, not a high-powered discovery test.

## Main result

The immediate expectation that immature neurons would provide the highest
PCDH19 detectability was **not** supported. PCDH19 probes were most detectable
in the transferred proliferating-progenitor state in every design group.

Values below are unweighted means of the three biological-sample percentages:

| Group | State | A+ | B+ | C+ | Any PCDH19 probe+ |
| --- | --- | ---: | ---: | ---: | ---: |
| WT-M | Progenitor | 12.26% | 5.89% | 8.61% | 22.12% |
| WT-M | Immature neuron | 6.78% | 3.15% | 4.24% | 12.04% |
| WT-F | Progenitor | 13.17% | 6.36% | 10.06% | 24.27% |
| WT-F | Immature neuron | 7.15% | 3.46% | 4.59% | 12.99% |
| HET-F | Progenitor | 5.21% | 9.95% | 15.19% | 25.89% |
| HET-F | Immature neuron | 2.60% | 5.03% | 7.30% | 13.33% |
| KO-M | Progenitor | 0.08% | 12.79% | 19.85% | 28.59% |
| KO-M | Immature neuron | 0.07% | 6.75% | 10.65% | 15.75% |

Three conclusions survive developmental-state stratification:

1. Probe A remains strongly WT-associated. In both progenitors and immature
   neurons, A detection is appreciable in WT-M/WT-F, nearly absent in KO-M,
   and intermediate in HET-F.
2. B and C remain enriched in KO-M relative to WT-M within both assigned
   states. HET-F again lies between the WT and KO profiles.
3. Progenitors provide roughly twice the A detection and substantially more
   any-probe detection than immature neurons. The progenitor compartment—not
   the immature-neuron compartment—is currently more informative for PCDH19
   probe evidence.

For the sex-matched sample-level contrast, progenitor KO-minus-WT effects were
-12.18 percentage points for A, +6.90 for B, +11.24 for C, and +6.47 for any
probe. Immature-neuron effects were -6.71, +3.60, +6.41, and +3.70 percentage
points. The directional consistency across samples is biologically informative
even though exact P values cannot be below 0.10 with this design.

## State composition

Broad-state composition varies substantially among samples. The unweighted
WT-M versus KO-M mean differences were +5.42 percentage points for progenitors,
-2.09 for immature neurons, and -3.32 for unassigned cells; exact P values were
0.60, 0.90, and 0.20. There is no sample-level evidence here for a reproducible
WT-M/KO-M composition shift.

Pooling cells is misleading: unequal sample sizes can change or reverse the
apparent group difference. Biological samples—not cells—remain the comparison
unit.

## Interpretation limits

- These are broad reference-transferred developmental states, not mature PV,
  SST, or other terminal interneuron identities.
- Raw `mat*` and `pro*` codes remain in the joined table but are not promoted
  to biological cell-type names or used for formal comparisons.
- `not_assigned_neural_state` is lack of a broad-state majority, not a novel
  cell type.
- `000` means no PCDH19 probe event was observed; it does not mean PCDH19 is
  biologically absent.
- HET patterns are descriptive and had no influence on the frozen classifier
  or thresholds.

## Reproducibility

```bash
./paper3_pcdh19/bin/run_step_08_pcdh19_developmental_state_probe_detectability.sh
sbatch paper3_pcdh19/slurm/step_08_pcdh19_developmental_state_probe_detectability.sbatch
```

Outputs are under:

```text
results/step_08_pcdh19_developmental_state_probe_detectability/
```

Primary files are the compressed exact all-cell join, per-sample state
composition, per-sample state/probe detectability, group summaries, the WT-M
versus KO-M exact comparisons, and three concise figures. Detailed reference-
transfer limitations remain in [MGE_REFERENCE_MAPPING.md](MGE_REFERENCE_MAPPING.md),
and the frozen classifier boundary remains in
[PCDH19_GENOTYPE_CLASSIFICATION_HANDOFF.md](PCDH19_GENOTYPE_CLASSIFICATION_HANDOFF.md).

The validated output-manifest SHA-256 is
`296a00caf2e11cc41ad81dddc154f77fb60d0b7b2887320a9a42ee87ba95b76b`.
