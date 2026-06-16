# Handoff: Jia S9 vs Shi Lineage Overlap

Date: 2026-06-16

## Purpose

Create a standalone workbook-only workflow to compare Jia et al. Science Data
S9 modules against Shi et al. marker tables, with the explicit biological
question:

```text
Are Jia's five MGE lineage modules mostly a higher-resolution split of Shi's
MGE fate map, and where does CRABP1/ANGPT2 remain contested?
```

This workflow is not label transfer and does not use expression matrices. It
compares marker/module gene lists, reports overlap statistics, and preserves a
curated biological crosswalk for interpretation.

## Biological Crosswalk To Preserve

| Collapsed lineage | Jia exact names | Jia module | Jia logic | Shi matching names |
|---|---|---:|---|---|
| Early subpallial cholinergic | LHX8/ISL1 | Module 1 | VZ RGC-derived, early, adult subpallial cholinergic; ZIC4/HMGA1-associated. | M5/M6, LHX8/ISL1/GBX2/ZIC1 cholinergic subpallial neurons |
| Early subpallial GABAergic | NR2F1/NR2F2 | Module 2 | VZ RGC-derived, early, subpallial GABAergic. | M4/M7, NR2F1/NR2F2/ZFHX3 subpallial GABAergic neurons |
| Cortical MGE interneuron output | EPHA5/MEF2C | Module 3 | SVZ RGC-derived, cortex-bound; maps broadly to adult cortical GABAergic types. | mostly M2, cortical interneuron branch |
| Cortical MGE interneuron output, chandelier-biased | LHX6/NFIA | Module 4 | SVZ RGC-derived, cortex-bound; Jia maps it especially to adult chandelier cells. | M2, possibly part of Shi's cortical branch |
| CRABP1 bridge / contested lineage | CRABP1/ANGPT2 | Module 5 | Jia calls it SVZ RGC-derived subpallial GABAergic, EPHA5-low/non-DEN; ETV1-associated. | M3, CRABP1/ANGPT2/ETV1/NFIA; Shi supports striatal plus cortical/PV-associated interpretation |

Collapsed model:

```text
Jia 5 -> Shi 4-ish

1. LHX8/ISL1 -> Shi M5/M6 -> subpallial cholinergic
2. NR2F1/NR2F2 -> Shi M4/M7 -> subpallial GABAergic
3. EPHA5/MEF2C + LHX6/NFIA -> Shi M2 -> cortical MGE interneurons
4. CRABP1/ANGPT2 -> Shi M3 -> CRABP1+ striatal/cortical-associated MGE lineage,
   not safely subpallial-only
```

## Inputs

Jia workbook:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/reference/science.adw1803_data_s9.xlsx
```

Default Shi workbook:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/reference/shi_2021_tables_s2_to_s9/science.abj6641_table_s6.xlsx
```

User-provided original Mac path for the same Shi Table S6 file:

```text
/Users/ecrespo/Downloads/science.abj6641_tables_s2_to_s9/science.abj6641_table_s6.xlsx
```

The file already exists on Great Lakes/NFS at the default reference path above,
so the first run does not require rsync.

Important caveat:

```text
Shi Table S6 contains GE progenitor subcluster marker sets named pC1-pC3,
pL1-pL3, and pM1-pM4. The available Shi S2-S9 Excel workbooks are present and
readable, but the narrative M2/M3/M5/M6/M7 names from the biological crosswalk
are not encoded as literal `cluster` values in the parsed Excel marker tables.
The workflow therefore computes overlap against the actual input clusters and
separately marks whether the expected crosswalk labels are directly encoded.
```

If a later source explicitly maps Shi numeric subclusters or table rows onto
the narrative M2/M3/M4/M5/M6/M7 lineage names, preserve that as a separate
crosswalk input instead of pretending those labels were already table clusters.

## Code Added

Standalone Python workflow:

```text
python_notebooks/scripts/jia_s9_shi_lineage_overlap.py
```

Slurm template:

```text
slurm_templates/45_jia_s9_shi_lineage_overlap.sbatch.template
```

Default run label:

```text
jia_s9_shi_table_s6_overlap_v1
```

Default output directory:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/jia_s9_shi_lineage_overlap/jia_s9_shi_table_s6_overlap_v1
```

Exhaustive S3-S9 run label:

```text
jia_s9_shi_s3_to_s9_exhaustive_v1
```

Exhaustive S3-S9 output directory:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/jia_s9_shi_lineage_overlap/jia_s9_shi_s3_to_s9_exhaustive_v1
```

Final executed Slurm run:

```text
job_id: 51811782
state: COMPLETED
exit_code: 0:0
stderr: empty
```

Final report:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/jia_s9_shi_lineage_overlap/jia_s9_shi_s3_to_s9_exhaustive_v1/reports/jia_s9_shi_lineage_overlap_report.md
```

Workbook inventory/accounting report:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/jia_s9_shi_lineage_overlap/jia_shi_workbook_inventory_v1/reports/jia_shi_workbook_inventory_report.md
```

Focused pM1-pM4 run label:

```text
jia_s9_shi_s6_pm_only_v1
```

Focused pM1-pM4 output directory:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/jia_s9_shi_lineage_overlap/jia_s9_shi_s6_pm_only_v1
```

Focused pM1-pM4 report:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/jia_s9_shi_lineage_overlap/jia_s9_shi_s6_pm_only_v1/reports/jia_s9_shi_lineage_overlap_report.md
```

## Focused Shi pM1-pM4 Result Snapshot

The biologically relevant direct Shi labels for this question are Table S6
`pM1`, `pM2`, `pM3`, and `pM4`. These labels are present in the Excel workbook
and were analyzed directly with:

```text
--shi-xlsx science.abj6641_table_s6.xlsx
--include-shi-cluster-regex '^pM[1-4]$'
```

Headline counts:

```text
Shi pM marker genes represented: 137 unique genes across 4 marker sets
Shared Jia/pM vocabulary: 48 genes
Jia-only vocabulary in this pM comparison: 1108 genes
Shi-pM-only vocabulary in this pM comparison: 89 genes
```

Table S6 pM gene-set sizes:

| Shi cluster | n genes |
|---|---:|
| pM1 | 48 |
| pM2 | 38 |
| pM3 | 48 |
| pM4 | 46 |

Best descriptive pM hit per Jia full module:

| Jia module | Jia exact name | Best pM cluster | Overlap genes | Fraction of Jia module | Fraction of pM set |
|---|---|---|---:|---:|---:|
| module_1 | LHX8/ISL1 | pM2 | 6 | 2.0% | 15.8% |
| module_2 | NR2F1/NR2F2 | pM3 | 7 | 5.9% | 14.6% |
| module_3 | EPHA5/MEF2C | pM1 | 11 | 4.8% | 22.9% |
| module_4 | LHX6/NFIA | pM1 | 14 | 3.9% | 29.2% |
| module_5 | CRABP1/ANGPT2 | pM2 | 4 | 1.1% | 10.5% |

pM-specific interpretation:

- The strongest descriptive pM signal is `pM1` for Jia cortical/MGE interneuron
  modules `EPHA5/MEF2C` and `LHX6/NFIA`.
- `NR2F1/NR2F2` points most toward `pM3` among the pM sets.
- `LHX8/ISL1` points most toward `pM2` among the pM sets, but only with a small
  absolute overlap.
- `CRABP1/ANGPT2` is weak against the pM marker sets; its best pM hit is `pM2`
  with four genes. This supports keeping it as a bridge/contested lineage
  rather than treating it as a clean pM-only class.
- All pM overlaps should be read descriptively because pM marker sets are small
  and the BH-adjusted hypergeometric q-values are not significant in this
  pM-only universe.

Collapse decision from the pM-only analysis:

| Question | Decision | Rationale |
|---|---|---|
| Collapse Shi `pM1-pM4` into one generic pM class? | No for analysis; yes only as a parent label. | The pM sets are distinct marker lists, and Jia modules point to different pM clusters. Keep `pM1-pM4` separate in figures/tables, with "Shi pM/MGE progenitor" only as a parent tier. |
| Collapse Jia `EPHA5/MEF2C` and `LHX6/NFIA`? | Yes for a reader-facing parent class; no for module-resolution interpretation. | Both best match `pM1`, with 11 and 14 overlapping genes, respectively. This supports a shared pM1-associated cortical/MGE interneuron-output tier, but Jia's split remains biologically useful. |
| Collapse Jia `LHX8/ISL1` and `NR2F1/NR2F2`? | No. | Their best pM hits differ: `LHX8/ISL1` points to `pM2`, while `NR2F1/NR2F2` points to `pM3`. Keep them separate, with "early MGE/subpallial output" only as a parent tier. |
| Collapse Jia `CRABP1/ANGPT2` into pM-only biology? | No. | Its pM overlap is weak: best hit `pM2`, 4 genes, 1.1% of the Jia module. Keep as bridge/contested. |

Best pM overlap genes per Jia full module:

| Jia module | Jia exact name | Best pM | n overlap | Overlapping genes |
|---|---|---:|---:|---|
| module_1 | LHX8/ISL1 | pM2 | 6 | ASCL1, CDK6, DLX1, HMGA1, NKX2-1, SMS |
| module_2 | NR2F1/NR2F2 | pM3 | 7 | ASCL1, GSX1, HES5, HMGA1, NKX2-1, OLIG1, SOX2 |
| module_3 | EPHA5/MEF2C | pM1 | 11 | BCL11A, DCX, DLX6, DLX6-AS1, HIST1H4C, LHX6, PBX1, PLXNA2, SOX11, SOX4, TIAM1 |
| module_4 | LHX6/NFIA | pM1 | 14 | BCL11A, BEST3, CDCA7, DCLK2, DLX6, LHX6, PDZRN4, PLS3, PLXNA4, RAB3IP, SOX11, SOX4, SP9, TCF4 |
| module_5 | CRABP1/ANGPT2 | pM2 | 4 | ASCL1, DLX1, HES6, NNAT |

All full-module pM overlaps:

| Jia module | pM cluster | n overlap | Overlapping genes |
|---|---|---:|---|
| module_1 LHX8/ISL1 | pM1 | 6 | DLX6, HMGA1, MEG3, NKX2-1, SOX11, SOX4 |
| module_1 LHX8/ISL1 | pM2 | 6 | ASCL1, CDK6, DLX1, HMGA1, NKX2-1, SMS |
| module_1 LHX8/ISL1 | pM3 | 6 | ASCL1, CDK6, HES5, HMGA1, NKX2-1, OLIG1 |
| module_1 LHX8/ISL1 | pM4 | 6 | HES5, HSPB1, MEG3, NES, NKX2-1, OLIG1 |
| module_2 NR2F1/NR2F2 | pM1 | 5 | HMGA1, NKX2-1, SOX11, SOX4, TAGLN3 |
| module_2 NR2F1/NR2F2 | pM2 | 6 | ASCL1, GADD45G, GSX1, HMGA1, NKX2-1, RGS16 |
| module_2 NR2F1/NR2F2 | pM3 | 7 | ASCL1, GSX1, HES5, HMGA1, NKX2-1, OLIG1, SOX2 |
| module_2 NR2F1/NR2F2 | pM4 | 5 | HES5, KCNQ1OT1, NKX2-1, OLIG1, SOX2 |
| module_3 EPHA5/MEF2C | pM1 | 11 | BCL11A, DCX, DLX6, DLX6-AS1, HIST1H4C, LHX6, PBX1, PLXNA2, SOX11, SOX4, TIAM1 |
| module_3 EPHA5/MEF2C | pM2 | 5 | ASCL1, DLX1, GSX1, NELL2, TFDP2 |
| module_3 EPHA5/MEF2C | pM3 | 7 | ASCL1, GSX1, NELL2, PID1, PLXNA2, SLC44A1, TFDP2 |
| module_3 EPHA5/MEF2C | pM4 | 3 | NTN4, SLC44A1, TFDP2 |
| module_4 LHX6/NFIA | pM1 | 14 | BCL11A, BEST3, CDCA7, DCLK2, DLX6, LHX6, PDZRN4, PLS3, PLXNA4, RAB3IP, SOX11, SOX4, SP9, TCF4 |
| module_4 LHX6/NFIA | pM2 | 3 | ASCL1, DLX1, GADD45G |
| module_4 LHX6/NFIA | pM3 | 1 | ASCL1 |
| module_4 LHX6/NFIA | pM4 | 1 | BCAN |
| module_5 CRABP1/ANGPT2 | pM1 | 4 | BCL11A, LHX6, SOX11, SOX4 |
| module_5 CRABP1/ANGPT2 | pM2 | 4 | ASCL1, DLX1, HES6, NNAT |
| module_5 CRABP1/ANGPT2 | pM3 | 2 | ASCL1, NNAT |
| module_5 CRABP1/ANGPT2 | pM4 | 4 | KCNQ1OT1, NTM, VEPH1, WWTR1 |

Genes shared with any pM marker set by Jia full module:

| Jia module | Shared with any pM / module genes | Genes |
|---|---:|---|
| module_1 LHX8/ISL1 | 14 / 298 | ASCL1, CDK6, DLX1, DLX6, HES5, HMGA1, HSPB1, MEG3, NES, NKX2-1, OLIG1, SMS, SOX11, SOX4 |
| module_2 NR2F1/NR2F2 | 13 / 118 | ASCL1, GADD45G, GSX1, HES5, HMGA1, KCNQ1OT1, NKX2-1, OLIG1, RGS16, SOX11, SOX2, SOX4, TAGLN3 |
| module_3 EPHA5/MEF2C | 19 / 228 | ASCL1, BCL11A, DCX, DLX1, DLX6, DLX6-AS1, GSX1, HIST1H4C, LHX6, NELL2, NTN4, PBX1, PID1, PLXNA2, SLC44A1, SOX11, SOX4, TFDP2, TIAM1 |
| module_4 LHX6/NFIA | 18 / 362 | ASCL1, BCAN, BCL11A, BEST3, CDCA7, DCLK2, DLX1, DLX6, GADD45G, LHX6, PDZRN4, PLS3, PLXNA4, RAB3IP, SOX11, SOX4, SP9, TCF4 |
| module_5 CRABP1/ANGPT2 | 12 / 351 | ASCL1, BCL11A, DLX1, HES6, KCNQ1OT1, LHX6, NNAT, NTM, SOX11, SOX4, VEPH1, WWTR1 |

## Jia S9 Values For DLX6, LHX8, And LHX6

These values come directly from the Jia Science Data S9 workbook. The
expression-like workbook field is `vst.mean`; this table does not recompute
expression from an expression matrix.

Machine-readable table:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/jia_s9_shi_lineage_overlap/jia_s9_shi_s6_pm_only_v1/tables/jia_s9_dlx6_lhx8_lhx6_values.tsv
```

Compact `vst.mean` matrix across Jia full modules:

| Gene | module_1 LHX8/ISL1 | module_2 NR2F1/NR2F2 | module_3 EPHA5/MEF2C | module_4 LHX6/NFIA | module_5 CRABP1/ANGPT2 |
|---|---:|---:|---:|---:|---:|
| DLX6 | 0.00614754 | absent | 0.0929846 | 0.132984 | absent |
| LHX8 | absent | 0.00301023 | absent | absent | absent |
| LHX6 | absent | absent | 0.0787723 | 0.156021 | 0.0594796 |

Detailed Jia S9 rows where each gene is present:

| Gene | Jia module | Module name | Sheet | Rank | pval | qval | vst.mean | vst.variance.standardized |
|---|---|---|---|---:|---:|---:|---:|---:|
| DLX6 | module_1 | LHX8/ISL1 | Module 1 | 115 | 0.398762 | 0.542677 | 0.00614754 | 1.28992 |
| DLX6 | module_3 | EPHA5/MEF2C | Module 3 | 77 | 6.34688e-61 | 1.24377e-60 | 0.0929846 | 1.73903 |
| DLX6 | module_4 | LHX6/NFIA | Module 4 | 87 | 1.59176e-19 | 4.29543e-19 | 0.132984 | 2.71327 |
| LHX8 | module_2 | NR2F1/NR2F2 | Module 2 | 29 | 0.0031442 | 0.00868982 | 0.00301023 | 1.29132 |
| LHX6 | module_3 | EPHA5/MEF2C | Module 3 | 53 | 1.66487e-92 | 4.08697e-92 | 0.0787723 | 2.13791 |
| LHX6 | module_4 | LHX6/NFIA | Module 4 | 71 | 2.44583e-20 | 6.82129e-20 | 0.156021 | 1.70988 |
| LHX6 | module_5 | CRABP1/ANGPT2 | Module 5 | 106 | 3.30219e-18 | 6.33465e-18 | 0.0594796 | 2.10795 |

All three genes are absent from the Jia S9 TF-only sheets in this workbook.

## Shi S6 pM Values For DLX6, LHX8, And LHX6

These values come directly from Shi Table S6, restricted to literal Shi
clusters `pM1-pM4`. Do not compare these columns directly to Jia `vst.mean`;
Shi Table S6 reports marker statistics for each pM cluster (`avg_logFC`,
`pct.1`, `pct.2`, `p_val_adj`).

Machine-readable table:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/jia_s9_shi_lineage_overlap/jia_s9_shi_s6_pm_only_v1/tables/shi_s6_pm_dlx6_lhx8_lhx6_values.tsv
```

Compact Shi pM marker-status matrix:

| Gene | pM1 | pM2 | pM3 | pM4 |
|---|---|---|---|---|
| DLX6 | present | absent | absent | absent |
| LHX8 | absent | absent | absent | absent |
| LHX6 | present | absent | absent | absent |

Detailed Shi Table S6 pM rows where each gene is present:

| Gene | Shi pM cluster | avg_logFC | pct.1 | pct.2 | p_val | p_val_adj |
|---|---|---:|---:|---:|---:|---:|
| DLX6 | pM1 | 0.323791 | 0.386 | 0.244 | 2.82215e-31 | 6.01399e-27 |
| LHX6 | pM1 | 0.731337 | 0.432 | 0.032 | 0 | 0 |

Interpretation for the pM comparison:

- In Shi Table S6, `DLX6` and `LHX6` are pM1 marker genes, not broad pM1-pM4
  marker genes.
- `LHX8` is not a Shi Table S6 pM1-pM4 marker in this workbook.
- This reinforces the earlier pM-specific conclusion: the Jia cortical/MGE
  interneuron modules, especially `LHX6/NFIA`, are the cleanest pM1-aligned
  signal; the `LHX8/ISL1` Jia module should not be interpreted as a direct
  Shi pM marker match based on LHX8 itself.

## Exhaustive S3-S9 Result Snapshot

The exhaustive run compared Jia Science Data S9 against Shi Tables S3-S9.

Headline counts:

```text
Jia unique module genes: 1156
Shi unique marker genes: 1396
Shi marker sets: 49
Shared Jia/Shi genes: 247
Jia-only genes in this comparison: 909
Shi-only genes in this comparison: 1149
Literal curated Shi M labels present as workbook cluster labels: 0/6
```

Full-module recovery against any Shi S3-S9 marker set:

| Jia module | Jia exact name | Shared / Jia genes | Fraction shared |
|---|---|---:|---:|
| module_1 | LHX8/ISL1 | 78 / 298 | 26.2% |
| module_2 | NR2F1/NR2F2 | 50 / 118 | 42.4% |
| module_3 | EPHA5/MEF2C | 76 / 228 | 33.3% |
| module_4 | LHX6/NFIA | 83 / 362 | 22.9% |
| module_5 | CRABP1/ANGPT2 | 71 / 351 | 20.2% |

Best observed full-module hits:

| Jia module | Best Shi workbook/context/cluster | Overlap | q-value |
|---|---|---:|---:|
| module_1 LHX8/ISL1 | S6 GE progenitor LGE `pL2` | 11 genes | 0.345 |
| module_2 NR2F1/NR2F2 | S6 GE progenitor LGE `pL2` | 8 genes | 0.0668 |
| module_3 EPHA5/MEF2C | S9 integrated human/mouse subcluster `5` | 13 genes | 0.0154 |
| module_4 LHX6/NFIA | S9 integrated human/mouse subcluster `2` | 15 genes | 0.0239 |
| module_5 CRABP1/ANGPT2 | S6 GE progenitor LGE `pL3` | 7 genes | 0.797 |

Interpretation:

- The Shi workbook labels and Jia module names are not directly the same naming
  system in the available S3-S9 Excel tables. Shi's curated `M2/M3/M4/M5/M6/M7`
  concepts from the notes are not literal cluster labels in these inputs, even
  though all Shi S2-S9 Excel files are accessible and readable.
- Gene-list overlap is present but partial. It supports a reader-facing
  conceptual crosswalk, not a direct one-to-one label merge.
- The strongest computed support for collapsing is the cortical pair:
  Jia `EPHA5/MEF2C` and `LHX6/NFIA` both point toward integrated human/mouse
  Shi subclusters, but their direct Jia-module overlap is modest
  (32 genes; Jaccard 0.0573), so keep the Jia split when discussing
  higher-resolution/chandelier-biased biology.
- Do not collapse `LHX8/ISL1` and `NR2F1/NR2F2` except as a parent
  "early subpallial output" tier; the biology is cholinergic vs GABAergic.
- Do not collapse `CRABP1/ANGPT2` into subpallial-only. Its weak best computed
  overlap and the Shi anatomical interpretation argue for keeping it as the
  CRABP1 bridge/contested lineage.

## What The Workflow Computes

1. Reads all Jia S9 full module sheets and TF-only sheets.
2. Reads one or more Shi Excel marker tables with automatic header-row
   detection.
3. Normalizes genes to uppercase symbols for overlap.
4. Computes pairwise Jia-module vs Shi-cluster overlap:
   - overlap gene count
   - Jaccard index
   - overlap coefficient
   - fraction of Jia module recovered in Shi cluster
   - fraction of Shi cluster recovered in Jia module
   - hypergeometric p-value and BH-adjusted q-value
5. Writes a curated crosswalk table that checks whether expected Shi labels are
   actually present in the input workbook.
6. Counts shared vs Jia-only vs Shi-only genes.
7. Audits anchor genes from the curated crosswalk.
8. Computes internal Jia module overlap to support collapse/not-collapse logic.
9. Writes a collapse-guidance table that separates reader-facing simplification
   from Jia-resolution lineage interpretation.
10. Plots full-module and TF-only overlap heatmaps/dotplots.

## Output Contract

Tables:

```text
tables/jia_s9_module_genes.tsv
tables/shi_marker_genes_long.tsv
tables/jia_s9_shi_pairwise_overlap.tsv
tables/jia_shi_crosswalk_marker_presence.tsv
tables/jia_s9_gene_set_summary.tsv
tables/shi_gene_set_summary.tsv
tables/jia_gene_membership_all_sources.tsv
tables/jia_module_shared_private_summary.tsv
tables/jia_shared_genes_by_module.tsv
tables/jia_private_genes_by_module.tsv
tables/shi_marker_set_recovery_by_jia.tsv
tables/jia_internal_module_overlap.tsv
tables/jia_shi_anchor_marker_audit.tsv
tables/jia_shi_naming_convention_audit.tsv
tables/jia_shi_collapse_guidance.tsv
tables/jia_s9_shi_overlap_run_manifest.tsv
```

Plots:

```text
plots/jia_full_vs_shi_overlap_coefficient_heatmap.png
plots/jia_full_vs_shi_overlap_coefficient_heatmap.pdf
plots/jia_tf_vs_shi_overlap_coefficient_heatmap.png
plots/jia_tf_vs_shi_overlap_coefficient_heatmap.pdf
plots/jia_full_vs_shi_overlap_dotplot.png
plots/jia_full_vs_shi_overlap_dotplot.pdf
plots/jia_tf_vs_shi_overlap_dotplot.png
plots/jia_tf_vs_shi_overlap_dotplot.pdf
```

Report:

```text
reports/jia_s9_shi_lineage_overlap_report.md
```

## Slurm Command

```bash
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
cp slurm_templates/45_jia_s9_shi_lineage_overlap.sbatch.template \
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/45_jia_s9_shi_lineage_overlap.sbatch
sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/45_jia_s9_shi_lineage_overlap.sbatch
```

To compare multiple Shi workbooks in one run, separate them with `:`:

```bash
SHI_XLSX="/path/to/shi_table_a.xlsx:/path/to/shi_table_b.xlsx" \
RUN_LABEL=jia_s9_shi_multi_table_overlap_v1 \
sbatch slurm_templates/45_jia_s9_shi_lineage_overlap.sbatch.template
```

Focused pM1-pM4 command:

```bash
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
RUN_LABEL=jia_s9_shi_s6_pm_only_v1 \
SHI_XLSX=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/reference/shi_2021_tables_s2_to_s9/science.abj6641_table_s6.xlsx \
INCLUDE_SHI_CLUSTER_REGEX='^pM[1-4]$' \
sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/45_jia_s9_shi_lineage_overlap_pm_only_v1.sbatch
```

Exhaustive S3-S9 run command used for the broad comparison:

```bash
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
cp slurm_templates/45_jia_s9_shi_lineage_overlap.sbatch.template \
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/45_jia_s9_shi_lineage_overlap_exhaustive_v1.sbatch
RUN_LABEL=jia_s9_shi_s3_to_s9_exhaustive_v1 \
SHI_XLSX="/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/reference/shi_2021_tables_s2_to_s9/science.abj6641_table_s3.xlsx:/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/reference/shi_2021_tables_s2_to_s9/science.abj6641_table_s4.xlsx:/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/reference/shi_2021_tables_s2_to_s9/science.abj6641_table_s5.xlsx:/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/reference/shi_2021_tables_s2_to_s9/science.abj6641_table_s6.xlsx:/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/reference/shi_2021_tables_s2_to_s9/science.abj6641_table_s7.xlsx:/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/reference/shi_2021_tables_s2_to_s9/science.abj6641_table_s8.xlsx:/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/reference/shi_2021_tables_s2_to_s9/science.abj6641_table_s9.xlsx" \
sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/45_jia_s9_shi_lineage_overlap_exhaustive_v1.sbatch
```

## Interpretation Rules

- Treat high overlap as gene-list support, not as proof of developmental
  equivalence.
- Keep Jia full modules and Jia TF-only modules separate.
- Do not force Module 5 CRABP1/ANGPT2 into "subpallial-only"; preserve it as a
  contested bridge until Shi anatomical evidence and Jia trajectory evidence
  are reconciled.
- If expected Shi M labels are missing from the input workbook, use the computed
  best-observed Shi cluster as a discovery result, not as the curated crosswalk.
