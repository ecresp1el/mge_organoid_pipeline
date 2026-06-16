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
pL1-pL3, and pM1-pM4. It does not directly expose the M2/M3/M5/M6/M7 names
from the biological crosswalk. The workflow therefore computes overlap against
the actual input clusters and separately marks expected crosswalk labels as
present or missing.
```

If a later Shi workbook with literal M2/M3/M4/M5/M6/M7 lineage labels is added,
pass it as an additional or replacement `--shi-xlsx`.

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
  concepts from the notes are not literal cluster labels in these inputs.
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
