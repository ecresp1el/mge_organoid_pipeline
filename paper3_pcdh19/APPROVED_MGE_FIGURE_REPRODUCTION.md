# Approved mouse MGE mosaic figure — 1 November 2026

This is the exact eight-panel figure approved in the conversation, titled
**Mouse MGE · E14.5 · n = 12**. Here n means 12 registered dissected MGE
samples, not the 72,000 displayed cells. It is mouse tissue, not organoids.
The final renderer is `scripts/render_HET_palette_expression.py`.

## Locations and workflow

Code, configuration and this handoff stay in Git. All inputs, intermediate
results, figures, submitted jobs and logs belong outside the checkout:

- Study runtime: `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19`
- Frozen approved package: `final_figures/approved_mge_mosaic_20261101/` under that runtime.
- Original exploratory run: `results/female_male_signaling_20261031/` under that runtime.
- Approved image: `approved/HET_palette_expression_figure.png` inside the frozen package; PDF and editable SVG sit beside it.
- Executed reproduction job: `jobs/approved_mge_mosaic_20261101/reproduce.sbatch`.
- Test log: `logs/approved_mge_reproduction-65177019.log`.

The frozen package contains `code/` (source/configuration snapshots),
`code_before_path_fix/` (historical source before moving runtime paths),
`inputs/` (all nine direct rendering inputs), `approved/` (reference exports),
`input_manifest.json`, `frozen_manifest.json`, `migration_manifest.json`,
`migration_summary.json`, and `historical_jobs_and_logs.tsv`.
Historical jobs/logs were preserved with their original run directories;
their embedded checkout paths describe past execution, not current entrypoints.
Do not submit those historical jobs unchanged. Use the reproduction job below.

## Exact reproduction

From the repository root, submit:

```bash
sbatch --wait paper3_pcdh19/slurm/reproduce_approved_mge_figure.sbatch
```

The wrapper verifies all frozen checksums, copies the nine inputs to a new
NFS reproduction directory, and runs the frozen renderer. It does not
recompute UMAPs, scores, differential expression or interaction selection.
It leaves the approved package untouched and writes a
`REPRODUCTION_VALIDATION.json` in the new runtime directory. Each invocation
uses a new UTC timestamp; `--output-root /exact/new/nfs/path` can be passed
to the sbatch script for a specific destination outside the checkout.

Verified on Slurm job **65177019**: the rendered PNG and gene-value table
both match the approved references **byte for byte**. Test directory:
`results/approved_mge_reproduction_20261101_155151/` under the study runtime.
PDF/SVG timestamps and generated element identifiers are not byte-compared.

Approved PNG SHA-256:
`7256cfe0433a61a595593ad36189edfcac9c3f7752738ad4008f01d3b65fd0c7`

Approved gene-value table SHA-256:
`3abe9efaa4944e4d2ce3589e2bd22958593dea1a1bd34a94310a74e656c236b9`

Environment: Python 3.6.8, NumPy 1.19.5, SciPy 1.5.4, matplotlib 3.3.4;
DejaVu Sans bundled with that matplotlib installation. The sbatch template
sets PYTHONPATH to the existing `software/gse94641_e15_5_label_transfer/python3.6-site-packages`
under the study runtime, loads R/4.4.3, and uses two CPUs and 8 GB RAM.
Canvas: 7.2 × 8.49 inches; PNG 450 dpi (3240 × 3820 pixels); PDF text type 42.

## Nine frozen direct inputs

Paths below are relative to `inputs/results/` inside the package. Every file
has a source path, byte count and SHA-256 in `input_manifest.json`.

| Panel | Input |
| --- | --- |
| A | `female_male_signaling_20261031/atlas_sex_genotype_display.tsv.gz` |
| B | `female_male_signaling_20261031/saved_MGE_LGE_CGE_scores.npz` |
| C | `neuronal_matched_selection_20261031/neuronal_exon1_signal.npz` |
| D | `mosaic_deep_dive_20261031/DE_C0006/within_HET_Aneg_BC1_vs_Apositive.tsv.gz` |
| E | `female_male_signaling_20261031/male_total_Pcdh19_display.npz` |
| F | `female_male_signaling_20261031/male_C0006_DE.tsv` |
| D/F | `female_male_signaling_20261031/matched_display_logCPM.tsv` |
| D/F | `female_male_signaling_20261031/matched_female_display_metadata.tsv` |
| G/H | `female_male_signaling_20261031/HET_replicated_interaction_selection.tsv` |

## Panel definitions and display choices

A reuses the saved joint UMAP, 6,000 cells per sample, seed 20261031.
The categorical mapping is provisional: RG=C0002; MGE IPC=C0003/C0008;
neurogenic progenitors=C0001/C0007/C0012/C0024; cycling neuroblasts=C0011;
immature inhibitory neurons=C0004/C0005/C0006/C0020; mixed neural/vascular=C0015;
other clusters gray. Atlas colors are distinct from the five sex/genotype
inset colors. The main renderer stores all hex values, geometry and fonts.

B uses the original saved MGE/LGE/CGE marker-program scores, joined to the
same cells. Score = mean ln(1+CPM) marker expression minus the mean of
expression-bin-matched controls. Scores are not 0–1 probabilities. Only
scores strictly >0.8 are colored, using Greens on a shared 0.8–1.5 scale;
values ≥1.5 saturate for display and remaining cells are gray. No refitting.
The original scorer/configuration are `scripts/hicat/validation_programs.py`
and `config/hicat_validation_programs.json`.

C/E both display ln(1 + feature UMI / total cell UMI × 10,000), with identical
0–3 axes and KDE bandwidth 0.18. Mean dots and interquartile strips use those
transformed values. Observations above 3 remain in density/statistic calculations.
C measures exon-1 probe A in C0006: WT females A-positive (1,992 cells),
HET WT cells A-positive (495), HET KO cells A=0 and B+C≥1 (2,180).
C's zero HET KO exon-1 signal follows directly from that selection.
E measures total Pcdh19=A+B+C in C0006: WT males are shown only if total>0
(2,094/20,981), whereas all 19,097 KO male cells are shown, including zeros.
This asymmetric positive selection is a plotting choice only; F uses all
eligible male cells. C and E therefore share units but use different Pcdh19
features and selections. Do not interpret their distributions as a matched
cross-sex knockout validation.

D/F share 14 genes, their hierarchy/order and scales. Joint TMM normalization
across 15 sample/subset pseudobulks gives logCPM (prior.count=2); each gene is
standardized across individual sample libraries within sex, then equally
weighted sample means are displayed. Average-linkage clustering uses the
concatenated female/male standardized profiles. Purple scale −1.5 to +1.5.
Bars show existing edgeR log2FC, range −1.6 to +1.6: **D is HET KO cells vs
HET WT cells, paired within the three HET females; the pure WT female column
is context, not the bar reference. F is KO males vs WT males, three samples
each.** Blue means lower and red higher in the KO group versus that WT
reference. Keys are at the far right. The male Pcdh19-positive display filter
is not applied to F. Candidate genes were explored/selected after screening;
these are not gene-wide FDR-significant discoveries.

G shows exactly the same five selected interactions as H; arrows indicate
ligand source to receptor receiver, loops within a population. They do not
prove physical cell contact. H compares KO→KO versus WT→WT within HET
females. Complex subunits are combined by geometric mean; expression
potential is sqrt(ligand CPM × receptor CPM). Every component must be ≥1 CPM
in all three HET samples in both states. All three changes must agree in
direction with absolute log2 change ≥0.2. Rank by the smallest absolute
replicate effect, keep unique interactions, at most two protocadherin pairs
and one per other predefined neural pathway family; display five.
Percent change = 100 × (2^median paired log2 change − 1).
G width = 0.6 + 2.2 × absolute percent change / maximum absolute change pt.
**No depth matching is used in the approved figure.** Archived matched-depth
audits remain diagnostic history, not the source of plotted values.
“WT cells” and “KO cells” within HET are display aliases for A-positive and
A-negative/downstream-positive probe states, not independently measured DNA
genotypes. These are inferred expression-based interaction candidates, not
measured signaling or a fitted CellChat communication-probability model.

## Upstream provenance and code

All paths below are relative to the study runtime unless stated otherwise.
The frozen inputs suffice for exact rendering; rebuilding from raw checkpoints
is a separate scientific analysis, not required for reproducing this figure.

- Raw counts and joint embedding:
  `results/primary_processing/06_technical_sample_batch_diagnostics/06_technical_sample_batch_diagnostics_20260903_114110_6087eeb/objects/pcdh19_step06_unintegrated_diagnostics.h5ad`.
- Consensus assignments:
  `results/hicat/03_full_data_consensus_benchmark/hicat_consensus_review_20260909_142348_efd9e13d/outputs/tables/all_cell_assignments.tsv.gz`.
- Original regional scores:
  `results/hicat/02_hierarchy_validation/hicat_validation_20260908_221625_d5d06a5d/outputs/full_data_projection/plot_program_scores.tsv.gz`.
- Raw probe A/B/C evidence:
  `results/pcdh19_probe_audit/per_sample/15662-JZ-N/pcdh19_probe_patterns.tsv` (N=1–12).
- Per-cell diagnostic cache: `results/mosaic_deep_dive_20261031/per_cell_marker_diagnostics.h5`.
- Sample metadata: repository `config/sample_key.csv`: JZ1–3 WT male,
  JZ4–6 WT female, JZ7–9 HET female, JZ10–12 KO male.
- Ligand/receptor resource: `results/paired_HET_signaling_20261031/expanded_signaling_pairs.tsv`
  and `signaling_gene_expression_support.tsv`; official CellChatDB.mouse
  interaction definitions expanded into complex subunits by `paired_HET_signaling.R`.

Source chain (all scripts are in the repository and frozen package):

1. `mosaic_deep_dive_prepare.py` / `mosaic_deep_dive_fit.R` produce probe-gated
   pseudobulks and female paired DE; associated report/finalize scripts preserve
   all tested definitions, including null results.
2. `extract_neuronal_exon1.py`, `extract_neuronal_violin_cells.py`,
   `prepare_neuronal_matched_cells.py` prepare female display caches.
3. `prepare_atlas_sex_genotype.py` and `prepare_saved_regional_scores.py`
   prepare A/B inputs without recomputing embeddings or program scores.
4. `prepare_male_pseudobulks.py`, `fit_male_pseudobulks.R`,
   `prepare_male_total_Pcdh19.py` prepare male DE and E.
5. `prepare_matched_figure_readouts.R` prepares the shared TMM display matrix.
6. `paired_HET_signaling.R`, `prepare_HET_crossstate_networks.py`,
   `audit_HET_network_differences.py`, and `render_HET_network_audit.py`
   supply interaction evidence and the frozen candidate selection. Audit code
   also computes depth-matched diagnostics; those diagnostics are not used here.
7. `render_HET_palette_expression.py` renders this approved figure.
   Older `render_*` scripts preserve iteration history and are not the final entrypoint.

Some historical validation/legend files retain earlier panel lettering or
superseded descriptions. This handoff describes the approved final layout;
the frozen renderer and nine inputs are authoritative for exact reproduction.

## Migration and Git audit

The reported 887 were untracked files, not commits. On 1 November 2026,
821 runtime files (483,698,254 bytes) across 13 run folders were copied to
`results/` on NFS, SHA-256 verified, then removed from the checkout. No
existing NFS run names collided and no Git history was rewritten.
`migration_manifest.json` records every original and destination path/hash.
55 new Python scripts now use `pcdh19_exploratory_paths.RUNTIME_ROOT` for
runtime results; their code/config lookups remain in the repository.
`PCDH19_RUNTIME_ROOT` can override the default with an external path; checkout
runtime roots are rejected. `.gitignore` excludes Paper 3 results/jobs/logs/final_figures.

## Display clarity revision (1 November 2026)

Subsequent requested display edits use `scripts/render_HET_display_clarity.py`.
Outputs: `results/female_male_signaling_20261031/HET_display_clarity_figure.{png,pdf,svg}`
on NFS. Source snapshot, exports and `PROVENANCE.json` are also saved under
`final_figures/mge_display_clarity_20261101/`. Job 65177438 completed successfully;
submitted job: `jobs/approved_mge_mosaic_20261101/display_clarity.sbatch`.
The nine direct scientific inputs and all D/F gene values are unchanged.
A now states 12 samples from 12 mice and 3 each WT female/HET female/WT male/KO male;
inset dot area is 0.045 rather than 0.19 pt². B removes cap/cutoff text while
retaining the >0.8 overlay rule and 1.5 color saturation. C uses a clear vertical
zero marker and “0 exon-1 signal” for the all-zero probe-gated HET KO group.
The previously approved frozen package and its exact reproduction command remain intact.

## Exploratory activity / IEG bottom row (1 November 2026)

`render_HET_activity_IEG_final.py` extends the display-clarity figure with a
full-width final-row panel I. Submit with
`sbatch --wait paper3_pcdh19/slurm/render_HET_activity_IEG_final.sbatch`.
The new package is `final_figures/mge_activity_IEG_20261101/` under the study
runtime; the combined exports are `HET_activity_IEG_final.{png,pdf,svg}`.
Prior approved and display-clarity exports are preserved. The base renderer
runs against copies of its seven OUT-directory inputs in the new package;
its other source paths remain unchanged. The package retains source snapshots,
paired IEG tables, a combined legend, and a SHA-256 provenance manifest.

I shows Egr1, Npas4 and Nr4a1 pseudobulk CPM, paired across three HET embryos,
in independently gated cycling-marker-negative inhibitory neurons. The
exploratory Npas4 elevation is prominently qualified: 0 WT / 9 KO UMIs from
16 WT / 151 KO cells, and the pattern weakens with a relaxed cycle gate.
WT/KO denote the existing expression-based probe states. These data do not
establish intracellular signaling activation, cortical residence, or a robust
IEG mechanism. Each gene has its own y scale; no significance test is added.

### Corrected four-row IEG dot layout

The user rejected the added fifth row. The current renderer keeps the original
7.2 × 8.49 inch canvas and panels A–F, and fits G/H/I across row four.
I is a compact Egr1/Npas4/Nr4a1 circle plot: area shows percent of cells
expressing; shared purple color shows mean log-normalized expression, including
zeros. Zero detection uses a hollow circle. WT/KO rows pool three HET samples
for descriptive display only. Source cell gate and sparse-count limitations
are unchanged. Updated exports live in `final_figures/mge_activity_IEG_dots_20261101/`
under the study runtime; `IEG_dot_values.tsv` retains the plotted values.
The earlier five-row package is a superseded draft.

### Typography and expanded neuronal IEG panel

The current four-row renderer standardizes panel letters to 11 pt, headings to
8 pt, body labels to 6–7 pt, and the dense atlas legend to 5.5 pt. Panel I now
uses nine gene rows and two probe-state columns (WT/KO): Npas4, Fos, Fosb, Junb,
Egr1, Egr2, Egr3, Arc and Nr4a1. All nine were already in the cached IEG screen;
none are removed for zero detection or selected by effect direction. The circle
size key is removed as requested; size still encodes detection percent, color
still encodes shared-scale mean log expression, and zero detection is hollow.
The complete values and definitions remain in the package TSV and legend.

Biological context for the gene panel: Madabhushi et al., Cell (2015),
https://pmc.ncbi.nlm.nih.gov/articles/PMC4886855/; The Neuroplasticity-Associated Arc Gene Is a Direct Transcriptional Target
of Early Growth Response (Egr) Transcription Factors, https://pmc.ncbi.nlm.nih.gov/articles/PMC1291244/;
Subregion-specific rules govern the distribution of neuronal immediate-early
gene induction,
https://pmc.ncbi.nlm.nih.gov/articles/PMC7519336/.
These references establish activity-associated neuronal IEG context, not
validation of activity or mechanism in this sparse E14.5 MGE subset.

### Explicit ligand–receptor labels and male IEG columns

G/H are now titled “Ligand–receptor network” and “Ligand–receptor changes”.
Subtitles identify the female HET comparison and source-to-receiver direction;
H labels its axis as change in inferred L–R strength. They retain the same
expression-based interaction candidates and values, including curated contact
interactions; this relabeling does not establish measured communication.

I removes the exploratory artwork footer and now shows four columns grouped
as HET female (WT cells, KO cells) and males (WT, KO). Known male samples use
the same neuronal identity/core/cycling-marker gate, with all eligible male
cells retained, without a Pcdh19 expression filter. The HET probe-state gates
are unchanged. This differs intentionally from the positive-only WT-male
Pcdh19 violin display in E, whose visualization filter is not used for IEGs.
`IEG_dot_per_sample_values.tsv` preserves all sample-level readouts alongside
the pooled display values. Scientific limitations remain in the full legend.

### Revised C/E wording and expression-scaled, ranked IEG circles

C now labels WT Females, HET WT Females and HET KO Females; the zero-signal
text is removed. E labels only WT Males and KO Males. These wording changes
do not change the documented probe states or E's positive-only WT display.

I no longer uses detection for circle size. For each gene and sample, compute
mean ln(1 + raw UMI / total UMI * 10000), including zeros; average the three
sample means equally per group, then divide by that gene's maximum over the
four groups. Color and area encode this 0–1 relative expression; positive area
is 6 + 30 * relative expression pt², with a 5 pt² hollow circle for zero.
Circle sizes are chosen for readability and do not represent cell fractions.
Each gene has its own scale, so absolute expression is not comparable across
rows; sparse genes may have a large relative circle. Raw values are preserved.

Rank rows by descending average of female-HET and male KO-minus-WT scaled
changes, weighting both comparisons equally. Alphabetical tie-breaking;
all-zero genes last. This is descriptive effect ranking, not a significance
ranking. All nine cached IEGs remain, including opposing changes and zeros.
Per-sample tables, display scores, relative values and row ranks are exported.

### Supplement expression with aligned KO-minus-WT circles

I retains its four expression columns and adds two signed change columns,
one for HET females and one for males, using the same gene rows and order.
Change is KO minus WT on the existing per-gene relative-expression scale,
not log fold change or statistical significance. Red means higher KO, blue
lower KO. Area encodes absolute change (6+30*abs(change) pt²), with hollow
zero references; color spans −1 to +1 consistently across both comparisons.
`IEG_KO_minus_WT_changes.tsv` preserves both scaled changes and unscaled
mean-log-expression differences. I is moved left and widened within row four;
H is compacted horizontally. Title: “Activity-dependent gene expression”.
The title identifies the gene class; there is no activity perturbation here.

### Final wording preference

User finalized the figure and requested singular sex labels everywhere.
All visible “Males/males/Females/females” are now singular. C uses exactly
“WT Female”, “HET Female WT cells”, “HET Female KO cells”; its longer group
labels use 6 pt text. Underlying probe gates, measurements, and results remain
unchanged. Preserve this singular-sex wording in subsequent figure revisions.
