# Final mouse MGE figure — 1 November 2026

Finalized by the user after the panel-I scale-spacing correction (job 65184824).
The final figure has four rows: A/B, C/D, E/F, G/H/I. Do not substitute the
older eight-panel approved package or the superseded five-row IEG draft.

## Saved package

NFS study runtime:
`/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19`

Frozen final package under that runtime:
`final_figures/mge_activity_IEG_FINAL_20261101/`

- `approved/HET_activity_IEG_final.{png,pdf,svg}`: final user-visible exports.
- `approved/HET_ACTIVITY_IEG_LEGEND.txt`: full scientific definitions.
- `approved/IEG_dot_values.tsv`, `IEG_dot_per_sample_values.tsv`,
  `IEG_KO_minus_WT_changes.tsv`: expression, sample and change readouts.
- `inputs/`: independent copies of all direct renderer inputs and copied IEG
  metadata; original NFS source paths and hashes are in `input_manifest.json`.
- `code/`: Python/R source snapshots, configuration, Slurm entrypoints,
  and handoffs. Code snapshots include files not yet tracked by Git.
- `generation/`: submitted final-render job, log, and Git status at freezing.
- `frozen_manifest.json`: sizes and SHA-256 hashes for every frozen artifact.

The working render package remains at
`final_figures/mge_activity_IEG_dots_20261101/`. The frozen package is the
reference for this finalized version; preserve it when rerunning analysis.

## Reproduce the finalized figure without rerunning analysis

From this checkout:

```bash
sbatch --wait paper3_pcdh19/slurm/reproduce_final_activity_IEG_figure.sbatch
```

The wrapper verifies every frozen hash, copies inputs into a fresh timestamped
NFS directory and executes the frozen renderer. It never writes into the
frozen package. It compares PNG and all three IEG data tables byte for byte
against the final references. PDF/SVG metadata can differ. Its
`REPRODUCTION_VALIDATION.json` records the environment, font and comparisons.
The frozen wrapper is also in `code/paper3_pcdh19/scripts/`.

## Definitions to retain when rebuilding the individual analyses

A/B reuse the saved atlas embedding and regional program scores. C/E show
Pcdh19 expression; D/F retain existing gene candidates and edgeR effects.
G/H show the same five curated expression-based ligand–receptor candidates,
without depth matching; these are not direct functional signaling measurements.
Earlier detailed source-chain and selection provenance remains in
`APPROVED_MGE_FIGURE_REPRODUCTION.md`, included in the code snapshot.

I uses existing inhibitory clusters C0004/5/6/20, Gad1 or Gad2 detection,
>=1 Snap25/Syt1/Rbfox3 detection, and zero Mki67/Top2a UMI. HET Female WT/KO
are A-positive versus A-negative/B+C>=1 probe states within each of three
samples. Male WT/KO are known genotypes, three separate samples each, with
no Pcdh19 expression filter. No direct Female-versus-Male inference is made.

I expression circles use equal-sample mean log-normalized expression, divided
by each gene's maximum across the four groups (0–1). Size and purple color
encode relative expression, not detection. Signed changes are KO minus WT on
that same gene-relative scale, separately in HET Female and Male: red higher
KO, blue lower KO. Size encodes absolute change. Rows rank by the average of
the two scaled changes; all-zero genes last. This ranking is descriptive,
not statistical significance. Scaling can emphasize extremely sparse genes.

All nine genes are retained: Npas4, Egr1, Fos, Egr3, Junb, Nr4a1, Arc, Fosb,
Egr2 (final display order). Preliminary takeaway: Npas4 is higher in both
comparisons, while the broader IEG response is mixed. HET Female counts are
particularly sparse and gate-sensitive. The activity-dependent title names
the gene class; there was no stimulation experiment or direct cleavage assay.

Preserve singular sex wording everywhere. C labels must be exactly
“WT Female”, “HET Female WT cells”, “HET Female KO cells”. E labels are
“WT Male”, “KO Male”. C has no zero-signal explanatory text. Panel I keeps
both expression and change views, with separated group labels and color scales.

When redoing upstream analyses, use a new runtime output directory and compare
results against this frozen version. Reproduction uses fixed inputs; rebuilding
raw processing, probe classification, DE, and interactions is a distinct task.

## Verified final reproduction

Slurm job 65185044 reproduced the final PNG and all three IEG tables byte
for byte from the frozen inputs and code. Environment: Python 3.6.8, NumPy
1.19.5, SciPy 1.5.4, matplotlib 3.3.4, bundled DejaVu Sans. Verification
report and log are archived in `generation/`. Reproduction directory:
`/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/final_IEG_reproduction_20261101_194555`.
