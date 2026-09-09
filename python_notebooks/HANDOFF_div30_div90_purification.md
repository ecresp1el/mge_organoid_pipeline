# DIV90 purification audit and separate DIV30 antecedent discovery

Resumed 2026-09-09 from the completed DIV90 Phase 2 PI figure package. This
continuation changes the optimization question; it preserves the earlier
scores, annotations, coordinates, gate results, reports and figure package.

## Authoritative starting interpretation

FGFR2/PTPRS retains 88.7% of PV-biased and 89.3% of dual-high cells, but changes
SST-biased composition only from 22.5% to 21.3%. This is weak purification.
Independent statistical significance is not established. It is a diagnostic
benchmark, not the preferred purification solution or an antibody recommendation.

## Separate output locations

All paths below are under
`/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/`:

- `div90_purification_audit_v1/`: reconstruction and reranking of all 2,885,531
  existing gate evaluations; additional signed and hierarchical gate searches.
- `div30_pv_antecedent_v1/`: independent DIV30 scoring, reference mapping,
  surface discovery and developmental/confound diagnostics.
- `div30_div90_purification_comparison_v1/`: original request copies, source
  preservation inventory and the final linked interpretation.

The two pasted instruction files are copied verbatim to the comparison
directory's `provenance/`. The original phase-one and phase-two directories
(including `pi_figure_package_v1/`) are read-only inputs to this work.

## Analysis contract

DIV90 searches target strict PV-biased cells, PV-biased plus an operational
high-confidence PV-supported subset of dual-high cells, and independently
rank-scaled high-PV/low-SST continuous targets. Dual-high quality filters are
not independent biological validation. Their thresholds and sensitivity
results must remain explicit. Every retained population includes unresolved
cells in its purity denominator.

Rank primarily by retained purity or PV:SST ratio at stated recovery floors,
not F1. Report total yield, all four original state fractions, conservative
dual-high fraction/recovery, target enrichment and SST depletion. Present
unconstrained maxima alongside useful recovery thresholds, because a tiny
pure fraction is not evidence of a practical sort. Three-marker searches are
bounded unless explicitly proven exhaustive; observed maxima are not a
theoretical upper bound on the entire surfaceome.

DIV30 uses every available cached cell without an ERBB4 or LHX6 entry
requirement. The frozen DIV90 PV program supplies the external biological
reference. Direct scoring, rank scoring, cell-level mapping and a sparse
supervised model must show their agreement or disagreement before a consensus
is interpreted. A top quantile exists even when absolute program evidence is
weak; examine gene support, tie handling, distance to reference, developmental
state and maturation dependence. Define a separate DIV30 surface feature set
from the annotation catalogue using DIV30 expression, excluding score genes.

Compare ages only after their searches are solved. A transcriptional
antecedent or developmental resemblance is a cross-sectional hypothesis;
neither nearest neighbors nor score agreement demonstrates lineage to later
PVALB protein. The later protein endpoint is supplied by the user and is not
remeasured by these computations.

## Preservation and runtime

Run code with `/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/python`,
`PYTHONDONTWRITEBYTECODE=1`, bounded BLAS threads and caches under the new
runtime directories. Do not rerun the full biological pipeline, embeddings,
identity recovery or the completed PI renderer.

`scripts/build_div30_div90_purification_comparison.py --stage initialize`
records the original request and frozen source inventory without overwriting
the before-manifest. `--stage verify` checks that the old directories remain
unchanged. Scientific artifacts are checked with SHA256; installed dependency
and cache files are tracked by file size and modification time.

## Initial independent verification

Before any new logic search, a separate scan of the two old HDF5 tables found
that their original retained fractions reach 36.94% strict PV-biased purity
only at 3.81% PV recovery. At a 50% recovery minimum, their maximum is 28.56%
purity (573 PV-biased cells among 2,006 retained). At a 70% minimum it is
27.84% (804 among 2,888). These checks do not include new positive or mixed-sign
gates. Exact source row IDs are saved in
`provenance/independent_legacy_purity_optima.tsv` in the comparison directory.

## Completion status

Both analyses are complete. Open
[the review index](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_div90_purification_comparison_v1/index.html),
[DIV90 report](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_purification_audit_v1/REPORT.md),
[DIV30 report](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_pv_antecedent_v1/REPORT.md),
and [linked interpretation](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_div90_purification_comparison_v1/REPORT.md).

DIV90 now covers 11,143,135 evaluated boundaries, including every legacy row.
At >=10% strict-PV recovery, `NTM <= 0 AND PCDH9 > 1.9864086000000001 AND
PTPRA <= 1.9811362250000002` retains 260 cells: 44.62% strict PV-biased,
8.08% SST-biased, and 10.78% PV recovery. At >=50% recovery, the maximum
observed strict-PV purity is 31.33%; at >=70%, it is 29.07%. Approximately
2x enrichment is reached at 9.94% PV recovery. The 65–75% purity solutions
preserve fewer than 2% of PV-biased cells. Three markers improve the >=10%
optimum by 2.71 percentage points relative to the best pair. These results
require both purity-oriented ranking and expanded signed/hierarchical logic.

DIV30 includes 90,631 cached cells from **six** samples. The historical
`9583-MW-N` / `9853-MW-N` sample alias is documented and reconciled without
changing source IDs. Metadata-only samples 7–9 (79B_old) are absent from the
frozen H5AD; they were not silently treated as analyzed. The candidate
top20 consensus population is 86.4% operational postmitotic, and requiring
ERBB4 RNA detection would lose 68.0% of it.

The DIV30 rule `GRIA2 > 1.1180085897445677 AND OPCML > 0.47731003165245056
AND NOTCH1 <= 0.88839703798294067` yields 58.36% top20 target composition,
2.92x enrichment and 50.81% recovery. Held-sample gate reselection at the 50%
training floor yields 53.49% composition but 39.73% actual recovery; it uses
different gates across folds. End-to-end gene-drop sensitivities and measured
developmental confound checks support a candidate transcriptional state, not
a confirmed future PV fate. FGFR2/PTPRS gives 1.56x enrichment of the separate
DIV30 target, primarily as a progenitor/maturation diagnostic.

The selected DIV30 and DIV90 >=50%-recovery rules each reduce target
composition when applied at their exact thresholds to the other age. This
supports investigating different stage-specific strategies, without proving
the developmental transitions. The schematic marks both transitions
TO BE VALIDATED and the later protein observation as user-reported.

Independent validations cover all legacy DIV90 counts, sampled and selected
DIV90 masks, all DIV30 archived count bounds, selected gates, random single
and pair masks, all 9,194 DIV30 triplets, held-sample masks, fixed direct-score
scaling, and cross-age applications. Preservation checks passed for 2,066
preexisting files, including SHA256 equality for 330 scientific artifacts.

New code consists of `div90_purification_{audit,validate}.py`,
`div30_antecedent_{analysis,surface}.py`,
`validate_div30_antecedent_surface.py`, the two `render_*` scripts, and the
comparison initialization/packaging scripts in `scripts/`. Runtime code
snapshots, source hashes and machine-readable validation remain with each
output. The compact `DIV30_DIV90_Purification_Review.zip` includes reports,
figures and summary tables; complete HDF5 archives and cell-level caches
remain on Turbo. Do not rerun discovery over these completed versions. Use a
new version for changed scientific definitions or additional data.
