# PCDH19 HiCAT technical pilot handoff

**Status: IN_REVIEW.** Computation completed on 2026-09-08; no approval inferred.

Review package: `hicat_pilot_20260908_202159_30835fdd`; SLURM job `60559849`.

## What the pilot answers

The Python adapter executes pooled expression clustering with recorded recursive
splits and DE-based merges, preserves the fitted models and cell identities, and
exports a reviewable AnnData/report package. It uses 100 randomly selected cells
from each of 12 samples: 1,200 cells and all 19,071 genes. Approved Step 02 supplies
raw counts; Step 06 supplies existing display coordinates only.

The baseline produced **4** final clusters, the second seed **3**, and the higher
separation-score threshold **2**. These counts show sensitivity even on the same
selected cells. **No full-dataset cluster count is established or locked.**

| Candidate | Seed | DE-score threshold | Root graph | After root merge | Recursive leaves | Final K |
|---|---:|---:|---:|---:|---:|---:|
| Baseline | 20260908 | 150 | 26 | 4 | 4 | 4 |
| Seed repeat | 20260909 | 150 | 27 | 4 | 4 | 3 |
| Higher separation requirement | 20260908 | 300 | 26 | 2 | 2 | 2 |

All final pairs passed the declared score/minimum-gene audit in their respective
candidate, but this does not establish stability. The baseline versus repeat-seed
adjusted Rand index was 0.758. Baseline groups C0002 and C0004 had best Jaccard
matches of 0.415 and 0.485: their distinction was less stable than the other groups.
The baseline cluster sizes were 745, 200, 49 and 206. Equal sample subsampling and
the small pilot limit conclusions about rare populations and population frequency.

## Open these assets

- [Three-page review PDF](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/00_method_and_technical_pilot/hicat_pilot_20260908_202159_30835fdd/outputs/figures/hicat_pilot_review.pdf): count progression,
  fixed UMAP, sizes, sample contributions, evidence genes, pairwise separation,
  recursive split tree, seed overlap and expression-similarity dendrogram.
- [Reusable pilot AnnData](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/00_method_and_technical_pilot/hicat_pilot_20260908_202159_30835fdd/outputs/pcdh19_hicat_pilot.h5ad): raw integer `.X`, exact
  normalized `.layers['log1p_cpm']`, all three final partitions and their recursive
  predecessors, root PCA/graph, and explicitly named existing display coordinates.
- [Exact `.uns` inventory](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/00_method_and_technical_pilot/hicat_pilot_20260908_202159_30835fdd/outputs/anndata_uns_inventory.json): structured metadata,
  provenance, resolved settings, candidate summaries and interpretation of slots.
  The only top-level key is `hicat`; configuration and candidate summaries are
  nested JSON strings within that dictionary for H5AD compatibility.
- [Count table](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/00_method_and_technical_pilot/hicat_pilot_20260908_202159_30835fdd/outputs/cluster_count_summary.tsv),
  [cell assignments](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/00_method_and_technical_pilot/hicat_pilot_20260908_202159_30835fdd/outputs/cell_assignments.tsv.gz), and
  [stability table](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/00_method_and_technical_pilot/hicat_pilot_20260908_202159_30835fdd/outputs/cluster_stability.tsv).
- [Complete output directory](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/00_method_and_technical_pilot/hicat_pilot_20260908_202159_30835fdd/outputs): each candidate's per-node PCA loadings,
  feature masks, centering, projections, graphs/indexes, memberships, DE scores
  and evidence genes, accepted merges, final merge model and endpoint audit.
- [Run directory](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/00_method_and_technical_pilot/hicat_pilot_20260908_202159_30835fdd): frozen Python adapter/upstream source, exact
  config/environment, selected cell IDs and input hashes, scheduler logs and
  printed/file-backed function progress. Keep this directory with the H5AD.

## What was not saved or inferred

No full-data assignments, final full-data K, biological annotations, bootstrap
consensus, new UMAP, correction/regression model, universal recursive PCA, or
new-cell classifier was produced. `.raw` is absent because raw counts are already
in `.X`. Per-node fitted models live outside the H5AD. The upstream DE-summary
API returns pair scores and evidence gene lists; full per-gene test statistics
were not saved. Gene names in the figure explain separation; they are not cell
type labels. This pilot does not validate the earlier Step 06 sample-structure
headline or change any preprocessing approval status.

## Validation

The final scheduler job completed successfully in 2 minutes 20 seconds, with
approximately 8.64 GiB peak batch memory (including input preparation).
All **19 saved-data checks passed**, including exact raw/normalized expression,
cell/gene identity, assignments, and structured `.uns` inventory equality after
reopening the H5AD. All **360 output manifest entries** matched their saved sizes
and SHA-256 hashes on independent verification. Five synthetic contract tests
passed, including the real upstream merge path and AnnData metadata export.
All three final PNGs are byte-identical to the visually reviewed previews.

The [HiCAT review ledger](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/APPROVAL_LEDGER.tsv)
records the terminal attempts and leaves successful packages IN_REVIEW. Use the
replacement run linked above for review; the earlier completed package retains
its original inventory-export limitation.

## Reproducibility and review boundary

The [method and tuning guide](HICAT_CLUSTER_DISCOVERY_PROTOCOL.md) documents
normalization units, every adjustable control, declared adapter differences,
asset schemas, plot panels, progress logging and a function-by-function index.
The [pilot config](config/hicat_pilot.json) is the central tuning surface.
Higher score requirements make evidence for retaining a split more demanding;
settings must be evaluated for stability rather than tuned toward a desired K.

Earlier attempts remain intact: `201306` failed on a gene-ID column assumption;
`201412` failed when upstream merging required list memberships. `201734`
completed clustering and its original 17 checks, but review found its external
`.uns` JSON inventory was a string representation. The replacement run preserves
structured mappings and checks the inventory against the reopened H5AD. None of
these historical packages was rewritten or silently promoted to approval.

The next checkpoint is review of this method and pilot, followed by full-data
resource/execution qualification and a stronger stability design. The dense
pilot execution path is deliberately capped at 1,200 cells; it must not simply
be pointed at all 446,349 cells. No full-data job has been submitted.
