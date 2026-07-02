# Handoff: Imaris CAREamics N2V Scale-Up

Date: 2026-07-01 15:40 EDT
Updated: 2026-07-02 corrected baseline complete; 8-bit display exports complete; QC review next

This handoff covers the non-scRNA imaging branch for two-channel Imaris `.ims`
confocal stacks. The raw `.ims` files are read only and are never modified.

## Scope

There are two complete F0-F3 movie sets, eight `.ims` files total:

| Set | Samples |
| --- | --- |
| `cl32_bive4_pv_reporter_40x` | `F0`, `F1`, `F2`, `F3` |
| `cl32_bive4_pv_reporter_40x_realbive4` | `F0`, `F1`, `F2`, `F3` |

The pilot/full run already completed for:

`cl32_bive4_pv_reporter_40x_realbive4_F1`

The other seven full denoising jobs were submitted as separate Slurm jobs and
all completed with Slurm exit code `0:0`. Post-run file validation found that
the `cl32_bive4_pv_reporter_40x_F0-F3` denoised outputs are invalid because
their denoised OME-TIFFs contain `NaN` values. The original
`cl32_bive4_pv_reporter_40x_realbive4_F0-F3` denoised outputs are finite image
data, but F0/F2/F3 showed visible restoration artifacts after visual QC. Those
finite-but-bad outputs should not be used as final clean data.

As of 2026-07-02 15:00 EDT, all eight corrected baseline outputs completed:

- `cl32_bive4_pv_reporter_40x_F0-F3`: trim empty Z, disable batch norm, predict
  from best checkpoint.
- `cl32_bive4_pv_reporter_40x_realbive4_F0-F3`: repredict from existing models
  using best checkpoint.

The next required step is QC review/comparison before any additional parameter
sweep or final movie generation.

As of 2026-07-02 15:28 EDT, all eight corrected baseline outputs also have
8-bit display/QC exports under `full_run/display_8bit/`. These are derived
copies only. The quantitative denoised float32 OME-TIFFs remain unchanged.

As of 2026-07-02, the BC43/realbive4 2x2 montage has a reproducible Fiji
Grid/Collection Stitching tile-configuration prep step. The parser reads the
TeraStitcher XML `Stack` entries and preserves the XML order rather than
guessing from F0-F3 filenames.

Generated layout file:

`/nfs/turbo/umms-parent/andor_micropscope_data_dump/exp17_pv_reporter_with_biver3and4/20x/TileConfiguration.txt`

Repo helpers:

```text
imaging/imaris_careamics/terastitcher_to_fiji_tile_config.py
imaging/imaris_careamics/fiji_grid_collection_stitching.ijm
slurm/prepare_fiji_stitching_tile_config.sbatch
slurm/run_fiji_grid_stitching_realbive4.sbatch
```

Submitted prep smoke job:

```text
52781973  fiji-tilecfg  COMPLETED, exit 0:0, elapsed 00:00:02
```

The layout generated from
`cl32_bive4_pv_reporter_40x_realbive4.xml` is:

```text
# Define the number of dimensions we are working on
dim = 3

# Define the image coordinates
cl32_bive4_pv_reporter_40x_realbive4_F3.ims; ; (0.0, 0.0, 0.0)
cl32_bive4_pv_reporter_40x_realbive4_F0.ims; ; (918.0, 0.0, 0.0)
cl32_bive4_pv_reporter_40x_realbive4_F2.ims; ; (0.0, 897.0, 0.0)
cl32_bive4_pv_reporter_40x_realbive4_F1.ims; ; (918.0, 897.0, 0.0)
```

Run Fiji from a Great Lakes desktop session, not from a headless Slurm batch:

```bash
module load fiji/1.5.4
/sw/pkgs/med/fiji/1.5.4/ImageJ-linux64 &
```

Fiji settings:

- `Plugins -> Stitching -> Grid/Collection Stitching`
- Type: `Positions from file`
- Layout file: generated `TileConfiguration.txt`
- Fusion method: `Linear Blending`
- Compute overlap: checked for the first test
- Regression threshold: `0.30`
- Max/avg displacement threshold: `2.50`
- Absolute displacement threshold: `3.50`
- Image output: `Fuse and display`

Important constraints: do not use the FusionStitcher stitched output, do not
manually guess F0-F3 order, trust the XML coordinates, and do not split channels
before stitching.

Actual batch stitch status:

```text
52782611  fiji-stitch-realbive4  CANCELLED before start; 180G request was too conservative for scheduling
52782628  fiji-stitch-realbive4  PENDING at last check; 96G, 8 CPU, 12h, Fiji heap 72g
```

The queued batch run uses `Image_output=[Write to disk]` instead of the
interactive `Fuse and display` setting so the Slurm run leaves durable output
files. It still uses Fiji Grid/Collection Stitching with `Positions from file`,
`Defined by TileConfiguration`, linear blending, compute-overlap enabled, and
the requested displacement thresholds.

Expected output directory:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/fiji_stitching/cl32_bive4_pv_reporter_40x_realbive4_xml_coords_compute_overlap/
```

Log paths:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/fiji-stitch-realbive4-52782628.out
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/fiji-stitch-realbive4-52782628.err
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/fiji_stitching/cl32_bive4_pv_reporter_40x_realbive4_xml_coords_compute_overlap/logs/run_fiji_grid_stitching_52782628.log
```

Status check:

```bash
squeue -j 52782628 -o '%.18i %.9P %.28j %.8u %.2t %.12M %.6D %.10m %R'
sacct -j 52782628 --format=JobID,JobName%28,State,ExitCode,Elapsed,MaxRSS
```

## Submitted Jobs

Submitted from:

`/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline`

| Job ID | Sample | Slurm status | Output validation | Elapsed | Output directory |
| --- | --- | --- | --- | --- | --- |
| `52693455` | `cl32_bive4_pv_reporter_40x_F0` | `COMPLETED`, exit `0:0` | Invalid: denoised outputs are `NaN` | `00:31:10` | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_F0` |
| `52693456` | `cl32_bive4_pv_reporter_40x_F1` | `COMPLETED`, exit `0:0` | Invalid: denoised outputs are `NaN` | `00:31:06` | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_F1` |
| `52693458` | `cl32_bive4_pv_reporter_40x_F2` | `COMPLETED`, exit `0:0` | Invalid: denoised outputs are `NaN` | `00:31:22` | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_F2` |
| `52693459` | `cl32_bive4_pv_reporter_40x_F3` | `COMPLETED`, exit `0:0` | Invalid: denoised outputs are `NaN` | `00:33:17` | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_F3` |
| `52693460` | `cl32_bive4_pv_reporter_40x_realbive4_F0` | `COMPLETED`, exit `0:0` | Valid finite denoised data | `00:51:39` | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_realbive4_F0` |
| `52693461` | `cl32_bive4_pv_reporter_40x_realbive4_F2` | `COMPLETED`, exit `0:0` | Valid finite denoised data | `00:51:38` | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_realbive4_F2` |
| `52693462` | `cl32_bive4_pv_reporter_40x_realbive4_F3` | `COMPLETED`, exit `0:0` | Valid finite denoised data | `00:51:31` | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_realbive4_F3` |

The logs contain Python multiprocessing/NFS temporary-directory cleanup warnings
(`OSError: [Errno 16] Device or resource busy`) after prediction, as seen in the
pilot run. The Slurm jobs still completed with exit code `0:0` and wrote the
expected outputs.

The first `cl32_bive4_pv_reporter_40x_F0-F3` set also contains 23 trailing
all-zero Z planes in the raw channel stacks. In the logs, N2V training loss turns
to `nan` during epoch 0, and the final denoised OME-TIFFs are all `NaN`. These
outputs should be treated as failed denoising runs despite Slurm success.

Post-run validation artifacts:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/scaleup_postrun_validation.tsv
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/scaleup_denoised_midz_contact_sheet.png
```

## Trim-Z Fix Rerun

On 2026-07-02, the failed `cl32_bive4_pv_reporter_40x_F0-F3` jobs were
resubmitted to separate output directories so the failed outputs remain intact.

Fix settings:

```text
TRIM_EMPTY_Z=true
DISABLE_BATCH_NORM=true
NORMALIZATION=mean_std
PREDICT_CHECKPOINT=best
TILE_OVERLAP_Z=8
TILE_OVERLAP_YX=32
```

The fix excludes leading/trailing all-zero Z planes from training/prediction,
then pads the denoised output back to the original Z size with zeros. For these
files, the nonzero Z range is `0-360`, with trailing zero planes `361-383`.

| Job ID | Sample | Status | Output directory |
| --- | --- | --- | --- |
| `52759673` | `cl32_bive4_pv_reporter_40x_F0` | Canceled after starting before best-checkpoint patch | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_F0_trimz_fix` |
| `52759674` | `cl32_bive4_pv_reporter_40x_F1` | Failed fast: checkpoint monitor config bug | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_F1_trimz_fix` |
| `52759675` | `cl32_bive4_pv_reporter_40x_F2` | Failed fast: checkpoint monitor config bug | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_F2_trimz_fix` |
| `52759676` | `cl32_bive4_pv_reporter_40x_F3` | Failed fast: checkpoint monitor config bug | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_F3_trimz_fix` |
| `52760969` | `cl32_bive4_pv_reporter_40x_F0` | Completed, exit `0:0`, elapsed `00:30:05` | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_F0_trimz_bestckpt_fix` |
| `52761073` | `cl32_bive4_pv_reporter_40x_F1` | Completed, exit `0:0`, elapsed `00:30:09` | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_F1_trimz_bestckpt_fix` |
| `52761074` | `cl32_bive4_pv_reporter_40x_F2` | Completed, exit `0:0`, elapsed `00:30:25` | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_F2_trimz_bestckpt_fix` |
| `52761075` | `cl32_bive4_pv_reporter_40x_F3` | Completed, exit `0:0`, elapsed `00:30:02` | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_F3_trimz_bestckpt_fix` |

## RealBiVe4 Visual QC Issue

After visual review on 2026-07-02, the original finite `realbive4` outputs were
not all acceptable. F1 looked usable, but F0/F2/F3 showed visible artifacts,
especially in green-channel denoising and max projections.

Diagnosis:

- The jobs trained one N2V model per field of view and predicted from the final
  epoch.
- Several field/channel models overtrained: their best validation loss occurred
  much earlier than epoch 49.
- Examples from the original runs:
  - `realbive4_F2` green: best saved checkpoint `0.3407` at epoch 9, final
    checkpoint `1.1929` at epoch 49.
  - `realbive4_F3` red: best saved checkpoint `0.4559` at epoch 9, final
    checkpoint `0.9536` at epoch 49.
  - `realbive4_F0` red: best saved checkpoint `0.2254` at epoch 29, final
    checkpoint `0.5287` at epoch 49.
- Because N2V has no ground-truth target, validation loss is a self-supervised
  blind-spot loss. A finite output can still be visually poor if the selected
  checkpoint is bad.

Pipeline patch:

- `denoise_ims_careamics.py` now defaults to `--predict-checkpoint best`.
- New training saves the top validation checkpoints instead of relying on the
  final epoch.
- `--model-source-dir` allows prediction-only reruns from existing trained
  models.
- `prediction_checkpoints.json` records the checkpoint used for each channel.

Best-checkpoint reprediction jobs were submitted without retraining and all
completed:

| Job ID | Sample | Status | Output directory |
| --- | --- | --- | --- |
| `52760972` | `cl32_bive4_pv_reporter_40x_realbive4_F0` | Completed, exit `0:0`, elapsed `00:40:28` | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_realbive4_F0_bestckpt_repredict` |
| `52760975` | `cl32_bive4_pv_reporter_40x_realbive4_F1` | Completed, exit `0:0`, elapsed `00:40:28` | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_realbive4_F1_bestckpt_repredict` |
| `52760976` | `cl32_bive4_pv_reporter_40x_realbive4_F2` | Completed, exit `0:0`, elapsed `00:41:15` | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_realbive4_F2_bestckpt_repredict` |
| `52760977` | `cl32_bive4_pv_reporter_40x_realbive4_F3` | Completed, exit `0:0`, elapsed `00:40:26` | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_realbive4_F3_bestckpt_repredict` |

Status command for corrected jobs:

```bash
squeue -j 52760969,52760972,52760975,52760976,52760977,52761073,52761074,52761075 \
  -o '%.18i %.9P %.35j %.8u %.2t %.12M %.6D %.12b %.10m %R'
```

At the time of this update, this queue is empty for the corrected job IDs.

Corrected output directories to review:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_F0_trimz_bestckpt_fix/full_run/
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_F1_trimz_bestckpt_fix/full_run/
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_F2_trimz_bestckpt_fix/full_run/
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_F3_trimz_bestckpt_fix/full_run/
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_realbive4_F0_bestckpt_repredict/full_run/
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_realbive4_F1_bestckpt_repredict/full_run/
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_realbive4_F2_bestckpt_repredict/full_run/
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_realbive4_F3_bestckpt_repredict/full_run/
```

## Parameter Optimization Plan

The next decision point is not "did Slurm finish?" but "which denoising
parameters preserve real biology while removing noise without creating shading,
grid, edge, or over-smoothing artifacts?" N2V is self-supervised, so there is no
ground-truth target. Final parameter agreement must combine model diagnostics
with visual QC.

### Immediate QC review

Before running any parameter sweep, review the corrected baseline outputs that
already exist. Open these first for each corrected sample:

| File | How to use it |
| --- | --- |
| `qc/raw_vs_denoised_example_slices.png` | Fast first look for over-smoothing, artifacts, and preserved structures |
| `qc/raw_vs_denoised_max_projections.png` | Check whether structures are preserved across Z |
| `merged_green_red_preview.ome.tif` | Fiji-friendly preview-scaled OME-TIFF |
| `green_denoised.ome.tif` and `red_denoised.ome.tif` | Full float32 denoised outputs; autoscale display in Fiji before judging |
| `max_projection_preview.mp4` | Quick preview only; not a final z-scroll movie |
| `prediction_checkpoints.json` | Confirms which checkpoint was used |
| `qc/denoised_finite_validation.txt` | Confirms no NaN/Inf values |

Recommended human review order:

1. Compare `realbive4_F1_bestckpt_repredict` with the earlier F1 output that
   looked acceptable.
2. Review `realbive4_F0`, `F2`, and `F3` best-checkpoint repredict outputs to
   see whether the horrible artifacts improved.
3. Review first-set `F0-F3_trimz_bestckpt_fix` to confirm the NaN/black-output
   failure is resolved.
4. Record pass/fail notes per field of view:
   - acceptable as baseline
   - improved but still artifacted
   - unusable / needs parameter sweep
5. Only after this review, choose whether to run the small parameter sweep below.

Do not judge final quality from the old invalid folders or from raw float32
OME-TIFFs opened in Fiji without display autoscaling. Prefer the QC PNGs and
`merged_green_red_preview.ome.tif` for the first pass.

For easier Fiji review, the corrected folders now also contain 8-bit display
exports in `full_run/display_8bit/`. These are useful for visual QC and sharing,
but they should not replace the quantitative denoised OME-TIFFs for analysis.

### Candidate knobs to compare

Start with a small controlled sweep before any broad rerun:

| Knob | Baseline | Candidate values | Why it matters |
| --- | --- | --- | --- |
| Prediction checkpoint | `best` | `best`, possibly selected earlier checkpoint if visual QC prefers it | Avoids overtrained final-epoch models |
| Epochs/checkpoint timing | `50` epochs, best checkpoint | `10`, `20`, `35`, `50` or existing saved checkpoints | N2V can look best before final epoch |
| Patch size | `16 x 64 x 64` | `16 x 64 x 64`, `24 x 96 x 96`, `32 x 128 x 128` if memory allows | Larger context may preserve 3D structures better |
| Prediction overlap | `8 x 32 x 32` | `8 x 32 x 32`, `12 x 48 x 48`, `16 x 64 x 64` | More overlap can reduce tile/grid artifacts |
| Batch norm | disabled for fragile runs | enabled vs disabled only if batch size supports it | Batch norm with batch size 1 can be unstable |
| Normalization | `mean_std` | `mean_std`, `min_max`, `quantile` | Weak/sparse green channel may be sensitive |
| Training scope | per-FOV model | per-FOV vs shared model across F0-F3 for each channel | Shared models may reduce field-specific failures |

### Recommended first sweep

Use one visually acceptable FOV and one visually bad FOV, not all eight files.
Recommended initial pair:

- acceptable/reference: `cl32_bive4_pv_reporter_40x_realbive4_F1`
- problematic: `cl32_bive4_pv_reporter_40x_realbive4_F2` or `F3`

Run small representative crops first. Do not use full-stack reruns to tune
parameters unless a candidate already passes crop-level QC.

Initial candidate set:

```text
A: best checkpoint, patch 16x64x64, overlap 8x32x32
B: best checkpoint, patch 24x96x96, overlap 12x48x48
C: best checkpoint, patch 32x128x128, overlap 16x64x64
D: shared F0-F3 channel model, patch 24x96x96, overlap 12x48x48
```

### Required QC assets for agreement

For every parameter candidate, generate the same fixed-scale review assets:

| Asset | Required content |
| --- | --- |
| Slice contact sheet | Raw, denoised, residual/difference, merged for selected Z slices |
| Max-projection contact sheet | Raw vs denoised vs residual for green and magenta/red |
| Z-scroll movie | 1x3 panel: `BiVe3-dTom` magenta, `PV-mNG` green, merged |
| Display metadata | Percentile limits fixed across the full stack/movie, not per slice |
| Quantitative summary | Finite check, intensity percentiles, residual percentiles, saturation report |
| Model summary | Checkpoint used, val loss, epochs, patch size, overlap, normalization, batch norm |

The residual/difference image is mandatory because it shows whether the model is
removing mostly noise or erasing real structures / inventing low-frequency
background.

### Agreement criteria

Before declaring final parameters, the chosen setting should satisfy all of the
following:

- No non-finite denoised pixels.
- No obvious tile/grid artifact in slices, max projections, or movies.
- No recurring autoscale changes in QC movies; display limits must be fixed.
- Green and red/magenta biological structures remain visible and spatially
  plausible.
- Background noise is reduced without flattening the whole green channel into a
  smooth field.
- Bright puncta/processes are not erased or strongly hallucinated.
- Residual images look noise-like, not structure-like.
- The same setting performs acceptably on both the good/reference FOV and a
  bad/problem FOV.

### Decision workflow

1. Let the currently submitted best-checkpoint and trim-Z jobs finish.
2. Build a QC comparison panel for the completed baseline outputs.
3. Select one good and one bad FOV for tuning.
4. Run the small candidate sweep on representative crops.
5. Generate fixed-scale QC panels and movies for all candidates.
6. Review side by side and choose the smallest set of acceptable parameters.
7. Run full-stack denoising only for the chosen parameter set.
8. Generate final publication/QC movies from the selected outputs.
9. Record the final agreed parameter set in this handoff and in the README.

Current status command:

```bash
squeue -j 52693455,52693456,52693458,52693459,52693460,52693461,52693462 \
  -o '%.18i %.9P %.35j %.8u %.2t %.12M %.6D %.12b %.10m %R'
```

Completion/resource summary:

```bash
sacct -j 52693455,52693456,52693458,52693459,52693460,52693461,52693462 \
  --format=JobID,JobName%30,Partition,State,Elapsed,ExitCode,NodeList,MaxRSS,ReqMem
```

## Expected Per-Sample Outputs

Each sample writes into:

`<OUTPUT_DIR>/full_run/`

Expected full-run files:

| Path | Meaning |
| --- | --- |
| `green_raw.ome.tif` | Raw channel 0 saved as OME-TIFF, unchanged except format conversion |
| `red_raw.ome.tif` | Raw channel 1 saved as OME-TIFF, unchanged except format conversion |
| `green_denoised.ome.tif` | CAREamics N2V denoised green channel |
| `red_denoised.ome.tif` | CAREamics N2V denoised red channel |
| `merged_green_red_preview.ome.tif` | Preview-scaled two-channel OME-TIFF |
| `max_projection_preview.mp4` | Quick max-projection preview movie |
| `models/green_n2v/` | Green-channel N2V model/checkpoints |
| `models/red_n2v/` | Red-channel N2V model/checkpoints |
| `prediction_checkpoints.json` | Checkpoint selected for green/red prediction |
| `qc/raw_vs_denoised_example_slices.png` | Slice-level raw vs denoised QC |
| `qc/raw_vs_denoised_max_projections.png` | Max-projection raw vs denoised QC |
| `qc/intensity_histograms.png` | Raw/denoised intensity distributions |
| `qc/clipping_saturation_warnings.txt` | Saturation/non-finite warning report |
| `run_metadata.json` | Input, internal dimensions, args, reader metadata |

Optional post-N2V 8-bit display exports live under:

`<OUTPUT_DIR>/full_run/display_8bit/`

| Path | Meaning |
| --- | --- |
| `green_denoised_8bit_display.ome.tif` | Green denoised channel converted to uint8 for display/QC |
| `red_denoised_8bit_display.ome.tif` | Red denoised channel converted to uint8 for display/QC |
| `merged_green_magenta_8bit_display.ome.tif` | RGB display stack with red channel shown as magenta and green as green |
| `denoised_8bit_green_max_projection.png` | Green 8-bit max-projection quicklook |
| `denoised_8bit_red_as_magenta_max_projection.png` | Red channel shown as magenta max-projection quicklook |
| `denoised_8bit_merged_green_magenta_max_projection.png` | Merged green/magenta max-projection quicklook |
| `display_scaling_8bit.json` | Reproducibility sidecar with input paths, dtype, axes, display limits, and clipping counts |

The 8-bit conversion is performed by:

`imaging/imaris_careamics/convert_denoised_to_8bit.py`

Default conversion command:

```bash
python imaging/imaris_careamics/convert_denoised_to_8bit.py \
  --input-dir /path/to/full_run \
  --overwrite
```

Default scaling is fixed per channel across a deterministic whole-stack sample:

```text
lower percentile: 0.5
upper percentile: 99.8
percentile z step: 2
percentile yx step: 2
```

The sampling avoids memory failures on the large `realbive4` stacks while still
using a fixed stack-wide display transform. Use `--percentile-z-step 1
--percentile-yx-step 1` only when exact all-pixel limits are required and enough
memory is available.

Great Lakes/Slurm wrapper:

```bash
FULL_RUN_DIR=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/<sample>/full_run \
sbatch slurm/convert_denoised_to_8bit.sbatch
```

Corrected baseline directories already converted to 8-bit display outputs:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_F0_trimz_bestckpt_fix/full_run/display_8bit/
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_F1_trimz_bestckpt_fix/full_run/display_8bit/
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_F2_trimz_bestckpt_fix/full_run/display_8bit/
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_F3_trimz_bestckpt_fix/full_run/display_8bit/
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_realbive4_F0_bestckpt_repredict/full_run/display_8bit/
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_realbive4_F1_bestckpt_repredict/full_run/display_8bit/
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_realbive4_F2_bestckpt_repredict/full_run/display_8bit/
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_realbive4_F3_bestckpt_repredict/full_run/display_8bit/
```

Expected logs:

| Path | Meaning |
| --- | --- |
| `<OUTPUT_DIR>/logs/run_full_denoising_<jobid>.log` | Main tee'd analysis log |
| `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/imaris-n2v-full-<jobid>.out` | Slurm stdout |
| `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/imaris-n2v-full-<jobid>.err` | Slurm stderr |

Optional panel movies are generated after denoising from the OME-TIFF outputs:

| Path | Meaning |
| --- | --- |
| `denoised_1x3_zstack_panel_native.mp4` | Native 1x3 Z-scroll panel: magenta, green, merged |
| `denoised_1x3_zstack_panel_fixedscale_prores.mov` | Artifact-resistant ProRes QC movie |
| `denoised_1x3_zstack_panel_fixedscale_prores.mov.display_limits.json` | Fixed display limits used for movie |
| `denoised_1x3_zstack_panel_fixedscale_prores_stills/*.png` | Lossless diagnostic stills |

## Important Parameters

Channels:

- channel 0: green / `PV-mNG`
- channel 1: red, displayed as magenta / `BiVe3-dTom`

CAREamics mode:

- one N2V model per channel
- 3D patches over `Z x Y x X`
- full-run patch size: `16 x 64 x 64`
- full-run prediction tile overlap: `8 x 32 x 32`
- full-run epochs: `50`
- full-run steps per epoch: `200`
- batch size: `1`
- prediction checkpoint: `best`

The prediction overlap was increased to half the patch size to reduce possible
patch-boundary/grid artifacts during tiled prediction. This does not modify raw
data or change training labels; it only changes how full-stack prediction chunks
are overlapped and stitched.

## Sample Configs

Sample configs live in:

`imaging/imaris_careamics/config/`

Configs for the eight expected samples:

```text
cl32_bive4_pv_reporter_40x_F0.env
cl32_bive4_pv_reporter_40x_F1.env
cl32_bive4_pv_reporter_40x_F2.env
cl32_bive4_pv_reporter_40x_F3.env
cl32_bive4_pv_reporter_40x_realbive4_F0.env
cl32_bive4_pv_reporter_40x_realbive4_F1.env
cl32_bive4_pv_reporter_40x_realbive4_F2.env
cl32_bive4_pv_reporter_40x_realbive4_F3.env
```

Submit a single full run:

```bash
SAMPLE_CONFIG=/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline/imaging/imaris_careamics/config/<sample>.env \
sbatch slurm/run_full_denoising.sbatch
```

Do not resubmit `cl32_bive4_pv_reporter_40x_realbive4_F1` unless intentionally
rerunning it, because its full outputs already exist.

## Completed Pilot Outputs

Completed sample:

`cl32_bive4_pv_reporter_40x_realbive4_F1`

Output directory:

`/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_realbive4_F1/full_run/`

Notable outputs:

```text
green_raw.ome.tif
red_raw.ome.tif
green_denoised.ome.tif
red_denoised.ome.tif
merged_green_red_preview.ome.tif
max_projection_preview.mp4
denoised_1x3_zstack_panel_native.mp4
denoised_1x3_zstack_panel_fixedscale_prores.mov
qc/raw_vs_denoised_example_slices.png
qc/raw_vs_denoised_max_projections.png
qc/intensity_histograms.png
qc/clipping_saturation_warnings.txt
```

## Quick Completion Check

```bash
for d in /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_{F0,F1,F2,F3} \
         /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_realbive4_{F0,F1,F2,F3}; do
  echo "== $d =="
  ls -lh "$d/full_run"/{green_denoised.ome.tif,red_denoised.ome.tif,max_projection_preview.mp4} 2>/dev/null || true
  test -f "$d/full_run/qc/clipping_saturation_warnings.txt" && tail -n 5 "$d/full_run/qc/clipping_saturation_warnings.txt"
done
```
