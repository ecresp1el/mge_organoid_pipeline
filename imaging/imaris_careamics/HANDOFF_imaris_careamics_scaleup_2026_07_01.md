# Handoff: Imaris CAREamics N2V Scale-Up

Date: 2026-07-01 15:40 EDT
Updated: 2026-07-02 post-run validation

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
their denoised OME-TIFFs contain `NaN` values. The
`cl32_bive4_pv_reporter_40x_realbive4_F0-F3` denoised outputs are valid finite
image data.

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
TILE_OVERLAP_Z=8
TILE_OVERLAP_YX=32
```

The fix excludes leading/trailing all-zero Z planes from training/prediction,
then pads the denoised output back to the original Z size with zeros. For these
files, the nonzero Z range is `0-360`, with trailing zero planes `361-383`.

| Job ID | Sample | Output directory |
| --- | --- | --- |
| `52759673` | `cl32_bive4_pv_reporter_40x_F0` | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_F0_trimz_fix` |
| `52759674` | `cl32_bive4_pv_reporter_40x_F1` | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_F1_trimz_fix` |
| `52759675` | `cl32_bive4_pv_reporter_40x_F2` | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_F2_trimz_fix` |
| `52759676` | `cl32_bive4_pv_reporter_40x_F3` | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_F3_trimz_fix` |

Status command:

```bash
squeue -j 52759673,52759674,52759675,52759676 \
  -o '%.18i %.9P %.35j %.8u %.2t %.12M %.6D %.12b %.10m %R'
```

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
| `qc/raw_vs_denoised_example_slices.png` | Slice-level raw vs denoised QC |
| `qc/raw_vs_denoised_max_projections.png` | Max-projection raw vs denoised QC |
| `qc/intensity_histograms.png` | Raw/denoised intensity distributions |
| `qc/clipping_saturation_warnings.txt` | Saturation/non-finite warning report |
| `run_metadata.json` | Input, internal dimensions, args, reader metadata |

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
