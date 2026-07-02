# Imaris CAREamics N2V Denoising

This is a standalone imaging subpipeline for two-channel Imaris `.ims` confocal files.
It is separate from the scRNA-seq workflow, but follows the same Great Lakes habits:
code stays in the repo, runtime artifacts go under `PROJECT_ROOT/results/...`, and
Slurm logs are written under the selected output directory plus the project log root.

## Inputs

- Imaris `.ims` file
- Channel 0: green
- Channel 1: red

The input `.ims` file is read only and is never overwritten or modified.

Current sample:

- Local mounted path: `/Volumes/umms-parent/andor_micropscope_data_dump/exp17_pv_reporter_with_biver3and4/20x/cl32_bive4_pv_reporter_40x_realbive4_F1.ims`
- Great Lakes path: `/nfs/turbo/umms-parent/andor_micropscope_data_dump/exp17_pv_reporter_with_biver3and4/20x/cl32_bive4_pv_reporter_40x_realbive4_F1.ims`
- Output root: `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_realbive4_F1`
- Sample config: `imaging/imaris_careamics/config/cl32_bive4_pv_reporter_40x_realbive4_F1.env`

Current scale-up handoff:

- `imaging/imaris_careamics/HANDOFF_imaris_careamics_scaleup_2026_07_01.md`

The current eight-movie batch has two complete F0-F3 sets:

| Set | Samples |
| --- | --- |
| `cl32_bive4_pv_reporter_40x` | `F0`, `F1`, `F2`, `F3` |
| `cl32_bive4_pv_reporter_40x_realbive4` | `F0`, `F1`, `F2`, `F3` |

`cl32_bive4_pv_reporter_40x_realbive4_F1` was the completed pilot/full run.
The other seven were submitted as separate Slurm jobs on 2026-07-01.

## Environment

Conda install:

```bash
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
bash imaging/imaris_careamics/setup_env.sh
conda activate imaris-careamics
python imaging/imaris_careamics/check_env.py
```

Pip fallback:

```bash
python3 -m venv .venv-imaris-careamics
source .venv-imaris-careamics/bin/activate
pip install -r imaging/imaris_careamics/requirements.txt
python imaging/imaris_careamics/check_env.py
```

For GPU runs, confirm that `check_env.py` prints `CUDA available: true` and a GPU name.
The default reader is `h5py`, which reads the `.ims` HDF5 structure directly and supports efficient smoke-test subsets.

## Interactive Commands

Environment check:

```bash
conda activate imaris-careamics
python imaging/imaris_careamics/check_env.py
```

Smoke test:

```bash
python imaging/imaris_careamics/denoise_ims_careamics.py \
  --input /path/to/sample.ims \
  --output-dir /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/sample_001 \
  --subset \
  --max-timepoints 2 \
  --max-z 16 \
  --crop-yx 512 \
  --epochs 3 \
  --batch-size 2 \
  --patch-size-z 16 \
  --patch-size-yx 64 \
  --tile-overlap-z 8 \
  --tile-overlap-yx 32 \
  --predict-checkpoint best \
  --require-gpu
```

Full run:

```bash
python imaging/imaris_careamics/denoise_ims_careamics.py \
  --input /path/to/sample.ims \
  --output-dir /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/sample_001 \
  --epochs 50 \
  --batch-size 2 \
  --patch-size-z 16 \
  --patch-size-yx 64 \
  --tile-overlap-z 8 \
  --tile-overlap-yx 32 \
  --require-gpu
```

Smoke outputs are written to `<output-dir>/smoke_test/`.
Full outputs are written to `<output-dir>/full_run/`.

By default, prediction uses the best validation-loss checkpoint
(`--predict-checkpoint best`) instead of the final epoch. This is important for
N2V runs because the model can overtrain even when outputs remain finite.

Prediction-only rerun from existing trained models:

```bash
python imaging/imaris_careamics/denoise_ims_careamics.py \
  --input /path/to/sample.ims \
  --output-dir /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/sample_001_bestckpt_repredict \
  --model-source-dir /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/sample_001/full_run/models \
  --predict-checkpoint best \
  --patch-size-z 16 \
  --patch-size-yx 64 \
  --tile-overlap-z 8 \
  --tile-overlap-yx 32 \
  --require-gpu
```

Native-resolution denoised panel movie:

```bash
python imaging/imaris_careamics/create_denoised_panel_movie.py \
  --green /path/to/full_run/green_denoised.ome.tif \
  --red /path/to/full_run/red_denoised.ome.tif \
  --output /path/to/full_run/denoised_1x3_zstack_panel_native.mp4
```

Artifact-resistant QC movie plus lossless stills:

```bash
python imaging/imaris_careamics/create_denoised_panel_movie.py \
  --green /path/to/full_run/green_denoised.ome.tif \
  --red /path/to/full_run/red_denoised.ome.tif \
  --output /path/to/full_run/denoised_1x3_zstack_panel_fixedscale_prores.mov \
  --codec-mode prores \
  --export-stills
```

8-bit display conversion after N2V:

```bash
python imaging/imaris_careamics/convert_denoised_to_8bit.py \
  --input-dir /path/to/full_run \
  --overwrite
```

This writes display/QC copies into `/path/to/full_run/display_8bit/`:

```text
green_denoised_8bit_display.ome.tif
red_denoised_8bit_display.ome.tif
merged_green_magenta_8bit_display.ome.tif
denoised_8bit_green_max_projection.png
denoised_8bit_red_as_magenta_max_projection.png
denoised_8bit_merged_green_magenta_max_projection.png
display_scaling_8bit.json
```

The 8-bit files are for Fiji/QC/display. They are not the quantitative N2V
outputs. The original `green_denoised.ome.tif` and `red_denoised.ome.tif` remain
unchanged. Display scaling uses fixed 0.5 to 99.8 percentile limits per channel
from a deterministic whole-stack sample by default
(`--percentile-z-step 2 --percentile-yx-step 2`), and the exact limits are saved
in `display_scaling_8bit.json`.

## Slurm Commands

Copy or submit from the repo root on Great Lakes:

```bash
mkdir -p /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs
sbatch slurm/check_gpu.sbatch
```

Sample-specific smoke test for `cl32_bive4_pv_reporter_40x_realbive4_F1`:

```bash
SAMPLE_CONFIG=/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline/imaging/imaris_careamics/config/cl32_bive4_pv_reporter_40x_realbive4_F1.env \
sbatch slurm/run_smoke_test.sbatch
```

Sample-specific full run:

```bash
SAMPLE_CONFIG=/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline/imaging/imaris_careamics/config/cl32_bive4_pv_reporter_40x_realbive4_F1.env \
sbatch slurm/run_full_denoising.sbatch
```

```bash
IMS_INPUT=/path/to/sample.ims \
OUTPUT_DIR=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/sample_001 \
sbatch slurm/run_smoke_test.sbatch
```

```bash
IMS_INPUT=/path/to/sample.ims \
OUTPUT_DIR=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/sample_001 \
sbatch slurm/run_full_denoising.sbatch
```

Useful Slurm environment overrides:

```bash
EPOCHS=100 BATCH_SIZE=4 PATCH_SIZE_Z=16 PATCH_SIZE_YX=96 PREDICT_CHECKPOINT=best \
IMS_INPUT=/path/to/sample.ims \
OUTPUT_DIR=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/sample_001 \
sbatch slurm/run_full_denoising.sbatch
```

Slurm prediction-only rerun from existing models:

```bash
IMS_INPUT=/path/to/sample.ims \
OUTPUT_DIR=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/sample_001_bestckpt_repredict \
MODEL_SOURCE_DIR=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/sample_001/full_run/models \
PREDICT_CHECKPOINT=best \
sbatch slurm/run_full_denoising.sbatch
```

Slurm 8-bit display conversion after a completed full run:

```bash
FULL_RUN_DIR=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/sample_001/full_run \
sbatch slurm/convert_denoised_to_8bit.sbatch
```

## Outputs

Each sample writes outputs under:

`/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/<sample>/full_run/`

Each full run directory contains:

- `green_raw.ome.tif`
- `red_raw.ome.tif`
- `green_denoised.ome.tif`
- `red_denoised.ome.tif`
- `merged_green_red_preview.ome.tif`
- `max_projection_preview.mp4`
- `models/green_n2v/`
- `models/red_n2v/`
- `qc/raw_vs_denoised_example_slices.png`
- `qc/raw_vs_denoised_max_projections.png`
- `qc/intensity_histograms.png`
- `qc/clipping_saturation_warnings.txt`
- `run_metadata.json`

Optional 8-bit display exports under `full_run/display_8bit/`:

- `green_denoised_8bit_display.ome.tif`
- `red_denoised_8bit_display.ome.tif`
- `merged_green_magenta_8bit_display.ome.tif`
- `denoised_8bit_green_max_projection.png`
- `denoised_8bit_red_as_magenta_max_projection.png`
- `denoised_8bit_merged_green_magenta_max_projection.png`
- `display_scaling_8bit.json`

Logs are written to:

- `<output-dir>/logs/run_full_denoising_<jobid>.log`
- `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/imaris-n2v-full-<jobid>.out`
- `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/imaris-n2v-full-<jobid>.err`

Preview scaling uses 0.5 to 99.8 percentiles per channel across the full loaded
movie or stack, not per slice. Raw and denoised OME-TIFF outputs are not preview-scaled.

The optional panel movie is a native-resolution Z-scroll with magenta `BiVe3-dTom`,
green `PV-mNG`, and merged views. The script writes a `.display_limits.json`
sidecar so the fixed display limits are recorded. Use `--codec-mode prores` or
`--codec-mode lossless_rgb` plus `--export-stills` when checking for denoising
grids, because standard H.264 can introduce block or chroma artifacts.

Prediction uses overlapping 3D chunks. The default overlap is half the patch size
(`8 x 32 x 32` for `16 x 64 x 64` patches), which is safer for avoiding visible
patch-boundary artifacts than a small overlap.

Quick status check for the current scale-up jobs:

```bash
squeue -j 52693455,52693456,52693458,52693459,52693460,52693461,52693462 \
  -o '%.18i %.9P %.35j %.8u %.2t %.12M %.6D %.12b %.10m %R'
```

Quick output-completion check:

```bash
for d in /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_{F0,F1,F2,F3} \
         /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/cl32_bive4_pv_reporter_40x_realbive4_{F0,F1,F2,F3}; do
  echo "== $d =="
  ls -lh "$d/full_run"/{green_denoised.ome.tif,red_denoised.ome.tif,max_projection_preview.mp4} 2>/dev/null || true
done
```

## Fiji Grid/Collection Stitching

Use the TeraStitcher XML coordinates to generate Fiji's
`TileConfiguration.txt`; do not infer the F0-F3 montage order by filename.

Current BC43/realbive4 montage:

```bash
python3 imaging/imaris_careamics/terastitcher_to_fiji_tile_config.py \
  --xml /nfs/turbo/umms-parent/andor_micropscope_data_dump/exp17_pv_reporter_with_biver3and4/20x/cl32_bive4_pv_reporter_40x_realbive4.xml \
  --output /nfs/turbo/umms-parent/andor_micropscope_data_dump/exp17_pv_reporter_with_biver3and4/20x/TileConfiguration.txt
```

Great Lakes Slurm prep job:

```bash
INPUT_DIR=/nfs/turbo/umms-parent/andor_micropscope_data_dump/exp17_pv_reporter_with_biver3and4/20x \
XML_NAME=cl32_bive4_pv_reporter_40x_realbive4.xml \
OUTPUT_NAME=TileConfiguration.txt \
sbatch slurm/prepare_fiji_stitching_tile_config.sbatch
```

Actual batch stitching job:

```bash
sbatch slurm/run_fiji_grid_stitching_realbive4.sbatch
```

The batch stitching job runs Fiji under `xvfb-run` and writes the fused image to
disk under:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/fiji_stitching/cl32_bive4_pv_reporter_40x_realbive4_xml_coords_compute_overlap/
```

For the batch run, `Image_output` is `Write to disk` so the result survives the
Slurm session. The interactive first-test setting below remains `Fuse and
display`.

The expected XML-derived layout is:

```text
# Define the number of dimensions we are working on
dim = 3

# Define the image coordinates
cl32_bive4_pv_reporter_40x_realbive4_F3.ims; ; (0.0, 0.0, 0.0)
cl32_bive4_pv_reporter_40x_realbive4_F0.ims; ; (918.0, 0.0, 0.0)
cl32_bive4_pv_reporter_40x_realbive4_F2.ims; ; (0.0, 897.0, 0.0)
cl32_bive4_pv_reporter_40x_realbive4_F1.ims; ; (918.0, 897.0, 0.0)
```

Launch Fiji from the Great Lakes desktop session:

```bash
module load fiji/1.5.4
/sw/pkgs/med/fiji/1.5.4/ImageJ-linux64 &
```

In Fiji, use `Plugins -> Stitching -> Grid/Collection Stitching`:

- Type: `Positions from file`
- Layout file: generated `TileConfiguration.txt`
- Fusion method: `Linear Blending`
- Compute overlap: checked for the first test
- Regression threshold: `0.30`
- Max/avg displacement threshold: `2.50`
- Absolute displacement threshold: `3.50`
- Image output: `Fuse and display`

Do not use the FusionStitcher stitched output, and do not split channels before
stitching.

## Workflow

```text
Raw .ims confocal stack
  |
  v
Read as T x C x Z x Y x X
  |
  v
Split green channel 0 and red channel 1
  |
  v
Train separate CAREamics N2V models with 3D ZYX patches
  |
  v
Predict every timepoint separately
  |
  v
Save OME-TIFF outputs, merged preview, MP4 max projection, models, and QC
```
