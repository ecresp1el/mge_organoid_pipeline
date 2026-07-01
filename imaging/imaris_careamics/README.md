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
  --require-gpu
```

Smoke outputs are written to `<output-dir>/smoke_test/`.
Full outputs are written to `<output-dir>/full_run/`.

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
EPOCHS=100 BATCH_SIZE=4 PATCH_SIZE_Z=16 PATCH_SIZE_YX=96 \
IMS_INPUT=/path/to/sample.ims \
OUTPUT_DIR=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/imaris_careamics/sample_001 \
sbatch slurm/run_full_denoising.sbatch
```

## Outputs

Each run directory contains:

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

Preview scaling uses 0.5 to 99.8 percentiles per channel across the full loaded
movie or stack, not per slice. Raw and denoised OME-TIFF outputs are not preview-scaled.

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
