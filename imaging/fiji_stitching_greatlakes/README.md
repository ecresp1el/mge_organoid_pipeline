# Portable Fiji Stitching for Great Lakes

This directory performs exactly one workflow:

```text
pre-stitch tiles -> Fiji Grid/Collection Stitching -> stitched OME-TIFF stacks
```

It does **not** make a movie, denoise, filter, rescale, or perform downstream
analysis. It does not use BigStitcher or FusionStitcher.

Each submission creates a new isolated run directory. Inputs are staged as
read-only symlinks and are never modified. The run saves its exact code and
configuration, Fiji's plane outputs, final channel OME-TIFF stacks, logs,
validation results, and checksums together.

## Great Lakes requirements

The pipeline uses centrally installed modules only:

- `fiji/1.5.4`
- `python3.11-anaconda/2024.02`

No personal Conda environment or repository path is required after submission.

## First prove Fiji stitching with synthetic data

```bash
cd fiji_stitching_greatlakes
./submit.sh configs/smoke.env
```

This makes a known two-channel 3D volume, cuts it into four overlapping tiles,
stitches them with the included Fiji runner, packages the planes as OME-TIFF,
and compares every output voxel with the known original. A valid smoke test
reports the following for both channels in `validation.json`:

```text
max_abs_difference: 0
exact_fraction: 1.0
```

Use the run path printed by `submit.sh` to inspect status:

```bash
./status.sh /path/to/the/run_directory
```

## Reproduce the corrected BC43/realbive4 stitch

Review `configs/realbive4.env`, especially `ACCOUNT`, `INPUT_DIR`, and
`RUN_PARENT`, then run:

```bash
./submit.sh configs/realbive4.env --dry-run
./submit.sh configs/realbive4.env
```

The supplied, visually reviewed tile layout is:

```text
F1 at (0,   0,   0)
F2 at (918, 0,   0)
F0 at (0,   896, 0)
F3 at (918, 896, 0)
```

Do not substitute the literal XML layout. The XML-derived job `52782973`
finished computationally but had the wrong montage orientation. The original
successful corrected stitch was job `52805392`.

## Outputs

```text
RUN_DIRECTORY/
├── code/                       exact bundle snapshot used by the job
├── config/                     submitted/resolved config and environment
├── inputs/                     copied layout and read-only tile symlinks
├── logs/                       Slurm stdout/stderr and pipeline log
├── work/classes/               compiled Java runner
├── results/planes/             Fiji fused planes: img_t*_z*_c*
├── results/stacks/             one stitched OME-TIFF per channel
├── validation.json             plane/stack integrity checks
├── manifest.json               inventory and SHA-256 checksums
└── SUCCESS.txt or FAILED.txt   unambiguous terminal status
```

The packaging step changes only the container: it streams Fiji's fused planes
into one `uint16`, `ZYX`, single-series OME-TIFF per channel. It does not change
pixel intensities.

## Adapt to another tile set

Copy `configs/realbive4.env` and change:

1. `INPUT_DIR` to the tile directory.
2. `LAYOUT_FILE` to a reviewed Fiji `TileConfiguration.txt`.
3. `RUN_PARENT`, `ACCOUNT`, and resource requests.
4. `OUTPUT_PREFIX` and `CHANNELS`.

The filenames in the layout must exactly match files in `INPUT_DIR`. Review the
layout geometry before a full-volume submission; coordinate-frame conventions
vary between microscopes.

## Share the code

```bash
./make_release.sh
```

This writes a `.tar.gz` archive and a SHA-256 file. The recipient can verify it:

```bash
sha256sum -c fiji_stitching_greatlakes_YYYYMMDD.tar.gz.sha256
```
