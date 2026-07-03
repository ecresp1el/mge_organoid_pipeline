# BC43/realbive4 Stitching Coordinate Audit

Date: 2026-07-03

This note documents why the first Fiji Grid/Collection stitch was marked
wrong-layout and why the corrected-stage layout was used for the final full
volume stitch.

## Inputs

Input directory:

```text
/nfs/turbo/umms-parent/andor_micropscope_data_dump/exp17_pv_reporter_with_biver3and4/20x
```

TeraStitcher-style XML:

```text
cl32_bive4_pv_reporter_40x_realbive4.xml
```

Tiles:

```text
cl32_bive4_pv_reporter_40x_realbive4_F0.ims
cl32_bive4_pv_reporter_40x_realbive4_F1.ims
cl32_bive4_pv_reporter_40x_realbive4_F2.ims
cl32_bive4_pv_reporter_40x_realbive4_F3.ims
```

## What The XML Said

The XML `Stack` entries had these positions:

```text
F3: ABS_H=0    ABS_V=0
F0: ABS_H=918  ABS_V=0
F2: ABS_H=0    ABS_V=897
F1: ABS_H=918  ABS_V=897
```

Passing those directly to Fiji placed:

```text
top-left:     F3
top-right:    F0
bottom-left:  F2
bottom-right: F1
```

That first full-volume stitch completed, but MIP review showed the organoid was
not spatially coherent. The output was therefore preserved as:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/fiji_stitching/cl32_bive4_pv_reporter_40x_realbive4_xml_coords_WRONG_LAYOUT_DO_NOT_USE_52782973/
```

## What The IMS Metadata Said

Each IMS file contains `DataSetInfo/Image` physical extents. The relevant values
were:

```text
tile  X     Y    Z    ExtMin0  ExtMin1  ExtMax0  ExtMax1
F0    1020  996  736  19553.9  1256.42  19873.0 1568.0
F1    1020  996  736  19553.9   975.997 19873.0 1287.58
F2    1020  996  736  19841.0   975.997 20160.1 1287.58
F3    1020  996  736  19841.0  1256.42  20160.1 1568.0
```

These extents are internally consistent:

- F0/F1 have lower `ExtMin0` than F3/F2, so they are the left column in the
  physical coordinate system.
- F1/F2 have lower `ExtMin1` than F0/F3.
- The physical spacing implied by the extents corresponds to about 918 pixels in
  X and 896 pixels in Y after conversion back to pixel units.
- Bio-Formats/Fiji opened each IMS as `1020x996x736 channels=2 frames=1`, which
  matches the IMS metadata `X/Y/Z`. The raw HDF5 arrays are padded to
  `1024x1024x736`, so the valid image size is the metadata/Bio-Formats size, not
  the padded HDF5 array size.

This means the IMS metadata was not simply corrupt. It gave useful tile spacing
and the left/right relationship. The problem was coordinate-frame interpretation
when turning XML/physical coordinates into Fiji display positions.

## MIP QC Comparison

To avoid debugging the whole 736-plane volume blind, the script below generated
per-tile max projections and stitched them using several candidate layouts:

```text
imaging/imaris_careamics/mip_stitch_ims_montage_qc.py
```

MIP QC output directory:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/fiji_stitching/cl32_bive4_pv_reporter_40x_realbive4_mip_qc/
```

Useful review files:

```text
mip_placement_three_way_contact_sheet.png
stitched_mip_corrected_stage_rgb_display.png
positions.tsv
positions_corrected_stage.tsv
```

The three-way MIP comparison showed:

- XML placement was wrong.
- IMS physical extent placement fixed left/right compared with the XML layout.
- Visual organoid continuity still required flipping the row placement for
  display. The corrected-stage layout best preserved the expected sphere-like
  organoid shape.

## Corrected Layout Used For Full Stitch

The corrected-stage layout used for the final full-volume Fiji stitch was:

```text
# Define the number of dimensions we are working on
dim = 3

# Define the image coordinates
cl32_bive4_pv_reporter_40x_realbive4_F1.ims; ; (0.0, 0.0, 0.0)
cl32_bive4_pv_reporter_40x_realbive4_F2.ims; ; (918.0, 0.0, 0.0)
cl32_bive4_pv_reporter_40x_realbive4_F0.ims; ; (0.0, 896.0, 0.0)
cl32_bive4_pv_reporter_40x_realbive4_F3.ims; ; (918.0, 896.0, 0.0)
```

The tile config file is:

```text
/nfs/turbo/umms-parent/andor_micropscope_data_dump/exp17_pv_reporter_with_biver3and4/20x/TileConfiguration.corrected_stage_mip_qc.txt
```

Final corrected full-volume stitch:

```text
52805392  fiji-stitch-realbive4  COMPLETED  exit 0:0
```

Final stack-writing job:

```text
52805393  fiji-stack-tiff  COMPLETED  exit 0:0
```

Final single-series uint16 OME-TIFF stacks:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/fiji_stitching/cl32_bive4_pv_reporter_40x_realbive4_corrected_stage_full/tiff_stacks/cl32_bive4_pv_reporter_40x_realbive4_corrected_stage_c1_stitched_stack.ome.tif
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/imaging/fiji_stitching/cl32_bive4_pv_reporter_40x_realbive4_corrected_stage_full/tiff_stacks/cl32_bive4_pv_reporter_40x_realbive4_corrected_stage_c2_stitched_stack.ome.tif
```

Verification of each final stack:

```text
pages=736
series=1
shape=(736, 1892, 1938)
axes=ZYX
dtype=uint16
OME=true
```

## Interpretation

Do not describe this as "IMS metadata was wrong" without qualification. The
evidence is more specific:

- The XML `ABS_H/ABS_V` layout was not directly usable as Fiji display
  coordinates for this montage.
- The IMS metadata extents were internally consistent and helped identify the
  left/right correction.
- The final row orientation was selected by visual MIP QC against the expected
  organoid morphology.

For future montages, first generate MIP QC candidates and inspect the montage
geometry before launching full-volume Fiji fusion.
