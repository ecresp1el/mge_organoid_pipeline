# Python Notebook Entry Point

This directory adds a Python-first workflow for converting selected Seurat `.rds`
objects to AnnData `.h5ad` files without changing the existing R/Slurm pipeline.

The repo remains code/config only. Runtime data and generated `.h5ad` files stay
under:

```bash
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder
```

## Fresh Login Test

These instructions are the ground truth for testing the Python entry point from a
fresh Great Lakes login node session.

### 1. Start From A Clean Login Shell

```bash
ssh elcrespo@gl-login1.arc-ts.umich.edu
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
```

### 2. Load Conda

```bash
module purge
module load Anaconda3
```

Module names can vary on Great Lakes. If `module load Anaconda3` fails, inspect
the available modules and load the available Anaconda, Miniconda, or Mambaforge
module:

```bash
module avail anaconda
module avail mamba
module avail miniconda
```

### 3. Create Or Update The Conda Environment

Create the env the first time:

```bash
conda env create -f python_notebooks/environment.yml
```

If the env already exists and you changed `environment.yml`, update it:

```bash
conda env update -n mge-organoid-python -f python_notebooks/environment.yml --prune
```

Activate it:

```bash
conda activate mge-organoid-python
```

### 4. Register The Notebook Kernel

```bash
python -m ipykernel install --user \
  --name mge-organoid-python \
  --display-name "Python (mge-organoid-python)"
```

### 5. Set Runtime Root

```bash
export PROJECT_ROOT=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder
```

### 6. Verify The Three Seurat Inputs

```bash
ls -lh \
  "$PROJECT_ROOT/results/shi_2019_paper_qc/shi_2019_seurat.rds" \
  "$PROJECT_ROOT/results/varela_this_paper/varela_this_paper_seurat.rds" \
  "/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/umap_props_output/clustered_day90_with_cluster_names_2.rds"
```

This should show three existing `.rds` files. Do not continue to conversion if
one is missing.

### 7. Verify Python Imports And Path Resolution

```bash
PYTHONPATH="$PWD/python_notebooks/src" python - <<'PY'
from mge_organoid_python import default_studies, resolve_project_root, validate_source_paths

project_root = resolve_project_root()
studies = default_studies()
print("PROJECT_ROOT =", project_root)
for study in studies:
    print(study.study_id, "=>", study.seurat_path)

missing = validate_source_paths(studies)
if missing:
    raise SystemExit("Missing inputs: " + repr(missing))
print("All default Seurat inputs are present.")
PY
```

### 8. Launch Jupyter

Use the login node for light notebook editing and path checks only. The Seurat
objects can be large, so full conversion should be done in a GUI/interactive
compute session when possible.

```bash
jupyter lab --no-browser
```

Open:

```text
python_notebooks/notebooks/01_seurat_to_anndata.ipynb
```

Select the kernel:

```text
Python (mge-organoid-python)
```

Run the notebook from top to bottom. The first cells only validate imports and
paths; the conversion cell writes `.h5ad` files.

## Quick Relaunch After The Env Exists

For later sessions, after the environment has already been created:

```bash
ssh elcrespo@gl-login1.arc-ts.umich.edu
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
module purge
module load Anaconda3
conda activate mge-organoid-python
export PROJECT_ROOT=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder
jupyter lab --no-browser
```

## Canonical Inputs

The notebook uses these Seurat objects:

```text
Shi 2019 paper QC:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_2019_paper_qc/shi_2019_seurat.rds

Varela DIV30:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/varela_this_paper/varela_this_paper_seurat.rds

Varela DIV90:
/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/umap_props_output/clustered_day90_with_cluster_names_2.rds
```

## Outputs

Converted AnnData files are cached under:

```text
$PROJECT_ROOT/results/python_anndata/
```

Expected files:

```text
shi_2019_paper_qc.h5ad
varela_div30.h5ad
varela_div90.h5ad
```

The notebook loads these files back into memory as `AnnData` objects after
conversion.

## Python Package

The notebook imports the repo-local package from:

```text
python_notebooks/src/mge_organoid_python/
```

Main API:

```python
from mge_organoid_python import (
    SeuratToAnnDataConverter,
    default_studies,
    resolve_project_root,
)

project_root = resolve_project_root()
converter = SeuratToAnnDataConverter(project_root=project_root)
adatas, reports = converter.convert_many(default_studies())
```
