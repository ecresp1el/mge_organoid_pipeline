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

### 2. Confirm Which Conda Is Available

Do this before creating or updating anything.

```bash
conda --version
```

Expected result:

```text
conda <version>
```

If you see a conda version, conda is already accessible. On the current Great
Lakes login node setup, this may come from your home Miniconda install rather
than from a module.

Confirm exactly which conda is being used:

```bash
type -a conda
which conda
conda info --base
```

Expected pattern for your current setup:

```text
/home/elcrespo/miniconda3
```

It is also okay if `type -a conda` reports that `conda` is a shell function;
that is normal when conda has been initialized in your shell. The important
check is that `conda info --base` points to the conda installation you intend to
use.

If `conda --version` says `conda: command not found`, then try loading a conda
module:

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

After loading any module, rerun:

```bash
conda --version
which conda
conda info --base
```

Do not continue until those commands work.

### 3. Check Whether The Environment Already Exists

List conda environments:

```bash
conda env list
```

Expected result if the env already exists:

```text
mge-organoid-python    <some path>
```

Expected result if the env does not exist:

```text
# conda environments:
base                  * /home/elcrespo/miniconda3
...
```

If `mge-organoid-python` is not listed, create it in step 4A.

If `mge-organoid-python` is listed, update it in step 4B.

### 4A. Create The Conda Environment

Run this only if `conda env list` did not show `mge-organoid-python`:

```bash
conda env create -f python_notebooks/environment.yml
```

Expected result:

```text
done
#
# To activate this environment, use
#
#     $ conda activate mge-organoid-python
```

After creation, confirm it exists:

```bash
conda env list
```

You should now see `mge-organoid-python`.

### 4B. Update The Conda Environment

Run this only if `conda env list` already showed `mge-organoid-python`:

```bash
conda env update -n mge-organoid-python -f python_notebooks/environment.yml --prune
```

Expected result:

```text
done
```

### 5. Activate The Conda Environment

```bash
conda activate mge-organoid-python
```

Confirm the active environment:

```bash
echo "$CONDA_DEFAULT_ENV"
which python
python --version
```

Expected result:

```text
mge-organoid-python
```

`which python` should point inside the `mge-organoid-python` environment, not to
system Python.

### 6. Register The Notebook Kernel

```bash
python -m ipykernel install --user \
  --name mge-organoid-python \
  --display-name "Python (mge-organoid-python)"
```

### 7. Set Runtime Root

```bash
export PROJECT_ROOT=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder
```

### 8. Verify The Three Seurat Inputs

```bash
ls -lh \
  "$PROJECT_ROOT/results/shi_2019_paper_qc/shi_2019_seurat.rds" \
  "$PROJECT_ROOT/results/varela_this_paper/varela_this_paper_seurat.rds" \
  "/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/umap_props_output/clustered_day90_with_cluster_names_2.rds"
```

This should show three existing `.rds` files. Do not continue to conversion if
one is missing.

### 9. Verify Python Imports And Path Resolution

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

### 10. Launch Jupyter

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
conda --version
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
