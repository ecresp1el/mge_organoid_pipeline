# Python Notebook Entry Point

This directory adds a Python-first workflow for converting selected Seurat `.rds`
objects to AnnData `.h5ad` files without changing the existing R/Slurm pipeline.

The repo remains code/config only. Runtime data and generated `.h5ad` files stay
under:

```bash
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder
```

## Resource Model

Use VS Code on the login node only for editing, setup, and light validation.
Great Lakes login nodes are not the place to convert large Seurat objects or
run memory-heavy AnnData work.

For real notebook work, request an interactive Slurm allocation, then connect VS
Code to the allocated compute node. The notebook interface stays in VS Code, but
the Python kernel runs on the compute node with the RAM/CPU resources requested
from Slurm.

Current Slurm facts confirmed from the login node:

```text
Available accounts: louisdan0, parent0
Default account for this workflow: parent0
Great Lakes user/account login: elcrespo
standard partition: roughly 180 GB RAM per node
largemem partition: roughly 1.5 TB RAM per node
```

GPU is not needed for the current Seurat `.rds` to AnnData `.h5ad` conversion.
Start with CPU/RAM on `standard`; use `largemem` only if conversion or loading
fails due to memory.

## VS Code On A Compute Node

The standard project helper for requesting notebook resources is:

```bash
python_notebooks/scripts/start_vscode_compute_job.sh
```

This uses the project default Slurm account:

```text
parent0
```

To override the account intentionally:

```bash
python_notebooks/scripts/start_vscode_compute_job.sh --account louisdan0
```

Default request:

```text
partition: standard
memory:    128GB
CPUs:      4
time:      04:00:00
```

For a larger memory test:

```bash
python_notebooks/scripts/start_vscode_compute_job.sh \
  --partition largemem \
  --mem 512GB \
  --time 04:00:00
```

Keep the terminal running after `salloc` grants the allocation. Note the
allocated compute node name from the Slurm output. Then connect VS Code to that
compute node using your Great Lakes SSH jump-host setup. Once VS Code is
connected to the compute node, use the same `mge-organoid-python` environment
and notebook kernel there.

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

Confirmed Great Lakes state for this project:

```text
conda --version      -> conda 25.3.1
conda info --base    -> /home/elcrespo/miniconda3
type -a conda        -> shell function plus /home/elcrespo/miniconda3/condabin/conda
```

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

Confirmed Great Lakes state for this project:

```text
mge-organoid-python    /home/elcrespo/miniconda3/envs/mge-organoid-python
```

Because this env already exists, use step 4B, not step 4A.

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

It is okay if conda also prints:

```text
WARNING: A newer version of conda exists.
```

Do not update base conda for this workflow. Continue to step 5.

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

Confirmed Great Lakes state for this project:

```text
echo "$CONDA_DEFAULT_ENV" -> mge-organoid-python
which python              -> ~/miniconda3/envs/mge-organoid-python/bin/python
python --version          -> Python 3.11.15
```

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
