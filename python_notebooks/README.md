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

## Current Step-By-Step Notebook Test

Use this section when starting from the confirmed state:

```text
login node: gl-login1
repo: /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
conda base: /home/elcrespo/miniconda3
env: mge-organoid-python
python: Python 3.11.15
default Slurm account: parent0
```

### A. Login-Node Setup Checks

These commands are safe on the login node because they only register the kernel
and validate paths/imports.

Goal: prove that the login session can see the repo, the conda env, the Jupyter
kernel, and the source files. This does not perform Seurat conversion.

```bash
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
conda activate mge-organoid-python
```

Register the kernel:

```bash
python -m ipykernel install --user \
  --name mge-organoid-python \
  --display-name "Python (mge-organoid-python)"
```

Expected output:

```text
Installed kernelspec mge-organoid-python in /home/elcrespo/.local/share/jupyter/kernels/mge-organoid-python
```

Why this matters: VS Code and Jupyter discover notebook environments through
kernel specs. Registering this once lets the notebook offer
`Python (mge-organoid-python)` as a selectable kernel.

Confirmed Great Lakes state for this project:

```text
Installed kernelspec mge-organoid-python in /home/elcrespo/.local/share/jupyter/kernels/mge-organoid-python
```

Set the runtime root:

```bash
export PROJECT_ROOT=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder
```

Confirm it:

```bash
echo "$PROJECT_ROOT"
```

Expected and confirmed value:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder
```

Why this matters: repo files live under `/home`, but data/results live under
`PROJECT_ROOT`. The Python code uses `PROJECT_ROOT` to keep generated `.h5ad`
files out of the git repo.

Confirm the three source objects:

```bash
ls -lh \
  "$PROJECT_ROOT/results/shi_2019_paper_qc/shi_2019_seurat.rds" \
  "$PROJECT_ROOT/results/varela_this_paper/varela_this_paper_seurat.rds" \
  "/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/umap_props_output/clustered_day90_with_cluster_names_2.rds"
```

Expected result: three `.rds` files are listed.

Why this matters: these are the exact Seurat inputs the notebook will convert.
Checking them before requesting compute resources catches path problems early.

Confirmed Great Lakes state for this project:

```text
Shi 2019 paper QC: 607M
Varela DIV30:      16G
Varela DIV90:      3.1G
```

Because the Varela objects are large, do not run conversion on the login node.

Confirm the Python package sees the same paths:

```bash
PYTHONPATH="$PWD/python_notebooks/src" python - <<'PY'
from mge_organoid_python import default_studies, resolve_project_root, validate_source_paths

print("PROJECT_ROOT =", resolve_project_root())
studies = default_studies()
for study in studies:
    print(study.study_id, "=>", study.seurat_path)

missing = validate_source_paths(studies)
if missing:
    raise SystemExit("Missing inputs: " + repr(missing))
print("All default Seurat inputs are present.")
PY
```

Stop here if any source object is missing.

Why this matters: this verifies the notebook-facing Python package imports
correctly and resolves the same canonical inputs that were checked with `ls`.

Confirmed Great Lakes state for this project:

```text
PROJECT_ROOT = /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder
shi_2019_paper_qc => /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_2019_paper_qc/shi_2019_seurat.rds
varela_div30 => /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/varela_this_paper/varela_this_paper_seurat.rds
varela_div90 => /nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/umap_props_output/clustered_day90_with_cluster_names_2.rds
All default Seurat inputs are present.
```

At this point, login-node setup validation is complete. Continue to resource
allocation before running notebook conversion.

### B. Request Notebook Resources

From the login-node VS Code terminal:

```bash
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
python_notebooks/scripts/start_vscode_compute_job.sh
```

This requests the default notebook allocation:

```text
account:   parent0
partition: standard
memory:    128GB
CPUs:      4
time:      04:00:00
```

Why this matters: the notebook conversion reads 607 MB, 3.1 GB, and 16 GB RDS
objects. The allocation moves the actual notebook kernel work off `gl-login1`
and onto a Slurm compute node with requested RAM/CPUs.

When Slurm grants the allocation, keep that terminal open and note the compute
node name.

Expected Slurm behavior:

```text
salloc: Granted job allocation <job_id>
salloc: Nodes <node_name> are ready for job
```

Confirmed Great Lakes state for this project:

```text
salloc: Pending job allocation 50282695
salloc: job 50282695 queued and waiting for resources
salloc: job 50282695 has been allocated resources
salloc: Granted job allocation 50282695
salloc: Nodes gl3103 are ready for job
```

In this confirmed example, the compute node is:

```text
gl3103
```

If the job remains queued, wait. If Slurm reports an account, partition, memory,
or time error, stop and fix the request before continuing.

### C. Connect VS Code To The Compute Node

Use VS Code Remote SSH to connect to the allocated compute node through the
Great Lakes login node jump host. The notebook must run from this compute-node
VS Code session, not from the `gl-login1` session.

Why this matters: opening the notebook from the login-node VS Code window would
start the kernel on `gl-login1`. Connecting VS Code to the allocated compute
node makes the notebook kernel run where the requested resources exist.

In the compute-node VS Code terminal, run:

```bash
hostname
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
conda activate mge-organoid-python
export PROJECT_ROOT=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder
```

`hostname` should show the allocated compute node, not `gl-login1`.

Confirm the env and runtime root on the compute node:

```bash
echo "$CONDA_DEFAULT_ENV"
which python
echo "$PROJECT_ROOT"
```

Expected output:

```text
mge-organoid-python
.../envs/mge-organoid-python/bin/python
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder
```

Confirmed Great Lakes state for this project:

```text
hostname                  -> gl3103.arc-ts.umich.edu
echo "$CONDA_DEFAULT_ENV" -> mge-organoid-python
which python              -> ~/miniconda3/envs/mge-organoid-python/bin/python
python --version          -> Python 3.11.15
echo "$PROJECT_ROOT"      -> /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder
```

This confirms the terminal is running on the allocated compute node, not on the
login node. It also confirms the notebook kernel environment and runtime data
root are available on the compute node.

Open:

```text
python_notebooks/notebooks/01_seurat_to_anndata.ipynb
```

Select:

```text
Python (mge-organoid-python)
```

If `Python (mge-organoid-python)` is not visible in VS Code, verify kernel
registration from the compute-node terminal:

```bash
hostname
conda activate mge-organoid-python
python -m jupyter kernelspec list
cat /home/elcrespo/.local/share/jupyter/kernels/mge-organoid-python/kernel.json
```

Expected meaning:

```text
hostname
```

should show the allocated compute node, not `gl-login1`.

```text
python -m jupyter kernelspec list
```

should include:

```text
mge-organoid-python    /home/elcrespo/.local/share/jupyter/kernels/mge-organoid-python
```

The `kernel.json` file should point to:

```text
/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/python
```

Confirmed Great Lakes state for this project:

```text
hostname -> gl3103.arc-ts.umich.edu

python -m jupyter kernelspec list includes:
mge-organoid-python    /home/elcrespo/.local/share/jupyter/kernels/mge-organoid-python

kernel.json argv includes:
/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/python

kernel.json display_name:
Python (mge-organoid-python)
```

If the kernelspec is missing, rerun:

```bash
python -m ipykernel install --user \
  --name mge-organoid-python \
  --display-name "Python (mge-organoid-python)"
```

Then close and reopen the notebook tab or use the VS Code command palette:

```text
Developer: Reload Window
```

If the kernel picker opens but does not show the env, choose:

```text
Select Another Kernel -> Python Environments -> Enter interpreter path
```

Then paste:

```text
/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/python
```

If clicking `Select Kernel` does not open anything at all, use the command
palette instead:

```text
Notebook: Select Notebook Kernel
```

If that still does not open a picker, reload the VS Code window and try again:

```text
Developer: Reload Window
```

Confirmed VS Code state for this project:

```text
Notebook kernel selected: mge-organoid-python (Python 3.11.15)
```

First notebook cell expected output:

```text
PosixPath('/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline')
```

Confirmed notebook state for this project:

```text
PosixPath('/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline')
```

Second notebook cell expected output:

```text
PROJECT_ROOT = /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder
shi_2019_paper_qc: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_2019_paper_qc/shi_2019_seurat.rds
varela_div30: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/varela_this_paper/varela_this_paper_seurat.rds
varela_div90: /nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/umap_props_output/clustered_day90_with_cluster_names_2.rds
```

If VS Code shows `Import "mge_organoid_python" could not be resolved` but the
cell runs and prints the paths above, Python is working correctly. That warning
comes from the editor language server, not the notebook runtime. The repo has a
`pyrightconfig.json` file that points VS Code/Pylance at `python_notebooks/src`;
reload the VS Code window if the warning does not clear.

Plain-language explanation of `pyrightconfig.json`:

```json
{
  "include": [
    "python_notebooks/src",
    "python_notebooks/notebooks"
  ],
  "extraPaths": [
    "python_notebooks/src"
  ]
}
```

`pyrightconfig.json` is for VS Code/Pylance static checking only. It does not
run the notebook and it does not change where Python looks at runtime.

`include` tells Pylance which repo folders to analyze. Here, it means:

```text
python_notebooks/src       -> analyze the custom Python package
python_notebooks/notebooks -> analyze notebook-side Python code
```

`extraPaths` tells Pylance where importable local packages live. Here, it means:

```text
python_notebooks/src
```

so the editor can understand:

```python
import mge_organoid_python
```

The notebook runtime still gets the same path from the first notebook cell,
which adds `python_notebooks/src` to `sys.path`. In short:

```text
pyrightconfig.json -> helps VS Code stop showing a false import warning
first notebook cell -> makes the package importable when the notebook runs
```

Run notebook cells in order until the smoke-test conversion cell. Do not use
`Run All` for the first test.

Before conversion, the notebook includes a resource diagnostic cell. Expected
output includes:

```text
Resource diagnostics before conversion
hostname: gl3121.arc-ts.umich.edu
python: /home/elcrespo/miniconda3/envs/mge-organoid-python/bin/python
PROJECT_ROOT: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder
Rscript: /home/elcrespo/miniconda3/envs/mge-organoid-python/bin/Rscript
SLURM_JOB_ID: ...
SLURM_CPUS_PER_TASK: 4
CPU count: ...
/proc/meminfo first lines:
MemTotal: ...
```

Why this matters: this proves the notebook is using the compute-node kernel and
shows the memory/Slurm context before any Seurat conversion starts.

Confirmed good notebook diagnostic for this project:

```text
Resource diagnostics before conversion
hostname: gl3121.arc-ts.umich.edu
python: /home/elcrespo/miniconda3/envs/mge-organoid-python/bin/python
PROJECT_ROOT: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder
Rscript: /home/elcrespo/miniconda3/envs/mge-organoid-python/bin/Rscript
SLURM_JOB_ID: 50282995
SLURM_JOB_NODELIST: gl3121
SLURM_CPUS_PER_TASK: 4
SLURM_MEM_PER_NODE: 131072
CPU count: 36
MemTotal: ~196 GB
MemAvailable: ~190 GB
```

This is the state required before running conversion.

The main conversion cell now starts with Shi only:

```python
study = studies_by_id["shi_2019_paper_qc"]
```

Why this matters: Shi is much smaller than the Varela objects. This first proves
that the R bridge, Seurat loading, UMAP transfer, H5AD writing, and AnnData
loading all work before attempting the 16 GB Varela object.

The smoke-test cell writes one output to:

```text
$PROJECT_ROOT/results/python_anndata/
```

The conversion cell should print progress before and during work. Expected
Python-side messages look like:

```text
Smoke test: converting/loading one small study first.
Cache directory: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/python_anndata
[YYYY-MM-DD HH:MM:SS] Study shi_2019_paper_qc: source=... target=... needs_conversion=True
[YYYY-MM-DD HH:MM:SS] Study shi_2019_paper_qc: starting RDS -> H5AD conversion
```

Expected R-side messages during conversion look like:

```text
[R YYYY-MM-DD HH:MM:SS] RETICULATE_PYTHON: /home/elcrespo/miniconda3/envs/mge-organoid-python/bin/python
[R YYYY-MM-DD HH:MM:SS] RETICULATE_AUTOCONFIGURE: FALSE
[R YYYY-MM-DD HH:MM:SS] Reading Seurat RDS: ...
[R YYYY-MM-DD HH:MM:SS] Loaded Seurat object with ... cells and ... features
[R YYYY-MM-DD HH:MM:SS] Converting Seurat object to SingleCellExperiment
[R YYYY-MM-DD HH:MM:SS] Transferring reduction to reducedDim X_umap: umap
[R YYYY-MM-DD HH:MM:SS] Writing H5AD: ...
[R YYYY-MM-DD HH:MM:SS] Finished writing H5AD: ...
```

Implementation note: conversion now runs in an external `Rscript` subprocess
instead of embedding R directly inside the notebook kernel. This matters because
native R, Seurat, or zellkonverter failures can crash an in-process `rpy2`
kernel without producing a normal Python traceback. With `Rscript`, the notebook
should keep running and show the R log/error if conversion fails.

`zellkonverter` uses `reticulate` internally. The converter forces reticulate to
use the existing notebook conda environment:

```text
RETICULATE_PYTHON=/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/python
RETICULATE_AUTOCONFIGURE=FALSE
```

Why this matters: without these variables, reticulate may try to download and
build its own Python through `pyenv`, which is slow and not appropriate for this
workflow. If you see output like `Installing pyenv` or `Installing
Python-3.14.0`, interrupt the notebook cell, restart the kernel, and rerun with
the updated converter.

Having 128 GB RAM is necessary but not sufficient to prevent every kernel crash:
kernel crashes can come from native-library segfaults, R/Python ABI conflicts,
or subprocesses being killed by memory pressure. The resource diagnostic plus
external `Rscript` path makes those cases easier to distinguish.

If the diagnostic shows this pattern, stop:

```text
hostname: gl-login3.arc-ts.umich.edu
PROJECT_ROOT: None
SLURM_JOB_ID: None
SLURM_JOB_NODELIST: None
SLURM_CPUS_PER_TASK: None
SLURM_MEM_PER_NODE: None
```

What it means:

```text
The notebook kernel is running on a login node, not on the allocated compute node.
The notebook kernel also did not inherit PROJECT_ROOT.
Do not run conversion in this state.
```

The notebook now raises an error if `hostname` starts with `gl-login` before any
conversion starts. The converter also refuses to convert on login nodes by
default.

Fix:

```text
1. Keep the salloc terminal open.
2. Use VS Code Remote SSH to connect to the allocated compute node, for example gl3121.
3. In that compute-node VS Code session, activate mge-organoid-python.
4. Export PROJECT_ROOT.
5. Reopen the notebook from that compute-node session.
6. Select mge-organoid-python (Python 3.11.15).
7. Rerun the diagnostic cell.
```

Alternative fix: connect VS Code to a Jupyter server started manually on the
compute node. From the compute-node terminal:

```bash
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
conda activate mge-organoid-python
export PROJECT_ROOT=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder
jupyter lab --no-browser --ip=0.0.0.0 --port=8899
```

Keep that terminal open. Jupyter prints a URL like:

```text
http://127.0.0.1:8899/lab?token=<token>
```

Do not paste the `127.0.0.1` URL into VS Code when VS Code is connected to a
login node. `127.0.0.1` means "this same machine", so from the login-node VS
Code session it points at the login node, not the compute node.

Replace `127.0.0.1` with the allocated compute node hostname. For example, if
the allocation is on `gl3121`, paste:

```text
http://gl3121.arc-ts.umich.edu:8899/lab?token=<token>
```

Use the same token printed by Jupyter.

In VS Code, connect the notebook to that existing server URL. If VS Code prompts
you to pick a Python environment instead of asking for a Jupyter server URL, you
are in the kernel/interpreter picker, not the existing-server picker. Cancel
that prompt.

From the notebook UI, use:

```text
Select Kernel -> Select Another Kernel -> Existing Jupyter Server -> Enter the URL of the running Jupyter server
```

Paste the full URL printed by the compute-node `jupyter lab` command.

If `Existing Jupyter Server` is not available in the notebook kernel selector,
open the command palette and search for one of these commands:

```text
Jupyter: Specify Jupyter Server for Connections
Jupyter: Select Jupyter Server
```

The correct path must ask for a server URL. If it only asks for a Python
environment, cancel and choose the existing-server path instead.

If a cached `.h5ad` already exists and is newer than the source `.rds`, the
cell should say it is using the existing cached H5AD instead of converting.

If the notebook was already open before this progress logging was added,
restart the notebook kernel and rerun the setup/import cells so it reloads the
updated `mge_organoid_python` module.

If the first notebook cell hangs after a previous kernel crash, the cell itself
is not doing conversion work. It only resolves the repo path and adds
`python_notebooks/src` to `sys.path`. A hang there usually means VS Code is
still attached to a wedged kernel process or has a stale notebook tab after the
file changed on disk.

The notebook file itself should have no saved running state. You can verify from
the repo root:

```bash
python3 - <<'PY'
import json
nb = json.load(open("python_notebooks/notebooks/01_seurat_to_anndata.ipynb"))
for i, cell in enumerate(nb["cells"][:5]):
    print(i, cell.get("cell_type"), "execution_count=", cell.get("execution_count"), "outputs=", len(cell.get("outputs", [])))
PY
```

Expected output for a clean notebook:

```text
0 markdown execution_count= None outputs= 0
1 markdown execution_count= None outputs= 0
2 code execution_count= None outputs= 0
3 code execution_count= None outputs= 0
4 code execution_count= None outputs= 0
```

If a cell starts running immediately when you open the notebook, check for stale
kernel processes:

```bash
ps -u "$USER" -f | rg 'jupyter|ipykernel|python.*kernel|ipython' || true
```

If you see an old `ipykernel_launcher` process for this notebook, shut down the
kernel from VS Code first. If the UI cannot stop it, terminate that specific
kernel PID from the terminal with `kill <PID>`, then reopen the notebook and
select `mge-organoid-python (Python 3.11.15)` again.

Only kill a PID that appears in `ps` from your current compute-node terminal. If
`kill <PID>` says `No such process`, that PID is not active in that terminal
context anymore; rerun the `ps` command and use the current PID, if any.

If `ps` prints no Jupyter/kernel process but the first cell appears to start
running immediately when the notebook opens, VS Code is showing stale notebook UI
state. The notebook file is not actually running. Use:

```text
Developer: Reload Window
```

Then reopen the notebook, select `mge-organoid-python (Python 3.11.15)`, and
run only the first code cell.

Recovery steps:

```text
1. Interrupt the notebook cell if it is still running.
2. Restart Kernel.
3. Close the notebook tab.
4. Reopen python_notebooks/notebooks/01_seurat_to_anndata.ipynb.
5. Select mge-organoid-python (Python 3.11.15).
6. Run only the first code cell again.
```

If it still hangs, reload VS Code:

```text
Developer: Reload Window
```

Then reconnect to the compute node if needed, reopen the notebook, select the
same kernel, and run only the first code cell.

Before reopening the notebook, you can confirm the kernel works from the
compute-node terminal:

```bash
python - <<'PY'
import sys
from pathlib import Path
print(sys.executable)
print(Path.cwd().resolve())
PY
```

Expected output includes:

```text
/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/python
/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
```

After the Shi smoke test succeeds, use the optional Varela cell. It is
intentionally written with commented lines:

```python
# varela_div30_path = converter.convert_file(studies_by_id["varela_div30"])
# print("Varela DIV30 cached at:", varela_div30_path)

# varela_div90_path = converter.convert_file(studies_by_id["varela_div90"])
# print("Varela DIV90 cached at:", varela_div90_path)
```

Uncomment and run one Varela conversion at a time. These use `convert_file`,
which writes `.h5ad` files but does not load the large `.h5ad` outputs back into
memory. This is safer for large objects.

Expected converted files:

```text
$PROJECT_ROOT/results/python_anndata/shi_2019_paper_qc.h5ad
$PROJECT_ROOT/results/python_anndata/varela_div30.h5ad
$PROJECT_ROOT/results/python_anndata/varela_div90.h5ad
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
