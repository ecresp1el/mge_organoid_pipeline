# Python Notebook Entry Point

This directory is the Python notebook entry point for working with cached
AnnData `.h5ad` files made from the project Seurat `.rds` objects.

The notebook workflow is load-only:

```text
open notebook -> load cached .h5ad files -> inspect/plot/analyze AnnData
```

Do not remake `.h5ad` files from the notebook during routine work. If a source
Seurat object changes, rebuild the cache with the Slurm batch job in
`Rebuild Cached H5AD Files`.

Runtime data stays outside the git repo under:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder
```

## Confirmed Project State

```text
repo:         /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
conda base:   /home/elcrespo/miniconda3
env:          mge-organoid-python
python:       Python 3.11.15
account:      parent0
PROJECT_ROOT: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder
```

Cached AnnData files were created successfully by Slurm job `50295259`:

```text
shi_2019_paper_qc.h5ad  n_obs=56,136  n_vars=21,191  has_umap=True
varela_div30.h5ad       n_obs=90,631  n_vars=18,082  has_umap=True
varela_div90.h5ad       n_obs=22,338  n_vars=18,082  has_umap=True
```

Expected cached files:

```text
$PROJECT_ROOT/results/python_anndata/shi_2019_paper_qc.h5ad
$PROJECT_ROOT/results/python_anndata/varela_div30.h5ad
$PROJECT_ROOT/results/python_anndata/varela_div90.h5ad
```

## Notebook 02 Monocle3 Batch Stage

Notebook `02_mgeo_rgc_ipc_progression.ipynb` expects Monocle3 pseudotime from
the Slurm stage below. The batch job exports the focused MGEO RGC/IPC subset,
runs Monocle3 in R, and writes outputs to:

```text
$PROJECT_ROOT/results/mgeo_rgc_ipc_monocle3/
```

Submit:

```bash
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
mkdir -p "$PROJECT_ROOT/jobs" "$PROJECT_ROOT/logs"
cp slurm_templates/09_mgeo_rgc_ipc_monocle3.sbatch.template \
  "$PROJECT_ROOT/jobs/09_mgeo_rgc_ipc_monocle3.sbatch"
sbatch "$PROJECT_ROOT/jobs/09_mgeo_rgc_ipc_monocle3.sbatch"
```

Key outputs:

```text
$PROJECT_ROOT/results/mgeo_rgc_ipc_monocle3/inputs/mgeo_rgc_ipc_input_manifest.tsv
$PROJECT_ROOT/results/mgeo_rgc_ipc_monocle3/mgeo_rgc_ipc_monocle3_pseudotime.csv
$PROJECT_ROOT/results/mgeo_rgc_ipc_monocle3/mgeo_rgc_ipc_monocle3_cds.rds
$PROJECT_ROOT/results/mgeo_rgc_ipc_monocle3/mgeo_rgc_ipc_monocle3_summary.tsv
```

The Slurm stage uses conda only for the Python AnnData export, then switches to
the Great Lakes Monocle3 module for R:

```text
module load Bioinformatics
module load Rmonocle3/1.3.7
```

After the job finishes, rerun notebook 02 through section 9. Section 9 loads
`mgeo_rgc_ipc_monocle3_pseudotime.csv` back into `mgeo.obs`.

For live Monocle3 troubleshooting in VS Code, request an interactive compute
node with:

```bash
python_notebooks/scripts/start_vscode_r_monocle_compute_job.sh
```

After Slurm grants the allocation, connect VS Code to the compute node and run
the commands printed by the script. Then open:

```text
exploration/09_mgeo_rgc_ipc_monocle3_interactive.R
```

That file uses `# %%` sections so you can run the Monocle3 workflow block by
block, inspect plots, and copy working choices back into the batch script.

## Fresh Login Setup

Use these steps from a fresh Great Lakes login-node session.

### 1. Open The Repo

```bash
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
```

Expected:

```bash
pwd
```

prints:

```text
/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
```

### 2. Confirm Conda

```bash
conda --version
type -a conda
conda info --base
```

Expected for this project:

```text
conda 25.3.1
/home/elcrespo/miniconda3
```

It is okay if `type -a conda` says `conda is a function`. The important check
is that `conda info --base` points to:

```text
/home/elcrespo/miniconda3
```

If `conda` is not found, stop and fix conda access before continuing.

### 3. Confirm Or Update The Environment

Check whether the environment exists:

```bash
conda env list
```

Expected:

```text
mge-organoid-python    /home/elcrespo/miniconda3/envs/mge-organoid-python
```

If the env exists, update it:

```bash
conda env update -n mge-organoid-python -f python_notebooks/environment.yml --prune
```

Expected:

```text
Solving environment: done
```

It is okay if conda prints a warning that a newer conda version exists. Do not
update base conda for this workflow.

If the env does not exist, create it:

```bash
conda env create -f python_notebooks/environment.yml
```

### 4. Activate The Environment

```bash
conda activate mge-organoid-python
```

Validate:

```bash
echo "$CONDA_DEFAULT_ENV"
which python
python --version
```

Expected:

```text
mge-organoid-python
~/miniconda3/envs/mge-organoid-python/bin/python
Python 3.11.15
```

### 5. Confirm Monocle3 For Notebook 02

Notebook `02_mgeo_rgc_ipc_progression.ipynb` runs Monocle3 from Python via
`Rscript`, so Monocle3 must be installed in the same active conda environment:

```bash
which Rscript
Rscript -e 'cat(requireNamespace("monocle3", quietly = TRUE), "\n")'
```

Expected:

```text
TRUE
```

If the check prints `FALSE`, install Monocle3 into this environment before
rerunning the Monocle3 notebook cell:

```bash
Rscript -e 'install.packages(c("remotes", "BiocManager"), repos = "https://cloud.r-project.org")'
Rscript -e 'remotes::install_github("cole-trapnell-lab/monocle3")'
```

### 5. Register The Notebook Kernel

```bash
python -m ipykernel install --user \
  --name mge-organoid-python \
  --display-name "Python (mge-organoid-python)"
```

Expected:

```text
Installed kernelspec mge-organoid-python in /home/elcrespo/.local/share/jupyter/kernels/mge-organoid-python
```

This is what lets VS Code show:

```text
Python (mge-organoid-python)
```

### 6. Set PROJECT_ROOT

```bash
export PROJECT_ROOT=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder
echo "$PROJECT_ROOT"
```

Expected:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder
```

### 7. Confirm Cached H5AD Files Exist

```bash
ls -lh \
  "$PROJECT_ROOT/results/python_anndata/shi_2019_paper_qc.h5ad" \
  "$PROJECT_ROOT/results/python_anndata/varela_div30.h5ad" \
  "$PROJECT_ROOT/results/python_anndata/varela_div90.h5ad"
```

Expected:

```text
shi_2019_paper_qc.h5ad
varela_div30.h5ad
varela_div90.h5ad
```

If any file is missing, go to `Rebuild Cached H5AD Files`.

### 8. Confirm Python Can Load The Cached AnnData Files

```bash
PYTHONPATH="$PWD/python_notebooks/src" python - <<'PY'
from mge_organoid_python import default_studies, load_cached_anndatas, resolve_project_root

project_root = resolve_project_root()
adatas, reports = load_cached_anndatas(default_studies(), project_root=project_root, backed="r")
for report in reports:
    print(report.study_id, report.n_obs, report.n_vars, report.has_umap)
for adata in adatas.values():
    adata.file.close()
PY
```

Expected:

```text
shi_2019_paper_qc 56136 21191 True
varela_div30 90631 18082 True
varela_div90 22338 18082 True
```

If this works, the cached AnnData files are valid.

### 9. Request Notebook Resources

Use the login node for setup checks only. For actually working in the notebook,
request a Slurm compute node so the notebook kernel has allocated RAM/CPUs.

From the login-node VS Code terminal:

```bash
python_notebooks/scripts/start_vscode_compute_job.sh
```

Default request:

```text
account:   parent0
partition: standard
memory:    128GB
CPUs:      4
time:      04:00:00
```

Expected Slurm output:

```text
salloc: Granted job allocation <job_id>
salloc: Nodes <node_name> are ready for job
```

Keep this terminal open. The allocation ends if this terminal exits.

### 10. Connect VS Code To The Compute Node

After Slurm grants the allocation, connect VS Code Remote SSH to the allocated
compute node. The node name is the `<node_name>` from the `salloc` output.

In the compute-node VS Code terminal:

```bash
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
conda activate mge-organoid-python
export PROJECT_ROOT=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder
hostname
echo "$CONDA_DEFAULT_ENV"
which python
echo "$PROJECT_ROOT"
```

Expected:

```text
<compute-node>.arc-ts.umich.edu
mge-organoid-python
~/miniconda3/envs/mge-organoid-python/bin/python
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder
```

Confirmed example:

```text
hostname                  -> gl3088.arc-ts.umich.edu
echo "$CONDA_DEFAULT_ENV" -> mge-organoid-python
which python              -> ~/miniconda3/envs/mge-organoid-python/bin/python
echo "$PROJECT_ROOT"      -> /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder
```

If `hostname` starts with `gl-login`, stop. The notebook is still attached to a
login node, not to your allocated compute node.

### 11. Start Jupyter On The Compute Node

In the compute-node terminal, start Jupyter with a URL that VS Code can connect
to:

```bash
jupyter lab --no-browser --ip=0.0.0.0 --port=8899
```

Keep this terminal open. Jupyter will print URLs like:

```text
http://<compute-node>.arc-ts.umich.edu:8899/lab?token=<token>
http://127.0.0.1:8899/lab?token=<token>
```

Use the URL with the compute-node hostname, not the `127.0.0.1` URL. For
example:

```text
http://gl3088.arc-ts.umich.edu:8899/lab?token=<token>
```

### 12. Open The Notebook And Connect It To That Jupyter Server

In VS Code, open the notebook first:

```text
python_notebooks/notebooks/01_seurat_to_anndata.ipynb
```

Then use the notebook kernel picker in the upper right of the notebook.

Click:

```text
Select Kernel
```

Then choose:

```text
Select Another Kernel
```

Then choose:

```text
Existing Jupyter Server
```

Then choose:

```text
Enter the URL of the running Jupyter Server
```

Paste the full compute-node URL printed by Jupyter. Use the `gl####` URL, not
the `127.0.0.1` URL:

```text
http://gl3088.arc-ts.umich.edu:8899/lab?token=<token>
```

After you paste the URL, VS Code may show the server host name, for example:

```text
gl3088.arc-ts.umich.edu
```

Press Enter to accept it.

Then select the kernel from that server:

```text
Python (mge-organoid-python)
```

Expected terminal evidence that VS Code reached the compute-node Jupyter server:

```text
302 GET / (@...)
302 GET /lab? (@...)
```

If VS Code asks you to pick a local Python environment, Conda env, or interpreter
before asking for a Jupyter server URL, cancel that prompt. You are in the wrong
picker. Go back to the notebook's `Select Kernel` button and choose:

```text
Select Another Kernel -> Existing Jupyter Server -> Enter the URL of the running Jupyter Server
```

Run the notebook cells in order. The notebook:

```text
1. adds python_notebooks/src to sys.path
2. imports the cached AnnData loader
3. prints runtime diagnostics
4. checks source and cached output paths
5. loads cached .h5ad files with backed="r"
6. displays validation reports
7. plots UMAPs from adata.obsm["X_umap"]
```

The notebook does not convert Seurat objects and does not overwrite `.h5ad`
files.

## Login Node Vs Compute Node

Light cached loading with `backed="r"` can work from a login-node notebook, but
the standard workflow for this project is to request a Slurm compute allocation
before opening the notebook for real analysis.

Use a Slurm compute node when you plan to do memory-heavy AnnData work, model
fitting, integration, large in-memory transformations, or anything that will
load large matrices into memory.

## Rebuild Cached H5AD Files

Only do this when a cached `.h5ad` file is missing or a source Seurat `.rds`
file has changed.

Submit the Slurm rebuild job:

```bash
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
cp slurm_templates/08_convert_python_anndata.sbatch.template \
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/08_convert_python_anndata.sbatch

sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/08_convert_python_anndata.sbatch
```

Default rebuild behavior:

```text
account:   parent0
partition: standard
memory:    160G
CPUs:      4
time:      12:00:00
studies:   all
overwrite: false
```

Because `overwrite=false`, existing current `.h5ad` files are skipped.

Monitor:

```bash
squeue -u elcrespo
tail -f /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/08_convert_python_anndata_<jobid>.log
```

Successful completion should show:

```text
All file conversions complete
shi_2019_paper_qc ... exists=True
varela_div30 ... exists=True
varela_div90 ... exists=True
```

## Canonical Inputs

```text
Shi 2019 paper QC:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_2019_paper_qc/shi_2019_seurat.rds

Varela DIV30:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/varela_this_paper/varela_this_paper_seurat.rds

Varela DIV90:
/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/umap_props_output/clustered_day90_with_cluster_names_2.rds
```

## Python API Used By The Notebook

```python
from mge_organoid_python import (
    default_studies,
    load_cached_anndatas,
    missing_cached_h5ads,
    resolve_project_root,
)

project_root = resolve_project_root()
studies = default_studies()
missing = missing_cached_h5ads(studies, project_root=project_root)
if missing:
    raise FileNotFoundError(missing)

adatas, reports = load_cached_anndatas(studies, project_root=project_root, backed="r")
```
