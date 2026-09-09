"""Run/resume the crude-first DIV90 Hypergate workflow from the existing assets.

Use the mge-organoid-python environment. Completed R checkpoints are reused;
incomplete checkpoints fail explicitly rather than being overwritten.
"""
from pathlib import Path
import os
import subprocess
import sys
import urllib.request

import pandas as pd
from div90_hypergate_paths import PROJECT_ROOT, REPO_ROOT

ROOT = REPO_ROOT
SCRIPTS = ROOT / "python_notebooks/scripts"
sys.path.insert(0, str(ROOT / "python_notebooks/src"))
from prepare_div90_hypergate import CORE, make_jobs

OUT = PROJECT_ROOT / "results/div90_hypergate_sst_pv"
SOURCE = PROJECT_ROOT
RSCRIPT = Path(sys.executable).parent / "Rscript"


def main():
    (OUT / "provenance").mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, PYTHONPATH=str(ROOT / "python_notebooks/src"),
               OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    subprocess.run([str(RSCRIPT), "-e", 'if(!requireNamespace("hypergate",quietly=TRUE)) install.packages("hypergate",repos="https://cloud.r-project.org")'], check=True, env=env)
    if not (OUT / "core_input.tsv.gz").exists():
        data = pd.read_csv(PROJECT_ROOT / "final_figures/div90_guidance_with_sst/tables/joined_with_sst.tsv.gz", sep="\t")
        data = data.loc[data.set_id.eq("cortical_only") & data.expression_available].copy()
        membership = pd.read_csv(SOURCE / "final_figures/fig_div90_loupe_recluster_annotations_v1/tables/div90_loupe_recluster_membership.tsv.gz", sep="\t",
                                 dtype={"cell_id": str, "cluster_id": str, "orig.ident": str})
        membership = membership.loc[membership.set_id.eq("cortical_only"), ["cell_id", "cluster_id", "orig.ident"]]
        data = data.merge(membership, on="cell_id", validate="one_to_one").rename(columns={"cluster_id": "cluster", "orig.ident": "sample"})
        data = data.loc[data.LHX6.gt(0) & data.ERBB4.gt(0) & ~data.cluster.astype(str).isin(["6", "7"])].copy()
        data["crude_label"] = data.SST.gt(0).map({True: "SST-like", False: "non-SST/PV-candidate"})
        assert len(data) == 4768 and data.SST.gt(0).sum() == 2919
        data.to_csv(OUT / "core_input.tsv.gz", sep="\t", index=False)
    if not (OUT / "core_jobs.tsv").exists():
        make_jobs(CORE, "crude", "core").to_csv(OUT / "core_jobs.tsv", sep="\t", index=False)
    if not (OUT / "core_fits.tsv").exists():
        subprocess.run([str(RSCRIPT), str(SCRIPTS / "run_div90_hypergate.R"), str(OUT / "core_input.tsv.gz"),
                        str(OUT / "core_jobs.tsv"), str(OUT / "core_fits.tsv")], check=True, env=env)
    assert len(pd.read_csv(OUT / "core_jobs.tsv", sep="\t")) == len(pd.read_csv(OUT / "core_fits.tsv", sep="\t", escapechar="\\"))
    surfaceome = OUT / "provenance/CSPA_human_surfaceome.xlsx"
    if not surfaceome.exists():
        url = "https://journals.plos.org/plosone/article/file?type=supplementary&id=10.1371/journal.pone.0121314.s003"
        with urllib.request.urlopen(url, timeout=60) as response:
            surfaceome.write_bytes(response.read())
    for script in ["prepare_div90_hypergate.py", "run_div90_hypergate_expanded.py", "summarize_div90_hypergate.py"]:
        subprocess.run([sys.executable, str(SCRIPTS / script)], check=True, env=env)


if __name__ == "__main__":
    main()
