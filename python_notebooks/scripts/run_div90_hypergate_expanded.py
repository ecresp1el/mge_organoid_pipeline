"""Run the two expanded searches sequentially, using four independent R workers."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import os
import subprocess
import pandas as pd
from div90_hypergate_paths import PROJECT_ROOT, REPO_ROOT

ROOT = REPO_ROOT
OUT = PROJECT_ROOT / "results/div90_hypergate_sst_pv"
RSCRIPT = "/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/Rscript"


def run_part(args):
    jobs, fits, log = args
    if fits.exists():
        n_expected = len(pd.read_csv(jobs, sep="\t"))
        n_observed = len(pd.read_csv(fits, sep="\t", escapechar="\\"))
        if n_expected == n_observed:
            return str(fits)
        raise RuntimeError(f"Incomplete checkpoint must be inspected: {fits}")
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    with log.open("w") as handle:
        subprocess.run([RSCRIPT, str(ROOT / "python_notebooks/scripts/run_div90_hypergate.R"),
                        str(OUT / "hypergate_input.tsv.gz"), str(jobs), str(fits)],
                       env=env, stdout=handle, stderr=subprocess.STDOUT, check=True)
    return str(fits)


def main():
    for definition in ["crude", "module"]:
        jobs = pd.read_csv(OUT / f"expanded_{definition}_jobs.tsv", sep="\t")
        tasks = []
        for worker in range(4):
            path = OUT / f"expanded_{definition}_jobs_{worker}.tsv"
            jobs.iloc[worker::4].to_csv(path, sep="\t", index=False)
            tasks.append((path, OUT / f"expanded_{definition}_fits_{worker}.tsv",
                          OUT / f"provenance/expanded_{definition}_{worker}.log"))
        print(f"Starting {definition}: {len(jobs)} Hypergate fits", flush=True)
        with ThreadPoolExecutor(max_workers=4) as pool:
            for result in pool.map(run_part, tasks):
                print(f"Finished {result}", flush=True)


if __name__ == "__main__":
    main()
