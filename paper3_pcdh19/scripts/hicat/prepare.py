"""Select the technical pilot from approved raw counts in the modern environment.

Only registered sample identity determines the balanced random subset.
Genotype, sex, previous clusters and provisional annotations never select or
fit cells. Step 06 PCA/UMAP is copied only as fixed display coordinates.
"""
from pathlib import Path
import argparse
import json
import h5py
import anndata as ad
import numpy as np
import pandas as pd
from .provenance import sha256, write_json


class PilotInputPreparer:
    """Create a small raw-count input while binding it to approved source bytes."""

    def __init__(self, run_dir, project_root):
        """Load frozen config and resolve exact Step 02 and Step 06 paths."""
        self.run_dir = Path(run_dir)
        self.project = Path(project_root)
        self.cfg = json.loads((self.run_dir/'config/hicat_pilot.json').read_text())

    def run(self):
        """Validate source identity, sample rows, and write a minimal pilot H5AD.

        Returns
        -------
        pathlib.Path
            Raw pilot input with fixed Step 06 display coordinates.

        Notes
        -----
        This reads the approved source and does not rewrite its metadata or
        status. All genes are kept. The sample-row ledger records exact IDs,
        source positions and the selection seed. No normalized matrix exists
        until the pilot's normalization phase.
        """
        cfg = self.cfg
        base = self.project/'results/primary_processing'
        source_run = base/'02_qc_filtering'/cfg['step02_run_id']
        source = source_run/'objects/pcdh19_step02_qc_filtered.h5ad'
        status = pd.read_csv(source_run/'STEP_STATUS.tsv', sep='\t')
        if len(status) != 1 or status.iloc[0]['status'] != 'APPROVED' or status.iloc[0]['run_id'] != cfg['step02_run_id']:
            raise ValueError('Pilot requires the exact approved Step 02 input')
        observed_hash = sha256(source)
        if observed_hash != cfg['step02_sha256']:
            raise ValueError('Approved Step 02 checksum mismatch')
        obj = ad.read_h5ad(source, backed='r')
        try:
            if obj.shape != (cfg['expected_cells'], cfg['expected_genes']):
                raise ValueError('Source shape mismatch')
            samples = obj.obs['technical_sample_id'].astype(str)
            if samples.nunique() != cfg['expected_samples']:
                raise ValueError('Sample count mismatch')
            rng = np.random.default_rng(cfg['selection_seed'])
            selected = []
            for sample in sorted(samples.unique()):
                eligible = np.flatnonzero(samples.to_numpy() == sample)
                selected.extend(rng.choice(eligible, min(cfg['pilot_cells_per_sample'],len(eligible)), replace=False))
            selected = np.sort(selected)
            fields = ['technical_sample_id','submitted_sample_name','genotype','sex','design_group',
                      'total_counts','n_genes_by_counts','pct_counts_mt']
            obs = obj.obs.iloc[selected][fields].copy()
            for field in fields[:5]:
                obs[field] = obs[field].astype(str)
            # Ensembl IDs are the AnnData index, not necessarily a separate column.
            var = pd.DataFrame({'gene_symbol': obj.var['gene_symbol'].astype(str)},
                               index=obj.var_names.copy())
            raw = obj[selected,:].X.copy()
            pilot = ad.AnnData(X=raw, obs=obs, var=var)
            ledger = pd.DataFrame({'cell_id': obs.index, 'source_row': selected,
                                   'sample': samples.iloc[selected].to_numpy(),
                                   'selection_seed': cfg['selection_seed']})
        finally:
            obj.file.close()
        display = base/'06_technical_sample_batch_diagnostics'/cfg['step06_display_run_id']/'objects/pcdh19_step06_unintegrated_diagnostics.h5ad'
        # Verify ID alignment rather than assuming matching row order across files.
        with h5py.File(display,'r') as handle:
            display_ids = handle['obs'][handle['obs'].attrs['_index']].asstr()[selected]
            if not np.array_equal(display_ids, pilot.obs_names.to_numpy()):
                raise ValueError('Step 06 display coordinates do not align with pilot cell IDs')
            pilot.obsm['X_umap_step06_display'] = handle['obsm/X_umap'][selected]
            pilot.obsm['X_pca_step06_display'] = handle['obsm/X_pca'][selected]
        output = self.run_dir/'inputs/pilot_raw_counts.h5ad'
        pilot.write_h5ad(output, compression='lzf')
        ledger.to_csv(self.run_dir/'inputs/pilot_cell_selection.tsv', sep='\t', index=False)
        provenance = dict(source=str(source), source_sha256=observed_hash, source_run=cfg['step02_run_id'],
                          display_source=str(display), display_sha256=sha256(display),
                          input_pilot_sha256=sha256(output), shape=list(pilot.shape),
                          selection='fixed random sample within each registered sample; all genes retained',
                          annotation_fields_used=False, sample_counts=ledger['sample'].value_counts().to_dict())
        write_json(self.run_dir/'inputs/input_identity.json',provenance)
        print(json.dumps(provenance),flush=True)
        return output


def main():
    """Parse run/project paths and prepare only the explicitly bounded pilot."""
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir',type=Path,required=True)
    parser.add_argument('--project-root',type=Path,required=True)
    args=parser.parse_args()
    PilotInputPreparer(args.run_dir,args.project_root).run()


if __name__ == '__main__':
    main()
