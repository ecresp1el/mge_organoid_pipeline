"""Execute Step 07 on the frozen E14.5 dissected mouse MGE HiCAT pilot.

The workflow verifies and reuses original counts, normalized expression,
partitions and display coordinates. It adds independent diagnostics and one
fixed-parent fine sensitivity partition. Full-data operations only calculate
canonical scores and display existing UMAP coordinates. Successful publication
is atomic and IN_REVIEW; it never authorizes the next scientific checkpoint.
"""
from pathlib import Path
import argparse
import gc
import json
import logging
import os
import shutil
import sys
import importlib.metadata
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from .provenance import Progress, write_json, sha256, manifest, json_value
from .validation_metrics import partition_overlap, sample_composition, phase_matched_identity
from .validation_sensitivity import FixedParentSensitivity, resolve_sensitivity_config
from .validation_programs import CanonicalProgramScorer
from .validation_full import FullDataProjector
from .validation_report import create_report


class ValidationWorkflow:
    """Own one new review package from input verification through saved checks."""

    def __init__(self, run_dir):
        """Load frozen config and create isolated staging/progress locations."""
        self.root = Path(run_dir)
        self.cfg = json.loads((self.root/'config/hicat_validation.json').read_text())
        self.cfg['programs_path'] = str(self.root/'config/hicat_validation_programs.json')
        self.cfg['hypotheses_path'] = str(self.root/'config/hicat_validation_hypotheses.json')
        self.cfg['published_output_dir'] = str(self.root/'outputs')
        self.workspace = self.root/'staging'
        self.output = self.workspace/'outputs'
        self.output.mkdir(parents=True)
        for name in ['figures', 'tables', 'marker_programs', 'stability',
                     'sample_composition', 'cell_cycle', 'full_data_projection',
                     'parameter_sensitivity', 'annotation_review', 'validation']:
            (self.output/name).mkdir()
        self.progress = Progress(self.root/'provenance')
        self.checks = []

    def check(self, name, condition, detail=''):
        """Record an explicit contract check; fail immediately on a violation."""
        self.checks.append(dict(check=name, status='PASS' if condition else 'FAIL', detail=detail))
        pd.DataFrame(self.checks).to_csv(self.output/'validation/checks.tsv', sep='\t', index=False)
        if not condition:
            raise ValueError(name+': '+detail)

    def load_input(self):
        """Verify exact source identity and expression before any fit or scoring."""
        cfg = self.cfg
        self.check('authorized_review_scope', all(cfg[key] is False for key in
                   ['full_data_clustering', 'batch_correction', 'regression', 'lock_annotations']))
        self.check('correct_biological_context', cfg['biological_context']['developmental_age'] == 'E14.5'
                   and cfg['biological_context']['organoid'] is False)
        source = Path(cfg['baseline_run'])/'pcdh19_hicat_coarse_fine.h5ad'
        self.check('pilot_sha256', sha256(source) == cfg['pilot_sha256'])
        self.check('full_source_sha256', sha256(cfg['full_data_source']) == cfg['full_data_sha256'])
        pilot = ad.read_h5ad(source)
        self.check('authorized_pilot_shape', pilot.shape == (cfg['pilot_cells'], cfg['genes']))
        self.check('unique_cell_and_gene_ids', pilot.obs_names.is_unique and pilot.var_names.is_unique)
        counts = pilot.obs.technical_sample_id.value_counts()
        self.check('exact_balanced_selection', len(counts) == cfg['samples'] and counts.eq(cfg['cells_per_sample']).all())
        self.check('nonnegative_sparse_raw_integer', sparse.issparse(pilot.X) and
                   pilot.X.dtype.kind in 'iu' and bool((pilot.X.data >= 0).all()))
        totals = np.asarray(pilot.X.sum(axis=1)).ravel()
        self.check('positive_cell_totals', bool((totals > 0).all()))
        expected = pilot.X.astype(np.float64).multiply((1e6/totals)[:, None]).tocsr()
        expected.data = np.log1p(expected.data)
        difference = expected - pilot.layers['log1p_cpm']
        self.check('saved_normalization_formula', difference.nnz == 0 or
                   float(np.max(np.abs(difference.data))) < 1e-10)
        del expected, difference
        baseline_config = json.loads(pilot.uns['hicat']['resolved_config_json'])
        frozen_baseline = json.loads((self.root/'config/hicat_coarse_fine.json').read_text())
        self.check('baseline_config_matches_original', baseline_config == frozen_baseline)
        self.check('baseline_counts', pilot.obs.hicat_coarse_baseline.nunique() == 4 and
                   pilot.obs.hicat_fine_baseline.nunique() == 38)
        self.check('repeat_counts', pilot.obs.hicat_coarse_seed_repeat.nunique() == 5 and
                   pilot.obs.hicat_fine_seed_repeat.nunique() == 39)
        write_json(self.root/'inputs/input_identity.json', dict(shape=list(pilot.shape),
            input_pilot=str(source), input_pilot_sha256=cfg['pilot_sha256'],
            full_display_source=cfg['full_data_source'], full_display_sha256=cfg['full_data_sha256'],
            biological_context=cfg['biological_context']))
        pilot.obs.to_csv(self.root/'inputs/original_pilot_assignments.tsv.gz', sep='\t', index_label='cell_id')
        # Copy plotted marker sources; link complete historical all-gene statistics.
        shutil.copy2(Path(cfg['baseline_run'])/'top20_genes_all_contexts.tsv',
                     self.output/'tables/top20_genes_all_contexts.tsv')
        write_json(self.output/'tables/existing_marker_statistics.json', dict(
            root=str(Path(cfg['baseline_run'])/'markers'),
            note='Complete all-gene statistics are preserved in the immutable input run; no top20 padding.'))
        return pilot, baseline_config

    def save_metrics(self, group, name, result):
        """Save each returned diagnostic table and its exact thresholds/summary."""
        target = self.output/group
        for key, value in result.items():
            if isinstance(value, (pd.Series, pd.DataFrame)):
                value.to_csv(target/(name+'_'+key+'.tsv'), sep='\t')
            elif isinstance(value, dict):
                write_json(target/(name+'_'+key+'.json'), value)

    def auxiliary_tables(self, pilot, program):
        """Save cluster score summaries, phase fractions, QC and descriptive eta².

        Eta² measures score variance associated with existing cluster labels;
        it does not identify the cause of a cluster or correct expression.
        Absolute means/medians and program-positive fractions are retained;
        any z scoring for heatmaps is a display-only operation in reporting.
        """
        scores = program['scores']
        pilot.obs.to_csv(self.output/'tables/pilot_cell_metadata.tsv.gz', sep='\t', index_label='cell_id')
        axes = scores.copy()
        for name in ['total_counts', 'n_genes_by_counts', 'pct_counts_mt']:
            axes['QC_'+name] = pd.to_numeric(pilot.obs[name], errors='coerce')
        records = []
        for level in ['coarse', 'fine']:
            labels = pilot.obs['hicat_'+level+'_baseline'].astype(str)
            pd.crosstab(labels, program['phase'], normalize='index').to_csv(
                self.output/'cell_cycle'/('phase_composition_'+level+'.tsv'), sep='\t')
            for statistic, frame in [('mean', scores.groupby(labels).mean()),
                                     ('median', scores.groupby(labels).median()),
                                     ('positive_fraction', program['positive'].groupby(labels).mean())]:
                frame.to_csv(self.output/'marker_programs'/(level+'_program_'+statistic+'.tsv'), sep='\t')
            labels.value_counts().sort_index().to_csv(self.output/'marker_programs'/(level+'_cluster_sizes.tsv'), sep='\t')
            for name in axes:
                values = axes[name]
                valid = values.notna()
                if not valid.any():
                    eta = np.nan
                else:
                    mean = values[valid].mean()
                    total = ((values[valid]-mean)**2).sum()
                    group_mean = values[valid].groupby(labels[valid]).mean()
                    group_size = values[valid].groupby(labels[valid]).size()
                    eta = float((((group_mean-mean)**2)*group_size).sum()/total) if total > 0 else 0.
                records.append(dict(level=level, axis=name, eta_squared=eta,
                                    interpretation='descriptive association; not causal or corrected'))
        pd.DataFrame(records).to_csv(self.output/'tables/axis_association.tsv', sep='\t', index=False)

    def reuse_completed_sensitivity(self, pilot, baseline_config):
        """Copy a hash-verified completed candidate from a failed reporting run.

        This recovery path avoids repeating the scientific comparison when a
        later display/I/O step failed. It requires identical pilot bytes,
        fitting code, upstream pin, complete engine controls and seed. Original
        run evidence is preserved, and the new package records its source.
        No missing parent, partial candidate or changed fit is accepted.
        """
        donor = Path(self.cfg['reuse_completed_sensitivity_run'])
        source = donor/'staging/outputs/parameter_sensitivity'
        donor_cfg = json.loads((donor/'config/hicat_validation.json').read_text())
        for key in ['pilot_sha256', 'baseline_run', 'pilot_cells', 'genes', 'seed',
                    'fine_changes', 'fixed_fine_score', 'fixed_baseline_parents', 'upstream_commit']:
            self.check('reused_sensitivity_'+key, donor_cfg[key] == self.cfg[key])
        resolved, _ = resolve_sensitivity_config(baseline_config, self.cfg)
        self.check('reused_sensitivity_complete_engine_config',
                   json.loads((source/'resolved_engine_config.json').read_text()) == resolved)
        for name in ['engine.py', 'validation_sensitivity.py']:
            self.check('reused_fitting_code_'+name,
                       sha256(donor/'code/hicat'/name) == sha256(self.root/'code/hicat'/name))
        identity = json.loads((donor/'inputs/input_identity.json').read_text())
        self.check('reused_candidate_input_identity', identity['input_pilot_sha256'] == self.cfg['pilot_sha256'])
        declared_path = donor/'provenance/completed_sensitivity_manifest.tsv'
        self.check('reuse_manifest_identity', sha256(declared_path) == self.cfg['reuse_sensitivity_manifest_sha256'])
        declared = pd.read_csv(declared_path, sep='\t')
        measured = manifest(source)
        self.check('reused_candidate_all_files_exact', measured.equals(declared))
        for parent in sorted(pilot.obs.hicat_coarse_baseline.astype(str).unique()):
            summary = json.loads((source/('fine_'+parent)/'candidate_summary.json').read_text())
            self.check('reused_parent_completed_'+parent,
                       summary['cells'] == int((pilot.obs.hicat_coarse_baseline.astype(str) == parent).sum()))
        target = self.output/'parameter_sensitivity'
        self.check('reuse_destination_empty', not any(target.iterdir()))
        shutil.copytree(source, target, dirs_exist_ok=True)
        alternative = pd.read_csv(target/'cell_assignments.tsv.gz', sep='\t', index_col=0).iloc[:, 0]
        parents = pd.read_csv(target/'parent_counts.tsv', sep='\t')
        self.check('reused_candidate_cell_order', alternative.index.equals(pilot.obs_names))
        self.check('reused_candidate_fixed_parents',
                   alternative.str.split('.').str[0].equals(pilot.obs.hicat_coarse_baseline.astype(str)))
        write_json(self.root/'provenance/reused_sensitivity.json', dict(
            donor_run=str(donor), candidate_assets=str(source), verified_files=len(declared),
            manifest_sha256=self.cfg['reuse_sensitivity_manifest_sha256'],
            reason='Recover completed controlled comparison after a downstream full-data metadata-reader failure.',
            scientific_comparison_refitted=False, baseline_labels_changed=False))
        return alternative, parents

    def save_checkpoint(self, pilot, original, program, report_summary, full_summary):
        """Write/round-trip the pilot object and inventory everything in uns."""
        pilot.uns['hicat_validation'] = dict(
            stage='02_hierarchy_validation', primary_processing_step='07',
            status='IN_REVIEW', run_id=self.root.name,
            tissue='dissected E14.5 mouse medial ganglionic eminence',
            annotation_locked=False, full_data_clustering=False,
            batch_correction=False, regression=False,
            input_pilot_sha256=self.cfg['pilot_sha256'],
            original_hicat_assets_root=self.cfg['baseline_run'],
            inherited_uns_context='uns.hicat describes the original expanded pilot; its relative external asset paths resolve under original_hicat_assets_root, not this validation package.',
            resolved_config_json=json.dumps(self.cfg, sort_keys=True),
            program_columns_json=json.dumps(program['scores'].columns.tolist()),
            program_model_asset='marker_programs/',
            full_display_assets='full_data_projection/',
            review_hypotheses_external='annotation_review/',
            original_labels='hicat_coarse_baseline,hicat_fine_baseline,hicat_coarse_seed_repeat,hicat_fine_seed_repeat',
            alternative_labels='hicat_fine_allen_reference; separate sensitivity partition, never adopted',
            saved_X='raw sparse integer counts; log1p_cpm retained unchanged in layers',
            readiness=str(report_summary.get('readiness', report_summary.get('decision', 'see summary.json'))))
        path = self.output/'pcdh19_hicat_hierarchy_validation.h5ad'
        pilot.write_h5ad(path, compression='lzf')
        reopened = ad.read_h5ad(path)
        self.check('saved_shape', reopened.shape == original.shape)
        self.check('saved_raw_counts_exact', (reopened.X != original.X).nnz == 0)
        self.check('saved_normalized_expression_exact',
                   (reopened.layers['log1p_cpm'] != original.layers['log1p_cpm']).nnz == 0)
        self.check('saved_cell_gene_order', reopened.obs_names.equals(original.obs_names) and reopened.var_names.equals(original.var_names))
        for column in ['hicat_coarse_baseline', 'hicat_fine_baseline',
                       'hicat_coarse_seed_repeat', 'hicat_fine_seed_repeat']:
            self.check('original_'+column+'_preserved',
                       reopened.obs[column].astype(str).equals(original.obs[column].astype(str)))
        self.check('original_umap_preserved', np.array_equal(reopened.obsm['X_umap_step06_display'], original.obsm['X_umap_step06_display']))
        self.check('candidate_nested_in_fixed_parents',
                   reopened.obs.hicat_fine_allen_reference.astype(str).str.split('.').str[0].equals(
                       reopened.obs.hicat_coarse_baseline.astype(str)))
        self.check('saved_program_scores', np.allclose(reopened.obsm['X_validation_program_scores'],
                   program['scores'].to_numpy(), equal_nan=True))
        self.check('no_locked_annotations', not reopened.uns['hicat_validation']['annotation_locked'])
        self.check('no_full_data_clustering', not reopened.uns['hicat_validation']['full_data_clustering'])
        write_json(self.output/'anndata_uns_inventory.json', reopened.uns)
        write_json(self.output/'anndata_slot_inventory.json', dict(
            shape=list(reopened.shape), X='raw integer counts', raw_present=reopened.raw is not None,
            layers=list(reopened.layers), obs=list(reopened.obs), var=list(reopened.var),
            obsm=list(reopened.obsm), obsp=list(reopened.obsp), uns=list(reopened.uns),
            external_full_data_sidecar=full_summary))
        del reopened

    def run(self):
        """Compute the authorized review, validate every output, publish IN_REVIEW."""
        logging.getLogger().setLevel(logging.INFO)
        logging.getLogger().addHandler(logging.FileHandler(self.root/'logs/upstream.log'))
        with self.progress.track('validation.verify_input'):
            pilot, baseline_cfg = self.load_input()
            original = pilot.copy()
        with self.progress.track('validation.canonical_programs', labels_used_for_fit=False):
            scorer = CanonicalProgramScorer(self.cfg['programs_path'])
            program = scorer.fit_transform(pilot, self.output/'marker_programs')
            self.auxiliary_tables(pilot, program)
        metrics = {}
        with self.progress.track('validation.existing_seed_and_sample_diagnostics'):
            for level in ['coarse', 'fine']:
                baseline = pilot.obs['hicat_'+level+'_baseline'].astype(str)
                repeat = pilot.obs['hicat_'+level+'_seed_repeat'].astype(str)
                metrics['seed_'+level] = partition_overlap(baseline, repeat, **self.cfg['stability_thresholds'])
                self.save_metrics('stability', level, metrics['seed_'+level])
                metrics['sample_'+level] = sample_composition(baseline, pilot.obs.technical_sample_id,
                                                             **self.cfg['sample_thresholds'])
                self.save_metrics('sample_composition', level, metrics['sample_'+level])
            identity = program['identity_scores'].dropna(axis=1, how='all')
            # Only developmental/regional programs participate, never cycle or QC.
            declared = json.loads(Path(self.cfg['programs_path']).read_text())['programs']
            identity_names = [item['name'] for item in declared if item['group'] in
                              ['development', 'maturation', 'region', 'non_neural'] and item['name'] in identity]
            self.check('identity_programs_available', len(identity_names) >= 5)
            self.check('regional_and_maturation_phase_comparisons',
                       all(key in identity_names for key in ['MGE', 'LGE', 'POA', 'later_neuronal_maturation']))
            metrics['phase_pairs'] = phase_matched_identity(
                pilot.obs.hicat_fine_baseline, pilot.obs.hicat_coarse_baseline,
                program['phase'], identity[identity_names], **self.cfg['phase_thresholds'])
            self.save_metrics('cell_cycle', 'matched_identity', metrics['phase_pairs'])
        with self.progress.track('validation.one_allen_comparison', fixed_coarse_parents=True):
            if self.cfg.get('reuse_completed_sensitivity_run'):
                alternative, parent_counts = self.reuse_completed_sensitivity(pilot, baseline_cfg)
            else:
                sensitivity = FixedParentSensitivity(self.output/'parameter_sensitivity', baseline_cfg, self.cfg, self.progress)
                alternative, parent_counts = sensitivity.run(pilot)
            pilot.obs['hicat_fine_allen_reference'] = pd.Categorical(alternative)
            metrics['allen_fine'] = partition_overlap(pilot.obs.hicat_fine_baseline, alternative,
                                                       **self.cfg['stability_thresholds'])
            metrics['parent_counts'] = parent_counts
            self.save_metrics('parameter_sensitivity', 'fine', metrics['allen_fine'])
        gc.collect()
        with self.progress.track('validation.full_data_marker_display', cells=self.cfg['full_cells'],
                                 clustering=False, existing_umap=True):
            full = FullDataProjector(self.cfg['full_data_source'], self.output/'full_data_projection',
                                     chunk_size=self.cfg['chunk_size']).run(pilot, scorer)
            self.check('full_display_cells', full['n_cells'] == self.cfg['full_cells'])
            self.check('full_display_pilot_cells', full['n_pilot_cells'] == self.cfg['pilot_cells'])
        with self.progress.track('validation.reports', fine_review_pages=38):
            report = create_report(self.workspace, pilot, program['scores'], program, metrics, self.cfg, full_context=full)
        with self.progress.track('validation.serialize_and_check'):
            # Reporting reads the staged sidecar; saved metadata points at its
            # final location after atomic publication.
            full['object_path'] = str(self.root/'outputs/full_data_projection/canonical_full_data.h5ad')
            write_json(self.output/'full_data_projection/full_data_projection_summary.json', full)
            self.save_checkpoint(pilot, original, program, report, full)
            write_json(self.output/'summary.json', report)
            pages = pd.read_csv(self.output/'annotation_review/fine_review_page_index.tsv', sep='\t')
            self.check('all_38_fine_review_pages', len(pages) == 38 and pages.cluster.nunique() == 38)
            fine_table = pd.read_csv(self.output/'annotation_review/cluster_validation_summary.tsv', sep='\t')
            self.check('fine_review_complete_partition', len(fine_table) == 38 and fine_table.n_cells.sum() == self.cfg['pilot_cells'])
            figure_index = pd.read_csv(self.output/'figures/figure_index.tsv', sep='\t')
            for row in figure_index.itertuples():
                for relative in [row.png, row.pdf] + str(row.source_tables).split(';'):
                    self.check('figure_asset_'+relative, (self.output/relative).is_file())
            for name in ['hicat_hierarchy_validation_main.pdf', 'hicat_hierarchy_validation_detailed.pdf']:
                self.check('report_'+name, (self.output/'figures'/name).is_file())
            self.check('review_readme_exists', (self.output/'REVIEW_README.md').is_file())
            write_json(self.output/'software_versions.json', dict(python=sys.version,
                packages={d.metadata['Name']:d.version for d in importlib.metadata.distributions()}))
            write_json(self.output/'STEP_STATUS.json', dict(status='IN_REVIEW', run_id=self.root.name,
                stage='02_hierarchy_validation', full_data_clustering_started=False, annotations_locked=False))
            for path in self.output.rglob('*.json'):
                text = path.read_text()
                if str(self.output) in text:
                    path.write_text(text.replace(str(self.output), str(self.root/'outputs')))
            manifest(self.output).to_csv(self.output/'output_manifest.tsv', sep='\t', index=False)
        destination = self.root/'outputs'
        if destination.exists():
            raise FileExistsError(str(destination))
        os.replace(self.output, destination)
        if (self.workspace/'REVIEW_README.md').exists():
            os.replace(self.workspace/'REVIEW_README.md', self.root/'REVIEW_README.md')
        self.progress.note('workflow', 'COMPLETE', {}, dict(outputs=str(destination), status='IN_REVIEW'))
        return destination


def main():
    """Execute only the explicitly named frozen Step 07 package."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', type=Path, required=True)
    args = parser.parse_args()
    print(ValidationWorkflow(args.run_dir).run(), flush=True)


if __name__ == '__main__':
    main()
