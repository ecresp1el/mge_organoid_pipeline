# DIV90 PI figure package

Created 2026-09-09 from the completed phase-two outputs. This is a presentation package, with no biological refitting, new gate discovery, UMAP calculation, or changes to original results.

## Runtime output

`/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_hypergate_sst_pv_phase2/pi_figure_package_v1/`

The share archive is `DIV90_Hypergate_PI_Figure_Package.zip`. On the Mac Turbo mount, replace `/nfs/turbo/umms-parent` with `/Volumes/umms-parent`.

The package contains exactly four main figures (PDF/SVG/600-dpi PNG), ten supplemental figures (one ten-page PDF and individual SVG/600-dpi PNG files), three standalone offline HTML pages, the one-page PI summary (PDF/PNG), and `PI_FIGURE_GUIDE.md`. Read the guide for each figure's conclusion, caveat and suggested 15-second explanation. `MAIN_FIGURE_CAPTIONS.md` and `supplement/SUPPLEMENTAL_CAPTIONS.md` contain display details; the latter gives full-precision RNA thresholds.

## Preserved interpretation

- Frozen entry: 4,768 cortical LHX6+/ERBB4+ cells; PVALB RNA detected in only three.
- Independent PV and SST scores define 1,076 PV-biased, 1,309 dual-high, 1,075 SST-biased and 1,308 unresolved cells. The source label `PV/SST hybrid` is displayed as dual-high; no underlying label was changed.
- The original FAT3/PTPRM retention rule preserves 79.1% of PV-biased and 47.9% of dual-high cells, losing 682 dual-high cells.
- Removing FGFR2-high OR PTPRS-high cells preserves 88.7% PV-biased, 89.3% dual-high and 84.4% of all entry cells. Retain the complement (at or below both frozen RNA thresholds). This is an RNA hypothesis for a prospective experiment, with modest SST depletion.
- Figure 4 proposes collecting retained and removed fractions plus an unsorted ERBB4+ control, identical maturation, and later protein, transcriptomic and physiological comparisons.
- CV/MW effects remain in the supplement and dedicated interactive page, separately by H9, 79B and 2E; 2E/CV has only 37 cells. Condition includes the reported approximately twofold glucose difference and is confounded with operator/culture differences.

Figure 1 uses the saved pre-entry local reclusters in `final_figures/div90_guidance_with_sst/tables/joined_with_sst.tsv.gz`. ERBB4 RNA is detected in 81.6% of 9,999 cortical profiles versus 18.4% of 6,121 subpallial profiles with available expression. Unavailable expression profiles are retained only as spatial context. Separate recluster distances are not comparable. ERBB4 is cortical-enriched, not exclusive; ERBB4-only live selection does not reproduce the recovered cortical/LHX6+ RNA criteria automatically.

## Code and checks

Orchestrator: [build_div90_hypergate_pi_package.py](scripts/build_div90_hypergate_pi_package.py). Separate modules render [main figures](scripts/div90_hypergate_pi_main.py), [supplement](scripts/div90_hypergate_pi_supplement.py), and [interactive pages](scripts/div90_hypergate_pi_interactive.py).

Use `/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/python` with `PYTHONDONTWRITEBYTECODE=1`; all rendering caches and runtime files belong on Turbo. The orchestrator supports `--stage entry|guide|all|package`. Do not rerun it over a published version; choose a new output version for changes. No rerender is necessary to inspect or share this completed package.

`PACKAGE_VALIDATION.json` checks file presence, page counts, main-figure vector assets, and 600-dpi PNGs. `provenance/interactive_validation.json` records offline Chromium checks of counts, filter controls, hover, zoom/reset, all cell-line comparisons and consistent percentage rounding. `provenance/source_preservation.json` records SHA256 equality for source artifacts. `FILE_MANIFEST.json` inventories packaged file hashes; `provenance/zip_validation.json` provides the archive checksum and CRC result. The browser runtime reuses the existing phase-two Turbo installation.
