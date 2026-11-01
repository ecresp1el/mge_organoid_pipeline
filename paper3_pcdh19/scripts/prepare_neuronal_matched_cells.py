from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import numpy as np,json
ROOT=Path(__file__).resolve().parents[1];SOURCE=RUNTIME_ROOT/'results/neuronal_mosaic_focus_20261031';OUT=RUNTIME_ROOT/'results/neuronal_matched_selection_20261031'
c=np.load(str(SOURCE/'neuronal_violin_cells.npz'));e=np.load(str(SOURCE/'neuronal_exon1_signal.npz'))
assert np.array_equal(c['group'],e['group'])
keep=(c['group']!=0)|(e['raw_exon1']>0)
np.savez_compressed(str(OUT/'neuronal_violin_cells.npz'),**{k:(c[k] if k=='genes' else c[k][keep]) for k in c.files})
np.savez_compressed(str(OUT/'neuronal_exon1_signal.npz'),**{k:e[k][keep] for k in e.files})
assert np.all(e['raw_exon1'][keep & (c['group']==0)]>0)
assert np.all(e['raw_exon1'][keep & (c['group']==1)]>0)
assert np.all(e['raw_exon1'][keep & (c['group']==2)]==0)
(OUT/'matched_selection_provenance.json').write_text(json.dumps({'selection':{'WT females':'A-positive inhibitory neurons','HET WT-like':'A-positive inhibitory neurons','HET KO-like':'A-negative, downstream-positive inhibitory neurons'},'cells_by_group':np.bincount(c['group'][keep]).tolist(),'group_1_and_2_same_criterion':True,'exon1_panel_is_selection_diagnostic':True,'expression_values_unchanged':True,'all_panels_use_matched_WT_reference':True},indent=2))
print('MATCHED_SELECTION_READY',np.bincount(c['group'][keep]).tolist())
