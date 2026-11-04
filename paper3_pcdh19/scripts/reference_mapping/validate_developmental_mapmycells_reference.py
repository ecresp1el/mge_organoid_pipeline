#!/usr/bin/env python3
"""Validate cached developmental labels with the actual MapMyCells taxonomy API.

Repeated labels such as Undefined need path-qualified machine IDs; the original
published labels remain unchanged in separate columns. This is metadata audit,
not reference selection, marker fitting, or query mapping.
"""
import argparse,json
from pathlib import Path
import h5py,pandas as pd
from anndata.io import read_elem
from cell_type_mapper.taxonomy.taxonomy_tree import TaxonomyTree


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audit-dir',type=Path,required=True);args=parser.parse_args()
    out=args.audit_dir;record=json.loads((out/'REFERENCE_VALIDATION.json').read_text())
    with h5py.File(record['LaManno']['reference'],'r') as h:obs=read_elem(h['obs'])
    published=['Class','Subclass','CellType'];working=obs[published].astype(str).copy()
    # Disambiguate repeated strings by their existing published parent path;
    # create no biological categories or manual cell assignments.
    columns=[]
    for depth in range(3):
        name='mmc_'+published[depth]+'_id';columns.append(name)
        working[name]=[json.dumps(list(row),ensure_ascii=False,separators=(',',':')) for row in working[published[:depth+1]].itertuples(index=False,name=None)]
    tree=TaxonomyTree.from_dataframe(working,columns)
    (out/'LaManno_path_qualified_taxonomy.json').write_text(tree.to_str())
    working.to_csv(out/'LaManno_published_labels_and_machine_IDs.tsv.gz',sep='\t',index_label='reference_cell_id',compression='gzip')
    record['LaManno']['checks']['path_qualified_MapMyCells_taxonomy_constructed']=True
    record['LaManno']['machine_hierarchy']=columns
    record['LaManno']['machine_ID_policy']='JSON-encoded published parent path; original labels retained verbatim; no biological regrouping'
    join=pd.read_csv(record['Bandler_CA301']['label_join'],sep='\t')
    assert len(join)==4481 and join.cell_id.is_unique and join.region.eq('MGE').all() and join.genotype.eq('WT').all()
    tree=TaxonomyTree.from_dataframe(join.set_index('cell_id'),['later_atlas_class','later_atlas_cluster'])
    (out/'Bandler_verified_MGE_taxonomy.json').write_text(tree.to_str())
    (out/'REFERENCE_VALIDATION.json').write_text(json.dumps(record,indent=2)+'\n')
    print('PASS: actual MapMyCells taxonomy API accepts both audited developmental hierarchies; no query mapping performed')


if __name__=='__main__':main()
