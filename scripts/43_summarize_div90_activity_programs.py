#!/usr/bin/env python3
"""Create complete, readable review tables from finished outputs; no new tests."""
import collections
import csv
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import time

os.environ['TZ']='America/Detroit';time.tzset()
run=Path(sys.argv[1]);tables=run/'tables'
now=dt.datetime.now().astimezone()
review=run/('review_{:%Y%m%d_%H%M%S_%Z}'.format(now));review.mkdir(exist_ok=False)
shutil.copy2(__file__,review/Path(__file__).name)
inputs=[]
def rows(name):
 p=tables/name;inputs.append(p)
 with p.open()as f:
  yield from csv.DictReader(f)
annotations=collections.defaultdict(set)
for r in rows('expanded_gene_program_annotations.csv'):annotations[r['gene']].update(r['categories'].split('|'))
present={r['gene']for r in rows('gene_coverage.csv')if r['present_in_object']=='TRUE'}
labels={'3':'CRABP1/ANGPT2','5':'NR2F1/NR2F2','9':'Shared neuronal internal branch','1':'LHX8/ISL1 state 1','2':'LHX8/ISL1 state 2','4':'LHX6/MEF2C-associated'}
branchpoint={'3':'BP11','5':'BP11','9':'BP11','1':'BP9','2':'BP9','4':'BP9'}
whole={}
for r in rows('individual_gene_branch_sample_expression.csv'):
 if r['sample']=='ALL' and r['segment']in labels:whole[(r['segment'],r['gene'])]=r
consistency={(r['segment'],r['gene']):r for r in rows('individual_gene_sample_consistency.csv')}
near=collections.defaultdict(dict)
for r in rows('consolidated_gene_evidence_all_early_tests.csv'):
 if r['gene']not in present:continue
 s=r['reported_daughter'];opp=r['daughter_B']if s==r['daughter_A']else r['daughter_A']
 near[(s,r['gene'])][opp]={'log2fc':float(r['log2fc']),'FDR':float(r['global_fdr']),'significant_elevated':r['candidate']=='TRUE'}
later=collections.defaultdict(list)
for r in rows('later_native_URD_all_tested_genes.csv'):
 if r['gene']not in present:continue
 for s,opp,sign in [(r['daughter_A'],r['daughter_B'],1),(r['daughter_B'],r['daughter_A'],-1)]:
  fc=sign*float(r['log2fc']);later[(s,r['gene'])].append({'opponent':opp,'log2fc':fc,'FDR':float(r['later_global_fdr']),'significant_elevated':fc>=.25 and float(r['later_global_fdr'])<=.05})
columns=['branchpoint','branch','branch_label','gene','programs','n_cells','n_expressing','percent_expressing','mean_normalized','mean_log2upx',
 'near_split_sibling_tests_available','near_split_siblings_significantly_lower','near_split_status','near_split_comparisons_JSON',
 'later_comparisons_JSON','samples_with_trend_support','samples_increasing','samples_decreasing','median_sample_first_to_last_log_expression_change']
lists=collections.defaultdict(list)
with (review/'ALL_2971_GENES_ALL_SIX_BRANCHES.csv').open('w')as f:
 w=csv.DictWriter(f,fieldnames=columns);w.writeheader()
 for s in ['3','5','9','1','2','4']:
  for g in sorted(present):
   key=(s,g);n=near[key];wins=sorted(k for k,v in n.items()if v['significant_elevated'])
   status='higher_than_both'if len(wins)==2 else 'higher_than_one'if len(wins)==1 else 'not_significantly_elevated_in_available_tests'
   if wins:lists[(s,'Near split: '+status)].append(g)
   if any(x['significant_elevated']for x in later[key]):lists[(s,'Later: elevated in the supported pair')].append(g)
   activity=any('activity'in x or 'calcium'in x for x in annotations[g])
   if activity and wins:lists[(s,'Activity/calcium near split: '+status)].append(g)
   a=whole.get(key,{});c=consistency.get(key,{})
   row={'branchpoint':branchpoint[s],'branch':s,'branch_label':labels[s],'gene':g,'programs':'|'.join(sorted(annotations[g])),
    'near_split_sibling_tests_available':len(n),'near_split_siblings_significantly_lower':'|'.join(wins),'near_split_status':status,
    'near_split_comparisons_JSON':json.dumps(n,sort_keys=True),'later_comparisons_JSON':json.dumps(later[key],sort_keys=True),
    'samples_with_trend_support':c.get('samples_with_support',''),'samples_increasing':c.get('samples_increasing',''),
    'samples_decreasing':c.get('samples_decreasing',''),'median_sample_first_to_last_log_expression_change':c.get('median_delta','')}
   for k in ['n_cells','n_expressing','percent_expressing','mean_normalized','mean_log2upx']:row[k]=a.get(k,'')
   w.writerow(row)
with (review/'COMPLETE_GENE_LISTS.txt').open('w')as f:
 f.write('Complete gene lists from the finished DIV90 activity-program analysis\nGenerated '+now.isoformat()+'\n\n')
 f.write('All '+str(len(present))+' profiled genes (alphabetical):\n'+', '.join(sorted(present))+'\n\n')
 for s in ['3','5','9','1','2','4']:
  f.write('BRANCH '+s+' '+labels[s]+'\n')
  for key in sorted(k for k in lists if k[0]==s):f.write(key[1]+' ('+str(len(lists[key]))+'):\n'+', '.join(lists[key])+'\n\n')
 f.write('Near-split and later tests have separate FDR families. One-sibling genes can have an untested other sibling; consult CSV eligibility.\nOnly BP11 5/3 and BP9 1/2 had supported later tests. No direct CRABP1/MEF2C DE test was performed.\nPseudotime trends are descriptive; late-only significance does not establish onset.\n')
manifest={'generated_at':now.isoformat(),'source_run':str(run),'profiled_genes':len(present),'branch_gene_rows':len(present)*6,'new_tests_performed':False,
 'inputs':{str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in dict.fromkeys(inputs)}}
(review/'provenance.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(review,flush=True)
