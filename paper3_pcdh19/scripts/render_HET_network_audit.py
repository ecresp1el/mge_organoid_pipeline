from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=RUNTIME_ROOT/'results/female_male_signaling_20261031'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':6,'pdf.fonttype':42,'svg.fonttype':'none','axes.linewidth':.5})
z=np.load(str(OUT/'HET_network_difference_matrices.npz'));raw=z['differences'][1:];depth=z['depth_balanced_differences'][:,1:].mean(0);bound=max(1,float(np.max(np.abs(raw))));depthbound=max(1,float(np.max(np.abs(depth))));short=['RG-like','MGE IPCs','Inhibitory','Cycling','Neurogenic'];titles=['WT-like → KO-like','KO-like → WT-like','KO-like → KO-like']
def draw_deltas(fig,boxes,matched=True,fs=5.2,scale=None):
 vals=depth if matched else raw;limit=depthbound if matched else bound
 if scale is not None:limit=scale
 for k,box in enumerate(boxes):
  ax=fig.add_axes(box);im=ax.imshow(vals[k],cmap='RdBu_r',vmin=-limit,vmax=limit,aspect='equal');ax.set_title(titles[k],fontsize=fs+.5,pad=3);ax.set_xticks(range(5));ax.set_xticklabels(short,fontsize=fs,rotation=50,ha='right');ax.set_yticks(range(5));ax.set_yticklabels(short if k==0 else ['']*5,fontsize=fs);ax.tick_params(length=0,pad=2)
  for spine in ax.spines.values():spine.set_visible(False)
  for i in range(5):
   for j in range(5):
    value=vals[k,i,j];ax.text(j,i,'{:+.0f}'.format(value) if abs(value)>=.5 else '0',ha='center',va='center',fontsize=fs-.4,color='white' if abs(value)>.65*limit else '#26292E')
  if k==0:ax.set_ylabel('Sender',fontsize=fs+.3)
 return im
fig=plt.figure(figsize=(7.2,4.8));fig.text(.035,.965,'Female HET · changes in supported pair counts',fontsize=10,weight='bold')
for yy,label,matched in [(.58,'Original depth',False),(.13,'Matched depth',True)]:
 fig.text(.035,yy+.32,label,fontsize=7,weight='bold');im=draw_deltas(fig,[[.12,yy,.22,.27],[.43,yy,.22,.27],[.74,yy,.22,.27]],matched,6,scale=bound)
cbax=fig.add_axes([.39,.025,.27,.015]);cb=fig.colorbar(im,cax=cbax,orientation='horizontal',ticks=[-bound,0,bound]);cb.ax.tick_params(length=0,labelsize=6);fig.text(.035,.03,'vs WT-like → WT-like',fontsize=6);cb.outline.set_visible(False)
for ext in ['png','pdf','svg']:fig.savefig(str(OUT/('HET_network_depth_audit.'+ext)),dpi=450,facecolor='white')
with open(str(OUT/'HET_network_continuous_pair_audit.tsv')) as f:allrows=list(csv.DictReader(f,delimiter='\t'))
family={'NRG':'Neuregulin','EPHA':'Ephrin','EPHB':'Ephrin','SLIT':'Slit','NCAM':'NCAM','PCDH':'Protocadherin'}
for x in ['SEMA3','SEMA4','SEMA5','SEMA6','SEMA7']:family[x]='Semaphorin'
eligible=[]
for r in allrows:
 if r['module'] not in family or r['consistent_direction']!='True':continue
 v=np.array([float(x) for x in r['sample_log2_changes'].split(';')]);r['replicated_effect_floor']=float(np.min(np.abs(v)))
 if r['replicated_effect_floor']<.2:continue
 eligible.append(r)
eligible.sort(key=lambda r:r['replicated_effect_floor'],reverse=True);chosen=[];seen=set();limits={}
for r in eligible:
 fam=family[r['module']]
 if r['interaction'] in seen or limits.get(fam,0)>=(2 if fam=='Protocadherin' else 1):continue
 chosen.append(r);seen.add(r['interaction']);limits[fam]=limits.get(fam,0)+1
 if len(chosen)==5:break
assert len(chosen)==5
shortname={'RG-like progenitors':'RG-like','MGE IPCs':'MGE IPCs','Immature inhibitory neurons':'Inhibitory neurons','Cycling neuroblasts':'Cycling neuroblasts','Neurogenic progenitors':'Neurogenic progenitors'}
f=plt.figure(figsize=(7.2,2.8));f.text(.035,.96,'HET KO-like vs WT-like · interaction candidates',fontsize=10,weight='bold');ax=f.add_axes([.69,.20,.27,.62]);effects=np.array([(2**float(r['median_log2_change'])-1)*100 for r in chosen]);ax.barh(range(5),effects,height=.58,color=['#326C9D' if v<0 else '#B64B50' for v in effects]);ax.axvline(0,color='#9FA5AD',lw=.5);ax.set_ylim(4.6,-.6);ax.set_xlim(-100,110);ax.set_xticks([-100,0,100]);ax.set_xticklabels(['−100%','0','+100%'],fontsize=6);ax.set_yticks([]);ax.set_xlabel('Median change in potential',fontsize=6);ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False);ax.spines['left'].set_visible(False)
lab=f.add_axes([.035,.20,.64,.62]);lab.set_xlim(0,1);lab.set_ylim(4.6,-.6);lab.axis('off')
for i,r in enumerate(chosen):
 l=r['ligand'];re=r['receptor'].replace('NRP1_PLXNA4','Nrp1/Plxna4');label=l+' contact' if l==re else l+'–'+re
 lab.text(0,i,label,fontsize=7,va='center');lab.text(.50,i,shortname[r['source_name']]+' → '+shortname[r['receiver_name']],fontsize=6,va='center')
 ax.text(effects[i]+(3 if effects[i]>0 else -3),i,'{:+.0f}%'.format(effects[i]),fontsize=6,ha='left' if effects[i]>0 else 'right',va='center')
f.text(.035,.86,'Same direction in 3/3 HET samples',fontsize=6)
for ext in ['png','pdf','svg']:f.savefig(str(OUT/('HET_replicated_interaction_changes.'+ext)),dpi=450,facecolor='white')
with open(str(OUT/'HET_replicated_interaction_selection.tsv'),'w') as ff:
 wr=csv.DictWriter(ff,fieldnames=list(chosen[0]),delimiter='\t');wr.writeheader();wr.writerows(chosen)
(OUT/'HET_NETWORK_INTERPRETATION.md').write_text('The original 100 raw edge counts were reconstructed exactly. Dense circular networks obscure differences; raw count increases largely disappear after matched-depth binomial thinning. The largest raw edge difference, MGE IPC to MGE IPC, changes from +256 to a mean +14.06 pairs; 50-thinning simulation range −13.775 to +50.775. Whole-network within-state difference falls from +2625 to a mean +111.7, simulation range −248.75 to +415.475. These simulation ranges are not biological confidence intervals. Main revised G shows mean depth-matched pair-count differences for each route versus WT-like to WT-like; supplementary audit shows original and matched depths on a common scale. Individual-pair candidates use continuous sqrt(ligand CPM × receptor CPM), geometric means for complexes, no zero pseudocounts. Every component is >=1 CPM in all three samples in both states. Candidates must show the same direction and absolute log2 change >=0.2 in each sample. Rank by the smallest absolute effect among the three samples, retain unique interactions, up to two protocadherin interactions and one per other predefined neural pathway family, display five. Selection is exploratory and does not establish FDR significance, functional signaling, or direct cell contact. Displayed median changes compare KO-like to KO-like with WT-like to WT-like. Count changes and continuous-pair effects answer different questions.\n')
print('REPLICATED_CANDIDATES',[(r['interaction'],r['source_name'],r['receiver_name'],round((2**float(r['median_log2_change'])-1)*100)) for r in chosen],flush=True)
