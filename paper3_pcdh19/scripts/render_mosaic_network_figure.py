from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import csv,gzip,json,copy
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Circle, Rectangle
from matplotlib.path import Path as MplPath
from scipy.cluster.hierarchy import linkage,dendrogram
ROOT=Path(__file__).resolve().parents[1];OUT=RUNTIME_ROOT/'results/female_male_signaling_20261031';FEMALE=RUNTIME_ROOT/'results/neuronal_matched_selection_20261031';DEEP=RUNTIME_ROOT/'results/mosaic_deep_dive_20261031';SIG=RUNTIME_ROOT/'results/paired_HET_signaling_20261031'
def read(p):
 with (gzip.open(str(p),'rt') if str(p).endswith('.gz') else open(str(p))) as f:return list(csv.DictReader(f,delimiter='\t'))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':6,'axes.titlesize':6.5,'axes.labelsize':6,'xtick.labelsize':5.5,'ytick.labelsize':5.5,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.5,'xtick.major.width':.5,'ytick.major.width':.5})
BLUE='#326C9D';RED='#B64B50';INK='#26292E';GRAY='#BFC4CC'
fig=plt.figure(figsize=(7.2,6.0))
rects={'a':[.055,.73,.62,.20],'b':[.13,.46,.21,.16],'c':[.42,.405,.555,.215],'d':[.75,.73,.225,.20],'e':[.13,.15,.21,.16],'f':[.42,.095,.555,.215]}

def title(k,text):
 x,y,w,h=rects[k];fig.text(x-.027,y+h+.020,k,fontsize=8,weight='bold',color=INK);fig.text(x,y+h+.020,text,fontsize=6.5,color=INK)
def clean(ax):ax.tick_params(length=2,pad=2)
# a: contextual atlas with broad, provenance-supported labels.
broad=[('RG-like progenitors',['C0002'],'#7186A6'),('MGE IPCs',['C0003','C0008'],'#CE9667'),('Neurogenic progenitors',['C0001','C0007','C0012','C0024'],'#B8A363'),('Cycling neuroblasts',['C0011'],'#9B789D'),('Immature inhibitory neurons',['C0004','C0005','C0006','C0020'],'#4D928A'),('Mixed neural / vascular',['C0015'],'#B77D86'),('Unresolved populations',[],GRAY)]
lookup={c:i for i,(_,cc,_) in enumerate(broad) for c in cc};cache=read(OUT/'atlas_sex_genotype_display.tsv.gz');xy=np.array([[float(r['UMAP1']),float(r['UMAP2'])] for r in cache]);category=np.array([lookup.get(r['cluster'],6) for r in cache]);sex=np.array([r['sex'] for r in cache]);geno=np.array([r['genotype'] for r in cache]);lim=[xy[:,0].min()-.4,xy[:,0].max()+.4,xy[:,1].min()-.4,xy[:,1].max()+.4]
po=np.random.RandomState(20261031).permutation(len(cache))
ax=fig.add_axes([.055,.745,.21,.17]);cols=np.array([broad[i][2] for i in category]);ax.scatter(xy[po,0],xy[po,1],s=.22,c=cols[po],alpha=.85,linewidths=0,rasterized=True);ax.set_xlim(lim[:2]);ax.set_ylim(lim[2:]);ax.set_aspect('equal');ax.axis('off');ax.set_title('Cell type',fontsize=6,pad=2)
insets=[('Female',sex=='Female','#0072B2',.315,.84),('Male',sex=='Male','#E69F00',.435,.84),('WT',geno=='WT','#222222',.315,.73),('HET',geno=='HET','#CC00AA',.435,.73),('KO',geno=='KO','#56B4E9',.555,.73)]
for name,mask,col,xx,yy in insets:
 ax=fig.add_axes([xx,yy,.105,.075]);ax.scatter(xy[po,0],xy[po,1],s=.10,c='#E6E7E9',linewidths=0,rasterized=True);ax.scatter(xy[mask,0],xy[mask,1],s=.19,c=col,alpha=.85,linewidths=0,rasterized=True);ax.set_xlim(lim[:2]);ax.set_ylim(lim[2:]);ax.set_aspect('equal');ax.axis('off');ax.set_title(name,fontsize=5.5,pad=1)
lg=fig.add_axes([.055,.653,.62,.065]);lg.axis('off');hh=[Line2D([0],[0],marker='o',ls='',color=col,label=name,markersize=2.5) for name,cc,col in broad];lg.legend(handles=hh,ncol=3,frameon=False,fontsize=5,loc='upper left',handletextpad=.35,handlelength=.5,columnspacing=.9,labelspacing=.25,borderaxespad=0)
title('a','')
# b/e: identical cell-level metric, axis range, violin settings and mean markers.
fcache=np.load(str(FEMALE/'neuronal_exon1_signal.npz'));mcache=np.load(str(OUT/'male_exon1_display.npz'));display_fractions={}
def violin(k,data,group,labels,text):
 ax=fig.add_axes(rects[k]);ds=[data[group==i] for i in range(len(labels))];non=[i for i,d in enumerate(ds) if np.ptp(d)>0];vp=ax.violinplot([ds[i] for i in non],positions=non,vert=False,widths=.64,showextrema=False,points=180,bw_method=.18)
 for body in vp['bodies']:body.set_facecolor('#C4C8CE');body.set_edgecolor('none');body.set_alpha(1);body.set_linewidth(0)
 for i,d in enumerate(ds):
  if i not in non:ax.scatter([0],[i],s=10,c=INK,marker='D',zorder=5,clip_on=False)
  q1,med,q3=np.percentile(d,[25,50,75]);ax.plot([q1,q3],[i,i],color=INK,lw=1.2,zorder=4);ax.scatter([d.mean()],[i],s=6,c=INK,zorder=5)
  ax.text(7.9,i-.29,'{:.1f}% detected'.format(100*(d>0).mean()),ha='right',va='bottom',fontsize=5,color='#666D75')
 ax.set_xlim(-.10,8.2);ax.set_xticks([0,4,8]);ax.set_ylim(len(ds)-.45,-.55);ax.set_yticks(range(len(ds)));ax.set_yticklabels(labels,fontsize=5.4);ax.set_xlabel('UMI / 10,000',fontsize=5.5,labelpad=2);clean(ax);ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0);title(k,text);display_fractions[k]=[float((d>8.2).mean()) for d in ds]
violin('b',fcache['exon1_CP10K'],fcache['group'],['WT female','HET WT-like','HET KO-like'],'Female · Pcdh19 exon 1')
violin('e',mcache['CP10K'],mcache['group'],['WT males','KO males'],'Male · Pcdh19 exon 1')
# A common pre-existing candidate panel makes sex comparisons direct.
fgenes=['Ank3','Kcnc2','Gnao1','Ptprd','Epha5','Snap47','Pcdha5','Pcdha10','Pcdhga5'];mgenes=['Eps8l1','Serpine2','Syt9','Ecel1','Lrfn5','Grm7','Cntn1','Sphkap','Nrp2'];genes=fgenes+mgenes
common_female={r['gene_symbol'] for r in read(DEEP/'DE_C0006/within_HET_Aneg_BC1_vs_Apositive.tsv.gz')};common_male={r['gene_symbol'] for r in read(OUT/'male_C0006_DE.tsv')};dropped=[g for g in genes if g not in common_female or g not in common_male];genes=[g for g in genes if g in common_female and g in common_male]
rows=read(OUT/'matched_display_logCPM.tsv');values={r['gene_symbol']:r for r in rows};md=read(OUT/'matched_female_display_metadata.tsv');go=[('WT_F','A_positive'),('HET_F','A_positive'),('HET_F','Aneg_BC1')];md=sorted(md,key=lambda r:(go.index((r['design_group'],r['gate'])),r['sample']));flibs=[r['library'] for r in md];mlibs=['15662-JZ-'+str(i)+'__C0006' for i in [1,2,3,10,11,12]]
fm=np.array([[np.mean([float(values[g][l]) for l in flibs[k:k+3]]) for k in [0,3,6]] for g in genes]);mm=np.array([[np.mean([float(values[g][l]) for l in mlibs[k:k+3]]) for k in [0,3]] for g in genes])
def standardize(a):return (a-a.mean(1)[:,None])/np.maximum(a.std(1)[:,None],1e-9)
fx=np.array([[float(values[g][l]) for l in flibs] for g in genes]);mx=np.array([[float(values[g][l]) for l in mlibs] for g in genes])
fz=(fm-fx.mean(1)[:,None])/np.maximum(fx.std(1)[:,None],1e-9);mz=(mm-mx.mean(1)[:,None])/np.maximum(mx.std(1)[:,None],1e-9);link=linkage(np.column_stack([fz,mz]),method='average');fd={r['gene_symbol']:r for r in read(DEEP/'DE_C0006/within_HET_Aneg_BC1_vs_Apositive.tsv.gz')};mdict={r['gene_symbol']:r for r in read(OUT/'male_C0006_DE.tsv')};fef=np.array([float(fd[g]['logFC']) if g in fd else np.nan for g in genes]);mef=np.array([float(mdict[g]['logFC']) if g in mdict else np.nan for g in genes]);orders={}
def matched_gene_block(k,z,effect,labels,text,contrast):
 x,y,w,h=rects[k];ax=fig.add_axes([x,y,w,h]);dd=dendrogram(link,no_plot=True);order=dd['leaves'];md=max(float(link[-1,2]),1e-9);expr_width=2.3;base=4.25;scale=.48
 # One coordinate system keeps every branch connected and places the leaves on the heatmap boundary.
 for ic,dc in zip(dd['icoord'],dd['dcoord']):ax.plot(-.8*np.asarray(dc)/md,(np.asarray(ic)-5)/10,color='#555B63',lw=.6,zorder=2)
 ax.imshow(z[order],extent=[0,expr_width,len(genes)-.5,-.5],aspect='auto',cmap='Purples',vmin=-1.5,vmax=1.5,zorder=1);ax.add_patch(Rectangle((0,-.5),expr_width,len(genes),fill=False,ec='#555B63',lw=.35,zorder=3))
 for i,j in enumerate(order):ax.text(expr_width+.10,i,genes[j],ha='left',va='center',fontsize=5.6)
 ev=effect[order];ax.barh(range(len(genes)),ev*scale,left=base,height=.60,color=[BLUE if q<0 else RED for q in ev],zorder=3);ax.plot([base,base],[-.5,len(genes)-.5],color='#A1A5AC',lw=.4)
 ax.set_xlim(-.85,5.10);ax.set_ylim(len(genes)-.5,-.5);ax.axis('off')
 for j,lab in enumerate(labels):ax.text((j+.5)*expr_width/len(labels),len(genes)-.15,lab.replace(' ','\n',1),ha='center',va='top',fontsize=5.2,clip_on=False)
 for tick in [-1,0,1]:
  xx=base+tick*scale;ax.plot([xx,xx],[len(genes)-.5,len(genes)-.35],color=INK,lw=.4,clip_on=False);ax.text(xx,len(genes)-.10,str(tick),ha='center',va='top',fontsize=5.2,clip_on=False)
 ax.text(base,len(genes)+.65,'log2FC',ha='center',va='top',fontsize=5.5,clip_on=False);ax.text(expr_width/2,-.85,'Relative expression',ha='center',va='bottom',fontsize=5.3,clip_on=False);ax.text(base,-.85,contrast,ha='center',va='bottom',fontsize=5.3,clip_on=False);title(k,text);orders[k]=[genes[i] for i in order]
matched_gene_block('c',fz,fef,['WT female','HET WT-like','HET KO-like'],'Female','KO-like / WT-like')
matched_gene_block('f',mz,mef,['WT males','KO males'],'Male','KO / WT')
assert orders['c']==orders['f']
# Individual interactions grouped by module; no family averaging.
ex={(r['stratum'],r['gene'],r['gate']):r for r in read(SIG/'signaling_gene_expression_support.tsv')}
state={'WT-like':'A_positive','KO-like':'Aneg_BC1'}
routekeys=[('WT-like','WT-like'),('WT-like','KO-like'),('KO-like','WT-like'),('KO-like','KO-like')]
families={'NRG':'Neuregulin','EPHA':'Ephrin','EPHB':'Ephrin','SLIT':'Slit','NCAM':'NCAM','PCDH':'Protocadherin'}
for fam in ['SEMA3','SEMA4','SEMA5','SEMA6','SEMA7']:families[fam]='Semaphorin'
def component(gs,gate):
 rr=[ex.get(('C0006',g,gate)) for g in gs]
 if any(r is None for r in rr):return None
 vals=np.array([[float(v) for v in r['sample_CPM'].split(';')] for r in rr])
 # Avoid near-zero denominators: every subunit expressed in all three samples in both states.
 if np.any(vals<1):return None
 return np.exp(np.mean(np.log(vals),axis=0))
audit=[];candidates=[]
for pair in read(SIG/'expanded_signaling_pairs.tsv'):
 if pair['pathway'] not in families:continue
 lg=pair['ligand_genes'].split(';');rg=pair['receptor_genes'].split(';')
 if 'Pcdh19' in lg+rg:continue
 ll={k:component(lg,v) for k,v in state.items()};rr={k:component(rg,v) for k,v in state.items()}
 if any(v is None for v in list(ll.values())+list(rr.values())):continue
 scores=np.array([np.sqrt(ll[x]*rr[y]) for x,y in routekeys]);means=scores.mean(axis=1);relative=means/means[0];paired=scores[3]/scores[0]
 consistent=bool(np.all(paired<1) or np.all(paired>1))
 rec=dict(module=families[pair['pathway']],interaction=pair['interaction'],ligand=pair['ligand'],receptor=pair['receptor'],WT_WT=means[0],WT_KO=relative[1],KO_WT=relative[2],KO_KO=relative[3],sample_ratios=';'.join(str(v) for v in paired),consistent=consistent)
 audit.append(rec)
 if consistent:candidates.append((rec,relative))
# Largest absolute within-state changes: three modules, up to two interactions each.
module_strength={}
for rec,rel in candidates:module_strength[rec['module']]=max(module_strength.get(rec['module'],0),abs(float(np.log2(rel[3]))))
selected_modules=sorted(module_strength,key=module_strength.get,reverse=True)[:3]
selected=[]
for module in selected_modules:
 eligible=sorted([x for x in candidates if x[0]['module']==module],key=lambda x:abs(float(np.log2(x[1][3]))),reverse=True)
 selected.extend(eligible[:2])
assert len(selected)>=4
names=[r['interaction'] for r,v in selected];matrix=np.log2(np.array([v for r,v in selected]));module_rows=[r for r,v in selected]
for rec in audit:rec['displayed']=rec['interaction'] in names
with open(str(OUT/'interaction_module_audit.tsv'),'w') as f:
 wr=csv.DictWriter(f,fieldnames=list(audit[0]),delimiter='\t');wr.writeheader();wr.writerows(audit)
def draw_interactions(target,rect,fs=5):
 ax=target.add_axes(rect);plotrows=[];positions=[];headers=[];last=None
 for i,(rec,rel) in enumerate(selected):
  if rec['module']!=last:
   headers.append((len(plotrows),rec['module']));plotrows.append([np.nan]*4);last=rec['module']
  positions.append(len(plotrows));plotrows.append(matrix[i])
 im=ax.imshow(np.array(plotrows),aspect='auto',cmap='RdBu_r',vmin=-1,vmax=1)
 # Module labels occupy their own rows in the margin; interaction labels remain concise.
 labels=[]
 for rec,rel in selected:
  lg=rec['ligand'];rg=rec['receptor'].replace('PLXNA1_KDR','Plxna1/Kdr');labels.append(lg+'–'+rg if lg!=rg else lg+' contact')
 ax.set_yticks(positions);ax.set_yticklabels(labels,fontsize=fs)
 ax.set_xticks(range(4));ax.set_xticklabels(['WT→WT','WT→KO','KO→WT','KO→KO'],rotation=50,ha='right',fontsize=fs);ax.tick_params(length=0,pad=2)
 for sp in ax.spines.values():sp.set_visible(False)
 for i in range(len(names)):
  for j in range(4):ax.text(j,positions[i],'{:.2f}'.format(2**matrix[i,j]),ha='center',va='center',fontsize=fs-.3,color='white' if abs(matrix[i,j])>.65 else INK)
 for yy,label in headers:ax.text(-.05,yy,label,transform=ax.get_yaxis_transform(),ha='right',va='center',fontsize=fs-.3,weight='bold',clip_on=False)
 return im
im=draw_interactions(fig,[.815,.745,.16,.185],4.8)
title('d','HET signaling')
cax=fig.add_axes([.835,.665,.12,.007]);cb=fig.colorbar(im,cax=cax,orientation='horizontal',ticks=[-1,0,1]);cb.ax.set_xticklabels(['½','1','2']);cb.ax.tick_params(length=0,pad=1,labelsize=4.8);cb.outline.set_visible(False)
# A separate larger D makes the interaction view available without the full figure.
detail=plt.figure(figsize=(4.5,3.2));di=draw_interactions(detail,[.49,.27,.47,.61],8)
detail.text(.06,.96,'Female HET · signaling interactions',fontsize=10,weight='bold')
dc=detail.add_axes([.61,.08,.28,.025]);cb2=detail.colorbar(di,cax=dc,orientation='horizontal',ticks=[-1,0,1]);cb2.ax.set_xticklabels(['½','1','2']);cb2.ax.tick_params(length=0,labelsize=7);cb2.set_label('Relative expression potential',fontsize=7,labelpad=2);cb2.outline.set_visible(False)
for ext in ['png','pdf','svg']:detail.savefig(str(OUT/('HET_signaling_interactions.'+ext)),dpi=450,facecolor='white')
# One compact shared key; black labels remain readable independently of colour.
cm=plt.get_cmap('Purples');legend=fig.add_axes([.59,.030,.12,.010]);gradient=np.tile(np.linspace(-1.5,1.5,100),(2,1));legend.imshow(gradient,aspect='auto',cmap=cm,vmin=-1.5,vmax=1.5);legend.axis('off');fig.text(.59,.018,'Low',fontsize=5);fig.text(.71,.018,'High',fontsize=5,ha='right');fig.text(.75,.033,'Expression scaled within sex',fontsize=5.5)
fig.add_artist(plt.Rectangle((.06,.031),.012,.010,transform=fig.transFigure,fc=BLUE,ec='none'));fig.text(.078,.031,'Lower in KO-like / KO',fontsize=5.3);fig.add_artist(plt.Rectangle((.29,.031),.012,.010,transform=fig.transFigure,fc=RED,ec='none'));fig.text(.308,.031,'Higher in KO-like / KO',fontsize=5.3)

# Add cell-type interaction-count networks below the existing figure.
network_rows=[r for r in read(OUT/'celltype_LR_network_counts.tsv') if float(r['CPM_threshold'])==1.]
cts=['C0002','C0003','C0006','C0011','C0001']
node_names=['RG-like','MGE IPCs','Immature\ninhibitory','Cycling\nneuroblasts','Neurogenic\nprogenitors']
node_cols=['#7186A6','#CE9667','#4D928A','#9B789D','#B8A363']
network_groups=['WT female','HET WT-like','HET KO-like']
max_count=max(int(r['pairs']) for r in network_rows);assert max_count>0
def network(target,box,group,fontsize=5.3):
 ax=target.add_axes(box);ax.set_xlim(-1.72,1.72);ax.set_ylim(-1.5,1.6);ax.set_aspect('equal');ax.axis('off')
 theta=np.pi/2-np.arange(5)*2*np.pi/5;pos=np.column_stack([np.cos(theta),np.sin(theta)])
 rr={(r['source'],r['receiver']):int(r['pairs']) for r in network_rows if r['group']==group}
 for i,source in enumerate(cts):
  for j,targetct in enumerate(cts):
   count=rr[(source,targetct)]
   if not count:continue
   width=2.8*count/max_count
   if i==j:
    xy=pos[i];out=xy/np.linalg.norm(xy);perp=np.array([-out[1],out[0]])
    pts=[xy+.10*perp,xy+.30*out+.22*perp,xy+.30*out-.22*perp,xy-.10*perp]
    path=MplPath(pts,[MplPath.MOVETO,MplPath.CURVE4,MplPath.CURVE4,MplPath.CURVE4]);patch=FancyArrowPatch(path=path,arrowstyle='-|>',mutation_scale=4.5,lw=width,color=node_cols[i],alpha=.42,zorder=1)
   else:patch=FancyArrowPatch(pos[i],pos[j],connectionstyle='arc3,rad=0.20',arrowstyle='-|>',mutation_scale=4.5,shrinkA=5,shrinkB=5,lw=width,color=node_cols[i],alpha=.38,zorder=1)
   ax.add_patch(patch)
 ax.scatter(pos[:,0],pos[:,1],s=30,c=node_cols,edgecolors='white',linewidths=.65,zorder=4)
 for i,xy in enumerate(pos):
  labpos=xy*1.34;ha='center' if abs(xy[0])<.1 else ('left' if xy[0]>0 else 'right')
  ax.text(labpos[0],labpos[1],node_names[i],ha=ha,va='center',fontsize=fontsize,zorder=5,bbox=dict(facecolor='white',edgecolor='none',pad=.25))
 ax.set_title(group,fontsize=fontsize+1,pad=1)
 return ax
# Preserve the physical dimensions of A-F when adding the new row.
original_height=6.;added_height=1.65;full_height=original_height+added_height
fig.set_size_inches(7.2,full_height)
for oldax in list(fig.axes):
 pp=oldax.get_position();oldax.set_position([pp.x0,(pp.y0*original_height+added_height)/full_height,pp.width,pp.height*original_height/full_height])
for tx in fig.texts:
 xx,yy=tx.get_position();tx.set_position((xx,(yy*original_height+added_height)/full_height))
for artist in fig.artists:
 if isinstance(artist,Rectangle):artist.set_y((artist.get_y()*original_height+added_height)/full_height);artist.set_height(artist.get_height()*original_height/full_height)
fig.text(.028,.198,'g',fontsize=8,weight='bold');fig.text(.055,.198,'Cell-type interaction networks',fontsize=6.5)
for xx,group in zip([.025,.345,.665],network_groups):network(fig,[xx,.024,.31,.161],group)
# Shared width key, in pair counts, for every network.
kax=fig.add_axes([.38,.004,.27,.018]);kax.axis('off')
for xx,val in zip([.02,.35,.68],[max(1,int(max_count/4)),max(1,int(max_count/2)),max_count]):
 kax.plot([xx,xx+.14],[.5,.5],lw=2.8*val/max_count,color='#555B63');kax.text(xx+.16,.5,str(val),ha='left',va='center',fontsize=5)
kax.set_xlim(0,1);kax.set_ylim(0,1)
fig.text(.22,.012,'Ligand–receptor pairs',fontsize=5.2)
# Larger comparison export uses exactly the same edges and width scale.
nfig=plt.figure(figsize=(9,3.1))
for xx,group in zip([.015,.345,.675],network_groups):network(nfig,[xx,.16,.31,.74],group,7)
nfig.text(.025,.96,'Expression-supported cell-type networks',fontsize=10,weight='bold')
nfig.text(.12,.055,'Ligand–receptor pairs',fontsize=7)
nk=nfig.add_axes([.40,.028,.43,.06]);nk.axis('off')
for xx,val in zip([.02,.35,.68],[max(1,int(max_count/4)),max(1,int(max_count/2)),max_count]):
 nk.plot([xx,xx+.14],[.5,.5],lw=2.8*val/max_count,color='#555B63');nk.text(xx+.16,.5,str(val),ha='left',va='center',fontsize=7)
nk.set_xlim(0,1);nk.set_ylim(0,1)
for ext in ['png','pdf','svg']:nfig.savefig(str(OUT/('celltype_interaction_networks.'+ext)),dpi=450,facecolor='white')
for ext in ['png','pdf','svg']:fig.savefig(str(OUT/('mosaic_network_figure.'+ext)),dpi=450,facecolor='white')
(OUT/'mosaic_network_validation.json').write_text(json.dumps({'figure_size_inches':[7.2,6.0],'atlas_views':['broad cell identity','Female inset','Male inset','WT inset','HET inset','KO inset'],'atlas_same_embedding':True,'atlas_render_order':'one fixed random interleaving, shared by all views; no sex/genotype drawn last as a block','atlas_cells':len(cache),'short_violin_height_inches':rects['b'][3]*6,'violin_orientation':'horizontal; no connecting mean lines; quartile strip and mean point; detection fraction displayed','tree_rendering':'one shared coordinate system; complete linkage segments meet the heatmap edge; gene labels to the right of the matrix','directional_panel':'curated neural signaling modules across all four HET state routes; expression potential normalized per module to WT-like within-state; not measured communication','module_names':selected_modules,'interaction_names':names,'selection':'all subunits >=1 CPM in every HET sample in both states; KO/WT within-state direction consistent in three samples; top three modules by largest absolute change, top two interactions per module','inset_palette':'separate from broad cell identity palette','genes_dropped_for_missing_female_contrast':dropped,'body_font_pt':[5.1,6.5],'B_E_same_measurement':'cell-level exon1 UMI / total cell UMI * 10000','B_E_same_x_range':[-.10,8.2],'B_E_same_KDE_bandwidth':.18,'B_E_above_range_fractions':display_fractions,'male_group_selection':'known genotype, all neurons; no probe gate','female_group_selection':'matched A-positive WT and HET WT-like; HET KO-like A-negative/downstream-positive','C_F_same_genes_and_order':orders['c'],'C_F_same_expression_scale':[-1.5,1.5],'C_F_same_effect_scale':[-1.6,1.6],'C_F_common_display_normalization':'joint TMM normalization of 15 sample/subset libraries; within-sex gene standardization across individual sample pseudobulks, then equal-sample group means','C_F_inferential_models':'unchanged paired HET and unpaired male edgeR fits','missing_effects':{'female':[g for g,q in zip(genes,fef) if not np.isfinite(q)],'male':[g for g,q in zip(genes,mef) if not np.isfinite(q)]},'exploratory_candidates_not_FDR_significant':True},indent=2))
with open(str(OUT/'mosaic_network_gene_values.tsv'),'w') as f:
 wr=csv.writer(f,delimiter='\t');wr.writerow(['gene','WT_female','HET_WT_like','HET_KO_like','WT_male','KO_male','female_paired_logFC','male_logFC'])
 for i,g in enumerate(genes):wr.writerow([g]+fm[i].tolist()+mm[i].tolist()+[fef[i],mef[i]])
(OUT/'MOSAIC_NETWORK_LEGEND.txt').write_text('A: broad cell annotations plus five single-colour insets on the same joint embedding: Female, Male, WT, HET and KO. Each inset has gray background context and only its named group highlighted. Display is balanced to 6000 QC cells per sample; all views use the same points and limits. Unresolved identities remain gray. B/E: same exon-1 readout (raw probe A UMI/total cell UMI*10000), horizontal linear 0-8.2 display, identical KDE bandwidth 0.18, quartile strips and unconnected mean points. Detected-cell fractions are annotated; all cells including zeros determine the violin shapes and statistics. Female WT and HET WT-like neurons are both A-positive; HET KO-like are A-negative/downstream-positive, so their zero exon1 signal is imposed by selection. Known male WT and KO neurons retain all cells, without probe gating. C/F: the same 14 genes, hierarchy, gene-row order, expression colour limits and effect-axis limits. Gene set is the intersection of genes tested in both contrasts from the union of nine previous female and nine previous male candidates; four female-untested genes are removed mechanically, without selecting effect directions. Joint TMM logCPM uses all 15 sample/subset libraries for display; group means weight samples equally and heatmaps standardize each gene across individual sample pseudobulks within sex before averaging by group. A single linkage clusters the concatenated five within-sex group profiles. Sidebars retain the existing paired HET and unpaired WT/KO male edgeR contrasts; they are not refitted by visual styling. D: module-based source-to-receiver expression potential; see the module methods below. Statistical findings remain exploratory and do not pass gene-wide FDR correction. The six-panel layout is prepared at 182.9 x 152.4 mm, with predominantly 5-7 pt editable text and native PDF/SVG exports.\n')
print('MATCHED_PUBLICATION_FIGURE_COMPLETE',flush=True)

with open(str(OUT/'MOSAIC_NETWORK_LEGEND.txt'),'a') as f:f.write('\nD supersedes the module averages: individual curated interactions, grouped by module. Every subunit must have >=1 CPM in all three HET samples in both probe states. Complex expression uses geometric means. Per-sample interaction potential=sqrt(source ligand CPM * receiver receptor CPM); displayed means weight the three samples equally, normalized separately for each interaction to its WT-like-to-WT-like mean. Candidates require consistent KO-like/WT-like within-state direction in all three samples. Display selects the three modules with largest absolute candidate log2 changes and up to two largest candidate interactions per module. This is exploratory effect-based selection, not a significance test; source audit contains every interaction passing expression support. WT/KO in D abbreviate WT-like/KO-like within female HET. No claim of functional communication.\n')
print('SELECTED_INTERACTIONS',names,flush=True)

with open(str(OUT/'MOSAIC_NETWORK_LEGEND.txt'),'a') as f:f.write('\nG: Directed cell-type networks for WT female (A-positive reference), HET WT-like and HET KO-like populations. Nodes are five provisional neural cell-type annotations from separate consensus clusters, with fixed size and position. Arrow colour denotes source type; width is the number of distinct curated ligand-receptor IDs whose ligand and receptor subunits are all >=1 CPM in the same >=2/3 samples. Self-loops denote within-cell-type potential. All groups use identical thresholds and linear width scales. Pcdh19 excluded. These are expression-supported candidate interactions, not significance-tested CellChat probabilities or demonstrated communication. Source counts, per-pair support, sensitivity thresholds and library metadata accompany the artwork.\n')

vp=OUT/'mosaic_network_validation.json';validation=json.loads(vp.read_text());validation.update(figure_size_inches=[7.2,7.65],network_groups=network_groups,network_cell_types=cts,network_max_pair_count=max_count,network_width_scale='common linear 0 to 2.8 pt',network_analysis='CELLTYPE_NETWORK_METHODS.json');vp.write_text(json.dumps(validation,indent=2))
