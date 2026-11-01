from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import csv,gzip,json,copy
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.text import Text
from matplotlib.patches import FancyArrowPatch, Circle, Rectangle
from matplotlib.path import Path as MplPath
from scipy.cluster.hierarchy import linkage,dendrogram
ROOT=Path(__file__).resolve().parents[1];OUT=RUNTIME_ROOT/'results/female_male_signaling_20261031';FEMALE=RUNTIME_ROOT/'results/neuronal_matched_selection_20261031';DEEP=RUNTIME_ROOT/'results/mosaic_deep_dive_20261031';SIG=RUNTIME_ROOT/'results/paired_HET_signaling_20261031'
def read(p):
 with (gzip.open(str(p),'rt') if str(p).endswith('.gz') else open(str(p))) as f:return list(csv.DictReader(f,delimiter='\t'))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':6,'axes.titlesize':6.5,'axes.labelsize':6,'xtick.labelsize':5.5,'ytick.labelsize':5.5,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.5,'xtick.major.width':.5,'ytick.major.width':.5})
BLUE='#326C9D';RED='#B64B50';INK='#26292E';GRAY='#BFC4CC'
fig=plt.figure(figsize=(7.2,6.0))
rects={'a':[.055,.73,.62,.20],'b':[.13,.46,.21,.16],'c':[.42,.385,.555,.235],'d':[.75,.73,.225,.20],'e':[.13,.15,.21,.16],'f':[.42,.075,.555,.235]}

def title(k,text):
 x,y,w,h=rects[k];offset=.005 if k in ['e','f'] else .020;fig.text(x-.027,y+h+offset,k,fontsize=8,weight='bold',color=INK);fig.text(x,y+h+offset,text,fontsize=6.5,color=INK)
def clean(ax):ax.tick_params(length=2,pad=2)
def enlarge_text(target):
 for artist in target.findobj(match=Text):artist.set_fontsize(max(6.5,artist.get_fontsize()*1.2))
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
fcache=np.load(str(FEMALE/'neuronal_exon1_signal.npz'));mcache=np.load(str(OUT/'male_total_Pcdh19_display.npz'));display_fractions={}
def violin(k,data,group,labels,text):
 upper=8.2 if k=='b' else max(8.2,float(np.percentile(data,99))*1.05)
 ax=fig.add_axes(rects[k]);ds=[data[group==i] for i in range(len(labels))];non=[i for i,d in enumerate(ds) if np.ptp(d)>0];vp=ax.violinplot([ds[i] for i in non],positions=non,vert=False,widths=.64,showextrema=False,points=180,bw_method=.18)
 for body in vp['bodies']:body.set_facecolor('#C4C8CE');body.set_edgecolor('none');body.set_alpha(1);body.set_linewidth(0)
 for i,d in enumerate(ds):
  if i not in non:ax.scatter([0],[i],s=10,c=INK,marker='D',zorder=5,clip_on=False)
  q1,med,q3=np.percentile(d,[25,50,75]);ax.plot([q1,q3],[i,i],color=INK,lw=1.2,zorder=4);ax.scatter([d.mean()],[i],s=6,c=INK,zorder=5)
  ax.text(upper*.96,i-.29,'{:.1f}% detected'.format(100*(d>0).mean()),ha='right',va='bottom',fontsize=5,color='#666D75')
 ax.set_xlim(-.10,upper);ax.set_xticks([0,4,8] if k=='b' else [0,round(upper/2,1),round(upper,1)]);ax.set_ylim(len(ds)-.45,-.55);ax.set_yticks(range(len(ds)));ax.set_yticklabels(labels,fontsize=5.4);ax.set_xlabel('UMI / 10,000',fontsize=5.5,labelpad=2);clean(ax);ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0);title(k,text);display_fractions[k]=[float((d>upper).mean()) for d in ds]
violin('b',fcache['exon1_CP10K'],fcache['group'],['WT female','HET WT-like','HET KO-like'],'Female · Pcdh19 exon 1')
violin('e',mcache['CP10K'],mcache['group'],['WT males','KO males'],'Male · total Pcdh19')
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
enlarge_text(detail)
for ext in ['png','pdf','svg']:detail.savefig(str(OUT/('HET_signaling_interactions.'+ext)),dpi=450,facecolor='white')
# One compact shared key; black labels remain readable independently of colour.
cm=plt.get_cmap('Purples');legend=fig.add_axes([.59,.020,.12,.010]);gradient=np.tile(np.linspace(-1.5,1.5,100),(2,1));legend.imshow(gradient,aspect='auto',cmap=cm,vmin=-1.5,vmax=1.5);legend.axis('off');fig.text(.59,.006,'Low',fontsize=5);fig.text(.71,.006,'High',fontsize=5,ha='right');fig.text(.75,.022,'Expression scaled within sex',fontsize=5.5)
fig.add_artist(plt.Rectangle((.06,.020),.012,.010,transform=fig.transFigure,fc=BLUE,ec='none'));fig.text(.078,.020,'Lower in KO-like / KO',fontsize=5.3);fig.add_artist(plt.Rectangle((.29,.020),.012,.010,transform=fig.transFigure,fc=RED,ec='none'));fig.text(.308,.020,'Higher in KO-like / KO',fontsize=5.3)

# G and H share exactly the same five selected interactions.
chosen=read(OUT/'HET_replicated_interaction_selection.tsv');assert len(chosen)==5
cts=['C0002','C0003','C0006','C0011'];network_groups=['KO-like vs WT-like'];max_count=len(chosen)
effects=np.array([(2**float(r['median_log2_change'])-1)*100 for r in chosen]);assert all(r['source'] in cts and r['receiver'] in cts for r in chosen)
original_height=6.;added_height=1.90;combined_height=original_height+added_height
fig.set_size_inches(7.2,combined_height)
for oldax in list(fig.axes):
 pp=oldax.get_position();oldax.set_position([pp.x0,(pp.y0*original_height+added_height)/combined_height,pp.width,pp.height*original_height/combined_height])
for tx in fig.texts:
 xx,yy=tx.get_position();tx.set_position((xx,(yy*original_height+added_height)/combined_height))
for artist in fig.artists:
 if isinstance(artist,Rectangle):artist.set_y((artist.get_y()*original_height+added_height)/combined_height);artist.set_height(artist.get_height()*original_height/combined_height)
def candidate_network(target,box,fs=5.4):
 ax=target.add_axes(box);ax.set_xlim(-1.8,1.8);ax.set_ylim(-1.65,1.65);ax.axis('off')
 positions={'C0002':(-.70,.60),'C0003':(.70,.60),'C0006':(-.70,-.65),'C0011':(.70,-.65)}
 colours={'C0002':'#7186A6','C0003':'#CE9667','C0006':'#4D928A','C0011':'#9B789D'}
 nodewords={'C0002':'RG-like','C0003':'MGE IPCs','C0006':'Inhibitory\nneurons','C0011':'Cycling\nneuroblasts'}
 paths={
 'PCDHB18_PCDHB18':[(-.83,.72),(-1.36,1.37),(.0,1.37),(-.57,.72)],
 'PCDHB21_PCDHB21':[(.87,-.52),(1.52,.0),(1.52,-1.27),(.87,-.78)],
 'SEMA4A_NRP1_PLXNA4':[(-.48,.60),(-.12,.60),(.12,.60),(.48,.60)],
 'EFNA4_EPHA2':[(-.48,-.65),(-.12,-.65),(.12,-.65),(.48,-.65)],
 'NCAM2_L1CAM':[(-.86,-.53),(-1.05,-.24),(-1.05,.24),(-.86,.48)]}
 textpositions={'PCDHB18_PCDHB18':(-.70,1.50),'PCDHB21_PCDHB21':(1.34,.02),'SEMA4A_NRP1_PLXNA4':(.48,1.05),'EFNA4_EPHA2':(0,-.35),'NCAM2_L1CAM':(-1.43,0)}
 edge_artists=[];network_text=[]
 for r,effect in zip(chosen,effects):
  ident=r['interaction'];col=BLUE if effect<0 else RED;pts=paths[ident];assert np.linalg.norm(np.array(pts[0])-positions[r['source']])<.3 and np.linalg.norm(np.array(pts[-1])-positions[r['receiver']])<.3
  path=MplPath(pts,[MplPath.MOVETO,MplPath.CURVE4,MplPath.CURVE4,MplPath.CURVE4]);arrow=FancyArrowPatch(path=path,arrowstyle='-|>',mutation_scale=8,lw=.6+2.2*abs(effect)/max(abs(effects)),color=col,zorder=2,capstyle='round',joinstyle='round');ax.add_patch(arrow);edge_artists.append(arrow)
  ligand=r['ligand'];receptor=r['receptor'].replace('NRP1_PLXNA4','Nrp1/Plxna4');label=ligand if ligand==receptor else ligand+'–'+receptor
  if ident=='NCAM2_L1CAM':label='Ncam2–\nL1cam'
  if ident=='SEMA4A_NRP1_PLXNA4':label='Sema4a–Nrp1/\nPlxna4'
  label+=('\n' if ident=='PCDHB21_PCDHB21' else ' ')+'{:+.0f}%'.format(effect)
  xx,yy=textpositions[ident];network_text.append(ax.text(xx,yy,label,ha='center',va='center',fontsize=fs-.3,color=INK,zorder=5))
 for ct,(xx,yy) in positions.items():
  ax.scatter([xx],[yy],s=55,c=colours[ct],edgecolors='white',linewidths=.7,zorder=4)
  label_y=.25 if yy>0 else -1.13
  network_text.append(ax.text(xx,label_y,nodewords[ct],ha='center',va='top',fontsize=fs,color=INK,zorder=5))
 ax._network_edges=edge_artists;ax._network_text=network_text
 return ax
fig.text(.028,1.72/combined_height,'g',fontsize=8,weight='bold');fig.text(.055,1.72/combined_height,'HET candidate network',fontsize=6.5)
main_network=candidate_network(fig,[.045,.10/combined_height,.37,1.50/combined_height])
fig.text(.435,1.72/combined_height,'h',fontsize=8,weight='bold');fig.text(.46,1.72/combined_height,'KO-like vs WT-like',fontsize=6.5)
labelsax=fig.add_axes([.46,.58/combined_height,.25,.92/combined_height]);labelsax.set_xlim(0,1);labelsax.set_ylim(4.6,-.6);labelsax.axis('off')
barax=fig.add_axes([.69,.58/combined_height,.285,.92/combined_height]);barax.barh(range(5),effects,height=.58,color=[BLUE if v<0 else RED for v in effects]);barax.axvline(0,color='#9FA5AD',lw=.5);barax.set_ylim(4.6,-.6);barax.set_xlim(-100,110);barax.set_xticks([-100,0,100]);barax.set_xticklabels(['−100%','0','+100%'],fontsize=5);barax.set_yticks([]);barax.set_xlabel('Change in inferred interaction strength (%)',fontsize=5.3,labelpad=2);barax.spines['left'].set_visible(False);clean(barax)
for i,r in enumerate(chosen):
 l=r['ligand'];re=r['receptor'].replace('NRP1_PLXNA4','Nrp1/Plxna4');label=l+' contact' if l==re else l+'–'+re
 labelsax.text(0,i,label,fontsize=5.8,va='center');barax.text(effects[i]+(3 if effects[i]>0 else -3),i,'{:+.0f}%'.format(effects[i]),fontsize=5,ha='left' if effects[i]>0 else 'right',va='center')
fig.add_artist(Rectangle((.46,.10/combined_height),.012,.075/combined_height,transform=fig.transFigure,fc=BLUE,ec='none'));fig.text(.480,.10/combined_height,'Lower in KO-like',fontsize=5.5)
fig.add_artist(Rectangle((.70,.10/combined_height),.012,.075/combined_height,transform=fig.transFigure,fc=RED,ec='none'));fig.text(.720,.10/combined_height,'Higher in KO-like',fontsize=5.5)
# A standalone version contains the identical graph and labels.
nfig=plt.figure(figsize=(4.2,3.1));nfig.text(.06,.95,'HET KO-like vs WT-like · candidate network',fontsize=9,weight='bold');candidate_network(nfig,[.04,.10,.92,.77],7);enlarge_text(nfig)
for ext in ['png','pdf','svg']:nfig.savefig(str(OUT/('HET_visual_clarity_network.'+ext)),dpi=450,facecolor='white')
with open(str(OUT/'HET_H_candidate_network_edges.tsv'),'w') as f:
 fields=['interaction','source','receiver','source_name','receiver_name','median_percent_change'];w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader()
 for r,v in zip(chosen,effects):w.writerow({key:(float(v) if key=='median_percent_change' else r[key]) for key in fields})
enlarge_text(fig)
for ext in ['png','pdf','svg']:fig.savefig(str(OUT/('HET_visual_clarity_figure.'+ext)),dpi=450,facecolor='white')
(OUT/'HET_visual_clarity_validation.json').write_text(json.dumps({'figure_size_inches':[7.2,6.0],'atlas_views':['broad cell identity','Female inset','Male inset','WT inset','HET inset','KO inset'],'atlas_same_embedding':True,'atlas_render_order':'one fixed random interleaving, shared by all views; no sex/genotype drawn last as a block','atlas_cells':len(cache),'short_violin_height_inches':rects['b'][3]*6,'violin_orientation':'horizontal; no connecting mean lines; quartile strip and mean point; detection fraction displayed','tree_rendering':'one shared coordinate system; complete linkage segments meet the heatmap edge; gene labels to the right of the matrix','directional_panel':'curated neural signaling modules across all four HET state routes; expression potential normalized per module to WT-like within-state; not measured communication','module_names':selected_modules,'interaction_names':names,'selection':'all subunits >=1 CPM in every HET sample in both states; KO/WT within-state direction consistent in three samples; top three modules by largest absolute change, top two interactions per module','inset_palette':'separate from broad cell identity palette','genes_dropped_for_missing_female_contrast':dropped,'body_font_pt':[5.1,6.5],'B_E_measurements':{'B':'exon1 UMI / total UMI * 10000','E':'total Pcdh19 UMI / total UMI * 10000, no probe gate'},'B_x_range':[-.10,8.2],'E_x_range':'99th percentile times 1.05, minimum 8.2','B_E_same_KDE_bandwidth':.18,'B_E_above_range_fractions':display_fractions,'male_group_selection':'known genotype, all neurons; no probe gate','female_group_selection':'matched A-positive WT and HET WT-like; HET KO-like A-negative/downstream-positive','C_F_same_genes_and_order':orders['c'],'C_F_same_expression_scale':[-1.5,1.5],'C_F_same_effect_scale':[-1.6,1.6],'C_F_common_display_normalization':'joint TMM normalization of 15 sample/subset libraries; within-sex gene standardization across individual sample pseudobulks, then equal-sample group means','C_F_inferential_models':'unchanged paired HET and unpaired male edgeR fits','missing_effects':{'female':[g for g,q in zip(genes,fef) if not np.isfinite(q)],'male':[g for g,q in zip(genes,mef) if not np.isfinite(q)]},'exploratory_candidates_not_FDR_significant':True},indent=2))
with open(str(OUT/'HET_visual_clarity_gene_values.tsv'),'w') as f:
 wr=csv.writer(f,delimiter='\t');wr.writerow(['gene','WT_female','HET_WT_like','HET_KO_like','WT_male','KO_male','female_paired_logFC','male_logFC'])
 for i,g in enumerate(genes):wr.writerow([g]+fm[i].tolist()+mm[i].tolist()+[fef[i],mef[i]])
(OUT/'HET_VISUAL_CLARITY_LEGEND.txt').write_text('A: broad cell annotations plus five single-colour insets on the same joint embedding: Female, Male, WT, HET and KO. Each inset has gray background context and only its named group highlighted. Display is balanced to 6000 QC cells per sample; all views use the same points and limits. Unresolved identities remain gray. B: exon-1 readout. E: TOTAL Pcdh19 (A+B+C) in ALL QC male C0006 neurons, including zeros, with no expression-based gating. Both use UMI/total cell UMI*10000, identical KDE bandwidth 0.18, quartile strips and unconnected mean points; E has its own 99th-percentile display limit. Detected-cell fractions are annotated; all cells including zeros determine the violin shapes and statistics. Female WT and HET WT-like neurons are both A-positive; HET KO-like are A-negative/downstream-positive, so their zero exon1 signal is imposed by selection. Known male WT and KO neurons retain all cells, without probe gating. C/F: the same 14 genes, hierarchy, gene-row order, expression colour limits and effect-axis limits. Gene set is the intersection of genes tested in both contrasts from the union of nine previous female and nine previous male candidates; four female-untested genes are removed mechanically, without selecting effect directions. Joint TMM logCPM uses all 15 sample/subset libraries for display; group means weight samples equally and heatmaps standardize each gene across individual sample pseudobulks within sex before averaging by group. A single linkage clusters the concatenated five within-sex group profiles. Sidebars retain the existing paired HET and unpaired WT/KO male edgeR contrasts; they are not refitted by visual styling. D: module-based source-to-receiver expression potential; see the module methods below. Statistical findings remain exploratory and do not pass gene-wide FDR correction. The six-panel layout is prepared at 182.9 x 152.4 mm, with predominantly 5-7 pt editable text and native PDF/SVG exports.\n')
print('MATCHED_PUBLICATION_FIGURE_COMPLETE',flush=True)

with open(str(OUT/'HET_VISUAL_CLARITY_LEGEND.txt'),'a') as f:f.write('\nD supersedes the module averages: individual curated interactions, grouped by module. Every subunit must have >=1 CPM in all three HET samples in both probe states. Complex expression uses geometric means. Per-sample interaction potential=sqrt(source ligand CPM * receiver receptor CPM); displayed means weight the three samples equally, normalized separately for each interaction to its WT-like-to-WT-like mean. Candidates require consistent KO-like/WT-like within-state direction in all three samples. Display selects the three modules with largest absolute candidate log2 changes and up to two largest candidate interactions per module. This is exploratory effect-based selection, not a significance test; source audit contains every interaction passing expression support. WT/KO in D abbreviate WT-like/KO-like within female HET. No claim of functional communication.\n')
print('SELECTED_INTERACTIONS',names,flush=True)

with open(str(OUT/'HET_VISUAL_CLARITY_LEGEND.txt'),'a') as f:f.write('\nG: One schematic network containing exactly the five interaction IDs displayed in H, without additional edges or cell types. Nodes denote the four source/receiver cell-type annotations represented in H. Each arrow is one selected ligand-receptor interaction; its direction denotes ligand source to receptor receiver. Within-type loops refer to populations, not proven autocrine signaling. Blue/red denote lower/higher KO-like-to-KO-like expression potential relative to WT-like-to-WT-like within female HET; arrow width is 0.6 + 2.2 * absolute median percent change / maximum absolute median percent change. Labels reproduce the same median percent changes as H. The graph is a schematic of selected expression-based candidates, not a network-wide interaction-count analysis or a significance-tested communication model. No depth matching.\n')

vp=OUT/'HET_visual_clarity_validation.json';validation=json.loads(vp.read_text());validation.update(figure_size_inches=[7.2,combined_height],network_groups=network_groups,network_cell_types=cts,network_interactions=[r['interaction'] for r in chosen],network_width_scale='0.6 + 2.2 * abs(median percent change) / max(abs(percent change)) pt',network_analysis='HET_H_candidate_network_edges.tsv',male_expression='total Pcdh19, all QC male C0006 neurons, no expression gating');vp.write_text(json.dumps(validation,indent=2))

vp=OUT/'HET_visual_clarity_validation.json';vv=json.loads(vp.read_text());vv.update(figure_size_inches=[7.2,combined_height],G='one network, exactly the five selected edges in H, no depth matching',H='five continuous interaction candidates from original CPMs, median KO-like within-state vs WT-like within-state change',H_same_direction_samples=3,H_no_depth_matching=True);vp.write_text(json.dumps(vv,indent=2))
with open(str(OUT/'HET_VISUAL_CLARITY_LEGEND.txt'),'a') as f:f.write('\nH: Individual curated interaction-potential candidates computed from original CPMs, without depth matching. Potential=sqrt(source ligand expression times receiver receptor expression), geometric means for complex subunits. Displayed percentages are median paired KO-like-to-KO-like versus WT-like-to-WT-like changes over three HET samples. All components >=1 CPM in all three samples in both states; consistent direction and absolute log2 change >=0.2 in every sample. Rank by the smallest absolute paired effect, retain unique interaction IDs, up to two protocadherin interactions and one per other predefined neural family, display five. Selected exploratory candidates, not FDR-confirmed signaling. G contains exactly these five interactions, with widths scaled to absolute median percent changes; no depth matching in this figure.\n')
print('NETWORK_AND_CANDIDATES_COMPLETE',flush=True)

vp=OUT/'HET_visual_clarity_validation.json';vv=json.loads(vp.read_text());vv.update(typography='all Text artists increased by 20%, minimum 6.5 pt',H_axis_label='Change in inferred interaction strength (%)',H_bar_axis_left=.69,G_text_backgrounds=False,G_arrow_style='rounded cubic curves, consistent arrowheads',H_color_key='blue lower in KO-like, red higher in KO-like');vp.write_text(json.dumps(vv,indent=2))

fig.canvas.draw();renderer=fig.canvas.get_renderer();collisions=[]
for ti,tx in enumerate(main_network._network_text):
 box=tx.get_window_extent(renderer).padded(2.5)
 for ei,edge in enumerate(main_network._network_edges):
  path=edge.get_path().transformed(edge.get_transform())
  if path.intersects_bbox(box,filled=False):collisions.append(dict(text=tx.get_text(),edge_index=ei))
vp=OUT/'HET_visual_clarity_validation.json';vv=json.loads(vp.read_text());vv['G_text_arrow_collisions_with_2_5px_clearance']=collisions;vp.write_text(json.dumps(vv,indent=2));print('NETWORK_LABEL_CLEARANCE',collisions,flush=True)
