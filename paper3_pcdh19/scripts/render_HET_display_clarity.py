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
rects={'a':[.055,.73,.62,.20],'b':[.13,.41,.21,.20],'c':[.385,.385,.445,.235],'d':[.75,.73,.225,.20],'e':[.13,.10,.21,.20],'f':[.385,.075,.445,.235]}

def figtext(row,*args,**kwargs):
 artist=fig.text(*args,**kwargs);artist._figure_row=row;return artist
def rowpatch(row,artist):
 artist._figure_row=row;fig.add_artist(artist);return artist
def panel_heading(row,letter,text):
 left=letter in ['A','C','E','G'];xx=.025 if left else .435
 art=figtext(row,xx,.95,letter,fontsize=11,weight='bold',color=INK);art._panel_heading=True
 art=figtext(row,.060 if left else .470,.95,text,fontsize=7.5,weight='bold',color=INK);art._panel_heading=True
def title(k,text):
 row={'a':0,'d':0,'b':1,'c':1,'e':2,'f':2}[k];letter={'a':'A','d':'B','b':'C','c':'D','e':'E','f':'F'}[k];panel_heading(row,letter,text)

def clean(ax):ax.tick_params(length=2,pad=2)
def enlarge_text(target):
 for ax in target.axes:
  for axis in [ax.xaxis,ax.yaxis]:
   ticks=axis.get_ticklabels();labels=[t.get_text() for t in ticks]
   if any(any(old in lab for old in ['KO-like','WT-like','RG-like']) for lab in labels):
    axis.set_ticklabels([lab.replace('KO-like','KO cells').replace('WT-like','WT cells').replace('RG-like progenitors','RG').replace('RG-like','RG') for lab in labels])
 for artist in target.findobj(match=Text):
  label=artist.get_text().replace('KO-like','KO cells').replace('WT-like','WT cells').replace('RG-like progenitors','RG').replace('RG-like','RG')
  if label=='KO / WT':label='KO males / WT males'
  artist.set_text(label);artist.set_fontsize(max(6.5,artist.get_fontsize()*1.2))
# a: contextual atlas with broad, provenance-supported labels.
title('a','Mouse MGE · E14.5')
figtext(0,.060,.91,'12 samples from 12 mice · dissected MGE',fontsize=6,color=INK)
figtext(0,.055,.535,'WT female: 3 · HET female: 3 · WT male: 3 · KO male: 3',fontsize=5,color=INK)
broad=[('RG-like progenitors',['C0002'],'#216869'),('MGE IPCs',['C0003','C0008'],'#B99A38'),('Neurogenic progenitors',['C0001','C0007','C0012','C0024'],'#C57750'),('Cycling neuroblasts',['C0011'],'#8F799E'),('Immature inhibitory neurons',['C0004','C0005','C0006','C0020'],'#7A883E'),('Mixed neural / vascular',['C0015'],'#8B6C5C'),('Unresolved populations',[],GRAY)]
lookup={c:i for i,(_,cc,_) in enumerate(broad) for c in cc};cache=read(OUT/'atlas_sex_genotype_display.tsv.gz');xy=np.array([[float(r['UMAP1']),float(r['UMAP2'])] for r in cache]);category=np.array([lookup.get(r['cluster'],6) for r in cache]);sex=np.array([r['sex'] for r in cache]);geno=np.array([r['genotype'] for r in cache]);lim=[xy[:,0].min()-.4,xy[:,0].max()+.4,xy[:,1].min()-.4,xy[:,1].max()+.4]
po=np.random.RandomState(20261031).permutation(len(cache))
ax=fig.add_axes([.055,.680,.22,.22]);ax._figure_row=0;cols=np.array([broad[i][2] for i in category]);ax.scatter(xy[po,0],xy[po,1],s=.22,c=cols[po],alpha=.85,linewidths=0,rasterized=True);ax.set_xlim(lim[:2]);ax.set_ylim(lim[2:]);ax.set_aspect('equal');ax.axis('off')
insets=[('Female',sex=='Female','#0072B2',.055,.570),('Male',sex=='Male','#E69F00',.135,.570),('WT',geno=='WT','#222222',.215,.570),('HET',geno=='HET','#CC00AA',.295,.570),('KO',geno=='KO','#56B4E9',.375,.570)]
for name,mask,col,xx,yy in insets:
 ax=fig.add_axes([xx,yy,.065,.075]);ax._figure_row=0;ax.scatter(xy[po,0],xy[po,1],s=.025,c='#E6E7E9',linewidths=0,rasterized=True);ax.scatter(xy[mask,0],xy[mask,1],s=.045,c=col,alpha=.70,linewidths=0,rasterized=True);ax.set_xlim(lim[:2]);ax.set_ylim(lim[2:]);ax.set_aspect('equal');ax.axis('off');ax.set_title(name,fontsize=5.5,pad=1)
lg=fig.add_axes([.280,.680,.16,.22]);lg._figure_row=0;lg.axis('off');hh=[Line2D([0],[0],marker='o',ls='',color=col,label=name.replace('Immature inhibitory neurons','Immature inhibitory\nneurons'),markersize=2.5) for name,cc,col in broad];lg.legend(handles=hh,ncol=1,frameon=False,fontsize=5,loc='upper left',handletextpad=.35,handlelength=.5,labelspacing=.15,borderaxespad=0)
title('a','')
# b/e: identical cell-level metric, axis range, violin settings and mean markers.
fcache=np.load(str(FEMALE/'neuronal_exon1_signal.npz'));mcache=np.load(str(OUT/'male_total_Pcdh19_display.npz'));display_fractions={}
def violin(k,data,group,labels,text):
 data=np.log1p(data);upper=3.0
 ax=fig.add_axes(rects[k]);ax._figure_row=1 if k=='b' else 2;ds=[data[group==i] for i in range(len(labels))];non=[i for i,d in enumerate(ds) if np.ptp(d)>0];vp=ax.violinplot([ds[i] for i in non],positions=non,vert=False,widths=.64,showextrema=False,points=180,bw_method=.18)
 for body in vp['bodies']:body.set_facecolor('#C4C8CE');body.set_edgecolor('none');body.set_alpha(1);body.set_linewidth(0)
 for i,d in enumerate(ds):
  if i not in non:
   ax.plot([0,0],[i-.18,i+.18],color=INK,lw=1.8,zorder=5,clip_on=False)
   ax.text(.16,i,'0 exon-1 signal' if k=='b' else '0 signal',ha='left',va='center',fontsize=5.3,color=INK)
  q1,med,q3=np.percentile(d,[25,50,75]);ax.plot([q1,q3],[i,i],color=INK,lw=1.2,zorder=4);ax.scatter([d.mean()],[i],s=6,c=INK,zorder=5)

 ax.set_xlim(-.15,upper);ax.set_xticks([0,1,2,3]);ax.set_ylim(len(ds)-.45,-.55);ax.set_yticks(range(len(ds)));ax.set_yticklabels(labels,fontsize=5.4);ax.set_xlabel('Log-normalized expression',fontsize=5.5,labelpad=2);clean(ax);ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0);title(k,text);display_fractions[k]=[float((d>upper).mean()) for d in ds]
violin('b',fcache['exon1_CP10K'],fcache['group'],['WT female','HET WT-like','HET KO-like'],'Female · Pcdh19 exon 1')
male_original_group=mcache['group'];male_positive=mcache['raw_counts']>0;male_display_mask=(male_original_group!=0)|male_positive
assert np.array_equal(male_positive,mcache['CP10K']>0)
violin('e',mcache['CP10K'][male_display_mask],male_original_group[male_display_mask],['WT males\nPcdh19+','KO males'],'Male · total Pcdh19')
# A common pre-existing candidate panel makes sex comparisons direct.
fgenes=['Ank3','Kcnc2','Gnao1','Ptprd','Epha5','Snap47','Pcdha5','Pcdha10','Pcdhga5'];mgenes=['Eps8l1','Serpine2','Syt9','Ecel1','Lrfn5','Grm7','Cntn1','Sphkap','Nrp2'];genes=fgenes+mgenes
common_female={r['gene_symbol'] for r in read(DEEP/'DE_C0006/within_HET_Aneg_BC1_vs_Apositive.tsv.gz')};common_male={r['gene_symbol'] for r in read(OUT/'male_C0006_DE.tsv')};dropped=[g for g in genes if g not in common_female or g not in common_male];genes=[g for g in genes if g in common_female and g in common_male]
rows=read(OUT/'matched_display_logCPM.tsv');values={r['gene_symbol']:r for r in rows};md=read(OUT/'matched_female_display_metadata.tsv');go=[('WT_F','A_positive'),('HET_F','A_positive'),('HET_F','Aneg_BC1')];md=sorted(md,key=lambda r:(go.index((r['design_group'],r['gate'])),r['sample']));flibs=[r['library'] for r in md];mlibs=['15662-JZ-'+str(i)+'__C0006' for i in [1,2,3,10,11,12]]
fm=np.array([[np.mean([float(values[g][l]) for l in flibs[k:k+3]]) for k in [0,3,6]] for g in genes]);mm=np.array([[np.mean([float(values[g][l]) for l in mlibs[k:k+3]]) for k in [0,3]] for g in genes])
def standardize(a):return (a-a.mean(1)[:,None])/np.maximum(a.std(1)[:,None],1e-9)
fx=np.array([[float(values[g][l]) for l in flibs] for g in genes]);mx=np.array([[float(values[g][l]) for l in mlibs] for g in genes])
fz=(fm-fx.mean(1)[:,None])/np.maximum(fx.std(1)[:,None],1e-9);mz=(mm-mx.mean(1)[:,None])/np.maximum(mx.std(1)[:,None],1e-9);link=linkage(np.column_stack([fz,mz]),method='average');fd={r['gene_symbol']:r for r in read(DEEP/'DE_C0006/within_HET_Aneg_BC1_vs_Apositive.tsv.gz')};mdict={r['gene_symbol']:r for r in read(OUT/'male_C0006_DE.tsv')};fef=np.array([float(fd[g]['logFC']) if g in fd else np.nan for g in genes]);mef=np.array([float(mdict[g]['logFC']) if g in mdict else np.nan for g in genes]);orders={}
def matched_gene_block(k,z,effect,labels,text,contrast):
 x,y,w,h=rects[k];ax=fig.add_axes([x,y,w,h]);ax._figure_row=1 if k=='c' else 2;dd=dendrogram(link,no_plot=True);order=dd['leaves'];md=max(float(link[-1,2]),1e-9);expr_width=2.3;base=4.25;scale=.48
 # One coordinate system keeps every branch connected and places the leaves on the heatmap boundary.
 for ic,dc in zip(dd['icoord'],dd['dcoord']):ax.plot(-.8*np.asarray(dc)/md,(np.asarray(ic)-5)/10,color='#555B63',lw=.6,zorder=2)
 ax.imshow(z[order],extent=[0,expr_width,len(genes)-.5,-.5],aspect='auto',cmap='Purples',vmin=-1.5,vmax=1.5,zorder=1);ax.add_patch(Rectangle((0,-.5),expr_width,len(genes),fill=False,ec='#555B63',lw=.35,zorder=3))
 for i,j in enumerate(order):ax.text(expr_width+.10,i,genes[j],ha='left',va='center',fontsize=5.6)
 ev=effect[order];ax.barh(range(len(genes)),ev*scale,left=base,height=.60,color=[BLUE if q<0 else RED for q in ev],zorder=3);ax.plot([base,base],[-.5,len(genes)-.5],color='#A1A5AC',lw=.4)
 ax.set_xlim(-.85,5.10);ax.set_ylim(len(genes)-.5,-.5);ax.axis('off')
 for j,lab in enumerate(labels):ax.text((j+.5)*expr_width/len(labels),len(genes)-.15,lab.replace(' ','\n',1),ha='center',va='top',fontsize=5.2,clip_on=False)
 for tick in [-1,0,1]:
  xx=base+tick*scale;ax.plot([xx,xx],[len(genes)-.5,len(genes)-.35],color=INK,lw=.4,clip_on=False);ax.text(xx,len(genes)-.10,str(tick),ha='center',va='top',fontsize=5.2,clip_on=False)
 ax.text(base,len(genes)+.65,'log2FC',ha='center',va='top',fontsize=5.5,clip_on=False);title(k,text);orders[k]=[genes[i] for i in order]
matched_gene_block('c',fz,fef,['WT female','HET WT-like','HET KO-like'],'Female HET · neuronal genes','KO-like / WT-like\nwithin HET females')
matched_gene_block('f',mz,mef,['WT males','KO males'],'Male WT vs KO · neuronal genes','KO / WT')
assert orders['c']==orders['f']
# B: reuse saved regional marker-program scores on the exact atlas embedding.
regional=np.load(str(OUT/'saved_MGE_LGE_CGE_scores.npz'));assert np.array_equal(regional['row'],np.array([int(r['row']) for r in cache]))
region_scores=regional['scores'];score_threshold=.8;score_limits=np.array([score_threshold,1.5]);highlighted=region_scores>score_threshold
from matplotlib.colors import LinearSegmentedColormap
regional_cmap=LinearSegmentedColormap.from_list('regional_positive',plt.get_cmap('Greens')(np.linspace(.30,.95,256)))
for j,xx in enumerate([.475,.635,.795]):
 ax=fig.add_axes([xx,.620,.145,.28]);ax._figure_row=0
 ax.scatter(xy[po,0],xy[po,1],s=.18,c='#D8DADD',linewidths=0,rasterized=True)
 selected=po[highlighted[po,j]]
 drawn=ax.scatter(xy[selected,0],xy[selected,1],s=.25,c=region_scores[selected,j],cmap=regional_cmap,vmin=score_limits[0],vmax=score_limits[1],linewidths=0,rasterized=True)
 ax.set_xlim(lim[:2]);ax.set_ylim(lim[2:]);ax.set_aspect('equal');ax.axis('off');ax.set_title(['MGE','LGE','CGE'][j],fontsize=6.5,pad=2)
scoreax=fig.add_axes([.944,.665,.005,.160]);scoreax._figure_row=0
cb=fig.colorbar(drawn,cax=scoreax,orientation='vertical',ticks=score_limits);cb.ax.set_yticklabels(['Low','High']);cb.ax.tick_params(length=0,pad=2,labelsize=5,labelcolor='#737B73');cb.outline.set_visible(False)
for yy in [0,1]:scoreax.plot([-.1,1.1],[yy,yy],transform=scoreax.transAxes,color='#737B73',lw=.35,clip_on=False)
title('d','Regional identity')
selected_modules=[];names=[]
# Expression colour key; comparison keys are placed locally beside C/F.
cm=plt.get_cmap('Purples');legend=fig.add_axes([.43,.020,.12,.010]);legend._figure_row=2;gradient=np.tile(np.linspace(-1.5,1.5,100),(2,1));legend.imshow(gradient,aspect='auto',cmap=cm,vmin=-1.5,vmax=1.5);legend.axis('off');figtext(2,.43,.006,'Low',fontsize=5);figtext(2,.55,.006,'High',fontsize=5,ha='right')
for panel in ['c','f']:
 row=1 if panel=='c' else 2
 labels=['Lower in\nHET KO cells\nvs HET WT cells','Higher in\nHET KO cells\nvs HET WT cells'] if panel=='c' else ['Lower in\nKO males\nvs WT males','Higher in\nKO males\nvs WT males']
 for offset,col,label in [(.165,BLUE,labels[0]),(.065,RED,labels[1])]:
  yy=rects[panel][1]+offset
  rowpatch(row,Rectangle((.852,yy-.004),.010,.008,transform=fig.transFigure,fc=col,ec='none'))
  figtext(row,.870,yy,label,fontsize=5,va='center',linespacing=1.25)

# G and H share exactly the same five selected interactions.
chosen=read(OUT/'HET_replicated_interaction_selection.tsv');assert len(chosen)==5
cts=['C0002','C0003','C0006','C0011'];network_groups=['KO-like vs WT-like'];max_count=len(chosen)
effects=np.array([(2**float(r['median_log2_change'])-1)*100 for r in chosen]);assert all(r['source'] in cts and r['receiver'] in cts for r in chosen)
original_height=6.;added_height=1.90;combined_height=original_height+added_height
fig.set_size_inches(7.2,combined_height)
for oldax in list(fig.axes):
 pp=oldax.get_position(original=True);oldax.set_position([pp.x0,(pp.y0*original_height+added_height)/combined_height,pp.width,pp.height*original_height/combined_height])
for tx in fig.texts:
 xx,yy=tx.get_position();tx.set_position((xx,(yy*original_height+added_height)/combined_height))
for artist in fig.artists:
 if isinstance(artist,Rectangle):artist.set_y((artist.get_y()*original_height+added_height)/combined_height);artist.set_height(artist.get_height()*original_height/combined_height)
def candidate_network(target,box,fs=5.4):
 ax=target.add_axes(box);ax._figure_row=3;ax.set_xlim(-1.8,1.8);ax.set_ylim(-1.65,1.65);ax.axis('off')
 positions={'C0002':(-.70,.60),'C0003':(.70,.60),'C0006':(-.70,-.65),'C0011':(.70,-.65)}
 colours={'C0002':'#216869','C0003':'#B99A38','C0006':'#7A883E','C0011':'#8F799E'}
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
panel_heading(3,'G','HET candidate network')
main_network=candidate_network(fig,[.045,.10/combined_height,.37,1.50/combined_height])
panel_heading(3,'H','KO-like vs WT-like within HET Female')
labelsax=fig.add_axes([.46,.58/combined_height,.25,.92/combined_height]);labelsax._figure_row=3;labelsax.set_xlim(0,1);labelsax.set_ylim(4.6,-.6);labelsax.axis('off')
barax=fig.add_axes([.69,.58/combined_height,.285,.92/combined_height]);barax._figure_row=3;barax.barh(range(5),effects,height=.58,color=[BLUE if v<0 else RED for v in effects]);barax.axvline(0,color='#9FA5AD',lw=.5);barax.set_ylim(4.6,-.6);barax.set_xlim(-100,110);barax.set_xticks([-100,0,100]);barax.set_xticklabels(['−100%','0','+100%'],fontsize=5);barax.set_yticks([]);barax.set_xlabel('Change in inferred interaction strength (%)',fontsize=5.3,labelpad=2);barax.spines['left'].set_visible(False);clean(barax)
for i,r in enumerate(chosen):
 l=r['ligand'];re=r['receptor'].replace('NRP1_PLXNA4','Nrp1/Plxna4');label=l+' contact' if l==re else l+'–'+re
 labelsax.text(0,i,label,fontsize=5.8,va='center');barax.text(effects[i]+(3 if effects[i]>0 else -3),i,'{:+.0f}%'.format(effects[i]),fontsize=5,ha='left' if effects[i]>0 else 'right',va='center')
rowpatch(3,Rectangle((.46,.10/combined_height),.012,.075/combined_height,transform=fig.transFigure,fc=BLUE,ec='none'));figtext(3,.480,.10/combined_height,'Lower in KO-like',fontsize=5.5)
rowpatch(3,Rectangle((.70,.10/combined_height),.012,.075/combined_height,transform=fig.transFigure,fc=RED,ec='none'));figtext(3,.720,.10/combined_height,'Higher in KO-like',fontsize=5.5)
# A standalone version contains the identical graph and labels.
nfig=plt.figure(figsize=(4.2,3.1));nfig.text(.06,.95,'HET KO-like vs WT-like · candidate network',fontsize=9,weight='bold');candidate_network(nfig,[.04,.10,.92,.77],7);enlarge_text(nfig)
for ext in ['png','pdf','svg']:nfig.savefig(str(OUT/('HET_display_clarity_network.'+ext)),dpi=450,facecolor='white')
with open(str(OUT/'HET_H_candidate_network_edges.tsv'),'w') as f:
 fields=['interaction','source','receiver','source_name','receiver_name','median_percent_change'];w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader()
 for r,v in zip(chosen,effects):w.writerow({key:(float(v) if key=='median_percent_change' else r[key]) for key in fields})
# Four independently spaced rows, with common left anchors and baselines.
old_height=combined_height;row_old={0:(5.26,7.78),1:(3.88,5.86),2:(1.912,3.94),3:(0.,1.85)}
row_new={3:(.12,1.85),2:(2.11,2.00),1:(4.25,1.98),0:(6.37,2.00)};combined_height=8.49
fig.set_size_inches(7.2,combined_height)
def new_y(row,y):
 lo,hi=row_old[row];bottom,height=row_new[row];return (bottom+(y*old_height-lo)*height/(hi-lo))/combined_height
for ax in fig.axes:
 row=ax._figure_row;pp=ax.get_position(original=True);lo,hi=row_old[row];bottom,height=row_new[row];ax.set_position([pp.x0,new_y(row,pp.y0),pp.width,pp.height*old_height*height/(hi-lo)/combined_height])
for tx in fig.texts:
 row=tx._figure_row;xx,yy=tx.get_position()
 if getattr(tx,'_panel_heading',False):tx.set_position((xx,(row_new[row][0]+row_new[row][1]-.15)/combined_height))
 else:tx.set_position((xx,new_y(row,yy)))
for artist in fig.artists:
 if isinstance(artist,Rectangle):
  row=artist._figure_row;lo,hi=row_old[row];bottom,height=row_new[row];artist.set_y(new_y(row,artist.get_y()));artist.set_height(artist.get_height()*old_height*height/(hi-lo)/combined_height)
enlarge_text(fig)
for ext in ['png','pdf','svg']:fig.savefig(str(OUT/('HET_display_clarity_figure.'+ext)),dpi=450,facecolor='white')
(OUT/'HET_display_clarity_validation.json').write_text(json.dumps({'figure_size_inches':[7.2,6.0],'atlas_views':['broad cell identity','Female inset','Male inset','WT inset','HET inset','KO inset'],'atlas_same_embedding':True,'atlas_render_order':'one fixed random interleaving, shared by all views; no sex/genotype drawn last as a block','atlas_cells':len(cache),'short_violin_height_inches':rects['b'][3]*6,'violin_orientation':'horizontal; no connecting mean lines; quartile strip and mean point; detection fraction displayed','tree_rendering':'one shared coordinate system; complete linkage segments meet the heatmap edge; gene labels to the right of the matrix','directional_panel':'curated neural signaling modules across all four HET state routes; expression potential normalized per module to WT-like within-state; not measured communication','module_names':selected_modules,'interaction_names':names,'selection':'all subunits >=1 CPM in every HET sample in both states; KO/WT within-state direction consistent in three samples; top three modules by largest absolute change, top two interactions per module','inset_palette':'separate from broad cell identity palette','genes_dropped_for_missing_female_contrast':dropped,'body_font_pt':[5.1,6.5],'B_E_measurements':{'B':'exon1 UMI / total UMI * 10000','E':'total Pcdh19 UMI / total UMI * 10000, no probe gate'},'B_x_range':[-.10,8.2],'E_x_range':'99th percentile times 1.05, minimum 8.2','B_E_same_KDE_bandwidth':.18,'B_E_above_range_fractions':display_fractions,'male_group_selection':'known genotype, all neurons; no probe gate','female_group_selection':'matched A-positive WT and HET WT-like; HET KO-like A-negative/downstream-positive','C_F_same_genes_and_order':orders['c'],'C_F_same_expression_scale':[-1.5,1.5],'C_F_same_effect_scale':[-1.6,1.6],'C_F_common_display_normalization':'joint TMM normalization of 15 sample/subset libraries; within-sex gene standardization across individual sample pseudobulks, then equal-sample group means','C_F_inferential_models':'unchanged paired HET and unpaired male edgeR fits','missing_effects':{'female':[g for g,q in zip(genes,fef) if not np.isfinite(q)],'male':[g for g,q in zip(genes,mef) if not np.isfinite(q)]},'exploratory_candidates_not_FDR_significant':True},indent=2))
with open(str(OUT/'HET_display_clarity_gene_values.tsv'),'w') as f:
 wr=csv.writer(f,delimiter='\t');wr.writerow(['gene','WT_female','HET_WT_like','HET_KO_like','WT_male','KO_male','female_paired_logFC','male_logFC'])
 for i,g in enumerate(genes):wr.writerow([g]+fm[i].tolist()+mm[i].tolist()+[fef[i],mef[i]])
(OUT/'HET_DISPLAY_CLARITY_LEGEND.txt').write_text('A: broad cell annotations plus five single-colour insets on the same joint embedding: Female, Male, WT, HET and KO. Each inset has gray background context and only its named group highlighted. Display is balanced to 6000 QC cells per sample; all views use the same points and limits. Unresolved identities remain gray. B: exon-1 readout. E: TOTAL Pcdh19 (A+B+C) in ALL QC male C0006 neurons, including zeros, with no expression-based gating. Both use UMI/total cell UMI*10000, identical KDE bandwidth 0.18, quartile strips and unconnected mean points; E has its own 99th-percentile display limit. Detected-cell fractions are annotated; all cells including zeros determine the violin shapes and statistics. Female WT and HET WT-like neurons are both A-positive; HET KO-like are A-negative/downstream-positive, so their zero exon1 signal is imposed by selection. Known male WT and KO neurons retain all cells, without probe gating. C/F: the same 14 genes, hierarchy, gene-row order, expression colour limits and effect-axis limits. Gene set is the intersection of genes tested in both contrasts from the union of nine previous female and nine previous male candidates; four female-untested genes are removed mechanically, without selecting effect directions. Joint TMM logCPM uses all 15 sample/subset libraries for display; group means weight samples equally and heatmaps standardize each gene across individual sample pseudobulks within sex before averaging by group. A single linkage clusters the concatenated five within-sex group profiles. Sidebars retain the existing paired HET and unpaired WT/KO male edgeR contrasts; they are not refitted by visual styling. D: module-based source-to-receiver expression potential; see the module methods below. Statistical findings remain exploratory and do not pass gene-wide FDR correction. The six-panel layout is prepared at 182.9 x 152.4 mm, with predominantly 5-7 pt editable text and native PDF/SVG exports.\n')
print('MATCHED_PUBLICATION_FIGURE_COMPLETE',flush=True)

with open(str(OUT/'HET_DISPLAY_CLARITY_LEGEND.txt'),'a') as f:f.write('\nD supersedes the module averages: individual curated interactions, grouped by module. Every subunit must have >=1 CPM in all three HET samples in both probe states. Complex expression uses geometric means. Per-sample interaction potential=sqrt(source ligand CPM * receiver receptor CPM); displayed means weight the three samples equally, normalized separately for each interaction to its WT-like-to-WT-like mean. Candidates require consistent KO-like/WT-like within-state direction in all three samples. Display selects the three modules with largest absolute candidate log2 changes and up to two largest candidate interactions per module. This is exploratory effect-based selection, not a significance test; source audit contains every interaction passing expression support. WT/KO in D abbreviate WT-like/KO-like within female HET. No claim of functional communication.\n')
print('SELECTED_INTERACTIONS',names,flush=True)

with open(str(OUT/'HET_DISPLAY_CLARITY_LEGEND.txt'),'a') as f:f.write('\nG: One schematic network containing exactly the five interaction IDs displayed in H, without additional edges or cell types. Nodes denote the four source/receiver cell-type annotations represented in H. Each arrow is one selected ligand-receptor interaction; its direction denotes ligand source to receptor receiver. Within-type loops refer to populations, not proven autocrine signaling. Blue/red denote lower/higher KO-like-to-KO-like expression potential relative to WT-like-to-WT-like within female HET; arrow width is 0.6 + 2.2 * absolute median percent change / maximum absolute median percent change. Labels reproduce the same median percent changes as H. The graph is a schematic of selected expression-based candidates, not a network-wide interaction-count analysis or a significance-tested communication model. No depth matching.\n')

vp=OUT/'HET_display_clarity_validation.json';validation=json.loads(vp.read_text());validation.update(figure_size_inches=[7.2,combined_height],network_groups=network_groups,network_cell_types=cts,network_interactions=[r['interaction'] for r in chosen],network_width_scale='0.6 + 2.2 * abs(median percent change) / max(abs(percent change)) pt',network_analysis='HET_H_candidate_network_edges.tsv',male_expression='total Pcdh19, all QC male C0006 neurons, no expression gating');vp.write_text(json.dumps(validation,indent=2))

vp=OUT/'HET_display_clarity_validation.json';vv=json.loads(vp.read_text());vv.update(figure_size_inches=[7.2,combined_height],G='one network, exactly the five selected edges in H, no depth matching',H='five continuous interaction candidates from original CPMs, median KO-like within-state vs WT-like within-state change',H_same_direction_samples=3,H_no_depth_matching=True);vp.write_text(json.dumps(vv,indent=2))
with open(str(OUT/'HET_DISPLAY_CLARITY_LEGEND.txt'),'a') as f:f.write('\nH: Individual curated interaction-potential candidates computed from original CPMs, without depth matching. Potential=sqrt(source ligand expression times receiver receptor expression), geometric means for complex subunits. Displayed percentages are median paired KO-like-to-KO-like versus WT-like-to-WT-like changes over three HET samples. All components >=1 CPM in all three samples in both states; consistent direction and absolute log2 change >=0.2 in every sample. Rank by the smallest absolute paired effect, retain unique interaction IDs, up to two protocadherin interactions and one per other predefined neural family, display five. Selected exploratory candidates, not FDR-confirmed signaling. G contains exactly these five interactions, with widths scaled to absolute median percent changes; no depth matching in this figure.\n')
print('NETWORK_AND_CANDIDATES_COMPLETE',flush=True)

vp=OUT/'HET_display_clarity_validation.json';vv=json.loads(vp.read_text());vv.update(typography='all Text artists increased by 20%, minimum 6.5 pt',H_axis_label='Change in inferred interaction strength (%)',H_bar_axis_left=.69,G_text_backgrounds=False,G_arrow_style='rounded cubic curves, consistent arrowheads',H_color_key='blue lower in KO-like, red higher in KO-like');vp.write_text(json.dumps(vv,indent=2))

fig.canvas.draw();renderer=fig.canvas.get_renderer();collisions=[]
for ti,tx in enumerate(main_network._network_text):
 box=tx.get_window_extent(renderer).padded(2.5)
 for ei,edge in enumerate(main_network._network_edges):
  path=edge.get_path().transformed(edge.get_transform())
  if path.intersects_bbox(box,filled=False):collisions.append(dict(text=tx.get_text(),edge_index=ei))
vp=OUT/'HET_display_clarity_validation.json';vv=json.loads(vp.read_text());vv['G_text_arrow_collisions_with_2_5px_clearance']=collisions;vp.write_text(json.dumps(vv,indent=2));print('NETWORK_LABEL_CLEARANCE',collisions,flush=True)

with open(str(OUT/'HET_DISPLAY_CLARITY_LEGEND.txt'),'a') as f:f.write('\nDisplay labels: WT cells and KO cells are shortened names for the previously defined WT-like and KO-like probe states within female HET tissue; relabeling does not establish cell-level DNA genotypes. RG is the shortened display label for the provisional RG-like annotation. C colour key refers to KO cells versus WT cells within female HET; F colour key refers to KO males versus WT males; G/H colour key refers to KO cells versus WT cells within female HET. Lower/higher describe gene expression in C/F and inferred interaction strength in G/H. Values, groups, definitions and fits unchanged.\n')

vp=OUT/'HET_display_clarity_validation.json';v=json.loads(vp.read_text());v.update(panel_rows=[['A','B'],['C','D'],['E','F'],['G','H']],panel_letter_font_pt=13.2,panel_heading_font_pt=9.0,row_gap_inches=.14,regional_scores='SAVED_REGIONAL_SCORES_PROVENANCE.json',removed_panel='former D signaling interaction heatmap',shared_panel_heading_baselines=True,figure_size_inches=[7.2,8.49]);vp.write_text(json.dumps(v,indent=2))
with open(str(OUT/'HET_DISPLAY_CLARITY_LEGEND.txt'),'a') as f:f.write('\nFINAL LAYOUT supersedes former letter assignments: A atlas; B saved MGE/LGE/CGE regional score overlays; C female exon1 violin; D female candidate heatmap; E male total Pcdh19 violin; F male candidate heatmap; G selected interaction network; H same interaction changes. The former signaling heatmap D is removed. B reuses original saved marker-program scores, joined by exact cell ID to the same 72000 displayed cells and original UMAP; pooled 2nd/98th-percentile score colour limits are shared across three regions. These scores provide descriptive regional identity support rather than proving dissection boundaries. Four rows separated by 0.14 inch; all panel letters use common left-column/right-column x anchors, shared row baselines, bold 13.2 pt letters and bold 9 pt titles.\n')

headings=[tx for tx in fig.texts if getattr(tx,'_panel_heading',False)]
letters={tx.get_text():tx for tx in headings if tx.get_text() in list('ABCDEFGH')}
assert len(letters)==8
for left,right in [('A','B'),('C','D'),('E','F'),('G','H')]:
 assert abs(letters[left].get_position()[1]-letters[right].get_position()[1])<1e-12
 assert abs(letters[left].get_position()[0]-.025)<1e-12 and abs(letters[right].get_position()[0]-.435)<1e-12
vp=OUT/'HET_display_clarity_validation.json';v=json.loads(vp.read_text());v['verified_panel_heading_positions']={k:list(tx.get_position()) for k,tx in letters.items()};vp.write_text(json.dumps(v,indent=2))
(OUT/'HET_DISPLAY_CLARITY_LEGEND.txt').write_text("""A: Broad provisional cell annotations, with separate Female, Male, WT, HET and KO insets on the same saved joint UMAP. The 72000 plotted QC cells comprise 6000 per sample across all 12 samples. RG is a shortened display label for the provisional RG-like population; no new annotation fit is performed.
B: Original saved MGE, LGE and CGE marker-program scores from the completed hierarchy-validation run, joined by exact cell ID to the same 72000 atlas display cells. No scores or UMAP coordinates are recomputed. All cells remain as gray background. Only saved regional scores strictly greater than 0.8 are overlaid in purple; scores <=0.8 stay gray. All three overlays share the same 0.8-to-pooled-98th-percentile colour scale, with upper clipping; the threshold is display-only and does not recalculate scores or filter the dataset. Scores provide descriptive regional-identity support rather than proving dissection boundaries. Full source and score summaries are in SAVED_REGIONAL_SCORES_PROVENANCE.json.
C: Female immature inhibitory neurons. Exon1 probe A UMI / total cell UMI * 10000. Pure WT female and HET WT cells use matched A-positive selections; HET KO cells are A-negative and downstream-positive. Their zero exon1 expression is therefore imposed by the selection. WT cells and KO cells are shortened display names for the previously defined WT-like and KO-like probe populations within female HET tissue, not confirmed cell-level DNA genotypes. Zero cells remain in their selected population. Violin mean points and quartile strips are unconnected.
D and F: Same 14 candidate genes, hierarchy, gene-row order, expression colour scale and effect limits as the preceding figure. Joint TMM logCPM display uses 15 sample/subset pseudobulk libraries; within-sex standardization precedes equal-sample group averaging. D retains the existing paired HET KO-cells versus WT-cells edgeR effects; F retains the unpaired KO-males versus WT-males effects. Blue/red indicate lower/higher gene expression for the explicitly labelled comparison. These candidate effects are exploratory and not gene-wide FDR significant.
E: Total Pcdh19 UMI (A+B+C) / total cell UMI * 10000. For this visualization ONLY, known WT males retain only immature inhibitory neurons with total Pcdh19 UMI >0. Known KO males retain ALL QC neurons in the same population, including zeros. The WT annotation reports the fraction retained from all eligible WT male neurons. This deliberately asymmetric display selection is not used for DE, network calculations, other panels or female populations. A+B+C can retain downstream signal in the known KO allele. Its display limit is based on the 99th percentile; all observations determine the violins, detection fractions and summaries.
G and H: Exactly the same five selected interaction candidates and effect values as the preceding figure, computed from original CPMs without depth matching. Potential=sqrt(source ligand expression * receiver receptor expression); complex subunits use geometric means. Each component is >=1 CPM in all three female HET samples in both states. Selection requires consistent direction and absolute log2 change >=0.2 in every sample; candidates rank by the smallest absolute sample effect, with unique interactions, up to two protocadherin interactions and one per other predefined neural family. Percentages are median KO-to-KO versus WT-to-WT changes across three paired HET samples. G arrows show source to receiver cell type; blue/red indicate lower/higher inferred strength in KO cells. Width scales to absolute median percent change. Within-type loops refer to populations, not demonstrated autocrine signaling. This is an expression-based exploratory schematic, not measured communication or an FDR-confirmed signaling model.
Layout: Former signaling heatmap D is removed. Panels read A/B, C/D, E/F, G/H. Four rows have 0.14 inch between their allocated blocks; panel letters share x anchors and row baselines, with bold 13.2 pt letters and bold 9 pt titles. Expression values, groups, definitions and fitted contrasts are unchanged by the layout edits.
""")
print('ALIGNED_PANEL_HEADINGS_VERIFIED',flush=True)

vp=OUT/'HET_display_clarity_validation.json';v=json.loads(vp.read_text());v.update(figure_size_inches=[7.2,8.49],regional_display_threshold=.8,regional_display_rule='strictly >0.8 coloured, all remaining cells grey',regional_highlighted_cells={name:int(highlighted[:,i].sum()) for i,name in enumerate(['MGE','LGE','CGE'])},regional_grey_cells={name:int((~highlighted[:,i]).sum()) for i,name in enumerate(['MGE','LGE','CGE'])},regional_shared_color_limits=score_limits.tolist(),A_layout='atlas and adjacent single-column biological legend above one aligned row of five equally sized sex/genotype insets',row_gap_inches=.14);vp.write_text(json.dumps(v,indent=2))
print('THRESHOLDED_REGIONAL_DISPLAY',v['regional_highlighted_cells'],flush=True)

vp=OUT/'HET_display_clarity_validation.json';v=json.loads(vp.read_text());v.update(E_WT_selection='plot only: total Pcdh19 A+B+C UMI >0',E_WT_original_cells=int((male_original_group==0).sum()),E_WT_displayed_cells=int(((male_original_group==0)&male_display_mask).sum()),E_KO_displayed_cells=int(((male_original_group==1)&male_display_mask).sum()),E_KO_original_cells=int((male_original_group==1).sum()),E_KO_selection='all eligible KO male neurons, including zeros',F_DE_selection_unchanged=True,regional_display_threshold=.8);assert v['E_KO_displayed_cells']==v['E_KO_original_cells'];vp.write_text(json.dumps(v,indent=2))
print('MALE_DISPLAY_SELECTION',v['E_WT_displayed_cells'],v['E_KO_displayed_cells'],flush=True)

vp=OUT/'HET_display_clarity_validation.json';v=json.loads(vp.read_text());v.update(figure_size_inches=[7.2,8.49],regional_shared_color_limits=[.8,1.5],regional_upper_visual_cap=1.5,regional_scale_orientation='vertical, right side; Low/High ticks; labelled cutoff 0.8 and cap 1.5',C_E_percentage_annotations=False,C_plot_rect=rects['b'],E_plot_rect=rects['e'],C_E_headers_and_neighboring_panels_fixed=True);vp.write_text(json.dumps(v,indent=2))
with open(str(OUT/'HET_DISPLAY_CLARITY_LEGEND.txt'),'a') as f:f.write('\nFINAL DISPLAY UPDATE: Regional UMAPs enlarged. Saved scores <=0.8 remain gray; scores >0.8 are colored with one shared 0.8–1.5 scale. Values >=1.5 receive the maximum color for visualization only; scores are not altered. The right-side vertical key has Low/High ticks, a Cutoff 0.8 marker and a Cap 1.5 marker. C/E percentage annotations are removed; only their plot axes are enlarged and lowered, while the headings and neighboring panels retain their positions. WT-male positive-only display selection remains unchanged.\n')
print('LARGE_UMAP_DISPLAY_COMPLETE',flush=True)

vp=OUT/'HET_display_clarity_validation.json';v=json.loads(vp.read_text());v.update(A_cell_type_title_removed=True,A_categorical_palette=['#216869','#B99A38','#C57750','#8F799E','#7A883E','#8B6C5C'],A_inset_colors_unchanged=True,B_colormap='truncated Greens',B_colorbar_width=.005,B_colorbar_height=.160,D_title='Female HET · neuronal genes',F_title='Male WT vs KO · neuronal genes',D_F_color_keys='Lower vs WT / Higher vs WT; D refers to WT cells in HET females, F refers to WT males',C_E_display_metric='natural log(1 + raw feature UMI / total cell UMI * 10000)',C_E_shared_display_limits=[0,3],C_E_population_selections_unchanged=True);vp.write_text(json.dumps(v,indent=2))
with open(str(OUT/'HET_DISPLAY_CLARITY_LEGEND.txt'),'a') as f:f.write('\nFINAL PALETTE AND EXPRESSION DISPLAY: A has a new categorical palette with no Cell type heading; the five sex/genotype inset colors are unchanged. G node colors follow the same cell-type palette as A. B uses green with an understated vertical scale, preserving the >0.8 display threshold and 1.5 visual cap. D is titled Female HET · neuronal genes; its bars compare KO cells versus WT cells within HET females. F is titled Male WT vs KO · neuronal genes; its bars compare KO males versus WT males. Their local blue/red keys state Lower vs WT / Higher vs WT, referring to gene expression relative to the specified WT reference. C/E are both plotted as ln(1 + feature UMI / total cell UMI * 10000), on the same 0–3 display range, with quartile strips and means calculated in those plotted units. C still measures exon1, E still measures total Pcdh19; all selections remain unchanged, including WT-male positive-only visualization. Values above 3 remain in density/statistic calculations but fall beyond the display axis; fractions are saved in validation. No DE or signaling calculation is refitted.\n')
print('PALETTE_AND_LOG_EXPRESSION_COMPLETE',flush=True)

vp=OUT/'HET_display_clarity_validation.json';v=json.loads(vp.read_text());v.update(A_title='Mouse MGE · E14.5 · n = 12',A_n_unit='dissected MGE samples',D_F_color_keys='D: lower/higher in HET KO cells vs HET WT cells; F: lower/higher in KO males vs WT males',D_F_keys_position='far right, stacked; duplicate bar headers removed',D_F_plot_rectangles={k:rects[k] for k in ['c','f']});vp.write_text(json.dumps(v,indent=2))

vp=OUT/'HET_display_clarity_validation.json';v=json.loads(vp.read_text());v.update(A_title='Mouse MGE · E14.5',A_subtitle='12 samples from 12 mice · dissected MGE',A_sample_breakdown={'WT female':3,'HET female':3,'WT male':3,'KO male':3},A_inset_dot_size=.045,B_cap_cutoff_text_removed=True,C_zero_label='0 exon-1 signal',C_HET_KO_cells=2180,C_zero_reason='A-negative/downstream-positive gate; no expression values changed');vp.write_text(json.dumps(v,indent=2))
