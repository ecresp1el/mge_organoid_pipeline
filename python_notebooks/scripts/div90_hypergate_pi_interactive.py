#!/usr/bin/env python3
"""Three self-contained, PI-facing views of the frozen DIV90 Phase 2 results.

No scoring, gate fitting, biological reanalysis, or source-file mutation occurs.
All figures use exact saved coordinates and full-precision RNA gate thresholds.
"""
from __future__ import annotations

import argparse
import html
import json
import os
from pathlib import Path
import sys

ROOT = Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_hypergate_sst_pv_phase2')
OUT = ROOT / 'pi_figure_package_v1'
STATES = ['PV-biased', 'PV/SST dual-high', 'SST-biased', 'unresolved']
RENAME = {'PV/SST hybrid': STATES[1], 'unresolved/immature': STATES[3]}
COLORS = ['#326DAB', '#39876F', '#CB733F', '#A6ACB6']
PAGES = [
    ('01_PV_SST_STATE_MAP.html', 'state',
     'DIV90 PV- and SST-associated developmental programs',
     'Hover over a cell to see its identity, sample, PV score, SST score and marker expression.',
     'Each dot is one of 4,768 recovered cortical LHX6+/ERBB4+ cells at day 90 in vitro (DIV90). '
     'The independent PV (parvalbumin)- and SST (somatostatin)-associated scores reveal developmental '
     'information before PVALB RNA is readily detectable. Dashed lines define operational score regions; '
     'dual-high means overlapping programs, not an established hybrid fate.'),
    ('02_VIRTUAL_FACS_GATE.html', 'gate',
     'Virtual FGFR2/PTPRS depletion gate',
     'This plot shows which DIV90 cells would be removed or retained by the proposed RNA-derived gate.',
     'The candidate removes cells with high FGFR2 RNA OR high PTPRS RNA, retaining the complement. '
     'Its main value is preserving PV (parvalbumin)-associated and dual-high cells. '
     'This is a prospective hypothesis for fluorescence-activated cell sorting (FACS), '
     'not a validated protein sorting protocol. SST means somatostatin; dual-high means overlapping PV/SST programs. RNA expression is not protein fluorescence.'),
    ('03_CULTURE_CONDITION_COMPARISON.html', 'condition',
     'Cell-line and culture-condition comparison',
     'Choose a cell line to compare CV and MW developmental-state composition.',
     'Compare culture conditions within the same cell line. CV is the higher-glucose condition (~2x); '
     'MW is the lower-glucose condition (~1x). PV means parvalbumin; SST means somatostatin. '
     'The four states describe continuous developmental programs, not proven fates. '
     'Condition is confounded with operator/culture differences; no glucose causality is claimed.')
]

TEMPLATE = r'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>__TITLE__</title><style>
*{box-sizing:border-box}body{margin:0;color:#243544;background:#F4F6F8;font:18px/1.5 system-ui,-apple-system,sans-serif}header{padding:28px 34px 20px;background:white;border-bottom:1px solid #DDE3E8}h1{font-size:34px;line-height:1.2;font-weight:680;margin:0 0 12px}h2{font-size:23px;line-height:1.25;margin:0 0 10px}h3{font-size:20px;margin:0 0 8px}p{margin:6px 0}main{max-width:1600px;margin:auto;padding:20px 28px}#instruction{font-size:21px}.explain{padding:18px 22px;border-left:5px solid #326DAB;background:#EAF0F6;border-radius:5px;margin-bottom:20px}.explain p{max-width:1250px}.controls{display:flex;gap:18px;align-items:end;flex-wrap:wrap;margin:0 0 18px;padding:18px;background:white;border:1px solid #DDE3E8;border-radius:8px}label{font-size:18px;font-weight:650;display:flex;flex-direction:column;gap:6px}select,button{font:18px system-ui;min-height:46px;border:1px solid #AAB8C4;border-radius:5px;background:white;padding:10px 12px;color:#243544}button{cursor:pointer;background:#EAF0F6;font-weight:650}select{min-width:200px}.cards{display:flex;flex-wrap:wrap;gap:12px;margin:0 0 18px}.card{flex:1;min-width:180px;background:white;border:1px solid #DDE3E8;border-radius:8px;padding:13px 17px;line-height:1.3}.card strong{font-size:28px;display:block;margin:7px 0}.card small{font-size:16px;color:#4F6170}.plotbox{display:grid;grid-template-columns:minmax(450px,1fr) 320px;gap:20px;background:white;padding:20px;border:1px solid #DDE3E8;border-radius:8px}canvas{width:100%;height:620px;display:block;touch-action:none}#conditionBars{height:420px}.landscapes{display:grid;grid-template-columns:1fr 1fr;gap:18px;margin:18px 0}.landscapes>div{background:white;padding:18px;border:1px solid #DDE3E8;border-radius:8px}.landscapes canvas{height:430px}.swatch{display:inline-block;width:15px;height:15px;border-radius:50%;margin-right:9px}#legend{font-size:18px;line-height:1.8;margin-bottom:18px}#hover{font-size:18px;line-height:1.5;overflow-wrap:anywhere;border-top:1px solid #DDE3E8;padding-top:16px}#hover b{font-weight:650}#summary{margin:0 0 18px}table{width:100%;border-collapse:collapse;background:white;font-size:18px}th,td{padding:10px 14px;text-align:left;border-bottom:1px solid #DDE3E8}th{background:#EAF0F6;font-weight:650}.tablewrap{overflow:auto}.note{color:#4F6170;margin:14px 0;font-size:17px}.emphasis{padding:12px 16px;background:#EDF4F0;border-radius:5px;margin-bottom:18px}.warning{background:#FFF2DE;padding:12px 16px;margin:12px 0;border-radius:5px}footer{padding:18px 0;color:#526472;font-size:16px}#floatingHover{position:fixed;right:24px;bottom:24px;max-width:440px;max-height:80vh;overflow:auto;padding:18px;background:white;border:1px solid #AAB8C4;border-radius:8px;box-shadow:0 5px 25px #24354433;z-index:5;font-size:18px;line-height:1.4;overflow-wrap:anywhere;pointer-events:none}.gradient{height:18px;background:linear-gradient(90deg,#EDF2F7,#326DAB);margin:8px 0}[hidden]{display:none!important}
@media(max-width:1000px){h1{font-size:29px}.plotbox{grid-template-columns:1fr}canvas{height:540px}.landscapes{grid-template-columns:1fr}main{padding:16px}header{padding:22px}aside{display:grid;grid-template-columns:1fr 1fr;gap:18px}#hover{border-top:0;padding-top:0}}
@media(max-width:620px){aside{display:block}.plotbox{padding:10px;display:block}canvas{height:450px}h1{font-size:26px}body{font-size:17px}th,td{font-size:16px;padding:8px}.controls{padding:12px}select{min-width:175px}.landscapes canvas{height:390px}}
</style></head><body>
<header><h1>__TITLE__</h1><p id="instruction">__INSTRUCTION__</p></header>
<main><section class="explain"><h2>What am I looking at?</h2><p>__EXPLANATION__</p></section>
<div id="controls" class="controls"></div><div id="message" class="emphasis"></div><div id="metrics" class="cards"></div><section id="summary"></section>
<div class="plotbox"><div><canvas id="plot" aria-label="Single-cell scatterplot"></canvas><canvas id="conditionBars" aria-label="Culture-condition developmental-state composition" hidden></canvas></div><aside><div id="legend"></div><div id="hover">Hover over a cell to see its identity, sample, scores and RNA expression.</div></aside></div>
<div class="landscapes" id="landscapes" hidden><div><h3 id="cvTitle">CV · higher glucose (~2x)</h3><canvas id="cvPlot" aria-label="CV developmental landscape"></canvas></div><div><h3 id="mwTitle">MW · lower glucose (~1x)</h3><canvas id="mwPlot" aria-label="MW developmental landscape"></canvas></div></div>
<p class="note" id="detail"></p><footer><p>Terms: DIV90 = day 90 in vitro; PV = parvalbumin; SST = somatostatin; PVALB = the parvalbumin gene; RNA = ribonucleic acid. Other capitalized gene symbols identify the measured genes.</p><p>Self-contained offline page: cell data and drawing code are embedded. Scroll over a scatterplot to zoom; use Reset view to restore the default. Dots use the saved coordinates without jitter. Dense overlaps can hide individual cells; all counts include every cell.</p></footer>
</main><div id="floatingHover" hidden aria-live="polite"></div><script type="application/json" id="payload">__PAYLOAD__</script><script>
'use strict';
const P=JSON.parse(document.getElementById('payload').textContent),rows=P.rows,idx=Object.fromEntries(P.columns.map((x,i)=>[x,i])),S=P.states,C=P.colors;
const A={shown:[],screens:{},bounds:{},stats:{},barHits:[]};window.PI_APP=A;
const E=(tag,props={})=>Object.assign(document.createElement(tag),props),byId=id=>document.getElementById(id),v=(i,k)=>rows[i][idx[k]],fmt=n=>n.toLocaleString('en-US');
// Round count ratios to one decimal with ties to even, matching the Python figures.
function pct(n,d){if(!d)return '—';const numerator=n*1000,base=Math.floor(numerator/d),remainder=numerator-base*d,rounded=base+(2*remainder>d||(2*remainder===d&&base%2===1)?1:0);return(rounded/10).toFixed(1)+'%';}
const esc=v=>String(v??'Unavailable').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const count=(ids,s)=>ids.filter(i=>v(i,'state')===s).length,all=rows.map((_,i)=>i),stateCounts=ids=>S.map(s=>count(ids,s));
function select(name,id,options,start){const label=E('label',{textContent:name}),node=E('select',{id});node.setAttribute('aria-label',name);for(const [value,label]of options)node.append(E('option',{value,textContent:label}));node.value=start;node.onchange=render;label.append(node);byId('controls').append(label);return node;}
let colorSelect=null,lineSelect,conditionSelect=null;
if(P.mode==='state')colorSelect=select('Color by','colorBy',[['state','Operational state'],['cell_line','Cell line'],['culture_operator','CV / MW culture condition'],['SST','SST RNA (somatostatin)'],['ERBB4','ERBB4 RNA'],['FGFR2','FGFR2 RNA'],['PTPRS','PTPRS RNA']],'state');
lineSelect=select('Cell line','cellLine',P.mode==='condition'?[['H9','H9'],['79B','79B'],['2E','2E']]:[['All','All cell lines'],['H9','H9'],['79B','79B'],['2E','2E']],P.mode==='condition'?'H9':'All');
if(P.mode!=='condition')conditionSelect=select('Culture condition','cultureCondition',[['All','CV and MW'],['CV','CV · higher glucose (~2x)'],['MW','MW · lower glucose (~1x)']],'All');
byId('controls').append(E('button',{id:'resetView',textContent:'Reset view',onclick:()=>{lineSelect.value=P.mode==='condition'?'H9':'All';if(conditionSelect)conditionSelect.value='All';if(colorSelect)colorSelect.value='state';A.bounds={};byId('floatingHover').hidden=true;byId('hover').textContent='Hover over a cell to see its identity, sample, scores and RNA expression.';render();}}));
const removed=all.map(i=>P.rule.rules.some(r=>v(i,r.gene)>r.threshold));
function card(name,number,detail,color){const div=E('div',{className:'card'});if(color)div.style.borderTop='5px solid '+color;div.append(E('div',{textContent:name}),E('strong',{textContent:number}),E('small',{textContent:detail}));byId('metrics').append(div);}
function table(headers,records){const wrap=E('div',{className:'tablewrap'}),t=E('table'),head=E('thead'),tr=E('tr'),body=E('tbody');headers.forEach(h=>tr.append(E('th',{textContent:h})));head.append(tr);for(const record of records){const row=E('tr');record.forEach(x=>row.append(E('td',{textContent:x})));body.append(row);}t.append(head,body);wrap.append(t);return wrap;}
function render(){A.shown=all.filter(i=>(lineSelect.value==='All'||v(i,'cell_line')===lineSelect.value)&&(!conditionSelect||conditionSelect.value==='All'||v(i,'culture_operator')===conditionSelect.value));A.stats={n:A.shown.length,stateCounts:stateCounts(A.shown)};byId('metrics').replaceChildren();byId('summary').replaceChildren();
if(P.mode==='state'){card('Entry cells displayed',fmt(A.shown.length),'4,768 recovered cortical LHX6+/ERBB4+ cells');S.forEach((s,j)=>card(s,pct(A.stats.stateCounts[j],A.shown.length),fmt(A.stats.stateCounts[j])+' cells',C[j]));byId('message').textContent='27.5% of all entry cells are dual-high / overlapping programs; coordinated hybrid identity is not yet established.';byId('detail').textContent='Global operational thresholds remain fixed when filtering. RNA colors show normalized, log-transformed expression (log1p of counts per 10,000); scores are independent standardized gene-program averages.';drawScatter('plot',A.shown,'sst_score','pv_score',colorSelect.value);}
if(P.mode==='gate'){const keep=A.shown.filter(i=>!removed[i]),drop=A.shown.filter(i=>removed[i]);A.stats.retained=keep.length;A.stats.removed=drop.length;A.stats.retainedCounts=stateCounts(keep);A.stats.removedCounts=stateCounts(drop);card('Entry cells displayed',fmt(A.shown.length),'Saved RNA-derived rule');card('Retained',fmt(keep.length),pct(keep.length,A.shown.length)+' total cell recovery','#39876F');card('Removed',fmt(drop.length),'Collect this fraction for validation','#CB733F');card('PV-biased recovery',pct(count(keep,S[0]),count(A.shown,S[0])),'Fraction of starting PV-biased cells',C[0]);card('Dual-high recovery',pct(count(keep,S[1]),count(A.shown,S[1])),'Fraction of starting dual-high cells',C[1]);byId('message').textContent='Remove FGFR2-high OR PTPRS-high; retain cells below or equal to BOTH saved RNA thresholds. Preservation is the main benefit; SST-biased composition changes only modestly.';byId('summary').append(table(['Composition within each fraction','Retained ('+fmt(keep.length)+' cells)','Removed ('+fmt(drop.length)+' cells)'],S.map(s=>[s,pct(count(keep,s),keep.length)+' · '+fmt(count(keep,s))+' cells',pct(count(drop,s),drop.length)+' · '+fmt(count(drop,s))+' cells'])));byId('detail').textContent='Exact saved RNA rule: remove FGFR2 > '+P.rule.rules[0].threshold+' OR PTPRS > '+P.rule.rules[1].threshold+'. Values are log1p of counts per 10,000 RNA molecules; they are not fluorescence cutoffs. Keep and mature both fractions plus an unsorted ERBB4+ control to test later PV/SST protein, transcriptomic identity and physiology. No protein separation is established by this view.';drawScatter('plot',A.shown,'FGFR2','PTPRS','gate');}
if(P.mode==='condition'){byId('plot').hidden=true;byId('conditionBars').hidden=false;byId('landscapes').hidden=false;const groups=['CV','MW'].map(k=>A.shown.filter(i=>v(i,'culture_operator')===k));A.stats.groups=Object.fromEntries(groups.map((g,j)=>[['CV','MW'][j],{n:g.length,counts:stateCounts(g)}]));groups.forEach((g,j)=>card(['CV · higher glucose (~2x)','MW · lower glucose (~1x)'][j],fmt(g.length)+' cells',lineSelect.value+' entry population'));byId('message').textContent='Compare CV versus MW within '+lineSelect.value+'. Both landscapes use exactly the same axes and operational thresholds.';const small=E('p',{className:lineSelect.value==='2E'?'warning':'note',textContent:'2E/CV has only 37 entry cells; its percentages are especially sensitive to small cell counts. Each line-condition combination is represented by one sample.'});byId('summary').append(small,table(['Developmental state','CV · higher glucose (~2x)','MW · lower glucose (~1x)'],S.map(s=>[s,...groups.map(g=>pct(count(g,s),g.length)+' · '+fmt(count(g,s))+' cells')])));byId('detail').textContent='CV versus MW culture condition includes an approximately 2-fold glucose difference. Condition is confounded with operator/culture differences; no glucose causality is claimed. Cells are observations within samples, not independent culture replicates.';drawBars(groups);drawScatter('cvPlot',groups[0],'sst_score','pv_score','state',true);drawScatter('mwPlot',groups[1],'sst_score','pv_score','state',true);byId('cvTitle').textContent='CV · higher glucose (~2x) · '+fmt(groups[0].length)+' cells';byId('mwTitle').textContent='MW · lower glucose (~1x) · '+fmt(groups[1].length)+' cells';legend('state',A.shown);}}
function scaleBounds(xkey,ykey){const xs=all.map(i=>v(i,xkey)),ys=all.map(i=>v(i,ykey)),a=Math.min(...xs),b=Math.max(...xs),c=Math.min(...ys),d=Math.max(...ys);return[a-(b-a)*.05,b+(b-a)*.05,c-(d-c)*.06,d+(d-c)*.06];}
function prepare(id){const canvas=byId(id),rect=canvas.getBoundingClientRect(),dpr=window.devicePixelRatio||1;canvas.width=Math.round(rect.width*dpr);canvas.height=Math.round(rect.height*dpr);const ctx=canvas.getContext('2d');ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,rect.width,rect.height);return{canvas,ctx,w:rect.width,h:rect.height};}
function colorFn(key){if(key==='state')return i=>C[S.indexOf(v(i,'state'))];if(key==='gate')return i=>removed[i]?'#CB733F':'#39876F';if(key==='cell_line')return i=>({'H9':'#326DAB','79B':'#A46C9F','2E':'#A98631'}[v(i,key)]);if(key==='culture_operator')return i=>v(i,key)==='CV'?'#8A5E97':'#31818E';const max=Math.max(...all.map(i=>v(i,key)));return i=>{const q=max?v(i,key)/max:0;return`rgb(${Math.round(237-187*q)},${Math.round(242-133*q)},${Math.round(247-76*q)})`;};}
function legend(key,ids){const root=byId('legend');root.replaceChildren(E('h3',{textContent:key==='gate'?'Proposed RNA gate':key==='state'?'Operational state':key==='cell_line'?'Cell line':key==='culture_operator'?'Culture condition':key+' RNA'}));const col=colorFn(key);if(['state','gate','cell_line','culture_operator'].includes(key)){const categories=key==='state'?S:key==='gate'?['Retained','Removed']:key==='cell_line'?['H9','79B','2E']:['CV','MW'];for(const category of categories){const div=E('div'),swatch=E('span',{className:'swatch'});let color=key==='state'?C[S.indexOf(category)]:key==='gate'?(category==='Retained'?'#39876F':'#CB733F'):key==='cell_line'?({'H9':'#326DAB','79B':'#A46C9F','2E':'#A98631'}[category]):category==='CV'?'#8A5E97':'#31818E';swatch.style.background=color;div.append(swatch,document.createTextNode(key==='culture_operator'?category+(category==='CV'?' · higher glucose (~2x)':' · lower glucose (~1x)'):category));root.append(div);}}else{const max=Math.max(...all.map(i=>v(i,key)));root.append(E('div',{className:'gradient'}),E('div',{textContent:'Undetected → higher RNA'}),E('p',{className:'note',textContent:'0 to '+max.toFixed(2)+' normalized, log-transformed expression. The color scale stays fixed across filters.'}));}}
function drawScatter(id,ids,xkey,ykey,color,paired=false){const{canvas,ctx,w,h}=prepare(id),m={left:74,right:19,top:36,bottom:70},pw=w-m.left-m.right,ph=h-m.top-m.bottom,bkey=paired?'paired':id;if(!A.bounds[bkey])A.bounds[bkey]=scaleBounds(xkey,ykey);const[a,b,c,d]=A.bounds[bkey],sx=x=>m.left+(x-a)/(b-a)*pw,sy=y=>m.top+ph-(y-c)/(d-c)*ph;ctx.font='15px system-ui';ctx.fillStyle='#526472';ctx.strokeStyle='#E5EAF0';ctx.lineWidth=1;for(let j=0;j<=4;j++){const x=a+(b-a)*j/4,y=c+(d-c)*j/4;ctx.beginPath();ctx.moveTo(sx(x),m.top);ctx.lineTo(sx(x),m.top+ph);ctx.moveTo(m.left,sy(y));ctx.lineTo(m.left+pw,sy(y));ctx.stroke();ctx.textAlign='center';ctx.fillText(x.toFixed(1),sx(x),h-43);ctx.textAlign='right';ctx.fillText(y.toFixed(1),m.left-10,sy(y)+5);}
const gate=xkey==='FGFR2',cutx=gate?P.rule.rules[0].threshold:P.cuts[0],cuty=gate?P.rule.rules[1].threshold:P.cuts[1];if(gate){ctx.fillStyle='#F8E8DE';ctx.fillRect(Math.max(m.left,sx(cutx)),m.top,Math.max(0,m.left+pw-Math.max(m.left,sx(cutx))),ph);ctx.fillRect(m.left,m.top,Math.min(pw,Math.max(0,sx(cutx)-m.left)),Math.max(0,Math.min(ph,sy(cuty)-m.top)));ctx.fillStyle='#EDF4F0';ctx.fillRect(m.left,Math.max(m.top,sy(cuty)),Math.min(pw,Math.max(0,sx(cutx)-m.left)),Math.max(0,m.top+ph-Math.max(m.top,sy(cuty))));}
ctx.save();ctx.beginPath();ctx.rect(m.left,m.top,pw,ph);ctx.clip();ctx.setLineDash([6,5]);ctx.strokeStyle='#536573';ctx.beginPath();ctx.moveTo(sx(cutx),m.top);ctx.lineTo(sx(cutx),m.top+ph);ctx.moveTo(m.left,sy(cuty));ctx.lineTo(m.left+pw,sy(cuty));ctx.stroke();ctx.setLineDash([]);const col=colorFn(color),numeric=!['state','gate','cell_line','culture_operator'].includes(color),ordered=numeric?[...ids].sort((i,j)=>v(i,color)-v(j,color)):ids;A.screens[id]=[];for(const i of ordered){const x=sx(v(i,xkey)),y=sy(v(i,ykey));ctx.globalAlpha=.72;ctx.fillStyle=col(i);ctx.beginPath();ctx.arc(x,y,paired?2.3:2.7,0,2*Math.PI);ctx.fill();if(x>=m.left&&x<=m.left+pw&&y>=m.top&&y<=m.top+ph)A.screens[id].push([x,y,i]);}ctx.restore();ctx.globalAlpha=1;ctx.fillStyle='#243544';ctx.textAlign='center';ctx.font='18px system-ui';ctx.fillText(gate?'FGFR2 RNA expression':'SST-associated program score',m.left+pw/2,h-12);ctx.save();ctx.translate(20,m.top+ph/2);ctx.rotate(-Math.PI/2);ctx.fillText(gate?'PTPRS RNA expression':'PV-associated program score',0,0);ctx.restore();ctx.font='16px system-ui';ctx.textAlign='left';ctx.fillText(gate?'Dashed lines: saved RNA thresholds':'PV-biased',m.left,m.top-13);ctx.textAlign='right';ctx.fillText(gate?'Remove if either marker is high':'Dual-high / overlapping programs',m.left+pw,m.top-13);if(gate){ctx.font='bold 17px system-ui';ctx.textAlign='left';ctx.fillStyle='#286850';ctx.fillText('RETAIN · both low',m.left+12,m.top+ph-12);ctx.textAlign='right';ctx.fillStyle='#A5592E';ctx.fillText('REMOVE · either high',m.left+pw-12,m.top+25);}else{ctx.font='16px system-ui';ctx.textAlign='left';ctx.fillStyle='#66717F';ctx.fillText('Unresolved',m.left+8,m.top+ph-9);ctx.textAlign='right';ctx.fillStyle='#AA5A2A';ctx.fillText('SST-biased',m.left+pw-8,m.top+ph-9);}if(!paired)legend(color,ids);canvas._params={xkey,ykey,color,paired};}
function drawBars(groups){const{ctx,w,h}=prepare('conditionBars'),m={left:70,right:25,top:20,bottom:72},ph=h-m.top-m.bottom,pw=w-m.left-m.right;ctx.font='17px system-ui';ctx.strokeStyle='#E5EAF0';ctx.fillStyle='#526472';for(let j=0;j<=4;j++){const y=m.top+ph-j/4*ph;ctx.beginPath();ctx.moveTo(m.left,y);ctx.lineTo(m.left+pw,y);ctx.stroke();ctx.textAlign='right';ctx.fillText((j*25)+'%',m.left-10,y+5);}A.barHits=[];groups.forEach((ids,g)=>{const x=m.left+pw*(g+.5)/2,bw=Math.min(165,pw*.30);let bottom=m.top+ph;S.forEach((s,j)=>{const n=count(ids,s),bh=ids.length?n/ids.length*ph:0;ctx.fillStyle=C[j];ctx.fillRect(x-bw/2,bottom-bh,bw,bh);ctx.fillStyle=j===3?'#243544':'white';ctx.font='bold 19px system-ui';ctx.textAlign='center';if(bh>25)ctx.fillText(pct(n,ids.length),x,bottom-bh/2+7);A.barHits.push({x:x-bw/2,y:bottom-bh,w:bw,h:bh,state:s,n,total:ids.length,condition:['CV','MW'][g]});bottom-=bh;});ctx.fillStyle='#243544';ctx.font='19px system-ui';ctx.fillText(['CV · ~2x glucose','MW · ~1x glucose'][g],x,h-39);ctx.font='17px system-ui';ctx.fillText(fmt(ids.length)+' cells',x,h-14);});}
function hoverCell(i,floating=false){const fields=[['Cell',v(i,'cell_id')],['Sample',v(i,'sample')],['Cell line',v(i,'cell_line')],['Culture condition',v(i,'condition')],['Operational state',v(i,'state')],['PV-associated score',v(i,'pv_score').toFixed(4)],['SST-associated score',v(i,'sst_score').toFixed(4)],...['SST','ERBB4','FGFR2','PTPRS','FAT3','PTPRM'].map(g=>[g+' RNA',v(i,g).toFixed(4)])];if(P.mode==='gate')fields.push(['Proposed fraction',removed[i]?'Removed':'Retained']);byId('hover').innerHTML='<h3>Cell details</h3>'+fields.map(([k,val])=>'<b>'+esc(k)+':</b> '+esc(val)).join('<br>');if(floating){byId('floatingHover').innerHTML=byId('hover').innerHTML;byId('floatingHover').hidden=false;}}
for(const id of ['plot','cvPlot','mwPlot']){const canvas=byId(id);canvas.addEventListener('mousemove',e=>{const r=canvas.getBoundingClientRect(),x=e.clientX-r.left,y=e.clientY-r.top;let best=null,dmin=144;for(const p of A.screens[id]||[]){const d=(p[0]-x)**2+(p[1]-y)**2;if(d<dmin){best=p;dmin=d;}}if(best)hoverCell(best[2],id!=='plot');});canvas.addEventListener('mouseleave',()=>{byId('floatingHover').hidden=true;});canvas.addEventListener('wheel',e=>{e.preventDefault();const p=canvas._params;if(!p)return;const key=p.paired?'paired':id,[a,b,c,d]=A.bounds[key],r=canvas.getBoundingClientRect(),fx=Math.max(0,Math.min(1,(e.clientX-r.left-74)/(r.width-93))),fy=1-Math.max(0,Math.min(1,(e.clientY-r.top-36)/(r.height-106))),cx=a+(b-a)*fx,cy=c+(d-c)*fy,k=e.deltaY>0?1.15:1/1.15;A.bounds[key]=[cx+(a-cx)*k,cx+(b-cx)*k,cy+(c-cy)*k,cy+(d-cy)*k];render();},{passive:false});}
byId('conditionBars').addEventListener('mousemove',e=>{const r=e.currentTarget.getBoundingClientRect(),x=e.clientX-r.left,y=e.clientY-r.top,hit=A.barHits.find(b=>x>=b.x&&x<=b.x+b.w&&y>=b.y&&y<=b.y+b.h);if(hit)byId('hover').innerHTML='<h3>'+esc(hit.condition+' · '+lineSelect.value)+'</h3><b>'+esc(hit.state)+'</b><br>'+fmt(hit.n)+' of '+fmt(hit.total)+' cells ('+pct(hit.n,hit.total)+')';});
window.addEventListener('resize',render);render();A.ready=true;
</script></body></html>'''


def run(cells, out, summary, params):
    """Write exactly three standalone interactive pages from saved cells/results."""
    import numpy as np
    import pandas as pd
    out = Path(out)
    if out.resolve() != OUT.resolve():
        raise ValueError(f'Interactive outputs must stay under {OUT}')
    df = pd.read_csv(cells, sep='\t', float_precision='round_trip') if isinstance(cells, (str, Path)) else cells.copy()
    summary = json.loads(Path(summary).read_text()) if isinstance(summary, (str, Path)) else summary
    params = json.loads(Path(params).read_text()) if isinstance(params, (str, Path)) else params
    assert len(df) == 4768
    df['state'] = df['state'].replace(RENAME)
    assert df['state'].value_counts().reindex(STATES).tolist() == [1076, 1309, 1075, 1308]
    assert len(df[(df.cell_line == '2E') & (df.culture_operator == 'CV')]) == 37
    rule = summary['rules'][summary['experimental_depletion']['gate_id']]
    assert rule['action'] == 'remove' and rule['logic'] == 'OR'
    assert [(r['gene'], r['op']) for r in rule['rules']] == [('FGFR2', '>'), ('PTPRS', '>')]
    removed = np.logical_or.reduce([df[r['gene']].to_numpy() > r['threshold'] for r in rule['rules']])
    assert int(removed.sum()) == 744 and int((~removed).sum()) == 4024
    columns = ['cell_id', 'sample', 'cell_line', 'culture_operator', 'condition', 'state',
               'pv_score', 'sst_score', 'SST', 'ERBB4', 'FGFR2', 'PTPRS', 'FAT3', 'PTPRM']
    values = df[columns].astype(object).where(pd.notnull(df[columns]), None).values.tolist()
    payload = {'columns': columns, 'rows': values, 'states': STATES, 'colors': COLORS, 'rule': rule,
               'cuts': [params['sst_threshold'], params['pv_threshold']]}
    directory = out / 'interactive'
    directory.mkdir(parents=True, exist_ok=True)
    manifest = []
    for filename, mode, title, instruction, explanation in PAGES:
        page = TEMPLATE
        replacements = {'__TITLE__': html.escape(title), '__INSTRUCTION__': html.escape(instruction),
                        '__EXPLANATION__': html.escape(explanation),
                        '__PAYLOAD__': json.dumps({**payload, 'mode': mode}, ensure_ascii=False, allow_nan=False, separators=(',', ':')).replace('</', '<\\/')}
        for key, value in replacements.items():
            page = page.replace(key, value)
        (directory / filename).write_text(page, encoding='utf-8')
        manifest.append({'file': f'interactive/{filename}', 'title': title, 'offline': True,
                         'cells_embedded': 4768, 'coordinates_jittered': False})
    return manifest


def verify_browser(out=OUT, cells=None):
    """Verify all three pages with Chromium offline and compare to saved data."""
    import pandas as pd
    out = Path(out)
    if out.resolve() != OUT.resolve():
        raise ValueError('Browser artifacts must remain under the PI package on Turbo')
    df = (pd.read_csv(ROOT / 'cells.tsv.gz', sep='\t', float_precision='round_trip')
          if cells is None else cells.copy())
    df['state'] = df['state'].replace(RENAME)
    qa = out / 'provenance'
    qa.mkdir(parents=True, exist_ok=True)
    tmp = out / 'cache/browser_tmp'
    tmp.mkdir(parents=True, exist_ok=True)
    os.environ['TMPDIR'] = str(tmp)
    import tempfile
    tempfile.tempdir = str(tmp)
    sys.path.insert(0, str(ROOT / 'runtime_deps'))
    os.environ['PLAYWRIGHT_BROWSERS_PATH'] = str(ROOT / 'cache/playwright')
    os.environ['NODE_OPTIONS'] = '--require=' + str(ROOT / 'provenance/playwright_autofs_shim.cjs')
    from playwright.sync_api import sync_playwright
    report = {'offline': True, 'all_network_blocked': True,
              'percentage_rounding': 'One decimal; exact count-ratio ties to even', 'assets': []}
    with sync_playwright() as pw:
        browser = pw.chromium.launch(headless=True, args=['--no-sandbox'])
        context = browser.new_context(offline=True, viewport={'width': 1500, 'height': 1200})
        for filename, mode, title, instruction, _ in PAGES:
            page = context.new_page()
            errors, external_requests = [], []
            page.on('pageerror', lambda error: errors.append(str(error)))
            page.on('request', lambda request: external_requests.append(request.url) if request.url.startswith(('http:', 'https:')) else None)
            page.goto((out / 'interactive' / filename).as_uri(), wait_until='load')
            page.wait_for_function('window.PI_APP && PI_APP.ready')
            assert page.locator('h1').inner_text() == title
            assert page.locator('#instruction').inner_text() == instruction
            assert page.get_by_role('heading', name='What am I looking at?').is_visible()
            assert page.evaluate('pct(14, 224)') == '6.2%'
            assert page.evaluate('pct(3, 16)') == '18.8%'
            assert page.evaluate('pct(0, 0)') == '—'
            stats = page.evaluate('PI_APP.stats')
            if mode != 'condition':
                assert stats['n'] == 4768 and stats['stateCounts'] == [1076, 1309, 1075, 1308]
            if mode == 'gate':
                assert (stats['removed'], stats['retained']) == (744, 4024)
                assert stats['retainedCounts'] == [954, 1169, 856, 1045]
                assert stats['removedCounts'] == [122, 140, 219, 263]
            if mode == 'state':
                for color in ['state', 'cell_line', 'culture_operator', 'SST', 'ERBB4', 'FGFR2', 'PTPRS']:
                    page.get_by_label('Color by', exact=True).select_option(color)
                    assert page.evaluate('PI_APP.screens.plot.length') == 4768
            condition_checks = {}
            if mode == 'condition':
                assert page.get_by_label('Cell line', exact=True).input_value() == 'H9'
                for line in ['H9', '79B', '2E']:
                    page.get_by_label('Cell line', exact=True).select_option(line)
                    actual = page.evaluate('PI_APP.stats.groups')
                    for condition in ['CV', 'MW']:
                        sub = df[(df.cell_line == line) & (df.culture_operator == condition)]
                        expected = {'n': len(sub), 'counts': sub.state.value_counts().reindex(STATES, fill_value=0).tolist()}
                        assert actual[condition] == expected, (line, condition, actual, expected)
                        percentages = page.evaluate('([counts, n]) => counts.map(c => pct(c, n))', [expected['counts'], expected['n']])
                        assert percentages == [f'{100 * n / expected["n"]:.1f}%' for n in expected['counts']]
                    condition_checks[line] = actual
                assert condition_checks['2E']['CV']['n'] == 37
            else:
                page.get_by_label('Cell line', exact=True).select_option('H9')
                page.get_by_label('Culture condition', exact=True).select_option('CV')
                sub = df[(df.cell_line == 'H9') & (df.culture_operator == 'CV')]
                assert page.evaluate('PI_APP.stats.n') == len(sub)
                if mode == 'gate':
                    expected = int(((sub.FGFR2 <= 1.2491153968) & (sub.PTPRS <= 3.40762345848)).sum())
                    assert page.evaluate('PI_APP.stats.retained') == expected
            canvas_id = 'cvPlot' if mode == 'condition' else 'plot'
            canvas = page.locator('#' + canvas_id)
            canvas.scroll_into_view_if_needed()
            point = page.evaluate(f'PI_APP.screens.{canvas_id}[0]')
            box = canvas.bounding_box()
            page.mouse.move(box['x'] + point[0], box['y'] + point[1])
            assert 'Cell details' in page.locator('#hover').inner_text()
            if mode == 'condition':
                assert page.locator('#floatingHover').is_visible()
            key = 'paired' if mode == 'condition' else 'plot'
            old_bounds = page.evaluate(f'PI_APP.bounds.{key}')
            page.mouse.wheel(0, -200)
            page.wait_for_timeout(150)
            assert page.evaluate(f'PI_APP.bounds.{key}') != old_bounds
            page.get_by_role('button', name='Reset view', exact=True).click()
            assert page.get_by_label('Cell line', exact=True).input_value() == ('H9' if mode == 'condition' else 'All')
            assert page.evaluate(f'PI_APP.bounds.{key}') == old_bounds
            if mode != 'condition':
                assert page.evaluate('PI_APP.stats.n') == 4768
            page.screenshot(path=str(qa / filename.replace('.html', '_preview.png')), full_page=True)
            assert not errors, errors
            assert not external_requests, external_requests
            report['assets'].append({'filename': filename, 'default_counts': stats, 'hover': True,
                                     'filtering': True, 'reset_and_zoom': True, 'javascript_errors': errors,
                                     'external_requests': external_requests, 'condition_counts': condition_checks})
            page.close()
        browser.close()
    (qa / 'interactive_validation.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    run(ROOT / 'cells.tsv.gz', OUT, ROOT / 'gate_summary.json', ROOT / 'state_parameters.json')
    if args.verify:
        print(json.dumps(verify_browser(), indent=2))
    else:
        print('Wrote exactly three standalone PI-facing interactive pages.')
