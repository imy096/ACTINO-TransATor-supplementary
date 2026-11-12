#!/usr/bin/env python3
"""Regenerate ACTINO Option-B figures using only archived audit tables/GenBank files.
No network, no synthetic points, no statistical inference. Python >=3.10.
"""
from pathlib import Path
from collections import Counter
import csv, re
import numpy as np
import pandas as pd
from Bio import SeqIO
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle, Patch

ROOT=Path(__file__).resolve().parent
S=ROOT/'source_data'; F=ROOT/'figures'; F.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':11,'axes.labelsize':9,'savefig.facecolor':'white','pdf.fonttype':42,'svg.fonttype':'none'})
COL={'PKS':'#234F7D','NRPS':'#398F85','RiPP':'#B97835','terpene':'#8569A7','other':'#7A8996','saccharide':'#C55E56'}

def save(fig,name):
    for ext in ['svg','png','pdf']:
        fig.savefig(F/f'{name}.{ext}',bbox_inches='tight',dpi=360)
    plt.close(fig)

# Figure 1. Audited per-genome protocluster category tallies.
g=pd.read_csv(S/'antismash_genome_category_totals.csv').sort_values('total_protoclusters',ascending=True)
fig,ax=plt.subplots(figsize=(9,8)); left=np.zeros(len(g));
for cat in ['PKS','NRPS','RiPP','terpene','other','saccharide']:
    vals=g[cat].to_numpy(); ax.barh(g['strain'],vals,left=left,color=COL[cat],height=.71,label=cat);left+=vals
for y,v in enumerate(g.total_protoclusters):ax.text(v+.6,y,str(v),va='center',fontsize=7)
ax.set_xlim(0,max(g.total_protoclusters)+8);ax.set_xlabel('Predicted antiSMASH protoclusters (count)')
ax.set_title('Predicted biosynthetic protocluster composition across 20 genome assemblies',loc='left',fontweight='bold',pad=14)
ax.legend(loc='upper center',bbox_to_anchor=(.5,-.09),ncol=6,frameon=False,fontsize=8)
ax.grid(axis='x',alpha=.16);ax.set_axisbelow(True);ax.spines[['top','right']].set_visible(False)
save(fig,'Figure_1_protocluster_distribution')

# Figure 2. antiSMASH hierarchical units and edge representation.
reg=pd.read_csv(S/'antismash_regions_674.csv'); pc=pd.read_csv(S/'antismash_all_protoclusters_878.csv'); kind=pd.read_csv(S/'candidate_cluster_kind_counts.csv')
fig,axs=plt.subplots(1,2,figsize=(10.3,4.1),gridspec_kw={'width_ratios':[1.2,1]})
units=['Region areas','Candidate-cluster\nobjects','Protoclusters']; counts=[len(reg),int(kind['count'].sum()),len(pc)]
for i,(nm,v) in enumerate(zip(units,counts)):
    axs[0].barh(i,v,color=['#315E7B','#73859C','#268A7C'][i],height=.55)
    axs[0].text(v+15,i,f'{v:,}',va='center',fontweight='bold')
axs[0].set_yticks(range(3),units);axs[0].invert_yaxis();axs[0].set_xlim(0,1020)
axs[0].set_xlabel('Reported feature objects (not additive)');axs[0].set_title('Annotation hierarchy',loc='left',fontweight='bold');axs[0].spines[['top','right']].set_visible(False)
for n,(label,d) in enumerate([('Regions',reg),('Protoclusters',pc)]):
    col='region_touches_contig_edge' if n==0 else 'protocluster_touches_contig_edge'
    val=int((d[col].astype(str)=='True').sum()); total=len(d)
    axs[1].barh(n,val/total*100,color='#BB794B',height=.45)
    axs[1].barh(n,100-val/total*100,left=val/total*100,color='#C7D5DE',height=.45)
    axs[1].text(3,n,f'Edge: {val}/{total} ({val/total*100:.1f}%)',va='center',fontsize=8,color='#222222')
axs[1].set_yticks([0,1],['Regions','Protoclusters']);axs[1].invert_yaxis();axs[1].set_xlim(0,100)
axs[1].set_xlabel('% of annotated units');axs[1].set_title('Contig-edge annotations',loc='left',fontweight='bold');axs[1].spines[['top','right']].set_visible(False)
fig.suptitle('Hierarchical detection counts and assembly-boundary context',fontweight='bold',y=1.02)
fig.tight_layout();save(fig,'Figure_2_annotation_hierarchy_edge')

# Figure 3. Genomic synteny from GenBank annotated protein CDS, target orientation.
entries=[('CCDA1',S/'CCDA1_antismash_region002.gbk','AC67HL_09695'),('Microbispora CL1-1',S/'Microbispora_CL1-1_61kb.gb','HCK04_RS19955'),('Streptomyces fimbriatus',S/'Streptomyces_fimbriatus_60kb.gb','ABD362_RS17570')]
labels=['DUF6039','MalonylT','DUF6039','VOC','DUF3224','PIG-L','SAM','Hydrolase','ACP','AMP','TARGET','HMG','Carrier','KS-like']
# category colors preserve distinct target and major biosynthetic groups.
colors={'TARGET':'#0C756A','MalonylT':'#B76A3B','ACP':'#548A9C','AMP':'#8362A1','HMG':'#AA8A39','Carrier':'#548A9C','KS-like':'#C46D5B','other':'#84939F'}
def fam(i):return colors.get(labels[i],'#84939F')
fig,ax=plt.subplots(figsize=(13,4.1))
for j,(name,gb,tag) in enumerate(entries):
    rec=SeqIO.read(gb,'genbank'); cds=sorted((f for f in rec.features if f.type=='CDS'),key=lambda f:int(f.location.start))
    idx=next((i for i,f in enumerate(cds) if f.qualifiers.get('locus_tag',[''])[0]==tag),None)
    assert idx is not None,tag
    tg=cds[idx]; strand=tg.location.strand
    ordered=cds if strand==1 else list(reversed(cds));tidx=ordered.index(tg)
    for rel in range(-10,4):
        f=ordered[tidx+rel]
        # relative position on the local transcription-oriented genome axis
        if strand==1:
            xa=(int(f.location.start)-int(tg.location.start))/1000;xb=(int(f.location.end)-int(tg.location.start))/1000
        else:
            xa=(int(tg.location.end)-int(f.location.end))/1000;xb=(int(tg.location.end)-int(f.location.start))/1000
        k=rel+10
        # arrow follows actual transcription direction, after coordinate reversal
        direct=int(f.location.strand)*strand
        y=2-j*.95;hh=.35
        tip=min(0.32,max(.075,(xb-xa)*.25))
        if xb-xa<.25: tip=(xb-xa)*.32
        if direct>0: poly=[(xa,y),(xb-tip,y),(xb,y+hh/2),(xb-tip,y+hh),(xa,y+hh)]
        else:poly=[(xb,y),(xa+tip,y),(xa,y+hh/2),(xa+tip,y+hh),(xb,y+hh)]
        ax.add_patch(Polygon(poly,facecolor=fam(k),edgecolor='white',linewidth=.5))
    ax.text(-17.7,2-j*.95+.17,name,ha='right',va='center',fontsize=9,fontweight='bold')
ax.axvline(0,color='#30485C',linewidth=.75,linestyle=':')
ax.set_xlim(-18,29);ax.set_ylim(-.5,3);ax.set_xlabel('Genomic position relative to target CDS start, in target-oriented coordinates (kb)');ax.set_yticks([])
ax.set_title('Conserved local gene organization across genera',fontweight='bold',loc='left',pad=14)
ax.spines[['left','right','top']].set_visible(False);ax.grid(axis='x',alpha=.14);ax.set_axisbelow(True)
ax.legend(handles=[Patch(facecolor=c,label=l) for l,c in [('Large target',colors['TARGET']),('ACP/carrier',colors['ACP']),('AMP-binding',colors['AMP']),('HMG-CoA synthase',colors['HMG']),('KS-like',colors['KS-like']),('Other neighbors',colors['other'])]],loc='lower center',bbox_to_anchor=(.5,-.40),ncol=3,frameon=False,fontsize=8)
fig.subplots_adjust(bottom=.24);save(fig,'Figure_3_cross_genus_gene_synteny')

# Figure 4. Selected domains in CCDA1, Microbispora and GnmT; threshold flagged ACP.
models=[('CCDA1','MIK3493667.1',7053,S/'reference_15_Pfam17_relaxed.domtblout'),('Microbispora','WP_168066089.1',7096,S/'reference_15_Pfam17_relaxed.domtblout'),('GnmT','ATY69569.1',7422,S/'GnmT_Pfam17_relaxed.domtblout')]
styles={'PF00109':('#315C7C','KS'),'PF08659':('#8E65A5','KR'),'PF00550':('#15847D','ACP'),'PF09924':('#B47835','PF09924'),'PF01212':('#C35E57','PF01212'),'PF00561':('#6A8498','PF00561')}
fig,ax=plt.subplots(figsize=(12.4,4.15))
for j,(label,acc,ln,path) in enumerate(models):
 y=2-j*.95
 ax.plot([0,ln],[y,y],color='#D3DDE3',lw=13,solid_capstyle='round',zorder=1)
 hits=[]
 for line in path.read_text().splitlines():
  if not line or line.startswith('#'):continue
  v=line.split(); fam=v[1].split('.')[0]
  if v[3]!=acc or fam not in styles:continue
  a,b=int(v[17]),int(v[18]);score=float(v[13])
  if fam=='PF00550' and (score<15 or b-a+1<45): continue
  if fam!='PF00550' and score<20:continue
  hits.append((a,b,score,fam))
 for a,b,score,fam in hits:
  color,_=styles[fam]
  ax.add_patch(Rectangle((a,y-.145),b-a+.5,.29,facecolor=color,edgecolor='white' if not (fam=='PF00550' and score<20.8) else '#94380D',linewidth=1 if score>=20.8 else 1.5,linestyle='-' if score>=20.8 else '--',zorder=5))
 ax.text(-200,y,label,ha='right',va='center',fontweight='bold',fontsize=10)
 ax.text(ln+70,y,f'{ln:,} aa',ha='left',va='center',fontsize=8,color='#596976')
ax.set_xlim(-600,7850);ax.set_ylim(-.5,2.5);ax.set_xticks(np.arange(0,8000,1000));ax.set_yticks([])
ax.set_xlabel('Amino-acid position within protein');ax.set_title('Profile-domain organization of large PKS-associated homologues',loc='left',fontweight='bold',pad=15)
ax.spines[['top','right','left']].set_visible(False);ax.grid(axis='x',alpha=.1);ax.set_axisbelow(True)
handles=[Patch(facecolor=v[0],label=v[1]) for v in styles.values()]
handles.append(Patch(facecolor=styles['PF00550'][0],edgecolor='#94380D',linestyle='--',label='ACP below Pfam GA score'))
ax.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,-.41),ncol=7,fontsize=8,frameon=False)
fig.subplots_adjust(bottom=.26);save(fig,'Figure_5_domain_architecture')

# Figure 5. Terminal ACP-like motif alignment 15 homologues + published GnmT.
sto=(S/'terminal_acp_alignment_15refs.sto').read_text().splitlines()
seqs=[]
for line in sto:
 t=line.split()
 if len(t)==2 and t[0].endswith('_terminal_ACP'):seqs.append((t[0].replace('_terminal_ACP',''),t[1]))
assert len(seqs)==15
for i,(lab,s) in enumerate(seqs):
 assert s[30]=='S',(lab,s[27:32])
seqs.append(('ATY69569.1 (GnmT)',next(line.split()[1] for line in (S/'terminal_acp_alignment_3refs.sto').read_text().splitlines() if line.startswith('GnmT_ACP6'))))
fig,ax=plt.subplots(figsize=(8.0,6.1))
for idx,(nm,al) in enumerate(seqs):
 motif=al[27:32]
 for j,c in enumerate(motif):
  color= '#D3E8E3' if c in 'ST' else '#F0DDBF' if c in 'GA' else '#E4E9ED'
  if j==3:color='#11756E'
  ax.add_patch(Rectangle((j,idx),.96,.92,facecolor=color,edgecolor='white'))
  ax.text(j+.48,idx+.46,c,ha='center',va='center',fontsize=10,fontweight='bold',color='white' if j==3 else '#203344')
ax.set_xlim(0,5);ax.set_ylim(0,len(seqs));ax.invert_yaxis()
ax.set_yticks(np.arange(len(seqs))+.45,[x[0] for x in seqs],fontsize=8)
ax.set_xticks(np.arange(5)+.48,['1','2','3','4*','5']);ax.set_xlabel('Five-residue window surrounding candidate serine; * = aligned S')
ax.set_title('Conserved noncanonical TSF context in terminal ACP-like regions',fontweight='bold',loc='left',pad=14)
for sp in ax.spines.values():sp.set_visible(False)
fig.subplots_adjust(left=.33,bottom=.12);save(fig,'Figure_6_terminal_acp_motifs')

# Figure 6. 18 loci synteny status / target colinear coverage, no statistical testing.
syn=pd.read_csv(S/'comparative_synteny_18.tsv',sep='\t');cov=pd.read_csv(S/'target_colinear_18.tsv',sep='\t')
d=syn.merge(cov[['Assembly','QueryCoveragePercent']],on='Assembly',validate='one_to_one')
rows=[]
for _,r in d.iterrows():
 vals=r['ObservedPattern'].split(',')
 assert len(vals)==14
 ref=['DUF6039','MalonylT','DUF6039','VOC','DUF3224','PIGL','SAM','Hydrolase','ACP','AMP','TARGET','HMG','PPbinding','KS_like']
 rows.append([0 if x=='MISSING' else (2 if x==y else 1) for x,y in zip(vals,ref)])
a=np.array(rows)
from matplotlib.colors import ListedColormap,BoundaryNorm
fig,(ax,bx)=plt.subplots(1,2,figsize=(12.3,6),gridspec_kw={'width_ratios':[4.8,1.55]},sharey=True)
cmap=ListedColormap(['#EAECEF','#D99A64','#237C72']);ax.imshow(a,cmap=cmap,aspect='auto',vmin=0,vmax=2)
for i in range(len(d)):
 for j in range(14):
  if a[i,j]!=2:ax.text(j,i,'–' if a[i,j]==0 else '!',ha='center',va='center',fontsize=7)
ax.set_xticks(range(14),[f'{i:+d}' for i in range(-10,4)],fontsize=8)
ax.set_yticks(range(len(d)),[x.replace('GCF_','').replace('GCA_','') for x in d.Assembly],fontsize=8)
ax.set_xlabel('Relative gene position, target = 0');ax.set_title('Annotation-level synteny (18 assembly loci)',loc='left',fontweight='bold',pad=12)
for sp in ax.spines.values():sp.set_visible(False)
ax.set_xticks(np.arange(-.5,14,1),minor=True);ax.set_yticks(np.arange(-.5,18,1),minor=True);ax.grid(which='minor',color='white',lw=.9);ax.tick_params(which='minor',bottom=False,left=False)
bx.barh(np.arange(len(d)),d.QueryCoveragePercent,color=['#CE8D53' if x<90 else '#2E8079' for x in d.QueryCoveragePercent],height=.68)
for k,x in enumerate(d.QueryCoveragePercent):bx.text(min(float(x)+2,102),k,f'{x:.1f}%',va='center',fontsize=7)
bx.set_xlim(0,115);bx.set_xlabel('Colinear query coverage (%)');bx.set_title('Protein coverage',loc='left',fontweight='bold',pad=12)
bx.spines[['top','right','left']].set_visible(False);bx.tick_params(axis='y',length=0)
fig.legend(handles=[Patch(facecolor='#237C72',label='Exact class match'),Patch(facecolor='#D99A64',label='Other annotated class'),Patch(facecolor='#EAECEF',label='Unavailable position')],loc='lower center',bbox_to_anchor=(.5,-.01),ncol=3,frameon=False,fontsize=8)
fig.subplots_adjust(wspace=.07,bottom=.13,left=.18,right=.96);save(fig,'Figure_4_synteny_homology_18_loci')

print('Reproduced figures:',sorted(x.name for x in F.glob('*.svg')))
