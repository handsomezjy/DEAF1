import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import style
import numpy as np

wt=pd.read_table("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/LoopStrength/WT.5kb.rescored_strength_coords.txt")
ko=pd.read_table("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/LoopStrength/KO.5kb.rescored_strength_coords.txt")
li=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']
df1=pd.merge(wt.loc[:,li+['strength']],ko.loc[:,li+['strength']],on=li, suffixes=('_WT', '_KO'))
df1['log2KO_WT']=np.log2((df1['strength_KO']+1e-6)/(df1['strength_WT']+1e-6))

wt=pd.read_table("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/LoopStrength/WT.10kb.rescored_strength_coords.txt")
ko=pd.read_table("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/LoopStrength/KO.10kb.rescored_strength_coords.txt")
df2=pd.merge(wt.loc[:,li+['strength']],ko.loc[:,li+['strength']],on=li, suffixes=('_WT', '_KO'))
df2['log2KO_WT']=np.log2((df2['strength_KO']+1e-6)/(df2['strength_WT']+1e-6))

df=pd.concat([df1,df2])
df=df[(df['strength_WT']>1)&(df['strength_KO']>1)&((df['strength_WT']<100))&(df['strength_KO']<100)]

style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
plt.rcParams['font.sans-serif']='Helvetica'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['axes.unicode_minus'] =False
fig,ax = plt.subplots(figsize=(3,3),constrained_layout=False)
tmpp=df[(df['log2KO_WT']<0.5)&(df['log2KO_WT']>=-0.5)]
plt.scatter(tmpp['strength_WT'],tmpp['strength_KO'],c='#DDDDDD',
            s=10,alpha=0.5,edgecolors='#DDDDDD', 
            rasterized=True)
plt.text(20,20,s='Unchanged\nn={}'.format(tmpp.shape[0]),c='grey')
tmp=df[(df['log2KO_WT']>=0.5)|(df['log2KO_WT']<=-0.5)]
tmp1=tmp[(tmp['log2KO_WT']>=0.5)]
plt.scatter(tmp1['strength_WT'],tmp1['strength_KO'],c='#F63049',
            s=10,alpha=0.5,edgecolors='#F63049', rasterized=True)
plt.text(5,20,s='Strengthened\nn={}'.format(tmp1.shape[0]),c='#F63049')
tmp2=tmp[(tmp['log2KO_WT']<=-0.5)]
plt.scatter(tmp2['strength_WT'],tmp2['strength_KO'],c='#008BFF',
            s=10,alpha=0.5,edgecolors='#008BFF', rasterized=True)
plt.text(20,5,s='Weakened\nn={}'.format(tmp2.shape[0]),c='#008BFF')
plt.plot([0,30],[0,30*2**0.5],linewidth=1,color='k')
plt.plot([0,30*2**0.5],[0,30],linewidth=1,color='k')
plt.ylim(0,30)
plt.xlim(0,30)
plt.xticks([0,10,20,30])
plt.yticks([0,10,20,30])
plt.xlabel('WT loop strength',fontsize=12)
plt.ylabel('DEAF1-mutant loop strength',fontsize=12)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_color('grey')
ax.spines['bottom'].set_color('grey')
ax.spines['left'].set_position(('outward', 5))
ax.spines['bottom'].set_position(('outward', 5))
plt.grid(axis='both', color='grey', linestyle=':',zorder=0)
# plt.legend(frameon=False,bbox_to_anchor=(0.6, 1), loc='upper left')
plt.savefig('/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/Pdf/LoopStrength.scatter.pdf',   ##
                bbox_inches = 'tight',
                facecolor='w')  

