import matplotlib.font_manager as font_manager
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd

up=pd.read_csv("/data/zhangjy/DEAF1/RNAseq_Analysis/DataSort/UpGenes.GOenrichment_BPresults.csv")
tmp=up['GeneRatio'].str.split('/',expand=True)
up['GeneRatio']=tmp[0].astype('int')/tmp[1].astype('int')
down=pd.read_csv("/data/zhangjy/DEAF1/RNAseq_Analysis/DataSort/DownGenes.GOenrichment_BPresults.csv")
tmp=down['GeneRatio'].str.split('/',expand=True)
down['GeneRatio']=tmp[0].astype('int')/tmp[1].astype('int')

style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
plt.rcParams['font.sans-serif']='Helvetica'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['axes.unicode_minus'] =False
fig = plt.figure(figsize=(2.2, 5))
gs = gridspec.GridSpec(4, 2,  height_ratios=[1,1,1,1],hspace=0.1,width_ratios=[1,0.1],wspace=0.3)
ax1=fig.add_subplot(gs[0:2,0])
colors1 = [(0, '#ED3500'), 
          (0.5,'#FFA673'),   # 红色，位置 0
          (1, '#EEEEEE')] 
colors2 = [(0, '#141E61'), 
          (0.5,'#6E85B2'),   # 红色，位置 0
          (1, '#EEEEEE')] 
ccmap=LinearSegmentedColormap.from_list('custom_cmap', colors1)
uptmp=up.loc[0:4,:].iloc[::-1]
sc1=plt.scatter(uptmp.loc[:,'GeneRatio'],np.arange(1,6),linewidths=0.5,edgecolors='k',
            s=uptmp.loc[:,'Count']*6,c=uptmp.loc[:,'p.adjust']
            ,cmap=ccmap,vmin=0.008,vmax=0.02)
plt.yticks(np.arange(1,6),uptmp.loc[:,'Description'])
plt.ylim(0.5,5.5)
plt.xlim(0,0.081)
plt.xticks(np.arange(0,0.1,0.02),['','','','',''])
plt.ylabel('Upregulated Genes',fontsize=12,color='#ED3500')
plt.grid(axis='both', color='grey', linestyle=':',zorder=0)
ax12=fig.add_subplot(gs[0,1])
cbar1 = fig.colorbar(sc1, cax=ax12,aspect=10,pad=30)
cbar1.set_label('p.adjust', fontsize=9)
cbar1.ax.tick_params(labelsize=8)
cbar1.set_ticks((0.008,0.015,0.02),(r'$2e^{-2}$',r'$15e^{-2}$',r'$8e^{-2}$'))#,(r'$2e^{-2}$',r'$15e^{-2}$',r'$8e^{-2}$')
cbar1.set_ticklabels([r'$8e^{-3}$',r'$1.5e^{-3}$',r'$2e^{-2}$'])

ccmap=LinearSegmentedColormap.from_list('custom_cmap', colors2)
ax2=fig.add_subplot(gs[2:,0])
downtmp=down.loc[[0,1,2,34,42],:].iloc[::-1]
sc2=plt.scatter(downtmp.loc[:,'GeneRatio'],np.arange(1,6),linewidths=0.5,edgecolors='k',
            s=downtmp.loc[:,'Count']*6,c=downtmp.loc[:,'p.adjust'],cmap=ccmap)
plt.yticks(np.arange(1,6),downtmp.loc[:,'Description'])
plt.ylim(0.5,5.5)
plt.xlim(0,0.081)
plt.xticks(np.arange(0,0.1,0.02))#,('0',r'$2e^{-2}$',r'$4e^{-2}$',r'$6e^{-2}$',r'$8e^{-2}$')
plt.grid(axis='both', color='grey', linestyle=':',zorder=0)
plt.xlabel('GeneRatio',fontsize=12)
plt.ylabel('Downregulated Genes',fontsize=12,color='#6E85B2')
ax21=fig.add_subplot(gs[2,1])
cbar2 = fig.colorbar(sc2, cax=ax21,aspect=10,pad=30)
cbar2.set_label('p.adjust', fontsize=9)
cbar2.ax.tick_params(labelsize=9)
cbar2.set_ticks((0.006,0.003,8e-5))
cbar2.set_ticklabels([r'$6e^{-3}$',r'$3e^{-3}$',r'$8e^{-5}$'])


all_counts = list(up.loc[0:4, 'Count']) + list(down.loc[0:4, 'Count'])
min_count, max_count = min(all_counts), max(all_counts)
count_values = [5,15,30]

handles = [
    plt.scatter([], [], s=count*6, color='w', alpha=0.7,linewidths=1,edgecolors='k')
    for count in count_values
]
ax_legend=fig.add_subplot(gs[3,1])
ax_legend.axis('off') 
ax_legend.legend(
    handles, [str(round(c, 2)) for c in count_values],
    loc='center',title='Count',
    ncol=1,frameon=False,
    fontsize=10,
    title_fontsize=10,
    labelspacing=1,
    handletextpad=0.2#,labelpad=2
)
plt.tight_layout(pad=0.8)
plt.savefig('/data/zhangjy/DEAF1/RNAseq_Analysis/Pdf/UpDowngenes.GO_BP.scatter.pdf',   ##
                bbox_inches = 'tight',
                facecolor='w')   
plt.show()
