import matplotlib.font_manager as font_manager
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import style
from matplotlib.colors import LinearSegmentedColormap

style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
plt.rcParams['font.sans-serif']='Helvetica'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['axes.unicode_minus'] =False
fig = plt.figure(figsize=(2.5, 2))
gs = gridspec.GridSpec(2, 3,  height_ratios=np.ones(2),
                       width_ratios=[1,1,0.1],hspace=0.05,wspace=0.05)
ma=["/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/APA_results/Difftypes/WT_loopStrenghthened/5000/gw/normedAPA.txt",
    "/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/APA_results/Difftypes/KO_loopStrenghthened/5000/gw/normedAPA.txt",
    "/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/APA_results/Difftypes/WT_loopWeakened/5000/gw/normedAPA.txt",
    "/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/APA_results/Difftypes/KO_loopWeakened/5000/gw/normedAPA.txt"]
k=0
colors = [(0, '#EEEEEE'), 
                (0.5,'#E76F2E'),
                (1, '#3E2C23')] 
ccmap=LinearSegmentedColormap.from_list('custom_cmap', colors)
xlabels=['WT','DEAF1-mutant']
ylabels=['Strenghthened','Weakened']
for i in range(2):
    for j in range(2):
        ax=fig.add_subplot(gs[i,j])
        data=read_apa(ma[k])
        k+=1
        heatmap=sns.heatmap(data,cmap=ccmap,rasterized=True,
                        vmin=0,vmax=3.1,
                        cbar=False,
                        xticklabels=False,
                        yticklabels=False,square=True)
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(1)
            spine.set_color('black')
        data=pd.DataFrame(data)
        center_mean = data.iloc[10,10]/np.mean(data.iloc[-5:,:5])
        val_str = f"{center_mean:.2f}"
        plt.text(5, 17, 'P2LL\n'+val_str, fontsize=10, color='k', ha='center', va='center')
        if i==0:
            plt.title(xlabels[j],fontsize=10)
        if j==0:
            plt.ylabel(ylabels[i],fontsize=10)
ax.set_xticks([0,10,20])
ax.set_xticklabels(['-50 kb','0','50 kb'], fontsize=9)
ax.tick_params(axis='x',bottom=True,top=False,labelbottom=True,length=2,width=0.8,pad=1)
cax = fig.add_axes([0.25, -0.05, 0.4, 0.03]) #left,bottom,width,height
cbar = plt.colorbar(
    heatmap.collections[0],
    cax=cax,
    orientation='horizontal'
)
cbar.outline.set_linewidth(0.8)
cbar.ax.tick_params(labelsize=9,width=0.8,length=2)
cbar.set_ticks([0,1,2,3])
cbar.ax.text(1.05, 0.5,'Obs/Exp',transform=cbar.ax.transAxes,ha='left',va='center',fontsize=9)

plt.savefig('/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/Pdf/StrengthededWeakened.loops.heatmap.pdf',   ##
        bbox_inches = 'tight',
        facecolor='w')   
