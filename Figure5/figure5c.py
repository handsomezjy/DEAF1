import matplotlib.font_manager as font_manager
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import style
from matplotlib.colors import LinearSegmentedColormap
from scipy.ndimage import gaussian_filter

sample_labels=["WT_CTCF","DEAF1_CTCF","WT_DEAF1","DEAF1_DEAF1"]
sample_boundaries=[0,100,200,300,400]
ma=pd.read_table("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/Domains/Cooltools/InsulationScores/heatmap_TADboundaries.CTCF_DEAF1.K562.GRCh38.gz",header=None,skiprows=1)

style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
plt.rcParams['font.sans-serif']='Helvetica'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['axes.unicode_minus'] =False
fig = plt.figure(figsize=(2.5, 3))
gs = gridspec.GridSpec(2, 3,  height_ratios=[1,3],width_ratios=[1,1,0.1],hspace=0.05,wspace=0.4)
colors1 = [(0, '#6EACDA'), 
          (0.5,'#EEEEEE'),   
          (1, '#ED3500')] 
ccmap=LinearSegmentedColormap.from_list('custom_cmap', colors1)
tmp1=ma.iloc[:,6:106]
tmp2=ma.iloc[:,106:206]
for i in range(2):
    ax = fig.add_subplot(gs[0,i])
    if i==0:
        tmp=tmp1.copy()#.sort_index(ascending=False)
        plt.plot(list(range(1,101)),tmp.mean(),color='#7F8487', linewidth = "1.5")
        plt.title("WT",x=0.5,rotation=0,y=1.0,fontsize=10,color='#7F8487')
    else:
        tmp=tmp2.copy()#.sort_index(ascending=False)
        plt.plot(list(range(1,101)),tmp.mean(),color='#8E1616', linewidth = "1.5")
        plt.title("DEAF1 mutant",x=0.5,rotation=0,y=1.0,fontsize=10,color='#8E1616')
        plt.yticks([0.7,1.2],[])
    plt.ylim(0.7,1.2)
    plt.yticks([0.7,1.2],fontsize=8)
    plt.xticks([1,50.5,100],[])
    plt.grid(axis='both', color='grey', linestyle=':',zorder=0)

    ax=fig.add_subplot(gs[1,i])
    smoothed_tmp = gaussian_filter(tmp, sigma=4)
    heatmap=sns.heatmap(np.array(smoothed_tmp),cmap=ccmap,rasterized=True,
                    vmin=0.2,vmax=1.5,
                    cbar=False,
                    xticklabels=False,
                    yticklabels=False)
    if i==0:
        plt.ylabel('TAD boundaries',labelpad=1,color='k')
    for _, spine in ax.spines.items():
        spine.set_visible(True)  # 显示边框
        spine.set_color('black')  # 设置边框颜色
        spine.set_linewidth(1)  # 设置边框宽度
    plt.xticks([1,50.5,100],['-0.5Mb','center','0.5Mb'],rotation=45)
im1 = heatmap.collections[0]
cax1 = fig.add_subplot(gs[1,2])
cbar1=fig.colorbar(im1, cax=cax1)  
cbar1.set_ticks([0.5,1,1.5])   
cbar1.set_label('CTCF ChIP-seq\nsignal(RPGC)')
plt.savefig('/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/Domains/Cooltools/InsulationScores/Pdf/TADboundaries.CTCF.heatmap.pdf',   ##
                bbox_inches = 'tight',
                facecolor='w')   
plt.show()
