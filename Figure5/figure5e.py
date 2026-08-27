import matplotlib.font_manager as font_manager
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import style
from matplotlib.colors import LinearSegmentedColormap
from scipy.ndimage import gaussian_filter

sample_labels=["WT.insulation_10kb","KO.insulation_10kb"]
sample_boundaries=[0,100,200]
ma=pd.read_table("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/Domains/Cooltools/InsulationScores/heatmap_TADboundaries.InsulationScores.K562.GRCh38.gz",header=None,skiprows=1)


style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
plt.rcParams['font.sans-serif']='Helvetica'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['axes.unicode_minus'] =False
fig = plt.figure(figsize=(2.5, 3))
gs = gridspec.GridSpec(2, 3,  height_ratios=[1,3],width_ratios=[1,1,0.1],hspace=0.05,wspace=0.4)
colors1 = [(0, '#F5824A'),
           (0.5, '#FEF3E2'),   
          (1, '#254F22')] 
ccmap=LinearSegmentedColormap.from_list('custom_cmap', colors1)
tmp1=ma.iloc[:,6:106]
tmp2=ma.iloc[:,106:]
for i in range(2):
    ax = fig.add_subplot(gs[0,i])
    if i==0:
        tmp=tmp1.copy().sort_index(ascending=False)
        plt.plot(list(range(1,101)),tmp.mean(),color='#7F8487', linewidth = "1.7")
        plt.title("WT",x=0.5,rotation=0,y=1.0,fontsize=10,color='#7F8487')
    else:
        tmp=tmp2.copy().sort_index(ascending=False)
        plt.plot(list(range(1,101)),tmp.mean(),color='#8E1616', linewidth = "1.7")
        plt.title("DEAF1 mutant",x=0.5,rotation=0,y=1.0,fontsize=10,color='#8E1616')
        plt.yticks([-0.4,-0.2,0],[])
    plt.ylim(-0.55,0.2)
    plt.yticks([-0.4,-0.2,0],fontsize=8)
    plt.xticks([1,50.5,100],[])
    plt.grid(axis='both', color='grey', linestyle=':',zorder=0)

    ax=fig.add_subplot(gs[1,i])
    smoothed_tmp = gaussian_filter(tmp, sigma=4)
    heatmap=sns.heatmap(np.array(smoothed_tmp),cmap=ccmap,rasterized=True,
                    vmin=-0.8,vmax=0.6,
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
cbar1.set_ticks([-0.8,-0.4,0,0.4])   
cbar1.set_label('Insulation scores')
plt.savefig('/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/Domains/Cooltools/InsulationScores/Pdf/TADboundaries.Insulation.heatmap.pdf',   ##
                bbox_inches = 'tight',
                facecolor='w')   
plt.show()

