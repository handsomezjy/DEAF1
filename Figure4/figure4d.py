import numpy as np
import pandas as pd
import scipy.stats as stats
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns

ma=pd.read_table("/data/zhangjy/DEAF1/RNAseq_Analysis/DataSort/Genes_tss.CTCF_DEAF1.matrix.gz",skiprows=1,header=None)
bed1=pd.read_table("/data/zhangjy/DEAF1/RNAseq_Analysis/DataSort/Up_genes_TSS.2kb.merged_CTCF_noDEAF1.bed",header=None)
bed2=pd.read_table("/data/zhangjy/DEAF1/RNAseq_Analysis/DataSort/Down_genes_TSS.2kb.merged_CTCF_noDEAF1.bed",header=None)
upctcf=pd.merge(ma,bed1,on=[0,1,2])
downctcf=pd.merge(ma,bed2,on=[0,1,2])

style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
plt.rcParams['font.sans-serif']='Helvetica'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['axes.unicode_minus'] = False
group_labels=['CTCF-bound\nUp-regulated\ngenes($\pm$ 2kb)','CTCF-bound\nDown-regulated\ngenes($\pm$ 2kb)']
fig = plt.figure(figsize=(4, 4))
main_gs = gridspec.GridSpec(1, 4, width_ratios=[1, 1, 1, 0.1], wspace=0.4)
left_gs = gridspec.GridSpecFromSubplotSpec(2,1 ,subplot_spec=main_gs[0, 0], 
                                        height_ratios=[210,255], hspace=0.08)
middle_gs = gridspec.GridSpecFromSubplotSpec(2,1 ,subplot_spec=main_gs[0, 1], 
                                        height_ratios=[210,255], hspace=0.08)
right_gs = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=main_gs[0, 2], 
                                        height_ratios=[210,255], hspace=0.08)
bar_gs = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=main_gs[0, 3], 
                                        height_ratios=[210,255], hspace=0.08)
colors1 = [(0, '#EEEEEE'), 
            (0.5,'#FFA673'),   # 红色，位置 0
            (1, '#ED3500')] 
colors2 = [(0, '#EEEEEE'), 
        (0.5,'#6E85B2'),   # 红色，位置 0
        (1, '#141E61')] 
colors=[colors1,colors2]
gs_heatmap=[left_gs,middle_gs]

for i in range(2):
    if i==0:
        tmp1=upctcf.iloc[:,6:406]
        tmp2=upctcf.iloc[:,806:1206]
    else:
        tmp1=downctcf.iloc[:,6:406]
        tmp2=downctcf.iloc[:,806:1206]
    tmp=[tmp1,tmp2]
    ims=[]
    for j in range(2):
        ccmap=LinearSegmentedColormap.from_list('custom_cmap', colors[j])
        ax=fig.add_subplot(gs_heatmap[i][j])
        smoothed_tmp = gaussian_filter(tmp[j], sigma=3)
        heatmap=sns.heatmap(np.array(smoothed_tmp),cmap=ccmap,rasterized=True,
                    vmin=0,vmax=0.2,
                    cbar=False,
                    xticklabels=False,
                    yticklabels=False)
        ims.append(heatmap.collections[0])
        if j==1:
            plt.xticks([0,199.5,399],['-2.0','TSS','2.0kb'],fontsize=9,rotation=0)
        for _, spine in ax.spines.items():
            spine.set_visible(True)  # 显示边框
            spine.set_color('black')  # 设置边框颜色
            spine.set_linewidth(0.8)  # 设置边框宽度
        if i==0:
            plt.ylabel(group_labels[j],labelpad=30,color=colors[j][2][1],rotation=0)
for j in range(2):
    ax=fig.add_subplot(right_gs[j])
    if j==0:
        tmp1=upctcf.iloc[:,6:406]
        tmp2=upctcf.iloc[:,806:1206]
    else:
        tmp1=downctcf.iloc[:,6:406]
        tmp2=downctcf.iloc[:,806:1206]
    plt.plot(list(range(1,401)),tmp1.mean(),color='#7F8487', linewidth = "1.5")
    plt.plot(list(range(1,401)),tmp2.mean(),color='#8E1616', linewidth = "1.5")
    plt.ylim(0,0.3)
    plt.yticks([0,0.1,0.2,0.3], fontsize=8)
    if j==1:
        plt.xticks([1,200.5,400],['-2.0','TSS','2.0kb'],rotation=0,fontsize=9)
    else:
        plt.xticks([])

cax1 = fig.add_subplot(bar_gs[0,0])  
cax2 = fig.add_subplot(bar_gs[1,0])  
cbar1=fig.colorbar(ims[0], cax=cax1)
cbar1.set_ticks([ 0, 0.2])
cbar2=fig.colorbar(ims[1], cax=cax2)
cbar2.set_ticks([ 0, 0.2])

bbox00 = left_gs[0,0].get_position(fig)
center_x1 = (bbox00.x0 + bbox00.x1) / 2
top_y= bbox00.y1+0.02
bbox01 = middle_gs[0,0].get_position(fig)
center_x2 = (bbox01.x0 + bbox01.x1) / 2
center_x3 = (bbox00.x0 + bbox01.x1) / 2
title_y=bbox00.y1+0.07
fig.text(center_x1, top_y, 'WT', ha='center', va='bottom', fontsize=12, fontweight='bold')
fig.text(center_x2, top_y, 'DEAF1 mutant', ha='center', va='bottom', fontsize=12, fontweight='bold',color='#8E1616')    
fig.text(center_x3, title_y, 'CTCF', ha='center', va='bottom', fontsize=14, fontweight='bold')  

plt.savefig('/data/zhangjy/DEAF1/RNAseq_Analysis/Pdf/DRE_TSS2kb.CTCF.heatmap.pdf',   ##
                    bbox_inches = 'tight',
                    facecolor='w')   
plt.show()
