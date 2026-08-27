import numpy as np
import pandas as pd
import scipy.stats as stats
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns

wt=pd.read_table("/data/zhangjy/DEAF1/RNAseq_Analysis/DataSort/Genes_tss.WT_histones.matrix.gz",skiprows=1,header=None)
sample_labels=["H3K27ac","H3K27me3","H3K4me1","H3K4me3","H3K9me3","RPB1"]
group_labels=["Upregulated\ngenes","Downregulated\ngenes"]
sample_boundaries=[0,400,800,1200,1600,2000,2400]
group_boundaries=[0,3820,9149]
ko=pd.read_table("/data/zhangjy/DEAF1/RNAseq_Analysis/DataSort/Genes_tss.KO_histones.matrix.gz",skiprows=1,header=None)

def histone_TSS(z):
    style.use('default')
    font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
    plt.rcParams['font.sans-serif']='Helvetica'
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['axes.unicode_minus'] = False

    fig = plt.figure(figsize=(4, 4))
    main_gs = gridspec.GridSpec(1, 4, width_ratios=[1, 1, 1, 0.1], wspace=0.4)
    left_gs = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=main_gs[0, 0], 
                                            height_ratios=[2,3], hspace=0.05)
    middle_gs = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=main_gs[0, 1], 
                                            height_ratios=[2,3], hspace=0.05)
    right_gs = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=main_gs[0, 2], 
                                            height_ratios=[1, 1, 1], hspace=0.05)
    bar_gs = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=main_gs[0, 3], 
                                            height_ratios=[1, 1, 1], hspace=0.05)
    gs_=[left_gs,middle_gs,right_gs]
    colors1 = [(0, '#EEEEEE'), 
            (0.5,'#FFA673'),   # 红色，位置 0
            (1, '#ED3500')] 
    colors2 = [(0, '#EEEEEE'), 
            (0.5,'#6E85B2'),   # 红色，位置 0
            (1, '#141E61')] 
    colors=[colors1,colors2]

    tmp1=wt.iloc[group_boundaries[0]:group_boundaries[1],(6+z*400):(6+(z+1)*400)]
    tmp2=ko.iloc[group_boundaries[0]:group_boundaries[1],(6+z*400):(6+(z+1)*400)]
    tmp3=wt.iloc[group_boundaries[1]:,(6+z*400):(6+(z+1)*400)]
    tmp4=ko.iloc[group_boundaries[1]:,(6+z*400):(6+(z+1)*400)]
    tmpdf=[tmp1,tmp3,tmp2,tmp4]
    ii=0
    ims=[]
    for i in range(2):
        if i==0:
            ax=fig.add_subplot(gs_[2][i,0])
        else:
            ax=fig.add_subplot(gs_[2][i+1,0])
        plt.plot(list(range(1,401)),tmpdf[i].mean(),color='#7F8487', linewidth = "1.5")
        plt.plot(list(range(1,401)),tmpdf[i+2].mean(),color='#8E1616', linewidth = "1.5")
        plt.ylim(-0.1,0.45)
        plt.yticks([0,0.2,0.4], fontsize=8)
        if i==1:
            plt.xticks([1,200.5,400],['-2.0','TSS','2.0kb'],rotation=0)#
        else:
            plt.xticks([1,200.5,400],['','',''],rotation=0)
        plt.grid(axis='both', color='grey', linestyle=':',zorder=0)
        for j in range(2):
            ax=fig.add_subplot(gs_[i][j,0])
            ccmap=LinearSegmentedColormap.from_list('custom_cmap', colors[j])
            print(f"数据类型: {tmpdf[ii].shape}")
            smoothed_tmp1 = gaussian_filter(tmpdf[ii], sigma=1.5)
            print(ii)
            # ii+=1
            heatmap=sns.heatmap(np.array(smoothed_tmp1),cmap=ccmap,rasterized=True,
                        vmin=-0.2,vmax=0.8,
                        cbar=False,
                        xticklabels=False,
                        yticklabels=False)
            if j==1:
                plt.xticks([0,199.5,399],['-2.0','TSS','2.0kb'],fontsize=10,rotation=0)
            for _, spine in ax.spines.items():
                spine.set_visible(True)  # 显示边框
                spine.set_color('black')  # 设置边框颜色
                spine.set_linewidth(0.8)  # 设置边框宽度
            if i==0:
                plt.ylabel(group_labels[j],labelpad=30,color=colors[j][2][1],rotation=0)
            ims.append(heatmap.collections[0])
            ii+=1

    cax1 = fig.add_subplot(bar_gs[0,0])  
    cax2 = fig.add_subplot(bar_gs[2,0])  
    cbar1=fig.colorbar(ims[0], cax=cax1)
    cbar1.set_ticks([-0.2, 0, 0.8])
    cbar2=fig.colorbar(ims[1], cax=cax2)
    cbar2.set_ticks([-0.2, 0, 0.8])

    bbox00 = left_gs[0,0].get_position(fig)
    center_x1 = (bbox00.x0 + bbox00.x1) / 2
    top_y= bbox00.y1+0.02
    bbox01 = middle_gs[0,0].get_position(fig)
    center_x2 = (bbox01.x0 + bbox01.x1) / 2
    center_x3 = (bbox00.x0 + bbox01.x1) / 2
    title_y=bbox00.y1+0.07
    fig.text(center_x1, top_y, 'WT', ha='center', va='bottom', fontsize=12, fontweight='bold')
    fig.text(center_x2, top_y, 'DEAF1 mutant', ha='center', va='bottom', fontsize=12, fontweight='bold',color='#8E1616')    
    fig.text(center_x3, title_y, sample_labels[z], ha='center', va='bottom', fontsize=14, fontweight='bold')    

    plt.savefig('/data/zhangjy/DEAF1/RNAseq_Analysis/Pdf/Heatmap_TSS.{x}.lineplot.pdf'.format(x=sample_labels[z]),   ##
                    bbox_inches = 'tight',
                    facecolor='w')   
    plt.show()

for z in range(len(sample_labels)):
    print(sample_labels[z])
    histone_TSS(z)
