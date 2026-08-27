import numpy as np
import pandas as pd
import scipy.stats as stats
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns

ma=pd.read_table("/data/zhangjy/DEAF1/ChiPseq_Analysis/ResultsSort/Deeptools/InputNew/Deeptools/mergedPeaks/Merged.CTCF_DEAF1.heatmapData.gz",header=None,skiprows=1)
sample_labels=["WT_CTCF","DEAF1_CTCF","WT_DEAF1","DEAF1_DEAF1"]
group_labels=["CTCF","DEAF1"]
sample_boundaries=[0,200,400,600,800]
group_boundaries=[0,69230,74359]
ctcfdeaf1=pd.read_table("/data/zhangjy/DEAF1/ChiPseq_Analysis/ResultsSort/Deeptools/InputNew/Deeptools/mergedPeaks/merged_CTCF_intersectDEAF1.bed",header=None)
ctcfnodeaf1=pd.read_table("/data/zhangjy/DEAF1/ChiPseq_Analysis/ResultsSort/Deeptools/InputNew/Deeptools/mergedPeaks/merged_CTCF_noDEAF1.bed",header=None)
deaf1noctcf=pd.read_table("/data/zhangjy/DEAF1/ChiPseq_Analysis/ResultsSort/Deeptools/InputNew/Deeptools/mergedPeaks/merged_DEAF1_noCTCF.bed",header=None)
deaf1ctcf=pd.read_table("/data/zhangjy/DEAF1/ChiPseq_Analysis/ResultsSort/Deeptools/InputNew/Deeptools/mergedPeaks/merged_DEAF1_intersectCTCF.bed",header=None)

style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
plt.rcParams['font.sans-serif']='Helvetica'
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Helvetica'
plt.rcParams['mathtext.it'] = 'Helvetica:italic'
plt.rcParams['mathtext.bf'] = 'Helvetica:bold'
plt.rcParams['axes.unicode_minus'] =False
plt.rcParams['pdf.fonttype'] = 42
# 高斯滤波平滑矩阵
fig = plt.figure(figsize=(4, 3))
gs = gridspec.GridSpec(3, 4,  height_ratios=[1,1,1],width_ratios=[1,1,1,0.1],hspace=0.05,wspace=0.4)
# colors1 = [(0, '#6EACDA'), 
#           (0.5,'#EEEEEE'),   
#           (1, '#ED3500')] 
colors1 = [(0, '#EEEEEE'), 
          (0.5,'#ED3500'),   
          (1, '#7F2020')] 
ccmap=LinearSegmentedColormap.from_list('custom_cmap', colors1)
gs_col1 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[:, 0], 
                                        #   height_ratios=[ctcfdeaf1.shape[0], ctcfnodeaf1.shape[0]],hspace=0.05)
                                        height_ratios=[1,10],hspace=0.05)
gs_col2 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[:, 1], 
                                        #   height_ratios=[ctcfdeaf1.shape[0], ctcfnodeaf1.shape[0]], hspace=0.05)
                                        height_ratios=[1,10],hspace=0.05)
gs_col3 = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs[:, 2], 
                                          height_ratios=[1, 1, 1], hspace=0.05)
gs_col4 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[:, 3], 
                                          height_ratios=[1, 2], hspace=0.05)
tmp1=pd.merge(ma,ctcfdeaf1,on=[0,1,2])
tmp2=pd.merge(ma,ctcfnodeaf1,on=[0,1,2])
tmpli=[tmp1,tmp2]
gscol_li=[gs_col1,gs_col2,gs_col3,gs_col4]
for i in range(2):
    gs_=gscol_li[i]
    for j in range(2):
        tmp=tmpli[j]
        ax = fig.add_subplot(gs_[j])
        df=tmp.iloc[:,(6+200*i):(6+(i+1)*200)]
        smoothed_tmp = gaussian_filter(df, sigma=5)
        heatmap=sns.heatmap(np.array(smoothed_tmp),cmap=ccmap,rasterized=True,
            vmin=0,vmax=15,
            cbar=False,
            xticklabels=False,
            yticklabels=False)
        for _, spine in ax.spines.items():
            spine.set_visible(True)  
            spine.set_color('black')  
            spine.set_linewidth(1)  
        if (i==0)&(j==0):
            ax.set_ylabel('$DEAF1^{+}$ $CTCF^{+}$\nn='+str(df.shape[0]),labelpad=1)
            plt.title('WT',x=0.5,rotation=0,y=1.01,fontsize=10,color='#7F8487')
            ax.yaxis.set_label_coords(-0.25, 0.5)
        elif (i==0)&(j==1):
            ax.set_ylabel('$DEAF1^{-}$ $CTCF^{+}$\nn='+str(tmp.shape[0]),labelpad=1)
            ax.yaxis.set_label_coords(-0.25, 0.5)
        elif (i==1)&(j==0):
            plt.title("DEAF1 mutant",x=0.5,rotation=0,y=1.01,fontsize=10,color='#8E1616')
        if j==1:
            plt.xticks([1,100.5,200],['-1.0','center','1.0kb'],rotation=45)
    df1=tmpli[i].iloc[:,6:206]
    df2=tmpli[i].iloc[:,206:406]
    if i>0:
        i+=1
    ax = fig.add_subplot(gs_col3[i])  
    plt.plot(list(range(1,201)),df1.mean(),color='#7F8487', linewidth = "1.5")
    plt.plot(list(range(1,201)),df2.mean(),color='#8E1616', linewidth = "1.5")
    plt.ylim(0,18)
    plt.yticks([0,5,10,15], fontsize=7)
    plt.xticks([1,100.5,200],['-1.0','center','1.0kb'],rotation=45)#
    # plt.grid(axis='both', color='grey', linestyle=':',zorder=0)

im1 = heatmap.collections[0]
cax1 = fig.add_subplot(gs_col4[1])
cbar1=fig.colorbar(im1, cax=cax1)
cbar1.set_ticks([0,5,10,15])

bbox00 = gs_col1[0,0].get_position(fig)
bbox01 = gs_col2[0,0].get_position(fig)
center_x1 = (bbox00.x0 + bbox01.x1) / 2
top_y = gs_col1[0,0].get_position(fig).y1
title_y = top_y + 0.1
fig.text(center_x1, title_y, 'CTCF ChIP-seq', ha='center', va='bottom', fontsize=12, fontweight='bold')
plt.savefig('/data/zhangjy/DEAF1/ChiPseq_Analysis/ResultsSort/Deeptools/InputNew/Deeptools/mergedPeaks/Pdf/Heatmap_CTCF.RPGC.heatmap.pdf',   ##
                bbox_inches = 'tight',
                facecolor='w')   
plt.show()

style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
plt.rcParams['font.sans-serif']='Helvetica'
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Helvetica'
plt.rcParams['mathtext.it'] = 'Helvetica:italic'
plt.rcParams['mathtext.bf'] = 'Helvetica:bold'
plt.rcParams['axes.unicode_minus'] =False
plt.rcParams['pdf.fonttype'] = 42
# 高斯滤波平滑矩阵
fig = plt.figure(figsize=(4, 3))
gs = gridspec.GridSpec(3, 4,  height_ratios=[1,1,1],width_ratios=[1,1,1,0.1],hspace=0.05,wspace=0.4)
# colors1 = [(0, '#6EACDA'), 
#           (0.5,'#EEEEEE'),   
#           (1, '#471396')] 
colors1 = [(0, '#EEEEEE'), 
          (0.5,'#294669'),   
          (1, '#471396')] 
ccmap=LinearSegmentedColormap.from_list('custom_cmap', colors1)
gs_col1 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[:, 0], 
                                          height_ratios=[deaf1ctcf.shape[0], deaf1noctcf.shape[0]],hspace=0.05)
gs_col2 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[:, 1], 
                                          height_ratios=[deaf1ctcf.shape[0], deaf1noctcf.shape[0]], hspace=0.05)
gs_col3 = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs[:, 2], 
                                          height_ratios=[1, 1, 1], hspace=0.05)
gs_col4 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[:, 3], 
                                          height_ratios=[1, 2], hspace=0.05)
tmp1=pd.merge(ma,deaf1ctcf,on=[0,1,2])
tmp2=pd.merge(ma,deaf1noctcf,on=[0,1,2])
tmpli=[tmp1,tmp2]
gscol_li=[gs_col1,gs_col2,gs_col3,gs_col4]
for i in range(2):
    gs_=gscol_li[i]
    for j in range(2):
        if j==0:
            tmp=tmpli[j]
            tmp=tmp[tmp[6]<15]
        else:
            tmp=tmpli[j]
        ax = fig.add_subplot(gs_[j])
        df=tmp.iloc[:,(6+200*i):(6+(i+1)*200)]
        smoothed_tmp = gaussian_filter(df, sigma=5)
        heatmap=sns.heatmap(np.array(smoothed_tmp),cmap=ccmap,rasterized=True,
            vmin=0,vmax=15,
            cbar=False,
            xticklabels=False,
            yticklabels=False)
        for _, spine in ax.spines.items():
            spine.set_visible(True)  
            spine.set_color('black')  
            spine.set_linewidth(1)  
        if (i==0)&(j==0):
            ax.set_ylabel('$DEAF1^{+}$ $CTCF^{+}$\nn='+str(2647),labelpad=1)
            plt.title('WT',x=0.5,rotation=0,y=1.01,fontsize=10,color='#7F8487')
            ax.yaxis.set_label_coords(-0.25, 0.5)
        elif (i==0)&(j==1):
            ax.set_ylabel('$DEAF1^{+}$ $CTCF^{-}$\nn='+str(tmp.shape[0]),labelpad=1)
            ax.yaxis.set_label_coords(-0.25, 0.5)
        elif (i==1)&(j==0):
            plt.title("DEAF1 mutant",x=0.5,rotation=0,y=1.01,fontsize=10,color='#8E1616')
        if j==1:
            plt.xticks([1,100.5,200],['-1.0','center','1.0kb'],rotation=45)
    df1=tmpli[i].iloc[:,6:206]
    df2=tmpli[i].iloc[:,206:406]
    if i>0:
        i+=1
    ax = fig.add_subplot(gs_col3[i])  
    plt.plot(list(range(1,201)),df1.mean(),color='#7F8487', linewidth = "1.5")
    plt.plot(list(range(1,201)),df2.mean(),color='#8E1616', linewidth = "1.5")
    plt.ylim(0,18)
    plt.yticks([0,5,10,15], fontsize=7)
    plt.xticks([1,100.5,200],['-1.0','center','1.0kb'],rotation=45)#
    # plt.grid(axis='both', color='grey', linestyle=':',zorder=0)

im1 = heatmap.collections[0]
cax1 = fig.add_subplot(gs_col4[1])
cbar1=fig.colorbar(im1, cax=cax1)
cbar1.set_ticks([0,5,10,15])

bbox00 = gs_col1[0,0].get_position(fig)
bbox01 = gs_col2[0,0].get_position(fig)
center_x1 = (bbox00.x0 + bbox01.x1) / 2
top_y = gs_col1[0,0].get_position(fig).y1
title_y = top_y + 0.1
fig.text(center_x1, title_y, 'DEAF1 ChIP-seq', ha='center', va='bottom', fontsize=12, fontweight='bold')
plt.savefig('/data/zhangjy/DEAF1/ChiPseq_Analysis/ResultsSort/Deeptools/InputNew/Deeptools/mergedPeaks/Pdf/Heatmap_DEAF1.RPGC.heatmap.pdf',   ##
                bbox_inches = 'tight',
                facecolor='w')   
plt.show()
