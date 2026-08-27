import h5py
import matplotlib.font_manager as font_manager
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib import style
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

f = h5py.File(
    "/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/Domains/Cooltools/ATA/Downsample/WT.boundary.downsample.pileup.npz",
    "r"
)
print(list(f.keys()))
mat = f['data'][:]
style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
plt.rcParams['font.sans-serif']='Helvetica'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['axes.unicode_minus'] =False
fig = plt.figure(figsize=(4/1.2, 3/1.2))
gs = gridspec.GridSpec(3, 2,  width_ratios=[4,0.2],height_ratios=[1,1,2],wspace=0.1)
colors1 = [(0, '#1C4D8D'),
           (0.5, '#FAF3E1'),   
          (1, '#C40C0C')] 
ccmap=LinearSegmentedColormap.from_list('custom_cmap', colors1)
ax = fig.add_subplot(gs[:,0])
heatmap=sns.heatmap(mat,cmap=ccmap,rasterized=True,
                        vmin=0.5,vmax=1.5,
                        cbar=False,
                        xticklabels=False,
                        yticklabels=False,square=True)
for spine in heatmap.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1)
        spine.set_color('black')
rect1 = Rectangle((0,17),13,3,
                  fill=False,
                  edgecolor='black',
                  linewidth=1,linestyle='--')
ax.add_patch(rect1)
rect2 = Rectangle((28,21),13,3,
                  fill=False,
                  edgecolor='black',
                  linewidth=1,linestyle='--')
ax.add_patch(rect2)
mean1 = np.mean(mat[17:20,0:13])
mean2 = np.mean(mat[21:24,28:41])
# add text
ax.text(6.5, 15.5,f'{mean1:.2f}',ha='center',va='center',fontsize=10,color='black')
ax.text(34.5, 26.5,f'{mean2:.2f}',ha='center',va='center',fontsize=10,color='black')
cax = fig.add_axes([0.25, -0.05, 0.4, 0.03]) #left,bottom,width,height
cbar = plt.colorbar(
    heatmap.collections[0],
    cax=cax,
    orientation='horizontal'
)
cbar.outline.set_linewidth(0.8)
cbar.ax.tick_params(labelsize=9,width=0.8,length=2)
cbar.set_ticks([0.5,1.0,1.5])
cbar.ax.text(1.05, 0.5,'Obs/Exp',transform=cbar.ax.transAxes,ha='left',va='center',fontsize=9)
ax.set_xticks([0,10,20,30,40])
ax.set_xticklabels(['-2Mb','-1Mb','0','1Mb','2Mb'], fontsize=8)
ax.set_yticks([0,10,20,30,40])
ax.set_yticklabels(['2Mb','1Mb','0','-1Mb','-2Mb'], fontsize=8)
ax.set_title('WT', fontsize=12)
fig.show()
fig.savefig('/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/Domains/Cooltools/ATA/Downsample/Pdf/WT_coolupBoundariesPileup.heatmap.pdf', 
                bbox_inches = 'tight',
                facecolor='w')


f = h5py.File(
    "/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/Domains/Cooltools/ATA/Downsample/KO.boundary.downsample.pileup.npz",
    "r"
)
print(list(f.keys()))
mat = f['data'][:]
style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
plt.rcParams['font.sans-serif']='Helvetica'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['axes.unicode_minus'] =False
fig = plt.figure(figsize=(4/1.2, 3/1.2))
gs = gridspec.GridSpec(3, 2,  width_ratios=[4,0.2],height_ratios=[1,1,2],wspace=0.1)
colors1 = [(0, '#1C4D8D'),
           (0.5, '#FAF3E1'),   
          (1, '#C40C0C')] 
ccmap=LinearSegmentedColormap.from_list('custom_cmap', colors1)
ax = fig.add_subplot(gs[:,0])
heatmap=sns.heatmap(mat,cmap=ccmap,rasterized=True,
                        vmin=0.5,vmax=1.5,
                        cbar=False,
                        xticklabels=False,
                        yticklabels=False,square=True)
for spine in heatmap.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1)
        spine.set_color('black')
rect1 = Rectangle((0,17),13,3,
                  fill=False,
                  edgecolor='black',
                  linewidth=1,linestyle='--')
ax.add_patch(rect1)
rect2 = Rectangle((28,21),13,3,
                  fill=False,
                  edgecolor='black',
                  linewidth=1,linestyle='--')
ax.add_patch(rect2)
mean1 = np.mean(mat[17:20,0:13])
mean2 = np.mean(mat[21:24,28:41])
# add text
ax.text(6.5, 15.5,f'{mean1:.2f}',ha='center',va='center',fontsize=10,color='black')
ax.text(34.5, 26.5,f'{mean2:.2f}',ha='center',va='center',fontsize=10,color='black')
cax = fig.add_axes([0.25, -0.05, 0.4, 0.03]) #left,bottom,width,height
cbar = plt.colorbar(
    heatmap.collections[0],
    cax=cax,
    orientation='horizontal'
)
cbar.outline.set_linewidth(0.8)
cbar.ax.tick_params(labelsize=9,width=0.8,length=2)
cbar.set_ticks([0.5,1.0,1.5])
cbar.ax.text(1.05, 0.5,'Obs/Exp',transform=cbar.ax.transAxes,ha='left',va='center',fontsize=9)
ax.set_xticks([0,10,20,30,40])
ax.set_xticklabels(['-2Mb','-1Mb','0','1Mb','2Mb'], fontsize=8)
ax.set_yticks([0,10,20,30,40])
ax.set_yticklabels(['2Mb','1Mb','0','-1Mb','-2Mb'], fontsize=8)
ax.set_title('DEAF1-mutant', fontsize=12)
fig.show()
fig.savefig('/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/Domains/Cooltools/ATA/Downsample/Pdf/KO_coolupBoundariesPileup.heatmap.pdf', 
                bbox_inches = 'tight',
                facecolor='w')
