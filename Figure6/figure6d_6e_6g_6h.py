import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import style
from scipy import stats

up_ep_k4=pd.read_table("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/MultiOmics_Analysis/matrix_files/up_EP__H3K4me3.mat.gz",skiprows=1,header=None)
down_ep_k4=pd.read_table("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/MultiOmics_Analysis/matrix_files/down_EP__H3K4me3.mat.gz",skiprows=1,header=None)
up_ep_k27=pd.read_table("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/MultiOmics_Analysis/matrix_files/up_EP__H3K27ac_both_anchors.mat.gz",skiprows=1,header=None)
down_ep_k27=pd.read_table("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/MultiOmics_Analysis/matrix_files/down_EP__H3K27ac_both_anchors.mat.gz",skiprows=1,header=None)

style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
plt.rcParams['font.sans-serif']='Helvetica'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['axes.unicode_minus'] =False
fig, ax = plt.subplots(figsize=(1.8, 1.5))

plt.plot(list(range(1,401)),down_ep_k4.iloc[:,6:406].mean(),color='#7F8487', linewidth = "1.5",label='WT')
plt.plot(list(range(1,401)),down_ep_k4.iloc[:,406:].mean(),color='#8E1616', linewidth = "1.5",label='DEAF1-mutant')
plt.ylim(0,30)
plt.yticks([0,10,20,30],fontsize=10)
plt.xticks([1,200.5,400],['-2.0kb','down-tss\n(within anchor)','2.0kb'],rotation=0)
plt.grid(axis='both', color='grey', linestyle=':',zorder=0)
plt.title('H3K4me3')
plt.legend(frameon=False,bbox_to_anchor=(0.5, 1.45), 
           loc='upper center',ncol=2)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_color('grey')
ax.spines['bottom'].set_color('grey')
ax.spines['left'].set_position(('outward', 5))
ax.spines['bottom'].set_position(('outward', 5))
plt.savefig("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/Pdf/H3K4me3.downEPtss2kb.lineplot.pdf",
                    bbox_inches = 'tight',
                    facecolor='w')

style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
plt.rcParams['font.sans-serif']='Helvetica'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['axes.unicode_minus'] =False
fig, ax = plt.subplots(figsize=(1.8, 1.5))

plt.plot(list(range(1,401)),up_ep_k4.iloc[:,6:406].mean(),color='#7F8487', linewidth = "1.5",label='WT')
plt.plot(list(range(1,401)),up_ep_k4.iloc[:,406:].mean(),color='#8E1616', linewidth = "1.5",label='DEAF1-mutant')
plt.ylim(0,32)
plt.yticks([0,10,20,30],fontsize=10)
plt.xticks([1,200.5,400],['-2.0kb','up-tss\n(within anchor)','2.0kb'],rotation=0)
plt.grid(axis='both', color='grey', linestyle=':',zorder=0)
plt.title('H3K4me3')
plt.legend(frameon=False,bbox_to_anchor=(0.5, 1.45), 
           loc='upper center',ncol=2)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_color('grey')
ax.spines['bottom'].set_color('grey')
ax.spines['left'].set_position(('outward', 5))
ax.spines['bottom'].set_position(('outward', 5))
plt.savefig("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/Pdf/H3K4me3.upEPtss2kb.lineplot.pdf",
                    bbox_inches = 'tight',
                    facecolor='w')


style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
fig, ax = plt.subplots(figsize=(1.25, 1.8))
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['axes.unicode_minus'] =False
plt.rcParams['font.sans-serif']='Helvetica'
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Helvetica'
plt.rcParams['mathtext.it'] = 'Helvetica:italic'
plt.rcParams['mathtext.bf'] = 'Helvetica:bold'
colors=['#7F8487','#8E1616']
data=[
    down_ep_k4.iloc[:,6+195:406-195].mean(axis=1),
    down_ep_k4.iloc[:,406+195:806-195].mean(axis=1)
]
positions = [0.9, 2.1]
box = ax.boxplot(data,positions=positions, patch_artist=True, widths=0.3,
                showmeans=False, showfliers=False) 
plt.ylim(-5,85)
plt.yticks((0,40,80))
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_color('grey')
ax.spines['bottom'].set_color('grey')
ax.spines['left'].set_position(('outward', 5))
ax.spines['bottom'].set_position(('outward', 5))
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.9)
for element in ['whiskers', 'caps', 'medians']:
    plt.setp(box[element], color='k', linewidth=0.8)
h_stat, p_kruskal = stats.kruskal(data[0],data[1])
ax.plot([0.7,2.3], [83,83], color='k', lw=1)
ax.text(1.5,83, f"P={p_kruskal:.3e}", 
        ha='center', va='bottom', color='k', fontsize=8)
plt.xticks((0.9,2.1),('WT','DEAF1-mutant'),fontsize=8)
plt.title('H3K4me3',y=1.05)
plt.xlabel('Down-tss(within anchor)')
plt.savefig("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/Pdf/H3K4me3.downEPtss2kb.boxplot.pdf",
                    bbox_inches = 'tight',
                    facecolor='w')

style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
fig, ax = plt.subplots(figsize=(1.2, 1.8))
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['axes.unicode_minus'] =False
plt.rcParams['font.sans-serif']='Helvetica'
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Helvetica'
plt.rcParams['mathtext.it'] = 'Helvetica:italic'
plt.rcParams['mathtext.bf'] = 'Helvetica:bold'
colors=['#7F8487','#8E1616']
data=[
    up_ep_k4.iloc[:,6+195:406-195].mean(axis=1),
    up_ep_k4.iloc[:,406+195:806-195].mean(axis=1)
]
positions = [0.9, 2.1]
box = ax.boxplot(data,positions=positions, patch_artist=True, widths=0.3,
                showmeans=False, showfliers=False) 
plt.ylim(0,95)
plt.yticks((0,40,80))
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_color('grey')
ax.spines['bottom'].set_color('grey')
ax.spines['left'].set_position(('outward', 5))
ax.spines['bottom'].set_position(('outward', 5))
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.9)
for element in ['whiskers', 'caps', 'medians']:
    plt.setp(box[element], color='k', linewidth=0.8)
h_stat, p_kruskal = stats.kruskal(data[0],data[1])
ax.plot([0.7,2.3], [94,94], color='k', lw=1)
ax.text(1.5,94, f"P={p_kruskal:.4f}", 
        ha='center', va='bottom', color='k', fontsize=8)
plt.xticks((0.9,2.1),('WT','DEAF1-mutant'),fontsize=8)
plt.title('H3K4me3',y=1.05)
plt.xlabel('Up-tss(within anchor)')
plt.savefig("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/Pdf/H3K4me3.upEPtss2kb.boxplot.pdf",
                    bbox_inches = 'tight',
                    facecolor='w')

style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
plt.rcParams['font.sans-serif']='Helvetica'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['axes.unicode_minus'] =False
fig, ax = plt.subplots(figsize=(1.8, 1.5))

plt.plot(list(range(1,401)),up_ep_k27.iloc[142:,6:406].mean(),color='#7F8487', linewidth = "1.5",label='WT')
plt.plot(list(range(1,401)),up_ep_k27.iloc[142:,406:].mean(),color='#8E1616', linewidth = "1.5",label='DEAF1-mutant')
plt.ylim(0,15)
plt.yticks([0,10],fontsize=10)
plt.xticks([1,200.5,400],['-2.0kb','distal E\n(within anchor)','2.0kb'],rotation=0)
plt.grid(axis='both', color='grey', linestyle=':',zorder=0)
plt.title('H3K27ac')
plt.legend(frameon=False,bbox_to_anchor=(0.5, 1.45), 
           loc='upper center',ncol=2)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_color('grey')
ax.spines['bottom'].set_color('grey')
ax.spines['left'].set_position(('outward', 5))
ax.spines['bottom'].set_position(('outward', 5))
plt.savefig("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/Pdf/H3K27ac.upEPtss2kb.lineplot.pdf",
                    bbox_inches = 'tight',
                    facecolor='w')


style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
plt.rcParams['font.sans-serif']='Helvetica'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['axes.unicode_minus'] =False
fig, ax = plt.subplots(figsize=(1.8, 1.5))

plt.plot(list(range(1,401)),down_ep_k27.iloc[222:,6:406].mean(),color='#7F8487', linewidth = "1.5",label='WT')
plt.plot(list(range(1,401)),down_ep_k27.iloc[222:,406:].mean(),color='#8E1616', linewidth = "1.5",label='DEAF1-mutant')
plt.ylim(0,13)
plt.yticks([0,10],fontsize=10)
plt.xticks([1,200.5,400],['-2.0kb','distal E\n(within anchor)','2.0kb'],rotation=0)
plt.grid(axis='both', color='grey', linestyle=':',zorder=0)
plt.title('H3K27ac')
plt.legend(frameon=False,bbox_to_anchor=(0.5, 1.45), 
           loc='upper center',ncol=2)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_color('grey')
ax.spines['bottom'].set_color('grey')
ax.spines['left'].set_position(('outward', 5))
ax.spines['bottom'].set_position(('outward', 5))
plt.savefig("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/Pdf/H3K27ac.downEPtss2kb.lineplot.pdf",
                    bbox_inches = 'tight',
                    facecolor='w')


style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
fig, ax = plt.subplots(figsize=(1.2, 1.8))
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['axes.unicode_minus'] =False
plt.rcParams['font.sans-serif']='Helvetica'
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Helvetica'
plt.rcParams['mathtext.it'] = 'Helvetica:italic'
plt.rcParams['mathtext.bf'] = 'Helvetica:bold'
colors=['#7F8487','#8E1616']
data=[
    down_ep_k27.iloc[222:,6+195:406-195].mean(axis=1),
    down_ep_k27.iloc[222:,406+195:806-195].mean(axis=1)
]
positions = [0.9, 2.1]
box = ax.boxplot(data,positions=positions, patch_artist=True, widths=0.3,
                showmeans=False, showfliers=False) 

plt.ylim(0,30)
plt.yticks((0,10,20,30))
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_color('grey')
ax.spines['bottom'].set_color('grey')
ax.spines['left'].set_position(('outward', 5))
ax.spines['bottom'].set_position(('outward', 5))
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.9)
for element in ['whiskers', 'caps', 'medians']:
    plt.setp(box[element], color='k', linewidth=0.8)
h_stat, p_kruskal = stats.kruskal(data[0],data[1])
ax.plot([0.7,2.3], [30,30], color='k', lw=1)
ax.text(1.5,30, f"P={p_kruskal:.4f}", 
        ha='center', va='bottom', color='k', fontsize=8)
plt.xticks((0.9,2.1),('WT','DEAF1-mutant'),fontsize=8)
plt.title('H3K27ac',y=1.05)
plt.xlabel('distal-E(within anchor)')
plt.savefig("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/Pdf/H3K27ac.downEPtss2kb.boxplot.pdf",
                    bbox_inches = 'tight',
                    facecolor='w')


style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
fig, ax = plt.subplots(figsize=(1.2, 1.8))
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['axes.unicode_minus'] =False
plt.rcParams['font.sans-serif']='Helvetica'
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Helvetica'
plt.rcParams['mathtext.it'] = 'Helvetica:italic'
plt.rcParams['mathtext.bf'] = 'Helvetica:bold'
colors=['#7F8487','#8E1616']
data=[
    up_ep_k27.iloc[142:,6+195:406-195].mean(axis=1),
    up_ep_k27.iloc[142:,406+195:806-195].mean(axis=1)
]
positions = [0.9, 2.1]
box = ax.boxplot(data,positions=positions, patch_artist=True, widths=0.3,
                showmeans=False, showfliers=False) 
plt.ylim(0,36)
plt.yticks((0,10,20,30))
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_color('grey')
ax.spines['bottom'].set_color('grey')
ax.spines['left'].set_position(('outward', 5))
ax.spines['bottom'].set_position(('outward', 5))
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.9)
for element in ['whiskers', 'caps', 'medians']:
    plt.setp(box[element], color='k', linewidth=0.8)
h_stat, p_kruskal = stats.kruskal(data[0],data[1])
ax.plot([0.7,2.3], [35,35], color='k', lw=1)
ax.text(1.5,35, f"P={p_kruskal:.3e}", 
        ha='center', va='bottom', color='k', fontsize=8)
plt.xticks((0.9,2.1),('WT','DEAF1-mutant'),fontsize=8)
plt.title('H3K27ac',y=1.05)
plt.xlabel('distal-E(within anchor)')
plt.savefig("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/Pdf/H3K27ac.upEPtss2kb.boxplot.pdf",
                    bbox_inches = 'tight',
                    facecolor='w')



