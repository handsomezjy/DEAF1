import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import style
from scipy import stats

w5=pd.read_table("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/LoopStrength/WT.5kb.rescored_strength_coords.txt")
w10=pd.read_table("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/LoopStrength/WT.10kb.rescored_strength_coords.txt")
k5=pd.read_table("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/LoopStrength/KO.5kb.rescored_strength_coords.txt")
k10=pd.read_table("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/LoopStrength/KO.10kb.rescored_strength_coords.txt")

ws=pd.concat([w5,w10]).sort_values(['chrom1','start1'])
ks=pd.concat([k5,k10]).sort_values(['chrom1','start1'])
ws_tmp=ws.iloc[:,[0,1,2,3,4,5,10]]
ws_tmp.columns=[0,1,2,3,4,5,'strengthWT']
ks_tmp=ks.iloc[:,[0,1,2,3,4,5,10]]
ks_tmp.columns=[0,1,2,3,4,5,'strengthKO']
down_ep=pd.read_table("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/APA_results/Difftypes/tmp_split_bedpe/down_anchor_EP_.bedpe",header=None)
down_pp=pd.read_table("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/APA_results/Difftypes/tmp_split_bedpe/down_anchor_PP_.bedpe",header=None)
down_none=pd.read_table("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/APA_results/Difftypes/tmp_split_bedpe/down_anchor_none_.bedpe",header=None)

down_ep_ls=pd.merge(down_ep,ws_tmp,on=[0,1,2,3,4,5])
down_ep_ls=pd.merge(down_ep_ls,ks_tmp,on=[0,1,2,3,4,5])

down_pp_ls=pd.merge(down_pp,ws_tmp,on=[0,1,2,3,4,5])
down_pp_ls=pd.merge(down_pp_ls,ks_tmp,on=[0,1,2,3,4,5])

down_none_ls=pd.merge(down_none,ws_tmp,on=[0,1,2,3,4,5])
down_none_ls=pd.merge(down_none_ls,ks_tmp,on=[0,1,2,3,4,5])


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
    down_ep_ls['strengthWT'],
    down_ep_ls['strengthKO']
]
positions = [0.9, 2.1]
box = ax.boxplot(data,positions=positions, patch_artist=True, widths=0.3,
                showmeans=False, showfliers=False) 

plt.ylim(-0.1,8)
plt.yticks((0,2,4,6))
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
ax.plot([0.7,2.3], [8,8], color='k', lw=1)
ax.text(1.5,8, f"P={p_kruskal:.4f}", 
        ha='center', va='bottom', color='k', fontsize=8)
plt.xticks((0.9,2.1),('WT','DEAF1-mutant'),fontsize=8)
plt.title('Loop strength',y=1.05)
plt.xlabel('E-P(down-tss)')
plt.savefig("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/Pdf/loopstrength.downEP.boxplot.pdf",
                    bbox_inches = 'tight',
                    facecolor='w')

up_ep=pd.read_table("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/APA_results/Difftypes/tmp_split_bedpe/up_anchor_EP_.bedpe",header=None)
up_pp=pd.read_table("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/APA_results/Difftypes/tmp_split_bedpe/up_anchor_PP_.bedpe",header=None)
up_none=pd.read_table("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/APA_results/Difftypes/tmp_split_bedpe/up_anchor_none_.bedpe",header=None)

up_ep_ls=pd.merge(up_ep,ws_tmp,on=[0,1,2,3,4,5])
up_ep_ls=pd.merge(up_ep_ls,ks_tmp,on=[0,1,2,3,4,5])

up_pp_ls=pd.merge(up_pp,ws_tmp,on=[0,1,2,3,4,5])
up_pp_ls=pd.merge(up_pp_ls,ks_tmp,on=[0,1,2,3,4,5])

up_none_ls=pd.merge(up_none,ws_tmp,on=[0,1,2,3,4,5])
up_none_ls=pd.merge(up_none_ls,ks_tmp,on=[0,1,2,3,4,5])


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
    up_ep_ls['strengthWT'],
    up_ep_ls['strengthKO']
]
positions = [0.9, 2.1]
box = ax.boxplot(data,positions=positions, patch_artist=True, widths=0.3,
                showmeans=False, showfliers=False) 
plt.ylim(0,8)
plt.yticks((0,2,4,6))
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
ax.plot([0.7,2.3], [8,8], color='k', lw=1)
ax.text(1.5,8, "n.s.", 
        ha='center', va='bottom', color='k', fontsize=8)
plt.xticks((0.9,2.1),('WT','DEAF1-mutant'),fontsize=8)
plt.title('Loop strength',y=1.05)
plt.xlabel('E-P(up-tss)')
plt.savefig("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/Pdf/loopstrength.upEP.boxplot.pdf",
                    bbox_inches = 'tight',
                    facecolor='w')
