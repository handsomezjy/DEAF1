import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt
from matplotlib import style
from matplotlib.patches import Patch
from scipy import stats
import pandas as pd

ctcfdeaf1=pd.read_table("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/Domains/Cooltools/Bedtools/CTCFDEAF1.TAD_boundaries_merged.10kb.downsample.PeakSignal.txt")
deaf1=pd.read_table("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/Domains/Cooltools/Bedtools/DEAF1alone.TAD_boundaries_merged.10kb.downsample.PeakSignal.txt")
ctcf=pd.read_table("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/Domains/Cooltools/Bedtools/CTCFalone.TAD_boundaries_merged.10kb.downsample.PeakSignal.txt")

style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
fig, ax = plt.subplots(figsize=(2.5, 3))
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['axes.unicode_minus'] =False
plt.rcParams['font.sans-serif']='Helvetica'
colors=['#7F8487','#8E1616']*2
data=[
    ctcfdeaf1["'E250604004_L01_WT_CTCF'"].to_list(),
    ctcfdeaf1["'E250609001_L01_DEAF1_CTCF'"].to_list(),
    ctcf["'E250604004_L01_WT_CTCF'"].to_list(),
    ctcf["'E250609001_L01_DEAF1_CTCF'"].to_list(),
]
positions = [1.75, 2.25, 3.25, 3.75]
box = ax.boxplot(data,positions=positions, patch_artist=True, widths=0.3,
                showmeans=False, showfliers=False) 
plt.ylim(0,2.5)
plt.ylabel('ChIP-seq Signal')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_color('grey')
ax.spines['bottom'].set_color('grey')
ax.spines['left'].set_position(('outward', 5))
ax.spines['bottom'].set_position(('outward', 5))
plt.xticks([1.9,3.6],[r'$DEAF1^{+} CTCF^{+}$',r'$DEAF1^{-} CTCF^{+}$'],
           rotation=0)
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.9)
for element in ['whiskers', 'caps', 'medians']:
    plt.setp(box[element], color='k', linewidth=0.8)
plt.title('CTCF at TAD boundaries')
legend_elements = [
    Patch(facecolor='#44444E', edgecolor='none', label='WT'),
    Patch(facecolor='#740A03', edgecolor='none', label='DEAF1-mutant')
]
ax.legend(handles=legend_elements,frameon=False,fontsize=9,loc='upper right',ncol=2,bbox_to_anchor=(1.1, 1.05),handletextpad=0.3)

h_stat, p_kruskal = stats.kruskal(data[0],data[1])
print(f"Kruskal-Wallis p-value: {p_kruskal:.4f}")
ax.plot([1.5,2.5], [2.1,2.1], color='k', lw=1)
ax.text(2, 2.1, 'n.s.' , 
        ha='center', va='bottom', color='k', fontsize=10)
h_stat, p_kruskal = stats.kruskal(data[2],data[3])
print(f"Kruskal-Wallis p-value: {p_kruskal}")
ax.plot([3,4], [2.1,2.1], color='k', lw=1)
ax.text(3.5, 2.1, f"P={p_kruskal:.2e}" , 
        ha='center', va='bottom', color='k', fontsize=10)
plt.savefig('/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/Domains/Cooltools/Bedtools/Pdf/CTCFsignal.TADboundaries.boxplot.pdf',   ##
                bbox_inches = 'tight',
                facecolor='w')  

