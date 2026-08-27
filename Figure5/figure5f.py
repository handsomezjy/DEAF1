import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt
from matplotlib import style
from matplotlib.patches import Patch
from scipy import stats
import pandas as pd

ctcfdeaf1=pd.read_table("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/Domains/Cooltools/Bedtools/CTCFDEAF1.TAD_boundaries_merged.10kb.downsample.InsulationScores.txt")
deaf1=pd.read_table("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/Domains/Cooltools/Bedtools/DEAF1alone.TAD_boundaries_merged.10kb.downsample.InsulationScores.txt")
ctcf=pd.read_table("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/Domains/Cooltools/Bedtools/CTCFalone.TAD_boundaries_merged.10kb.downsample.InsulationScores.txt")
ctcfdeaf1=ctcfdeaf1.dropna()
deaf1=deaf1.dropna()
ctcf=ctcf.dropna()

style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
fig, ax = plt.subplots(figsize=(2.5, 3))
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['axes.unicode_minus'] =False
plt.rcParams['font.sans-serif']='Helvetica'
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Helvetica'
plt.rcParams['mathtext.it'] = 'Helvetica:italic'
plt.rcParams['mathtext.bf'] = 'Helvetica:bold'
colors=['#7F8487','#8E1616']*3
data=[
    ctcfdeaf1["'WT.insulation_10kb'"].to_list(),
    ctcfdeaf1["'KO.insulation_10kb.downsample'"].to_list(),
    deaf1["'WT.insulation_10kb'"].to_list(),
    deaf1["'KO.insulation_10kb.downsample'"].to_list(),
    ctcf["'WT.insulation_10kb'"].to_list(),
    ctcf["'KO.insulation_10kb.downsample'"].to_list()
]
positions = [0.82,1.15,1.85, 2.15, 2.85, 3.15]
box = ax.boxplot(data,positions=positions, patch_artist=True, widths=0.2,
                showmeans=False, showfliers=False) 
plt.ylim(-2,0.5)
plt.yticks((-2,-1,0))
plt.ylabel('Insulation scores')
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
plt.xticks((1,2,3),('$DEAF1^{+}$ $CTCF^{+}$','$DEAF1^{+}$ $CTCF^{-}$','$DEAF1^{-}$ $CTCF^{+}$'),rotation=40)
h_stat, p_kruskal = stats.kruskal(data[0],data[1])
print(f"Kruskal-Wallis p-value: {p_kruskal:.4f}")
ax.plot([0.7,1.3], [0.1,0.1], color='k', lw=1)
ax.text(1,0.1, f"P={p_kruskal:.4f}", 
        ha='center', va='bottom', color='k', fontsize=9)
h_stat, p_kruskal = stats.kruskal(data[2],data[3])
print(f"Kruskal-Wallis p-value: {p_kruskal:.4f}")
ax.plot([1.7,2.3], [0.3,0.3], color='k', lw=1)
ax.text(2,0.3, f"P={p_kruskal:.4f}", ha='center', va='bottom', color='k', fontsize=9)
h_stat, p_kruskal = stats.kruskal(data[4],data[5])
print(f"Kruskal-Wallis p-value: {p_kruskal:.4e}")
ax.plot([2.7,3.3], [.4,.4], color='k', lw=1)
ax.text(3,.4, f"P={p_kruskal:.2e}", ha='center', va='bottom', color='k', fontsize=9)
plt.title('TAD-boundary peaks',y=1.05)
legend_elements = [
    Patch(facecolor=colors[0], edgecolor='k', label='WT'),
    Patch(facecolor=colors[1], edgecolor='k', label='DEAF1-mutant')
]
plt.legend(loc='upper right',bbox_to_anchor=(1.35,0.2),handles=legend_elements, frameon=False,ncol=1)
plt.savefig('/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/Domains/Cooltools/Bedtools/Pdf/CTCFDEAF1.Insulationscores.TADboundaries.boxplot.pdf',   ##
                bbox_inches = 'tight',
                facecolor='w')  
