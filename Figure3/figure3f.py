import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

ctcfdeaf1=pd.read_table("/data/zhangjy/DEAF1/ChiPseq_Analysis/ResultsSort/Peaksignal/merged_CTCF_intersectDEAF1.PeakSignal.txt")
deaf1=pd.read_table("/data/zhangjy/DEAF1/ChiPseq_Analysis/ResultsSort/Peaksignal/merged_DEAF1_noCTCF.PeakSignal.txt")
ctcf=pd.read_table("/data/zhangjy/DEAF1/ChiPseq_Analysis/ResultsSort/Peaksignal/merged_CTCF_noDEAF1.PeakSignal.txt")
deaf1=deaf1[~(deaf1.iloc[:, -4:] > 100).any(axis=1)]
ctcf=ctcf[~(ctcf.iloc[:, -4:] > 100).any(axis=1)]
ctcfdeaf1=ctcfdeaf1[~(ctcfdeaf1.iloc[:, -4:] > 100).any(axis=1)]

style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
fig, ax = plt.subplots(figsize=(2., 2.7))
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['axes.unicode_minus'] =False
plt.rcParams['font.sans-serif']='Helvetica'
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Helvetica'
plt.rcParams['mathtext.it'] = 'Helvetica:italic'
plt.rcParams['mathtext.bf'] = 'Helvetica:bold'
colors=['#7F8487','#8E1616']*2
data=[
    np.log2(ctcfdeaf1["'KO_deaf1'"]/ctcfdeaf1["'WT_deaf1'"]).to_list(),
    np.log2(deaf1["'KO_deaf1'"]/deaf1["'WT_deaf1'"]).to_list(),
    np.log2(ctcfdeaf1["'KO_ctcf'"]/ctcfdeaf1["'WT_ctcf'"]).to_list(),
    np.log2(ctcf["'KO_ctcf'"]/ctcf["'WT_ctcf'"]).to_list()
]
positions = [0.82,1.15,1.85, 2.15]
box = ax.boxplot(data,positions=positions, patch_artist=True, widths=0.2,
                showmeans=False, showfliers=False) 
plt.ylim(-1.8,2.1)
plt.xlim(0.5,2.5)
plt.yticks((-1,0,1,2))
plt.ylabel('log2FC(ChIP-seq signal)')
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
plt.xticks((0.82,1.15,1.85, 2.15),
           ('$DEAF1^{+}$ $CTCF^{+}$','$DEAF1^{+}$ $CTCF^{-}$','$DEAF1^{+}$ $CTCF^{+}$','$DEAF1^{-}$ $CTCF^{+}$'),rotation=60)
h_stat, p_kruskal = stats.kruskal(data[0],data[1])
print(f"Kruskal-Wallis p-value: {p_kruskal:.4e}")
ax.plot([0.7,1.3], [1.8,1.8], color='k', lw=1)
ax.text(1,1.8, f"DEAF1\nP={p_kruskal:.2e}", ha='center', va='bottom', color='k', fontsize=9)

h_stat, p_kruskal = stats.kruskal(data[1],data[2])
print(f"Kruskal-Wallis p-value: {p_kruskal:.4e}")
ax.plot([1.7,2.3], [2.1,2.1], color='k', lw=1)
ax.text(2,2.1, f"CTCF\nP={p_kruskal:.2e}", ha='center', va='bottom', color='k', fontsize=9)
plt.savefig('/data/zhangjy/DEAF1/ChiPseq_Analysis/ResultsSort/Deeptools/InputNew/Deeptools/mergedPeaks/Pdf/CTCFDEAF1.chipseqSignalRPGC.boxpplot.pdf',   ##
                bbox_inches = 'tight',
                facecolor='w')  
