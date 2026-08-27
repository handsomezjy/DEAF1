import matplotlib.font_manager as font_manager
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import scipy.stats as stats

fimocd=pd.read_table("/data/zhangjy/DEAF1/ChiPseq_Analysis/ResultsSort/Peaksignal/merged_CTCF_intersectDEAF1/fimo.tsv")
fimoc=pd.read_table("/data/zhangjy/DEAF1/ChiPseq_Analysis/ResultsSort/Peaksignal/merged_CTCF_noDEAF1/fimo.tsv")
fimod=pd.read_table("/data/zhangjy/DEAF1/ChiPseq_Analysis/ResultsSort/Peaksignal/merged_DEAF1_noCTCF/fimo.tsv")
# cba=pd.read_table("/data/zhangjy/DEAF1/ChiPseq_Analysis/ResultsSort/Deeptools/InputNew/BedFiles/CTCF_before_alone.bed",header=None)
li1=[2647,2534,66446]
li2=[len(fimocd['sequence_name'].unique()),len(fimod['sequence_name'].unique()),len(fimoc['sequence_name'].unique())]
li3=[x-y for x, y in zip(li1, li2)]
li4=[-np.log2(fimocd.groupby('sequence_name')['p-value'].min()),
     -np.log2(fimod.groupby('sequence_name')['p-value'].min()),
     -np.log2(fimoc.groupby('sequence_name')['p-value'].min())]


style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
plt.rcParams['font.sans-serif']='Helvetica'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['axes.unicode_minus'] =False
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Helvetica'
plt.rcParams['mathtext.it'] = 'Helvetica:italic'
plt.rcParams['mathtext.bf'] = 'Helvetica:bold'
fig = plt.figure(figsize=(2.2,3))
gs = gridspec.GridSpec(2,1 ,height_ratios=[1,1],hspace=0.1)
colors=['#444444','#FFD369']
ax1=fig.add_subplot(gs[0])
for i in range(3):   
    plt.bar(i,li3[i]/li1[i],color=colors[0],width=0.5)
    plt.bar(i,li2[i]/li1[i],color=colors[1],width=0.5,bottom=li3[i]/li1[i])
    plt.text(i-0.25,li3[i]/li1[i]/2,s=li3[i])
    plt.text(i-0.25,li2[i]/li1[i]/2+li3[i]/li1[i],s=li2[i])
plt.ylim(0,1)
# plt.xlim(0.5,3.5)
plt.yticks((0,0.5,1),(0,'50%',1))
legend_elements = [
    Patch(facecolor=colors[0], edgecolor='k', label='CTCF motif-'),
    Patch(facecolor=colors[1], edgecolor='k', label='CTCF motif+')
]
plt.legend(loc='upper center',bbox_to_anchor=(0.5,1.3),handles=legend_elements, frameon=False,ncol=2,columnspacing=0.8,handletextpad=0.2)

ax2=fig.add_subplot(gs[1])
colors=['#BBBBBB']*3
for ax in [ax1,ax2]:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('grey')
    ax.spines['bottom'].set_color('grey')
    ax.spines['left'].set_position(('outward', 5))
    ax.spines['bottom'].set_position(('outward', 5))
positions = list(range(3))
box = ax.boxplot(li4,positions=positions, patch_artist=True, widths=0.3,
                showmeans=False, showfliers=False) 
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.9)
for element in ['whiskers', 'caps', 'medians']:
    plt.setp(box[element], color='k', linewidth=0.8)
plt.xticks((0,1,2),
           ('$DEAF1^{+}$\n$CTCF^{+}$','$DEAF1^{+}$\n$CTCF^{-}$','$DEAF1^{-}$\n$CTCF^{+}$'),rotation=0)
plt.ylim(12,38)
plt.yticks((15,20,25))
plt.ylabel('-log2(p-value)\nCTCF motif scores')
h_stat, p_kruskal = stats.kruskal(li4[0],li4[1])
print(f"Kruskal-Wallis p-value: {p_kruskal:.2e}")
ax2.plot([-0.1,1.1], [24,24], color='k', lw=1)
ax2.text(0.5,24, f"P={p_kruskal:.2e}", 
        ha='center', va='bottom', color='k', fontsize=10)
h_stat, p_kruskal = stats.kruskal(li4[0],li4[2])
print(f"Kruskal-Wallis p-value: {p_kruskal:.2e}")
ax2.plot([-0.1,2.1], [34,34], color='k', lw=1)
ax2.text(1,34, f"P={p_kruskal:.2e}", 
        ha='center', va='bottom', color='k', fontsize=10)
h_stat, p_kruskal = stats.kruskal(li4[1],li4[2])
print(f"Kruskal-Wallis p-value: {p_kruskal:.2e}")
ax2.plot([0.9,2.1], [29,29], color='k', lw=1)
ax2.text(1.5,29, f"P={p_kruskal:.2e}", 
        ha='center', va='bottom', color='k', fontsize=10)
plt.savefig('/data/zhangjy/DEAF1/ChiPseq_Analysis/ResultsSort/Pdf/CTCFDEAF1.ctcfmotifscores.pdf',   ##
                bbox_inches = 'tight',
                facecolor='w')  
