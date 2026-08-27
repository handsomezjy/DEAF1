import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import scikit_posthocs as sp

ctcf=np.load("/data/zhangjy/DEAF1/Histone_Analysis/Datasort/Bedtools/CTCFbeforeAlone.histones.npz")
deaf1=np.load("/data/zhangjy/DEAF1/Histone_Analysis/Datasort/Bedtools/DEAF1beforeAlone.histones.npz")
cd=np.load("/data/zhangjy/DEAF1/Histone_Analysis/Datasort/Bedtools/CTCFbeforewithDEAF1.histones.npz")

def remove_outliers_iqr(data):
    if len(data) == 0:
        return data
    q1 = np.percentile(data, 25)
    q3 = np.percentile(data, 75)
    iqr = q3 - q1
    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr
    return data[(data >= lower_bound) & (data <= upper_bound)]
for i in range(5):   
    style.use('default')
    font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
    plt.rcParams['font.sans-serif']='Helvetica'
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['axes.unicode_minus'] = False
    plt.rcParams['mathtext.rm'] = 'Helvetica'
    plt.rcParams['mathtext.it'] = 'Helvetica:italic'
    plt.rcParams['mathtext.bf'] = 'Helvetica:bold'

    fig, ax = plt.subplots(figsize=(2, 2.6))
    g1_raw = cd['matrix'][:, i]
    g2_raw = ctcf['matrix'][:, i]
    g3_raw = deaf1['matrix'][:, i]
    
    g1_raw = g1_raw[~np.isnan(g1_raw)]
    g2_raw = g2_raw[~np.isnan(g2_raw)]
    g3_raw = g3_raw[~np.isnan(g3_raw)]

    group1 = remove_outliers_iqr(g1_raw)
    group2 = remove_outliers_iqr(g2_raw)
    group3 = remove_outliers_iqr(g3_raw)
    box = ax.boxplot(
        [group1, group2, group3],
        patch_artist=True,     # 允许填充颜色
        showfliers=False,      # 不显示离群点
        widths=0.3,
        medianprops=dict(color='black', linewidth=1),
        boxprops=dict(linewidth=1),
        whiskerprops=dict(linewidth=1),
        capprops=dict(linewidth=1)
    )
    labels = [
    'DEAF1$^{+}$ CTCF$^{+}$',
    'DEAF1$^{-}$ CTCF$^{+}$',
    'DEAF1$^{+}$ CTCF$^{-}$'
    ]
    colors = ['#FFDBB5','#825B32','#6CBEC7']
    for patch, color in zip(box['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.9)
    for element in ['whiskers', 'caps', 'medians']:
        plt.setp(box[element], color='k', linewidth=0.8)
    
    ax.set_ylabel('Normalized ChIP-seq\nsignal(RPGC)')
    ax.set_title(ctcf['labels'][i].split('.')[0].split('-')[1])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('grey')
    ax.spines['bottom'].set_color('grey')
    ax.spines['left'].set_position(('outward', 5))
    ax.spines['bottom'].set_position(('outward', 5))
    ax.set_xticklabels(labels,rotation=35,fontsize=8)

    data_list = [group1, group2, group3]
    stat_kw, p_kw = stats.kruskal(*data_list)
    if p_kw < 0.05:
        p_matrix = sp.posthoc_dunn(data_list, p_adjust='bonferroni')
        comparison_pairs = [(1, 2), (2, 3), (1, 3)]
        caps_y = [cap.get_ydata()[0] for cap in box['caps']]
        max_y = max(caps_y) if caps_y else max([max(g) for g in data_list])
        
        y_shift = max_y * 0.08
        start_y = max_y + y_shift
        
        for idx, (g1, g2) in enumerate(comparison_pairs):
            p_val = p_matrix.loc[g1, g2]
    
            if p_val < 0.01:
                p_text = f"P = {p_val:.1e}" if p_val < 0.001 else f"P = {p_val:.3f}"
            else:
                p_text = "n.s."
                
            current_y = start_y + idx * y_shift * 1.5
            h = y_shift * 0.4 
            ax.plot([g1, g1, g2, g2], [current_y - h, current_y, current_y, current_y - h], 
                    color='black', linewidth=0.6)
            ax.text((g1 + g2) / 2, current_y + h * 0.2, p_text, 
                    ha='center', va='bottom', fontsize=7.5, color='black')
        
        ax.set_ylim(top=start_y + len(comparison_pairs) * y_shift * 1.6)

    plt.tight_layout()
    plt.savefig("/data/zhangjy/DEAF1/Histone_Analysis/Pdf/"+ctcf['labels'][i]+'.lineplot.pdf',
                        bbox_inches = 'tight',
                        facecolor='w')
