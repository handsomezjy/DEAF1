import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

# 数据定义
tot1 = 40250
y1 = np.array([15.3888199, 3.1652174, 2.6782609, 0.4173913, 3.3664596,
               2.5416149, 5.4509317, 13.6844720, 27.1677019, 0.1515528, 25.9875776])

tot2 = 3657
y2 = np.array([8.55892808, 0.41017227, 0.32813782, 0.08203445, 0.08203445,
               0.65627564, 0.21875855, 2.43368882, 5.52365327, 0, 81.70631665])

tmp = pd.DataFrame([y1, y2])

# 样式与字体设置
plt.style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
plt.rcParams['font.sans-serif'] = 'Helvetica'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['pdf.fonttype'] = 42

# 颜色与图例标签
colors = [
    '#347928', '#C0EBA6', '#E7FBE6', '#1A4870', '#81DAE3',
    '#AD49E1', '#EBD3F8', '#821131', '#EE66A6', '#F6FB7A', '#FCCD2A'
]
labels = [
    'Promoter (<=1kb)', 'Promoter (1-2kb)', 'Promoter (2-3kb)',
    "5' UTR", "3' UTR", "1st Exon",
    'Other Exon', '1st Intron', 'Other Intron',
    'Downstream (<=300)', 'Distal Intergenic'
]

# 绘图
fig, ax = plt.subplots(figsize=(4/2, 1/2), dpi=150)
plt.grid(axis='x', color='#F6F1E9', linestyle=':', zorder=0)

for i in range(tmp.shape[1]):
    if i == 0:
        plt.barh(tmp.index, tmp[i], color=colors[i], edgecolor=None, label=labels[i])
    else:
        plt.barh(tmp.index, tmp[i], left=tmp.loc[:, :(i-1)].sum(axis=1), color=colors[i], edgecolor=None, label=labels[i])

plt.xticks([0, 25, 50, 75, 100], [0, '25%', '50%', '75%', '100%'], fontsize=8)
plt.yticks([0, 1], ['CTCF', 'DEAF1'], fontsize=8)

legend = plt.legend(
    bbox_to_anchor=(0.5, 3),
    fontsize=6,
    loc='upper center',
    ncol=3,
    frameon=False,
    columnspacing=0.4
)
legend.texts[0].set_weight('bold')

fig.show()
fig.savefig(
    '/data/zhangjy/DEAF1/Figure1_Analysis/Pdf/Pdf_new/DEAF1_CTCF.allpeaks.chipseeker.pdf',
    bbox_inches='tight',
    facecolor='w'
)

