import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

# 数据界限与矩阵读取
group_boundaries = [0, 38027, 39403, 41367]
matrix = pd.read_table(
    "/data/zhangjy/DEAF1/Figure1_Analysis/Deeptools/MNase_ATAC/Newinput/matrix_ATAC.gz",
    skiprows=1,
    header=None
)

# 字体与样式配置
plt.style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
plt.rcParams['font.sans-serif'] = 'Helvetica'

# 绘图初始化
fig, ax = plt.subplots(figsize=(2.5, 2), dpi=100)
colors = ['#825B32', '#6CBEC7', '#FFDBB5']
group_labels = [r'$DEAF1^{-} CTCF^{+}$', r'$DEAF1^{+} CTCF^{-}$', r'$DEAF1^{+} CTCF^{+}$']

# 循环计算各组 95% 置信区间并绘制折线与阴影
for i in range(3):
    sub_df1 = matrix.iloc[group_boundaries[i]:group_boundaries[i+1], 6:]
    ci_lower1, ci_upper1 = stats.t.interval(
        0.95, df=len(sub_df1) - 1, loc=sub_df1.mean(axis=0), scale=sub_df1.sem(axis=0)
    )
    (p1,) = plt.plot(list(range(400)), sub_df1.mean(), color=colors[i], label=group_labels[i], lw=1.5)
    ax.fill_between(list(range(400)), ci_lower1, ci_upper1, color=colors[i], alpha=0.4)

# 坐标轴与标签美化
plt.xticks([0, 199, 399], ['-2.0', 'center', '+2.0'], fontsize=10)
plt.title('ATAC-seq', fontdict={'weight': 'normal', 'size': 13}, x=0.5, y=0.99)
plt.yticks((0, 100, 200))

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel('distance to peak center(kb)', fontsize=12)

# 偏移坐标轴起点位置
ax.spines['left'].set_position(('data', -5))
plt.xlim(0, 399)
plt.ylim(0, 200)
ax.spines['bottom'].set_position(('data', -10))

fig.show()
fig.savefig(
    '/data/zhangjy/DEAF1/Figure1_Analysis/Deeptools/MNase_ATAC/Newinput/DEAF1_CTCFaloneShared.ATACline.95CI.pdf', 
    bbox_inches='tight',
    facecolor='w'
)
