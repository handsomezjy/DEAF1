import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import matplotlib.patches as mpatches

# 读取 deeptools 导出的矩阵文件
ma = pd.read_table(
    "/data/zhangjy/DEAF1/ChiPseq_Analysis/ResultsSort/Deeptools/InputNew/Test/DEAF1before_allResults.RPGC.gz",
    skiprows=1,
    header=None
)

# 计算 95% 置信区间 (CI)
sub_df1 = ma.iloc[0:2102, 6:206]
ci_lower1, ci_upper1 = stats.t.interval(
    0.95, df=len(sub_df1) - 1, loc=sub_df1.mean(axis=0), scale=sub_df1.sem(axis=0)
)

sub_df2 = ma.iloc[2102:, 6:206]
ci_lower2, ci_upper2 = stats.t.interval(
    0.95, df=len(sub_df2) - 1, loc=sub_df2.mean(axis=0), scale=sub_df2.sem(axis=0)
)

# 样式与字体设置
plt.style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
plt.rcParams['font.sans-serif'] = 'Helvetica'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['pdf.fonttype'] = 42

# 绘图
fig, ax = plt.subplots(figsize=(3.5/2, 3/2), constrained_layout=False)

(p1,) = plt.plot(list(range(1, 201)), ma.iloc[0:2102, 6:206].mean(), color='#D2C1B6', linewidth="1.7")
(p2,) = plt.plot(list(range(1, 201)), ma.iloc[2102:, 6:206].mean(), color='#6CBEC7', linewidth="1.7")

ax.fill_between(list(range(1, 201)), ci_lower1, ci_upper1, color='#D2C1B6', alpha=0.4, label='Interval')
ax.fill_between(list(range(1, 201)), ci_lower2, ci_upper2, color='#6CBEC7', alpha=0.4, label='Interval')

patch1 = mpatches.Patch(color='#D2C1B6', alpha=0.35, linewidth=0)
patch2 = mpatches.Patch(color='#6CBEC7', alpha=0.35, linewidth=0)

ax.legend(
    bbox_to_anchor=(1.4, 1.25),
    handles=[(p1, patch1), (p2, patch2)],
    labels=["DEAF1-CTCF", "DEAF1-only"],
    frameon=False,
    fontsize=8,
    handletextpad=0.4,
    borderpad=0.2,
    loc="upper right",
)

ax.text(0.65, 0.91, "Line: mean\nShade: 95% CI", transform=ax.transAxes, fontsize=8, color="#555555", va="top", ha="left")

plt.ylim(0, 25)
plt.xticks([1, 100.5, 200], ['-1.0', 0, '1.0'])
plt.yticks([0, 10, 20])
plt.xlabel('distance to peak center(kb)')

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

plt.title('DEAF1 ChIPseq signal', y=1.2)

plt.savefig(
    '/data/zhangjy/DEAF1/ChiPseq_Analysis/ResultsSort/Deeptools/InputNew/Pdf/lineplot.beforeDEAF1.95CI.pdf',
    bbox_inches='tight',
    facecolor='w'
)
plt.show()

ma=pd.read_table("/data/zhangjy/DEAF1/ChiPseq_Analysis/ResultsSort/Deeptools/InputNew/Test/CTCFbefore_allResults.RPGC.gz",skiprows=1,header=None)
sub_df1 = ma.iloc[0:2102, 6:206]
ci_lower1, ci_upper1 = stats.t.interval(
    0.95, df=len(sub_df1) - 1, loc=sub_df1.mean(axis=0), scale=sub_df1.sem(axis=0)
)
sub_df2 = ma.iloc[2102:,6:206]
ci_lower2, ci_upper2 = stats.t.interval(
    0.95, df=len(sub_df2) - 1, loc=sub_df2.mean(axis=0), scale=sub_df2.sem(axis=0)
)
style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
plt.rcParams['font.sans-serif']='Helvetica'
plt.rcParams['pdf.fonttype'] = 42
# 高斯滤波平滑矩阵
plt.rcParams['axes.unicode_minus'] = False
fig,ax = plt.subplots(figsize=(3.5/2,3/2),constrained_layout=False)
(p1,) =plt.plot(list(range(1,201)),ma.iloc[0:2102,6:206].mean(),color='#D2C1B6', linewidth = "1.7")
(p2,) =plt.plot(list(range(1,201)),ma.iloc[2102:,6:206].mean(),color='#825B32', linewidth = "1.7")

ax.fill_between(list(range(1,201)), ci_lower1, ci_upper1, color='#D2C1B6', alpha=0.4)
ax.fill_between(list(range(1,201)), ci_lower2, ci_upper2, color='#825B32', alpha=0.4)
patch1 = mpatches.Patch(color='#D2C1B6', alpha=0.35, linewidth=0)
patch2 = mpatches.Patch(color='#825B32', alpha=0.35, linewidth=0)
ax.legend(bbox_to_anchor=(1.4, 1.),
    handles=[(p1, patch1), (p2, patch2)],
    labels=["DEAF1-CTCF", "CTCF-only"],  
    frameon=False,  
    fontsize=8,  
    handletextpad=0.4,
    borderpad=0.2,
    loc="upper right",  
)
ax.text(0.65,0.65,"Line: mean\nShade: 95% CI",transform=ax.transAxes,fontsize=8,color="#555555",va="top",ha="left",)

plt.ylim(0,25)
plt.yticks([0,10,20])
plt.xticks([1,100.5,200],['-1.0',0,'1.0'])#
plt.xlabel('distance to peak center(kb)')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.title('CTCF ChIPseq signal',y=1.1)
plt.savefig('/data/zhangjy/DEAF1/ChiPseq_Analysis/ResultsSort/Deeptools/InputNew/Pdf/lineplot.beforeCTCF.95CI.pdf',   ##
                bbox_inches = 'tight',
                facecolor='w')  
plt.show()
