import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

# 数据与样式配置
plt.style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
plt.rcParams['font.sans-serif'] = 'Helvetica'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['pdf.fonttype'] = 42

# 绘图初始化
fig, ax = plt.subplots(figsize=(3/1.8, 3/1.8), constrained_layout=False)

# 绘制柱状图
plt.bar([1, 2], [1196/1771, 34394/38479], color=['#FB9E3A', '#FB9E3A'], width=0.4, zorder=3)

# 刻度与范围配置
plt.yticks((0, 0.25, 0.5, 0.75, 1), (0, '25', '50', '75', 100))
plt.ylim(0, 1)
plt.xlim(0.5, 2.5)
plt.xticks([1, 2], [r'$DEAF1^+$' + '\n' + r'$CTCF^+$', r'$DEAF1^+$' + '\n' + r'$CTCF^-$'], rotation=0)

# 坐标轴与边框美化
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_color('grey')
ax.spines['bottom'].set_color('grey')
ax.tick_params(axis='both', color='grey') 

plt.ylabel('CTCF retained(%)')

# 保存与展示
fig.show()
fig.savefig(
    '/data/zhangjy/DEAF1/Figure1_Analysis/Pdf/Pdf_new/CTCFretained_bar.pdf', 
    bbox_inches='tight',
    facecolor='w'
)

style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
plt.rcParams['font.sans-serif']='Helvetica'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['pdf.fonttype'] = 42
fig,ax = plt.subplots(figsize=(3.5/2,3/2),constrained_layout=False)

plt.bar([1,2],[1768/2102,697/1477],color=['#406093','#406093'],width=0.4,zorder=3)
plt.yticks((0,0.25,0.5,0.75,1),(0,'25','50','75',100))
plt.ylim(0,1)
plt.xlim(0.5,2.5)
plt.xticks([1,2],[r'$DEAF1^+$' + '\n' + r'$CTCF^+$', r'$DEAF1^-$' + '\n' + r'$CTCF^+$'],rotation=0)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_color('grey')
ax.spines['bottom'].set_color('grey')
ax.tick_params(axis='both', color='grey') 
plt.ylabel('DEAF1 retained(%)')
fig.show()
fig.savefig('/data/zhangjy/DEAF1/Figure1_Analysis/Pdf/Pdf_new/DEAF1retained_bar.pdf', 
                bbox_inches = 'tight',
                facecolor='w')
