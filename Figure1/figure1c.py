import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from matplotlib_venn import venn2, venn2_circles

# 设置默认样式与自定义字体
style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
plt.rcParams['font.sans-serif'] = 'Helvetica'

# 颜色与画布配置
colors = ['#6CBEC7', '#825B32']
fig, ax = plt.subplots(figsize=(4, 4))

# 绘制韦恩图
v = venn2(subsets=(1, 1, 1), set_labels=(' ', ' '))

for area in ['01', '10', '11']:
    v.get_patch_by_id(area).set_alpha(0.9)
    v.get_label_by_id(area).set_size(13)

v.get_patch_by_id('10').set_color(colors[1])
v.get_patch_by_id('01').set_color(colors[0])

# 设置文本内容与颜色
v.get_label_by_id('10').set_text(str(40250-2102) + '\nCTCF')
v.get_label_by_id('01').set_color('#0D5EA6') 
v.get_label_by_id('11').set_color('#0D5EA6') 
v.get_label_by_id('01').set_text(str(3657-2102) + '\nDEAF1')
v.get_label_by_id('11').set_text('2102\n(' + f'{2102/3657:.2%}' + ')')

v.get_patch_by_id('11').set_alpha(0.3)
c = venn2_circles(subsets=(1, 1, 1), color="#F3F2EC", linewidth=2)

fig.show()
fig.savefig(
    '/data/zhangjy/DEAF1/ChiPseq_Analysis/ResultsSort/Deeptools/InputNew/Pdf/WT_DEAF1_CTCF.pie.pdf', 
    bbox_inches='tight',
    facecolor='w'
)
