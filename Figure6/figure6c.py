import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt
from matplotlib import style

loops=pd.read_csv("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/Bedtools/loop_anchors_feature_stats.csv")
df=loops[loops['Strengthened']==1]
tmp=df[((df['A1_TSS']==1)&((df['A2_H3K27ac']==1)))|((df['A2_TSS']==1)&((df['A1_H3K27ac']==1)))].iloc[:,0:6]
tmp1=df.drop(index=tmp.index)
tmp1=tmp1[((tmp1['A1_TSS']==1)&((tmp1['A2_TSS']==1)))|((tmp1['A2_TSS']==1)&((tmp1['A1_TSS']==1)))].iloc[:,0:6]

style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
plt.rcParams['font.sans-serif']='Helvetica'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['axes.unicode_minus'] =False
# fig, ax = plt.subplots(figsize=(1, 3))
values = [df.shape[0]-tmp.shape[0]-tmp1.shape[0],tmp.shape[0],tmp1.shape[0]]  
labels = ['None','EP','PP']
colors = ['#E0D9D9','#E4D329','#0D0B61']
fig, ax = plt.subplots(figsize=(2,2))

wedges, texts, autotexts = ax.pie(
    values,
    # labels=labels,
    colors=colors[:len(values)],
    autopct='%1.1f%%',
    startangle=0,
    counterclock=False,
    # labeldistance=0.5
)
x, y = autotexts[2].get_position()
autotexts[2].set_position((x*2.3, y*1.8))
ax.set_title(f'Strengthened loops\nn=4244', fontsize=10,y=0.9)
ax.legend(wedges,labels,loc='center left',bbox_to_anchor=(0.9, 0.8),
    frameon=False,fontsize=8,handlelength=1,handletextpad=0.5)
fig.show()
fig.savefig('/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/Pdf/Strengthened.Ep_PP.pie.pdf', 
                bbox_inches = 'tight',
                facecolor='w')


df=loops[loops['Weakened']==1]
tmp=df[((df['A1_TSS']==1)&((df['A2_H3K27ac']==1)))|((df['A2_TSS']==1)&((df['A1_H3K27ac']==1)))].iloc[:,0:6]
tmp1=df.drop(index=tmp.index)
tmp1=tmp1[((tmp1['A1_TSS']==1)&((tmp1['A2_TSS']==1)))|((tmp1['A2_TSS']==1)&((tmp1['A1_TSS']==1)))].iloc[:,0:6]

style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
plt.rcParams['font.sans-serif']='Helvetica'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['axes.unicode_minus'] =False
# fig, ax = plt.subplots(figsize=(1, 3))
values = [df.shape[0]-tmp.shape[0]-tmp1.shape[0],tmp.shape[0],tmp1.shape[0]]  
labels = ['None','EP','PP']
colors = ['#E0D9D9','#E4D329','#0D0B61']
fig, ax = plt.subplots(figsize=(2,2))

wedges, texts, autotexts = ax.pie(
    values,
    # labels=labels,
    colors=colors[:len(values)],
    autopct='%1.1f%%',
    startangle=0,
    counterclock=False,
    # labeldistance=0.5
)
x, y = autotexts[2].get_position()
autotexts[2].set_position((x*2.3, y*1.8))
ax.set_title(f'Weakened loops\nn=6487', fontsize=10,y=0.9)
ax.legend(wedges,labels,loc='center left',bbox_to_anchor=(0.9, 0.8),
    frameon=False,fontsize=8,handlelength=1,handletextpad=0.5)
fig.show()
fig.savefig('/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/Pdf/Weakened.Ep_PP.pie.pdf', 
                bbox_inches = 'tight',
                facecolor='w')
