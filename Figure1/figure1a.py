import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import style

def calculate_ci_z(row):
    n = len(row)                    # 样本量
    mean = row.mean()               # 均值
    std = row.std(ddof=1)           # 样本标准差（ddof=1表示无偏估计）
    sem = std / np.sqrt(n)          # 标准误
    z = 1.96                        # 95%置信水平的Z值
    margin = z * sem                # 置信区间半宽
    return pd.Series({
        'mean': mean,
        'ci_lower': mean - margin,
        'ci_upper': mean + margin
    })

## 应用函数并生成结果
ma=pd.read_table("/data/zhangjy/CTCF/TRN/Cells_bedtools/BedtoolsK562/TF_intersectTF/tf_overlap_matrix.tsv")
ma.set_index('TF',inplace=True)
df=ma.apply(calculate_ci_z, axis=1)
df['CTCF']=ma['CTCF']


df=df.sort_values(by='ci_upper',ascending=False)

style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
plt.rcParams['font.sans-serif']='Helvetica'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['axes.unicode_minus'] = False

fig, ax = plt.subplots(figsize=(4, 2.5))
plt.plot(np.arange(1,df.shape[0]+1,1),df['CTCF'],color="#DDA853",linewidth=1.5,label='CTCF')
plt.plot(np.arange(1,df.shape[0]+1,1),df['ci_upper'],color="#B2B2B2",linewidth=1,label='95% CI')
#confidence interval
plt.fill_between(
    np.arange(1,df.shape[0]+1,1),df['ci_upper'], y2=0, color='#B2B2B2', alpha=0.8
)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_position(('data',-800))
ax.spines['left'].set_position(('data',-10))

plt.ylim(0,df['CTCF'].max()+1000)
plt.xlim(0,df.shape[0])
plt.xticks([])
plt.ylabel('Intersection Num')
plt.xlabel('Transcription Factors')
plt.hlines(int(df.loc['CTCF','ci_upper']),0,df.shape[0]
           ,colors='k',linestyles='dotted')
plt.text(320,
         int(df.loc['CTCF','ci_upper'])+1000,
         'CTCF 95% CI : {x}'.format(x=int(df.loc['CTCF','ci_upper'])),
         fontsize=9)
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
tmp=df.sort_values(by='CTCF',ascending=False)#.iloc[20:30]
tmpli=['RAD21','REST','SMC3','MAZ','ZNF143','MAX','DEAF1','MYC','CTCFL','YY1','RB1']
tmp=tmp.loc[tmpli,:]
for i,j in tmp.iterrows():
    plt.scatter(df.index.get_loc(i),
                j['CTCF'],
                color='k',
                s=5,zorder=2)
    plt.text(df.index.get_loc(i),
            j['CTCF']+1000,
            i,
            fontsize=9)
    
plt.legend(loc='upper right',frameon=False) 
plt.savefig('/data/zhangjy/CTCF/TRN/Cells_bedtools/BedtoolsK562/TF_intersectTF/CTCF_TFs526.95CI.lineplot.pdf',
                bbox_inches = 'tight',
                facecolor='w')
