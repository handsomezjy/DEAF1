import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import style

wt=pd.read_table("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/cooltools/CompartmentAnalysis/WT.MicroC.100kb.cis.vecs.tsv").dropna()
ko=pd.read_table("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/cooltools/CompartmentAnalysis/KO.MicroC.100kb.cis.vecs.tsv").dropna()
df=pd.merge(wt.iloc[:,[0,1,2,4]],ko.iloc[:,[0,1,2,4]],on=['chrom','start','end'])

style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
plt.rcParams['font.sans-serif']='Helvetica'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['axes.unicode_minus'] =False

fig,ax = plt.subplots(figsize=(3,3),constrained_layout=False)
plt.scatter(df['E1_x'],df['E1_y'],s=10,c='#456882',alpha=0.7,rasterized=True)
plt.xlim(-2,8)
plt.ylim(-2,8)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_color('grey')
ax.spines['bottom'].set_color('grey')
corr = np.corrcoef(df['E1_x'],df['E1_y'])[0, 1]
print(corr)
ax.plot([-2,7.7], [-2,7.7], color='k', lw=1,linestyle='--')
ax.spines['left'].set_position(('data',-2.3))
ax.spines['bottom'].set_position(('data',-2.3))

plt.xlabel('WT compartment scores',fontsize=11)
plt.ylabel('DEAF1 mutant compartment scores',fontsize=11)
ax.text(x=0.1,y=0.8,s='Pearson r={:.3f}'.format(corr), transform=ax.transAxes,fontsize=12)
plt.savefig('/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/Pdf/CompartmentScores.pearson.scatter.pdf',   ##
                bbox_inches = 'tight',
                facecolor='w')  
