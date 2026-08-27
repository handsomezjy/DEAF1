import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

deseq=pd.read_table("/data/zhangjy/DEAF1/RNAseq_Analysis/Step3_RPKM/total.gene_counts.Deseq2results.csv")
deseq=deseq.dropna()
deseq['-log10(p-value)']=-np.log10(deseq['pvalue'])

max_col=deseq['-log10(p-value)'].replace([np.inf, -np.inf], np.nan).max()
deseq['-log10(p-value)'] = deseq['-log10(p-value)'].replace({np.inf: max_col, -np.inf: max_col})

tmp=pd.DataFrame(deseq[(deseq['log2FoldChange']>1)&(deseq['pvalue']<0.05)].index)

rnasequp=deseq[(deseq['log2FoldChange']>1)&(deseq['pvalue']<0.05)][['log2FoldChange','-log10(p-value)']].reset_index()

rnaseqdown=deseq[(deseq['log2FoldChange']<-1)&(deseq['pvalue']<0.05)][['log2FoldChange','-log10(p-value)']].reset_index()

style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
plt.rcParams['font.sans-serif']='Helvetica'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['axes.unicode_minus'] =False

fig,ax = plt.subplots(figsize=(3,4),constrained_layout=False)
tmp1=deseq[(deseq['log2FoldChange']>1)&(deseq['pvalue']<0.05)]
plt.scatter(x=tmp1['log2FoldChange'],y=tmp1['-log10(p-value)'],c='#DC3C22',s=1,
            label='Up ('+str(tmp1.shape[0])+')',rasterized=True)

tmp2=deseq[(deseq['log2FoldChange']<-1)&(deseq['pvalue']<0.05)]
plt.scatter(x=tmp2['log2FoldChange'],y=tmp2['-log10(p-value)'],c='#3D74B6',s=1,
            label='Down ('+str(tmp2.shape[0])+')',rasterized=True)

li=list(set(deseq.index.to_list())-set(tmp1.index.to_list()+tmp2.index.to_list()))
tmp3=deseq.loc[li,:]
plt.scatter(x=tmp3['log2FoldChange'],y=tmp3['-log10(p-value)'],c='#7A7A73',s=1,label='NoSignifi ('+str(tmp3.shape[0])+')')
plt.vlines(x=[-1,1],ymin=-20,ymax=300,colors='k',linestyles='--',linewidth=1)
plt.hlines(y=-np.log10(0.05),xmin=-15,xmax=15,colors='k',linestyles='--',linewidth=1)
plt.ylim(-20,300)
plt.xlim(-15,15)
plt.xticks([-10,0,10])#
plt.yticks([0,100,200,300])
plt.xlabel('distance to peak center(kb)')
plt.legend(frameon=False, loc='upper left', fontsize=10,
           bbox_to_anchor=(0.55, 0.95),markerscale=5,
           labelspacing=0.5,handletextpad=0.2)# 
plt.ylabel(r'$-\log_{10}(\mathrm{p-value})$', fontsize=12)
plt.xlabel(r'$\log_2(\mathrm{fold\ change})$', fontsize=12)
plt.grid(axis='both', color='grey', linestyle=':',zorder=0)

plt.savefig('/data/zhangjy/DEAF1/RNAseq_Analysis/Pdf/Volcano_DEAF1mutant.pdf',   ##
                bbox_inches = 'tight',
                facecolor='w')  
plt.show()
