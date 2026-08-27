import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import matplotlib.gridspec as gridspec

path='/data/zhangjy/DEAF1/RNAseq_Analysis/Step3_RPKM/'
li=['ko_1','ko_2','wt_1','wt_2']
files=os.listdir(path)

fin=pd.DataFrame()
fin2=pd.DataFrame()
def FPKM(df):
    return df.iloc[:,-1]/((df.iloc[:,-2]/1e3)*(df.iloc[:,-1].sum()/1e6))
ii=0
for i in files:
    if i.endswith('gene_counts.txt'):
        df=pd.read_table(path+li[ii]+".gene_counts.txt",skiprows=1)
        df[li[ii]]=FPKM(df)
        if ii==0:
            fin=df[['Geneid',li[ii]]]
            fin2=df.iloc[:,[0,6]]
        else:
            fin=pd.merge(df[['Geneid',li[ii]]],fin,on='Geneid')
            fin2=pd.merge(df.iloc[:,[0,6]],fin2,on='Geneid')
        ii+=1

gtf=pd.read_table("/data/zhangjy/Reference/Human/hg38/gencode.v48.annotation.gtf",skiprows=5,header=None)
gene=gtf[gtf[2]=='gene']
tmp=gene[8].str.split('"',expand=True)[[1,5]]
tmp.columns=['Geneid','Genename']
fin=pd.merge(tmp,fin,on=['Geneid'])

style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
plt.rcParams['font.sans-serif']='Helvetica'
plt.rcParams['pdf.fonttype'] = 42
fig = plt.figure(figsize=(5, 2))
gs = gridspec.GridSpec(1, 2,  wspace=0.3)
li=['wt_1','wt_2','ko_1','ko_2']
li2=['CTCF','DEAF1']
colors=['#3C3D37','#3C3D37','#B8001F','#B8001F']
for i in range(2):
    ax1=fig.add_subplot(gs[0,i])
    tmpli=fin[fin['Genename']==li2[i]].loc[:,li].iloc[0].values.tolist()
    y_values=[np.log2(i+0.01) for i in tmpli]
    plt.bar([1,2,3,4], y_values,color=colors,width=0.6,zorder=3)
    plt.grid(axis='y', color='grey', linestyle=':',zorder=0)
# plt.yticks((0,0.25,0.5,0.75,1),(0,'25','50','75',100))
    plt.ylim(0,6)
    plt.xticks((1,2,3,4),('Rep 1','Rep 2','Rep 1','Rep 2'),fontsize=9,rotation=0)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['left'].set_color('grey')
    ax1.spines['bottom'].set_color('grey')
    ax1.tick_params(axis='both', color='grey') 
    if i==0:
        plt.ylabel('log2(FPKM+0.01)')
        plt.hlines(y=3.8,xmin=1,xmax=2,colors='k')
        plt.text(x=1.2,y=3.95,s='WT',fontsize=10)
        plt.hlines(y=3.8,xmin=3,xmax=4,colors='k')
        plt.text(x=2.8,y=3.95,s='DEAF1 mutant',fontsize=10)
    else:
        plt.hlines(y=5.8,xmin=1,xmax=2,colors='k')
        plt.text(x=1.2,y=5.95,s='WT',fontsize=10)
        plt.hlines(y=5.8,xmin=3,xmax=4,colors='k')
        plt.text(x=2.8,y=5.95,s='DEAF1 mutant',fontsize=10)
    plt.title(li2[i]+' expression',y=1.05)


fig.show()
fig.savefig('/data/zhangjy/DEAF1/RNAseq_Analysis/Pdf/CTCF_DEAF1.expression_bar.pdf', 
                bbox_inches = 'tight',
                facecolor='w')
