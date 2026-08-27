import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt
from matplotlib import style

style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
plt.rcParams['font.sans-serif']='Helvetica'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['axes.unicode_minus'] =False
fig, ax = plt.subplots(figsize=(2, 3))

colors=['#825B32','#6CBEC7','#FFDBB5','#7F8487']
labels=[r'$DEAF1^{-} CTCF^{+}$',r'$DEAF1^{+} CTCF^{-}$',r'$DEAF1^{+} CTCF^{+}$','None']
plt.bar(1,210/469,color=colors[0],width=0.5,label=labels[0])
plt.text(0.8,210/469/2,s='210')
plt.bar(1,9/469,color=colors[1],width=0.5,bottom=210/469,label=labels[1])
plt.text(0.8,0.42,s='9')
plt.bar(1,3/469,color=colors[2],width=0.5,bottom=(210+9)/469,label=labels[2])
plt.text(0.8,0.48,s='3')
plt.bar(1,(469-(210+9+3))/469,color=colors[3],width=0.5,bottom=(210+9+3)/469,label=labels[3])
plt.text(0.8,0.8,s=str(469-(210+9+3)))

plt.bar(2,255/608,color=colors[0],width=0.5)
plt.text(1.8,255/608/2,s='255')
plt.bar(2,8/608,color=colors[1],width=0.5,bottom=255/608)
plt.text(1.8,0.38,s='8')
plt.bar(2,1/608,color=colors[2],width=0.5,bottom=(255+8)/608)
plt.text(1.8,0.43,s='1')
plt.bar(2,(608-(255+8+1))/608,color=colors[3],width=0.5,bottom=(255+8+1)/608)
plt.text(1.8,0.8,s=str(608-(255+8+1)))

locs, labels =plt.xticks([1,2],['Up-regulated\ngenes','Down-regulated\ngenes'])
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
labels[0].set_color('#F63049') 
labels[1].set_color('#008BFF')
plt.yticks([0,0.25,0.5,0.75,1],[0,'25%','50%','75%','100%'])
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(
    handles[::-1], 
    labels[::-1],
    loc='lower center',          
    bbox_to_anchor=(1.5, 0.5), 
    ncol=1,                       
    frameon=False,                
    fontsize=10
)
plt.ylim(0,1)
plt.text(2.5,0,s='TSS $\pm$ 2kb')
plt.savefig("/data/zhangjy/DEAF1/RNAseq_Analysis/Pdf/DRE_TSS.DEAF1_CTCF.binding.bar.pdf",
            bbox_inches = 'tight',
            facecolor='w')
