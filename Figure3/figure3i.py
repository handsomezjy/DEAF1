import matplotlib.font_manager as font_manager
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import fanc.plot as fancplot
from fanc.architecture.regions import GenomicRegion

def FancTrackLinePlot(filename, color, label, fill=True, ylabel=None, ylim=None, ax=None):
    track = fanc.load(filename)
    hp = fancplot.LinePlot(track, fill=fill, style='mid', colors=[color], labels=[label], ylabel=ylabel, ylim=ylim, n_yticks=2, ax=ax)
    return hp

style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
plt.rcParams['font.sans-serif']='Helvetica'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['axes.unicode_minus'] =False

region = GenomicRegion('chr1', 860000, 920000)
fanc_region = f'{region.chromosome}:{(region.start/1000000):.2f}mb-{(region.end/1000000):.2f}mb'
fig = plt.figure(figsize=(3, 3.5))
gs = gridspec.GridSpec(5, 1,  height_ratios=[1,0.8,0.8,0.8,0.8,],hspace=0.05)#list(np.ones(9)),wspace=0.4

ax1=fig.add_subplot(gs[0,0])
filename = "/data/zhangjy/Reference/Human/hg38/gencode.v29.primary_assembly.annotation_UCSC_names.gtf"
p = fancplot.GenePlot(filename, group_by='gene_id', squash=True, label_field='gene_name',
                        color_forward='#999999', color_reverse='#999999', show_arrows=True,
                        arrow_size=4, line_width=1, box_height=1,
                        # relative_marker_step=0.004*2.54, 
                       font_size=7,ax=ax1)
p.plot(fanc_region)
plt.xticks((700000,800000,900000),('Chr11:0.7mb','0.8mb','0.9mb'))
ax1.xaxis.tick_top()
ax1.xaxis.set_ticks_position('top')###刻度方向
ax1.xaxis.set_label_position('top') 
ax1.spines['bottom'].set_visible(False)
ax1.spines['top'].set_visible(True)
ax1.minorticks_off()

ax2=fig.add_subplot(gs[1,0])
filename="/data/zhangjy/DEAF1/ChiPseq_Analysis/Step5_combineRep/Bigwig/E250604004_L01_WT_DEAF1.bw"
p = FancTrackLinePlot(filename, '#34699A', 'WT DEAF1', True, "WT DEAF1", ylim=(0,2), ax=ax2)
p.plot(fanc_region)
ax2.set_ylabel('WT DEAF1                    ', rotation=0, x=-0.9, y=0.35,color='#34699A')
ax2.set_xticks([])
ax2.set_yticks([])
ax2.spines[:].set_visible(False)
plt.text(x=910000,y=0.9,s='[0,2]',color='#34699A')

filename="/data/zhangjy/DEAF1/ChiPseq_Analysis/Step5_combineRep/Bigwig/E250609001_L01_DEAF1_DEAF1.bw"
ax3=fig.add_subplot(gs[2,0])
p = FancTrackLinePlot(filename, '#FBA834', 'WT CTCF', True, "WT DEAF1", ylim=(0,2), ax=ax3)
p.plot(fanc_region)
ax3.set_ylabel('DEAF1(DEAF1 mutant)                    ', rotation=0, x=-1, y=0.35,color='#FBA834')
ax3.set_xticks([])
ax3.set_yticks([])
ax3.spines[:].set_visible(False)
plt.text(x=910000,y=0.9,s='[0,2]',color='#FBA834')

filename="/data/zhangjy/DEAF1/ChiPseq_Analysis/Step5_combineRep/Bigwig/E250604004_L01_WT_CTCF.bw"
ax4=fig.add_subplot(gs[3,0])
p = FancTrackLinePlot(filename, '#077A7D', 'WT CTCF', True, "WT DEAF1", ylim=(0,2), ax=ax4)
p.plot(fanc_region)
ax4.set_ylabel('WT CTCF                     ', rotation=0, x=-1, y=0.35,color='#077A7D')
ax4.set_xticks([])
ax4.set_yticks([])
ax4.spines[:].set_visible(False)
plt.text(x=910000,y=0.9,s='[0,2]',color='#077A7D')

filename="/data/zhangjy/DEAF1/ChiPseq_Analysis/Step5_combineRep/Bigwig/E250609001_L01_DEAF1_CTCF.bw"
ax5=fig.add_subplot(gs[4,0])
p = FancTrackLinePlot(filename, '#BED754', 'CTCF(DEAF1 mutant)', True, "CTCF(DEAF1 mutant)", ylim=(0,2), ax=ax5)
p.plot(fanc_region)
ax5.set_ylabel('CTCF(DEAF1 mutant)                    ', rotation=0, x=-1, y=0.35,color='#BED754')
ax5.set_xticks([])
ax5.set_yticks([])
ax5.spines[:].set_visible(False)
plt.text(x=910000,y=0.9,s='[0,2]',color='#BED754')

fig.savefig('/data/zhangjy/DEAF1/Figure1_Analysis/Pdf/IGV.CTCFDEAF1.sitesLostExam1.genomictracks.pdf', 
            bbox_inches = 'tight',
            facecolor='w')
