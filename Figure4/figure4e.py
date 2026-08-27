import math
import matplotlib.font_manager as font_manager
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import pyBigWig

from matplotlib import style

import fanc.plot as fancplot
from fanc.architecture.regions import GenomicRegion

def FancTrackLinePlot(filename, color, label, fill=True, ylabel=None, ylim=None, ax=None):
    track = fanc.load(filename)
    hp = fancplot.LinePlot(track, fill=fill, style='mid', colors=[color], labels=[label], ylabel=None, ylim=ylim, n_yticks=2, ax=ax)
    return hp
def exampleGene(chrom,sta,end,gene):
    style.use('default')
    font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
    plt.rcParams['font.sans-serif']='Helvetica'
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['axes.unicode_minus'] =False

    region = GenomicRegion(chrom,sta,end)
    fanc_region = f'{region.chromosome}:{(region.start/1000000):.2f}mb-{(region.end/1000000):.2f}mb'
    fig = plt.figure(figsize=(3, 4.5))
    gs = gridspec.GridSpec(7, 1,  height_ratios=[1,0.8,0.8,0.8,0.8,0.8,0.8],hspace=0.05)#list(np.ones(9)),wspace=0.4

    ax1=fig.add_subplot(gs[0,0])
    filename = "/data/Public_Data/References/human/hg38/gencode.v38.protein_coding.annotation.gtf"
    p = fancplot.GenePlot(filename, group_by='gene_name', squash=True, label_field='gene_name',
                            color_forward='#999999', color_reverse='#999999', show_arrows=True,
                            arrow_size=4, line_width=1, box_height=1,
                            # relative_marker_step=0.004*2.54, 
                        font_size=7,ax=ax1)
    p.plot(fanc_region)
    plt.xticks((sta, end),(f"{chrom}:{sta/1e6:.2f} Mb", f"{end/1e6:.2f} Mb"))
    ax1.xaxis.tick_top()
    ax1.xaxis.set_ticks_position('top')###刻度方向
    ax1.xaxis.set_label_position('top') 
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['top'].set_visible(True)
    ax1.minorticks_off()

    ax2=fig.add_subplot(gs[1,0])
    filename="/data/zhangjy/DEAF1/RNAseq_Analysis/Step4_BW/Combined/wt_RPKM.bw"
    p=fancplot.LinePlot(filename, fill=True, style='mid', colors='#3C3D37', 
                        labels=None, ylabel='WT DEAF1', ylim=(0,50), n_yticks=2, ax=ax2)
    p.plot(fanc_region)
    ax2.set_ylabel('WT RNA-seq                    ', rotation=0, x=-0.9, y=0.35,color='#3C3D37')
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.spines[:].set_visible(False)
    ax2.text(0.98, 0.95,"0,50",transform=ax2.transAxes,ha='right',va='top',fontsize=9)

    ax3=fig.add_subplot(gs[2,0])
    filename="/data/zhangjy/DEAF1/RNAseq_Analysis/Step4_BW/Combined/ko_RPKM.bw"
    p=fancplot.LinePlot(filename, fill=True, style='mid', colors='#B8001F', 
                        labels=None, ylabel='WT DEAF1', ylim=(0,50), n_yticks=2, ax=ax3)
    p.plot(fanc_region)
    ax3.set_ylabel('RNA-seq(DEAF1 mutant)                    ', rotation=0, x=-0.9, y=0.35,color='#B8001F')
    ax3.set_xticks([])
    ax3.set_yticks([])
    ax3.spines[:].set_visible(False)
    ax3.text(0.98, 0.95,"0,50",transform=ax3.transAxes,ha='right',va='top',fontsize=9)

    ax4=fig.add_subplot(gs[3,0])
    filename="/data/zhangjy/DEAF1/ChiPseq_Analysis/Step5_combineRep/Bigwig/BW_RPGC/E250604004_L01_WT_DEAF1_RPGC.bw"
    x_max, max_v = get_region_max_pos(filename, chrom, sta, end)
    # ctcflim=math.ceil(max_v / 10) * 10
    ctcflim=5
    p = FancTrackLinePlot(filename, '#34699A', 'WT CTCF', True, "WT DEAF1", ylim=(0,ctcflim), ax=ax4)
    p.plot(fanc_region)
    ax4.set_xticks([])
    ax4.set_yticks([])
    ax4.spines['left'].set_visible(False)
    ax4.text(0.98, 0.95,f"0,{ctcflim}",transform=ax4.transAxes,ha='right',va='top',fontsize=9)

    filename="/data/zhangjy/DEAF1/ChiPseq_Analysis/Step5_combineRep/Bigwig/BW_RPGC/E250609001_L01_DEAF1_DEAF1_RPGC.bw"
    x_max, max_v = get_region_max_pos(filename, chrom, sta, end)
    # ctcflim=math.ceil(max_v / 10) * 10
    ax5=fig.add_subplot(gs[4,0])
    p = FancTrackLinePlot(filename, '#FBA834', 'WT CTCF', True, "WT DEAF1", ylim=(0,ctcflim), ax=ax5)
    p.plot(fanc_region)
    ax5.set_xticks([])
    ax5.set_yticks([])
    ax5.spines['left'].set_visible(False)
    ax5.text(0.98, 0.95,f"0,{ctcflim}",transform=ax5.transAxes,ha='right',va='top',fontsize=9)

    filename="/data/zhangjy/DEAF1/ChiPseq_Analysis/Step5_combineRep/Bigwig/BW_RPGC/E250604004_L01_WT_CTCF_RPGC.bw"
    x_max, max_v = get_region_max_pos(filename, chrom, sta, end)
    ctcflim=8
    # ctcflim=math.ceil(max_v / 10) * 10
    ax6=fig.add_subplot(gs[5,0])
    p = FancTrackLinePlot(filename, '#077A7D', 'WT CTCF', True, "WT DEAF1", ylim=(0,ctcflim), ax=ax6)
    p.plot(fanc_region)
    ax6.set_xticks([])
    ax6.set_yticks([])
    ax6.spines['left'].set_visible(False)
    ax6.text(0.98, 0.95,f"0,{ctcflim}",transform=ax6.transAxes,ha='right',va='top',fontsize=9)

    filename="/data/zhangjy/DEAF1/ChiPseq_Analysis/Step5_combineRep/Bigwig/BW_RPGC/E250609001_L01_DEAF1_CTCF_RPGC.bw"
    x_max, max_v = get_region_max_pos(filename, chrom, sta, end)
    # ctcflim=math.ceil(max_v / 10) * 10
    ax7=fig.add_subplot(gs[6,0])
    p = FancTrackLinePlot(filename, '#BED754', 'WT CTCF', True, "WT DEAF1", ylim=(0,ctcflim), ax=ax7)
    p.plot(fanc_region)
    ax7.set_xticks([])
    ax7.set_yticks([])
    ax7.spines['left'].set_visible(False)
    ax7.text(0.98, 0.95,f"0,{ctcflim}",transform=ax7.transAxes,ha='right',va='top',fontsize=9)

    li1=[ax1,ax2,ax3,ax4,ax5,ax6,ax6,ax7]#
    for i  in li1:
        for line in i.lines:
            line.set_rasterized(True)
        for col in i.collections:
            col.set_rasterized(True)

    fig.savefig(f'/data/zhangjy/DEAF1/RNAseq_Analysis/Pdf/{gene}.genomictracks.pdf', 
                bbox_inches = 'tight',
                facecolor='w')

exampleGene('chr9', 76380000, 76800000,"Downgene_GCNT1")
exampleGene('chr1', 155600000, 155900000,"Upgene_SYT11")
