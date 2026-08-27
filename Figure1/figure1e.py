mport math
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import matplotlib.gridspec as gridspec

import fanc.newplot as fancplot
from fanc.architecture.regions import GenomicRegion

def plotExample(chrom,sta,end):
    style.use('default')
    font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
    plt.rcParams['font.sans-serif']='Helvetica'
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['axes.unicode_minus'] =False
    region = GenomicRegion(chrom,sta,end)
    fanc_region = f'{region.chromosome}:{(region.start/1000000):.2f}mb-{(region.end/1000000):.2f}mb'
    fig = plt.figure(figsize=(3, 2))
    gs = gridspec.GridSpec(3, 1,  height_ratios=[1,1,1],hspace=0.01)#list(np.ones(9)),wspace=0.4
    ax1=fig.add_subplot(gs[0,0])
    filename = "/data/Public_Data/References/human/hg38/gencode.v38.protein_coding.annotation.gtf"
    p = fancplot.GenePlot(filename, group_by='gene_id', squash=True, label_field='gene_name',
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
    filename="/data/zhangjy/DEAF1/ChiPseq_Analysis/Step5_combineRep/Bigwig/BW_RPGC/E250604004_L01_WT_DEAF1_RPGC.bw"
    x_max, max_v = get_region_max_pos(filename, chrom, sta, end)
    deaf1lim=math.ceil(max_v / 10) * 10
    p = FancTrackLinePlot(filename, '#34699A', 'WT DEAF1', True, "WT DEAF1", ylim=(0,deaf1lim), ax=ax2)
    p.plot(fanc_region)
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.spines['left'].set_visible(False)
    ax2.text(0.98, 0.95,f"0,{deaf1lim}",transform=ax2.transAxes,ha='right',va='top',fontsize=9)

    filename="/data/zhangjy/DEAF1/ChiPseq_Analysis/Step5_combineRep/Bigwig/BW_RPGC/E250604004_L01_WT_CTCF_RPGC.bw"
    x_max, max_v = get_region_max_pos(filename, chrom, sta, end)
    ctcflim=math.ceil(max_v / 10) * 10
    ax4=fig.add_subplot(gs[2,0])
    p = FancTrackLinePlot(filename, '#077A7D', 'WT CTCF', True, "WT DEAF1", ylim=(0,ctcflim), ax=ax4)
    p.plot(fanc_region)
    ax4.set_xticks([])
    ax4.set_yticks([])
    ax4.spines['left'].set_visible(False)
    ax4.text(0.98, 0.95,f"0,{ctcflim}",transform=ax4.transAxes,ha='right',va='top',fontsize=9)

    pos_ax1 = ax2.get_position()
    pos_ax2 = ax4.get_position()
    mid_y = (pos_ax1.y0 + pos_ax1.y1) / 2
    label_x = pos_ax1.x0-0.1  
    fig.text(label_x, mid_y, 'DEAF1\n(RPGC)', rotation=0,ha='center', va='center', fontsize=10,color='#34699A')
    mid_y = (pos_ax2.y0 + pos_ax2.y1) / 2
    label_x = pos_ax2.x0-0.1  
    fig.text(label_x, mid_y, 'CTCF\n(RPGC)', rotation=0,ha='center', va='center', fontsize=10,color='#077A7D')
    li1=[ax1,ax2,ax4]#
    for i  in li1:
        for line in i.lines:
            line.set_rasterized(True)
        for col in i.collections:
            col.set_rasterized(True)

    fig.savefig('/data/zhangjy/DEAF1/Figure1_Analysis/Pdf/Pdf_new/Example.CTCF_DEAF1_cobinding.genomictracks_'+fanc_region+'.RPGC.pdf', 
                bbox_inches = 'tight',
                facecolor='w')

plotExample("chr1",2150000,2680000)
