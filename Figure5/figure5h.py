import math
import os
import tempfile
import bioframe
import cooler
import fanc
import fanc.plotting as fancplot
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig
from genomic_regions import GenomicRegion, as_region
from matplotlib import style
from matplotlib.colors import LinearSegmentedColormap

def Mb(x):
    mb=x/1e6
    return f"{mb:.2f} Mb"

hic_file1="/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/HicFiles/WT.contact_map.KR.hic"
hic_file2="/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/HicFiles/Downsampling/KO.contact_map.downsample.KR.hic"

insulation='/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/Domains/Cooltools/'
wt_insulation_table=pd.read_table(insulation+'WT.insulation_10kb.bedGraph',header=None,names=['chrom','start','end','log2_insulation_score'])
c_insulation_table=pd.read_table(insulation+'KO.insulation_10kb.downsample.bedGraph',header=None,names=['chrom','start','end','log2_insulation_score'])

def get_region_max_pos(bw_path, chrom, sta, end):
    bw = pyBigWig.open(bw_path)
    vals = np.array(bw.values(chrom, sta, end, numpy=True))
    bw.close()

    if np.all(np.isnan(vals)):
        return None, 0

    idx = np.nanargmax(vals)
    x_max = sta + idx
    max_v = vals[idx]
    return x_max, max_v

def FancTrackLinePlot(filename, color, label, fill=True, ylabel=None, ylim=None, ax=None):
    track = fanc.load(filename)
    hp = fancplot.LinePlot(track, fill=fill, style='mid', colors=[color], labels=[label], ylabel=ylabel, ylim=ylim, n_yticks=2, ax=ax)
    return hp
    

def FancTrackLinePlot_safe(filename, color, label, fill=True, ylabel=None, ylim=None, ax=None, region=None, log_transform=False):
    track = fanc.load(filename)
    hp = fancplot.LinePlot(track, fill=fill, style='mid', colors=[color], labels=[label],
                           ylabel=ylabel, ylim=ylim, n_yticks=2, ax=ax)
    
    orig_line_values = hp._line_values
    def _line_values_safe(self, plot_region):
        has_valid_data = False
        try:
            for i, x, y in orig_line_values(plot_region):
                if x is None or len(x) == 0:
                    continue
                has_valid_data = True
                
                # 【核心逻辑】：对 Y 值进行 Log2(y + 1) 转换
                if log_transform:
                    y = np.log2(y + 1)
                    
                yield i, x, y
        except Exception as e:
            print(f"Warning: original _line_values failed for {filename}: {e}")
        if not has_valid_data and region is not None:
            try:
                chrom, start_end = region.split(":")
                start, end = map(int, start_end.split("-"))
            except:
                start, end = 0, 1000
                
            x_mock = np.arange(start, end)
            y_mock = np.zeros(end - start)
            yield 0, x_mock, y_mock
    hp._line_values = _line_values_safe.__get__(hp)
    return hp

def get_OEmatrix(chr,sta,end):
    wt_clr = cooler.Cooler("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/Mcool/WT.MicroC.5kb.cool")
    ko_clr = cooler.Cooler("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/HicFiles/Downsampling/coolFiles/KO.contact_map.downsample.KR.5kb.cool")
    region = chr+":"+str(sta)+'-'+str(end)

    sta=int(region.split('-')[0].split(':')[1])
    end=int(region.split('-')[1])
    wt_mat = wt_clr.matrix(balance=True).fetch(region)
    ko_mat = ko_clr.matrix(balance=True).fetch(region)
    wt_exp = pd.read_csv("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/cooltools/WT.MicroC.5kb.expected.smoothed.tsv", sep="\t")
    ko_exp = pd.read_csv("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/HicFiles/Downsampling/KO.contact_map.downsample.KR.5kb.expected.smoothed.tsv", sep="\t")
    chrom = region.split(":")[0]

    def build_expected_array(exp_df):
        max_dist = exp_df['dist'].max()
        exp_arr = np.full(max_dist + 1, np.nan)
        for _, row in exp_df.iterrows():
            d = int(row['dist'])
            exp_arr[d] = row['balanced.avg.smoothed']
        return exp_arr

    wt_exp_chr = wt_exp[wt_exp['region1'] == chrom]
    ko_exp_chr = ko_exp[ko_exp['region1'] == chrom]
    wt_expected = build_expected_array(wt_exp_chr)
    ko_expected = build_expected_array(ko_exp_chr)

    def compute_oe_fast(mat, expected):
        n = mat.shape[0]
        i, j = np.indices((n, n))
        dist = np.abs(i - j)
        dist[dist >= len(expected)] = len(expected) - 1
        exp_vals = expected[dist]
        oe = mat / exp_vals
        oe[(exp_vals <= 0) | np.isnan(exp_vals)] = np.nan
        return oe

    wt_oe = compute_oe_fast(wt_mat, wt_expected)
    ko_oe = compute_oe_fast(ko_mat, ko_expected)
    wt_oe = np.nan_to_num(wt_oe, nan=0)
    ko_oe = np.nan_to_num(ko_oe, nan=0)
    # wt_oe = np.log2(wt_oe)
    # ko_oe = np.log2(ko_oe)
    n = wt_oe.shape[0]

    mask_upper = np.triu(np.ones((n, n)), k=1)
    mask_lower = np.tril(np.ones((n, n)))

    merged = wt_oe * mask_upper + ko_oe * mask_lower
    # colors = [(0, '#EEEEEE'), 
    #                 (0.5,'#E76F2E'),
    #                 (1, '#3E2C23')] 
    colors = [(0, '#EEEEEE'), 
                    (1,'#E76F2E')] 
    colors = [(0, '#EEEEEE'), 
                    (0.5,'#E76F2E'),
                    (1, '#B12C00')] 
    ccmap=LinearSegmentedColormap.from_list('custom_cmap', colors)
    return merged

def plot_TAD(chrom,sta,end):
    style.use('default')
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['axes.unicode_minus'] =False
    plt.rcParams['axes.unicode_minus'] =False
    region = GenomicRegion(chrom, sta, end)

    wt_insul_region = bioframe.select(wt_insulation_table, (region.chromosome,region.start,region.end))
    c_insul_region = bioframe.select(c_insulation_table, (region.chromosome,region.start,region.end))

    fanc_region = f'{region.chromosome}:{(region.start/1000000):.2f}mb-{(region.end/1000000):.2f}mb'
    fig = plt.figure(figsize=(3.5, 4.7))
    gs = gridspec.GridSpec(9, 1,  height_ratios=[3,3]+list(np.ones(7)),hspace=0.1)
    colors1 = [(0, '#FAF3E1'),   
            (1, '#C40C0C')] 
    ccmap=LinearSegmentedColormap.from_list('custom_cmap', colors1)
    ax1=fig.add_subplot(gs[0,0])
    WT_filename = f'{hic_file1}@5kb@KR'
    WT_hic = fanc.load(WT_filename)
    WT_HiC_plot = fancplot.TriangularMatrixPlot(WT_hic, rasterized=True,
                                                vmin=0, vmax=15, 
                                                # oe=True,log=True,vmin=-2, vmax=2,
                                                max_dist='0.3mb', colormap=ccmap, show_colorbar=False, ax=ax1)
    WT_HiC_plot.plot(fanc_region)
    ax1.set_xticks([])
    # plt.title(chrom+':'+str(sta)+'-'+str(end))

    ax2=fig.add_subplot(gs[1,0])
    KD_filename = f'{hic_file2}@5kb@KR'
    KD_hic = fanc.load(KD_filename)
    KD_HiC_plot = fancplot.TriangularMatrixPlot(KD_hic,  rasterized=True,
                                                vmin=0, vmax=15,
                                                # oe=True,log=True,vmin=-2, vmax=2,
                                                max_dist='0.3mb', colormap=ccmap, show_colorbar=False, ax=ax2) #norm="lin", 
    KD_HiC_plot.plot(fanc_region)
    ax2.set_xticks([])
    colors=['#44444E','#740A03']
    ax3=fig.add_subplot(gs[2,0])
    ax3.plot(wt_insul_region[['start', 'end']].mean(axis=1),
                wt_insul_region['log2_insulation_score'],
                color=colors[0],linewidth=0.5)

    ax4=fig.add_subplot(gs[3,0])
    ax4.plot(c_insul_region[['start', 'end']].mean(axis=1),
                c_insul_region['log2_insulation_score'],
                color=colors[1],linewidth=0.5)
    minis=min(c_insul_region['log2_insulation_score'])
    minis=math.floor(minis * 2) / 2
    for i in [ax3,ax4]:
        i.spines['right'].set_visible(False)
        i.spines['top'].set_visible(False)
        i.set_xticks([])
        i.set_ylim(minis,1)
        i.set_yticks([minis,1],[])
        i.spines['bottom'].set_position(('data',0))
        i.text(1.1, 0.95,f"{minis},1",transform=i.transAxes,ha='right',va='top',fontsize=9)

    ax5=fig.add_subplot(gs[4,0])
    pos1 = ax5.get_position()
    pos2 = ax5.get_position()
    ax5.set_position([pos1.x0, pos2.y0, pos1.width, pos2.height])
    filename="/data/zhangjy/DEAF1/ChiPseq_Analysis/Step5_combineRep/Bigwig/BW_RPGC/E250604004_L01_WT_DEAF1_RPGC.bw"
    p = FancTrackLinePlot(filename, colors[0], 'WT DEAF1', True,  ylim=(0,10), ax=ax5)
    p.plot(region)
    ax5.text(1.1, 0.95,"0,10",transform=ax5.transAxes,ha='right',va='top',fontsize=9)

    ax6=fig.add_subplot(gs[5,0], sharex=ax1)
    pos1 = ax1.get_position()
    pos2 = ax6.get_position()
    ax6.set_position([pos1.x0, pos2.y0, pos1.width, pos2.height])
    filename="/data/zhangjy/DEAF1/ChiPseq_Analysis/Step5_combineRep/Bigwig/BW_RPGC/E250609001_L01_DEAF1_DEAF1_RPGC.bw"
    p = FancTrackLinePlot(filename, colors[1], 'DEAF1(DEAF1 mutant)', True,  ylim=(0,10), ax=ax6)
    p.plot(region)
    ax6.text(1.1, 0.95,"0,10",transform=ax6.transAxes,ha='right',va='top',fontsize=9)

    filename1="/data/zhangjy/DEAF1/ChiPseq_Analysis/Step5_combineRep/Bigwig/BW_RPGC/E250604004_L01_WT_CTCF_RPGC.bw"
    x_max11, max_v1 = get_region_max_pos(filename1, chrom, sta, end)
    filename2="/data/zhangjy/DEAF1/ChiPseq_Analysis/Step5_combineRep/Bigwig/BW_RPGC/E250609001_L01_DEAF1_CTCF_RPGC.bw"
    x_max2, max_v2 = get_region_max_pos(filename2, chrom, sta, end)
    ctcflim=math.ceil(max(max_v1,max_v2) / 10) * 10
    ax7=fig.add_subplot(gs[6,0], sharex=ax1)
    pos2 = ax7.get_position()
    ax7.set_position([pos1.x0, pos2.y0, pos1.width, pos2.height])
    p = FancTrackLinePlot(filename1, colors[0], 'WT CTCF', True,  ylim=(0,ctcflim), ax=ax7)
    p.plot(region)
    ax7.text(1.1, 0.95,f"0,{ctcflim}",transform=ax7.transAxes,ha='right',va='top',fontsize=9)

    ax8=fig.add_subplot(gs[7,0], sharex=ax1)
    pos2 = ax8.get_position()
    ax8.set_position([pos1.x0, pos2.y0, pos1.width, pos2.height])
    p = FancTrackLinePlot(filename2, colors[1], 'CTCF(DEAF1 mutant)', True,  ylim=(0,ctcflim), ax=ax8)
    p.plot(region)
    ax8.text(1.1, 0.95,f"0,{ctcflim}",transform=ax8.transAxes,ha='right',va='top',fontsize=9)

    ax9=fig.add_subplot(gs[8,0])
    filename = "/data/zhangjy/Reference/Human/hg38/gencode.v48.annotation.protein_coding.gtf"
    with tempfile.NamedTemporaryFile(mode="w", suffix=".gtf", delete=False) as tmp_gtf:
        with open(filename, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9:
                    continue
                attrs = fields[8]
                if attrs.startswith('"') and attrs.endswith('"'):
                    attrs = attrs[1:-1]
                attrs = attrs.replace('""', '"')
                fields[8] = attrs
                tmp_gtf.write("\t".join(fields) + "\n")
    tmp_gtf_path = tmp_gtf.name
    p = fancplot.GenePlot(tmp_gtf_path, group_by='gene_id', squash=True, label_field='gene_name',
                            color_forward='#999999', color_reverse='#999999', show_arrows=True,
                            arrow_size=4, line_width=1, box_height=1,
                        font_size=7,ax=ax9)
    p.plot(region)
    ax9.set_xticks([sta, end])
    ax9.set_xticklabels([f"{chrom}:{sta/1e6:.3f} Mb",f"{end/1e6:.3f} Mb"])
    ax9.spines[:].set_visible(False)
    ax9.minorticks_off()

    pos_ax3 = ax3.get_position()
    pos_ax4 = ax4.get_position()
    mid_y = (pos_ax3.y0 + pos_ax4.y1) / 2
    label_x = pos1.x0-0.15  
    fig.text(label_x, mid_y, 'log2(IS)', rotation=0,ha='center', va='center', fontsize=10)
    pos_ax5 = ax5.get_position()
    pos_ax6 = ax6.get_position()
    mid_y = (pos_ax5.y0 + pos_ax6.y1) / 2
    fig.text(label_x, mid_y, 'DEAF1\n(RPGC)', rotation=0,ha='center', va='center', fontsize=10)
    pos_ax7 = ax7.get_position()
    pos_ax8 = ax8.get_position()
    mid_y = (pos_ax7.y0 + pos_ax8.y1) / 2
    label_x = pos1.x0-0.15  
    fig.text(label_x, mid_y, 'CTCF\n(RPGC)', rotation=0,ha='center', va='center', fontsize=10)


    li1=[ax5,ax6,ax7,ax8]#
    for i  in li1:
        i.set_xticks([])
        i.set_yticks([])
        i.set_xlim(sta,end)
        i.spines[:].set_visible(False)
        for line in i.lines:
            line.set_rasterized(True)
        for col in i.collections:
            col.set_rasterized(True)
    plt.savefig(f"/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/Domains/Cooltools/Bedtools/Pdf/"+chrom+":"+str(sta)+'-'+str(end)+'.TADboundariesIS.pdf',   ##
            bbox_inches = 'tight',dpi=150,
            facecolor='w')  


plot_TAD("chr7",71540000,71950000)
plot_TAD("chr1",13740000,14050000)
