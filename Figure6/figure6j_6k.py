import math
import tempfile
import cooler
import fanc
import fanc.plotting as fancplot
import matplotlib.font_manager as font_manager
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig
from matplotlib import style
from matplotlib.colors import LinearSegmentedColormap
import os

w5=pd.read_table("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/LoopStrength/WT.5kb.rescored_strength_coords.txt")
w10=pd.read_table("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/LoopStrength/WT.10kb.rescored_strength_coords.txt")
k5=pd.read_table("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/LoopStrength/KO.5kb.rescored_strength_coords.txt")
k10=pd.read_table("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/LoopStrength/KO.10kb.rescored_strength_coords.txt")

ws=pd.concat([w5,w10]).sort_values(['chrom1','start1'])
ks=pd.concat([k5,k10]).sort_values(['chrom1','start1'])
ws_tmp=ws.iloc[:,[0,1,2,3,4,5,10]]
ws_tmp.columns=[0,1,2,3,4,5,'strengthWT']
ks_tmp=ks.iloc[:,[0,1,2,3,4,5,10]]
ks_tmp.columns=[0,1,2,3,4,5,'strengthKO']
loopstrength=pd.merge(ws_tmp,ks_tmp,on=[0,1,2,3,4,5])

def plot_loops_track(chr_name, sta, end, loop_df, ax_wt,ax_ko):
    valid_loops = loop_df[(loop_df.iloc[:, 0] == chr_name) & (loop_df.iloc[:, 1] >= sta) & (loop_df.iloc[:, 4] <= end)]
    max_strength_in_window = max(valid_loops.iloc[:, 6].max(), valid_loops.iloc[:, 7].max())
    if max_strength_in_window == 0:
            max_strength_in_window = 1.0
    MIN_LW = 0.3
    MAX_LW = 3

    for _, row in valid_loops.iterrows():
        anchor1_mid = (row.iloc[1] + row.iloc[2]) / 2
        anchor2_mid = (row.iloc[4] + row.iloc[5]) / 2
        
        wt_strength = row.iloc[6]
        ko_strength = row.iloc[7]
        
        lw_wt = MIN_LW + (MAX_LW - MIN_LW) * (wt_strength / max_strength_in_window)
        lw_ko = MIN_LW + (MAX_LW - MIN_LW) * (ko_strength / max_strength_in_window)
        
        loop_span = anchor2_mid - anchor1_mid
        arc_height = loop_span * 0.15 
        
        path_data_wt = [
            (Path.MOVETO, (anchor1_mid, 0)),
            (Path.CURVE3, ( (anchor1_mid + anchor2_mid)/2, arc_height )),
            (Path.CURVE3, (anchor2_mid, 0))
        ]
        codes_wt, verts_wt = zip(*path_data_wt)
        path_wt = Path(verts_wt, codes_wt)
        patch_wt = patches.PathPatch(path_wt, edgecolor='#7F8487', facecolor='none', 
                                    lw=lw_wt, alpha=0.8, zorder=3)
        ax_wt.add_patch(patch_wt)

        path_data_ko = [
            (Path.MOVETO, (anchor1_mid, 0)),
            (Path.CURVE3, ( (anchor1_mid + anchor2_mid)/2, arc_height )), # 改为正值
            (Path.CURVE3, (anchor2_mid, 0))
        ]
        codes_ko, verts_ko = zip(*path_data_ko)
        path_ko = Path(verts_ko, codes_ko)
        patch_ko = patches.PathPatch(path_ko, edgecolor='#8E1616', facecolor='none', 
                                    lw=lw_ko, alpha=0.8, zorder=3)
        ax_ko.add_patch(patch_ko)

    max_span = (end - sta) * 0.18  
    sta_mb_str = chr_name+":"+f"{sta / 1e6:.3f} Mb"
    end_mb_str = f"{end / 1e6:.3f} Mb"
    for ax in [ax_wt, ax_ko]:
        ax.set_xlim(sta, end+1)
        ax.set_ylim(0, max_span) 
        ax.axis('off')
        ax.hlines(y=0, xmin=sta, xmax=end, colors='#2C3947', linewidth=0.8, zorder=2)
        tick_h = max_span * 0.1 # 动态计算刻度线高度
        ax.vlines(x=[sta, end], ymin=0, ymax=tick_h, colors='#2C3947', linewidth=0.8, zorder=2)
    ax_wt.text(x=sta-(end-sta)*0.1, y=-tick_h*1.2, s=sta_mb_str, ha='left', va='top', fontsize=8, color='#555555')
    ax_wt.text(x=end+(end-sta)*0.1, y=-tick_h*1.2, s=end_mb_str, ha='right', va='top', fontsize=8, color='#555555')
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

def Gene_anchor(region,gene,path):
    sta=int(region.split(":")[1].split('-')[0])
    end=int(region.split(":")[1].split('-')[1])
    distancex=(end-sta)/16
    chrom=region.split(":")[0]
    colors=['#44444E','#740A03']
    style.use('default')
    font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
    plt.rcParams['font.sans-serif']='Helvetica'
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['axes.unicode_minus'] =False
    fig = plt.figure(figsize=(3.5, 4.4))
    fig.canvas.draw()
    gs = gridspec.GridSpec(13, 1,  height_ratios=[2,2]+list(np.ones(11)),hspace=0.1)
    loop_df=loopstrength
    ax1 = fig.add_subplot(gs[0, 0])  # 第一行：WT
    ax2 = fig.add_subplot(gs[1, 0])  # 第二行：KO/DEAF1-mutant
    plot_loops_track(chrom, sta, end, loop_df, ax1, ax2)

    ax3=fig.add_subplot(gs[2,0])
    pos1 = ax3.get_position()
    pos2 = ax3.get_position()
    ax3.set_position([pos1.x0, pos2.y0, pos1.width, pos2.height])
    filename="/data/zhangjy/DEAF1/ChiPseq_Analysis/Step5_combineRep/Bigwig/BW_RPGC/E250604004_L01_WT_DEAF1_RPGC.bw"
    p = FancTrackLinePlot(filename, colors[0], 'WT DEAF1', True,  ylim=(0,10), ax=ax3)
    p.plot(region)
    ax3.text(1.1, 0.95,"0,10",transform=ax3.transAxes,ha='right',va='top',fontsize=9)

    ax4=fig.add_subplot(gs[3,0], sharex=ax1)
    pos1 = ax1.get_position()
    pos2 = ax4.get_position()
    ax4.set_position([pos1.x0, pos2.y0, pos1.width, pos2.height])
    filename="/data/zhangjy/DEAF1/ChiPseq_Analysis/Step5_combineRep/Bigwig/BW_RPGC/E250609001_L01_DEAF1_DEAF1_RPGC.bw"
    p = FancTrackLinePlot(filename, colors[1], 'DEAF1(DEAF1 mutant)', True,  ylim=(0,10), ax=ax4)
    p.plot(region)
    ax4.text(1.1, 0.95,"0,10",transform=ax4.transAxes,ha='right',va='top',fontsize=9)

    filename1="/data/zhangjy/DEAF1/ChiPseq_Analysis/Step5_combineRep/Bigwig/BW_RPGC/E250604004_L01_WT_CTCF_RPGC.bw"
    x_max11, max_v1 = get_region_max_pos(filename1, chrom, sta, end)
    filename2="/data/zhangjy/DEAF1/ChiPseq_Analysis/Step5_combineRep/Bigwig/BW_RPGC/E250609001_L01_DEAF1_CTCF_RPGC.bw"
    x_max2, max_v2 = get_region_max_pos(filename2, chrom, sta, end)
    ctcflim=math.ceil(max(max_v1,max_v2) / 10) * 10
    ax5=fig.add_subplot(gs[4,0], sharex=ax1)
    pos2 = ax5.get_position()
    ax5.set_position([pos1.x0, pos2.y0, pos1.width, pos2.height])
    p = FancTrackLinePlot(filename1, colors[0], 'WT CTCF', True,  ylim=(0,ctcflim), ax=ax5)
    p.plot(region)
    ax5.text(1.1, 0.95,f"0,{ctcflim}",transform=ax5.transAxes,ha='right',va='top',fontsize=9)

    ax6=fig.add_subplot(gs[5,0], sharex=ax1)
    pos2 = ax6.get_position()
    ax6.set_position([pos1.x0, pos2.y0, pos1.width, pos2.height])
    p = FancTrackLinePlot(filename2, colors[1], 'CTCF(DEAF1 mutant)', True,  ylim=(0,ctcflim), ax=ax6)
    p.plot(region)
    ax6.text(1.1, 0.95,f"0,{ctcflim}",transform=ax6.transAxes,ha='right',va='top',fontsize=9)

    ax7=fig.add_subplot(gs[6,0])
    pos1 = ax7.get_position()
    pos2 = ax7.get_position()
    ax7.set_position([pos1.x0, pos2.y0, pos1.width, pos2.height])
    filename1="/data/zhangjy/DEAF1/Histone_Analysis/Step4_CombinedRep/BW_combined/RPGC/WT-H3K27ac.RPGC.bw"
    filename2="/data/zhangjy/DEAF1/Histone_Analysis/Step4_CombinedRep/BW_combined/RPGC/KO-H3K4me1.RPGC.bw"
    x_max1, max_v1 = get_region_max_pos(filename1, chrom, sta, end)
    x_max2, max_v2 = get_region_max_pos(filename2, chrom, sta, end)
    h3k27aclim=math.ceil(max(max_v1,max_v2) / 10) * 10
    p = FancTrackLinePlot_safe(filename1, colors[0], 'WT', True,  ylim=(0,h3k27aclim), ax=ax7)
    p.plot(region)
    ax7.text(1.1, 0.95,f"0,{h3k27aclim}",transform=ax7.transAxes,ha='right',va='top',fontsize=9)

    ax8=fig.add_subplot(gs[7,0])
    pos1 = ax1.get_position()
    pos2 = ax8.get_position()
    ax8.set_position([pos1.x0, pos2.y0, pos1.width, pos2.height])
    p = FancTrackLinePlot_safe(filename2, colors[1], 'WT', True,  ylim=(0,h3k27aclim), ax=ax8)
    p.plot(region)
    ax8.text(1.1, 0.95,f"0,{h3k27aclim}",transform=ax8.transAxes,ha='right',va='top',fontsize=9)

    ax9=fig.add_subplot(gs[8,0])
    pos1 = ax9.get_position()
    pos2 = ax9.get_position()
    ax9.set_position([pos1.x0, pos2.y0, pos1.width, pos2.height])
    filename1="/data/zhangjy/DEAF1/Histone_Analysis/Step4_CombinedRep/BW_combined/RPGC/WT-H3K4me3.RPGC.bw"
    filename2="/data/zhangjy/DEAF1/Histone_Analysis/Step4_CombinedRep/BW_combined/RPGC/KO-H3K4me3.RPGC.bw"
    x_max1, max_v1 = get_region_max_pos(filename1, chrom, sta, end)
    x_max2, max_v2 = get_region_max_pos(filename2, chrom, sta, end)
    h3k27aclim=math.ceil(max(max_v1,max_v2) / 10) * 10
    p = FancTrackLinePlot_safe(filename1, colors[0], 'WT', True,  ylim=(0,h3k27aclim), ax=ax9)
    p.plot(region)
    ax9.text(1.1, 0.95,f"0,{h3k27aclim}",transform=ax9.transAxes,ha='right',va='top',fontsize=9)

    ax10=fig.add_subplot(gs[9,0])
    pos1 = ax1.get_position()
    pos2 = ax10.get_position()
    ax10.set_position([pos1.x0, pos2.y0, pos1.width, pos2.height])
    p = FancTrackLinePlot_safe(filename2, colors[1], 'WT', True,  ylim=(0,h3k27aclim), ax=ax10)
    p.plot(region)
    ax10.text(1.1, 0.95,f"0,{h3k27aclim}",transform=ax10.transAxes,ha='right',va='top',fontsize=9)

    ax11=fig.add_subplot(gs[10,0])
    pos1 = ax11.get_position()
    pos2 = ax11.get_position()
    ax11.set_position([pos1.x0, pos2.y0, pos1.width, pos2.height])
    filename1="/data/zhangjy/DEAF1/RNAseq_Analysis/Step4_BW/Combined/wt_RPKM.bw"
    filename2="/data/zhangjy/DEAF1/RNAseq_Analysis/Step4_BW/Combined/ko_RPKM.bw"
    x_max1, max_v1 = get_region_max_pos(filename1, chrom, sta, end)
    x_max2, max_v2 = get_region_max_pos(filename2, chrom, sta, end)
    # rnaseqlim=math.ceil(max(max_v1,max_v2) / 10) * 10
    rnaseqlim=200
    p = FancTrackLinePlot_safe(filename1, colors[0], 'WT', True,  ylim=(0,rnaseqlim), ax=ax11)
    p.plot(region)
    ax11.text(1.1, 0.95,f"0,{rnaseqlim}",transform=ax11.transAxes,ha='right',va='top',fontsize=9)

    ax12=fig.add_subplot(gs[11,0])
    pos1 = ax1.get_position()
    pos2 = ax12.get_position()
    ax12.set_position([pos1.x0, pos2.y0, pos1.width, pos2.height])
    p = FancTrackLinePlot_safe(filename2, colors[1], 'WT', True,  ylim=(0,rnaseqlim), ax=ax12)
    p.plot(region)
    ax12.text(1.1, 0.95,f"0,{rnaseqlim}",transform=ax12.transAxes,ha='right',va='top',fontsize=9)

    ax13=fig.add_subplot(gs[12,0])
    filename = "/data/zhangjy/Reference/Human/hg38/gencode.v48.annotation.gtf"
    target_genes = [gene]  
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

                if any(f'gene_name "{gene}";' in attrs for gene in target_genes):
                    fields[8] = attrs
                    tmp_gtf.write("\t".join(fields) + "\n")
    tmp_gtf_path = tmp_gtf.name
    p = fancplot.GenePlot(tmp_gtf_path, group_by='gene_id', squash=True, label_field='gene_name',
                            color_forward='#999999', color_reverse='#999999', show_arrows=True,
                            arrow_size=4, line_width=1, box_height=1,
                        font_size=7,ax=ax13)
    p.plot(region)
    ax13.set_xticks([sta, end])
    ax13.set_xticklabels([f"{chrom}:{sta/1e6:.3f} Mb",f"{end/1e6:.3f} Mb"])
    ax13.spines[:].set_visible(False)
    ax13.text(x_max11, -0.1, x_max11,transform=ax13.get_xaxis_transform(),ha='center', va='top',fontsize=8, color=colors[0])
    ax13.minorticks_off()

    pos_ax3 = ax3.get_position()
    pos_ax4 = ax4.get_position()
    mid_y = (pos_ax3.y0 + pos_ax4.y1) / 2
    label_x = pos1.x0-0.15  
    fig.text(label_x, mid_y, 'DEAF1\n(RPGC)', rotation=0,ha='center', va='center', fontsize=10, fontweight='bold')
    pos_ax5 = ax5.get_position()
    pos_ax6 = ax6.get_position()
    mid_y = (pos_ax5.y0 + pos_ax6.y1) / 2
    fig.text(label_x, mid_y, 'CTCF\n(RPGC)', rotation=0,ha='center', va='center', fontsize=10, fontweight='bold')
    pos_ax7 = ax7.get_position()
    pos_ax9 = ax9.get_position()
    mid_y = (pos_ax7.y0 + pos_ax9.y1) / 2
    fig.text(label_x, mid_y, 'H3K4me3\n(RPGC)', rotation=0,ha='center', va='center', fontsize=10, fontweight='bold')
    pos_ax9 = ax9.get_position()
    pos_ax10 = ax10.get_position()
    mid_y = (pos_ax9.y0 + pos_ax10.y1) / 2
    fig.text(label_x, mid_y, 'H3K27ac\n(RPGC)', rotation=0,ha='center', va='center', fontsize=10, fontweight='bold')
    pos_ax11 = ax11.get_position()
    pos_ax12 = ax12.get_position()
    mid_y = (pos_ax11.y0 + pos_ax12.y1) / 2
    fig.text(label_x, mid_y, 'RNA-seq\n(RPGC)', rotation=0,ha='center', va='center', fontsize=10, fontweight='bold')

    li1=[ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10,ax11,ax12]#
    for i  in li1:
        i.set_xticks([])
        i.set_yticks([])
        i.set_xlim(sta,end)
        i.spines[:].set_visible(False)
        for line in i.lines:
            line.set_rasterized(True)
        for col in i.collections:
            col.set_rasterized(True)
    plt.savefig(f"/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/mustache/Downsampling/Pdf/{path}/"+gene+"_"+chrom+":"+str(sta)+'-'+str(end)+'.plusRNAseq.pdf',   ##
            bbox_inches = 'tight',dpi=150,
            facecolor='w')  

Gene_anchor("chr5:178200000-178900000","COL23A1","EPexample")
Gene_anchor("chr17:68750000-68965000","ABCA8","EPexample")