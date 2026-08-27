import matplotlib.font_manager as font_manager
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig
from matplotlib import style

boundary_bed = "/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/Domains/Cooltools/InsulationScores/WT_KO.TAD_boundaries_merged.10kb.bed"
ctcf_peak_bed = "/data/zhangjy/DEAF1/ChiPseq_Analysis/Step5_combineRep/Narrowpeaks_combineDEAF1input/merged_CTCFpeaks_signal_classified.bed"
wt_ins_bw = "/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/Domains/Cooltools/WT.insulation_10kb.bw"
ko_ins_bw = "/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/Domains/Cooltools/KO.insulation_10kb.bw"
out_prefix = "/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/Domains/Cooltools/InsulationScores/TAD_boundary_CTCF_peak_insulation_correlation"


def mean_bw_signal(bw, chrom, start, end):
    if chrom not in bw.chroms():
        return np.nan
    chrom_len = bw.chroms()[chrom]
    start = max(0, int(start))
    end = min(int(end), chrom_len)
    if start >= end:
        return np.nan
    val = bw.stats(chrom, start, end, type="mean")[0]
    return np.nan if val is None else val
def add_insulation_signal(boundary_bed):
    boundaries = pd.read_csv(
        boundary_bed,
        sep="\t",
        header=None,
        names=["chrom", "start", "end"]
    )

    wt_bw = pyBigWig.open(wt_ins_bw)
    ko_bw = pyBigWig.open(ko_ins_bw)

    records = []

    for _, row in boundaries.iterrows():
        chrom = row["chrom"]
        start = int(row["start"])
        end = int(row["end"])
        center = int((start + end) / 2)

        wt_ins = mean_bw_signal(wt_bw, chrom, start, end)
        ko_ins = mean_bw_signal(ko_bw, chrom, start, end)

        records.append({
            "chrom": chrom,
            "boundary_start": start,
            "boundary_end": end,
            "center": center,
            "WT_insulation": wt_ins,
            "KO_insulation": ko_ins,
        })

    wt_bw.close()
    ko_bw.close()

    df = pd.DataFrame(records)

    df["insulation_change"] = df["KO_insulation"] - df["WT_insulation"]
    df["boundary_insulation_gain"] = df["WT_insulation"] - df["KO_insulation"]

    return df

def add_ctcf_peak_signal(boundary_df, ctcf_peak_bed):
    peaks = pd.read_csv(
        ctcf_peak_bed,
        sep="\t",
        header=None,
        names=["chrom", "peak_start", "peak_end", "WT_CTCF_peak_signal", "KO_CTCF_peak_signal", "signal_ratio", "class"]
    )

    records = []

    for _, b in boundary_df.iterrows():
        chrom = b["chrom"]
        start = int(b["boundary_start"])
        end = int(b["boundary_end"])

        sub = peaks[
            (peaks["chrom"] == chrom) &
            (peaks["peak_end"] > start) &
            (peaks["peak_start"] < end)
        ]

        rec = b.to_dict()

        rec["CTCF_peak_n"] = len(sub)

        if len(sub) > 0:
            rec["WT_CTCF_peak_mean"] = sub["WT_CTCF_peak_signal"].mean()
            rec["KO_CTCF_peak_mean"] = sub["KO_CTCF_peak_signal"].mean()
            rec["WT_CTCF_peak_max"] = sub["WT_CTCF_peak_signal"].max()
            rec["KO_CTCF_peak_max"] = sub["KO_CTCF_peak_signal"].max()
        else:
            rec["WT_CTCF_peak_mean"] = np.nan
            rec["KO_CTCF_peak_mean"] = np.nan
            rec["WT_CTCF_peak_max"] = np.nan
            rec["KO_CTCF_peak_max"] = np.nan

        records.append(rec)

    df = pd.DataFrame(records)

    pseudocount = 0.01

    df["CTCF_peak_delta_mean"] = df["KO_CTCF_peak_mean"] - df["WT_CTCF_peak_mean"]
    df["CTCF_peak_log2FC_mean"] = np.log2(
        (df["KO_CTCF_peak_mean"] + pseudocount) /
        (df["WT_CTCF_peak_mean"] + pseudocount)
    )

    df["CTCF_peak_delta_max"] = df["KO_CTCF_peak_max"] - df["WT_CTCF_peak_max"]
    df["CTCF_peak_log2FC_max"] = np.log2(
        (df["KO_CTCF_peak_max"] + pseudocount) /
        (df["WT_CTCF_peak_max"] + pseudocount)
    )

    df = df.replace([np.inf, -np.inf], np.nan)

    return df
boundary_df = add_insulation_signal(boundary_bed)
df = add_ctcf_peak_signal(boundary_df, ctcf_peak_bed)
df2=df[df['CTCF_peak_n']!=0]
df2=df2.fillna(0)
ma=pd.read_table("/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/Domains/Cooltools/InsulationScores/heatmap_TADboundaries.CTCF_DEAF1.K562.GRCh38.gz",header=None,skiprows=1)
tmp1=df2[df2['insulation_change']>0.25].rename(columns={'chrom':0,'boundary_start':1,'boundary_end':2})
tmp2=df2[(df2['insulation_change']<=0.25)&(df2['insulation_change']>=-0.25)].rename(columns={'chrom':0,'boundary_start':1,'boundary_end':2})
tmp3=df2[df2['insulation_change']<-0.25].rename(columns={'chrom':0,'boundary_start':1,'boundary_end':2})
ma1=pd.merge(tmp1[[0,1,2]],ma,on=[0,1,2])
ma2=pd.merge(tmp2[[0,1,2]],ma,on=[0,1,2])
ma3=pd.merge(tmp3[[0,1,2]],ma,on=[0,1,2])
style.use('default')
font_manager.fontManager.addfont("/data/zhangjy/anaconda3/envs/vscodeJuypter/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/HELVETICA.ttf")
fig = plt.figure(figsize=(1.5,4.5))
gs = gridspec.GridSpec(3, 1,height_ratios=[1,1,1],hspace=0.2)
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['axes.unicode_minus'] =False
plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['pdf.fonttype'] = 42
colors = ['#7F8487','#8E1616']
li=[ma1,ma2,ma3]
li2=[f'log2(IS FC)>0.25\nn={ma1.shape[0]}',f'|log2(IS FC)| ≤ 0.25\nn = {ma2.shape[0]}',f'log2(IS FC)<-0.25\nn={ma3.shape[0]}']
for i in range(3):
    ax1 = fig.add_subplot(gs[i,0])
    plt.plot(list(range(1,51)),li[i].iloc[:,6+25:106-25].mean(),color=colors[0], linewidth = "1.7", alpha=0.85, zorder=3)
    plt.plot(list(range(1,51)),li[i].iloc[:,106+25:206-25].mean(),color=colors[1], linewidth = "1.7", alpha=0.90, zorder=2)
    # plt.yticks([0,-0.5])
    plt.ylim(0.4,2.2)
    if i<2:
        plt.xticks([1,25.5,50],[])
    else:
        plt.xticks([1,25.5,50],['-250','0','250'])
        plt.xlabel('Distance to TAD boundaries(Kb)')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['left'].set_position(('outward', 5))
    ax1.spines['bottom'].set_position(('outward', 2))
    ax1.text(1.4, 0.95,li2[i],transform=ax1.transAxes,ha='right',va='top',fontsize=9)
fig.text(-0.1, 0.5, 'log2(CTCF FC)', rotation=90,ha='center', va='center', fontsize=10)
plt.savefig('/data/zhangjy/DEAF1/MicroC_Analysis/Combined_Analysis/FinalCombined/Domains/Cooltools/InsulationScores/Pdf/TADboundaries_IS.CTCFsignal.pdf',   ##
                bbox_inches = 'tight',
                facecolor='w')  
