#modified bases by AT or CG in 1kb
import pandas as pd
import pyranges as pr
import seaborn as sns
import matplotlib.pyplot as plt

from utils import make_palette

WINDOW_SIZE = 1000

def overlap_with_windows(df_bedmethyl, df_windows):

    win = df_windows.rename(columns={"Start": "Start", "End": "End"})
    bed = df_bedmethyl.rename(columns={
        "Start_chrom_pos": "Start",
        "End_chrom_pos": "End",
        "sample" : "sample"
    })
    
    win_pr = pr.PyRanges(win)
    bed_pr = pr.PyRanges(bed)
    overlaps = bed_pr.join(win_pr)

    return overlaps.df

def summarize_cg(overlap_df):

    df = overlap_df.copy()

    df = df[
        df["mod_code"].isin(["h", "m"]) &
        (df["mod"] == 1)
    ]

    summary = (
        df.groupby(
            ["Chromosome", "Start", "End", "strand", "sample"],
            as_index=False
        )
        .size()
        .rename(columns={"size": "mod"})
    )

    summary.loc[summary["strand"] == "-", "mod"] *= -1

    return summary

def summarize_at(overlap_df):

    df = overlap_df.copy()

    df = df[
        df["mod_code"].isin(["a"]) &
        (df["mod"] == 1)
    ]

    summary = (
        df.groupby(
            ["Chromosome", "Start", "End", "strand", "sample"],
            as_index=False
        )
        .size()
        .rename(columns={"size": "mod"})
    )

    summary.loc[summary["strand"] == "-", "mod"] *= -1

    return summary

def plot_windows(df, outpath, ylab):

    df = df.copy()

    chrom_order = list(dict.fromkeys(df["Chromosome"]))
    df["Chromosome"] = pd.Categorical(
        df["Chromosome"],
        categories=chrom_order,
        ordered=True
    )

    g = sns.FacetGrid(
        df,
        row="sample",
        col="Chromosome",
        sharex=False,
        sharey=True,
        height=3
    )

def draw(data, **kwargs):
    ax = plt.gca()

    for strand, sub in data.groupby("strand"):

    color = "red" if strand == "+" else "blue"

    ax.bar(
        sub["Start"],
        sub["mod"],
        width=WINDOW_SIZE,
        align="edge",
        color=color
    )

    ax.axhline(0, color="black", linewidth=0.8)
    ax.tick_params(axis="x", labelbottom=False)

g.map_dataframe(draw)

g.set_axis_labels("Genomic Position (kb)", ylab)
g.set_titles(col_template="{col_name}", row_template="")

plt.tight_layout()
plt.savefig(outpath, dpi=300)
plt.close()
