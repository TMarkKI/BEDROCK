#modified bases by AT or CG in 1kb
import pandas as pd
import pyranges as pr
import seaborn as sns
import matplotlib.pyplot as plt

from utils import make_palette

WINDOW_SIZE = 1000

def make_windows(df_windows):
    return pr.PyRanges(
        df_windows.rename(columns={
            "Start": "Start",
            "End": "End"
        })
    )

def bed_to_ranges(df_bed):
    return pr.PyRanges(
        df_bed.rename(columns={
            "Start_chrom_pos": "Start",
            "End_chrom_pos": "End"
        })
    )

def summarize_modifications(bed_df, window_df, mod_codes):
    bed = bed_to_ranges(bed_df)
    windows = make_windows(window_df)

    overlaps = bed.join(windows)
    df = overlaps.df

    df = df[
        df["mod_code"].isin(mod_codes) &
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

def plot_mod_windows(df, outpath, ylab):

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
