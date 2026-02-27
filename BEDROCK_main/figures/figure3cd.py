#modified bases by AT or CG in 1kb
import pandas as pd
import numpy as np
import pyranges as pr
import seaborn as sns
import matplotlib.pyplot as plt

from utils import make_palette

WINDOW_SIZE = 1000

def get_chromosome_sizes(bed_df):

    chrom_sizes = (
        bed_df
        .groupby("Chromosome")["End_chrom_pos"]
        .max()
        .reset_index()
        .rename(columns={"End_chrom_pos": "chrom_length"})
    )
    
    return chrom_sizes

def make_1kb_windows(chrom_sizes, window_size = 1000):
    windows = []
    for _, row in chrom_sizes.iterrows():
        chrom = row["Chromosome"]
        length = row["chrom_length"]

        starts = np.arrange(0, length, window_size)
        ends = starts + window_size
        ends[ends > length] = length

        df = pd.DataFrame({
            "Chromosome": chrom,
            "Start": starts,
            "End": ends
        })

        windows.append(df)

    windows_df = pd.concat(windows, ignore_index=True)

    return pr.Pyranges(windows_df)

def bed_to_ranges(df_bed):
    df = df_bed.rename(columns={
        "Start_chrom_pos": "Start",
        "End_chrom_pos": "End"
    })
    
    return pr.PyRanges(df)


def summarize_modifications(bed_df, window_pr, mod_codes):

    filtered = bed_df[(bed_df["mod_code"].isin(mod_codes)) & (bed_df["mod"] == 1)].copy()
    bed = bed_to_ranges(filtered)
    overlaps = bed.join(window_pr)
    df = overlaps.df

    summary = (
        df.groupby(
            ["Chromosome", "Start_b", "End_b", "strand", "sample"]
        )
        .size()
        .reset_index(name="mod_count")
        .rename(columns={
            "Start_b": "window_start",
            "End_b": "window_end"
        })
    )
    
    all_windows = window_pr.df.rename(
        columns={"Start": "window_start", "End": "window_end"}
        )

    samples = bed_df["sample"].unique()
    strands = ["+", "-"]

    template = (
        all_windows
        .assign(key=1)
        .merge(
            pd.DataFrame({
                "sample": samples,
                "strand": strands,
                "key": 1
            }),
            on="key"
        )
        .drop("key", axis=1)
    )

    mod_count_summary_full = template.merge(
        summary,
        on=["Chromosome", "window_start", "window_end", "strand", "sample"],
        how="left"
    )
    mod_count_summary_full["mod_count"] = mod_count_summary_full["mod_count"].fillna(0)

    mod_count_summary_full.loc[mod_count_summary_full["strand"] == "-", "mod_count"] *= -1

    
    return mod_count_summary_full
    
    
def plot_mod_windows(df, outpath, ylab):

    print(df.columns)
    print(df["mod_count"].abs().max())

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
                sub["window_start"],
                sub["mod_count"],
                width=WINDOW_SIZE,
                align="edge",
                color=color
            )

        ax.axhline(0, color="black", linewidth=0.8)
        ax.set_xticks(ax.get_xticks())
        ax.set_xticklabels((ax.get_xticks() / 1000).astype(int))

    g.map_dataframe(draw)

    g.set_axis_labels("Genomic Position (kb)", ylab)
    g.set_titles(col_template="{col_name}", row_template="")

    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()
