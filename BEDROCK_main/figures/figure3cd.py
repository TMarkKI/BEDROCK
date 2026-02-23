#modified bases by AT or CG in 1kb
import pandas as pd
import pyranges as pr
import seaborn as sns
import matplotlib.pyplot as plt

from utils import make_palette

WINDOW_SIZE = 1000


def make_windows(df_windows):
    
    assert {"Chromosome", "start", "end"}.issubset(df_windows.columns)
    
    return pr.PyRanges(
        df_windows.rename(
            columns={"start": "Start", "end": "End"}
        )
        
    )


def bed_to_ranges(df):
    return pr.PyRanges(
        df.rename(
            columns={
                "Chromosome": "Chromosome",
                "Start_chrom_pos": "Start",
                "End_chrom_pos": "End",
            }
        )
    )


def summarize_modifications(
    bed_df,
    window_df,
    mod_codes,
):
    bed = bed_to_ranges(bed_df)
    windows = make_windows(window_df)

    overlaps = bed.join(windows)

    df = overlaps.df
    df = df[df["mod_code"].isin(mod_codes) & (df["mod"] == 1)]

    summary = (
        df.groupby(
            ["Chromosome", "Start", "End", "strand", "sample_name"]
        )
        .size()
        .reset_index(name="mod")
    )

    summary.loc[summary.strand == "-", "mod"] *= -1

    return summary


def plot_mod_windows(df, outpath, ylab):
    
    df["strand"] = pd.Categorical(
        df["strand"],
        categories=["+", "-"],
        ordered=True
    )

    chrom_order = list(dict.fromkeys(df["Chromosome"]))

    df["Chromosome"] = pd.Categorical(
        df["Chromosome"],
        categories=chrom_order,
        ordered=True
    )

    palette = make_palette(df["strand"].unique(), palette="Set1")

    g = sns.FacetGrid(
        df,
        row="sample_name",
        col="Chromosome",
        sharex=False,
        height=3,
    )

    g.map_dataframe(
        sns.barplot,
        x="Start",
        y="mod",
        hue="strand",
        palette=palette,
        dodge=False,
    )

    g.add_legend()
    handles, labels = g.axes[0, 0].get_legend_handles_labels()
    g.legend.remove()

    plt.tight_layout()
    
    g.fig.legend(
        handles=handles,
        labels=labels,
        title="Strand",
        bbox_to_anchor=(1.02, 0.5),
        loc="center left",
        borderaxespad=0
    )
    
    for ax in g.axes.flatten():
        ax.tick_params(axis="x", labelbottom=False)

    plt.savefig(outpath, dpi=300)
    plt.close()
