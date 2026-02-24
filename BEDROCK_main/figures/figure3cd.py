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

    df = df.copy()
    df["signed_mod"] = df["Modified Base Count"]
    df.loc[df["strand"] == "-", "signed_mod"] *= -1

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
        x="Genomic Position(kb)",
        y="Modified Base Count",
        hue="strand",
        palette=palette,
        dodge=False,
    )

    g.set_titles(col_template="{col_name}", row_template="")
    
    for ax, sample in zip(g.axes[:, -1], g.row_names):
        ax.text(
            1.02, 0.5,
            sample,
            transform=ax.transAxes,
            rotation=-90,
            va="center",
            ha="left",
            fontsize=10
        )

    g.fig.legend(
        handles=handles,
        labels=labels,
        title="Strand",
        bbox_to_anchor=(1.02, 0.5),
        loc="center left",
        borderaxespad=0
    )
    
    for ax in g.axes.flatten():
        ax.axhline(0, color="black", linewidth=0.8)
        ax.set_ylabel(ylab)
        ax.text(
            0.99, 0.75, "+",
            transform=ax.transAxes,
            ha="right",
            va="center",
            fontsize=12,
            alpha=0.7
        )
        ax.text(
            0.99, 0.25, "−",
            transform=ax.transAxes,
            ha="right",
            va="center",
            fontsize=12,
            alpha=0.7
        )

    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()
