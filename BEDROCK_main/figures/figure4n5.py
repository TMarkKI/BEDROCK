##genes and methylation
import pyranges as pr
import seaborn as sns
import matplotlib.pyplot as plt


def assign_genes(df_bed, genes_pr):

    mods = pr.PyRanges(
        df_bed.rename(
            columns={
                "Start_chrom_pos": "Start",
                "End_chrom_pos": "End",
            }
        )
    )

    overlaps = mods.join(genes_pr)

    return overlaps.df.rename(columns={"Name": "gene"})

def summarise_gene_methylation(
    df,
    mod_codes=None,
    stat="mean",
):
    if mod_codes is not None:
        df = df[df["mod_code"].isin(mod_codes)]

    agg = {"mean": "mean", "median": "median"}[stat]

    summary = (
        df.groupby(["gene", "sample_name"])["percent_mod"]
        .agg(agg)
        .reset_index(name=f"{stat}_methylation")
    )

    return summary

def to_heatmap_matrix(df, value_col):
    return (
        df.pivot(index="gene", columns="sample_name", values=value_col)
        .fillna(0)
    )

def plot_heatmap(
    matrix,
    outfile,
    scale=None,
    cmap="vlag",
    figsize=(12, 20),
    vmin=None,
    vmax=None,
    title=None,
):
    data = matrix.copy()

    if scale == "row":
        data = (data - data.mean(axis=1).values[:, None]) / (
            data.std(axis=1).values[:, None] + 1e-6
        )

    sns.clustermap(
        data,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        metric="euclidean",
        method="complete",
        yticklabels=False,
    )
    plt.savefig(outfile, dpi=300)
    plt.close()

def run_figure4(df_bed, genes_pr, outdir):
    df_bed = df_bed[df_bed["mod_score"] >= 50] #threshold for depth
    df = assign_genes(df_bed, genes_pr)
    df = df.dropna(subset=["gene"])

    configs = [
        ("4a", None, "mean", None),
        ("4b", None, "median", "row"),
        ("4c", ["a"], "mean", None),
        ("4d", ["a"], "median", "row"),
        ("4e", ["h"], "mean", None),
        ("4f", ["h"], "median", "row"),
        ("4g", ["m"], "mean", "row"),
        ("4h", ["m"], "median", "row"),
    ]

    for label, mods, stat, scale in configs:
        summary = summarise_gene_methylation(df, mods, stat)
        mat = to_heatmap_matrix(summary, f"{stat}_methylation")
        plot_heatmap(
            mat,
            outdir / f"figure{label}.tiff",
            scale=scale,
            cmap="vlag",
        )

#figure 5
##average methylation per site within gene
def run_figure5(df_bed, genes_pr, outdir):
    df_bed = df_bed[df_bed["mod_score"] >= 50]
    df = assign_genes(df_bed, genes_pr)
    df = df.dropna(subset=["gene"])

    configs = [
        ("5a", None),
        ("5b", ["a"]),
        ("5c", ["m"]),
        ("5d", ["h"]),
    ]

    for label, mods in configs:
        summary = summarise_gene_methylation(df, mods, stat="mean")
        mat = to_heatmap_matrix(summary, "mean_methylation")

        plot_heatmap(
            mat,
            outdir / f"figure{label}.pdf",
            cmap="vlag",
            vmin=0,
            vmax=1,
            figsize=(15, 25),
        )
