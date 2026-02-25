##genes and methylation
import pyranges as pr
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

chrom_map = {
    "Pf3D7_01_v3": "1",
    "Pf3D7_02_v3": "2",
    "Pf3D7_03_v3": "3",
    "Pf3D7_04_v3": "4",
    "Pf3D7_05_v3": "5",
    "Pf3D7_06_v3": "6",
    "Pf3D7_07_v3": "7",
    "Pf3D7_08_v3": "8",
    "Pf3D7_09_v3": "9",
    "Pf3D7_10_v3": "10",
    "Pf3D7_11_v3": "11",
    "Pf3D7_12_v3": "12",
    "Pf3D7_13_v3": "13",
    "Pf3D7_14_v3": "14",
    "Pf3D7_API_v3": "API",
    "Pf_M76611": "MIT",
}

def assign_genes(df_bed, genes_pr):

    mods = pr.PyRanges(
        chromosomes=df_bed["Chromosome"],
        starts=df_bed["Start_chrom_pos"],
        ends=df_bed["End_chrom_pos"],
    )

    overlaps = mods.join(genes_pr)
    odf = overlaps.df

    if "Name" in odf.columns:
        gene_col = "Name"
    elif "gene_id" in odf.columns:
        gene_col = "gene_id"
    elif "ID" in odf.columns:
        gene_col = "ID"
    else:
        raise RuntimeError("No gene identifier column found in GFF (expected Name, gene_id, or ID)")

    odf = odf[
        ["Chromosome", "Start", "End", gene_col]
    ].rename(columns={gene_col: "gene"})

    df_bed = df_bed.merge(
        odf,
        left_on=["Chromosome", "Start_chrom_pos", "End_chrom_pos"],
        right_on=["Chromosome", "Start", "End"],
        how="left",
    )

    df_bed.drop(columns=["Start", "End"], inplace=True)
    
    return df_bed


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

def select_top_500_genes(matrix, n=500):
    variances = matrix.var(axis=1)
    top_genes = variances.sort_values(ascending=False).head(n).index
    return matrix.loc[top_genes]

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
        row_cluster=False,
        column_cluster=True,
        yticklabels=False,
    )
    plt.savefig(outfile, dpi=300)
    plt.close()

def remap_chromosomes(genes_pr, chrom_map):
    df = genes_pr.df.copy()

    df["Chromosome"] = df["Chromosome"].astype(str).replace(chrom_map)

    return pr.PyRanges(df)

def run_figure4(df_bed, genes_pr, outdir):
    genes_pr = remap_chromosomes(genes_pr, chrom_map)
    
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

        top_500 = select_top_500_genes(mat, n=500)

        plot_heatmap(
            top_500,
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
