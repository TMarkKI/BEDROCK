# %%
##python project file for methylation analysis from nanopore data
##python main.py --spreadsheet samples.csv --chr-list chr_map.csv --fai ref.fasta.fai --ref ref.fasta --annotation genes.gff3 --coverage-thresholds 50 100 150 200 250 300
##optional:  --outdir results --no-log-depth

##Project layout
##methylation_analysis/
##  args.py
##  config.py
##  io.py
##  utils.py
##  data_loading.py
##  figures/
##      fig1.py
##      fig2.py
##      fig3.py
##      fig4.py
##      fig5.py
##      fig6.py
##  main.py
##requirements.txt

# %%
# args.py
import argparse

def get_args():
    parser = argparse.ArgumentParser(
        description="Methylation analysis from bedmethyl generated from Nanopore data"
    )

    parser.add_argument(
        "--spreadsheet",
        required=True,
        help="CSV with #1;sample_name #2;sample_type #3;bedmethyl_path #4;depth_path"
    )

    parser.add_argument(
        "--chr-list",
        required=True,
        help="Chr list should have #1=ref_chr;ref_chr name #2=assigned_chr;assigned chr number/name"
    )

    parser.add_argument(
        "--fai",
        required=True,
        help="Reference FASTA index (.fai)"
    )

    parser.add_argument(
        "--ref",
        required=True,
        help="Reference FASTA file (.fasta)"
    )

    parser.add_argument(
        "--annotation",
        required=True,
        help="GFF3 annotation file (.gff)"
    )

    parser.add_argument(
        "--outdir",
        default="results",
        help="Output directory (default: results)"
    )

    parser.add_argument(
        "--coverage-thresholds",
        type=int,
        nargs="+",
        default=[50, 100, 150, 200, 250, 300],
        help="Coverage thresholds for QC and composition plots"
    )

    parser.add_argument(
        "--no-log-depth",
        action="store_true",
        help="Disable log10 scaling for depth plots"
    )

    parser.add_argument(
        "--help",
        help="python main.py --spreadsheet samples.csv --chr-list chr_map.csv --fai ref.fasta.fai --ref ref.fasta --annotation genes.gff3 --coverage-thresholds (default: 50 100 150 200 250 300) OPTIONAL:  --outdir results --no-log-depth"
    )

    return parser.parse_args()

# %%
#config.py
##read all datafiles and set thresholds

BEDMETHYL_COLUMNS = [
    "Chromosome",
    "Start_chrom_pos",
    "End_chrom_pos",
    "mod_code",
    "mod_score",
    "strand",
    "strand_pos_start",
    "strand_pos_end",
    "color",
    "mod_cov",
    "percent_mod",
    "mod",
    "canonical_mod",
    "other_mod",
    "delete_mod",
    "fail_mod",
    "variant_mod",
    "no_call",
]

DEPTH_COLUMNS = [
    "Chromosome",
    "Position",
    "Coverage_depth",
]

DROP_BEDMETHYL_COLUMNS = {
    "color",
    "mod_cov",
}


SAMPLE_PALETTE = "colorblind"        # samples
CHROM_PALETTE = "tab20"         # chromosomes
BASE_CALL_PALETTE = "Set2"      # mod / canonical / etc
MODIFICATION_PALETTE = "Dark2"  # 5mC / 5hmC / 6mA
EXTENDED_PALETTE = "hls"

# %%
#io.py
#input/output

# io.py
from pathlib import Path
import pandas as pd
import pyranges as pr

def load_chr_map(chr_list_path):
    df = pd.read_csv(chr_list_path, sep="\t", header=0)
    return dict(zip(df[0], df[1]))

def load_fai(fai_path, chr_map):
    df = pd.read_csv(
        fai_path,
        sep="\t",
        header=None,
        names=["chrom", "length", "offset", "linebases", "linewidth"]
    )
    df["chrom"] = df["chrom"].map(chr_map)

    if df["chrom"].isna().any():
        missing = df.loc[df["chrom"].isna(), "chrom"].unique()
        raise ValueError(f"Unmapped chromosomes in FAI: {missing}")

    return df



def load_samples(spreadsheet, chr_map):
    sample_list = pd.read_csv(spreadsheet)

    required_cols = {
            "sample_name", "sample_type",
            "bedmethyl_path", "depth_path"
    }

    missing = required_cols - set(sample_list.columns)
    if missing:
        raise ValueError(f"Spreadsheet missing columns: {missing}")

    samples = {}

    for _, row in sample_list.iterrows():
        name = row["sample_name"]

        bed = pd.read_csv(row["bedmethyl_path"], sep="\t")
        depth = pd.read_csv(row["depth_path"], sep="\t")

        bed = bed.rename(columns={"chrom": "Chromosome"})
        depth = depth.rename(columns={"chrom": "Chromosome"})

        bed["Chromosome"] = bed["Chromosome"].map(chr_map)
        if bed["Chromosome"].isna().any():
            raise ValueError(f"Unmapped chromosomes in BED for sample {name}")

        depth["Chromosome"] = depth["Chromosome"].map(chr_map)
        if depth["Chromosome"].isna().any():
            raise ValueError(f"Unmapped chromosomes in BED for sample {name}")

        samples[name] = {
            "bed": bed,
            "depth": depth,
            "type": row["sample_type"],
        }
        
    return samples

def load_genes_as_pyranges(gff_path):
    df = pd.read_csv(
        gff_path,
        sep="\t",
        comment="#",
        header=None,
        names=[
            "seqnames", "source", "feature",
            "start", "end", "score",
            "strand", "frame", "attributes"
        ]
    )
    
    df = df.rename(columns={"seqnames": "Chromosome"})

    df = df[df["feature"] == "gene"]

    df["gene_id"] = (
        df["attributes"]
        .str.extract(r"ID=([^;]+)")
    )

    df["Start"] = df["start"] - 1
    df["End"] = df["end"]

    return pr.PyRanges(df.drop(columns=["start", "end"]))

# %%
#utils.py

import seaborn as sns

def make_palette(categories, palette="colorblind", fallback = "hls"):
    categories = sorted(set(categories))
    n = len(categories)
    
    try:
        colours = sns.color_palette(palette, n)
    except ValueError:
        colours = sns.color_palette(fallback, n)

    return dict(zip(categories, colours))

# %%
#data_loading.py
##loading and preprocessing of spreadsheets and datafiles

import pandas as pd
from Bio import SeqIO
from config import BEDMETHYL_COLUMNS, DROP_BEDMETHYL_COLUMNS

def load_chr_map(chr_list_file):
    df = pd.read_csv(chr_list_file)
    return dict(zip(df["ref_chr"], df["assigned_chr"]))

def load_fai(fai_file, chr_map):
    df = pd.read_csv(
        fai_file, sep="\t", header=None,
        names=["Chromosome", "length", "offset", "linebases", "linewidth"]
    )
    df["Chromosome"] = df["Chromosome"].replace(chr_map)
    return df

def load_samples(spreadsheet, chr_map):
    meta = pd.read_csv(spreadsheet)
    samples = {}

    for _, row in meta.iterrows():
        bed = pd.read_csv(row.bedmethyl_path, sep="\t", header=None)

        n_expected = len(BEDMETHYL_COLUMNS)
        n_found = bed.shape[1]

        if n_found > n_expected:
            print(
                f"[WARN] {row.sample}: BEDMethyl file has {n_found} columns; "
                f"only the first {n_expected} will be used. "
                f"Extra columns ignored: {n_found - n_expected}"
            )
        elif n_found < n_expected:
            raise ValueError(
                f"[ERROR] {row.sample}: BEDMethyl file has {n_found} columns, "
                f"but {n_expected} expected. File format may be incompatible. "
                f"Recheck if in Bedmethyl format (remember only ModKit bedmethyl format is compatible with this BEDROCK version)."
            )

        bed = bed.iloc[:, :len(BEDMETHYL_COLUMNS)]
        bed.columns = BEDMETHYL_COLUMNS

        print(f"[INFO] {row.sample}: BEDMethyl column mapping:")

        for i, col in enumerate(BEDMETHYL_COLUMNS):
            print(f"  column {i:02d} → {col}")


        bed = bed.drop(columns=[c for c in DROP_BEDMETHYL_COLUMNS if c in bed.columns])
        bed["Chromosome"] = bed["Chromosome"].replace(chr_map)
        unmapped = bed["Chromosome"].isna().sum()

        if unmapped > 0:
            print(
                f"[WARN] {row.sample}: {unmapped} BED entries could not be "
                "mapped via chr_map and were set to NA"
            )

        depth = pd.read_csv(row.depth_path, sep="\t", header=None)
        depth.columns = ["Chromosome", "Position", "Coverage_depth"]
        depth["Chromosome"] = depth["Chromosome"].replace(chr_map)

        samples[row.sample] = {
            "type": row.type,
            "bed": bed,
            "depth": depth,
            "mean_depth": depth.Coverage_depth.mean()
        }

        print(f"Loaded {row.sample} | Mean depth {depth.Coverage_depth.mean():.2f} | Depth file {row.depth_path}")

    return samples


# %%
#figure 1
##create figure 1 for QC = depth/coverage

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from utils import make_palette
from config import SAMPLE_PALETTE

def figure_1a(samples, outdir):
    df = pd.concat(
        [
            s["depth"].assign(sample=name)
            for name, s in samples.items()
        ]
    )

    sample_palette = make_palette(
        df["sample"].unique(),
        palette=SAMPLE_PALETTE
    )

    summary = (
        df.groupby(["Chromosome", "sample"])
          .Coverage_depth.mean()
          .reset_index()
    )

    plt.figure(figsize=(8, 6))
    sns.barplot(
        data=summary,
        x="Chromosome",
        y="Coverage_depth",
        hue="sample",
        palette=sample_palette
    )
    plt.yscale("log")
    plt.ylabel("Mean Coverage Depth (log10)")
    plt.tight_layout()
    plt.savefig(f"{outdir}/figure_1a.tiff", dpi=300)
    plt.close()


def figure_1b(df_depth_all, fai, thresholds):
    out = []

    for t in thresholds:
        tmp = (
            df_depth_all[df_depth_all.Coverage_depth >= t]
            .groupby(["Chromosome", "sample"])
            .size()
            .reset_index(name="covered")
            .merge(fai, on="Chromosome")
        )
        tmp["percent"] = tmp.covered / tmp.length * 100
        tmp["threshold"] = f">={t}x"
        out.append(tmp)

    df = pd.concat(out)

    sample_palette = make_palette(
        df["sample"].unique(),
        palette=SAMPLE_PALETTE
    )


    g = sns.catplot(
        data=df,
        x="Chromosome",
        y="percent",
        hue="sample",
        col="threshold",
        kind="bar",
        col_wrap=3,
        height=4,
        palette=sample_palette
    )
    g.set_axis_labels("Chromosome", "Percent Covered")
    g.savefig(f"{outdir}/figure_1b.tiff", dpi=300)


# %%
#figure 2
##bases composition (failed/mod/canonical/variant)

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from utils import make_palette
from config import BASE_CALL_PALETTE


def base_composition(df, threshold):
    df = df[df.mod_score >= threshold]

    summary = (
        df.groupby(["Chromosome", "sample"])
        [["mod", "canonical_mod", "other_mod",
          "delete_mod", "fail_mod", "variant_mod", "no_call"]]
        .sum()
        .reset_index()
    )

    long = summary.melt(
        id_vars=["Chromosome", "sample"],
        var_name="base_call_type",
        value_name="n_sites"
    )

    long["percent"] = (
        long.groupby(["Chromosome", "sample"])
            .n_sites.transform(lambda x: x / x.sum() * 100)
    )

    long["threshold"] = f">={threshold}x"
    return long

def plot_base_composition(df, outdir):
    order = [
        "canonical_mod",
        "mod",
        "variant_mod",
        "other_mod",
        "fail_mod",
        "delete_mod",
        "no_call",
    ]

    df["base_call_type"] = pd.Categorical(
        df["base_call_type"],
        categories=order,
        ordered=True
    )        
    
    base_palette = make_palette(
        df["base_call_type"].unique(),
        palette=BASE_CALL_PALETTE
    )



    g = sns.catplot(
        data=df,
        x="Chromosome",
        y="percent",
        hue="base_call_type",
        row="threshold",
        col="sample",
        kind="bar",
        multiple="stack",
        height=4,
        palette=base_palette
    )
    g.set_xticklabels(rotation=45)
    g.set_axis_labels("Chromosome", "Percentage of Bases Called")
    
    plt.tight_layout()
    g.savefig(f"{outdir}/figure2a.tiff", dpi=300)
    plt.close()

# %%
#figure 3a and 3b
##AT or CG pairs per 1kb over forward and reverse strands
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO

from utils import make_palette


WINDOW_SIZE = 1000

##count bases in FR strands - fixed window size
def window_base_counts(fasta_path, chrom_map, bases):

    records = []

    for record in SeqIO.parse(fasta_path, "fasta"):
        chrom = chrom_map.get(record.id, record.id)
        seq = str(record.seq).upper()
        rev = str(record.seq.reverse_complement()).upper()

        for start in range(0, len(seq), WINDOW_SIZE):
            end = min(start + WINDOW_SIZE, len(seq))

            for base, strand_seq, strand in [
                (bases[0], seq, "+"),
                (bases[1], rev, "-"),
            ]:
                count = strand_seq[start:end].count(base)
                records.append(
                    dict(
                        Chromosome=chrom,
                        start=start,
                        end=end,
                        strand=strand,
                        Count=count if strand == "+" else -count,
                    )
                )

    return pd.DataFrame(records)


def plot_window_counts(df, outpath, ylab):

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
        col="Chromosome",
        col_wrap=4,
        sharex=False,
        sharey=True,
        height=3,
    )

    g.map_dataframe(
        sns.barplot,
        x="start",
        y="Count",
        hue="strand",
        palette=palette,
        dodge=False,
    )

    g.add_legend()
    g.set_axis_labels("Genomic position (bp)", ylab)

    for ax in g.axes.flatten():
        ax.tick_params(axis="x", labelbottom=False)

    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()


# %%
#figure 3c and 3d
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
            ["Chromosome", "start", "end", "strand", "sample"]
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
        row="sample",
        col="Chromosome",
        sharex=False,
        height=3,
    )

    g.map_dataframe(
        sns.barplot,
        x="start",
        y="mod",
        hue="strand",
        palette=palette,
        dodge=False,
    )

    g.add_legend()
    g.set_axis_labels("Genomic position (bp)", ylab)

    for ax in g.axes.flatten():
        ax.tick_params(axis="x", labelbottom=False)

    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()


# %%
#merge figure 3

from figures.fig3_reference_windows import window_base_counts, plot_window_counts
from figures.fig3_modification_windows import summarize_modifications, plot_mod_windows


def run_figure3(ref_fasta, bed_all, chr_map, outdir):

    # ---- FIGURE 3A/B ----
    cg_df = window_base_counts(ref_fasta, chr_map, bases=("C", "G"))
    plot_window_counts(cg_df, outdir / "figure3a.tiff", "CG Count")

    at_df = window_base_counts(ref_fasta, chr_map, bases=("A", "T"))
    plot_window_counts(at_df, outdir / "figure3b.tiff", "AT Count")

    # ---- FIGURE 3C/D ----
    cg_mods = summarize_modifications(bed_all, cg_df, mod_codes=["h", "m"])
    plot_mod_windows(cg_mods, outdir / "figure3c.tiff", "Modified Cytosine Count")

    at_mods = summarize_modifications(bed_all, at_df, mod_codes=["a"])
    plot_mod_windows(at_mods, outdir / "figure3d.tiff", "Modified Adenosine Count")


# %%
#figure 4&5
##genes and methylation
import pyranges as pr
import seaborn as sns
import matplotlib.pyplot as plt


def assign_genes(df_bed, genes_pr):
    mods = pr.PyRanges(
        chromosomes=df_bed["Chromosome"],
        starts=df_bed["Start_chrom_pos"],
        ends=df_bed["End_chrom_pos"],
        sample=df_bed["sample"],
        mod_code=df_bed["mod_code"],
        mod_score=df_bed["mod_score"],
        percent_mod=df_bed["percent_mod"],
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
        df.groupby(["gene", "sample"])["percent_mod"]
        .agg(agg)
        .reset_index(name=f"{stat}_methylation")
    )

    return summary

def to_heatmap_matrix(df, value_col):
    return (
        df.pivot(index="gene", columns="sample", values=value_col)
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

# %%
#figure 6
##differential methylation
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np

def compute_gene_methylation(mod_summary, threshold):
    df = mod_summary.copy()

    df["methylated"] = df["mod_percentage"] >= threshold

    summary = (
        df.groupby(["gene", "sample"])
        .agg(
            total_sites=("mod_percentage", "size"),
            methylated_sites=("methylated", "sum"),
        )
        .reset_index()
    )

    summary["proportion_methylated"] = (
        summary["methylated_sites"] / summary["total_sites"]
    )

    return summary

def attach_gene_positions(methylation_df, gene_positions):
    gene_positions_unique = gene_positions.drop_duplicates("gene_id")

    out = methylation_df.merge(
        gene_positions_unique,
        left_on="gene",
        right_on="gene_id",
        how="left",
    )

    out["chr"] = out["seqnames"].astype(str)
    out["pos"] = out["midpoint"]

    return out

def plot_manhattan_by_position(df, outpath, title):
    g = sns.FacetGrid(
        df,
        row="sample",
        col="chr",
        sharex=False,
        sharey=True,
        height=2.5,
        aspect=1.2,
    )

    g.map_dataframe(
        sns.scatterplot,
        x="pos",
        y="proportion_methylated",
        alpha=0.7,
        s=15,
    )

    g.set_titles(row_template="{row_name}", col_template="Chr {col_name}")
    g.set_axis_labels("Genomic position", "Proportion methylated")

    for ax in g.axes.flat:
        ax.set_xticks([])
        ax.grid(False, axis="x")

    g.fig.suptitle(title, y=1.02)
    g.savefig(outpath, dpi=300)
    plt.close()

def plot_gene_manhattan(df, outpath, title):
    plt.figure(figsize=(30, 10))
    sns.scatterplot(
        data=df,
        x="gene",
        y="proportion_methylated",
        hue="sample",
        alpha=0.7,
        s=20,
    )

    plt.xticks(rotation=90, fontsize=5)
    plt.xlabel("Gene")
    plt.ylabel("Proportion methylated")
    plt.title(title)
    plt.grid(False, axis="x")
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()

def top_genes_overall(df, prop=0.05):
    n = int(len(df["gene"].unique()) * prop)

    top = (
        df.groupby("gene")["proportion_methylated"]
        .max()
        .nlargest(n)
        .index
    )

    return df[df["gene"].isin(top)]

def top_genes_per_sample(df, prop=0.05):
    return (
        df.sort_values("proportion_methylated", ascending=False)
        .groupby("sample")
        .head(lambda x: int(np.ceil(len(x) * prop)))
    )


def plot_top_genes(df, outpath, title):
    plt.figure(figsize=(15, 10))
    sns.scatterplot(
        data=df,
        x="gene",
        y="proportion_methylated",
        hue="sample",
        s=25,
        alpha=0.8,
    )
    plt.xticks(rotation=90, fontsize=6)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()


def plot_top_genes_by_sample(df, outpath, title):
    g = sns.FacetGrid(
        df,
        col="sample",
        col_wrap=2,
        sharex=False,
        height=4,
    )

    g.map_dataframe(
        sns.scatterplot,
        x="gene",
        y="proportion_methylated",
        alpha=0.7,
        s=20,
    )

    for ax in g.axes.flat:
        ax.tick_params(axis="x", rotation=90, labelsize=6)
        ax.grid(False, axis="x")

    g.fig.suptitle(title, y=1.02)
    g.savefig(outpath)
    plt.close()

def run_figure6(mod_summary, gene_positions, outdir):
    outdir = Path(outdir)

    for threshold, label in [(10, "10_percent"), (5, "5_percent"), (1, "1_percent")]:
        methyl = compute_gene_methylation(mod_summary, threshold)
        methyl_pos = attach_gene_positions(methyl, gene_positions)

        # 6a / 6f
        plot_manhattan_by_position(
            methyl_pos,
            outdir / f"figure6a_{label}.tiff",
            f"Manhattan plot (≥{threshold}%)",
        )

        # 6b / 6g
        plot_gene_manhattan(
            methyl,
            outdir / f"figure6b_{label}.pdf",
            f"Gene-level methylation (≥{threshold}%)",
        )

        # top 5% overall
        top_overall = top_genes_overall(methyl)
        plot_top_genes(
            top_overall,
            outdir / f"figure6c_{label}.pdf",
            "Top 5% methylated genes",
        )

        # top 5% per sample
        top_sample = top_genes_per_sample(methyl)
        plot_top_genes_by_sample(
            top_sample,
            outdir / f"figure6d_{label}.pdf",
            "Top 5% genes per sample",
        )

        # export table
        top_overall.to_csv(
            outdir / f"top5_percent_genes_methylated_{label}.csv",
            index=False,
        )


# %%
# main.py
from pathlib import Path
import pandas as pd
import sys
import json

from args import get_args
from io import load_chr_map, load_fai, load_samples, load_genes_as_pyranges
from figures.fig1_qc import figure_1a, figure_1b
from figures.fig2_base_comp import base_composition, plot_base_composition
from figures.fig3 import run_figure3
from figures.fig4_gene_heatmaps import run_figure4, assign_genes
from figures.fig5_site_heatmaps import run_figure5
from figures.fig6_differential_methylation import run_figure6
#from config import SAMPLE_COLOURS, BASE_COLOURS

def main():
    print(f"BEDROCK - Bedmethyl Visualisation Tool; V0.1; used Modkit output only")
    
    args = get_args()

    # ---- validate argument inputs ----

    required_files = {
        "spreadsheet": args.spreadsheet,
        "chr_list": args.chr_list,
        "fai": args.fai,
        "ref": args.ref,
        "annotation": args.annotation,
    }

    for name, path in required_files.items():
        p = Path(path)
        if not p.exists():
            sys.exit(f"[ERROR] {name} file not found: {p}")
        if p.stat().st_size == 0:
            sys.exit(f"[ERROR] {name} file is empty: {p}")
    

    # ---- create output directory ----
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    with open(outdir / "run_args.json", "w") as fh:
        json.dump(vars(args), fh, indent=2)


    # ---- load data ----
    chr_map = load_chr_map(args.chr_list)
    fai = load_fai(args.fai, chr_map)
    thresholds = args.coverage_thresholds

    samples = load_samples(args.spreadsheet, chr_map)

    # flatten depth once
    df_depth_all = pd.concat(
        [s["depth"].assign(sample=name) for name, s in samples.items()]
    )

    bed_all = pd.concat(
        [s["bed"].assign(sample=name) for name, s in samples.items()]
    )

    # ---- FIGURE 1 ----
    figure_1a(
        samples,
        outdir
    )

    figure_1b(
        df_depth_all,
        fai,
        thresholds,
        outdir=outdir,
    )

    # ---- FIGURE 2 ----

    df_base = pd.concat(
        base_composition(bed_all, t) for t in thresholds
    )

    plot_base_composition(
        df_base,
        outdir
    )

    # ---- FIGURE 3A/B ----
    run_figure3(
        ref_fasta=args.ref,
        bed_all=bed_all,
        chr_map=chr_map,
        outdir=outdir,
    )

    # ---- FIGURE 4&5 ----
    genes_pr = load_genes_as_pyranges(args.annotation)

    run_figure4(
        bed_all, 
        genes_pr, 
        outdir,
    )

    run_figure5(
        bed_all, 
        genes_pr, 
        outdir,
    )

    # ---- FIGURE 6 ----

    df_gene = assign_genes(bed_all, genes_pr)
    df_gene = df_gene.dropna(subset=["gene"])

    mod_summary = (
        df_gene[["gene", "sample", "percent_mod"]]
        .rename(columns={"percent_mod": "mod_percentage"})
    )

    gene_positions = (
        genes_pr.df
        .assign(midpoint=lambda x: (x.Start + x.End) // 2)
        [["gene_id", "seqnames", "midpoint"]]
    )

    run_figure6(
        mod_summary,
        gene_positions,
        outdir,
    )


if __name__ == "__main__":
    main()


