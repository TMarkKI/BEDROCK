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
