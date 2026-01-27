import pandas as pd
import pyranges as pr
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

def load_genes_as_pyranges(gff_path):
    df = pd.read_csv(
        gff_path,
        sep="\t",
        comment="#",
        header=None,
        names=[
            "Chromosome", "source", "feature",
            "start", "end", "score",
            "strand", "frame", "attributes"
        ]
    )

    df = df[df["feature"] == "gene"].copy()

    if df.empty:
        raise ValueError("No 'gene' features found in GFF file")

    df["gene_id"] = df["attributes"].str.extract(r"ID=([^;]+)")

    if df["gene_id"].isna().any():
        n_missing = df["gene_id"].isna().sum()
        print(
            f"[WARN] {n_missing} gene entries missing gene_id; "
            "loaded but without name"
        )

    df["Start"] = df["start"] - 1
    df["End"] = df["end"]

    df = df[[
        "Chromosome",
        "Start",
        "End",
        "strand",
        "gene_id"
    ]]

    return pr.PyRanges(df)    
