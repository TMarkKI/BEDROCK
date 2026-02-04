from pathlib import Path
import pandas as pd
import sys
import json

from args import get_args
from data_loading import load_chr_map, load_fai, load_samples, load_genes_as_pyranges
from check_install import check_install
from figures.figure1 import figure_1a, figure_1b
from figures.figure2 import base_composition, plot_base_composition
from figures.merge_figure3 import run_figure3
from figures.figure4n5 import run_figure4, assign_genes, run_figure5
from figures.figure6 import run_figure6

def main():
    print(f"BEDROCK - Bedmethyl Visualisation Tool; V0.1; used Modkit output only")

    # ---- check packages ----
    check_install()
    
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
        [s["depth"].assign(sample_name=name) for name, s in samples.items()]
    )

    bed_all = pd.concat(
        [s["bed"].assign(sample_name=name) for name, s in samples.items()]
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
