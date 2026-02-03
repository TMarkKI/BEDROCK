import argparse

def get_args():
    parser = argparse.ArgumentParser(
        prog="BEDROCK",
        description="BEDROCK – Bedmethyl visualisation tool (Modkit output only)",
        epilog=(
            "Example usage:\n"
            "  python main.py "
            "--spreadsheet samples.csv "
            "--chr-list chr_map.csv "
            "--fai ref.fasta.fai "
            "--ref ref.fasta "
            "--annotation genes.gff3 "
            "--coverage-thresholds 50 100 150 200 250 300 "
            "--outdir results\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    
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

    return parser.parse_args()
