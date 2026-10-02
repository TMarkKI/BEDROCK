# genome-wide methylation summary (5mC / 5hmC / 6mA vs reference)
import pandas as pd
from Bio import SeqIO
 
MOD_BASE_MAP = {"m": "C", "h": "C", "a": "A"}
BASE_MOD_CODES = {"C": ["m", "h"], "A": ["a"]}
MOD_NAME_MAP = {"m": "5mC", "h": "5hmC", "a": "6mA"}

MOD_COL = "mod"
COV_COL = "mod_score"
POS_COLS = ["Chromosome", "Start_chrom_pos", "strand"]
 
 
def count_reference_bases(fasta_path, chrom_map=None):
    chrom_map = chrom_map or {}
    records = []
    genome_totals = {"A": 0, "C": 0}
    
    for record in SeqIO.parse(fasta_path, "fasta"):
        chrom = chrom_map.get(record.id, record.id)
        seq = str(record.seq).upper()
        a, t, c, g = (seq.count(b) for b in "ATCG")
        genome_totals["A"] += a + t
        genome_totals["C"] += c + g
        records.append({"Chromosome": chrom, "A": a, "T" : t, "C": c, "G" : g})
    
    return genome_totals, pd.DataFrame(records)

def _to_positions(df):
    return (
        df.groupby(POS_COLS, as_index=False)
        .agg(mod=(MOD_COL, "sum"), cov=(COV_COL, "max"))
    )

def _metrics(pos, ref_total, min_mod_reads):
    n_mod_reads = pos["mod"].sum()
    total_reads = pos["cov"].sum()
    n_pos_cov = len(pos)
    n_pos_mod = int((pos["mod"] >= min_mod_reads).sum())
    return {
        "n_modified_reads": n_mod_reads,
        "total_reads": total_reads,
        "percent_modified_reads": n_mod_reads / total_reads * 100 if total_reads else float("nan"),
        "n_positions_covered": n_pos_cov,
        "n_positions_modified": n_pos_mod,
        "total_ref_positions": ref_total,
        "percent_ref_positions_modified": n_pos_mod / ref_total * 100,
        "percent_covered_positions_modified": n_pos_mod / n_pos_cov *100 if n_pos_cov else float("nan"),
    }

def base_totals_summary(samples, reference_counts):
    rows = []
    for sample_name, sample in samples.items():
        bed = sample["bed"]
        for ref_base, codes in BASE_MOD_CODES.items():
            sub = bed[bed["mod_code"].isin(codes)]
            if sub.empty:
                continue
            pos = _to_positions(sub)
            n_mod = pos["mod"].sum()
            total_calls = pos["cov"].sum()
            rows.append({
                "sample_name": sample_name,
                "ref_base": ref_base,
                "total_ref_bases": reference_counts[ref_base],
                "n_positions_covered": len(pos),
                "total_calls": total_calls,
                "n_modified_calls": n_mod,
                "n_unmodified_calls": total_calls - n_mod,
                "percent_modified_calls": n_mod / total_calls * 100 if total_calls else float("nan"),
                "mean_depth_per_covered_position": total_calls / len(pos),
            })
    return pd.DataFrame(rows)
 
def methylation_summary(samples, fasta_path, chrom_map, outdir, min_mod_reads=1):
    reference_counts, _ = count_reference_bases(fasta_path, chrom_map)

    rows = []
    for sample_name, sample in samples.items():
        bed = sample["bed"]

        for code in bed["mod_code"].unique():
            ref_base = MOD_BASE_MAP[code]
            pos = _to_positions(bed[bed["mod_code"] == code])
            rows.append({
                "sample_name": sample_name,
                "mod_code": code,
                "modification": MOD_NAME_MAP[code],
                "ref_base": ref_base,
                **_metrics(pos, reference_counts[ref_base], min_mod_reads),
            })

        c_bed = bed[bed["mod_code"].isin(["m", "h"])]
        if not c_bed.empty:
            pos = _to_positions(c_bed)
            rows.append({
                "sample_name": sample_name,
                "mod_code": "m+h",
                "modification": "5mC+5hmC combined",
                "ref_base": "C",
                **_metrics(pos, reference_counts["C"], min_mod_reads),
            })

    out = pd.DataFrame(rows)

    mod_order = ["5mC", "5hmC", "5mC+5hmC combined", "6mA"]
    out["modification"] = pd.Categorical(
        out["modification"],
        base_totals = base_totals_summary(samples, reference_counts),
        base_totals.to_csv(f"{outdir}/base_totals_summary.tsv", sep="\t", index=False),
        categories=mod_order + [m for m in out["modification"].unique() if m not in mod_order],
        ordered=True,
    )

    txt_path = f"{outdir}/methylation_summary.txt"
    with open(txt_path, "w") as f:
        f.write("Genome-wide methylation summary (all calls, no mod_score filter)\n")
        f.write(f"Position called modified if >= {min_mod_reads} modified read(s)\n")
        f.write("=" * 65 + "\n")
        f.write("Total base calls (C = 5mC+5hmC, A = 6mA)\n")
        for _, r in base_totals.iterrows():
            f.write(
                f"\nSample: {r['sample_name']}  base: {r['ref_base']}\n"
                f"  calls:     {r['n_modified_calls']:>14,.0f} modified "
                f"/ {r['total_calls']:>14,.0f} total = {r['percent_modified_calls']:.4f}%\n"
                f"  coverage:  {r['n_positions_covered']:>14,.0f} of {r['total_ref_bases']:,.0f} "
                f"reference {r['ref_base']} positions covered, "
                f"mean depth {r['mean_depth_per_covered_position']:.1f}\n"
            )

    print(out)
    return out
