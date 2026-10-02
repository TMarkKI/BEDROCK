# genome-wide methylation summary (5mC / 5hmC / 6mA vs reference)
import pandas as pd
from Bio import SeqIO
 
MOD_BASE_MAP = {"m": "C", "h": "C", "a": "A"}
MOD_NAME_MAP = {"m": "5mC", "h": "5hmC", "a": "6mA"}

MOD_COL = "mod_score"
COV_COL = "mod_cov"
POS_COLS = ["Chromosome", "start_chrom_pos", "strand"]
 
 
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
 
def methylation_summary(samples, fasta_path, chrom_map, outdir, min_mod_reads=1):
    reference_counts, _ = count_reference_bases(fasta_path, chrom_map)
 
    rows = []
    for sample_name, sample in samples.items():
        bed = sample["bed"]
 
        for code in bed["mod_code"].unique():
            ref_base = MOD_BASE_MAP[code]
            pos = _to_positions(bed[bed["mod_code"] == code])
            row.append({
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
                 **_metrics(pos, reference_counts["C"]. min_mod_reads),
             })
          
    out = pd.DataFrame(rows)

    mod_order = ["5mC", "5hmC", "5mC+5hmC combined", "6mA"]
    out["modification"] = pd.Categorical(
        out["modification"],
        categories=mod_order + [m for m in out["modification"].unique() if m not in mod_order],
        ordered=True,
    )

    txt_path = f"{outdir}/methylation_summary.txt"
    with open(txt_path, "w") as f:
         f.write("Genome-wide methylation summary (all calls, no mod_score filter)\n")
        f.write(f"Position called modified if >= {min_mod_reads} modified read(s)\n")
        f.write("=" * 65 + "\n")
        for sample_name in out["sample_name"].unique():
            f.write(f"\nSample: {sample_name}\n")
            sub = out[out["sample_name"] == sample_name].sort_values("modification")
            for _, r in sub.iterrows():
                f.write(f"  {r['modification']}\n")
                f.write(
                    f"    reads:     {r['n_modified_reads']:>14,.0f} modified "
                    f"/ {r['total_reads']:>14,.0f} total = {r['percent_modified_reads']:.4f}%\n"
                )
                f.write(
                    f"    positions: {r['n_positions_modified']:>14,.0f} modified "
                    f"/ {r['total_ref_positions']:>14,.0f} reference {r['ref_base']} "
                    f"= {r['percent_ref_positions_modified']:.4f}% "
                    f"({r['percent_covered_positions_modified']:.4f}% of covered)\n"
                )

    print(out)
    return out
