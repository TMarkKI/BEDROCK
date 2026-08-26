# genome-wide methylation summary (5mC / 5hmC / 6mA vs reference)
import pandas as pd
from Bio import SeqIO
 
MOD_BASE_MAP = {"m": "C", "h": "C", "a": "A"}
MOD_NAME_MAP = {"m": "5mC", "h": "5hmC", "a": "6mA"}
 
 
def count_reference_bases(fasta_path, chrom_map=None):
    chrom_map = chrom_map or {}
    records = []
    genome_totals = {"A": 0, "C": 0}
    
    for record in SeqIO.parse(fasta_path, "fasta"):
        chrom = chrom_map.get(record.id, record.id)
        seq = str(record.seq).upper()
        a_count = seq.count("A")
        t_count = seq.count("T")
        c_count = seq.count("C")
        g_count = seq.count("G")
        genome_totals["A"] += a_count
        genome_totals["A"] += t_count
        genome_totals["C"] += c_count
        genome_totals["G"] += g_count
        records.append({"Chromosome": chrom, "A": a_count, "C": c_count})
    per_chrom = pd.DataFrame(records)
    
    return genome_totals, per_chrom
 
 
def methylation_summary(samples, fasta_path, chrom_map, outdir):
    reference_counts, _ = count_reference_bases(fasta_path, chrom_map)
 
    rows = []
    for sample_name, sample in samples.items():
        bed = sample["bed"]
 
        per_mod = (
            bed.groupby("mod_code")["mod"]
            .sum()
            .reset_index()
            .rename(columns={"mod": "n_modified"})
        )
        per_mod["ref_base"] = per_mod["mod_code"].map(MOD_BASE_MAP)
        per_mod["modification"] = per_mod["mod_code"].map(MOD_NAME_MAP)
        per_mod["total_ref_bases"] = per_mod["ref_base"].map(reference_counts)
        per_mod["percent_of_reference"] = (
            per_mod["n_modified"] / per_mod["total_ref_bases"] * 100
        )
        per_mod["sample_name"] = sample_name
        rows.append(per_mod)
 
        c_sub = per_mod[per_mod.ref_base == "C"]
        if not c_sub.empty:
            n_mod_c = c_sub["n_modified"].sum()
            total_c = reference_counts["C"]
            rows.append(pd.DataFrame([{
                "mod_code": "m+h",
                "n_modified": n_mod_c,
                "ref_base": "C",
                "modification": "5mC+5hmC combined",
                "total_ref_bases": total_c,
                "percent_of_reference": n_mod_c / total_c * 100,
                "sample_name": sample_name,
            }]))
 
    out = pd.concat(rows, ignore_index=True)
    
    mod_order = ["5mC", "5hmC", "5mC+5hmC combined", "6mA"]
    out["modification"] = pd.Categorical(
        out["modification"],
        categories=mod_order + [m for m in out["modification"].unique() if m not in mod_order],
        ordered=True,
    )
 
    txt_path = f"{outdir}/methylation_summary.txt"
    with open(txt_path, "w") as f:
        f.write("Genome-wide methylation summary (all calls, no mod_score filter)\n")
        f.write("=" * 65 + "\n")
        for sample_name in out["sample_name"].unique():
            f.write(f"\nSample: {sample_name}\n")
            sub = out[out["sample_name"] == sample_name].sort_values("modification")
            for _, row in sub.iterrows():
                f.write(
                    f"  {row['modification']:<20s} "
                    f"{row['n_modified']:>12,.0f} modified {row['ref_base']} "
                    f"/ {row['total_ref_bases']:>14,.0f} reference {row['ref_base']} "
                    f"= {row['percent_of_reference']:.4f}%\n"
                )
 
    print(out)
    
    return out
