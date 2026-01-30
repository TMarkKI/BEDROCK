from figures.figure3ab import window_base_counts, plot_window_counts
from figures.figure3cd import summarize_modifications, plot_mod_windows


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
