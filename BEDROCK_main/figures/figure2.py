##bases composition (failed/mod/canonical/variant)
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from utils import make_palette
from config import BASE_CALL_PALETTE


def base_composition(df, threshold):
    df = df[df.mod_score >= threshold]

    summary = (
        df.groupby(["Chromosome", "sample_name"])
        [["mod", "canonical_mod", "other_mod",
          "delete_mod", "fail_mod", "variant_mod", "no_call"]]
        .sum()
        .reset_index()
    )

    long = summary.melt(
        id_vars=["Chromosome", "sample_name"],
        var_name="base_call_type",
        value_name="n_sites"
    )

    long["percent"] = (
        long.groupby(["Chromosome", "sample_name"])
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



    g = sns.histplot(
        data=df,
        x="Chromosome",
        y="percent",
        hue="base_call_type",
        row="threshold",
        col="sample_name",
        multiple="stacked",
        kind="bar",
        height=4,
        palette=base_palette
    )
    g.set_xticklabels(rotation=45)
    g.set_axis_labels("Chromosome", "Percentage of Bases Called")
    
    plt.tight_layout()
    g.savefig(f"{outdir}/figure2a.tiff", dpi=300)
    plt.close()
