##bases composition (failed/mod/canonical/variant)
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

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

    samples = df["sample_name"].unique()
    thresholds = df["threshold"].unique()
    
    base_palette = make_palette(order,palette=BASE_CALL_PALETTE)

    fig, axes = plt.subplots(
        nrows=len(samples),
        ncols=1,
        figsize=(6, 4 * len(samples)),
        sharex=True
    )


    if len(samples) == 1:
        axes = [axes]

    for ax, sample in zip(axes, samples):
        sdf = df[df["sample_name"] == sample]
        
        pivot = (
            sdf
            .pivot_table(
                index="threshold",
                columns="base_call_type",
                values="percent",
                aggfunc="sum",
                fill_value=0
            )
            .reindex(thresholds)
        )
        
        bottom = np.zeros(len(pivot))

        for base in order:
            if base not in pivot.columns:
                continue

            
            ax.bar(
                pivot.index,
                pivot[bc],
                bottom=bottom,
                color=base_palette.get(bc),
                label=bc
            )
            bottom += pivot[bc].values
            
            
        ax.set_title(sample)
        ax.set_ylabel("Percentage of Base-Type Called")
        ax.set_ylim(0, 100)

    axes[-1].set_xlabel("Depth Threshold")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        title="Type of Base",
        bbox_to_anchor=(1.02, 0.5),
        loc="center left"
    )

    plt.tight_layout()
    fig.savefig(f"{outdir}/figure2.tiff", dpi=300)
    plt.close()
