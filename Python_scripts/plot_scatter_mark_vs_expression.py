#!/usr/bin/env python3
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

# 1) load the matrix
m = pd.read_csv(snakemake.input.matrix, sep="\t", index_col=0)

# 2) decide which columns to plot
marks = ["H3K27ac", "H3K9me3"]

# 3) compute the 99th percentile cap for RNA
rna_cap = np.nanpercentile(m["RNA"], 99)

for mark, out_png in zip(marks, snakemake.output):
    x = m[mark]
    # 4) cap RNA at 99th pct, drop log
    y = m["RNA"].clip(upper=rna_cap)

    plt.figure(figsize=(5,5))
    sns.scatterplot(x=x, y=y, alpha=0.4)



    # 6) Spearman correlation
    rho, _ = spearmanr(x, y, nan_policy="omit")
    plt.text(
        0.05, 0.95,
        f"Spearman ? = {rho:.2f}",
        transform=plt.gca().transAxes,
        va="top", ha="left",
        fontsize=10, color="darkred"
    )

    plt.xlabel(mark)
    plt.ylabel(f"RNA (capped at 99th pct = {rna_cap:.1f})")
    plt.title(f"{snakemake.wildcards.tissue}: {mark} vs RNA")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()
