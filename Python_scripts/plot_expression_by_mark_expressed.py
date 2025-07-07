#!/usr/bin/env python3
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

# 1) load
df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)

# 2) choose your flank marks
h3a = next(c for c in df if c.startswith("H3K27ac") and c.endswith("_flank"))
h3r = next(c for c in df if c.startswith("H3K9me3") and c.endswith("_flank"))

# 3) classify by chromatin context
df["group"] = np.where(df[h3a] > df[h3r], "active", "inactive")

# — NEW: 4) filter to expressed TEs only
df = df[df["RNA"] > 1].copy()

# 5) log2-transform for plotting
df["RNA_log2"] = np.log2(df["RNA"] + 1)

# 6) counts for each category
counts = df["group"].value_counts()

# 7) plot
fig, ax = plt.subplots(figsize=(5,4))
sns.violinplot(
    x="group", y="RNA_log2", data=df,hue="group", legend=False,
    palette={"active":"#fc8d62","inactive":"#66c2a5"},
    inner="quartile", cut=0, ax=ax
)
sns.stripplot(
    x="group", y="RNA_log2", data=df,
    color="k", size=2, alpha=0.4, jitter=0.2, ax=ax
)

# 8) add n-counts below each violin
ymin = ax.get_ylim()[0]
for i, g in enumerate(["active","inactive"]):
    ax.text(i, ymin - 0.5, f"n={counts.get(g,0)}",
            ha="center", va="top")

# 9) labels & title
ax.set(xlabel="", ylabel="log2(RNA+1)",
       title="Expression (RPKM>1) by chromatin context")

# 10) stats on raw RNA
act = df.loc[df["group"]=="active", "RNA"].dropna()
ina = df.loc[df["group"]=="inactive", "RNA"].dropna()
_, p = mannwhitneyu(act, ina, alternative="greater")
ax.text(0.5, 0.95, f"p = {p:.2e}",
        ha="center", va="center", transform=ax.transAxes)

plt.tight_layout()
plt.savefig(snakemake.output[0], dpi=150)
plt.close(fig)
