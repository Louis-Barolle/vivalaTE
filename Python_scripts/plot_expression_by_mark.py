#!/usr/bin/env python3
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)

# choose flank columns
h3a = next(c for c in df if c.startswith("H3K27ac") and c.endswith("_flank"))
h3r = next(c for c in df if c.startswith("H3K9me3") and c.endswith("_flank"))

# classify
df["group"] = df.apply(lambda r: "active" if r[h3a] > r[h3r] else "inactive", axis=1)

# log transform
df["RNA_log2"] = np.log2(df["RNA"].fillna(0) + 1)

# counts
counts = df["group"].value_counts()

# plot
fig, ax = plt.subplots(figsize=(5,4))
sns.violinplot(x="group", y="RNA_log2", data=df,hue="group",legend=False, palette={"active":"#fc8d62","inactive":"#66c2a5"}, inner="quartile", cut=0, ax=ax)
sns.stripplot(x="group", y="RNA_log2", data=df, color="k", size=2, alpha=0.4, jitter=0.2, ax=ax)
for i, g in enumerate(["active","inactive"]):
    ax.text(i, ax.get_ylim()[0] - 1, f"n={counts.get(g,0)}", ha="center")
ax.set(xlabel="", ylabel="log2(RNA+1)", title="Expression by chromatin context - gut")

# stats on raw RNA
act = df[df["group"]=="active"]["RNA"].dropna()
ina = df[df["group"]=="inactive"]["RNA"].dropna()
_, p = mannwhitneyu(act, ina, alternative="greater")
ax.text(0.5, 0.95, f"p = {p:.2e}", ha="center", va="center", transform=ax.transAxes)

plt.tight_layout()
plt.subplots_adjust(bottom=0.25)
plt.savefig(snakemake.output[0], dpi=150, bbox_inches="tight")
plt.close(fig)
