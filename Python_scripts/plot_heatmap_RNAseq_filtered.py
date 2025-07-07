#!/usr/bin/env python3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# load & coerce
m = (pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
       .apply(pd.to_numeric, errors="coerce"))

# filter by RNA ≥ 1
f = m[m["RNA"] >= 1]


# **exclude TEs beginning with "micropia"**
f = f[~f.index.str.startswith("Micropia")]

# sort rows so highest‐RNA TEs are at the top
f = f.sort_values("RNA", ascending=False)
# sort rows so highest‐RNA TEs are at the top
f = f.sort_values("RNA", ascending=False)

# save TE list
f.index.to_series().to_csv(snakemake.output.te_list, header=False)

# split ChIP vs RNA
chips   = [c for c in f.columns if c != "RNA"]
chip_mat = f[chips]
rna_mat  = f[["RNA"]]

# sizing
h = max(6, len(f) * 0.3)
w = max(6, len(chips) * 0.5 + 2)

fig, (ax1, ax2) = plt.subplots(
    1, 2, figsize=(w, h), sharey=True,
    gridspec_kw={'width_ratios':[len(chips),1], 'wspace':0.4},
    constrained_layout=True
)

# ChIP as before
sns.heatmap(chip_mat, cmap="viridis", ax=ax1, cbar=True)
ax1.set(title="TEs with RNA≥1 – ChIP", xlabel="Samples", ylabel="TE Names")
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, ha="right")

# RNA with log scale
sns.heatmap(
    rna_mat,
    cmap="Reds",
    ax=ax2,
    cbar=True,
    norm=LogNorm(vmin=max(rna_mat.min().min(), 0.1),  # avoid log(0)
                 vmax=rna_mat.max().max())
)
ax2.set(title="TEs with RNA≥1 – RNA (log scale)", xlabel="RNA", ylabel="")
ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45, ha="right")

plt.savefig(snakemake.output.heatmap, dpi=150)
plt.close(fig)