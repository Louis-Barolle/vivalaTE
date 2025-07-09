#!/usr/bin/env python3
import pandas as pd, seaborn as sns, matplotlib.pyplot as plt

# load & coerce
m = (pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
       .apply(pd.to_numeric, errors="coerce")
       .dropna(how="all"))

# split ChIP vs RNA
chips  = [c for c in m.columns if c != "RNA"]
chip   = m[chips]
rna    = m[["RNA"]]

# size
h = max(6, m.shape[0] * 0.2)
w = max(6, len(chips) * 1.0 + 2)

fig, (ax1, ax2) = plt.subplots(
    1, 2, figsize=(w, h), sharey=True,
    gridspec_kw={'width_ratios':[len(chips),1], 'wspace':0.4},
    constrained_layout=True
)

# ChIP panel
sns.heatmap(chip, cmap="viridis", ax=ax1, cbar=True, cbar_kws={'pad':0.02})
ax1.set(title="ChIP signals", xlabel="Samples", ylabel="TE Names")
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, ha="right")

# RNA panel
vmax = rna["RNA"].max()
sns.heatmap(rna, cmap="Reds", ax=ax2, cbar=True,
            vmin=0, vmax=vmax,
            cbar_kws={'ticks':[0, vmax/2, vmax], 'pad':0.02})
ax2.set(title="RNA-seq", xlabel="RNA signal", ylabel="")

fig.savefig(snakemake.output[0])
plt.close(fig)
