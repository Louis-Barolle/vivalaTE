#!/usr/bin/env python3
import os, pandas as pd, seaborn as sns, matplotlib.pyplot as plt

# load & coerce
m = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)\
      .apply(pd.to_numeric, errors="coerce")

chips = [c for c in m.columns if c != "RNA"]
fams  = m.index.to_series().str.split("_").str[0]
out   = snakemake.output[0]; os.makedirs(out, exist_ok=True)

for fam in sorted(fams.unique()):
    sub = m[fams == fam]
    if sub.empty: continue
    cs, rs = sub[chips], sub[["RNA"]]
    h = max(6, sub.shape[0] * 0.2)
    w = max(6, len(chips) * 1.0 + 2)
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(w, h), sharey=True,
        gridspec_kw={'width_ratios':[len(chips),1],'wspace':0.4},
        constrained_layout=True
    )
    sns.heatmap(cs, cmap="viridis", ax=ax1, cbar=True, cbar_kws={'pad':0.02})
    ax1.set(title=f"{fam} – ChIP upstream", xlabel="Samples")
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, ha="right")
    ax1.set_ylabel("TE Names")
    vmax = rs["RNA"].max()
    sns.heatmap(rs, cmap="Reds", ax=ax2, cbar=True,
                vmin=0, vmax=vmax,
                cbar_kws={'ticks':[0, vmax/2, vmax],'pad':0.02})
    ax2.set(title=f"{fam} – RNA-seq inside of TE", xlabel="RNA signal")
    fig.savefig(f"{out}/heatmap_{fam}.png")
    plt.close(fig)
