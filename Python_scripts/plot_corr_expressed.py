#!/usr/bin/env python3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Load matrix
m = pd.read_csv(snakemake.input.matrix, sep="\t", index_col=0)

# (Optional) report how many loci have both RNA==0 and H3K27ac==0
n_zero_zero = ((m["RNA"] == 0) & (m["H3K27ac"] == 0)).sum()
print(f"Number of zero–zero pairs: {n_zero_zero}")

# Set up figure with two panels
fig = plt.figure(figsize=(12,6))
gs  = gridspec.GridSpec(1, 2, width_ratios=[1,1], wspace=0.3)

# Panel 1: all TEs
ax0 = fig.add_subplot(gs[0])
corr_all = m.corr(method="spearman")
sns.heatmap(
    corr_all,
    cmap="RdBu_r",
    center=0,
    vmin=-1, vmax=1,
    cbar_kws={"shrink": .6},
    ax=ax0
)
ax0.set_title("Spearman ρ: All TEs")

# Panel 2: expressed only
ax1 = fig.add_subplot(gs[1])
m_expr = m[m["RNA"] >= 1]
corr_expr = m_expr.corr(method="spearman")
sns.heatmap(
    corr_expr,
    cmap="RdBu_r",
    center=0,
    vmin=-1, vmax=1,
    cbar_kws={"shrink": .6},
    ax=ax1
)
ax1.set_title("Spearman ρ: Expressed TEs (RNA ≥ 1)")

# Save
plt.tight_layout()
plt.savefig(snakemake.output.png, dpi=300)
