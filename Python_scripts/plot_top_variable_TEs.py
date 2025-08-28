#!/usr/bin/env python3
import pandas as pd, seaborn as sns, matplotlib.pyplot as plt

mat = pd.read_csv(snakemake.input.mat, sep="\t", index_col=0)

N = 50
vars = mat.var(axis=1)
top = vars.nlargest(N).index
sub = mat.loc[top]

# z-score (avoid div-by-zero with .replace(0,1))
z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, 1), axis=0)

# NaN-proof for clustering
z = z.dropna(how="all", axis=0).dropna(how="all", axis=1).fillna(0)

h = max(6, 0.2 * z.shape[0])
g = sns.clustermap(
    z, cmap="vlag", row_cluster=True, col_cluster=True,
    figsize=(10, h), xticklabels=True, yticklabels=True
)
g.fig.suptitle(f"Top {N} most variable TE insertions", fontsize=12, y=1.02)
plt.savefig(snakemake.output[0], dpi=150, bbox_inches="tight")
plt.close()
