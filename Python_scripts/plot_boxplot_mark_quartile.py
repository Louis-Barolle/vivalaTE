#!/usr/bin/env python3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# 1) load the combined matrix
m = pd.read_csv(snakemake.input.matrix, sep="\t", index_col=0)

# 2) pick your mark column
mark_col = "H3K27ac"
if mark_col not in m.columns:
    raise ValueError(f"Expected column {mark_col!r} in {snakemake.input.matrix}")

# 3) drop the top 1% of values to avoid extreme outliers
p99 = m[mark_col].quantile(0.99)
m = m[m[mark_col] <= p99]

# 4) cut into quartiles (on the filtered data)
m["quartile"] = pd.qcut(m[mark_col], 4, labels=["Q1","Q2","Q3","Q4"])

# 5) compute percent expressed (RNA ≥ 1) in each quartile
expr_frac = (m["RNA"] >= 1).groupby(m["quartile"]).mean() * 100

# 6) plot
plt.figure(figsize=(6,4))
sns.barplot(x=expr_frac.index, y=expr_frac.values, palette="Blues_d")
plt.ylabel("Percent expressed (RNA ≥ 1)")
plt.xlabel(f"{mark_col} upstream quartile (≤99th pct)")
plt.ylim(0,100)
plt.title("TE expression rate by H3K27ac quartile\n(top 1% of values removed)")
plt.tight_layout()

# 7) save
plt.savefig(snakemake.output[0], dpi=300)
