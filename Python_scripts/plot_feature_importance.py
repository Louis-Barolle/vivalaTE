#!/usr/bin/env python3
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
import matplotlib.pyplot as plt

try:
    matrix_fp = snakemake.input.matrix
except AttributeError:
    matrix_fp = snakemake.input[0]

m = pd.read_csv(matrix_fp, sep="\t", index_col=0)

X = m.drop(columns=["RNA"], errors="ignore")
y = m["RNA"]

if X.shape[1] == 0:
    fig, ax = plt.subplots(figsize=(4,2))
    ax.text(0.5, 0.5, "No ChIP/RNA marks found\nskipping", 
            ha="center", va="center", fontsize=12)
    ax.axis("off")
    out_fp = getattr(snakemake.output, "png", snakemake.output[0])
    plt.savefig(out_fp, dpi=150, bbox_inches="tight")
    exit(0)


rf = RandomForestRegressor(n_estimators=100, random_state=0)
rf.fit(X, y)
imp = pd.Series(rf.feature_importances_, index=X.columns)
top10 = imp.sort_values(ascending=False).head(10)

fig, ax = plt.subplots(figsize=(6,4))
top10.sort_values().plot.barh(color="steelblue", ax=ax)
ax.set_xlabel("Feature importance")
tissue = getattr(snakemake.wildcards, "tissue", "")
ax.set_title(f"Top 10 predictors of RNA — {tissue}")
plt.tight_layout()


out_fp = getattr(snakemake.output, "png", snakemake.output[0])
plt.savefig(out_fp, dpi=150, bbox_inches="tight")
