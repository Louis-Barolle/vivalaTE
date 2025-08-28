#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt

mat = pd.read_csv(snakemake.input.mat, sep="\t", index_col=0)
te = snakemake.params.TE_ID

if te not in mat.index:
    raise ValueError(f"TE {te} not found in matrix")

vals = mat.loc[te]
# clean up column labels (optional)
vals.index = vals.index.str.replace("_flank", "")

plt.figure(figsize=(6,4))
vals.plot.bar(color="C2")
plt.ylabel("Mean signal / RNA")
plt.xticks(rotation=45, ha="right")
plt.title(f"TE copy {te} across genomes & tissues")
plt.tight_layout()
plt.savefig(snakemake.output[0], dpi=150)
