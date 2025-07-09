#!/usr/bin/env python3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

m = pd.read_csv(snakemake.input.matrix, sep="\t", index_col=0)

# log-transform to compress dynamic range
m["logRNA"] = np.log2(m["RNA"] + 1)
m["logChIP"] = np.log2(m["H3K27ac"] + 1)

plt.figure(figsize=(5,5))
sns.regplot(x="logChIP", y="logRNA", data=m,
            scatter_kws={"s":10, "alpha":0.5}, line_kws={"color":"red"})
plt.xlabel("log2(mean H3K27ac + 1)")
plt.ylabel("log2(RNA + 1)")
plt.title(f"RNA vs. H3K27ac — {snakemake.wildcards.tissue}")
plt.tight_layout()
plt.savefig(snakemake.output.png, dpi=150)
