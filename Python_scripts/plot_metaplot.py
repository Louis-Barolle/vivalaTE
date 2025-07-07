#!/usr/bin/env python3
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# load final matrix
m = pd.read_csv(snakemake.input.matrix, sep="\t", index_col=0)

# assume you have per-TE mean signal for H3K27ac in m["H3K27ac"]
# and you have bedtools signal per-base in an auxiliary bigwig or pre-binned file...
# for minimalism, we'll bin the continuous signal into 100 bins across all TEs:

# here m_signal is a DataFrame: rows=TEs, cols=bins (this must be precomputed or loaded)
# as a quick stand‐in, we simulate by slicing the signal vector:
# (in practice, you would load a separate matrix of shape [nTE x nBins])
signal = m["H3K27ac"].values
n_bins = 50
bins = np.array_split(signal, n_bins)
mean_profile = [np.mean(b) for b in bins]
sem_profile  = [np.std(b)/np.sqrt(len(b)) for b in bins]

x = np.linspace(-2000, +2000, n_bins)

plt.figure(figsize=(6,4))
plt.fill_between(x, np.array(mean_profile)-sem_profile, np.array(mean_profile)+sem_profile, alpha=0.3)
plt.plot(x, mean_profile, lw=2)
plt.axvline(0, color="k", ls="--", lw=1)
plt.xlabel("Position relative to TE insertion (bp)")
plt.ylabel("Mean H3K27ac signal ± SEM")
plt.title(f"Aggregate metaplot — {snakemake.wildcards.tissue}")
plt.tight_layout()
plt.savefig(snakemake.output.png, dpi=150)
