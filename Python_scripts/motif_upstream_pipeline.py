#!/usr/bin/env python3
import os
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# step 1: extract upstream FASTA
os.makedirs(os.path.dirname(snakemake.output.fimo), exist_ok=True)
subprocess.run(
    f"bedtools getfasta -fi {snakemake.input.genome} "
    f"-bed {snakemake.input.bed} -fo {snakemake.output.fasta}",
    shell=True, check=True
)

# step 2: run FIMO
outdir = os.path.dirname(snakemake.output.fimo)
subprocess.run(
    f"fimo --thresh 1e-4 --oc {outdir} "
    f"{snakemake.input.motifs} {snakemake.output.fasta}",
    shell=True, check=True
)
os.replace(f"{outdir}/fimo.tsv", snakemake.output.fimo)

# step 3: plot motif‐hit density
df = pd.read_csv(snakemake.output.fimo, sep="\t", comment="#")
region_len, n_bins = 5000, 100
edges = np.linspace(1, region_len, n_bins + 1)
counts, _ = np.histogram(df["start"], bins=edges)
density = counts / df["sequence_name"].nunique()
centers = (edges[:-1] + edges[1:]) / 2

plt.figure(figsize=(6,4))
sns.lineplot(x=centers, y=density, color="steelblue")
plt.xlabel("Position upstream of TE (bp)")
plt.ylabel("Avg motif hits per TE")
plt.title("Motif density across 5 kb upstream")
plt.tight_layout()
plt.savefig(snakemake.output.density, dpi=150)
