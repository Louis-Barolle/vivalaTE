#!/usr/bin/env python3
import os, sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# io
mat = pd.read_csv(snakemake.input.mat, sep="\t", index_col=0)

# gather unique TE IDs (skip missing/empty files silently)
uniq_ids = []
for u in snakemake.input.uniq:
    try:
        if os.path.getsize(u) > 0:
            uniq_ids.extend(pd.read_csv(u, header=None).iloc[:,0].astype(str).tolist())
    except FileNotFoundError:
        pass
uniq_ids = pd.Index(pd.unique(uniq_ids))

genome = snakemake.wildcards.genome
out_png = snakemake.output[0]

def placeholder(msg):
    fig, ax = plt.subplots(figsize=(6, 3))
    ax.axis("off")
    ax.text(0.5, 0.5, msg, ha="center", va="center", fontsize=12)
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    sys.exit(0)

if len(uniq_ids) == 0:
    placeholder(f"No unique TE insertions detected in {genome}")

sub = mat.loc[mat.index.intersection(uniq_ids)]
if sub.empty:
    placeholder(f"No unique TE insertions detected in {genome}")

# split columns (keep original order)
def is_rna(c):
    return c.endswith("_RNA") or c.split("_")[-1] == "RNA"

rna_cols  = [c for c in sub.columns if is_rna(c)]
chip_cols = [c for c in sub.columns if c not in rna_cols]

# z-score by row (avoid inf/nan)
def zscore_rows(df):
    if df.empty: return df
    z = df.sub(df.mean(axis=1), axis=0)
    denom = df.std(axis=1).replace(0, np.nan)
    z = z.div(denom, axis=0)
    return z.replace([np.inf, -np.inf], np.nan).fillna(0)

zr = zscore_rows(sub[rna_cols])   if rna_cols  else pd.DataFrame(index=sub.index)
zc = zscore_rows(sub[chip_cols])  if chip_cols else pd.DataFrame(index=sub.index)



# rename rows to include origin + shortened ID
import re
def short_te_id(te):
    m = re.search(r"_(\d+)-\d+$", te)
    return te[:m.start()] + f"_{m.group(1)}" if m else te

sub = sub.copy()
sub.index = [f"{short_te_id(i)} [{wildcards.genome}•unique]" for i in sub.index]
# simple plotting (no clustering), same row order in both panels
MAX_ROWS_LABELS = 400
HEIGHT_PER_ROW  = 0.16

nrows = sub.shape[0]
h = min(22, max(6, HEIGHT_PER_ROW * nrows))
show_rows = nrows <= MAX_ROWS_LABELS

fig, axes = plt.subplots(1, 2, figsize=(14, h), gridspec_kw={"width_ratios":[1,1]})
fig.suptitle(f"Unique TE insertions in {genome} — RNA vs chromatin (row-Z, no clustering)", y=0.995, fontsize=12)

# RNA
ax = axes[0]
if zr.empty:
    ax.axis("off"); ax.set_title("RNA (none)")
else:
    vmax_r = float(np.nanpercentile(np.abs(zr.values), 98)) or 1.0
    sns.heatmap(zr, ax=ax, cmap="vlag", center=0, vmin=-vmax_r, vmax=vmax_r,
                cbar_kws={"shrink":0.6}, yticklabels=show_rows, xticklabels=True)
    ax.set_title("RNA (row-Z)")
    ax.set_xlabel("")
    ax.set_ylabel("TE copy" if show_rows else "")
    if not show_rows: ax.set_yticklabels([])
    ax.set_xticklabels(zr.columns, rotation=90, fontsize=8)

# Chromatin
ax = axes[1]
if zc.empty:
    ax.axis("off"); ax.set_title("Chromatin marks (none)")
else:
    vmax_c = float(np.nanpercentile(np.abs(zc.values), 98)) or 1.0
    sns.heatmap(zc, ax=ax, cmap="viridis", center=0, vmin=-vmax_c, vmax=vmax_c,
                cbar_kws={"shrink":0.6}, yticklabels=False, xticklabels=True)
    ax.set_title("Chromatin marks (row-Z)")
    ax.set_xlabel("")
    ax.set_xticklabels(zc.columns, rotation=90, fontsize=8)

plt.tight_layout(rect=[0,0,1,0.98])
fig.savefig(out_png, dpi=150, bbox_inches="tight")
plt.close(fig)
