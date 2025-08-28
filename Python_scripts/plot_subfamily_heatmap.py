#!/usr/bin/env python3
import pandas as pd, numpy as np
import seaborn as sns, matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

mat = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)

# split columns
rna_cols  = [c for c in mat.columns if c.endswith("RNA")]
chip_cols = [c for c in mat.columns if c not in rna_cols]

def zscore_rows(df):
    z = df.sub(df.mean(axis=1), axis=0)
    denom = df.std(axis=1).replace(0, np.nan)
    z = z.div(denom, axis=0)
    return z.replace([np.inf, -np.inf], np.nan).fillna(0)

def order_for_panels(z_rna, z_chip):
    # common row order from joint clustering (keeps panels aligned)
    both = pd.concat([z for z in [z_rna, z_chip] if not z.empty], axis=1)
    if both.shape[0] >= 2:
        row_order = leaves_list(linkage(pdist(both.values, metric="euclidean"),
                                        method="average"))
    else:
        row_order = np.arange(both.shape[0])

    def col_order(z):
        if z.shape[1] >= 2:
            return leaves_list(linkage(pdist(z.T.values, metric="euclidean"),
                                       method="average"))
        else:
            return np.arange(z.shape[1])

    return row_order, col_order(z_rna), col_order(z_chip)

def draw_two_panel(m, out, title, active_only=False):
    if active_only and rna_cols:
        mask = m[rna_cols].gt(1).any(axis=1)
        m = m.loc[mask]

    rna = m[rna_cols]  if rna_cols  else pd.DataFrame(index=m.index)
    chip = m[chip_cols] if chip_cols else pd.DataFrame(index=m.index)

    if rna.empty and chip.empty:
        fig, ax = plt.subplots(figsize=(6, 2))
        ax.axis("off"); ax.text(0.5, 0.5, "No data", ha="center", va="center")
        fig.suptitle(title)
        fig.savefig(out, dpi=150, bbox_inches="tight"); plt.close()
        return

    zr, zc = zscore_rows(rna) if not rna.empty else rna, \
             zscore_rows(chip) if not chip.empty else chip

    # align row order across panels
    if zr.empty:
        rows = np.arange(zc.shape[0])
    elif zc.empty:
        rows = np.arange(zr.shape[0])
    else:
        rows, cr, cc = order_for_panels(zr, zc)
        zr = zr.iloc[rows, :].iloc[:, cr] if not zr.empty else zr
        zc = zc.iloc[rows, :].iloc[:, cc] if not zc.empty else zc

    # sizes
    nrows = zr.shape[0] if not zr.empty else zc.shape[0]
    h = max(6, 0.22 * nrows)
    fig, axes = plt.subplots(1, 2, figsize=(14, h), gridspec_kw={"width_ratios":[1,1]})
    fig.suptitle(title, y=0.995, fontsize=12)

    # RNA panel
    ax = axes[0]
    if zr.empty:
        ax.axis("off"); ax.set_title("RNA (no columns)")
    else:
        vmax = np.nanpercentile(np.abs(zr.values), 98) or 1.0
        sns.heatmap(zr, ax=ax, cmap="vlag", center=0, vmin=-vmax, vmax=vmax,
                    cbar_kws={"shrink":0.6})
        ax.set_title("RNA (row-Z)")
        ax.set_xlabel("")
        ax.set_ylabel("TE subfamily")
        ax.set_xticklabels(zr.columns, rotation=90, fontsize=8)

    # Histone panel
    ax = axes[1]
    if zc.empty:
        ax.axis("off"); ax.set_title("Chromatin marks (no columns)")
    else:
        vmax = np.nanpercentile(np.abs(zc.values), 98) or 1.0
        sns.heatmap(zc, ax=ax, cmap="viridis", center=0, vmin=-vmax, vmax=vmax,
                    cbar_kws={"shrink":0.6})
        ax.set_title("Chromatin marks (row-Z)")
        ax.set_xlabel("")
        ax.set_yticklabels([])  # rows already labeled on the left
        ax.set_xticklabels(zc.columns, rotation=90, fontsize=8)

    plt.tight_layout(rect=[0,0,1,0.98])
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)

# (1) all rows
draw_two_panel(mat, snakemake.output[0],
               "TE subfamilies across genomes/tissues — RNA vs histone (separate scales)")

# (2) active only (RNA>1 anywhere)
draw_two_panel(mat, snakemake.output[1],
               "Expressed TE subfamilies (RNA>1) — RNA vs histone (separate scales)",
               active_only=True)
