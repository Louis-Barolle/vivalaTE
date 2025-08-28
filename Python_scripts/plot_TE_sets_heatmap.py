#!/usr/bin/env python3
import re, numpy as np, pandas as pd
import seaborn as sns, matplotlib.pyplot as plt

# --- tiny helpers ---
def base_ids(s):
    # "TEID::chr:start-end" → "TEID"
    return s.astype(str).str.replace(r"::.*$", "", regex=True)

def short_te_id(te):
    # remove trailing "-end" in "..._start-end" → "..._start"
    te = str(te).split("::", 1)[0]
    m = re.search(r"_(\d+)-\d+$", te)
    return te[:m.start()] + f"_{m.group(1)}" if m else te

def zscore_cols(df):
    # column-wise Z: normalize each sample/column (mean=0, SD=1)
    if df.empty: return df
    mu = df.mean(axis=0)
    sd = df.std(axis=0).replace(0, 1)
    z = df.sub(mu, axis=1).div(sd, axis=1)
    return z.replace([np.inf, -np.inf], 0).fillna(0)

def draw_two_panel(mat_sub, out_png, title, rna_cols, chip_cols):
    rna  = mat_sub.loc[:, [c for c in mat_sub.columns if c in rna_cols]]
    chip = mat_sub.loc[:, [c for c in mat_sub.columns if c in chip_cols]]
    zr, zc = zscore_cols(rna), zscore_cols(chip)

    if zr.empty and zc.empty:
        fig, ax = plt.subplots(figsize=(6, 2)); ax.axis("off")
        ax.text(0.5, 0.5, "no rows after filtering", ha="center", va="center")
        fig.suptitle(title, y=0.98); fig.savefig(out_png, dpi=150, bbox_inches="tight")
        plt.close(fig); return

    nrows = (zr.shape[0] if not zr.empty else zc.shape[0])
    h = max(5, 0.18 * nrows)
    fig, axes = plt.subplots(1, 2, figsize=(14, h), gridspec_kw={"width_ratios":[1,1]})
    fig.suptitle(title, y=0.995, fontsize=12)

    # RNA panel
    ax = axes[0]
    if zr.empty:
        ax.axis("off"); ax.set_title("RNA (none)")
    else:
        vmax = np.nanpercentile(np.abs(zr.values), 98) or 1.0
        sns.heatmap(zr, ax=ax, cmap="vlag", center=0, vmin=-vmax, vmax=vmax,
                    cbar_kws={"shrink":0.6})
        ax.set_title("RNA (col-Z)"); ax.set_xlabel(""); ax.set_ylabel("TE copy")
        ax.set_xticklabels(zr.columns, rotation=90, fontsize=8)

    # ChIP panel
    ax = axes[1]
    if zc.empty:
        ax.axis("off"); ax.set_title("ChIP (none)")
    else:
        vmax = np.nanpercentile(np.abs(zc.values), 98) or 1.0
        sns.heatmap(zc, ax=ax, cmap="viridis", center=0, vmin=-vmax, vmax=vmax,
                    cbar_kws={"shrink":0.6})
        ax.set_title("Chromatin marks (col-Z)"); ax.set_xlabel(""); ax.set_yticklabels([])
        ax.set_xticklabels(zc.columns, rotation=90, fontsize=8)

    plt.tight_layout(rect=[0,0,1,0.98])
    fig.savefig(out_png, dpi=150, bbox_inches="tight"); plt.close(fig)

# --- load matrix (keep only the two genomes’ columns) ---
g1 = snakemake.wildcards.genome1
g2 = snakemake.wildcards.genome2

mat = pd.read_csv(snakemake.input.mat, sep="\t", index_col=0)
keep_cols = [c for c in mat.columns if c.startswith(g1 + "_") or c.startswith(g2 + "_")]
mat = mat.loc[:, keep_cols]

rna_cols  = [c for c in keep_cols if c.endswith("_RNA")]
chip_cols = [c for c in keep_cols if not c.endswith("_RNA")]

# --- build label map {TE_ID -> "short_id [ORIGIN•shared|unique]"} ---
label = {}

# RBH (shared): first col = g1 copy, second col = g2 copy
rbh = pd.read_csv(snakemake.input.rbh, sep="\t", header=None, usecols=[0,1])
for te in base_ids(rbh.iloc[:,0]):
    s = short_te_id(te); label[te] = f"{s} [{g1}•shared]"
for te in base_ids(rbh.iloc[:,1]):
    s = short_te_id(te); label[te] = f"{s} [{g2}•shared]"

# uniques
u1 = base_ids(pd.read_csv(snakemake.input.uniq1, sep="\t", header=None).iloc[:,0])
u2 = base_ids(pd.read_csv(snakemake.input.uniq2, sep="\t", header=None).iloc[:,0])
for te in u1:
    s = short_te_id(te); label.setdefault(te, f"{s} [{g1}•unique]")
for te in u2:
    s = short_te_id(te); label.setdefault(te, f"{s} [{g2}•unique]")

# --- subsets, apply labels by renaming index, plot ---
shared_ids = pd.Index(u for col in [rbh.iloc[:,0], rbh.iloc[:,1]] for u in base_ids(col)).unique()

sub_shared = mat.loc[mat.index.isin(shared_ids)].copy()
sub_shared.index = [label.get(i, short_te_id(i)) for i in sub_shared.index]
draw_two_panel(sub_shared, snakemake.output.shared,
               f"Shared TE copies (RBH) — {g1} ↔ {g2}", rna_cols, chip_cols)

sub_u1 = mat.loc[mat.index.isin(u1)].copy()
sub_u1.index = [label.get(i, short_te_id(i)) for i in sub_u1.index]
draw_two_panel(sub_u1, snakemake.output.unique1,
               f"{g1}-specific TE copies (vs {g2})", rna_cols, chip_cols)

sub_u2 = mat.loc[mat.index.isin(u2)].copy()
sub_u2.index = [label.get(i, short_te_id(i)) for i in sub_u2.index]
draw_two_panel(sub_u2, snakemake.output.unique2,
               f"{g2}-specific TE copies (vs {g1})", rna_cols, chip_cols)
