#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns

# 1) load & annotate
df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
df["Expressed"] = df["RNA"] > 1
df["Location"] = df["context"]
df["State"]    = (df["H3K27ac_flank"] > df["H3K9me3_flank"]) \
                    .map({True: "H3K27ac", False: "H3K9me3"})

# 2) dynamic categories & counts
cats_loc   = df["Location"].value_counts().index.tolist()
cats_state = df["State"].value_counts().index.tolist()

ctx_all   = df["Location"].value_counts().reindex(cats_loc,   fill_value=0)
ctx_expr  = df[df["Expressed"]]["Location"]\
                 .value_counts().reindex(cats_loc, fill_value=0)

st_all    = df["State"].value_counts().reindex(cats_state,   fill_value=0)
st_expr   = df[df["Expressed"]]["State"]\
                 .value_counts().reindex(cats_state, fill_value=0)

pal_loc   = sns.color_palette("pastel",   len(cats_loc))
pal_state = sns.color_palette("muted",    len(cats_state))

# 3) plot 2×2 pies
fig, axes = plt.subplots(2, 2, figsize=(5,5))
for ax, data, cats, pal, title in [
    (axes[0,0], ctx_all,   cats_loc,   pal_loc,   "Context (all)"),
    (axes[0,1], ctx_expr,  cats_loc,   pal_loc,   "Context (expr)"),
    (axes[1,0], st_all,    cats_state, pal_state, "Chromatin (all)"),
    (axes[1,1], st_expr,   cats_state, pal_state, "Chromatin (expr)"),
]:
    ax.pie(data, labels=cats, colors=pal,
           startangle=90, autopct="%1.1f%%", pctdistance=0.7,
           textprops={"fontsize":8})
    ax.set_title(title, fontsize=10)
    ax.axis("equal")

# 4) legends
ctx_patches   = [Patch(color=pal_loc[i],   label=c) for i,c in enumerate(cats_loc)]
state_patches = [Patch(color=pal_state[i], label=c) for i,c in enumerate(cats_state)]

fig.legend(handles=ctx_patches,   title="Genomic context",
           loc="center right", bbox_to_anchor=(1.3,0.75), fontsize=8)
fig.legend(handles=state_patches, title="Dominant upstream mark",
           loc="center right", bbox_to_anchor=(1.3,0.35), fontsize=8)

plt.tight_layout()
plt.savefig(snakemake.output[0], dpi=150, bbox_inches="tight")
