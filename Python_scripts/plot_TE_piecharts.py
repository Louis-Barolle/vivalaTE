#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# 1) load & annotate
df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
df["Expressed"] = df["RNA"] > 1
df["Location"] = df["context"]               # "genic" / "intergenic"
df["State"] = (df["H3K27ac_flank"] > df["H3K9me3_flank"]) \
                 .map({True: "H3K27ac", False: "H3K9me3"})

# 2) counts (ensure fixed order)
ctx_all  = df["Location"].value_counts().reindex(["genic","intergenic"], fill_value=0)
ctx_expr = df[df["Expressed"]]["Location"] \
              .value_counts().reindex(["genic","intergenic"], fill_value=0)
st_all   = df["State"].value_counts().reindex(["H3K27ac","H3K9me3"], fill_value=0)
st_expr  = df[df["Expressed"]]["State"] \
              .value_counts().reindex(["H3K27ac","H3K9me3"], fill_value=0)

# 3) plot grid
fig, axes = plt.subplots(2,2, figsize=(5,5))
pal_ctx   = ["#1f78b4", "#a6cee3"]
pal_state = ["#fc8d62", "#66c2a5"]

# draw pies with percentages inside
for ax, data, colors, title in [
    (axes[0,0], ctx_all,  pal_ctx,   "Context (all)"),
    (axes[0,1], ctx_expr, pal_ctx,   "Context (expr)"),
    (axes[1,0], st_all,   pal_state, "Chromatin (all)"),
    (axes[1,1], st_expr,  pal_state, "Chromatin (expr)"),
]:
    ax.pie(
        data,
        colors=colors,
        startangle=90,
        autopct="%1.1f%%",
        pctdistance=0.7,
        textprops={"fontsize": 8},
    )
    ax.set_title(title, fontsize=10)
    ax.axis("equal")

# 4) legends
ctx_labels   = ["genic", "intergenic"]
state_labels = ["H3K27ac", "H3K9me3"]

ctx_patches   = [Patch(color=pal_ctx[i],   label=ctx_labels[i])   for i in range(2)]
state_patches = [Patch(color=pal_state[i], label=state_labels[i]) for i in range(2)]

# place two legends to the right
fig.legend(
    handles=ctx_patches,
    title="Genomic context",
    loc="center right",
    bbox_to_anchor=(1.3, 0.75),
    fontsize=8,
)
fig.legend(
    handles=state_patches,
    title="Dominant upstream mark",
    loc="center right",
    bbox_to_anchor=(1.3, 0.35),
    fontsize=8,
)

plt.tight_layout()
plt.savefig(snakemake.output[0], dpi=150, bbox_inches="tight")
