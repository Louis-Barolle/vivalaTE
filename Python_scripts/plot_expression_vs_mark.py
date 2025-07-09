#!/usr/bin/env python3
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, fisher_exact

# 1) load summary table
df = pd.read_csv(snakemake.input.summary, sep="\t", index_col=0)

# 2) pick upstream‐flank columns by name
act_col = next(c for c in df.columns if "H3K27ac" in c and c.endswith("_flank"))
rep_col = next(c for c in df.columns if "H3K9me3" in c and c.endswith("_flank"))

# 3) classify each TE
df["expressed"] = df["RNA"] > 1
# make a string category so palette keys match
df["expr_cat"] = df["expressed"].map({False: "not expr.", True: "expr."})
df["has_act"]  = df[act_col] > df[rep_col]
df["has_rep"]  = df[rep_col] > df[act_col]
cats = ["not expr.", "expr."]
# 4) stats
u_act, p_act_cont = mannwhitneyu(
    df.loc[df.expressed, act_col],
    df.loc[~df.expressed, act_col],
    alternative="greater"
)
tbl_act = pd.crosstab(df.expressed, df.has_act)
odds_act, p_act_bin = fisher_exact(tbl_act, alternative="greater")

u_rep, p_rep_cont = mannwhitneyu(
    df.loc[df.expressed, rep_col],
    df.loc[~df.expressed, rep_col],
    alternative="less"
)
tbl_rep = pd.crosstab(df.expressed, df.has_rep)
odds_rep, p_rep_bin = fisher_exact(tbl_rep, alternative="less")

# 5) set up 2×2 figure
fig, axes = plt.subplots(2, 2, figsize=(10, 8))
(ax1, ax2), (ax3, ax4) = axes

# helper to draw violin + strip + n
def violin_with_strip(ax, col, title, pval, palette):
    sns.violinplot(
        x="expr_cat", y=col, data=df, hue="expr_cat", legend=False,
        palette=palette, inner="quartile", cut=0, ax=ax
    )
    if ax.get_legend():
        ax.get_legend().remove()
    subs = df.groupby("expr_cat", group_keys=False).apply(
        lambda g: g.sample(min(len(g), 200), random_state=0)
    )
    sns.stripplot(
        x="expr_cat", y=col, data=subs,
        color="k", size=2, alpha=0.4, jitter=0.2, ax=ax
    )
    counts = df["expr_cat"].value_counts()
    y0 = ax.get_ylim()[0]
    for i, cat in enumerate(["not expr.", "expr."]):
        ax.text(i, y0 - (abs(y0)*0.05), f"n={counts.get(cat,0)}",
                ha="center", va="top")
    ax.set_xticks([0,1])
    ax.set_xticklabels(["not expr.","expr."])
    ax.set_ylabel(col)
    ax.set_title(f"{title}\n{pval}")

# colors keyed by our new strings:
pal = {"not expr.":"#66c2a5", "expr.":"#fc8d62"}

# top‐left: H3K27ac intensity
violin_with_strip(
    ax1, act_col,
    title="Active mark (H3K27ac)",
    pval=f"M-W p={p_act_cont:.1e}",
    palette=pal 
)

# top‐right: H3K27ac % w/ mark
freq_act = df.groupby("expr_cat")["has_act"].mean() * 100
vals_act = [freq_act.get(c, 0) for c in cats]
ax2.bar([0,1], vals_act,
        color=[pal[c] for c in cats],
        width=0.6)
for i, pct in enumerate(vals_act):
    ax2.text(i, pct + 2, f"{pct:.1f}%", ha="center")
ax2.set_xticks([0,1])
ax2.set_xticklabels(cats)
ax2.set_ylabel("% w/ H3K27ac")
ax2.set_title(f"Active mark freq\nOR={odds_act:.2f}, p={p_act_bin:.1e}")
# bottom‐left: H3K9me3 intensity
violin_with_strip(
    ax3, rep_col,
    title="Repressive mark (H3K9me3)",
    pval=f"M-W p={p_rep_cont:.1e}",
    palette=pal
)

# bottom‐right: H3K9me3 % w/ mark
freq_rep = df.groupby("expr_cat")["has_rep"].mean() * 100
vals_rep = [freq_rep.get(c, 0) for c in cats]
ax4.bar([0,1], vals_rep,
        color=[pal[c] for c in cats],
        width=0.6)
for i, pct in enumerate(vals_rep):
    ax4.text(i, pct + 2, f"{pct:.1f}%", ha="center")
ax4.set_xticks([0,1])
ax4.set_xticklabels(cats)
ax4.set_ylabel("% w/ H3K9me3")
ax4.set_title(f"Repressive mark freq\nOR={odds_rep:.2f}, p={p_rep_bin:.1e}")
plt.tight_layout()
plt.savefig(snakemake.output.png, dpi=150)
plt.close(fig)
