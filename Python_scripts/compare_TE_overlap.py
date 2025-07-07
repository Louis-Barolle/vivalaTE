#!/usr/bin/env python3
import os, pandas as pd, matplotlib.pyplot as plt, seaborn as sns
from matplotlib_venn import venn2

os.makedirs(os.path.dirname(snakemake.output.overlap), exist_ok=True)

gut  = pd.read_csv(snakemake.input.gut,  sep="\t", index_col=0)
head = pd.read_csv(snakemake.input.head, sep="\t", index_col=0)
g = set(gut.index[gut["RNA"] >= 1])
h = set(head.index[head["RNA"] >= 1])
common = sorted(g & h)
pd.Series(common).to_csv(snakemake.output.overlap, index=False, header=False)

# Venn
fig, ax = plt.subplots(figsize=(4,4))
venn2([g,h], ("Gut","Head"), ax=ax)
plt.savefig(snakemake.output.venn, dpi=150)
plt.close(fig)

# Composition
dfc = pd.read_csv(snakemake.input.cls, sep="\t",
                  usecols=["full_name","repeatmasker_type","repeatmasker_subtype"]
                 ).rename(columns={"full_name":"Family","repeatmasker_type":"Class","repeatmasker_subtype":"Subclass"})
df = (pd.DataFrame({"TE":common})
      .assign(Family=lambda d: d.TE.str.split(pat="_", n=1).str[0])
      .merge(dfc, on="Family", how="left")
      .fillna("Unknown")
)
comps = {lvl: df[lvl].value_counts(normalize=True) for lvl in ["Family","Class","Subclass"]}

fig, axes = plt.subplots(1,3,figsize=(15,5))
for ax,lvl in zip(axes, comps):
    sns.barplot(x=comps[lvl].values, y=comps[lvl].index, ax=ax, orient="h", palette="tab20")
    ax.set_title(lvl); ax.set_xlabel("Prop"); ax.xaxis.set_major_formatter(lambda x, _:f"{x:.0%}")
plt.tight_layout()
plt.savefig(snakemake.output.composition, dpi=150)
