#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 1. Load expression matrix and filter for expressed TEs
m = pd.read_csv(snakemake.input.matrix, sep="\t", index_col=0)
m_expr = m[m["RNA"] >= 1].reset_index().rename(columns={"index": "TE"})

# 2. Derive Family from TE name (assumes "Family_chr-start-end" format)
m_expr["Family"] = m_expr["TE"].str.split(pat="_", n=1).str[0]

# 3. Load your class/subclass table, selecting only the needed columns
df_classes = pd.read_csv(
    snakemake.input.classes,
    sep="\t",
    usecols=["full_name", "repeatmasker_type", "repeatmasker_subtype"],
    dtype=str
).rename(columns={
    "full_name": "Family",
    "repeatmasker_type": "Class",
    "repeatmasker_subtype": "Subclass"
})

# 4. Merge annotations on Family
df = m_expr.merge(df_classes, on="Family", how="left")



# Fill any missing with 'Unknown' to guarantee plotting
df["Class"]    = df["Class"].fillna("Unknown")
df["Subclass"] = df["Subclass"].fillna("Unknown")


# 5. Compute proportional counts
levels = ["Family", "Class", "Subclass"]
counts = {
    lvl: df[lvl].value_counts(normalize=True).sort_values(ascending=False)
    for lvl in levels
}

# 6. Plot
fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharex=False)
for ax, lvl in zip(axes, levels):
    sns.barplot(
        x=counts[lvl].values,
        y=counts[lvl].index,
        ax=ax,
        palette="tab20",
        orient="h"
    )
    ax.set_title(f"Expressed TEs by {lvl}")
    ax.set_xlabel("Proportion")
    ax.set_ylabel("")
    ax.xaxis.set_major_formatter(lambda x, pos: f"{x:.0%}")

plt.suptitle(f"TE Composition in {snakemake.wildcards.tissue}", fontsize=16, y=1.02)
plt.tight_layout()
plt.savefig(snakemake.output.png, dpi=150)
