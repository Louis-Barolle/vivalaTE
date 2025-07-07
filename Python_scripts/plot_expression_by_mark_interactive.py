#!/usr/bin/env python3
import pandas as pd
import numpy as np
import plotly.express as px

# 1) load
df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
df.index.name = "TE"  # ensure the index is named so we can hover it

# 2) pick flank marks
h3a = next(c for c in df if c.startswith("H3K27ac") and c.endswith("_flank"))
h3r = next(c for c in df if c.startswith("H3K9me3")  and c.endswith("_flank"))

# 3) classify
df["context"] = np.where(df[h3a] > df[h3r], "active", "inactive")

# 4) log‐transform RNA
df["RNA_log2"] = np.log2(df["RNA"].fillna(0) + 1)

# 5) build a combined DataFrame for plotting
plot_df = df.reset_index()  # brings TE names into a column

# 6) create the interactive violin + strip/jitter
fig = px.violin(
    plot_df, x="context", y="RNA_log2", color="context",
    category_orders={"context":["active","inactive"]},
    color_discrete_map={"active":"#fc8d62","inactive":"#66c2a5"},
    hover_data=["TE","RNA"], box=True, points="all",  # box=True draws the quartiles, points="all" adds jittered dots
    labels={"RNA_log2":"log₂(RNA+1)","context":"Chromatin context"},
    title="Expression by chromatin context"
)

# 7) annotate counts in the subtitle
counts = plot_df["context"].value_counts()
fig.update_layout(
    title={
      "text": f"Expression by chromatin context<br>"
              f"n(active)={counts['active']}, n(inactive)={counts['inactive']}",
      "x":0.5
    }
)

# 8) save to standalone HTML
fig.write_html(snakemake.output[0])
