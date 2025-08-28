#!/usr/bin/env python3
import os, re, gzip
import numpy as np
import pandas as pd
import seaborn as sns, matplotlib.pyplot as plt
import pyBigWig
from matplotlib.ticker import FuncFormatter

# parse genome & tissue from your RNA bigWig paths
RX = re.compile(r".*/(?P<genome>[^/]+)/genomic_data_(?P<tissue>[^/]+)/RNA-seq/.*\.bw$")

def parse_sample(path):
    m = RX.match(path)
    if not m:
        raise ValueError(f"Unexpected RNA bigWig layout: {path}")
    return m.group("genome"), m.group("tissue")

def read_gff_intervals(gff_path, require_feature=None):
    """Return list[(chrom, start0, end)] from (gz)GFF; convert 1-based inclusive → 0-based half-open."""
    opn = gzip.open if gff_path.endswith(".gz") else open
    out = []
    with opn(gff_path, "rt") as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            chrom, feat, start, end = f[0], f[2], int(f[3]), int(f[4])
            if require_feature and feat != require_feature:
                continue
            out.append((chrom, start - 1, end))  # 0-based half-open
    return out

def merge_by_chrom(intervals):
    """Merge overlaps per-chrom. intervals: list[(chrom, s, e)] → dict chrom -> list[(s,e)]."""
    byc = {}
    for c, s, e in intervals:
        byc.setdefault(c, []).append((s, e))
    merged = {}
    for c, lst in byc.items():
        lst.sort()
        cur = []
        for s, e in lst:
            if not cur or s > cur[-1][1]:
                cur.append([s, e])
            else:
                cur[-1][1] = max(cur[-1][1], e)
        merged[c] = [(s, e) for s, e in cur]
    return merged

def sum_bw_over_intervals(bw, merged):
    total = 0.0
    bw_chroms = bw.chroms()
    for chrom, spans in merged.items():
        if chrom not in bw_chroms:
            continue
        for s, e in spans:
            v = bw.stats(chrom, s, e, type="sum")[0]
            if v is not None:
                total += float(v)
    return total

# map GFFs to genomes (by filename prefix)
te_gff_by_genome   = {os.path.basename(p).split(".")[0]: p for p in snakemake.input.te_gff}
gene_gff_by_genome = {os.path.basename(p).split(".")[0]: p for p in snakemake.input.gene_gff}

# cache merged intervals per genome
_te_cache, _gene_cache = {}, {}

def get_te_intervals(genome):
    if genome not in _te_cache:
        te_gff = te_gff_by_genome.get(genome)
        if not te_gff:
            raise ValueError(f"No TE GFF for genome {genome}")
        _te_cache[genome] = merge_by_chrom(read_gff_intervals(te_gff, require_feature=None))
    return _te_cache[genome]

def get_gene_intervals(genome):
    if genome not in _gene_cache:
        ggff = gene_gff_by_genome.get(genome)
        if not ggff:
            raise ValueError(f"No gene GFF for genome {genome}")
        # keep it simple: use exons as the genic territory
        exons = read_gff_intervals(ggff, require_feature="exon")
        if not exons:
            raise ValueError(f"No 'exon' features in {ggff}")
        _gene_cache[genome] = merge_by_chrom(exons)
    return _gene_cache[genome]

# --- compute per bigWig, then collapse replicates to (genome,tissue) ---
rows = []
for bw_path in snakemake.input.rna_bw:
    g, t = parse_sample(bw_path)
    with pyBigWig.open(bw_path) as bw:
        te_sum   = sum_bw_over_intervals(bw,   get_te_intervals(g))
        gene_sum = sum_bw_over_intervals(bw, get_gene_intervals(g))
    rows.append({"genome": g, "tissue": t, "te_sum": te_sum, "gene_sum": gene_sum})

df_raw = pd.DataFrame(rows)

# average replicates per genome×tissue
grp = (
    df_raw.groupby(["genome", "tissue"], as_index=False)[["te_sum", "gene_sum"]]
          .mean()
)

# build plotting dataframe with a 'sample' label and %TE
df = grp.copy()
df["sample"] = df["genome"].astype(str) + "_" + df["tissue"].astype(str)
df["total_signal"] = df["te_sum"] + df["gene_sum"]
# avoid divide-by-zero → NaN; fill with 0 for plotting
df["pct"] = (df["te_sum"] / df["total_signal"]) * 100.0
df["pct"] = df["pct"].replace([np.inf, -np.inf], np.nan).fillna(0.0)
df["library"] = "RNA-seq"  # informative label for the table
df = df.sort_values(["genome", "tissue", "sample"])

def nice_cap(x):
    for cap in [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100]:
        if x <= cap:
            return cap
    return 100

ymax = nice_cap(float(np.nanmax(df["pct"])) if len(df) else 1.0)

# compact summary table above the plot
tab = {
    "Sample":        df["sample"].astype(str).tolist(),
    "Library":       df["library"].astype(str).tolist(),
    "Total signal":  df["total_signal"].round(0).astype("Int64").astype(str).tolist(),
    "TE signal":     df["te_sum"].round(0).astype("Int64").astype(str).tolist(),
    "TE %":         [f"{v:.2f}%" if np.isfinite(v) else "—" for v in df["pct"]],
}

# --- figure: table + barplot ---
sns.set_context("talk")
fig = plt.figure(figsize=(max(12, 0.35 * len(df)), 6))
gs  = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[1, 3])

# table panel
ax0 = fig.add_subplot(gs[0]); ax0.axis("off")
tbl = ax0.table(
    cellText=list(map(list, zip(*tab.values()))),  # rows
    colLabels=list(tab.keys()),
    loc="center", cellLoc="center"
)
tbl.auto_set_font_size(False); tbl.set_fontsize(8); tbl.scale(1, 1.15)

# barplot panel
ax = fig.add_subplot(gs[1])
sns.barplot(data=df, x="sample", y="pct", color="#4055A8", ax=ax)
ax.set_ylim(0, ymax)
ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:g}%"))
ax.set_xlabel("")
ax.set_ylabel("% of total RNA assigned to TEs")
ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha="center", fontsize=8)
ax.grid(axis="y", color="0.85", linewidth=0.8)

# on-bar labels
for p in ax.patches:
    h = p.get_height()
    if np.isfinite(h) and h > 0:
        ax.annotate(f"{h:.2f}%", (p.get_x() + p.get_width()/2, h),
                    ha="center", va="bottom", fontsize=7,
                    xytext=(0, 2), textcoords="offset points")

fig.suptitle("TE fraction of transcriptome per sample (replicates averaged)", y=0.99, fontsize=12)
plt.tight_layout(rect=[0, 0, 1, 0.98])
plt.savefig(snakemake.output[0], dpi=150, bbox_inches="tight")
plt.close(fig)
