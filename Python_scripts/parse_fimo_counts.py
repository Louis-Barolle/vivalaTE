#!/usr/bin/env python3
import sys
from collections import Counter

# usage: parse_fimo_counts.py fimo.tsv bed.bed out.tsv
fimo_tsv, bed_bed, out_tsv = sys.argv[1], sys.argv[2], sys.argv[3]

# 1) count up all FIMO hits per sequence name
counts = Counter()
with open(fimo_tsv) as fh:
    # first line of fimo.tsv is header
    header = fh.readline().rstrip().split("\t")
    seq_idx = header.index("sequence_name")
    for line in fh:
        parts = line.rstrip().split("\t")
        counts[parts[seq_idx]] += 1

# 2) read the BED to get your full list of TEs (col 4)
names = []
with open(bed_bed) as fh:
    for line in fh:
        if not line.startswith("#") and line.strip():
            name = line.split("\t")[3]
            names.append(name)

# 3) write out one row per TE (zero if no hits)
with open(out_tsv, "w") as out:
    out.write("TE\tTFBS_count\n")
    for n in names:
        out.write(f"{n}\t{counts.get(n,0)}\n")
