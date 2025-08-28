configfile: "config.yaml"
import os
import time, datetime as dt

# --- tiny pretty banners ---
START = time.time()
def banner(title, sub=None, w=70):
    w = max(w, len(title)+4, len(sub)+4 if sub else 0)
    print(f"\n╔{'═'*w}╗")
    print(f"║ {title.center(w-2)} ║")
    if sub: print(f"║ {sub.center(w-2)} ║")
    print(f"╚{'═'*w}╝\n")

onstart:
    banner("Louis Barolle's epic code now starting B) : Comparing different genomes",
           f"Started: {dt.datetime.now():%Y-%m-%d %H:%M:%S}")

onsuccess:
    banner("Workflow finished :D !!!", f"Elapsed: {(time.time()-START)/60:.1f} min")

onerror:
    banner("Despite my immense intellect, the workflow failed :(",
           "Check .snakemake/log/* and rule logs")
# --- end banners ---




# config
GENOMES = config["genomes"]
TISSUES = config["tissues"]
BASE = config["base_dir"]

def tpl(key, genome=None, tissue=None):
    return config[f"{key}_template"].format(
        base_dir=BASE, genome=genome, tissue=tissue
    )

RNA_DIRS = { (g,t): tpl("rna_dir", genome=g, tissue=t)
             for g in GENOMES for t in TISSUES }

from glob import glob
RNA_REPS = { (g,t): sorted(
                glob(os.path.join(RNA_DIRS[(g,t)], "*.bw")) +
                glob(os.path.join(RNA_DIRS[(g,t)], "*.bigWig"))
            )
            for g in GENOMES for t in TISSUES }


RNA_BW_ALL = [bw for v in RNA_REPS.values() for bw in v]


# only inter-genome, ordered pairs (no self/self, no duplicates)
PAIRS = [(g1, g2) for g1 in GENOMES for g2 in GENOMES if g1 < g2]

# per-pair outputs live under compare/{g1}_vs_{g2}/
RBH       = [f"compare/{g1}_vs_{g2}/rbh.tsv"                           for g1, g2 in PAIRS]
UNIQ      = [f"compare/{g1}_vs_{g2}/{g1}_unique_to_{g2}.tsv"           for g1, g2 in PAIRS] + \
            [f"compare/{g1}_vs_{g2}/{g2}_unique_to_{g1}.tsv"           for g1, g2 in PAIRS]
SHARED_HP = [f"compare/{g1}_vs_{g2}/shared_heatmap.png"                for g1, g2 in PAIRS]
UNIQ_HP   = [f"compare/{g1}_vs_{g2}/{g1}_unique_to_{g2}.png"           for g1, g2 in PAIRS] + \
            [f"compare/{g1}_vs_{g2}/{g2}_unique_to_{g1}.png"           for g1, g2 in PAIRS]

# per-genome “all-uniques” heatmap 
UNIQ_GENOME_HP = [f"compare/{g}_unique_heatmap.png" for g in GENOMES]

rule all:
    input:
        "compare/subfamily_matrix.tsv",
        "compare/subfamily_heatmap_all.png",
        "compare/subfamily_heatmap_active.png",
        *RBH,
        *UNIQ,
        "compare/TE_matrix.tsv",
        "compare/variable_TE_heatmap.png",
        *SHARED_HP,
        *UNIQ_HP,
        *UNIQ_GENOME_HP,
        "compare/TE_transcript_percentage.png"

rule aggregate_subfamily_matrix:
    input:
        expand("{genome}/results/{tissue}/TE_summary.tsv",
               genome=GENOMES, tissue=TISSUES)
    output: "compare/subfamily_matrix.tsv"
    run:
        import pandas as pd, os
        dfs = []
        for fn in input:
            g,_,t,_ = fn.split("/")[:4]
            sig = pd.read_csv(fn, sep="\t", index_col=0)
            fam = sig.index.to_series().str.split("_").str[0]
            for c in ("H3K27ac_flank","H3K9me3_flank","RNA"):
                if c in sig:
                    s = sig[c].groupby(fam).mean()
                    s.name = f"{g}_{t}_{c}"
                    dfs.append(s)
        mat = pd.concat(dfs, axis=1).fillna(0)
        os.makedirs("compare", exist_ok=True)
        mat.to_csv(output[0], sep="\t")

rule plot_subfamily_heatmap:
    input: "compare/subfamily_matrix.tsv"
    output:
        "compare/subfamily_heatmap_all.png",
        "compare/subfamily_heatmap_active.png"
    script: "Python_scripts/plot_subfamily_heatmap.py"

rule extract_TE_flank:
    input:
        bed = "{genome}/bed/{genome}_TEs.sorted.bed",
        fa  = "{genome}/Genome/{genome}.fasta",
        fai = "{genome}/Genome/{genome}.fasta.fai"
    output: "compare/{genome}_TEs_1000bpFlank.fa"
    params: flank = 1000
    run:
        import subprocess
        g = wildcards.genome
        cs = f"{g}/bed/{g}.chrom.sizes"
        os.makedirs(os.path.dirname(cs), exist_ok=True)
        # write chrom.sizes from .fai
        with open(input.fai) as fin, open(cs, "w") as fout:
            for line in fin:
                chrom, length = line.split()[:2]
                fout.write(f"{chrom}\t{length}\n")
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        cmd = (
            f"bedtools slop -i {input.bed} -g {cs} -b {params.flank} "
            f"| bedtools getfasta -fi {input.fa} -bed - -fo {output[0]} -name"
        )
        subprocess.run(cmd, shell=True, check=True)
        os.remove(cs)


# Build BLAST DB once per genome (avoids concurrent rebuilds)
rule format_blastdb:
    input:
        fa = "compare/{genome}_TEs_1000bpFlank.fa"
    output:
        touch("compare/db/{genome}.dbdone")
    params:
        prefix = "compare/db/{genome}"
    shell:
        r"""
        mkdir -p compare/db
        makeblastdb -in {input.fa} -dbtype nucl -out {params.prefix}
        touch {output}
        """


rule reciprocal_best_hits:
    input:
        fa1 = "compare/{genome1}_TEs_1000bpFlank.fa",
        fa2 = "compare/{genome2}_TEs_1000bpFlank.fa",
        db1 = "compare/db/{genome1}.dbdone",
        db2 = "compare/db/{genome2}.dbdone"
    output:
        rbh   = "compare/{genome1}_vs_{genome2}/rbh.tsv",
        uniq1 = "compare/{genome1}_vs_{genome2}/{genome1}_unique_to_{genome2}.tsv",
        uniq2 = "compare/{genome1}_vs_{genome2}/{genome2}_unique_to_{genome1}.tsv"
    shell:
        r"""
        pairdir="compare/{wildcards.genome1}_vs_{wildcards.genome2}"
        mkdir -p "$pairdir"

        # BLAST using prebuilt DB prefixes
        blastn -query {input.fa1} -db compare/db/{wildcards.genome2} \
               -qcov_hsp_perc 90 -outfmt 6 \
        | sort -k1,1 -k12,12nr | sort -u -k1,1 \
        > "$pairdir/{wildcards.genome1}_to_{wildcards.genome2}.tsv"

        blastn -query {input.fa2} -db compare/db/{wildcards.genome1} \
               -qcov_hsp_perc 90 -outfmt 6 \
        | sort -k1,1 -k12,12nr | sort -u -k1,1 \
        > "$pairdir/{wildcards.genome2}_to_{wildcards.genome1}.tsv"

        # RBH join (A.query col1 ↔ B.subject col2)
        join -t $'\t' -1 1 -2 2 \
             <(sort -k1,1 "$pairdir/{wildcards.genome1}_to_{wildcards.genome2}.tsv") \
             <(sort -k2,2 "$pairdir/{wildcards.genome2}_to_{wildcards.genome1}.tsv") \
        > {output.rbh}

        # uniques = all headers in FASTA minus RBH IDs
        comm -23 \
          <(grep '^>' {input.fa1} | sed 's/^>//' | sort -u) \
          <(cut -f1 {output.rbh} | sort -u) \
          > {output.uniq1}

        comm -23 \
          <(grep '^>' {input.fa2} | sed 's/^>//' | sort -u) \
          <(cut -f2 {output.rbh} | sort -u) \
          > {output.uniq2}

        rm -f "$pairdir/{wildcards.genome1}_to_{wildcards.genome2}.tsv" \
              "$pairdir/{wildcards.genome2}_to_{wildcards.genome1}.tsv"
        """


rule aggregate_TE_matrix:
    input:
        expand("{genome}/results/{tissue}/TE_summary.tsv",
               genome=GENOMES, tissue=TISSUES)
    output: "compare/TE_matrix.tsv"
    run:
        import pandas as pd, os
        dfs = []
        for fn in input:
            parts = fn.split("/")
            g, t = parts[0], parts[2]
            sig = pd.read_csv(fn, sep="\t", index_col=0)
            for col in ("H3K27ac_flank","H3K9me3_flank","RNA"):
                if col in sig.columns:
                    s = sig[col]
                    s.name = f"{g}_{t}_{col}"
                    dfs.append(s)
        mat = pd.concat(dfs, axis=1)
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        mat.to_csv(output[0], sep="\t")

rule plot_top_variable_TEs:
    input: mat="compare/TE_matrix.tsv"
    output: "compare/variable_TE_heatmap.png"
    script: "Python_scripts/plot_top_variable_TEs.py"

rule plot_single_TE:
    input:  mat="compare/TE_matrix.tsv"
    params: TE_ID=lambda wc: wc.TE_ID
    output: "compare/single_TE_{TE_ID}.png"
    script: "Python_scripts/plot_single_TE.py"

# per-pair heatmaps go inside the pair dir now
rule plot_shared_vs_unique_TEs:
    input:
        mat   = "compare/TE_matrix.tsv",
        rbh   = "compare/{genome1}_vs_{genome2}/rbh.tsv",
        uniq1 = "compare/{genome1}_vs_{genome2}/{genome1}_unique_to_{genome2}.tsv",
        uniq2 = "compare/{genome1}_vs_{genome2}/{genome2}_unique_to_{genome1}.tsv"
    output:
        shared  = "compare/{genome1}_vs_{genome2}/shared_heatmap.png",
        unique1 = "compare/{genome1}_vs_{genome2}/{genome1}_unique_to_{genome2}.png",
        unique2 = "compare/{genome1}_vs_{genome2}/{genome2}_unique_to_{genome1}.png"
    script: "Python_scripts/plot_TE_sets_heatmap.py"

# per-genome “all uniques across others” heatmap (still in compare/)
rule plot_unique_genome_heatmap:
    input:
        mat="compare/TE_matrix.tsv",
        uniq=lambda wc: [
            f"compare/{min(wc.genome,o)}_vs_{max(wc.genome,o)}/{wc.genome}_unique_to_{o}.tsv"
            for o in GENOMES if o != wc.genome
        ]
    output:
        "compare/{genome}_unique_heatmap.png"
    script:
        "Python_scripts/plot_unique_genome_heatmap.py"
        
        
rule plot_TE_transcript_percentage:
    input:
        rna_bw = RNA_BW_ALL,
        te_gff = expand(config["te_gff_template"],
                        base_dir=BASE, genome=GENOMES),
        gene_gff = expand(config["gene_gff_template"],
                          base_dir=BASE, genome=GENOMES)
    output:
        "compare/TE_transcript_percentage.png"
    script:
        "Python_scripts/plot_TE_transcript_percentage.py"
