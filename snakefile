import os, glob, re
configfile: "config.yaml"

# 1) project-level config
BASE    = config["base_dir"]       
GENOMES = config["genomes"]
TISSUES = config["tissues"]

# where your helper scripts live
SCRIPT_DIR = os.path.join(BASE, "Python_scripts")
def S(script_name):
    return os.path.join(SCRIPT_DIR, script_name)

# helper to fill in path templates from config.yaml
def tpl(key, genome=None, tissue=None):
    return config[f"{key}_template"].format(
        base_dir=BASE,
        genome=genome,
        tissue=tissue
    )

# 2) per-genome × per-tissue input dirs
CHIP_DIRS = { (g,t): tpl("chip_dir", genome=g, tissue=t)  for g in GENOMES for t in TISSUES }
RNA_DIRS  = { (g,t): tpl("rna_dir",  genome=g, tissue=t)  for g in GENOMES for t in TISSUES }

# 3) discover samples
def list_bw(d):
    return [os.path.splitext(os.path.basename(p))[0]
            for p in glob.glob(os.path.join(d, "*.bw"))]

CHIP_SAMPLES = { (g,t): list_bw(CHIP_DIRS[(g,t)]) for g in GENOMES for t in TISSUES }
RNA_SAMPLES  = { (g,t): list_bw(RNA_DIRS[(g,t)])  for g in GENOMES for t in TISSUES }

# split ChIP by mark
ACT_SAMPLES = { (g,t): [s for s in CHIP_SAMPLES[(g,t)] if "H3K27ac" in s]
                for g in GENOMES for t in TISSUES }
REP_SAMPLES = { (g,t): [s for s in CHIP_SAMPLES[(g,t)] if "H3K9me3" in s]
                for g in GENOMES for t in TISSUES }

# 4) static resources
TE_CLASSES = config["te_classes_tsv"].format(base_dir=BASE)

rule all:
    input:
        # per-genome, per-tissue outputs
        expand("{genome}/results/{tissue}/final_matrix.tsv",                genome=GENOMES, tissue=TISSUES),
        expand("{genome}/results/{tissue}/family_heatmaps",                 genome=GENOMES, tissue=TISSUES),
        expand("{genome}/results/{tissue}/heatmap_RNAseq_filtered.png",     genome=GENOMES, tissue=TISSUES),
        expand("{genome}/results/{tissue}/TE_list_rna_ge1.txt",             genome=GENOMES, tissue=TISSUES),
        expand("{genome}/results/{tissue}/expression_by_mark.png",          genome=GENOMES, tissue=TISSUES),
        expand("{genome}/results/{tissue}/expression_by_mark_expressed.png",genome=GENOMES, tissue=TISSUES),
        expand("{genome}/results/{tissue}/expression_by_mark.html",         genome=GENOMES, tissue=TISSUES),
        expand("{genome}/results/{tissue}/expression_vs_mark.png",          genome=GENOMES, tissue=TISSUES),
        expand("{genome}/results/{tissue}/TE_piecharts.png",                genome=GENOMES, tissue=TISSUES),
        expand("{genome}/results/{tissue}/scatter_mark_vs_expr_H3K27ac.png",genome=GENOMES, tissue=TISSUES),
        expand("{genome}/results/{tissue}/scatter_mark_vs_expr_H3K9me3.png",genome=GENOMES, tissue=TISSUES),
        expand("{genome}/results/{tissue}/boxplot_quartiles.png",           genome=GENOMES, tissue=TISSUES),
        expand("{genome}/results/{tissue}/corr_panels.png",                 genome=GENOMES, tissue=TISSUES),
        expand("{genome}/results/{tissue}/metaplot_signal.png",             genome=GENOMES, tissue=TISSUES),
        expand("{genome}/results/{tissue}/scatter_rna_chip.png",            genome=GENOMES, tissue=TISSUES),
        expand("{genome}/results/{tissue}/feature_importance.png",          genome=GENOMES, tissue=TISSUES),
        expand("{genome}/results/{tissue}/TE_composition.png",              genome=GENOMES, tissue=TISSUES),

        # per-genome static files
        expand("{genome}/bed/{genome}_genes.bed",      genome=GENOMES),
        expand("{genome}/results/TE_context.tsv",      genome=GENOMES),
        expand("{genome}/results/{tissue}/TE_summary.tsv", genome=GENOMES, tissue=TISSUES),

        # compare within each genome
        expand("{genome}/results/compare/gut_head_overlap.txt",             genome=GENOMES),
        expand("{genome}/results/compare/gut_head_venn.png",                genome=GENOMES),
        expand("{genome}/results/compare/gut_head_overlap_composition.png", genome=GENOMES)

##### BED preparation #####
rule gff_to_bed:
    input:
        te_gff=lambda wc: tpl("te_gff", genome=wc.genome)
    output:
        bed="{genome}/bed/{genome}_TEs.bed"
    run:
        import gzip, pandas as pd
        rec=[]
        with gzip.open(input.te_gff,"rt") as f:
            for L in f:
                if L.startswith("#"): continue
                chrom,_,_,s,e,_,strand,_,attrs = L.rstrip().split("\t")
                m = re.search(r"Family=([^;]+)", attrs)
                if not m: continue
                fam, start, end = m.group(1), int(s)-1, int(e)
                name = f"{fam}_{chrom}_{start}-{end}"
                rec.append((chrom, start, end, name, strand))
        pd.DataFrame(rec, columns=["chr","start","end","name","strand"]) \
          .to_csv(output.bed, sep="\t", header=False, index=False)

rule gff3_to_gene_bed:
    input:
        gff3=lambda wc: tpl("gene_gff", genome=wc.genome)
    output:
        gene_bed="{genome}/bed/{genome}_genes.bed"
    run:
        import gzip, pandas as pd
        recs=[]
        with gzip.open(input.gff3, "rt") as f:
            for L in f:
                if L.startswith("#"): continue
                chrom, _, tp, s, e, _, _, _, attrs = L.rstrip().split("\t")
                if tp != "gene": continue
                start, end = int(s)-1, int(e)
                gid = next((kv.split("=",1)[1] for kv in attrs.split(";") if kv.startswith("ID=")),
                           f"{chrom}:{start}-{end}")
                recs.append((chrom, start, end, gid))
        pd.DataFrame(recs, columns=["chr","start","end","gene_id"]) \
          .to_csv(output.gene_bed, sep="\t", header=False, index=False)

rule extract_upstream:
    input:
        bed="{genome}/bed/{genome}_TEs.bed"
    output:
        up1k="{genome}/bed/{genome}_TEs.upstream1kb.bed"
    run:
        import pandas as pd
        df=pd.read_csv(input.bed, sep="\t", header=None,
                       names=["chr","start","end","name","strand"])
        df2=df.assign(
            start=df.start.sub(5000).clip(lower=0).where(df.strand=="+", df.end),
            end  =df.start.where(df.strand=="+", df.end.add(5000))
        )[['chr','start','end','name']]
        df2.to_csv(output.up1k, sep="\t", header=False, index=False)

rule sort_bed:
    input:
        up1k="{genome}/bed/{genome}_TEs.upstream1kb.bed"
    output:
        up1k_sorted="{genome}/bed/{genome}_TEs.upstream1kb.sorted.bed"
    shell:
        "bedtools sort -i {input.up1k} > {output.up1k_sorted}"

rule sort_te_bed:
    input:
        te="{genome}/bed/{genome}_TEs.bed"
    output:
        te_sorted="{genome}/bed/{genome}_TEs.sorted.bed"
    shell:
        "bedtools sort -i {input.te} > {output.te_sorted}"

rule make_TE_start_windows:
    input:
        bed="{genome}/bed/{genome}_TEs.bed"
    output:
        start125="{genome}/bed/{genome}_TEs.start125bp.bed"
    run:
        import pandas as pd
        df=pd.read_csv(input.bed, sep="\t", header=None,
                       names=["chr","start","end","TE","strand"])
        rec=[]; half=125
        for _,r in df.iterrows():
            if r.strand=="+":
                s=max(0,r.start-half); e=r.start+half
            else:
                s=max(0,r.end-half);   e=r.end+half
            rec.append((r.chr,s,e,r.TE))
        pd.DataFrame(rec, columns=["chr","start","end","TE"]) \
          .to_csv(output.start125, sep="\t", header=False, index=False)

rule sort_start_windows:
    input:
        start125="{genome}/bed/{genome}_TEs.start125bp.bed"
    output:
        start125_sorted="{genome}/bed/{genome}_TEs.start125bp.sorted.bed"
    shell:
        "bedtools sort -i {input.start125} > {output.start125_sorted}"


##### Signal mapping #####
rule map_chip:
    input:
        bed="{genome}/bed/{genome}_TEs.upstream1kb.sorted.bed",
        bw=lambda wc: os.path.join(CHIP_DIRS[(wc.genome,wc.tissue)], f"{wc.sample}.bw")
    output:
        "{genome}/results/{tissue}/{sample}_chip_signal.tsv"
    run:
        import pyBigWig, pandas as pd
        bwf=pyBigWig.open(input.bw)
        df=pd.read_csv(input.bed, sep="\t", header=None, names=["chr","start","end","TE"])
        vals=[(r.TE, bwf.stats(r.chr,r.start,r.end,type="mean")[0] or 0)
              for _,r in df.iterrows()]
        bwf.close()
        pd.DataFrame(vals, columns=["TE", wildcards.sample]) \
          .to_csv(output[0], sep="\t", index=False)

rule map_rna:
    input:
        bed="{genome}/bed/{genome}_TEs.sorted.bed",
        bw=lambda wc: os.path.join(RNA_DIRS[(wc.genome,wc.tissue)], f"{wc.sample}.bw")
    output:
        "{genome}/results/{tissue}/{sample}_RNA_signal.tsv"
    run:
        import pyBigWig, pandas as pd
        bw=pyBigWig.open(input.bw)
        df=pd.read_csv(input.bed, sep="\t", header=None,
                       usecols=[0,1,2,3], names=["chr","start","end","TE"])
        vals=[(r.TE, bw.stats(r.chr,r.start,r.end,type="mean")[0] or 0)
              for _,r in df.iterrows()]
        bw.close()
        pd.DataFrame(vals, columns=["TE", wildcards.sample]) \
          .to_csv(output[0], sep="\t", index=False)


##### Combine & summarise #####
rule combine_all:
    input:
        chip=lambda wc: [f"{wc.genome}/results/{wc.tissue}/{s}_chip_signal.tsv"
                         for s in CHIP_SAMPLES[(wc.genome,wc.tissue)]],
        rna =lambda wc: [f"{wc.genome}/results/{wc.tissue}/{s}_RNA_signal.tsv"
                         for s in RNA_SAMPLES[(wc.genome,wc.tissue)]]
    output:
        "{genome}/results/{tissue}/final_matrix.tsv"
    run:
        import os, pandas as pd
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        dfs_c=[pd.read_csv(x, sep="\t", index_col=0) for x in input.chip]
        dfs_r=[pd.read_csv(x, sep="\t", index_col=0) for x in input.rna]
        df_c = pd.concat(dfs_c,axis=1).groupby(lambda c: c.rsplit("_",1)[-1],axis=1).mean()
        df_r=pd.concat(dfs_r, axis=1).mean(axis=1).to_frame("RNA")
        pd.concat([df_c, df_r], axis=1).to_csv(output[0], sep="\t")

rule classify_context:
    input:
        bed   = "{genome}/bed/{genome}_TEs.sorted.bed",
        genes = "{genome}/bed/{genome}_genes.bed"
    output:
        "{genome}/results/TE_context.tsv"
    run:
        import os, pandas as pd, subprocess
        res_dir = os.path.dirname(output[0])
        os.makedirs(res_dir, exist_ok=True)
        tmp = f"{wildcards.genome}/bed/tmp_TE.bed"
        subprocess.run(f"cp {input.bed} {tmp}", shell=True, check=True)
        genic_file = f"{wildcards.genome}/results/_genic.txt"
        cmd = (
            f"bedtools intersect -wa -a {tmp} -b {input.genes} "
            f"| cut -f4 | sort | uniq > {genic_file}"
        )
        subprocess.run(cmd, shell=True, check=True)
        all_te = pd.read_csv(input.bed, sep="\t", header=None, usecols=[3], names=["TE"])
        genic  = pd.read_csv(genic_file, header=None, names=["TE"])
        all_te["context"] = all_te.TE.isin(genic.TE).map({True: "genic", False: "intergenic"})
        all_te.to_csv(output[0], sep="\t", index=False)
        os.remove(tmp)
        os.remove(genic_file)


rule summarize_TE_signals:
    input:
        chip=lambda wc: [f"{wc.genome}/results/{wc.tissue}/{s}_chip_signal.tsv"
                         for s in CHIP_SAMPLES[(wc.genome,wc.tissue)]],
        rna =lambda wc: [f"{wc.genome}/results/{wc.tissue}/{s}_RNA_signal.tsv"
                         for s in RNA_SAMPLES[(wc.genome,wc.tissue)]],
        cont="{genome}/results/TE_context.tsv"
    output:
        "{genome}/results/{tissue}/TE_summary.tsv"
    run:
        import os, pandas as pd
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        df_ctx=pd.read_csv(input.cont, sep="\t", index_col=0)
        df_act=pd.concat([pd.read_csv(x, sep="\t", index_col=0)
                          for x in input.chip if "H3K27ac" in x], axis=1).mean(axis=1)
        df_rep=pd.concat([pd.read_csv(x, sep="\t", index_col=0)
                          for x in input.chip if "H3K9me3" in x], axis=1).mean(axis=1)
        df_rna=pd.concat([pd.read_csv(x, sep="\t", index_col=0)
                          for x in input.rna], axis=1).mean(axis=1).to_frame("RNA")
        summary=pd.DataFrame({"context":df_ctx.context,
                              "H3K27ac_flank":df_act,
                              "H3K9me3_flank":df_rep}).join(df_rna).fillna(0)
        summary.to_csv(output[0], sep="\t")


##### Visualization #####
rule heatmap:
    input:
        "{genome}/results/{tissue}/final_matrix.tsv"
    output:
        "{genome}/results/{tissue}/heatmap.png"
    script:
        S("plot_heatmap2.py")

rule heatmap_by_family:
    input:
        "{genome}/results/{tissue}/final_matrix.tsv"
    output:
        directory("{genome}/results/{tissue}/family_heatmaps")
    script:
        S("plot_heatmap_by_family.py")

rule heatmap_rna_filtered:
    input:
        "{genome}/results/{tissue}/final_matrix.tsv"
    output:
        heatmap="{genome}/results/{tissue}/heatmap_RNAseq_filtered.png",
        te_list="{genome}/results/{tissue}/TE_list_rna_ge1.txt"
    script:
        S("plot_heatmap_RNAseq_filtered.py")

rule plot_expression_by_mark:
    input:
        "{genome}/results/{tissue}/TE_summary.tsv"
    output:
        "{genome}/results/{tissue}/expression_by_mark.png"
    script:
        S("plot_expression_by_mark.py")

rule plot_expression_by_mark_expressed:
    input:
        "{genome}/results/{tissue}/TE_summary.tsv"
    output:
        "{genome}/results/{tissue}/expression_by_mark_expressed.png"
    script:
        S("plot_expression_by_mark_expressed.py")

rule plot_expression_by_mark_interactive:
    input:
        "{genome}/results/{tissue}/TE_summary.tsv"
    output:
        "{genome}/results/{tissue}/expression_by_mark.html"
    script:
        S("plot_expression_by_mark_interactive.py")

rule plot_expression_vs_mark:
    input:
        summary="{genome}/results/{tissue}/TE_summary.tsv"
    output:
        png="{genome}/results/{tissue}/expression_vs_mark.png"
    script:
        S("plot_expression_vs_mark.py")

rule plot_TE_pies:
    input:
        "{genome}/results/{tissue}/TE_summary.tsv"
    output:
        "{genome}/results/{tissue}/TE_piecharts.png"
    script:
        S("plot_TE_piecharts.py")

rule plot_corr_panels:
    input:
        matrix="{genome}/results/{tissue}/final_matrix.tsv"
    output:
        png="{genome}/results/{tissue}/corr_panels.png"
    script:
        S("plot_corr_expressed.py")

rule scatter_mark_vs_expression:
    input:
        matrix="{genome}/results/{tissue}/final_matrix.tsv"
    output:
        "{genome}/results/{tissue}/scatter_mark_vs_expr_H3K27ac.png",
        "{genome}/results/{tissue}/scatter_mark_vs_expr_H3K9me3.png"
    script:
        S("plot_scatter_mark_vs_expression.py")

rule boxplot_quartiles:
    input:
        matrix="{genome}/results/{tissue}/final_matrix.tsv"
    output:
        "{genome}/results/{tissue}/boxplot_quartiles.png"
    script:
        S("plot_boxplot_mark_quartile.py")

rule metaplot:
    input:
        matrix="{genome}/results/{tissue}/final_matrix.tsv"
    output:
        png="{genome}/results/{tissue}/metaplot_signal.png"
    script:
        S("plot_metaplot.py")

rule scatter_rna_chip:
    input:
        matrix="{genome}/results/{tissue}/final_matrix.tsv"
    output:
        png="{genome}/results/{tissue}/scatter_rna_chip.png"
    script:
        S("plot_rna_chip_scatter.py")

rule feature_importance:
    input:
        matrix="{genome}/results/{tissue}/final_matrix.tsv"
    output:
        "{genome}/results/{tissue}/feature_importance.png"
    script:
        S("plot_feature_importance.py")

rule plot_te_composition:
    input:
        matrix="{genome}/results/{tissue}/final_matrix.tsv",
        classes=TE_CLASSES
    output:
        png="{genome}/results/{tissue}/TE_composition.png"
    script:
        S("plot_TE_composition.py")

rule compare_gut_head:
    input:
        gut = "{genome}/results/gut/final_matrix.tsv",
        head= "{genome}/results/head/final_matrix.tsv",
        cls = TE_CLASSES
    output:
        overlap     = "{genome}/results/compare/gut_head_overlap.txt",
        venn        = "{genome}/results/compare/gut_head_venn.png",
        composition = "{genome}/results/compare/gut_head_overlap_composition.png"
    script:
        S("compare_TE_overlap.py")
