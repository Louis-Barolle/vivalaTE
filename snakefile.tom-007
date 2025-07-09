import os, glob, re
configfile: "config.yaml"

# 1) tissues to process in parallel
TISSUES    = config["tissues"]
TE_GFF     = config["te_gff"]
GENE_GFF   = config["gene_gff"]
GENOME    = config["genome_fasta"]
MOTIFS    = config["motif_meme"]

# 2) per‐tissue directories
CHIP_DIRS  = {t: config[f"chip_dir_{t}"] for t in TISSUES}
RNA_DIRS   = {t: config[f"rna_dir_{t}"]  for t in TISSUES}

# 3) sample lists per tissue
CHIP_SAMPLES = {
    t: [os.path.splitext(os.path.basename(p))[0]
        for p in glob.glob(os.path.join(CHIP_DIRS[t],"*.bw"))]
    for t in TISSUES
}
RNA_SAMPLES = {
    t: [os.path.splitext(os.path.basename(p))[0]
        for p in glob.glob(os.path.join(RNA_DIRS[t],"*.bw"))]
    for t in TISSUES
}
 
# 4) split ChIP replicates by mark
ACT_SAMPLES = {t: [s for s in CHIP_SAMPLES[t] if "H3K27ac" in s] for t in TISSUES}
REP_SAMPLES = {t: [s for s in CHIP_SAMPLES[t] if "H3K9me3" in s] for t in TISSUES}

rule all:
    input:
        expand("results/{tissue}/final_matrix.tsv",                tissue=TISSUES),
        expand("results/{tissue}/family_heatmaps",                 tissue=TISSUES),
        expand("results/{tissue}/heatmap_RNAseq_filtered.png",     tissue=TISSUES),
        expand("results/{tissue}/TE_list_rna_ge1.txt",             tissue=TISSUES),
        expand("results/{tissue}/expression_by_mark.png",          tissue=TISSUES),
        expand("results/{tissue}/expression_by_mark_expressed.png",tissue=TISSUES),
        expand("results/{tissue}/expression_by_mark.html",         tissue=TISSUES),
        expand("results/{tissue}/expression_vs_mark.png",          tissue=TISSUES),
        expand("results/{tissue}/TE_piecharts.png",                tissue=TISSUES),
        expand("results/{tissue}/scatter_mark_vs_expr_H3K27ac.png",tissue=TISSUES),
        expand("results/{tissue}/scatter_mark_vs_expr_H3K9me3.png",tissue=TISSUES),
        expand("results/{tissue}/boxplot_quartiles.png",           tissue=TISSUES),
        expand("results/{tissue}/corr_panels.png",                 tissue=TISSUES),
        expand("results/{tissue}/metaplot_signal.png",            tissue=TISSUES),
        expand("results/{tissue}/scatter_rna_chip.png",           tissue=TISSUES),
        expand("results/{tissue}/feature_importance.png",         tissue=TISSUES),
        expand("results/{tissue}/TE_composition.png",             tissue=TISSUES),
        "bed/TOM-007_genes.bed",
        "results/TE_context.tsv",
        expand("results/{tissue}/TE_summary.tsv",                 tissue=TISSUES),
        "results/compare/gut_head_overlap.txt",
        "results/compare/gut_head_overlap_composition.png",
        "results/compare/gut_head_venn.png",
        "results/TE_upstream1kb.fa",
        "results/fimo_upstream1kb/fimo.tsv",
        "results/plots/fimo_upstream1kb_density.png"
        
        

##### BED preparation #####
rule gff_to_bed:
    input:  te_gff=TE_GFF
    output: bed="bed/TOM-007_TEs.bed"
    run:
        import gzip,pandas as pd
        rec=[]
        with gzip.open(input.te_gff,"rt") as f:
            for L in f:
                if L.startswith("#"): continue
                chrom,_,_,s,e,_,strand,_,attrs=L.rstrip().split("\t")
                m=re.search(r"Family=([^;]+)",attrs)
                if not m: continue
                fam, start, end = m.group(1), int(s)-1, int(e)
                name = f"{fam}_{chrom}_{start}-{end}"
                rec.append((chrom,start,end,name,strand))
        pd.DataFrame(rec,columns=["chr","start","end","name","strand"]) \
          .to_csv(output.bed,sep="\t",header=False,index=False)

rule extract_upstream:
    input:  "bed/TOM-007_TEs.bed"
    output: "bed/TOM-007_TEs.upstream1kb.bed"
    run:
        import pandas as pd
        df=pd.read_csv(input[0],sep="\t",header=None,
                       names=["chr","start","end","name","strand"])
        df2 = df.assign(
            start=df.start.sub(5000).clip(lower=0).where(df.strand=="+",df.end),
            end  =df.start.where(df.strand=="+",df.end.add(5000))
        )[['chr','start','end','name']]
        df2.to_csv(output[0],sep="\t",header=False,index=False)

rule sort_bed:
    input:  "bed/TOM-007_TEs.upstream1kb.bed"
    output: "bed/TOM-007_TEs.upstream1kb.sorted.bed"
    shell: "bedtools sort -i {input} > {output}"

rule sort_te_bed:
    input:  "bed/TOM-007_TEs.bed"
    output: "bed/TOM-007_TEs.sorted.bed"
    shell: "bedtools sort -i {input} > {output}"

rule make_TE_start_windows:
    input:  bed="bed/TOM-007_TEs.bed"
    output: bed="bed/TOM-007_TEs.start125bp.bed"
    run:
        import pandas as pd
        df=pd.read_csv(input.bed,sep="\t",header=None,
                       names=["chr","start","end","TE","strand"])
        rec=[]; half=125
        for _,r in df.iterrows():
            if r.strand=="+":
                s=max(0,r.start-half); e=r.start+half
            else:
                s=max(0,r.end-half);   e=r.end+half
            rec.append((r.chr,s,e,r.TE))
        pd.DataFrame(rec,columns=["chr","start","end","TE"]) \
          .to_csv(output.bed,sep="\t",header=False,index=False)

rule sort_start_windows:
    input:  "bed/TOM-007_TEs.start125bp.bed"
    output: "bed/TOM-007_TEs.start125bp.sorted.bed"
    shell: "bedtools sort -i {input} > {output}"

##### Signal mapping #####
rule map_chip:
    input:
        bed="bed/TOM-007_TEs.upstream1kb.sorted.bed",
        bw=lambda wc: os.path.join(CHIP_DIRS[wc.tissue], f"{wc.sample}.bw")
    output: "results/{tissue}/{sample}_chip_signal.tsv"
    run:
        import pyBigWig,pandas as pd
        bwf=pyBigWig.open(input.bw)
        df=pd.read_csv(input.bed,sep="\t",header=None,
                       names=["chr","start","end","TE"])
        vals=[(r.TE, bwf.stats(r.chr,r.start,r.end,type="mean")[0] or 0)
              for _,r in df.iterrows()]
        bwf.close()
        pd.DataFrame(vals,columns=["TE",wildcards.sample]) \
          .to_csv(output[0],sep="\t",index=False)

rule map_rna:
    input:
        bed="bed/TOM-007_TEs.sorted.bed",
        bw=lambda wc: os.path.join(RNA_DIRS[wc.tissue], f"{wc.sample}.bw")
    output: "results/{tissue}/{sample}_RNA_signal.tsv"
    run:
        import pyBigWig,pandas as pd
        bw=pyBigWig.open(input.bw)
        df=pd.read_csv(input.bed,sep="\t",header=None,
                       usecols=[0,1,2,3],
                       names=["chr","start","end","TE"])
        vals=[(r.TE, bw.stats(r.chr,r.start,r.end,type="mean")[0] or 0)
              for _,r in df.iterrows()]
        bw.close()
        pd.DataFrame(vals,columns=["TE",wildcards.sample]) \
          .to_csv(output[0],sep="\t",index=False)

##### Combine & summarise #####
rule combine_all:
    input:
        chip=lambda wc: [f"results/{wc.tissue}/{s}_chip_signal.tsv" for s in CHIP_SAMPLES[wc.tissue]],
        rna =lambda wc: [f"results/{wc.tissue}/{s}_RNA_signal.tsv"    for s in RNA_SAMPLES[wc.tissue]]
    output: "results/{tissue}/final_matrix.tsv"
    run:
        import os,pandas as pd
        os.makedirs(os.path.dirname(output[0]),exist_ok=True)
        dfs_c=[pd.read_csv(x,sep="\t",index_col=0) for x in input.chip]
        dfs_r=[pd.read_csv(x,sep="\t",index_col=0) for x in input.rna]
        df_c=pd.concat(dfs_c,axis=1).groupby(lambda c:c.split("_")[0],axis=1).mean()
        df_r=pd.concat(dfs_r,axis=1).mean(axis=1).to_frame("RNA")
        pd.concat([df_c,df_r],axis=1).to_csv(output[0],sep="\t")

rule classify_context:
    input:
        bed="bed/TOM-007_TEs.sorted.bed",
        genes="bed/TOM-007_genes.bed"
    output: "results/TE_context.tsv"
    run:
        import pandas as pd,subprocess
        tmp="bed/tmp_TE.bed"
        subprocess.run(f"cp {input.bed} {tmp}",shell=True,check=True)
        subprocess.run(f"bedtools intersect -wa -a {tmp} -b {input.genes} "
                       "| cut -f4 | sort | uniq > results/_genic.txt",
                       shell=True,check=True)
        all_te=pd.read_csv(input.bed,sep="\t",header=None,usecols=[3],names=["TE"])
        genic=pd.read_csv("results/_genic.txt",header=None,names=["TE"])
        all_te["context"]=all_te.TE.isin(genic.TE).map({True:"genic",False:"intergenic"})
        all_te.to_csv(output[0],sep="\t",index=False)
        subprocess.run("rm bed/tmp_TE.bed results/_genic.txt",shell=True)

rule summarize_TE_signals:
    input:
        chip=lambda wc: expand("results/{tissue}/{sample}_chip_signal.tsv",
                               tissue=wc.tissue, sample=CHIP_SAMPLES[wc.tissue]),
        rna =lambda wc: expand("results/{tissue}/{sample}_RNA_signal.tsv",
                               tissue=wc.tissue, sample=RNA_SAMPLES[wc.tissue]),
        cont="results/TE_context.tsv"
    output: "results/{tissue}/TE_summary.tsv"
    run:
        import os,pandas as pd
        os.makedirs(os.path.dirname(output[0]),exist_ok=True)
        df_ctx=pd.read_csv(input.cont,sep="\t",index_col=0)
        df_act=pd.concat([pd.read_csv(x,sep="\t",index_col=0)
                          for x in input.chip if "H3K27ac" in x],axis=1).mean(axis=1)
        df_rep=pd.concat([pd.read_csv(x,sep="\t",index_col=0)
                          for x in input.chip if "H3K9me3" in x],axis=1).mean(axis=1)
        df_rna=pd.concat([pd.read_csv(x,sep="\t",index_col=0)
                          for x in input.rna],axis=1).mean(axis=1).to_frame("RNA")
        summary=pd.DataFrame({"context":df_ctx.context,
                              "H3K27ac_flank":df_act,
                              "H3K9me3_flank":df_rep}).join(df_rna).fillna(0)
        summary.to_csv(output[0],sep="\t")

##### Visualization #####
rule heatmap:
    input:  "results/{tissue}/final_matrix.tsv"
    output: "results/{tissue}/heatmap.png"
    script: "Python_scripts/plot_heatmap2.py"

rule heatmap_by_family:
    input:  "results/{tissue}/final_matrix.tsv"
    output: directory("results/{tissue}/family_heatmaps")
    script: "Python_scripts/plot_heatmap_by_family.py"

rule heatmap_rna_filtered:
    input:  "results/{tissue}/final_matrix.tsv"
    output:
        heatmap="results/{tissue}/heatmap_RNAseq_filtered.png",
        te_list="results/{tissue}/TE_list_rna_ge1.txt"
    script: "Python_scripts/plot_heatmap_RNAseq_filtered.py"

rule plot_expression_by_mark:
    input:  "results/{tissue}/TE_summary.tsv"
    output: "results/{tissue}/expression_by_mark.png"
    script: "Python_scripts/plot_expression_by_mark.py"

rule plot_expression_by_mark_expressed:
    input:  "results/{tissue}/TE_summary.tsv"
    output: "results/{tissue}/expression_by_mark_expressed.png"
    script: "Python_scripts/plot_expression_by_mark_expressed.py"

rule plot_expression_by_mark_interactive:
    input:  "results/{tissue}/TE_summary.tsv"
    output: "results/{tissue}/expression_by_mark.html"
    script: "Python_scripts/plot_expression_by_mark_interactive.py"

rule plot_expression_vs_mark:
    input:  summary="results/{tissue}/TE_summary.tsv"
    output: png="results/{tissue}/expression_vs_mark.png"
    script: "Python_scripts/plot_expression_vs_mark.py"

rule plot_TE_pies:
    input:  "results/{tissue}/TE_summary.tsv"
    output: "results/{tissue}/TE_piecharts.png"
    script: "Python_scripts/plot_TE_piecharts.py"

rule plot_corr_panels:
    input:  matrix="results/{tissue}/final_matrix.tsv"
    output: png="results/{tissue}/corr_panels.png"
    script: "Python_scripts/plot_corr_expressed.py"

rule scatter_mark_vs_expression:
    input:  matrix="results/{tissue}/final_matrix.tsv"
    output:
        "results/{tissue}/scatter_mark_vs_expr_H3K27ac.png",
        "results/{tissue}/scatter_mark_vs_expr_H3K9me3.png"
    script: "Python_scripts/plot_scatter_mark_vs_expression.py"

rule boxplot_quartiles:
    input:  matrix="results/{tissue}/final_matrix.tsv"
    output: "results/{tissue}/boxplot_quartiles.png"
    script: "Python_scripts/plot_boxplot_mark_quartile.py"

rule metaplot:
    input:  matrix="results/{tissue}/final_matrix.tsv"
    output: png="results/{tissue}/metaplot_signal.png"
    script: "Python_scripts/plot_metaplot.py"

rule scatter_rna_chip:
    input:  matrix="results/{tissue}/final_matrix.tsv"
    output: png="results/{tissue}/scatter_rna_chip.png"
    script: "Python_scripts/plot_rna_chip_scatter.py"

rule feature_importance:
    input:  matrix="results/{tissue}/final_matrix.tsv"
    output: png="results/{tissue}/feature_importance.png"
    script: "Python_scripts/plot_feature_importance.py"

rule plot_te_composition:
    input:
        matrix="results/{tissue}/final_matrix.tsv",
        classes="/store/EQUIPES/REPEAT/MEMBERS/louis.barolle/Exp_23/ressources/TEClasses.tsv"
    output: png="results/{tissue}/TE_composition.png"
    script: "Python_scripts/plot_TE_composition.py"


        
rule compare_gut_head:
    input:
        gut         = "results/gut/final_matrix.tsv",
        head        = "results/head/final_matrix.tsv",
        cls         = config["te_classes_tsv"]
    output:
        overlap     = "results/compare/gut_head_overlap.txt",
        venn        = "results/compare/gut_head_venn.png",
        composition = "results/compare/gut_head_overlap_composition.png"
    script:
        "Python_scripts/compare_TE_overlap.py"


rule motif_upstream_pipeline:
    input:
        bed    = "bed/TOM-007_TEs.upstream1kb.sorted.bed",
        genome = "/store/EQUIPES/REPEAT/MEMBERS/louis.barolle/Exp_23/TOM-007/Genome/TOM-007.fasta",
        motifs = "/store/EQUIPES/REPEAT/MEMBERS/louis.barolle/Exp_23/ressources/JASPAR/20250626162532_JASPAR2024_combined_matrices_1185657_meme.txt"
    output:
        fasta   = "results/TE_upstream1kb.fa",
        fimo    = "results/fimo_upstream1kb/fimo.tsv",
        density = "results/plots/fimo_upstream1kb_density.png"
    script:
        "Python_scripts/motif_upstream_pipeline.py"


