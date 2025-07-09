# vivalaTE

A Snakemake pipeline for integrated analysis of transposable elements (TEs), chromatin marks, and RNA-seq across multiple genomes and tissues on a compute cluster.

## Quickstart

```bash
git clone https://github.com/Louis-Barolle/vivalaTE.git
cd vivalaTE
snakemake --cores 4
```

All input files (genome BED/GFF, bigWig signals, config.yaml) live in the same directory as the Snakefile. Results will be written under <genome>/results/<tissue>/….

Configuration
Edit config.yaml in the project root:

yaml
base_dir: "/path/to/vivalaTE"
genomes: ["TOM-007", "AKA-017", …]
tissues: ["gut", "head", …]
chip_dir_template: "{base_dir}/{genome}/genomic_data_{tissue}/CHIP-seq"
rna_dir_template: "{base_dir}/{genome}/genomic_data_{tissue}/RNA-seq"
te_classes_tsv: "{base_dir}/ressources/TEClasses.tsv"
Ensure your .bw files are named <sample>.bw and located under the paths defined above.

Workflow
BED preparation: GFF ? TE and gene BEDs; create upstream windows.

Signal mapping: compute mean ChIP/RNA signal per TE.

Combine & summarize: build final_matrix.tsv and TE_summary.tsv.

Visualization: heatmaps, violins/boxplots, scatterplots, pie charts, metaplots, feature-importance, overlap analyses.

Scripts live in Python_scripts/ and are invoked via Snakemake’s script: directive.

Contact
Louis Barolle
louis.barolle@i2bc.paris-saclay.fr

License
MIT License — see LICENSE for details.