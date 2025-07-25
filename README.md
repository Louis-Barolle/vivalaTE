

````markdown
# vivalaTE

A Snakemake pipeline for integrated analysis of transposable elements (TEs), chromatin marks, and RNA-seq across multiple genomes and tissues on a compute cluster.

## Quickstart

```bash
git clone https://github.com/Louis-Barolle/vivalaTE.git
cd vivalaTE
snakemake --cores 4
````

*All input files (genome BED/GFF, bigWig signals, config.yaml) live in the same directory as the Snakefile. Results will be written under `results/<genome>/<tissue>/…`.*

## Configuration

All base data (Genome,TE,Gene annotations as well as CHIP seq and RNA-seq) were downloaded from https://www.biologiaevolutiva.org/gonzalez_lab/drosomics/DATA/ 
All 3 replicates were downloaded per experiment, and for CHIP-seq data, the file names were manually annotated to contain *_H3k9me3.bw or *_H3k27ac.bw

Edit `config.yaml` in the project root to suit your paths:

```yaml
base_dir: "/path/to/vivalaTE"
genomes:
  - "TOM-007"
  - "AKA-017"
  # …
tissues:
  - "gut"
  - "head"
  # …
chip_dir_template: "{base_dir}/{genome}/genomic_data_{tissue}/CHIP-seq"
rna_dir_template:  "{base_dir}/{genome}/genomic_data_{tissue}/RNA-seq"

te_gff_template:    "{base_dir}/{genome}/TE/{genome}.TE.gff3.gz"
gene_gff_template:  "{base_dir}/{genome}/Genes/{genome}.gff3.gz"
genome_fasta_template: "{base_dir}/{genome}/Genome/{genome}.fasta"

te_classes_tsv:     "{base_dir}/ressources/TEClasses.tsv"
```

Ensure your `.bw` files are named `<sample>.bw` and located under the directories defined above.

## Workflow

1. **BED preparation**

   * Convert GFF ? TE and gene BEDs
   * Extract upstream windows

2. **Signal mapping**

   * Compute mean ChIP/RNA signal per TE

3. **Combine & summarize**

   * Build `final_matrix.tsv` (all marks + RNA)
   * Build `TE_summary.tsv` (context + average signals)

4. **Visualization**

   * Heatmaps, violin/boxplots, scatterplots, pie charts, metaplots
   * Feature-importance (Random Forest)
   * Overlap analyses (venn diagrams, barplots)

All helper scripts live in `Python_scripts/` and are invoked automatically via Snakemake’s `script:` directives.

## Contact

Louis Barolle <i>[louis.barolle@i2bc.paris-saclay.fr](mailto:louis.barolle@i2bc.paris-saclay.fr)</i>

## License

This project is released under the MIT License. See [LICENSE](LICENSE) for details.

```
```
