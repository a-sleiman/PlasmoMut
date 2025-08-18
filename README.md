# PlasmoMut: A Variant Discovery Pipeline
A bioinformatics pipeline to identify candidate variants associated with treatment failure in *Plasmodium falciparum*. From raw sequencing data to variant filtering, this project documents every step needed to reproduce the full analysis and uncover candidate mutations.

## Project Goal

To identify mutations exclusively present in *P. falciparum* isolates from treatment failure cases, but absent in treatment success cases. The ultimate aim is to generate a list of candidate variants potentially linked to resistance against antimalarial drugs (e.g., lumefantrine).

---

## Folder Structure

| Folder         | Description |
|----------------|-------------|
| `input_data/`  | Raw FASTQ, reference genome, sample metadata |
| `bash_pipeline/` | Bash scripts for read mapping & VCF generation |
| `R_scripts/`   | R scripts for filtering and visualization |
| `output_data/` | Merged VCFs, filtered candidate variants, plots |
| `reference_data/` | Necessary data, reference genome, project design, workflow diagrams, notes |

---

## Steps Overview

- Process raw FASTQ files to produce a merged VCF (ASEQ1.sh): fastq -> .bam + .vcf + .vcf_merged.vcf
- Identify candidate variants exclusive to treatment failure samples (PlasmoMutFinder.R): .vcf_merged.vcf -> candidate variants
- Visualize read mapping quality and variant distributions (SeqQualPlot15.R): .bam -> plot.html + plot.jpg

---

##  Pipeline Overview

1. **Sequencing Output**  
   - Illumina paired-end FASTQ files

2. **Read Mapping & Variant Calling** (`ASEQ1.sh`)  
   - Align reads to *P. falciparum* reference genome  
   - Generate merged VCF of all samples

3. **Variant Filtering** (`Mutfind.R`)  
   - Separate samples into Success vs. Failure groups  
   - Identify variants exclusive to the failure group  
   - Output filtered VCF and candidate table

4. **Quality Visualization** (`SeqQualPlot15.R`)  
   - Sankey diagram of read processing  
   - Mapping proportions per sample  
   - Distribution of reads over target genes

---

## Required Files

- `input_data/*.fastq/`: Raw sequencing reads  
- `input_data/Pfalciparum.genome.fasta/`: Reference genome FASTA + index files  
- `input_data/labels.xlsx`: Sample metadata with `Sample` and `Group` columns

---

##  Dependencies

### Bash & Bioinformatics Tools

- `bwa`, `samtools`, `GATK`, `bcftools` (for ASEQ1.sh)

### R Packages

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("VariantAnnotation", "vcfR", "rtracklayer", "Rsamtools", "GenomicAlignments"))
install.packages("openxlsx")
install.packages(c("openxlsx", "networkD3", "htmlwidgets", "webshot2", "tidyr", "tibble", "ggplot2", "dplyr", "SankeyDiagram", "stringr"))

---

## Disclaimer

This project was developed as part of a research internship. The VCF and metadata files were provided by the hosting laboratory. This pipeline is released for educational and reproducibility purposes, but does not include raw sequencing data due to privacy and size constraints.

---

## Citation

If you use or adapt this pipeline, please credit:
In Preparation

---
