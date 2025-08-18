# PlasmoMut: A Variant Discovery Pipeline
A bioinformatics pipeline to identify candidate variants associated with treatment failure in *Plasmodium falciparum*. From raw sequencing data to variant filtering, this project documents every step needed to reproduce the full analysis and uncover candidate mutations.

## Project Goal

To identify mutations exclusively present in *P. falciparum* isolates from treatment failure cases, but absent in treatment success cases. The ultimate aim is to generate a list of candidate variants potentially linked to resistance against antimalarial drugs.

---

## Folder Structure

| Folder         | Description |
|----------------|-------------|
| `R_scripts/`   | R scripts for filtering and visualization |
| `annotations/` | Gene coordinates and annotations |
| `bash_scripts/` | Bash scripts for read mapping & VCF generation |
| `input_data/`  | Raw FASTQ, reference genome, sample metadata |
| `known_sites/` | VCF files of known variants used in variant recalibration |
| `output_data/` | Merged VCFs, filtered candidate variants, plots |

---

## Steps Overview

- Process raw FASTQ files to produce a merged VCF (ASEQ1.sh)
- Identify candidate variants exclusive to treatment failure samples (PlasmoMutFinder.R)
- Visualize read mapping quality and variant distributions (SeqQualPlot15.R)

---

##  Pipeline Overview

1. **Sequencing Output**  
   - Illumina paired-end FASTQ files

2. **Read Mapping & Variant Calling** (`ASEQ1.sh`)  
   - Align reads to *P. falciparum* reference genome  
   - Generate BAM and merged VCF files of all samples

3. **Variant Filtering** (`Mutfind.R`)  
   - Identify variants exclusive to the failure group  
   - Output filtered VCF and candidate table

4. **Quality Visualization** (`SeqQualPlot15.R`)  
   - Sankey diagram of read processing  
   - Mapping proportions per sample  
   - Distribution of reads over target genes

---

## Required Files

Before running the pipeline, make sure you have the following files organized in your project directory and your scripts:
- `input_data/*.fastq/`: Raw sequencing reads (ILLUMINA pared-end)
- `input_data/Pfalciparum.genome.fasta/`: Plasmodium falciparum reference genome in FASTA format
- - `input_data/annotations/`: gene annotation reference files (e.g., GFF3) used to map variants to genes and funtional regions
- `input_data/known_sites/`: reference variant databases (e.g., dbSNP) used for base quality score recalibration (BQSR) and variant filtering
- `input_data/labels.xlsx` – Excel file with sample information:  
  - `Sample` column: sample names matching .BAM bash-output files.  
  - `Group` column: `Succes` or `Failure` (treatment outcome labels).  
- `input_data/Newtargets15/`: BED file containing genomic target regions for on-target read analysis
- Any index files generated during preprocessing (`.bai`, `.vcf`, `.bam`).  
- Any intermediate files needed for mapping, filtering, or plotting steps.  

---

##  Dependencies

### Bash & Bioinformatics Tools

- `bwa`, `samtools`, `GATK`, `bcftools` (for ASEQ1.sh)

### R Packages

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c(
    "VariantAnnotation", 
    "vcfR", 
    "rtracklayer", 
    "Rsamtools", 
    "GenomicAlignments"
))
install.packages(c(
    "openxlsx", 
    "networkD3", 
    "htmlwidgets", 
    "webshot2", 
    "tidyr", 
    "tibble", 
    "ggplot2", 
    "dplyr", 
    "SankeyDiagram", 
    "stringr"
))
```

---

## Disclaimer

This project was developed as part of a research internship. The VCF and metadata files were provided by the hosting lab `ESCAPE – EpidémioSurveillance et circulation des parasites dans les environnements – UR 7510`. This pipeline is released for educational and reproducibility purposes, but does not include raw sequencing data due to privacy and size constraints.

---

## Citation

If you use or adapt this pipeline, please credit:
In Preparation

---
