# PlasmoMut: A Variant Discovery Pipeline
**PlamoMut** is a reproducible bioinformatics workflow to identify candidate variants associated with Artemisinin-based Combination Therapy (ACT) treatment failure in *Plasmodium falciparum*. From raw sequencing data to variant filtering, this project documents every step needed to reproduce the full analysis and uncover candidate mutations. The pipeline includes both a Bash workflow for variant discovery and an R script for filtering, summarizing, and visualizing missense variants.

---

## Project Goal

To identify mutations found exclusively in ACT treatment failure *P. falciparum* isolates, but absent in treatment success isolates. The ultimate aim is to generate a list of candidate variants potentially linked to antimalarial resistance.

---

## Features

- Process raw sequencing data to generate high-confidence SNPs
- Filter variants exclusive to treatment failure cases
- Focus on missense variants for functional analysis
- Generate summary reports and interactive UpSet plots
- Fully reproducible workflow with modular scripts

---

## Repository Structure

```
PlasmoMut/
├─ input/
│  └─ metadata/labels.xlsx
├─ output/
│  ├─ vcf/
│  └─ R/
├─ Scripts/
│  ├─ Bash/PlasmoVar.sh
│  └─ R/PlasmoMutFinder.R
└─ docs/
   └─ PlasmoMutFinder_User_Guide.md
```

## Workflow

1. **Input Data**
   - Raw FASTQ sequencing reads (paired-end, Illumina).
   - Reference genome (Pfalciparum.genome.fasta + index files).
   - Sample metadata (labels.xlsx).
   - Gene annotations (Pfalciparum.gff3 + resistance gene lists).
   - Known variant VCFs (known_sites/).
   - Target regions (*.bed).
2. **Mapping & Variant Calling (PlasmoVar.sh)**
   - Runs bwa, samtools, and GATK to:
     1. Map FASTQ → BAM.
     2. Mark duplicates & sort.
     3. Jointly call variants → merged VCF.
3. **Candidate Variant Filtering (PlasmoMutFinder.R)**
   - Reads merged VCF.
   - Uses labels.xlsx to split samples into Failure vs Success.
   - Retains variants exclusive to the failure group.
   - Exports missense_variant.vcf + results_summary.tsv.
4. **Quality Control Visualization (SeqQualPlot15.R)**
   - Generates plots:
     1. Sankey diagram of read filtering steps.
     2. Barplots of mapping proportions.
     3. Read distribution over targets.
   - Outputs sankey_plot.html/png, barplot.jpeg, scatterplot.jpeg.

---

## Quick Start

1. Clone the repository:

```bash
git clone https://github.com/yourusername/PlasmoMut.git
cd PlasmoMut
```

2. Generate filtered VCFs using the Bash script:

```bash
bash Scripts/Bash/PlasmoVar.sh
```

3. Run the R filtering and plotting pipeline:

```r
setwd("C:/Stages/StageESCAPE/PlasmoMut")  # adjust path
source("Scripts/R/PlasmoMutFinder.R")
```

4. Check output in `output/R/`:

* `variants_Fail_uniques.vcf`
* `missense_variant.vcf`
* `summary_report.txt`
* `upset_plot.jpeg`

---

## User Guide

Detailed instructions, workflow explanation, and tips are in the [PlasmoMutFinder User Guide](user_guide/01_Run_PlasmoVar.md).

---

## Disclaimer

This project was developed as part of a research internship. The VCF and metadata files were provided by the hosting lab `ESCAPE – EpidémioSurveillance et circulation des parasites dans les environnements – UR 7510`. This pipeline is released for educational and reproducibility purposes.
Raw FASTQ sequencing files are not included in this repository due to size.
Publicly available sequencing data for this study will be deposited at the European Nucleotide Archive (ENA).
Once available, the ENA project accession will be provided here for direct download.

---

## Contribution & Support

* Feel free to open issues for bugs, questions, or feature requests
* Pull requests are welcome for improvements

## Citation

In Preparation

---
