# PlasmoMut: A Variant Discovery Pipeline
A bioinformatics pipeline to identify candidate variants associated with treatment failure in *Plasmodium falciparum*. From raw sequencing data to variant filtering, this project documents every step needed to reproduce the full analysis and uncover candidate mutations.

## Project Goal

To identify mutations exclusively present in *P. falciparum* isolates from treatment failure cases, but absent in treatment success cases. The ultimate aim is to generate a list of candidate variants potentially linked to resistance against antimalarial drugs (e.g., lumefantrine).

---

## Pipeline Overview

1. **Sequencing Output**  
   - Illumina paired-end FASTQ files from 200+ isolates

2. **Read Mapping & Variant Calling**  
   - Align reads to the *P. falciparum* reference genome  
   - GATK best practices to call variants (VCFs)

3. **Sample Labeling**  
   - Metadata table assigning each sample as `Success` or `Failure`

4. **Variant Filtering (R)**  
   - Identify variants exclusively present in the `Failure` group  
   - Export candidate list and modified VCF

---

## Disclaimer

This project was developed as part of a research internship. The VCF and metadata files were provided by the hosting laboratory. This pipeline is released for educational and reproducibility purposes, but does not include raw sequencing data due to privacy and size constraints.

---

## Citation

If you use or adapt this pipeline, please credit:
In Preparation
---
