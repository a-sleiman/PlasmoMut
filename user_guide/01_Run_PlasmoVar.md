# PlasmoVar User Guide

This guide explains how to use the PlasmoVar pipeline step by step.

## 1. Overview

PlasmoVar is a reproducible pipeline for variant discovery in *Plasmodium falciparum* commonly known as Malaria.
It takes raw paired-end FASTQ files as input and produces high-confidence, annotated SNPs and merged VCF. These outputs are then used for downstream analysis with PlasmoMutFinder.R to identify candidate mutations potentially associated with ACT treatment failure.
This guide explains how to install dependencies, run the pipeline, and interpret outputs.

---

## 2. Requirements

* **Operating System**: Linux (this script was tested and created on Ubuntu 20.04+, but it may work on previous versions as well)
* **Dependencies** (must be in PATH - those are tools and packages we used in order to get the needed files and results and files):

  * `bwa`
  * `snpeff`
  * `samtools`
  * `picard`
  * `gatk`
  * `bcftools`
  * `fastqc`

* **Note on PICARD and SnpEff installation** In this script, PICARD and SnpEff are used in their .jar form (Java executables). This means that the script calls them via:

```
java -jar /path/to/picard.jar ...
java -jar /path/to/snpEff.jar ...
```

If you prefer to use the Bioconda versions of these tools, you can install them with:

```
conda install -c bioconda picard snpeff
```

After installing through Bioconda, you do not need to use java -jar anymore. The commands will instead look like this:

```
picard MarkDuplicates ...
snpEff ann ...
```

If you switch to Bioconda installations, remember to remove or edit the java -jar ... parts in the script. Both approaches will give the same results. Choose whichever fits your environment best.

* **Check installation with:**

```bash
bwa --version
samtools --version
gatk --version
```

---

## 3. Input Files

The pipeline expects the following structure inside `PlasmoMut/input/`:

- **fastq/** → Paired-end FASTQ files (`*_R1.fastq`, `*_R2.fastq`)
- **reference/** → Reference genome files:
  - `Pfalciparum.genome.fasta` (only required input; indexes are auto-generated)
  - `annotations/` (annotation files for filtering and barcode extraction)
    - `snpEff.config` → SnpEff configuration file
    - `regions-20130225.onebased.txt.gz` → Core genome regions for SNP filtering
    - `regions.hdr` → Header file for the regions file
    - `global_barcode_tidy.txt.gz` → Known barcode SNP positions
    - `global_barcode.hdr` → Header file for barcode SNP annotations
- **known_sites/** → Known variant sites used for recalibration
- **metadata/** → metadata for downstream R analysis

---

## 4. Running the Pipeline

Run the pipeline with:

```bash
bash PlasmoVar.sh
```

The script will after indexing the reference genome and quality check:

1. Perform read alignment (`bwa mem`)
2. Sort and index BAM files (`samtools sort`)
3. Mark duplicates (`picard MarkDuplicates`)
4. Base quality recalibration (`gatk BaseRecalibrator` / `ApplyBQSR`)
5. Call variants (`gatk HaplotypeCaller`)
6. Joint genotyping (`gatk CombineGVCFs` / `GenotypeGVCFs`)
7. Apply hard filters (`vcftools`)
8. Extract SNPs in core regions
9. Annotate global barcode SNPs

---

## 5. Output Files

All results are written to `PlasmoMut/output/`.

- **bam/** → Processed BAM files (sorted, deduplicated, recalibrated)
- **vcf/** → Variant calls
  - `<sample>.g.vcf` → Individual sample gVCFs
  - `combine.g.vcf` → Combined gVCF from all samples
  - `calling_GVCF.vcf` → Jointly genotyped variants
  - `SNP_filtered2.vcf` → Final filtered SNPs (main input for PlasmoMutFinder.R)
  - `pass_SNPs.vcf` → Strictly PASS-only variants
  - `coreregions.vcf` → Variants in core genome regions
  - `barcodecoreregions.vcf` → Final SNP barcode calls
- **qc/** → Quality control reports (`fastqc` HTML + TXT)
- **stats/** → Alignment statistics (`samtools flagstat`)

---

## 6. Example Workflow

```bash
# Place sequencer FASTQ files into /input/fastq
ls PlasmoMut/input/fastq

# Run pipeline
bash PlasmoVar.sh

# Final output needed for further analysis
less PlasmoMut/output/vcf_files/SNP_filtered2.vcf
```

---

## 7. Troubleshooting

If an error occurs, it will be displayed clearly in the terminal with [ERROR] or [WARNING] messages. Common issues include:
- No FASTQ files found → Make sure all input files end with _R1.fastq and _R2.fastq and are placed in PlasmoMut/input/fastq/.
- Reference genome missing → Check that $REF_GEN points to the correct .fasta file.
- Missing annotations or known sites → Ensure all files in PlasmoMut/input/reference/annotations/ and PlasmoMut/input/known_sites/ exist.
- Java tools errors (Picard/SnpEff) → Confirm the java -jar paths are correct, or use Bioconda-installed versions.
If all inputs are in place and paths are correct, the pipeline should run without issues.

---

## 8. Next Step

Once you have `PlasmoMut/output/vcf/SNP_filtered2.vcf`, proceed to `PlasmoMut/Scripts/R/PlasmoMutFinder.R` for candidate mutation analysis.

---

## 9. Support

For issues or questions, contact me or open an issue on the GitHub repository.