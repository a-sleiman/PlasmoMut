#!/bin/bash
set -e

### ============================================================
# PlamoVar: A Variant Discovery Pipeline for Plasmodium falciparum
# From FASTQ to high-confidence annotated SNPs
### ============================================================

### ============================================================
# Original Version Author: Romain Coppee
# Current Script Author: Anthony Sleiman
# Creation date: 24/07/2025
# Last modification: 20/08/2025
### ============================================================

### ============================================================
### Prepare your working environment
### ============================================================
# Recommended: create a dedicated conda environment with all dependencies
# If you don’t already have one:
# conda create -n ASEQv1 -c bioconda -c conda-forge \
#     bcftools \
#     bwa \
#     samtools \
#     snpeff \
#     picard
# conda activate ASEQv1
#
# If the environment exists but you are missing a tool, install it with:
# conda install -c bioconda -c conda-forge <toolname>
#
# Tools documentation:
# - bcftools: https://samtools.github.io/bcftools/
# - BWA: http://bio-bwa.sourceforge.net/
# - samtools: http://www.htslib.org/
# - SnpEff: https://pcingola.github.io/SnpEff/
# - Picard: https://broadinstitute.github.io/picard/

### ============================================================
### Tool locations
### ============================================================
# If using conda: tools are in PATH, no need to specify full paths
# If using jar files manually downloaded, specify their paths here:

PICARD="/mnt/c/Stages/StageESCAPE/Tools/picard.jar"
SNPEFF="/mnt/c/Stages/StageESCAPE/Tools/snpEff/snpEff.jar"

### ============================================================
### External data
### ============================================================
# Reference genome and annotations (adjust paths as needed)

REF_GEN="/mnt/c/Stages/StageESCAPE/PlasmoMut/input/reference/Pfalciparum.genome.fasta"
GEN_CROSS="/mnt/c/Stages/StageESCAPE/PlasmoMut/input/known_sites"
FILES_ANNOT="/mnt/c/Stages/StageESCAPE/PlasmoMut/input/reference/annotations"
FASTQ_DIR="/mnt/c/Stages/StageESCAPE/PlasmoMut/input/fastq"

# Create output directories

QC_DIR="/mnt/c/Stages/StageESCAPE/PlasmoMut/output/qc"
BAM_DIR="/mnt/c/Stages/StageESCAPE/PlasmoMut/output/bam"
STATS_DIR="/mnt/c/Stages/StageESCAPE/PlasmoMut/output/stats"
FASTQ_DIR="/mnt/c/Stages/StageESCAPE/PlasmoMut/input/fastq"
VCF_DIR="/mnt/c/Stages/StageESCAPE/PlasmoMut/output/vcf"

mkdir -p "$QC_DIR" "$BAM_DIR" "$STATS_DIR" "$VCF_DIR"

echo "[INFO] Environment and paths set up successfully."

### ============================================================
### Set pipeline parameters
### ============================================================

# ----- Number of threads to use
read -p "How many threads do you want to use? (default: 4) " threads
threads=${threads:-4}
echo "[INFO] Using $threads threads."

# Choose the reference genome
# Check if the file exists
if [ ! -f "$REF_GEN" ]; then
    echo "[ERROR] Reference genome not found at $REF_GEN"
    exit 1
fi

echo "[INFO] Using reference genome: $REF_GEN"

### ============================================================
### Index the reference genome (only if missing)
### ============================================================

# BWA index
if [ ! -f "$REF_GEN.bwt" ]; then
    echo "[INFO] Indexing reference genome with BWA..."
    bwa index "$REF_GEN"
else
    echo "[INFO] BWA index already exists, skipping."
fi

# FASTA index
if [ ! -f "$REF_GEN.fai" ]; then
    echo "[INFO] Creating FASTA index with samtools..."
    samtools faidx "$REF_GEN"
else
    echo "[INFO] FASTA index already exists, skipping."
fi

# Sequence dictionary
if [ ! -f "${REF_GEN%.fasta}.dict" ]; then
    echo "[INFO] Creating sequence dictionary with GATK..."
    gatk CreateSequenceDictionary -R "$REF_GEN" -O "${REF_GEN%.fasta}.dict"
else
    echo "[INFO] Sequence dictionary already exists, skipping."
fi

### ============================================================
### FASTQ processing: QC + Mapping
### ============================================================

# Count FASTQ pairs
count=$(ls -1 "$FASTQ_DIR"/*_R1.fastq | wc -l)

echo "[INFO] Found $count FASTQ R1 files in the directory."

if [ "$count" -eq 0 ]; then
    echo "[ERROR] No FASTQ files found. Aborting."
    exit 1
fi

# Loop over all R1 FASTQ files
for fq1 in "$FASTQ_DIR"/*_R1.fastq; do
    fq2="${fq1/_R1.fastq/_R2.fastq}"
    sample=$(basename "$fq1" _R1.fastq)

    echo "[INFO] Processing sample: $sample"

    # Check R2 exists
    if [ ! -f "$fq2" ]; then
        echo "[WARNING] Missing R2 file for $fq1 — skipping sample."
        continue
    fi

    ### Quality control
    echo "[INFO] Running FastQC for $sample..."
    fastqc -t "$threads" -o "$QC_DIR" "$fq1" "$fq2"

    ### Alignment to reference genome
    echo "[INFO] Mapping reads for $sample with BWA..."
    bwa mem -t "$threads" -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA\tLB:WGA\tPU:${sample}_unit" "$REF_GEN" "$fq1" "$fq2" | \
        samtools view -b - | \
        samtools sort -@ "$threads" -o "$BAM_DIR/${sample}.sorted.bam"

    ### Index BAM and generate stats
    echo "[INFO] Indexing BAM for $sample..."
    samtools index "$BAM_DIR/${sample}.sorted.bam"

    echo "[INFO] Generating flagstat for $sample..."
    samtools flagstat "$BAM_DIR/${sample}.sorted.bam" > "$STATS_DIR/${sample}.flagstat.txt"

    echo "[INFO] Done with sample: $sample"
    echo "--------------------------------------------"
done

echo "[INFO] All samples processed."

### ============================================================
### Post-processing: Fix BAM headers + MarkDuplicates
### ============================================================

for bam in "$BAM_DIR"/*.sorted.bam; do
    filename=$(basename "$bam")            # e.g. 112313321571.sorted.bam
    base="${filename%.bam}"         # e.g. 112313321571.sorted
    sample="${base%.sorted}"      # e.g. 112313321571

    echo "[INFO] Checking read groups for: $filename"
    if samtools view -H "$bam" | grep -q "^@RG"; then
        echo "[INFO] Read group already present for $sample"
    else
        echo "[WARNING] No read group found, adding with Picard..."
        tmp_bam="$BAM_DIR/${base}.rg.bam"
        java -jar "$PICARD" AddOrReplaceReadGroups \
            -I "$bam" \
            -O "$tmp_bam" \
            -LB WGA \
            -PL ILLUMINA \
            -PU NA \
            -SM "$sample" \
            -CREATE_INDEX true
        mv -f "$tmp_bam" "$bam"
        echo "[INFO] Added read group to $filename"
    fi

    echo "[INFO] Running MarkDuplicates for $bam..."
    dedup_bam="$BAM_DIR/${base}.dedup.bam"
    java -jar "$PICARD" MarkDuplicates \
        -INPUT "$bam" \
        -OUTPUT "$dedup_bam" \
        -METRICS_FILE "$BAM_DIR/${base}.metrics.txt" \
        -REMOVE_DUPLICATES true \
        -VALIDATION_STRINGENCY LENIENT \
        -CREATE_INDEX true

    echo "[INFO] Deduplication complete for ${base}.bam: $dedup_bam"
done
rm -f "$BAM_DIR"/*.sorted.bam
rm -f "$BAM_DIR"/*.sorted.bam.bai

### ============================================================
### Base Quality Score Recalibration (BQSR)
### ============================================================

# Process only deduplicated BAM files
for b in "$BAM_DIR"/*.dedup.bam; do
    sample=$(basename "$b")
    echo "[INFO] Running BaseRecalibrator for $sample..."
    gatk BaseRecalibrator \
        -R "$REF_GEN" \
        -I "$b" \
        --known-sites "$GEN_CROSS/3d7_hb3.combined.final.vcf.gz" \
        --known-sites "$GEN_CROSS/hb3_dd2.combined.final.vcf.gz" \
        --known-sites "$GEN_CROSS/7g8_gb4.combined.final.vcf.gz" \
        -O "$BAM_DIR/${sample}.pass1.table"
    echo "[INFO] BaseRecalibrator finished for $sample"

done

### ============================================================
### Apply BQSR to deduplicated BAMs
### ============================================================

for b in "$BAM_DIR"/*.dedup.bam; do
    sample=$(basename "$b")
    echo "[INFO] Applying BQSR for $sample..."
    gatk ApplyBQSR \
        -R "$REF_GEN" \
        -I "$b" \
        -bqsr "$BAM_DIR/${sample}.pass1.table" \
        -O "$BAM_DIR/${sample}.bam.fix"
    echo "[INFO] BQSR applied for $sample"
done

### ============================================================
### Replace old dedup BAMs with recalibrated versions
### ============================================================

echo "[INFO] Replacing old dedup BAMs with BQSR-applied BAMs..."
rm -f "$BAM_DIR"/*.dedup.bam
rm -f "$BAM_DIR"/*.dedup.bai

for fixfile in "$BAM_DIR"/*.fix; do
    mv -- "$fixfile" "${fixfile%.fix}"
    echo "[INFO] Renamed $fixfile to ${fixfile%.fix}"
done

for baifile in "$BAM_DIR"/*.bam.fix.bai; do
    newname="${baifile/.bam.fix.bai/.bam.bai}"
    mv -- "$baifile" "$newname"
    echo "[INFO] Renamed $baifile to $newname"
done

### ============================================================
### Variant Calling with HaplotypeCaller
### ============================================================

for b in "$BAM_DIR"/*.bam; do
    sample=$(basename "$b")
    echo "[INFO] Running HaplotypeCaller for $sample..."
    gatk HaplotypeCaller \
        -R "$REF_GEN" \
        -I "$b" \
        -ERC GVCF \
        --native-pair-hmm-threads "$threads" \
        --max-alternate-alleles 6 \
        -O "$VCF_DIR/${sample}.g.vcf"
    echo "[INFO] HaplotypeCaller finished for $sample"
done

### ============================================================
### Filter VCFs by mean depth (min-meanDP 5)
### ============================================================

for f in "$VCF_DIR"/*.g.vcf; do
    echo "[INFO] Processing $f..."
    base=$(basename "$f" .g.vcf)
    vcftools --vcf "$f" \
        --min-meanDP 5 \
        --recode \
        --recode-INFO-all \
        --out "$VCF_DIR/${base}.fix"
    echo "[INFO] Filtered file created: $VCF_DIR/${base}.fix.recode"
    rm "$f"
done

### ============================================================
### Rename filtered VCFs to standard names
### ============================================================

for f in "$VCF_DIR"/*.fix.recode.vcf; do
    newname="${f%.fix.recode.vcf}.g.vcf"
    mv -- "$f" "$newname"
    echo "[INFO] $f renamed: $newname"
done

### ============================================================
### Combine GVCFs and Joint Genotyping
### ============================================================

VCF_FILES=()

# Collect all single-sample GVCFs
for file in "$VCF_DIR"/*.g.vcf; do
    VCF_FILES+=(--variant "$file")
    echo "[INFO] Added to combine list: $file"
done

# Combine GVCFs into a single multi-sample GVCF
echo "[INFO] Running CombineGVCFs..."
gatk CombineGVCFs \
    -R "$REF_GEN" \
    "${VCF_FILES[@]}" \
    -O "$VCF_DIR/combine.g.vcf"
echo "[INFO] CombineGVCFs completed: $VCF_DIR/combine.g.vcf"

# Joint genotyping to produce final multi-sample VCF
echo "[INFO] Running GenotypeGVCFs..."
gatk GenotypeGVCFs \
    -R "$REF_GEN" \
    -V "$VCF_DIR/combine.g.vcf" \
    --max-alternate-alleles 6 \
    -O "$VCF_DIR/calling_GVCF.vcf"
echo "[INFO] GenotypeGVCFs completed: $VCF_DIR/calling_GVCF.vcf"

## ============================================================
## Build recalibration model for SNPs
## ============================================================

echo "[INFO] Running VariantRecalibrator on calling_GVCF.vcf..."

gatk VariantRecalibrator \
    -R "$REF_GEN" \
    -V "$VCF_DIR/calling_GVCF.vcf" \
    --resource:cross1,known=false,training=true,truth=true,prior=15.0 "$GEN_CROSS/3d7_hb3.combined.final.vcf.gz" \
    --resource:cross2,known=false,training=true,truth=true,prior=15.0 "$GEN_CROSS/hb3_dd2.combined.final.vcf.gz" \
    --resource:cross3,known=false,training=true,truth=true,prior=15.0 "$GEN_CROSS/7g8_gb4.combined.final.vcf.gz" \
    -mode SNP \
    -an QD \
    -an FS \
    -an SOR \
    -an DP \
    --max-gaussians 8 \
    -mq-cap 70 \
    -O "$VCF_DIR/recal_GVCF.recal" \
    --tranches-file "$VCF_DIR/recal_GVCF.tranches"

echo "[INFO] VariantRecalibrator finished for calling_GVCF.vcf"

### ============================================================
### Apply VQSR to filter SNPs
### ============================================================

echo "[INFO] Applying VQSR to calling_GVCF.vcf..."

gatk ApplyVQSR \
    -R "$REF_GEN" \
    -V "$VCF_DIR/calling_GVCF.vcf" \
    -mode SNP \
    --truth-sensitivity-filter-level 99.0 \
    --recal-file "$VCF_DIR/recal_GVCF.recal" \
    --tranches-file "$VCF_DIR/recal_GVCF.tranches" \
    -O "$VCF_DIR/applied_recal_GVCF.vcf"

echo "[INFO] ApplyVQSR finished for calling_GVCF.vcf"

### ============================================================
### Annotate variants using snpEff
### ============================================================

echo "[INFO] Annotating variants in applied_recal_GVCF.vcf with snpEff..."

java -jar "$SNPEFF" \
    -c "$FILES_ANNOT/snpEff.config" \
    -no-downstream \
    -no-upstream \
    -onlyProtein \
    Pf3D7v3 \
    "$VCF_DIR/applied_recal_GVCF.vcf" > "$VCF_DIR/annotated.vcf"

echo "[INFO] snpEff annotation completed: annotated.vcf"

### ============================================================
### Annotate variants with core regions
### ============================================================

echo "[INFO] Annotating core regions on annotated.vcf..."

bcftools annotate \
    -a "$FILES_ANNOT/regions-20130225.onebased.txt.gz" \
    -h "$FILES_ANNOT/regions.hdr" \
    -Ov \
    -o "$VCF_DIR/coreregions.vcf" \
    -c CHROM,FROM,TO,RegionType "$VCF_DIR/annotated.vcf"

echo "[INFO] Core region annotation completed: $VCF_DIR/coreregions.vcf"

### ============================================================
### Annotate global barcode SNPs
### ============================================================

echo "[INFO] Annotating global barcode SNPs on coreregions.vcf..."

bcftools annotate \
    -a "$FILES_ANNOT/global_barcode_tidy.txt.gz" \
    -h "$FILES_ANNOT/global_barcode.hdr" \
    -Ov \
    -o "$VCF_DIR/barcodecoreregions.vcf" \
    -c CHROM,FROM,TO,GlobalBarcode "$VCF_DIR/coreregions.vcf"

echo "[INFO] Global barcode annotation completed: $VCF_DIR/barcodecoreregions.vcf"

### ============================================================
### Select only biallelic SNPs
### ============================================================

echo "[INFO] Selecting biallelic SNPs from barcodecoreregions.vcf..."

gatk SelectVariants \
    -R "$REF_GEN" \
    -V "$VCF_DIR/barcodecoreregions.vcf" \
    -select-type SNP \
    --restrict-alleles-to BIALLELIC \
    -O "$VCF_DIR/SNPs.vcf"

echo "[INFO] SelectVariants completed: $VCF_DIR/SNPs.vcf"

### ============================================================
### Apply additional variant-level filters
### ============================================================

echo "[INFO] Applying variant filtration on SNPs.vcf..."

gatk VariantFiltration \
    -R "$REF_GEN" \
    --filter-name LowQualVQ -filter "VQSLOD <= 0.0" \
    --filter-name NotCore -filter "RegionType != 'Core'" \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -V "$VCF_DIR/SNPs.vcf" \
    -O "$VCF_DIR/SNP_filtered.vcf"

echo "[INFO] Variant filtration completed: $VCF_DIR/SNP_filtered.vcf"

### ============================================================
### Keep only PASS variants
### ============================================================

echo "[INFO] Selecting only non-filtered variants..."

gatk SelectVariants \
    -R "$REF_GEN" \
    -V "$VCF_DIR/SNP_filtered.vcf" \
    -O "$VCF_DIR/SNP_filtered2.vcf" \
    -select 'vc.isNotFiltered()'

echo "[INFO] Non-filtered variant selection completed: $VCF_DIR/SNP_filtered2.vcf."
echo "[INFO] This VCF is ready for PlasmoMutFinder.R analysis."
echo "[INFO] You can continue with additional filtering if desired, but stopping here is recommended — this ensures you work with high-confidence SNPs."

### ============================================================
### Remove all genotypes flagged by filters
### ============================================================

echo "[INFO] Removing filtered genotypes to create final PASS SNPs VCF..."

vcftools --vcf "$VCF_DIR/SNP_filtered2.vcf" \
         --remove-filtered-all \
         --recode \
         --stdout > "$VCF_DIR/pass_SNPs.vcf"

echo "[INFO] Final filtered SNPs ready: $VCF_DIR/pass_SNPs.vcf"
echo "[INFO] ALL DONE — your high-quality SNPs VCF is ready for analysis."
