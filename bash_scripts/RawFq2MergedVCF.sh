#!/bin/bash
set -e

# Original Version Author: Romain Coppee
# Script Anthony Sleiman
# Creation date: 24/07/2025
# Last modification: 06/08/2025

### Prepare your working environment
# # ---------- First, Install the required tools if not already in your environment:
# # ----- If you don't have an environment yet, please create one as such:
# conda create -n ASEQv1 -c bioconda -c conda-forge \
#     bcftools multiqc bwa samtools snpeff picard
# conda activate ASEQv1
# # ----- if you have the environment but miss one or multiple tool, install the missing tools as such:
# conda install -c bioconda -c conda-forge \
#     bcftools \
#     multiqc \
#     bwa \
#     samtools \
#     snpeff \
#     picard
# echo "Tools installed in current active environment"

PICARD="/mnt/c/Stages/StageESCAPE/Tools/picard.jar"
SNPEFF=/mnt/c/Stages/StageESCAPE/Tools/snpEff/snpEff.jar

### Now that your environment is ready you can start analysing your genomes
#---------- Precise the needed external data (reference genome etc.) and their placment on your machine as such:
REF_GEN=/mnt/c/Stages/StageESCAPE/plasmodium_project/Pfalciparum.genome.fasta
GEN_CROSS=/mnt/c/Stages/StageESCAPE/PlasmoMutGIT/known_sites
FILES_ANNOT=/mnt/c/Stages/StageESCAPE/plasmodium_project/annotations











#---------- Set some parameters:
#----- Number of threads to use:
read -p "How many threads do you want to use? (default: 4) " threads
threads=${threads:-4}
#----- Choose the reference genome to use:
refgen=$(ls *.fasta | head -n1)
# Check if any .fasta file was found
if [ -z "$refgen" ]; then
    echo "No .fasta reference genome found in the directory: " $(pwd)
    exit 1
fi
echo "Found reference genome: $refgen"
# Ask if to use the found file
read -p "Use this reference? (y/n): " confirm
if [[ "$confirm" != "y" && "$confirm" != "Y" ]]; then
    echo "Available .fasta files:"
    ls *.fasta
    read -p "Enter the name of the reference genome file you want to use: " refgen
    # Check if the file exists
    if [ ! -f "$refgen" ]; then
        echo "File '$refgen' not found. Aborting."
        exit 1
    fi
fi
echo "Using reference genome: $refgen"
#----- Index the reference genome (only once, only if the indexing files are missing)
if [ ! -f "$refgen.bwt" ]; then
    echo "Indexing reference genome..."
    bwa index "$refgen"
else
    echo "BWA index already exists, skipping."
fi
if [ ! -f "$refgen.fai" ]; then
    samtools faidx "$refgen"
else
    echo "FASTA index already exists, skipping."
fi
if [ ! -f "${refgen%.fasta}.dict" ]; then
    echo "Creating sequence dictionary..."
    gatk CreateSequenceDictionary -R "$refgen" -O "${refgen%.fasta}.dict"
else
    echo "Sequence dictionary already exists, skipping."
fi
#---------- Now let's find the FASTQ pairs and treat them:
count=$(ls -1 *_R1.fastq | wc -l)
echo "Found $count FASTQ file pairs in the directory."
# Loop over all _R1.fastq files
mkdir -p quality_control
mkdir -p bam_files
mkdir -p Stats
for fq1 in *_R1.fastq;
do
    fq2="${fq1/_R1.fastq/_R2.fastq}"
    sample="${fq1%%_R1.fastq}"
    echo "Processing sample: $sample"
    if [ ! -f "$fq2" ]; then
        echo "Missing R2 file for $fq1 — skipping sample."
        continue
    fi
    #----- Quality control:
    fastqc -t "$threads" -o quality_control "$fq1" "$fq2"
    #----- Map our files to the reference genome:
    bwa mem -t "$threads" -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA\tLB:WGA\tPU:${sample}_unit" "$refgen" "$fq1" "$fq2" | \
        samtools view -b - | \
        samtools sort -@ "$threads" -o "bam_files/${sample}.sorted.bam"
    samtools index "bam_files/${sample}.sorted.bam"
    samtools flagstat "bam_files/${sample}.sorted.bam" > "Stats/${sample}.flagstat.txt"
    echo "Done with sample: $sample"
    echo "--------------------------------------------"
done
echo "All samples processed."










# ### Done, now next step is fixing the BAM files we got, let's go...
#---------- Complete Header to fix BAM files, then index them
#----- Check if the bam file has a ReadGroup Header, if not add one using PICARD tool:
for bam in bam_files/*.bam; do
    filename=$(basename "$bam")                  # e.g. 210307573.sorted.bam
    echo "121 file name: $filename"
    base="${filename%.bam}"                      # e.g. 210307573.sorted
    echo "123 base: $base"
    sample="${base%.sorted}"                     # e.g. 210307573
    echo "125 sample: $sample"
    echo "Checking read group in: $filename"
    if samtools view -H "$bam" | grep -q "^@RG"; then
        echo "Read group found. No action needed."
    else
        echo "No read group found. Rewriting BAM with RG info..."
        tmp_bam="bam_files/${base}.tmp.bam"
        echo "132 tmp_bam: $tmp_bam"
        java -jar "$PICARD" AddOrReplaceReadGroups \
            -I "$bam" \
            -O "$tmp_bam" \
            -LB WGA \
            -PL illumina \
            -PU NA \
            -SM "$sample"
        mv -f "$tmp_bam" "$bam"
        echo "Overwrote $filename with RG-added version"
    fi
    echo "Running MarkDuplicates on: $bam"
    dedup_bam="bam_files/${base}.dedup.bam"
    echo "146 dedup_bam: $dedup_bam"
    bai_dedup="${dedup_bam%.bam}.bai"
    echo "148 bai_dedup: $bai_dedup"
    java -jar "$PICARD" MarkDuplicates \
        -INPUT "$bam" \
        -OUTPUT "$dedup_bam" \
        -METRICS_FILE "bam_files/${base}.metrics.txt" \
        -REMOVE_DUPLICATES true \
        -VALIDATION_STRINGENCY LENIENT \
        -CREATE_INDEX true
    mv -f "$dedup_bam" "$bam"
    mv -f "$bai_dedup" "${bam}.bai"
    echo "Dedup done for: ${base}.bam"
    echo "------------------------------------------"
done
#----- Recalibrate base qualities in a .bam file so that quality metrics match actual observed error rates:
for b in bam_files/*.bam;
do
    echo "$b"
    sample=$(basename "$b")
    gatk BaseRecalibrator \
        -R $REF_GEN \
        -I "$b" \
        --known-sites "$GEN_CROSS/3d7_hb3.combined.final.vcf.gz" \
        --known-sites "$GEN_CROSS/hb3_dd2.combined.final.vcf.gz" \
        --known-sites "$GEN_CROSS/7g8_gb4.combined.final.vcf.gz" \
        -O "bam_files/${sample}.pass1.table"
    echo "BaseRecalibrator PROCESSED $b"
done
#----- Apply BQSR to input BAM file:
for f in bam_files/*.bam;
do
    echo "Input BAM: [$f]"
    echo "Reference Genome: [$REF_GEN]"
    echo "f is $f"
    file=$(basename "$f")
    echo "file is $file"
    gatk ApplyBQSR \
        -R $REF_GEN \
        -I "$f" \
        -bqsr "bam_files/${file}.pass1.table" \
        -O "bam_files/${file}.bam.fix"
done
#----- Replace the old bams with the fresh new generated ones:
echo "Removing .bam and .bam.bai files"
rm bam_files/*.bam
rm bam_files/*.bam.bai
echo ".bam and .bam.bai files removed"
for fixfile in bam_files/*.fix;
do
    mv -- "$fixfile" "${fixfile%.fix}"
    echo "Renamed $fixfile to ${fixfile%.fix}"
done
for baifile in bam_files/*.bam.fix.bai;
do
    newname="${baifile/.bam.fix.bai/.bam.bai}"
    mv -- "$baifile" "$newname"
    echo "Renamed $baifile to $newname"
done
# ----- Generate a complete variant calling from each fixed bam file:
for f in bam_files/*.bam;
do
    echo "f is $f"
    gatk HaplotypeCaller \
        -R $REF_GEN \
        -I "$f" \
        -ERC GVCF \
        --native-pair-hmm-threads 4 \
        --max-alternate-alleles 6 \
        -O "$f.g.vcf"
    echo "HaplotypeCaller PROCESSED $f"
done










# ### Great, now that we have then VCF files generated and placed in the bam_files folder, let's combine them
# #---------- Filter VCF file by mean depth
for f in bam_files/*.vcf;
do
    echo "f is $f"
    base=$(basename "$f" .vcf)
    echo "base is $base"
    vcftools --vcf "$f" \
        --min-meanDP 5 \
        --recode \
        --recode-INFO-all \
        --out "bam_files/${base}.fix"
    echo "file created: bam_files/${base}.fix"
    rm "$f"
done

for f in bam_files/*.vcf;
do 
    mv -- "$f" "${f%.fix.recode.vcf}.vcf"
    echo "fixing VCF PROCESSED $f"
done

VCF_FILES=()

for file in bam_files/*.g.vcf; do
    VCF_FILES+=(--variant "$file")
done

gatk CombineGVCFs \
    -R "$REF_GEN" \
    "${VCF_FILES[@]}" \
    -O bam_files/combine.vcf
echo "CombineGVCF Processed"
 
# Perform joint genotyping on one or more samples pre-called with HaplotypeCaller
gatk GenotypeGVCFs \
    -R $REF_GEN \
    -V bam_files/combine.vcf \
    --max-alternate-alleles 6 \
    -O bam_files/calling_GVCF.vcf
echo "GenotypeGVCFs PROCESSED"










# ### Perfect, let's now filter. Let's produce VCF files containing only high-quality SNPs from whole genome sequencing data
#---------- Build a recalibration model to score variant quality for filtering purposes
gatk VariantRecalibrator \
    -R $REF_GEN \
    -V bam_files/calling_GVCF.vcf \
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
    -O bam_files/recal_GVCF.vcf \
    --tranches-file bam_files/recal_GVCF.tranches \
    --rscript-file bam_files/recal_GVCF.plots.R
echo "VariantRecalibrator PROCESSED"
#---------- Apply a score cutoff to filter variants based on a recalibration table
gatk ApplyVQSR \
    -R $REF_GEN \
    -V bam_files/calling_GVCF.vcf \
    -mode SNP \
    --truth-sensitivity-filter-level 99.0 \
    --recal-file bam_files/recal_GVCF.vcf \
    --tranches-file bam_files/recal_GVCF.tranches \
    -O bam_files/applied_recal_GVCF.vcf
echo "ApplyRecalibration PROCESSED"
# ---------- Filtration: Hard-filter SNPs based on annotation thresholds - Added by GPT cause idk further testings needed!!!!
gatk VariantFiltration \
    -R $REF_GEN \
    -V bam_files/calling_GVCF.vcf \
    --filter-expression "QD < 2.0" --filter-name "LowQD" \
    --filter-expression "FS > 60.0" --filter-name "HighFS" \
    --filter-expression "SOR > 3.0" --filter-name "HighSOR" \
    --filter-expression "DP < 10" --filter-name "LowDP" \
    -O bam_files/filtered_GVCF.vcf
echo "VariantFiltration PROCESSED"
#---------- Annotate variants using snpEff"
java -jar $SNPEFF \
    -c $FILES_ANNOT/snpEff.config \
    -no-downstream \
    -no-upstream \
    -onlyProtein \
    Pf3D7v3 \
    bam_files/applied_recal_GVCF.vcf > bam_files/annotated.vcf   
echo "snpEff done"
#---------- Include core regions from Pf genetic crosses version 1
bcftools annotate \
    -a "annotations/regions-20130225.onebased.txt.gz" \
    -h "annotations/regions.hdr" \
    -Ov \
    -o bam_files/coreregions.vcf \
    -c CHROM,FROM,TO,RegionType bam_files/annotated.vcf
echo "Done"
#---------- Annotate global barcode SNPs from Neafsey et al., 2008"
bcftools annotate \
    -a annotations/global_barcode_tidy.txt.gz \
    -h annotations/global_barcode.hdr \
    -Ov \
    -o bam_files/barcodecoreregions.vcf \
    -c CHROM,FROM,TO,GlobalBarcode bam_files/coreregions.vcf
echo "Done"
#---------- Select only biallelic SNPs
gatk SelectVariants \
    -R $REF_GEN \
    -V bam_files/barcodecoreregions.vcf \
    -select-type SNP \
    --restrict-alleles-to BIALLELIC \
    -O bam_files/SNPs.vcf
echo "SelectVariants PROCESSED"
#---------- Annotate VCF file with additional filters at the variant level
gatk VariantFiltration \
    -R $REF_GEN \
    --filter-name LowQualVQ -filter "VQSLOD <= 0.0" \
    --filter-name NotCore -filter "RegionType != 'Core'" \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -V bam_files/SNPs.vcf \
    -O bam_files/SNP_filtered.vcf
echo "Done"
#---------- Clean VCF for analysis ready
gatk SelectVariants \
    -R $REF_GEN \
    -V bam_files/SNP_filtered.vcf \
    -O bam_files/SNP_filtered2.vcf \
    -select 'vc.isNotFiltered()'
echo "Done"
#---------- Exclude all genotypes with a filter flag not equal to "."(a missing value) or PASS
vcftools --vcf bam_files/SNP_filtered2.vcf --remove-filtered-all --recode --stdout > bam_files/pass_SNPs.vcf
echo "FIltering PROCESSED"

echo "ALL DONE, CONGRATS, now let's filter using R script... Good Luck ;)"
