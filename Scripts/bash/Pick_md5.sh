#!/bin/bash

### ============================================================
# Pick_md5.sh: A code snippet to get md5 names
# Author: Anthony Sleiman
# Creation date: 19/08/2025
# Last modification: 21/08/2025
### ============================================================

folder="/mnt/c/Stages/StageESCAPE/PlasmoMut/input/fastq"
output="/mnt/c/Stages/StageESCAPE/PlasmoMut/output/md5_summary.tsv"

echo -e "sample\tforward_file_name\tforward_md5_name\treverse_file_name\treverse_md5_name" > "$output"

for r1 in "$folder"/*_R1.fastq; do

    r2="${r1/_R1.fastq/_R2.fastq}"

    sample=$(basename "$r1" _R1.fastq)

    md5_r1=$(md5sum "$r1" | awk '{print $1}')
    md5_r2=$(md5sum "$r2" | awk '{print $1}')

    echo -e "${sample}\t$(basename "$r1")\t${md5_r1}\t$(basename "$r2")\t${md5_r2}" >> "$output"
done
