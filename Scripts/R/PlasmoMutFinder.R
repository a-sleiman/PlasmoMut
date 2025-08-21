# --- --- --- --- --- --- --- -- Resistance Variant Finder --- --- --- --- --- -

# Created by: Anthony Sleiman
# Creation date: 30/07/2025
# Last modifications: 12/08/2025

# ------------------------------ Project Description ---------------------------

# R script to identify variants exclusive to treatment failure cases from a VCF

# ------------------------------ Required Libraries ----------------------------

# Install if missing:
# BiocManager::install("VariantAnnotation")

library(VariantAnnotation)
library(openxlsx)

# ------------------------------ Load log --------------------------------------

start_time <- Sys.time()
log_path <- "summary_report.txt"
log <- file(log_path, open = "wt")
sink(log, split = TRUE)

# ------------------------------ Load Input Data -------------------------------

cat("Loading input files...\n")

# Load sample classification table
label <- read.xlsx("input_data/labels.xlsx")

cat("Label file loaded. You have ", nrow(label), " samples.\n")

# Load merged multi-sample VCF file
vcf <- readVcf("output_data/bam_files/SNP_filtered.vcf")
nb_raw_vcf <- nrow(vcf)
cat("VCF file loaded.", nb_raw_vcf, "total variants.\n")

# Extract genotype data (GT field)
gt <- geno(vcf)$GT

# Extract sample names
samples <- colnames(gt)

# Match classification info to VCF sample order
label_vector <- label$Group[match(samples, label$Sample)]

# Verify label file
if (any(is.na(label_vector))) {
  stop("Some samples in the VCF not found in labels file.")
}

# Rename genotype columns to include group classification
colnames(gt) <- paste0(samples, "_", label_vector)

# ------------------------------ Encode Genotype Patterns ----------------------

# Identify columns corresponding to treatment groups
cols_Fail <- grep("_Fail$", colnames(gt))
cols_Success <- grep("_Success$", colnames(gt))

# Initialize genotype summary vectors
code_Fail <- character(nrow(gt))
code_Success <- character(nrow(gt))

cat("Encoding allele patterns by group...\n")

# For each variant, extract the second allele (third character) and concatenate
# the pattern across all samples for each group
for (i in 1:nrow(gt)) {
  chars_Fail <- character(length(cols_Fail))
  chars_Success <- character(length(cols_Success))
  
  # Extract characters from treatment failure samples
  for (j in seq_along(cols_Fail)) {
    val <- gt[i, cols_Fail[j]]
    chars_Fail[j] <- substr(val, 3, 3)
  }
  code_Fail[i] <- paste(chars_Fail, collapse = "")
  
  # Extract characters from treatment Successs samples
  for (j in seq_along(cols_Success)) {
    val <- gt[i, cols_Success[j]]
    chars_Success[j] <- substr(val, 3, 3)
  }
  code_Success[i] <- paste(chars_Success, collapse = "")
}

cat("Encoding complete for", nb_raw_vcf, "variants.\n")

# Append summary codes to genotype matrix
gt <- cbind(gt, code_Fail, code_Success)

# ------------------------------ Filter Candidate Variants ---------------------

cat("Filtering variants unique to the 'Fail' group...\n")

# Retain only variants with at least one alternate allele in the failure group
# and none in the Successs group
keep_rows <- grepl("1", gt[, "code_Fail"]) & !grepl("1", gt[, "code_Success"])
gt_filtered <- gt[keep_rows, ]

nb_gt_filtered <- nrow(gt_filtered)
cat("Done.", nb_gt_filtered, "variants retained with ALT alleles in 'Fail' group only.\n")
cat((nb_raw_vcf - nb_gt_filtered), "variants discarded at this step.\n")
# 
# # Filter by min 2 mutations in fail
# count_ones_fail <- integer(nrow(gt_filtered))  # empty vector to store counts
# for (i in 1:nrow(gt_filtered)) {
#   chars <- strsplit(gt_filtered[i, "code_Fail"], split = "")[[1]]
#   count_ones_fail[i] <- sum(chars == "1")
# }
# has_ones_success <- grepl("1", gt_filtered[, "code_Success"])
# keep_rows <- (count_ones_fail >= 2) & (!has_ones_success)
# gt_filtered_2 <- gt_filtered[keep_rows, ]

# ------------------------------ Regenerate Filtered VCF -----------------------

# Identify variant IDs to retain
variants_to_keep <- rownames(gt_filtered)
# variants_to_keep <- rownames(gt_filtered_2)

# Subset the original VCF to retain only candidate variants
vcf_filtered <- vcf[rownames(vcf) %in% variants_to_keep, ]

# Write the filtered VCF for downstream analysis
writeVcf(vcf_filtered, filename = "variants_Fail_uniques.vcf")

cat("Filtered VCF saved.\n")

# ------------------------------ Keep Only Missense Variants -------------------

cat("Filtering for missense variants...\n")

# Check the new vcf if annotated
if (!"ANN" %in% names(info(vcf_filtered))) {
  stop("VCF with no annotation.")
}

# load the annotation alone
ann <- info(vcf_filtered)$ANN

# load a false logical list as long as the annotations
missense_idx <- logical(length(ann))

# change the false to TURE if the annotation is a missense_variant mutation
for (i in seq_along(ann)) {
  if (any(grepl("missense_variant", ann[[i]]))) {
    missense_idx[i] <- TRUE
  }
}

# Filter VCF by keeping only the TRUE = missense_variant
vcf_missense <- vcf_filtered[missense_idx]

# Extract genotype data of this new vcf
gt_missense <- geno(vcf_missense)$GT

nb_miss_idx <- sum(missense_idx)
cat(nb_miss_idx, "missense variants retained.\n")
cat(nb_gt_filtered - nb_miss_idx, "variants discarded (non-missense).\n")

# Save VCF as a new .vcf
writeVcf(vcf_missense, filename = "missense_variant.vcf")

cat("Missense variant VCF saved.\n")

# ------------------------------ Summary Report --------------------------------

end_time <- Sys.time()
elapsed <- round(difftime(end_time, start_time, units = "secs"), 1)

cat("\n ----------------- Summary Report --------------------\n")
cat("Total variants in original VCF:        ", nb_raw_vcf, "\n")
cat("Variants after Fail-only filtering:    ", nb_gt_filtered, "(", round(100 * (nb_gt_filtered / nb_raw_vcf), 1), "% )\n")
cat("Variants after missense filtering:     ", nb_miss_idx, "(", round(100 * (nb_miss_idx / nb_raw_vcf), 1), "% )\n")
cat("Total variants lost (in 2 steps):      ", nb_raw_vcf - nb_miss_idx, "(", round(100 * (nb_raw_vcf - nb_miss_idx) / nb_raw_vcf, 1), "% )\n")
cat("Execution time:                        ", elapsed, "seconds\n")
cat("------------------------------------------------------\n")
cat("Done!\n")

sink()  # stop logging

# ------------------------------ Plot Mutations --------------------------------

### IL RESTE A METTRE VRAIMENT LE NOM DE LA MUTAITON DONC LA POSITION LE GENE ETC.... ON TROUVE CES INFO DANS LE VCF

only_mut <- readVcf("missense_variant.vcf")
only_mut_gt <- geno(only_mut)$GT
presence_absence <- function(gt) ifelse(gt == "0/0" | gt == "./.", 0, 1)
binary_matrix <- apply(only_mut_gt, c(1,2), presence_absence)
binary_df <- as.data.frame(t(binary_matrix))
upset(
  binary_df,
  nsets = ncol(binary_df),      # number of mutations
  sets.bar.color = "grey20",     # couleur des barres de sets à gauche
  main.bar.color = "steelblue",  # couleur des intersections
  keep.order = TRUE
)
