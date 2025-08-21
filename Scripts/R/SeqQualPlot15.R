# --- --- --- --- --- --- --- -- Sequencing Quality Plotter --- --- --- --- --- 

# Author: Anthony Sleiman
# Creation date: 07/08/2025
# Last modification: 18/08/2025

# ------------------------------ Project Description ---------------------------

# R script for processing sequencing data and visualizing read mapping and quality.
# It generates 1) Sankey diagram of read processing steps. 2) Sample proportions of
# mapped and unmapped reads. 3) Distribution of reads mapped to specific target genes.

# ------------------------------ Required Libraries ----------------------------

library(rtracklayer) # read and 1 base transform .bed files
library(Rsamtools) # scan and read .bam files (scanBam) -> To count total reads
library(GenomicAlignments) # read only mapped reads from .bam files
library(networkD3) # create Sankey diagram
library(htmlwidgets) # save plot as html
library(webshot2) # capture html to save as png
library(tidyr) # pivot_longer function
library(tibble) # rownames_to_column
library(dplyr) # pipe, easier and cleanier code
library(ggplot2) # trace plots
library(stringr) # text on plots

# ############################## PROCESS DATA ##################################

## ----------------------------- Import Data -----------------------------------

# Import BED (0 based, will be automatically 1 based with import)
bed_gr <- import("New15Targets.bed")
bed_gr_df <- as.data.frame(bed_gr)
n_amplicon <- nrow(bed_gr_df) # Count amplicon number
# Load .bam files name
bam_files <- list.files("bam_files/", pattern="\\.sorted.bam.bam$", full.names=TRUE)
# Create empty dataframe to stock results bases
results <- data.frame(File = character(), 
                      Reads = integer(), 
                      Mapped = integer(), 
                      Unmapped = integer(), 
                      OnTarget = integer(),
                      OffTarget = integer(),
                      stringsAsFactors = FALSE)
# Counts table initialised with gene_id column
all_counts_df <- data.frame(gene_id = 1:n_amplicon)

## ----------------------------- Treat Data ------------------------------------

for (file in bam_files) {
  # Read .bam file
  scan <- scanBam(BamFile(file))
  # Read only mapped reads form .bam
  bam_gr <- granges(readGAlignments(file))
  bam_gr_df <- as.data.frame(bam_gr)
  
  # Count total read count
  Reads <- length(row.names(as.data.frame(scan)))
  # Count Mapped reads
  Mapped <- length(bam_gr)
  # Calculate Unmapped reads
  Unmapped <- Reads - Mapped
  
  # Find overlaps position .bam file with target BED
  hits <- findOverlaps(bam_gr, bed_gr)
  # Finds the actual regions overlapping (1|3 ; 2|14...) 
  ov <- pintersect(bam_gr[queryHits(hits)], bed_gr[subjectHits(hits)])
  # Get the length of each mapped read, this way if 1 read overlaps 2 regions, we chose the longest
  ov_width <- width(ov)
  # Create a data frame with the info per mapped read (multimapped included)
  hits_df <- data.frame(read_id = queryHits(hits), gene_id = subjectHits(hits),
                        ov_len  = ov_width)  
  # For each read, keep the row with the largest overlap
  hits_df <- hits_df[order(hits_df$read_id, -hits_df$ov_len), ]
  hits_df <- hits_df[!duplicated(hits_df$read_id), ]
  # Count number or reads mapped to each gene region
  counts <- table(factor(hits_df$gene_id, levels = 1:n_amplicon))
  
  # Count OnTarget reads
  OnTarget <- length(unique(hits_df$read_id))
  # Calculate OffTarget reads
  OffTarget <- Mapped - OnTarget
  
  # Save results
  results <- rbind(results, data.frame(File = basename(file), Reads = Reads, 
                                       Mapped = Mapped, Unmapped = Unmapped, 
                                       OnTarget = OnTarget, OffTarget = OffTarget))
  # Add this sample's counts as a column in all_counts_df
  all_counts_df[[basename(file)]] <- as.integer(counts)
}

## ----------------------------- Review and Export Results ---------------------

# Transpose results
results <- t(results)
colnames(results) <- results[1, ]
results <- results[-1, ]

# Change col names to get the actual amplicons name
rownames(all_counts_df) <- bed_gr_df$name
all_counts_df <- all_counts_df[ , -1]

# Merge rsults and gene counts into the results data frame
results <- rbind(results, all_counts_df)

# Extract a data frame with only the amplicons info
df_with_gene <- tibble::rownames_to_column(all_counts_df, "gene")
# Initialize a new group column with NA
df_with_gene$group <- NA_character_
# Assign groups based on gene names
df_with_gene$group[grepl("^Mdr1", df_with_gene$gene)] <- "Mdr1"
df_with_gene$group[grepl("^Dhps", df_with_gene$gene)] <- "Dhps"
df_with_gene$group[grepl("^Dhfr", df_with_gene$gene)] <- "Dhfr"
df_with_gene$group[grepl("^Cytb", df_with_gene$gene)] <- "Cytb"
df_with_gene$group[grepl("^Cor", df_with_gene$gene)]  <- "Cor"
df_with_gene$group[grepl("^Crt", df_with_gene$gene)]  <- "Crt"
df_with_gene$group[grepl("^K13", df_with_gene$gene)]  <- "K13"

# Sum within each group
numeric_cols <- sapply(df_with_gene, is.numeric)
summed_df <- aggregate(
  df_with_gene[, numeric_cols, drop = FALSE],
  by = list(group = df_with_gene$group),
  FUN = sum)
rownames(summed_df) <- summed_df$group
summed_df <- summed_df[, -1]

# Merge the only Genes data frame with the previous results data frame
results <- rbind(results ,summed_df)
# Remove ".sorted.bam" from all column names
colnames(results) <- sub("\\.sorted\\.bam$", "", colnames(results))

# Export to TSV
write.table(results, "results_summary.tsv", sep = "\t", quote = FALSE, row.names = TRUE)

# ############################## Sankey Diagram ################################

## ----------------------------- Load and Prepare data -------------------------

# Read data from .tsv
tsv <- read.table("results_summary.tsv", sep = "\t", header = TRUE, row.names = 1)
# Calculate total reads of all Samples
tsv_sum <- rowSums(tsv)

# Prepare data frame for sankey links
stainley <- data.frame(
  source = c("Reads", "Reads", "Mapped", "Mapped",
             "OnTarget", "OnTarget", "OnTarget", "OnTarget", "OnTarget", "OnTarget", "OnTarget",
             "Cor",
             "Crt", "Crt", "Crt", "Crt",
             "Cytb",
             "Dhfr", "Dhfr",
             "Dhps", "Dhps",
             "K13",
             "Mdr1", "Mdr1", "Mdr1", "Mdr1"),
  target = c("Mapped", "Unmapped", "OnTarget", "OffTarget",
             "Cor", "Crt", "Cytb", "Dhfr", "Dhps", "K13", "Mdr1",
             "Cor_Fragment3_G50E_R100K",
             "Crt_FragmentJeromeAsaap_C72R_M74I_N75D_E_K76I_T93S_H97Y_L_Q", "Crt_Fragment4_F145I", "Crt_Fragment6_Q271E", "Crt_Fragment10_R371I",
             "Cytb_Fragment4_I258M_Y268C_S_N_272_275",
             "Dhfr_Fragment1_pos16_C50R_N51I_C59R", "Dhfr_Fragment2_S108N_I164L",
             "Dhps_Fragment1_S436A_S_G437A", "Dhps_Fragment2_K540E_A581G_A613S",
             "K13_Fragment2_A675V_C580Y_P574L_V568G_R561H_R622I",
             "Mdr1_Fragment1_N86Y", "Mdr1_Fragment2_Y184F", "Mdr1_Fragment3_S1034C_R_N1042D", "Mdr1_Fragment4_D1246Y"),
  value  = c(tsv_sum["Mapped"], tsv_sum["Unmapped"], tsv_sum["OnTarget"], tsv_sum["OffTarget"],
             tsv_sum["Cor"], tsv_sum["Crt"], tsv_sum["Cytb"], tsv_sum["Dhfr"], tsv_sum["Dhps"], tsv_sum["K13"], tsv_sum["Mdr1"],
             tsv_sum["Cor_Fragment3_G50E_R100K"],
             tsv_sum["Crt_FragmentJeromeAsaap_C72R_M74I_N75D_E_K76I_T93S_H97Y_L_Q"], tsv_sum["Crt_Fragment4_F145I"], tsv_sum["Crt_Fragment6_Q271E"], tsv_sum["Crt_Fragment10_R371I"],
             tsv_sum["Cytb_Fragment4_I258M_Y268C_S_N_272_275"],
             tsv_sum["Dhfr_Fragment1_pos16_C50R_N51I_C59R"], tsv_sum["Dhfr_Fragment2_S108N_I164L"],
             tsv_sum["Dhps_Fragment1_S436A_S_G437A"], tsv_sum["Dhps_Fragment2_K540E_A581G_A613S"],
             tsv_sum["K13_Fragment2_A675V_C580Y_P574L_V568G_R561H_R622I"],
             tsv_sum["Mdr1_Fragment1_N86Y"], tsv_sum["Mdr1_Fragment2_Y184F"], tsv_sum["Mdr1_Fragment3_S1034C_R_N1042D"], tsv_sum["Mdr1_Fragment4_D1246Y"]
             ))

# Create nodes data frame - unique nodes
stainley_nodes <- data.frame(
  name = unique(c(as.character(stainley$source), as.character(stainley$target)))
)

# Match source/target names to node IDs
stainley$IDsource <- match(stainley$source, stainley_nodes$name) - 1
stainley$IDtarget <- match(stainley$target, stainley_nodes$name) - 1

my_color <- 'd3.scaleOrdinal()
  .domain(["Reads","Mapped","Unmapped","OnTarget","OffTarget","Cor","Crt","Cytb","Dhfr","Dhps","Hrp2","Hrp3","K13","Mdr1","Pm", "Dhfr_Fragment1_pos16_C50R_N51I_C59R", "Dhfr_Fragment2_S108N_I164L", "Mdr1_Fragment1_N86Y", "Mdr1_Fragment2_Y184F", "Mdr1_Fragment3_S1034C_R_N1042D", "Mdr1_Fragment4_D1246Y", "Crt_FragmentJeromeAsaap_C72R_M74I_N75D_E_K76I_T93S_H97Y_L_Q", "Crt_Fragment4_F145I", "Crt_Fragment6_Q271E", "Crt_Fragment10_R371I", "Dhps_Fragment1_S436A_S_G437A", "Dhps_Fragment2_K540E_A581G_A613S", "Cor_Fragment3_G50E_R100K", "K13_Fragment2_A675V_C580Y_P574L_V568G_R561H_R622I", "Cytb_Fragment4_I258M_Y268C_S_N_272_275"])
  .range(["black","gray","silver","steelblue","tomato","chocolate","coral","crimson","darkcyan","darksalmon","darkseagreen","darkslateblue","greenyellow","indigo","lightsteelblue", "darkcyan", "darkcyan", "indigo", "indigo", "indigo", "indigo", "coral", "coral", "coral", "coral", "darksalmon", "darksalmon", "chocolate", "greenyellow", "crimson"])'

## ----------------------------- Let's Sankey ----------------------------------

# Build Sankey network
p <- sankeyNetwork(
  Links = stainley, 
  Nodes = stainley_nodes,
  Source = "IDsource", 
  Target = "IDtarget",
  Value = "value", 
  NodeGroup = "name",
  colourScale = my_color,
  sinksRight = FALSE, 
  nodeWidth = 20, 
  fontSize = 14, 
  fontFamily = "Arial"
)
p

# Save Sankey as an HTML file
saveWidget(p, "sankey_plot.html", selfcontained = TRUE)
# Save Sankey as .png
webshot("sankey_plot.html", "sankey_plot.png", vwidth = 1200, vheight = 800)

# ############################## Bar Plot ######################################

## ----------------------------- Load and Prepare data -------------------------

# Read data from .tsv
tsv <- read.table("results_summary.tsv", sep = "\t", header = TRUE, row.names = 1)

# Prepare bar data
df_bar <- pivot_longer(
  rownames_to_column(tsv, var = "Category"),
  cols = -Category,
  names_to = "Sample",
  values_to = "Count"
)
bar_data <- df_bar %>%
  filter(Category %in% c("Mapped", "Unmapped")) %>%
  pivot_wider(names_from = Category, values_from = Count) %>%
  mutate(
    Total_reads = Mapped + Unmapped,
    Mapped_pct = Mapped / Total_reads * 100,
    Unmapped_pct = Unmapped / Total_reads * 100
  ) %>%
  pivot_longer(cols = c(Mapped, Unmapped), names_to = "Category", values_to = "Count") %>%
  mutate(Count_k = Count / 1000)  # convert to thousands for plotting

## ----------------------------- Let's Bar Plot --------------------------------

# Plot
p1 <- ggplot(bar_data, aes(x = Sample, y = Count_k, fill = Category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  geom_text(
    data = bar_data %>% filter(Category == "Mapped"),
    aes(label = paste0(round(Mapped_pct, 2), "%")),
    position = position_stack(vjust = 0.5),
    color = "white",
    size = 3
  ) +
  scale_fill_manual(values = c("Mapped" = "skyblue3", "Unmapped" = "grey")) +
  ylab(expression("No. reads (x10"^3*")")) +
  xlab("") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
  )
p1

# Save
ggsave("barplot.jpeg", plot = p1, width = 8, height = 6, dpi = 300)

# ############################## Scatter plot ###################################

## ----------------------------- Load and Prepare data -------------------------

# Prepare Scatter data
gene_data <- df_bar[!(df_bar$Category %in% c("Reads", "Mapped", "Unmapped", "OnTarget", "OffTarget")), ]
first_cor_row <- which(gene_data[[1]] == "Cor")[1]
gene_data <- gene_data[-(1:first_cor_row-1), ]

## ----------------------------- Let's Scatter Plot ----------------------------

# Plot
p2 <- ggplot(gene_data, aes(x = Count, y = Sample, color = Category)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 1000, linetype = "dotted", color = "black") +
  xlab(expression("No. reads overlapping target")) +
  ylab("") +
  scale_x_log10() +
  scale_color_manual(
    values = c(
      "Cor"   = "chocolate",
      "Crt"   = "coral",
      "Cytb"  = "#DC143C",
      "Dhfr"  = "darkcyan",
      "Dhps"  = "darksalmon",
      "K13"   = "greenyellow",
      "Mdr1"  = "#4B0082"
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
  )
p2

# Save 
ggsave("scatterplot.jpeg", plot = p2, width = 8, height = 6, dpi = 300)