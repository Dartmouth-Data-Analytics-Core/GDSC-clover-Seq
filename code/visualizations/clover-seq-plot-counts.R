#~~~~~~~~~~~~~~~~~~~~~~~~ README ~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 
# Title: clover-seq-plot-counts.R
# Description: Plot exploratory plots for smRNAs and tRNAs
#
# Author: Mike Martinez
# Lab: GDSC
# Project: Clover-Seq
# Date created: 05/18/25
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# LOAD LIBRARIES AND SET PATHS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#----- Libraries
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(cowplot))
suppressMessages(library(stats))

#----- Set command line args
args <- commandArgs(trailingOnly = TRUE)
metadata <- args[1]
refLevel <- args[2]

if (length(args) < 2 | length(args) > 2) {
  stop("Usage: RScript clover-seq-plot-count.R <Sample_list_SE.txt> <ref level>")
}

#----- Check input metadata exists
if (!file.exists(metadata)) {
  stop(paste0("Input metadata file: ", metadata, " does not exists!"))
}

#----- Set input directories
trnaDir <- "03_tRNA_counts/"
smrnaDir <- "04_smRNA_counts/"
normalizeDir <- "05_normalized/"

#----- Set output directory and create
opDir <- "08_plots/"
if (!dir.exists(opDir)) {
  dir.create(opDir)
}

#----- Check all input directories exist
message("--------------------------------------------------")
message(paste("Input directories:", trnaDir, smrnaDir, normalizeDir, sep = "\n\t"))
message("Checking that all input directories exist...")
dirsToCheck <- c(trnaDir, smrnaDir, normalizeDir)
for (i in dirsToCheck) {
  if (!dir.exists(i)) {
    stop(paste(i, " Does not exist or is empty!\n"))
  } 
}
message("All directories found. Starting script...")
message("--------------------------------------------------")

#----- Function to safely read in CSVs
read_file_safe <- function(filepath, sep, row.names) {
  result <- tryCatch(
    {
      data <- read.csv(filepath, sep = sep, row.names = row.names)
      message(paste0(filepath, " loaded successfully\n"))
      return(data)
    },
    error = function(e) {
      message(paste0("Error: Failed to load ", filepath, ": ", e$message))
      return(NULL)
    }
  )
  return(result)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# READ IN THE METADATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
meta <- read_file_safe(metadata, sep = ",")
colnames(meta)[1] <- c("Sample")

#----- Ensure all necessary columns are present before continuing
checkCols <- c("Sample", "Group")
for (i in checkCols) {
  if (! i %in% colnames(meta)) {
    stop(paste0(i, " column is missing from Metadata!"))
  }
}

#----- Ensure reference level is in metadata
uniqueLevels <- unique(meta$Group)
if (! refLevel %in% uniqueLevels) {
  stop(paste0(refLevel, " is not a level in your metadata!"))
}

#----- Ensure reference level is first factor
meta$Group <- factor(meta$Group)
meta$Group <- relevel(meta$Group, ref = refLevel)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# READ IN THE DATA - tRNA ISOTYPE COUNTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#trna <- read.csv(paste0(normalizeDir, "normalized_tRNA_isotype_counts.csv"), sep = ",")
trna <- read_file_safe(paste0(normalizeDir, "normalized_tRNA_isotype_counts.csv"), sep = ",")
rownames(trna) <- trna[,1]
colnames(trna)[1] <- c("tRNA")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# FORMAT DATA FOR PLOTTING
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Pivot the data
trnaLong <- trna %>%
  tidyr::pivot_longer(-c(tRNA), names_to = "Sample", values_to = "Count")

#----- Add column for isoacceptor information
trnaLong$Isoacceptor <- sub("^[^-]+-([^-]+)-.*$", "\\1", trnaLong$tRNA)

#----- Add column for codon information
trnaLong$Codon <- sub("^[^-]+-[^-]+-([^-]+)-.*$", "\\1", trnaLong$tRNA)

#----- Pivot the data back to wide format
trnaWide <- trnaLong %>%
  tidyr::unite("Isoacceptor_Codon", Isoacceptor, Codon, sep = "_") %>%
  tidyr::pivot_wider(names_from = Sample, values_from = Count, values_fill = 0)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# DETAILED FIGURE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Convert to a long format for ggplot2
plot_data <- trnaWide %>%
  tidyr::pivot_longer(-c(Isoacceptor_Codon, tRNA), names_to = "Sample", values_to = "Count") %>%
  tidyr::separate(Isoacceptor_Codon, into = c("Isoacceptor", "Codon"), sep = "_")

#----- Factor the Isoacceptors
plot_data$Isoacceptor <- factor(plot_data$Isoacceptor, 
                                   levels = c(
                                     "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", 
                                     "Glu", "Gly", "His", "Ile", "iMet", "Leu",
                                     "Lys", "Met", "Phe", "Pro", "SeC", "Ser",
                                     "Sup", "Thr", "Trp", "Tyr", "Val", "Und"))


#----- Plot detailed plot
isotype_plot <- ggplot2::ggplot(plot_data, aes(x = Sample, y = Count, fill = Sample)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  
  geom_jitter(aes(color = Isoacceptor), width = 0.2, size = 1.5, alpha = 0.8) +  
  facet_wrap(~ Codon, scales = "free", ) +  
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right", 
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold", hjust = 0.5)) +
  labs(title = "Grouped Boxplot of Normalized tRNA Isotype Counts by Sample & Anticodon",
    x = "Sample",
    y = "Count",
    color = "Isoacceptor",
    fill = "Sample")
ggplot2::ggsave(paste0(opDir, "Grouped_boxplot_norm_tRNA_isotypes_by_Sample_and_Anticodon.png"), isotype_plot, width = 14, height = 12)
message("\tPlotted Grouped_boxplot_norm_tRNA_isotypes_by_Sample_and_Anticodon.png\n")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# NORMALIZED COUNTS ABSOLUTE ABUNDANCES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Plot general plot normalized absolute abundances
normalizedCounts <- ggplot2::ggplot(plot_data, aes(x = Isoacceptor, y = Count, fill = Sample)) +
  geom_bar(stat = "identity") +
  labs(y = "Normalized Counts",
       x = "") +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "right", 
    legend.background = element_rect(fill = "white", color = NA),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA))
ggplot2::ggsave(paste0(opDir, "Isoacceptor_counts_by_sample_normalized.png"),
                normalizedCounts, 
                width = 12, height = 10)
message("\tPlotted Isoacceptor_counts_by_sample_normalized.png\n")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# NORMALIZED ISOACCEPTOR ABSOLUTE ABUNDANCES PER SAMPLE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Calculate total counts for each Iso-acceptor
plot_data_aa <- suppressMessages(plot_data %>%
  group_by(Sample, Isoacceptor) %>%
  summarize(total_AA = sum(Count))
)

#----- Add group information to plot data
plot_data_aa$Group <- meta$Group[match(plot_data_aa$Sample, meta$Sample)]

#----- Plot isoacceptor absolute abundances per samples
isoacceptors <- ggplot2::ggplot(plot_data_aa, aes(x = Sample, y = total_AA, fill = Isoacceptor)) +
  geom_bar(stat = "identity") +
  labs(y = "Normalized Counts",
       x = "") +
  theme_minimal(base_size = 16) +
  facet_grid(~Group, scales = "free") +
  theme(
    legend.position = "right", 
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA))
ggplot2::ggsave(paste0(opDir, "Isoacceptor_counts_normalized.png"),
                isoacceptors, width = 12, height = 10)
message("\tPlotted Isoacceptor_counts_normalized.png\n")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# PLOT RELATIVE ABUNDANCES OF CCA ENDS (NON-NORMALIZED)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
message("--------------------------------------------------")
message("Calculating relative abundance of CCA ends with non-normalize data...\n")

#----- Read in the data
#ccaNN <- read.csv(paste0(trnaDir, "tRNA_ends_counts.txt"), sep = "\t", row.names = NULL)
ccaNN <- read_file_safe(paste0(trnaDir, "tRNA_ends_counts.txt"), sep = "\t", row.names = NULL)
colnames(ccaNN)[1] <- c("tRNA")

#----- Get data in long format
ccaNNLong <- suppressMessages(ccaNN %>%
  tidyr::pivot_longer(-c(tRNA, end), names_to = "Sample", values_to = "Count") %>%
  group_by(Sample) %>%
  mutate(RelAbund = Count / sum(Count))
)

#----- Append metadata
ccaNNLong$Group <- meta$Group[match(ccaNNLong$Sample, meta$Sample)]

#----- Plot
endRelative <- ggplot2::ggplot(ccaNNLong, aes(x = Sample, y = RelAbund, fill = end)) +
  geom_bar(stat = "identity") +
  labs(x = "",
       y = "Absolute Count",
       fill = "End") +
  theme_minimal(base_size = 16) +
  facet_grid(~Group, scales = "free") +
  theme(
    legend.position = "right", 
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA))
ggplot2::ggsave(paste0(opDir, "CCA_ends_Relative_Abundances.png"),
                endRelative, width = 12, height = 10)
message("\tPlotted CCA_ends_Relative_Abundances.png\n")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# NORMALIZE DATA WITH SIZE FACTORS FOR ABSOLUTE ABUNDANCES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
message("--------------------------------------------------")

#----- Read in size factors
message("Calculating normalization using 03_tRNA_counts/tRNA_isotype_counts_size_factors.csv\n")
trna_sf <- read.table(paste0(normalizeDir, "tRNA_isotype_counts_size_factors.csv"), 
                      sep = ",", 
                      header = TRUE, 
                      check.names = FALSE, 
                      quote = "\"")

#----- Remove the first column
trna_sf[,1] <- NULL

#----- Get size factors as a named vector
sizeFactors <- as.numeric(unlist(trna_sf[1,]))

#----- Remove the first element (X, which is just 1)
#names(sizeFactors) <- colnames(trna_sf)
names(sizeFactors) <- meta$Sample
message("Size Factors: ")
print(sizeFactors)

#----- Check the number of columns in size factors matches the number of samples in metadata
if (length(sizeFactors) != length(unique(meta$Sample))) {
  stop("Number of size factor columns does not equal number of samples in metadata!")
} 

#----- Prepare the data (normalization)
message("\tNormalizing CCA end points with tRNA size factors...\n")
ccaNorm <- ccaNN
ccaNorm[,1:2] <- NULL

ccaNorm <- tryCatch(
  {
    result <- sweep(ccaNorm, 2, sizeFactors, "/")
    message("\tNormalization successful. Saving data to 05_normalized/CCA_ends_normalized.csv\n")
    result
  },
  error = function(e) {
    stop("Normalization failed: ", e$message)
  }
)

#----- Re-append the first 2 columns
ccaNorm$tRNA <- ccaNN$tRNA
ccaNorm$end <- ccaNN$end

#----- Write as a csv to the normalized data folder
write.csv(ccaNorm, file = "05_normalized/CCA_ends_normalized.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# PREPARE DATA AND PLOT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Get data in long format
ccaNormLong <- suppressMessages(ccaNorm %>%
  tidyr::pivot_longer(-c(tRNA, end), names_to = "Sample", values_to = "Count") %>%
  group_by(end, Sample) %>%
  summarize(total = sum(Count))
)

#----- Append metadata
ccaNormLong$Group <- meta$Group[match(ccaNormLong$Sample, meta$Sample)]

#----- Plot
endTypePlot <- ggplot2::ggplot(ccaNormLong, aes(x = Sample, y = total, fill = end)) +
  geom_bar(stat = "identity") +
  labs(x = "",
       y = "Absolute Count",
       fill = "End") +
  theme_minimal(base_size = 16) +
  facet_grid(~Group, scales = "free") +
  theme(
    legend.position = "right", 
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA))
ggplot2::ggsave(paste0(opDir, "CCA_ends_normalized_absolute_abundances.png"),
                endTypePlot, width = 12, height = 10)
message("\tPlotted CCA_ends_normalized_absolute_abundances.png\n")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SMRNA DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Read in the data
smRNA <- read_file_safe(paste0(smrnaDir, "smRNA_raw_counts_by_sample.txt"), sep = "\t", row.names = 1)
smRNA$Type <- rownames(smRNA)

#----- Calculate relative abundances and format
smRNAlong <- smRNA %>%
  tidyr::pivot_longer(-c(Type), names_to = "Sample", values_to = "Count") %>%
  dplyr::group_by(Sample) %>%
  dplyr::mutate(RelAbund = Count / sum(Count)) %>%
  dplyr::mutate(Type2 = dplyr::recode(Type, 
                        "miRNA" = "smRNA",
                        "misc_RNA" = "Misc RNA",
                        "Mt_rRNA" = "Mito RNA",
                        "Mt_tRNA" = "Mito RNA",
                        "other" = "Misc RNA",
                        "pretRNA" = "tRNA",
                        "ribozyme" = "Ribozyme",
                        "rRNA" = "rRNA",
                        "rRNA_pseudogene" = "rRNA",
                        "scaRNA" = "smRNA",
                        "snoRNA" = "smRNA",
                        "snRNA" = "smRNA",
                        "sRNA" = "smRNA",
                        "tRNA" = "tRNA",
                        "tRNA_antisense" = "tRNA")) %>%
  dplyr::mutate(Type2 = factor(Type2, levels = c("Mito RNA", "Ribozyme", "rRNA", "smRNA", "tRNA", "Misc RNA")))

#----- Add in Group information
smRNAlong$Group <- meta$Group[match(smRNAlong$Sample, meta$Sample)]

#----- Plot for higher-level relative abundances
a <- ggplot2::ggplot(smRNAlong, aes(x = Sample, y = RelAbund, fill = Type2)) +
  geom_bar(stat = "identity") +
  labs(y = "Relative Abundance",
       x = "",
       fill = "RNA-Subgroup") +
  facet_grid(~Group, scales = "free") +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "right", 
    legend.background = element_rect(fill = "white", color = NA),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA))


#----- Subset for just smRNA subtype and calculate relative abundance of granular classifications
smRNAtypes <- smRNAlong %>%
  dplyr::filter(Type2 == "smRNA") %>%
  dplyr::group_by(Sample) %>%
  dplyr::mutate(RelAbund = Count / sum(Count))

#----- Set color palette
pal <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00")

#----- Plot granular relative abundances for smRNA subtype
x <- ggplot2::ggplot(smRNAtypes, aes(x = Sample, y = RelAbund, fill = Type)) +
  geom_bar(stat = "identity") +
  geom_bar(stat = "identity") +
  labs(y = "Relative Abundance",
       x = "",
       fill = "",
       title = "smRNA Subgroup: Granular") +
  facet_grid(~Group, scales = "free") +
  scale_fill_manual(values = pal[1:5]) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "bottom", 
    legend.background = element_rect(fill = "white", color = NA),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA))

#----- Subset for just tRNA subtypes and calculate granular relative abundances
tRNAtypes <- smRNAlong %>%
  dplyr::filter(Type2 == "tRNA") %>%
  dplyr::group_by(Sample) %>%
  dplyr::mutate(RelAbund = Count / sum(Count))

#----- Plot granular relative abundances for tRNA subtypes
y <- ggplot2::ggplot(tRNAtypes, aes(x = Sample, y = RelAbund, fill = Type)) +
  geom_bar(stat = "identity") +
  geom_bar(stat = "identity") +
  labs(y = "Relative Abundance",
       x = "",
       fill = "",
       title = "tRNA Subgroup: Granular") +
  facet_grid(~Group, scales = "free") +
  scale_fill_manual(values = pal[6:8]) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "bottom", 
    legend.background = element_rect(fill = "white", color = NA),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA))
  
subTypePlot <- cowplot::plot_grid(x, y, ncol = 2)
full <- cowplot::plot_grid(a, subTypePlot, nrow = 2)
ggplot2::ggsave(paste0(opDir, "smRNA_Relative_Abundances.png"), full, width = 12, height = 10)
message("\tPlotted smRNA_Relative_Abundances.png")

message("--------------------------------------------------")
message("Count plotting finished!")
message("--------------------------------------------------")
