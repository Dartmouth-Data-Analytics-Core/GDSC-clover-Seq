#~~~~~~~~~~~~~~~~~~~~~~~~ README ~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 
# Title: testing mismatch code
# Description: Test code to produce mismatch heatmap
#
# Author: Mike Martinez
# Lab: Orellana
# Project: clover-Seq development
# Date created: 6/12/25
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# LOAD LIBRARIES AND SET PATHS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#----- Libraries
library(tidyr)
library(dplyr)
library(ggplot2)
#library(ggtext)

#----- Set command line args
args <- commandArgs(trailingOnly = TRUE)

#----- Check that all arguments are supplied
if (length(args) < 3 | length(args) > 3) {
  stop("Usage: RScript <mature_tRNA_mismatches.txt> <Sample_list_SE.txt> <reference level>")
}

#----- Set variables based on command line args
mismatchData <- args[1]
metadata <- args[2]
refLevel <- args[3]

#----- Set input directories
trnaDir <- "03_tRNA_counts/"

#----- Set output directories
mismatchDir <- "03_tRNA_counts/mismatch_figures/"

#---- Create output directories
if (!dir.exists(mismatchDir)) {
  dir.create(mismatchDir)
}

#----- Check all input directories exist
message("--------------------------------------------------")
message(paste("Input directories:", trnaDir, sep = "\n\t"))
message("Checking that all input directories exist...")
if (!dir.exists(trnaDir)) {
  stop(paste(trnaDir, " Does not exist or is empty!\n"))
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
# READ IN THE DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
mis <- read_file_safe(mismatchData, sep = "\t", row.names = NULL)
meta <- read_file_safe(metadata, sep = ",", row.names = NULL)

#----- Ensure reference level is in metadata
uniqueLevels <- unique(meta$Group)
if (! refLevel %in% uniqueLevels) {
  stop(paste0(refLevel, " is not a level in your metadata!"))
}

#----- Ensure reference level is first factor
meta$Group <- factor(meta$Group)
meta$Group <- relevel(meta$Group, ref = refLevel)

#----- Get a vector of unique groups in the data
groups <- unique(meta$Group)

#----- Group 1 is reference level
group1 <- meta[meta$Group == groups[1],]$Sample_ID
group2 <- meta[meta$Group == groups[2],]$Sample_ID

#----- Extract just the group1 values
misGroup1 <- mis[mis$Sample %in% group1,]
misGroup2 <- mis[mis$Sample %in% group2,]

#----- Group by feature and position and sum mismatches in misGroup1
summaryGroup1 <- misGroup1 %>%
  dplyr::group_by(Feature, position) %>%
  dplyr::summarise(Total_Mismatches = mean(mismatchedbases), .groups = "drop")

#----- Group by feature and position and sum mismatches in misGroups
summaryGroup2 <- misGroup2 %>%
  dplyr::group_by(Feature, position) %>%
  dplyr::summarise(Total_Mismatches = mean(mismatchedbases), .groups = "drop")

# Join the two summaries
mismatch_comparison <- dplyr::inner_join(
  summaryGroup1, summaryGroup2,
  by = c("Feature", "position"),
  suffix = c("_group1", "_group2")
)

# Subtract group2 from group1
mismatches <- mismatch_comparison %>%
  dplyr::mutate(Relative_Mismatches = Total_Mismatches_group2 - Total_Mismatches_group1) %>%
  as.data.frame()

#----- Summed mismatches relative to one another
mismatches <- as.data.frame(mismatches$Relative_Mismatches)
colnames(mismatches) <- c("Relative_Mismatches")
mismatches$Feature <- summaryGroup1$Feature
mismatches$Position <- summaryGroup1$position

#----- pivot wider the data
heatmap_matrix <- mismatches %>%
  tidyr::pivot_wider(
    names_from = Position,
    values_from = Relative_Mismatches,
    values_fill = list(Relative_Mismatches = 0))

#----- Convert to dataframe
heatmap_matrix <- as.data.frame(heatmap_matrix)
rownames(heatmap_matrix) <- heatmap_matrix$Feature
heatmap_matrix$Feature <- NULL

#----- Get annotation column
annotations <- as.data.frame(rownames(heatmap_matrix))
colnames(annotations) <- c("Feature")

#----- Set isoacceptors
annotations$Isoacceptor <- ifelse(
  grepl("^tRX-", annotations$Feature),
  "tRX",
  sub("^tRNA-([^-]+)-.*", "\\1", annotations$Feature))

#----- Clean annotations
rownames(annotations) <- annotations$Feature
annotations$Feature <- NULL

# Isoacceptor color palette
colors <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
         "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928",
         "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC",
         "#E5D8BD", "#FDDAEC", "#F2F2F2", "#1B9E77", "#D95F02", "#7570B3", "#E7298A")

#----- Assign colors to isoacceptors
unique_isoacceptors <- sort(unique(annotations$Isoacceptor))
iso_colors <- setNames(colors[seq_along(unique_isoacceptors)], unique_isoacceptors)
annoColors <- list(Isoacceptor = iso_colors)

#----- Plot pheatmap
summaryFile <- c("Z_scored_summed_mismatches_relative_to_reference.png")
title <- paste0("Z-scaled summed mismatches relative to ", refLevel)
pheatmap::pheatmap(as.matrix(heatmap_matrix),
         scale = "row",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         file = paste0(trnaDir, summaryFile),
         width = 12, 
         height = 10,
         color = colorRampPalette(c("dark blue", "white", "dark red"))(100), 
         annotation_colors = annoColors,
         annotation_row = annotations,
         main = title)
message(paste0("\tPlotted ", summaryFile))


#----- Prepare for loop (set isoacceptors in heatmap matrix)
heatmap_matrix$Isoacceptor <- ifelse(
  grepl("^tRX-", rownames(heatmap_matrix)),
  "tRX",
  sub("^tRNA-([^-]+)-.*", "\\1", rownames(heatmap_matrix)))

#----- Iterate through
isoAcceptors <- unique(heatmap_matrix$Isoacceptor)
for (i in isoAcceptors) {
  #----- Subset for just the i-th isoacceptor
  subset <- heatmap_matrix[heatmap_matrix$Isoacceptor == i,]
  subset$Isoacceptor <- NULL
  message(paste0("\t", i, "..."))
  
  #----- Get annotations
  anno <- as.data.frame(rownames(subset))
  colnames(anno) <- c("Feature")
  anno$Isoacceptor <- i
  rownames(anno) <- anno$Feature
  anno$Feature <- NULL
  
  #----- Set filename
  name <- paste0(i, "_heatmap.png")
  heatmapFile <- paste0(mismatchDir, name)

  #----- Determine clustering logic
  do_cluster <- nrow(subset) > 2
  
  #----- Plot pheatmap
  pheatmap::pheatmap(as.matrix(subset),
           scale = "row",
           cluster_rows = do_cluster,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           annotation_row = anno,
           filename = heatmapFile,
           color = colorRampPalette(c("dark blue", "white", "dark red"))(100), 
           width = 14, 
           height = 10)
  message(paste0("\tPlotted ", name))

}

message("Done!")
































