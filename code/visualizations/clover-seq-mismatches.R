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
library(tidyr)
library(dplyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: RScript <mature_tRNA_mismatches.txt> <Sample_list_SE.txt> <reference level>")
}

mismatchData <- args[1]
metadata     <- args[2]
refLevel     <- args[3]

trnaDir     <- "03_tRNA_counts/"
mismatchDir <- "03_tRNA_counts/mismatch_figures/"

if (!dir.exists(mismatchDir)) dir.create(mismatchDir, recursive = TRUE)

message("--------------------------------------------------")
message(paste("Input directories:", trnaDir, sep = "\n\t"))
message("Checking that all input directories exist...")
if (!dir.exists(trnaDir)) stop(paste(trnaDir, "does not exist or is empty!\n"))
message("All directories found. Starting script...")
message("--------------------------------------------------")

read_file_safe <- function(filepath, sep, row.names) {
  tryCatch({
    data <- read.csv(filepath, sep = sep, row.names = row.names)
    message(paste0(filepath, " loaded successfully\n"))
    data
  }, error = function(e) {
    message(paste0("Error: Failed to load ", filepath, ": ", e$message))
    NULL
  })
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SPRINZL POSITION ORDER
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
sprinzl_order <- c(
  '-1','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16',
  '17','17a','18','19','20','20a','20b','21','22','23','24','25','26','27','28',
  '29','30','31','32','33','34','35','36','37','38','39','40','41','42','43','44',
  '45','e1','e2','e3','e4','e5','e6','e7','e8','e9','e10','e11','e12','e13','e14',
  'e15','e16','e17','e18','e19','46','47','48','49','50','51','52','53','54','55',
  '56','57','58','59','60','61','62','63','64','65','66','67','68','69','70','71',
  '72','73','74','75','76'
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# READ IN THE DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
mis  <- read_file_safe(mismatchData, sep = "\t", row.names = NULL)
meta <- read_file_safe(metadata,     sep = ",",  row.names = NULL)

uniqueLevels <- unique(meta$Group)
if (!refLevel %in% uniqueLevels) {
  stop(paste0(refLevel, " is not a level in your metadata!"))
}

meta$Group <- relevel(factor(meta$Group), ref = refLevel)
groups <- levels(meta$Group)

group1 <- meta[meta$Group == groups[1], ]$Sample_ID
group2 <- meta[meta$Group == groups[2], ]$Sample_ID

misGroup1 <- mis[mis$Sample %in% group1, ]
misGroup2 <- mis[mis$Sample %in% group2, ]

summaryGroup1 <- misGroup1 %>%
  dplyr::group_by(Feature, position) %>%
  dplyr::summarise(Total_Mismatches = mean(mismatchedbases), .groups = "drop")

summaryGroup2 <- misGroup2 %>%
  dplyr::group_by(Feature, position) %>%
  dplyr::summarise(Total_Mismatches = mean(mismatchedbases), .groups = "drop")

mismatch_comparison <- dplyr::inner_join(
  summaryGroup1, summaryGroup2,
  by = c("Feature", "position"),
  suffix = c("_group1", "_group2")
)

mismatches <- mismatch_comparison %>%
  dplyr::mutate(Relative_Mismatches = Total_Mismatches_group2 - Total_Mismatches_group1) %>%
  dplyr::select(Feature, Position = position, Relative_Mismatches) %>%
  as.data.frame()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# PIVOT AND ORDER COLUMNS BY SPRINZL POSITION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
heatmap_matrix <- mismatches %>%
  tidyr::pivot_wider(
    names_from  = Position,
    values_from = Relative_Mismatches,
    values_fill = 0
  ) %>%
  as.data.frame()

rownames(heatmap_matrix) <- heatmap_matrix$Feature
heatmap_matrix$Feature <- NULL

# Add any missing Sprinzl positions as zero-filled columns, then reorder
missing_cols <- setdiff(sprinzl_order, colnames(heatmap_matrix))
heatmap_matrix[, missing_cols] <- 0
heatmap_matrix <- heatmap_matrix[, sprinzl_order, drop = FALSE]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ANNOTATIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
annotations <- data.frame(
  Isoacceptor = ifelse(
    grepl("^tRX-", rownames(heatmap_matrix)),
    "tRX",
    sub("^tRNA-([^-]+)-.*", "\\1", rownames(heatmap_matrix))
  ),
  row.names = rownames(heatmap_matrix)
)

colors <- c(
  "#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C",
  "#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928",
  "#FBB4AE","#B3CDE3","#CCEBC5","#DECBE4","#FED9A6","#FFFFCC",
  "#E5D8BD","#FDDAEC","#F2F2F2","#1B9E77","#D95F02","#7570B3","#E7298A"
)
unique_isoacceptors <- sort(unique(annotations$Isoacceptor))
iso_colors  <- setNames(colors[seq_along(unique_isoacceptors)], unique_isoacceptors)
annoColors  <- list(Isoacceptor = iso_colors)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SUMMARY HEATMAP (ALL tRNAs)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
summaryFile <- paste0(mismatchDir, "Z_scored_summed_mismatches_relative_to_reference.png")
title       <- paste0("Z-scaled summed mismatches relative to ", refLevel)

pheatmap::pheatmap(
  as.matrix(heatmap_matrix),
  scale          = "row",
  cluster_rows   = FALSE,
  cluster_cols   = FALSE,
  show_rownames  = FALSE,
  filename       = summaryFile,
  width          = 12,
  height         = 10,
  color          = colorRampPalette(c("dark blue", "white", "dark red"))(100),
  annotation_colors = annoColors,
  annotation_row = annotations,
  main           = title
)
message(paste0("\tPlotted ", basename(summaryFile)))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# PER-ISOACCEPTOR HEATMAPS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
heatmap_matrix$Isoacceptor <- ifelse(
  grepl("^tRX-", rownames(heatmap_matrix)),
  "tRX",
  sub("^tRNA-([^-]+)-.*", "\\1", rownames(heatmap_matrix))
)

for (i in unique(heatmap_matrix$Isoacceptor)) {
  message(paste0("\t", i, "..."))

  subset <- heatmap_matrix[heatmap_matrix$Isoacceptor == i, ]
  subset$Isoacceptor <- NULL

  anno <- data.frame(Isoacceptor = rep(i, nrow(subset)))
  rownames(anno) <- rownames(subset)

  heatmapFile <- paste0(mismatchDir, i, "_heatmap.png")
  do_cluster  <- nrow(subset) > 2

  pheatmap::pheatmap(
    as.matrix(subset),
    scale          = "row",
    cluster_rows   = do_cluster,
    cluster_cols   = FALSE,
    show_rownames  = TRUE,
    annotation_row = anno,
    filename       = heatmapFile,
    color          = colorRampPalette(c("dark blue", "white", "dark red"))(100),
    width          = 14,
    height         = 10
  )
  message(paste0("\tPlotted ", basename(heatmapFile)))
}

message("Done!")
