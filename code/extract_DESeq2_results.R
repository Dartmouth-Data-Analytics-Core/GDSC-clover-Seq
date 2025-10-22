#~~~~~~~~~~~~~~~~~~~~~~~~ README ~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 
# Title: extract_DESeq2_results.R
# Description: Extract DESeq2 results from Rds files
#
# Author: Mike Martinez
# Lab: GDSC
# Project: Clover-Seq
# Date created: 05/27/25
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# LOAD LIBRARIES AND SET PATHS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#----- Libraries
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggblend))
suppressMessages(library(ggrepel))
suppressMessages(library(cowplot))
suppressMessages(library(DESeq2))
suppressMessages(library(ashr))
suppressMessages(library(stats))
suppressMessages(library(SummarizedExperiment))

#----- Source visualization code
source("code/visualizations/clover-seq-DEGs.R")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# READ IN THE DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Arg 1: Path to RDS file directory

#----- Set command line args
args <- commandArgs(trailingOnly = TRUE)

#----- Check that all arguments are supplied
if (length(args) < 5 | length(args) > 5) {
  stop("Usage: RScript PCA.R <Path to Rds directory> <metadata> <reference level> <log2FC magnitude threshold> <padj threshold>")
}

#----- Set variables based on command line args
rdsDir <- args[1]
metadata <- args[2]
refLevel <- args[3]
log2FC_thresh <- as.numeric(args[4])
padj_thresh <- as.numeric(args[5])

#----- Set FC thresholds from argument 4
denBound <- abs(log2FC_thresh) * -1
numBound <- abs(log2FC_thresh)


#----- Set output directories
opDir <- "10_differential_expression/"
figDir <- paste0(opDir, "Figures/")

#---- Create output directories
if (!dir.exists(opDir)) {
  dir.create(opDir)
}

if(!dir.exists(figDir)) {
  dir.create(figDir)
}

#----- Check all input directories exist
message("--------------------------------------------------")
message(paste("Input directories:", rdsDir, sep = "\n\t"))
message("Checking that all input directories exist...")
if (!dir.exists(rdsDir)) {
    stop(paste(rdsDir, " Does not exist or is empty!\n"))
}
message("All directories found. Starting script...\n")
message(paste0("Log2 FC magnitude threshold set as: ", log2FC_thresh))
message(paste0("padj threshold set as: ", padj_thresh))
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

#----- Function to safely read in Rds files
read_rds_safe <- function(filepath) {
  result <- tryCatch(
    {
      data <- readRDS(filepath)
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
# READ IN METADATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
meta <- read_file_safe(metadata, sep = ",", row.names = NULL)
rownames(meta) <- meta$Sample_ID

#----- Ensure reference level is in metadata
uniqueLevels <- unique(meta$Group)
if (! refLevel %in% uniqueLevels) {
  stop(paste0(refLevel, " is not a level in your metadata!"))
}

#----- Ensure reference level is first factor
meta$Group <- factor(meta$Group)
meta$Group <- relevel(meta$Group, ref = refLevel)

#----- vector of factored group levels
sampleLevels <- levels(meta$Group)

#----- Assign non-refernece level
numerator <- sampleLevels[2]
message(paste0("\tNumerator: ", numerator))
message(paste0("\tDenominator: ", refLevel))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# READ IN DATA FILES FROM INPUT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
message("\nReading RDS files...\n")

#----- Read in the gene level Rds file
geneLevel <- read_rds_safe(paste0(rdsDir, "gene_level_DESeq2_object.Rds"))

#----- Read in the tRNA-isotype Rds file
isotypes <- read_rds_safe(paste0(rdsDir, "tRNA_isotype_DESeq2_object.Rds"))

#----- Function to assign description column to results
addDescription <- function(results, num, den, numThresh, denThresh, padj_thresh) {
    results$Group <- ifelse(results$log2FoldChange < denThresh & results$padj <= padj_thresh, paste0("Enriched in ", refLevel),
                        ifelse(results$log2FoldChange > numThresh & results$padj <= padj_thresh, paste0("Enriched in ", num), "ns"))

    results$Group <- factor(results$Group, levels = c(
      paste0("Enriched in ", refLevel),
      paste0("Enriched in ", num),
      "ns"
    ))

    return(results)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# GENE LEVEL RESULTS: SHRINKAGE 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Extract results from Rds
geneLevelCoef<- DESeq2::resultsNames(geneLevel)[2]

#----- shrink Results
#geneLevelData <- as.data.frame(results(geneLevel, name = geneLevelCoef))
geneLevelData <- DESeq2::lfcShrink(geneLevel,
                                      coef = geneLevelCoef,
                                      type = "ashr")


#----- Extract normalized counts
geneLevelCounts <- as.data.frame(counts(geneLevel, normalized = TRUE))

#----- Combine the shrink results and counts
geneLevelResults <- cbind(geneLevelData, geneLevelCounts)
geneLevelResults <- na.omit(geneLevelResults)

#----- Format description column
geneLevelResults <- addDescription(geneLevelResults, numerator, refLevel, numBound, denBound, padj_thresh)

#----- Save as csv
write.csv(geneLevelResults, file = paste0(opDir, "gene_level_DESeq2_shrink_Results.csv"))
message(paste0("Saved gene level results to ", paste0(opDir, "gene_level_DESeq2_shrink_Results.csv")))

#----- Plot MA
MAplot(geneLevelResults, "Gene Level MA Plot", "gene_level_MA_plot.png")
VolcanoPlot(geneLevelResults, "Gene Level Volcano Plot", "gene_level_volcano_plot.png")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ISOTYPE LEVEL RESULTS 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Extract results from Rds
isoLevelCoef <- DESeq2::resultsNames(isotypes)[2]

#----- shrink Results
isoLevelData <- DESeq2::lfcShrink(isotypes,
                                   coef = isoLevelCoef,
                                   type = "ashr")

#----- Extract normalized counts
isoLevelCounts <- as.data.frame(counts(isotypes, normalized = TRUE))

#----- Combine the shrink results and counts
isoLevelResults <- cbind(isoLevelData, isoLevelCounts)
isoLevelResults <- na.omit(isoLevelResults)

#----- Format description column
isoLevelResults <- addDescription(isoLevelResults, numerator, refLevel, numBound, denBound, padj_thresh)

#----- Save as csv
write.csv(isoLevelResults, file = paste0(opDir, "tRNA_isotype_DESeq2_shrink_Results.csv"))
message(paste0("Saved isotype level results to ", paste0(opDir, "tRNA_isotype_DESeq2_shrink_Results.csv")))

#----- Plot MA
MAplot(isoLevelResults, "tRNA Isotype MA Plot", "tRNA_isotype_MA_plot.png")
VolcanoPlot(geneLevelResults, "tRNA Isotype Volcano Plot", "tRNA_isotype_volcano_plot.png")





