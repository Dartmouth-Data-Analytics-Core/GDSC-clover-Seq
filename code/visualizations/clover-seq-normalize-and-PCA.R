#~~~~~~~~~~~~~~~~~~~~~~~~ README ~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 
# Title: Full_tRNA_Genome_PCA.R
# Description: Generate PCA for full genome + tRNA counts
#
# Author: Mike Martinez
# Lab: GDSC
# Project: Clover-Seq
# Date created: 05/16/25
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
suppressMessages(library(DESeq2))
suppressMessages(library(stats))
suppressMessages(library(SummarizedExperiment))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# READ IN THE DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Arg 1: sample sheet
# Arg 2: reference level


#----- Set command line args
args <- commandArgs(trailingOnly = TRUE)

#----- Check that all arguments are supplied
if (length(args) < 2 | length(args) > 2) {
  stop("Usage: RScript PCA.R <Sample_list_SE.txt> <reference level>")
}

#----- Set variables based on command line args
metadata <- args[1]
refLevel <- args[2]

#----- Set input directories
trnaDir <- "03_tRNA_counts/"

#----- Set output directories
normalizedDir <- "05_normalized/"
pcaDir <- "06_PCA/"
rdsDir <- "07_rds_files/"

#---- Create output directories
if (!dir.exists(normalizedDir)) {
  dir.create(normalizedDir)
}

if (!dir.exists(pcaDir)) {
  dir.create(pcaDir)
}

if (!dir.exists(rdsDir)) {
  dir.create(rdsDir)
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
# READ IN DATA FILES FROM INPUT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Read in netadata
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

#----- Read in the data files
data <- read_file_safe(paste0(trnaDir, "gene_level_counts_collapsed.txt"), sep = "\t", row.names = 1)
trna <- read_file_safe(paste0(trnaDir, "tRNA_isotype_counts.txt"), sep = "\t", row.names = 1)


#----- Check all samples are shared between data and metadata
check1 <- all(rownames(meta) %in% colnames(data))
check2 <- all(rownames(meta) %in% colnames(trna))

#----- Stop if not
if (check1 != TRUE | check2 != TRUE){
    stop("Metadata samples do not match sample information in data!")
} else {
  message("All samples shared between counts and metadata!")
}

#----- Set colors
myColors <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2",
  "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7",
  "#9C755F", "#BAB0AC")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# CREATE DESEQ2 OBJECT FOR FULL TRNA + SMRNA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
message("--------------------------------------------------")
message("Building DESeq2 dataset for gene-level counts...\n")
ddsFull <- DESeq2::DESeqDataSetFromMatrix(countData = data,
                                  colData = meta,
                                  design = ~Group)
colData(ddsFull)$Group <- stats::relevel(colData(ddsFull)$Group, ref = refLevel)

#----- Run DESeq2 on full genome + tRNA counts
ddsFull <- DESeq2::DESeq(ddsFull)

#----- Save rds
saveRDS(ddsFull, file = paste0(rdsDir, "gene_level_DESeq2_object.Rds"))
message(paste0("\n\tSaved RDS to ", paste0(rdsDir, "gene_level_DESeq2_object.Rds")))

#----- Extract size factors
full_sizeFactors <- DESeq2::sizeFactors(ddsFull)
full_sizeFactors <- t(full_sizeFactors)
colnames(full_sizeFactors) <- NULL
full_sizeFactors <- as.vector(full_sizeFactors)
write.table(rbind(rownames(meta), full_sizeFactors), file = paste0(normalizedDir, "gene_level_counts_size_factors.csv"), row.names=FALSE, col.names=FALSE)
message(paste0("\tSize factors saved to ", paste0(normalizedDir, "gene_level_counts_size_factors.csv")))

#----- Extract normalized counts
normData <- DESeq2::counts(ddsFull, normalized = TRUE) %>%
  as.data.frame()
write.csv(normData, file = paste0(normalizedDir, "normalized_gene_level_counts.csv"))
message(paste0("\tNormalized counts saved to ", paste0(normalizedDir, "normalized_gene_level_counts.csv")))

#----- Perform rlog transformation
ddsFullRlog <- rlog(ddsFull) 
ddsFullRlogMat <- SummarizedExperiment::assay(ddsFullRlog) %>%
  as.matrix()

#----- Explore variance
varianceFull <- apply(ddsFullRlogMat, 1, var)
fullVars <- sort(varianceFull, decreasing = TRUE)

#----- Convert to dataframe
varianceFull <- data.frame(
  Features = seq_along(fullVars),
  Variance = fullVars
)

# Save the plot as a PNG
fullVar <- ggplot2::ggplot(varianceFull, aes(x = Features, y = Variance)) +
  geom_point() +
  theme_classic() +
  labs(
    title = "tRNA + smRNA Gene Variance Plot",
    x = "Number of Features (tRNAs + smRNAs)",
    y = "Variance") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12))
ggplot2::ggsave(fullVar, filename = paste0(pcaDir, "gene_level_variance_plot.png"), width = 6, height = 6)
message("\tPlotted gene_level_variance_plot.png\n")

#----- Generate PCs
generatePCs <- function(MAT, VARS, NFEATURES) {
  #----- Set a variable for the number of genes (features) to be used for PCA and clustering
  var_features_n <- NFEATURES
  #----- Order variance and select the rows (genes) with the most variance
  select <- order(VARS, decreasing = TRUE)[1:var_features_n]
  #----- Subset vsd values for genes by top variance ranks
  vsd_sub <- MAT[select,]
  #----- Transpose the matrix
  vsd_sub <- t(vsd_sub)
  #----- Run principal component analysis
  message(paste0("Running PCA on ", NFEATURES, " most variable features..."))
  pca <- prcomp(vsd_sub)
  #----- extract the variance explained by each PC
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  #----- subset for first 5 elements
  percentVar <- percentVar[1:5]
  #----- add names to the percentVar vector
  message("Percent variations:")
  percentVar <- paste(round(percentVar*100,2), "%", sep = "")
  names(percentVar) <- c("PC1", "PC2", "PC3", "PC4", "PC5")
  message(paste0(names(percentVar), sep = " "))
  message(paste0(percentVar, " "))
  cat("\n")
  #----- Construct a data frame with PC loadings and add sample labels
  pca_res <- list()
  pca_df <- as.data.frame(pca$x)
  pca_eigenvecs <- as.data.frame(pca$rotation[,1:3])
  pca_res[["Loadings"]] <- pca_df
  pca_res[["Eigenvectors"]] <- pca_eigenvecs
  pca_res[["percent_var"]] <- percentVar
  #----- Return PCA values
  return(pca_res)
}

#----- Calculate PCs and extract loadings
PCsFull <- generatePCs(ddsFullRlogMat, fullVars, 500)
loadingsFull <- PCsFull[[1]]
loadingsFull$Sample <- rownames(loadingsFull)

#----- Ensure data is in same order as metadata
loadingsFull <- loadingsFull[match(rownames(meta), rownames(loadingsFull)),]
loadingsFull$Group <- meta$Group

#----- Extract the rotations
write.csv(loadingsFull, file = paste0(pcaDir, "gene_level_loadings.csv"))
message(paste0("\tPCA loadings saved to ", paste0(pcaDir, "gene_level_loadings.csv")))

#----- Plot PCA
PCAplotFull <- ggplot2::ggplot(loadingsFull, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 5) +
  geom_text_repel() +
  labs(title = "tRNA + smRNA PCA", 
       x = paste0("PC1: ", PCsFull[[3]][1]),
       y = paste0("PC2: ", PCsFull[[3]][2])) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        legend.position = "right",
        panel.background = element_rect(fill = "white", color = NA),  
        plot.background = element_rect(fill = "white", color = NA),
        title = element_text(face = "bold"))
ggplot2::ggsave(PCAplotFull, file = paste0(pcaDir, "gene_level_PCA.png"), width = 6, height = 6)
message("\tPlotted gene_level_PCA.png")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# TRNA ONLY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
message("--------------------------------------------------")
message("Building DESeq2 dataset for tRNA-isoform counts...\n")
ddstrna <- DESeq2::DESeqDataSetFromMatrix(countData = trna,
                                          colData = meta,
                                          design = ~Group)
colData(ddstrna)$Group <- stats::relevel(colData(ddstrna)$Group, ref = refLevel)

#----- Run DESeq2 on trna genome + tRNA counts
ddstrna <- DESeq2::DESeq(ddstrna)
saveRDS(ddstrna, file = paste0(rdsDir, "tRNA_isotype_DESeq2_object.Rds"))
message(paste0("\n\tSaved RDS to ", paste0(rdsDir, "tRNA_isotype_DESeq2_object.Rds")))


#----- Extract size factors
trna_sizeFactors <- DESeq2::sizeFactors(ddstrna)
trna_sizeFactors <- t(trna_sizeFactors)
colnames(trna_sizeFactors) <- NULL
trna_SizeFactors <- as.vector(trna_sizeFactors)
write.table(rbind(rownames(meta), trna_sizeFactors), file = paste0(normalizedDir, "trna_isotype_counts_size_factors.csv"), row.names=FALSE, col.names=FALSE)
write.csv(trna_sizeFactors, file = paste0(normalizedDir, "tRNA_isotype_counts_size_factors.csv"))
message(paste0("\tSize factors saved to ", paste0(normalizedDir, "tRNA_isotype_counts_size_factors.csv")))

#----- Extract normalized counts
normData <- DESeq2::counts(ddstrna, normalized = TRUE) %>%
  as.data.frame()
write.csv(normData, file = paste0(normalizedDir, "normalized_tRNA_isotype_counts.csv"))
message(paste0("\tNormalized counts saved to ", paste0(normalizedDir, "normalized_tRNA_isotype_counts.csv")))

#----- Perform rlog transformation
ddstrnaRlog <- rlog(ddstrna) 
ddstrnaRlogMat <- SummarizedExperiment::assay(ddstrnaRlog) %>%
  as.matrix()

#----- Explore variance
variancetrna <- apply(ddstrnaRlogMat, 1, var)
trnaVars <- sort(variancetrna, decreasing = TRUE)

#----- Convert to dataframe
variancetrna <- data.frame(
  Features = seq_along(trnaVars),
  Variance = trnaVars
)

# Save the plot as a PNG
trnaVar <- ggplot2::ggplot(variancetrna, aes(x = Features, y = Variance)) +
  geom_point() +
theme_classic() +
  labs(
    title = "tRNA Isodecoder Variance Plot",
    x = "Number of Features (tRNAs + smRNAs)",
    y = "Variance") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12))
ggplot2::ggsave(trnaVar, filename = paste0(pcaDir, "tRNA_isotype_variance_plot.png"), width = 6, height = 6)
message("\tPlotted tRNA_isotype_variance_plot.png\n")

#----- Calculate PCs and extract loadings
message(paste0("Length Variance tRNA ", length(variancetrna)))

PCstrna <- generatePCs(ddstrnaRlogMat, trnaVars, nrow(variancetrna))
loadingstrna <- PCstrna[[1]]
loadingstrna$Sample <- rownames(loadingstrna)

#----- Ensure data is in same order as metadata
loadingstrna <- loadingstrna[match(rownames(meta), rownames(loadingstrna)),]
loadingstrna$Group <- meta$Group

#----- Extract the rotations
write.csv(loadingstrna, file = paste0(pcaDir, "tRNA_isotype_loadings.csv"))
message(paste0("PCA loadings saved to ", paste0(pcaDir, "tRNA_isotype_loadings.csv")))

#----- Plot PCA
PCAplottrna <- ggplot2::ggplot(loadingstrna, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 5) +
  geom_text_repel() +
  labs(title = "tRNA Isodecoder PCA", 
       x = paste0("PC1: ", PCstrna[[3]][1]),
       y = paste0("PC2: ", PCstrna[[3]][2])) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        legend.position = "right",
        panel.background = element_rect(fill = "white", color = NA),  
        plot.background = element_rect(fill = "white", color = NA),
        title = element_text(face = "bold"))
ggplot2::ggsave(PCAplottrna, file = paste0(pcaDir, "tRNA_isotype_PCA.png"), width = 6, height = 6)
message("\tPlotted tRNA_isotype_PCA.png")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# COMBINED PCA PLOT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Extract common legend for plots and remove from individual plots
commonLegend <- suppressMessages(cowplot::get_legend(PCAplotFull))
PCAplotFull <- PCAplotFull + theme(legend.position = "none")
PCAplottrna <- PCAplottrna + theme(legend.position = "none")

#----- Create a combined plot
combinedPCA <- cowplot::plot_grid(
  cowplot::plot_grid(PCAplotFull, PCAplottrna, ncol = 1, align = "v"),
  commonLegend,
  ncol = 2,
  rel_widths = c(1, 0.3)
)
ggplot2::ggsave(paste0(pcaDir, "PCA_Analysis_Summary.png"), width = 8, height = 8)
message("\tPlotted PCA_Analysis_Summary.png")


#----- Try and get wild
loadingsFull$Analysis <- c("tRNAs + smRNAs")
loadingstrna$Analysis <- c("tRNA Isodecoders")
loadingsCombined <- rbind(loadingsFull, loadingstrna)
PCAplotBoth <- ggplot2::ggplot(loadingsCombined, aes(x = PC1, y = PC2, color = Group, label = Sample, shape = Analysis)) +
  geom_point(size = 5) +
  geom_text_repel() +
  geom_line(aes(group = Sample), linetype = "dotted", color = "darkgray") +
  scale_shape_manual(values = c("tRNAs + smRNAs" = 16, "tRNA Isodecoders" = 2)) +
  labs(title = "", 
       x = paste0("PC1 tRNA: ", PCstrna[[3]][1], "\nPC1 tRNA + smRNA: ", PCsFull[[3]][1]),
       y = paste0("PC2 tRNA: ", PCstrna[[3]][2], "\nPC2 tRNA + smRNA: ", PCsFull[[3]][2])) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        legend.position = "right",
        panel.background = element_rect(fill = "white", color = NA),  
        plot.background = element_rect(fill = "white", color = NA))
ggplot2::ggsave(PCAplotBoth, file = paste0(pcaDir, "PCA_Direct_Comparison.png"), width = 8, height = 8)
message("\tPlotted PCA_Direct_Comparison.png")




