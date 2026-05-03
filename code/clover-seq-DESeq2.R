#~~~~~~~~~~~~~~~~~~~~~~~~ README ~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Title: clover-seq-DESeq2.R
# Description: Normalize tRNA counts, generate PCA, and run
#              differential expression using tRAX-style DESeq2
#              logic (estimateSizeFactors separately, sweep
#              normalization, betaPrior=TRUE, volcano plots).
#
# Author: Mike Martinez
# Lab: GDSC
# Project: Clover-Seq
# Date created: 04/30/26
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# LOAD LIBRARIES AND SET PATHS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(cowplot))
suppressMessages(library(DESeq2))
suppressMessages(library(stats))
suppressMessages(library(scales))
suppressMessages(library(SummarizedExperiment))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# PARSE ARGUMENTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Arg 1: sample sheet
# Arg 2: reference level
# Arg 3: (optional) comparisons file — two-column tab-delimited,
#         each row is one comparison (group1, group2).
#         If omitted, all pairwise combinations are run.

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2 | length(args) > 3) {
  stop("Usage: RScript clover-seq-DESeq2.R <Sample_list_SE.txt> <reference level> [comparisons.txt]")
}

metadata <- args[1]
refLevel <- args[2]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SET DIRECTORIES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
trnaDir       <- "03_Raw_Quant/tRNA_counts/"
normalizedDir <- "04_Expression/"
pcaDir        <- "07_Plots/PCA/"
rdsDir        <- "04_Expression/"
deDir         <- "04_Expression/Differential_Expression/"

if (!dir.exists(normalizedDir)) dir.create(normalizedDir, recursive = TRUE)
if (!dir.exists(pcaDir))        dir.create(pcaDir,        recursive = TRUE)
if (!dir.exists(rdsDir))        dir.create(rdsDir,        recursive = TRUE)
if (!dir.exists(deDir))         dir.create(deDir,         recursive = TRUE)

message("--------------------------------------------------")
message(paste("Input directories:", trnaDir, sep = "\n\t"))
message("Checking that all input directories exist...")
if (!dir.exists(trnaDir)) {
  stop(paste(trnaDir, " Does not exist or is empty!\n"))
}
message("All directories found. Starting script...")
message("--------------------------------------------------")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# HELPER FUNCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

read_file_safe <- function(filepath, sep, row.names) {
  tryCatch(
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
}

generatePCs <- function(MAT, VARS, NFEATURES) {
  select  <- order(VARS, decreasing = TRUE)[1:NFEATURES]
  vsd_sub <- t(MAT[select, ])
  message(paste0("Running PCA on ", NFEATURES, " most variable features..."))
  pca        <- prcomp(vsd_sub)
  percentVar <- pca$sdev^2 / sum(pca$sdev^2)
  percentVar <- percentVar[1:5]
  percentVar <- paste(round(percentVar * 100, 2), "%", sep = "")
  names(percentVar) <- c("PC1", "PC2", "PC3", "PC4", "PC5")
  message("Percent variations:")
  message(paste0(names(percentVar), sep = " "))
  message(paste0(percentVar, " "))
  cat("\n")
  list(
    Loadings     = as.data.frame(pca$x),
    Eigenvectors = as.data.frame(pca$rotation[, 1:3]),
    percent_var  = percentVar
  )
}

#----- Volcano plot y-axis helpers (from tRAX)
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv   <- function(x) base^(-x)
  scales::trans_new(
    paste0("reverselog-", format(base)), trans, inv,
    scales::log_breaks(base = base),
    domain = c(1e-100, Inf)
  )
}

reverselog_breaks <- function(minprob) {
  minprobadj <- min(c(minprob, 0.001))
  minfloor   <- 5 * (10^floor(log10(minprobadj)))
  minlog     <- log10(minfloor)
  10^(seq(minlog, 1))
}

#----- Differential expression runner (tRAX-style logic)
runDE <- function(dds, normCounts, meta, comparisons, typename, deDir) {

  message("--------------------------------------------------")
  message(paste0("Running differential expression for: ", typename))

  dashinterc <- 1.5
  outputFmt  <- ".pdf"

  #----- Results per comparison: extract log2FC (col 2) and padj (col 6)
  compareresults <- lapply(comparisons, function(currcompare) {
    res  <- DESeq2::results(dds,
                            contrast    = c("Group", currcompare[1], currcompare[2]),
                            cooksCutoff = TRUE)
    name <- paste(currcompare[1], currcompare[2], sep = "_")
    resdf <- as.data.frame(res)

    logcol        <- resdf[, 2, drop = FALSE]
    colnames(logcol) <- name
    padjcol       <- resdf[, 6, drop = FALSE]
    colnames(padjcol) <- name

    list(name = name, logFC = logcol, padj = padjcol)
  })

  alllogvals <- do.call(cbind, lapply(compareresults, `[[`, "logFC"))
  allprobs   <- do.call(cbind, lapply(compareresults, `[[`, "padj"))

  #----- Write dispersions
  write.table(
    cbind(rownames(normCounts), DESeq2::dispersions(dds)),
    file      = paste0(deDir, typename, "_dispersions.txt"),
    quote     = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  #----- Write padj and logval tables
  write.table(allprobs,   paste0(deDir, typename, "_padjs.txt"),   sep = "\t")
  write.table(alllogvals, paste0(deDir, typename, "_logvals.txt"), sep = "\t")
  message(paste0("\tWrote padjs and logvals for ", typename))

  #----- Volcano plot per comparison
  for (currcomp in comparisons) {
    currpair <- paste(currcomp[1], currcomp[2], sep = "_")
    pairname <- sub(":", "_", currpair, fixed = TRUE)

    currlogval <- alllogvals[, currpair]
    currpadj   <- allprobs[,   currpair]
    genename   <- rownames(allprobs)

    currsampledata <- data.frame(genename, currlogval, currpadj,
                                 stringsAsFactors = FALSE)

    validpadj    <- currsampledata$currpadj[!is.na(currsampledata$currpadj)]
    pvalcutoff   <- if (length(validpadj) >= 10) sort(validpadj)[10] else max(validpadj)
    displayfeats <- ifelse(
      abs(currsampledata$currlogval) > dashinterc &
        !is.na(currsampledata$currpadj) &
        currsampledata$currpadj < pvalcutoff,
      as.character(currsampledata$genename), ""
    )

    minpadj <- min(validpadj)
    if (!is.finite(minpadj) || minpadj <= 0) minpadj <- 1e-4

    currplot <- ggplot2::ggplot(currsampledata, aes(x = currlogval, y = currpadj)) +
      geom_point() +
      scale_x_continuous() +
      geom_text_repel(label = displayfeats,
                      min.segment.length = unit(0, "lines"),
                      segment.color = "red") +
      scale_y_continuous(trans   = reverselog_trans(10),
                         breaks  = reverselog_breaks(minpadj),
                         labels  = scales::scientific) +
      geom_hline(yintercept = 0.05,  linetype = 2) +
      geom_hline(yintercept = 0.005, linetype = 2) +
      geom_vline(xintercept =  dashinterc, linetype = 2) +
      geom_vline(xintercept = -dashinterc, linetype = 2) +
      theme_bw() +
      xlab("Log2-Fold Change") +
      ylab("Adjusted P-value") +
      ggtitle(paste(currpair, typename)) +
      theme(legend.box    = "horizontal",
            aspect.ratio  = 1,
            axis.text.x   = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      labs(caption = c(currcomp[1], currcomp[2])) +
      theme(plot.caption = element_text(size = 16, hjust = c(1, 0)))

    outfile <- paste0(deDir, typename, "_", pairname, "_volcano", outputFmt)
    ggplot2::ggsave(outfile, currplot)
    message(paste0("\tPlotted: ", basename(outfile)))
  }

  #----- Median normalized counts per group (tRAX-style)
  groupnames <- as.character(unique(meta$Group))
  medcounts  <- list()
  for (grp in groupnames) {
    cols <- rownames(meta)[as.character(meta$Group) == grp]
    if (length(cols) > 1) {
      medcounts[[grp]] <- apply(normCounts[, cols, drop = FALSE], 1, median)
    } else {
      medcounts[[grp]] <- normCounts[, cols]
    }
  }
  medcountmat <- do.call(cbind, medcounts)
  colnames(medcountmat) <- names(medcounts)
  write.table(medcountmat, paste0(deDir, typename, "_medians.txt"), quote = FALSE)
  message(paste0("\tWrote medians for ", typename))

  #----- Combined table: log2_* + padj_* + group medians
  logcols  <- alllogvals
  padjcols <- allprobs
  colnames(logcols)  <- paste("log2", colnames(logcols),  sep = "_")
  colnames(padjcols) <- paste("padj", colnames(padjcols), sep = "_")
  allcombinevals <- as.matrix(cbind(logcols, padjcols, medcountmat))
  write.table(allcombinevals, paste0(deDir, typename, "_combine.txt"), col.names = NA)
  message(paste0("\tWrote combined table for ", typename))
  message("--------------------------------------------------")
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# READ IN DATA FILES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
meta <- read_file_safe(metadata, sep = ",", row.names = NULL)
rownames(meta) <- meta$Sample_ID

uniqueLevels <- unique(meta$Group)
if (!refLevel %in% uniqueLevels) {
  stop(paste0(refLevel, " is not a level in your metadata!"))
}

meta$Group <- factor(meta$Group)
meta$Group <- relevel(meta$Group, ref = refLevel)

data <- read_file_safe(paste0(trnaDir, "gene_level_counts_collapsed.txt"), sep = "\t", row.names = 1)
trna <- read_file_safe(paste0(trnaDir, "tRNA_isotype_counts.txt"),         sep = "\t", row.names = 1)

check1 <- all(rownames(meta) %in% colnames(data))
check2 <- all(rownames(meta) %in% colnames(trna))
if (!check1 | !check2) {
  stop("Metadata samples do not match sample information in data!")
} else {
  message("All samples shared between counts and metadata!")
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# BUILD COMPARISONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
if (length(args) == 3) {
  pairtable  <- read.table(args[3], stringsAsFactors = FALSE)
  pairreduce <- pairtable[pairtable[, 1] %in% uniqueLevels &
                            pairtable[, 2] %in% uniqueLevels, ]
  comparisons <- apply(pairreduce, 1, list)
  comparisons <- lapply(comparisons, unlist)
} else {
  comparisons <- combn(as.character(unique(meta$Group)), 2, simplify = FALSE)
}

message(paste0("Comparisons to run: ", length(comparisons)))
for (comp in comparisons) message(paste0("\t", comp[1], " vs ", comp[2]))

#----- Set colors
myColors <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2",
              "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7",
              "#9C755F", "#BAB0AC")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# GENE-LEVEL (tRNA + smRNA)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
message("--------------------------------------------------")
message("Building DESeq2 dataset for gene-level counts...\n")

ddsFull <- DESeq2::DESeqDataSetFromMatrix(
  countData = data,
  colData   = meta,
  design    = ~Group
)
colData(ddsFull)$Group <- stats::relevel(colData(ddsFull)$Group, ref = refLevel)

#----- Estimate size factors separately (tRAX-style)
ddsFull <- DESeq2::estimateSizeFactors(ddsFull)

#----- Normalize via sweep (tRAX-style)
full_normalizedCounts <- sweep(data, 2, sizeFactors(ddsFull), "/")

#----- Write size factors (two-row: sample names, then values)
write.table(
  rbind(rownames(meta), sizeFactors(ddsFull)),
  file      = paste0(normalizedDir, "gene_level_counts_size_factors.csv"),
  row.names = FALSE,
  col.names = FALSE
)
message(paste0("\tSize factors saved to ", normalizedDir, "gene_level_counts_size_factors.csv"))

#----- Write normalized counts
write.csv(full_normalizedCounts, file = paste0(normalizedDir, "normalized_gene_level_counts.csv"))
message(paste0("\tNormalized counts saved to ", normalizedDir, "normalized_gene_level_counts.csv"))

#----- Run DESeq2 (tRAX-style: betaPrior=TRUE)
ddsFull <- DESeq2::DESeq(ddsFull, betaPrior = TRUE)

saveRDS(ddsFull, file = paste0(rdsDir, "gene_level_DESeq2_object.Rds"))
message(paste0("\n\tSaved RDS to ", rdsDir, "gene_level_DESeq2_object.Rds"))

#----- Differential expression
runDE(ddsFull, full_normalizedCounts, meta, comparisons, "gene_level", deDir)

#----- rlog for PCA
ddsFullRlog    <- DESeq2::rlog(ddsFull)
ddsFullRlogMat <- SummarizedExperiment::assay(ddsFullRlog) %>% as.matrix()

varianceFull <- apply(ddsFullRlogMat, 1, var)
fullVars     <- sort(varianceFull, decreasing = TRUE)

varianceFull <- data.frame(Features = seq_along(fullVars), Variance = fullVars)

fullVar <- ggplot2::ggplot(varianceFull, aes(x = Features, y = Variance)) +
  geom_point() +
  theme_classic() +
  labs(title = "tRNA + smRNA Gene Variance Plot",
       x     = "Number of Features (tRNAs + smRNAs)",
       y     = "Variance") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text  = element_text(size = 12))
ggplot2::ggsave(fullVar, filename = paste0(pcaDir, "gene_level_variance_plot.png"), width = 6, height = 6)
message("\tPlotted gene_level_variance_plot.png\n")

PCsFull      <- generatePCs(ddsFullRlogMat, fullVars, 500)
loadingsFull <- PCsFull[[1]]
loadingsFull$Sample <- rownames(loadingsFull)
loadingsFull <- loadingsFull[match(rownames(meta), rownames(loadingsFull)), ]
loadingsFull$Group <- meta$Group

write.csv(loadingsFull, file = paste0(pcaDir, "gene_level_loadings.csv"))
message(paste0("\tPCA loadings saved to ", pcaDir, "gene_level_loadings.csv"))

PCAplotFull <- ggplot2::ggplot(loadingsFull, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 5) +
  geom_text_repel() +
  labs(title = "tRNA + smRNA PCA",
       x     = paste0("PC1: ", PCsFull[[3]][1]),
       y     = paste0("PC2: ", PCsFull[[3]][2])) +
  theme_classic(base_size = 16) +
  theme(axis.text        = element_text(size = 14, face = "bold"),
        axis.title       = element_text(size = 16, face = "bold"),
        legend.position  = "right",
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA),
        title            = element_text(face = "bold"))
ggplot2::ggsave(PCAplotFull, file = paste0(pcaDir, "gene_level_PCA.png"), width = 6, height = 6)
message("\tPlotted gene_level_PCA.png")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# TRNA ISOTYPE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
message("--------------------------------------------------")
message("Building DESeq2 dataset for tRNA-isoform counts...\n")

ddstrna <- DESeq2::DESeqDataSetFromMatrix(
  countData = trna,
  colData   = meta,
  design    = ~Group
)
colData(ddstrna)$Group <- stats::relevel(colData(ddstrna)$Group, ref = refLevel)

#----- Estimate size factors separately (tRAX-style)
ddstrna <- DESeq2::estimateSizeFactors(ddstrna)

#----- Normalize via sweep (tRAX-style)
trna_normalizedCounts <- sweep(trna, 2, sizeFactors(ddstrna), "/")

#----- Write size factors (two-row: sample names, then values)
write.table(
  rbind(rownames(meta), sizeFactors(ddstrna)),
  file      = paste0(normalizedDir, "tRNA_isotype_counts_size_factors.csv"),
  row.names = FALSE,
  col.names = FALSE
)
message(paste0("\tSize factors saved to ", normalizedDir, "tRNA_isotype_counts_size_factors.csv"))

#----- Write normalized counts
write.csv(trna_normalizedCounts, file = paste0(normalizedDir, "normalized_tRNA_isotype_counts.csv"))
message(paste0("\tNormalized counts saved to ", normalizedDir, "normalized_tRNA_isotype_counts.csv"))

#----- Run DESeq2 (tRAX-style: betaPrior=TRUE)
ddstrna <- DESeq2::DESeq(ddstrna, betaPrior = TRUE)

saveRDS(ddstrna, file = paste0(rdsDir, "tRNA_isotype_DESeq2_object.Rds"))
message(paste0("\n\tSaved RDS to ", rdsDir, "tRNA_isotype_DESeq2_object.Rds"))

#----- Differential expression
runDE(ddstrna, trna_normalizedCounts, meta, comparisons, "tRNA_isotype", deDir)

#----- rlog for PCA
ddstrnaRlog    <- DESeq2::rlog(ddstrna)
ddstrnaRlogMat <- SummarizedExperiment::assay(ddstrnaRlog) %>% as.matrix()

variancetrna <- apply(ddstrnaRlogMat, 1, var)
trnaVars     <- sort(variancetrna, decreasing = TRUE)

variancetrna <- data.frame(Features = seq_along(trnaVars), Variance = trnaVars)

trnaVar <- ggplot2::ggplot(variancetrna, aes(x = Features, y = Variance)) +
  geom_point() +
  theme_classic() +
  labs(title = "tRNA Isodecoder Variance Plot",
       x     = "Number of Features (tRNAs)",
       y     = "Variance") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text  = element_text(size = 12))
ggplot2::ggsave(trnaVar, filename = paste0(pcaDir, "tRNA_isotype_variance_plot.png"), width = 6, height = 6)
message("\tPlotted tRNA_isotype_variance_plot.png\n")

message(paste0("Length Variance tRNA ", length(variancetrna)))

PCstrna      <- generatePCs(ddstrnaRlogMat, trnaVars, nrow(variancetrna))
loadingstrna <- PCstrna[[1]]
loadingstrna$Sample <- rownames(loadingstrna)
loadingstrna <- loadingstrna[match(rownames(meta), rownames(loadingstrna)), ]
loadingstrna$Group <- meta$Group

write.csv(loadingstrna, file = paste0(pcaDir, "tRNA_isotype_loadings.csv"))
message(paste0("PCA loadings saved to ", pcaDir, "tRNA_isotype_loadings.csv"))

PCAplottrna <- ggplot2::ggplot(loadingstrna, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 5) +
  geom_text_repel() +
  labs(title = "tRNA Isodecoder PCA",
       x     = paste0("PC1: ", PCstrna[[3]][1]),
       y     = paste0("PC2: ", PCstrna[[3]][2])) +
  theme_classic(base_size = 16) +
  theme(axis.text        = element_text(size = 14, face = "bold"),
        axis.title       = element_text(size = 16, face = "bold"),
        legend.position  = "right",
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA),
        title            = element_text(face = "bold"))
ggplot2::ggsave(PCAplottrna, file = paste0(pcaDir, "tRNA_isotype_PCA.png"), width = 6, height = 6)
message("\tPlotted tRNA_isotype_PCA.png")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# COMBINED PCA PLOT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
commonLegend <- suppressMessages(cowplot::get_legend(PCAplotFull))
PCAplotFull  <- PCAplotFull + theme(legend.position = "none")
PCAplottrna  <- PCAplottrna + theme(legend.position = "none")

combinedPCA <- cowplot::plot_grid(
  cowplot::plot_grid(PCAplotFull, PCAplottrna, ncol = 1, align = "v"),
  commonLegend,
  ncol       = 2,
  rel_widths = c(1, 0.3)
)
ggplot2::ggsave(paste0(pcaDir, "PCA_Analysis_Summary.png"), width = 8, height = 8)
message("\tPlotted PCA_Analysis_Summary.png")

loadingsFull$Analysis <- c("tRNAs + smRNAs")
loadingstrna$Analysis <- c("tRNA Isodecoders")
loadingsCombined <- rbind(loadingsFull, loadingstrna)

PCAplotBoth <- ggplot2::ggplot(
    loadingsCombined,
    aes(x = PC1, y = PC2, color = Group, label = Sample, shape = Analysis)
  ) +
  geom_point(size = 5) +
  geom_text_repel() +
  geom_line(aes(group = Sample), linetype = "dotted", color = "darkgray") +
  scale_shape_manual(values = c("tRNAs + smRNAs" = 16, "tRNA Isodecoders" = 2)) +
  labs(title = "",
       x     = paste0("PC1 tRNA: ", PCstrna[[3]][1], "\nPC1 tRNA + smRNA: ", PCsFull[[3]][1]),
       y     = paste0("PC2 tRNA: ", PCstrna[[3]][2], "\nPC2 tRNA + smRNA: ", PCsFull[[3]][2])) +
  theme_classic(base_size = 16) +
  theme(axis.text        = element_text(size = 14, face = "bold"),
        axis.title       = element_text(size = 16, face = "bold"),
        legend.position  = "right",
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA))
ggplot2::ggsave(PCAplotBoth, file = paste0(pcaDir, "PCA_Direct_Comparison.png"), width = 8, height = 8)
message("\tPlotted PCA_Direct_Comparison.png")
