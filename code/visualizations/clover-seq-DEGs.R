#~~~~~~~~~~~~~~~~~~~~~~~~ README ~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 
# Title: clover-seq-DEGs.R
# Description: Plotting functions for differential expression
#
# Author: Mike Martinez
# Lab: GDSC
# Project: Clover-Seq
# Date created: 05/28/25
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 

#----- Load librarires
suppressMessages(library(ggplot2))
suppressMessages(library(ggblend))
suppressMessages(library(ggrepel))

#----- Function to plot MA plots
MAplot <- function(results, title, figName) {
  #----- Set colors
  colors <- c("firebrick", 
              "steelblue", 
              "#b0b0b0") 

  #----- Filter to just rows with top 15% base mean values (for labelling)
  threshold <- quantile(results$baseMean, 0.85, na.rm = TRUE)
  filtered <- subset(results, baseMean > threshold) 

  #----- Get the top 5 downregulated and top 5 upregulated genes based on l2FC
  top_hitsB <- filtered[order(filtered$log2FoldChange, decreasing = FALSE), ][1:5, ]
  top_hitsA <- filtered[order(filtered$log2FoldChange, decreasing = TRUE), ][1:5, ] 

  #----- Combine the hits
  top_hits <- rbind(top_hitsB, top_hitsA)
  top_hits$feature <- rownames(top_hits)

  #----- Trace hits
  trace_filter <- abs(results$log2FoldChange) > log2FC_thresh & results$padj <= padj_thresh

  #----- Plot
  MAFigure <- ggplot(results, aes(x = log(baseMean), y = log2FoldChange, fill = Group, size = -log10(padj))) +
  geom_point(data = results, aes(fill = Group), 
             color = "white", fill = "grey85", shape = 21, alpha = 0.4, stroke = 0.3) +
  geom_point(data = results[trace_filter, ], 
             aes(fill = Group), color = "black", shape = 21, alpha = 0.75, stroke = 0.8) +
  geom_text_repel(data = top_hits, aes(label = feature), size = 4, max.overlaps = 110) +
  annotate("text", x = -Inf, y = 1.2, label = paste0("Enriched in ", numerator), hjust = -0.1, vjust = 0, fontface = 4) +
  annotate("text", x = -Inf, y = -1.2, label = paste0("Enriched in ", refLevel), hjust = -0.1, vjust = 1, fontface = 4) +
  geom_hline(yintercept = -1, linetype = "dashed") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_hue(direction = 1) +
  labs(
    y = "Log2 Fold Change",
    x = "Log10 Base Mean",
    size = "-log10 padj",
    color = "",
    fill = ""
  ) +
  scale_size(range = c(1, 10)) +
  theme_classic(base_size = 16) +
  theme(axis.title = element_text(face = "bold")) +
  guides(fill = "none")

  ggsave(paste0(figDir, figName), MAFigure, width = 8, height = 8)
  message(paste0("Plotted ", paste0(figDir, figName)))
}

#----- Function to plot volcano plots
VolcanoPlot <- function(results, title, figName) {
  #----- Set colors
  colors <- c("firebrick", 
              "steelblue", 
              "#b0b0b0") 

  #----- Filter to just rows with top 15% base mean values (for labelling)
  threshold <- quantile(results$baseMean, 0.85, na.rm = TRUE)
  filtered <- subset(results, baseMean > threshold) 

  #----- Get the top 5 downregulated and top 5 upregulated genes based on l2FC
  top_hitsB <- filtered[order(filtered$log2FoldChange, decreasing = FALSE), ][1:5, ]
  top_hitsA <- filtered[order(filtered$log2FoldChange, decreasing = TRUE), ][1:5, ] 

  #----- Combine the hits
  top_hits <- rbind(top_hitsB, top_hitsA)
  top_hits$feature <- rownames(top_hits)

  #----- Plot
  VolFigure <- ggplot(results, aes(x = log2FoldChange, y = -log10(padj), color = Group)) +
  geom_point(alpha = 0.5) +  
  geom_text_repel(data = top_hits, aes(label = feature), size = 4, max.overlaps = 100) +
  geom_vline(xintercept = denBound, linetype = "dashed") +
  geom_vline(xintercept = numBound, linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed") +
  labs(y = "-log10 padj",
       x = "Log2 Fold Change",
       title = title) +
  theme_minimal(base_size = 16) +
  theme(axis.title = element_text(face = "bold"),
        panel.background = element_rect(fill = "white", color = NA),  
        plot.background = element_rect(fill = "white", color = NA),
        title = element_text(face = "bold"))
  ggsave(paste0(figDir, figName), VolFigure, width = 8, height = 8)
  message(paste0("Plotted ", paste0(figDir, figName)))
}

#----- Function to plot heatmap


