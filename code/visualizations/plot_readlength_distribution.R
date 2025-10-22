#~~~~~~~~~~~~~~~~~~~~~~~~ README ~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 
# Title: plot_readlength_distribution.R
# Description: Visualize read length distributions per sample
#
# Author: Mike Martinez
# Lab: Genomic Data Science Core
# Project: CloverSeq-Pipeline
# Date created: 5/14/25
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# LOAD LIBRARIES AND SET PATHS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#----- Libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(magrittr)
library(tidyr)
library(RColorBrewer)
library(ggrepel)
library(cowplot)

#----- Set command line args
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: RScript plot_readlength_distribution.R <input read length dist> <output_dir>")
}

input_tsv <- args[1]
outputDir <- args[2]

if (!dir.exists(outputDir)) {
  dir.create(outputDir)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# READ IN THE READ LENGTHS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
dist <- read.csv(input_tsv, sep = "\t")

rld_plot <- ggplot2::ggplot(dist, aes(x = Length, y = Count, fill = Sample)) +
  geom_col(position = "dodge", width = 0.8, color = "black") +
  facet_wrap(~Sample, scales = "free_y") +
  scale_x_continuous(breaks = seq(0, 100, by = 10)) +
  theme_classic(base_size = 16) +
  labs(title = "Read Length Distribution by Sample",
       x = "Read Length",
       y = "Count") +
  theme(panel.spacing = unit(1, "lines"),
        legend.position = "none",
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"))
ggplot2::ggsave(rld_plot, file = paste0(outputDir, "read_length_distribution_by_sample.png"), width = 10, height = 10)
