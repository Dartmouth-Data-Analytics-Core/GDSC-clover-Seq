#~~~~~~~~~~~~~~~~~~~~~~~~ README ~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 
# Title: clover-seq-plot-coverages.R
# Description: Plot coverage plots for mature tRNA features
#
# Author: Mike Martinez
# Lab: GDSC
# Project: Clover-Seq
# Date created: 05/30/25
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
covFile <- args[1]

if (length(args) < 1 | length(args) > 1) {
  stop("Usage: RScript clover-seq-plot-coverages.R <coverage file>")
}

#----- Check input metadata exists
if (!file.exists(covFile)) {
  stop(paste0("Input coverage file: ", covFile, " does not exists!"))
}

#----- Set input directories
trnaDir <- "03_tRNA_counts/"

#----- Set output directory and create
opDir <- "08_plots/mature_tRNA_coverages/"
if (!dir.exists(opDir)) {
  dir.create(opDir)
}

#----- Check all input directories exist
message("--------------------------------------------------")
message(paste("Input directories:", trnaDir, sep = "\n\t"))
message("Checking that all input directories exist...")
dirsToCheck <- c(trnaDir)
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
# READ IN THE COVERAGE DATA AND FORMAT FOR PLOTTING
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
covs <- read_file_safe(paste0(trnaDir, "mature_tRNA_coverages.txt"), sep = "\t", row.names = NULL)

#----- Set explicit ordering for Sprinzl positions
positions <- c(
  "-1","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17",
  "17a","18","19","20","20a","20b","21","22","23","24","25","26","27","28",
  "29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45",
  'e1','e2','e3','e4','e5','e6','e7','e8','e9','e10','e11',
  'e12','e13','e14','e15','e16','e17','e18','e19',"46","47","48",
  "49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66",
  "67","68","69","70","71","72","73","74","75","76")

#----- Pivot to long format for plotting
covsLong <- tidyr::pivot_longer(data = covs,
                                -c(Sample, Feature),
                                names_to = "Position",
                                values_to = "normalized_coverage")

#----- Clean the positions and factor
covsLong$Position <- gsub("^X", "", covsLong$Position)
covsLong$Position <- factor(covsLong$Position, levels = positions)

#----- Get a vector of unique features to iterate over
maturetRNAs <- unique(covsLong$Feature)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ITERATE OVER THE DATA AND PLOT COVERAGE PLOTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
message(paste0("Plotting normalized coverages for...\n"))
for (i in maturetRNAs) {
    #----- Subset the data
    subset <- covsLong[covsLong$Feature == i,]

    #----- Debugging
    message(paste0("\t", i))

    #----- Set plot filename
    plotName <- paste0(i, "_normalized_coverage.png")
    fileName <- paste0(opDir, plotName)

    #----- Plot
    covPlot <- ggplot2::ggplot(subset, aes(x = Position, y = normalized_coverage, fill = Sample, group = Sample)) +
        geom_col(width = 1) +
        labs(title = i,
            y = "Normalized Coverage",
            x = "") +
        theme_classic() +
        facet_grid(Sample~.) +
        theme(
            axis.text.x = element_text(angle = 90, size = 9),
            legend.position = "none",
            axis.title.y = element_text(face = "bold", size = 16),
            strip.text = element_text(face = "bold", size = 16),
            title = element_text(face = 4, size = 16))
    ggplot2::ggsave(fileName, covPlot, width = 12, height = 10)
}

message("--------------------------------------------------")
message("Coverage plotting finished!")
message("--------------------------------------------------")
