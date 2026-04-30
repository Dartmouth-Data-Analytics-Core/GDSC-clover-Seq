#~~~~~~~~~~~~~~~~~~~~~~~~ README ~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Title: tRAX mismatch plots
# Description: Generates mismatch heatmaps and supporting
#              diagnostic figures from tRAX coverage output
#
# Author: tRAX / Mike Martinez
# Lab: Orellana
# Project: clover-Seq development
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# LOAD LIBRARIES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(scales)
library(plyr)
library(gridExtra)
library(getopt)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# PARSE ARGUMENTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
spec <- matrix(c(
  'mismatch'    , 'm', 1, "character", "mismatches file",
  'samples'     , 's', 1, "character", "sample CSV file (required)",
  'trna'        , 't', 1, "character", "trna table file (required)",
  'comparisons' , 'c', 1, "character", "comparisons file",
  'directory'   , 'd', 1, "character", "output directory (required)",
  'help'        , 'h', 0, "logical",   "this help"
), ncol = 5, byrow = TRUE)

opt <- getopt(spec)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SET DIRECTORIES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
directory  <- opt$directory
dumpDir    <- paste0(directory, "/Misc/")
outputFmt  <- ".pdf"

#----- Create output directories
if (!dir.exists(directory)) dir.create(directory, recursive = TRUE)
if (!dir.exists(dumpDir))   dir.create(dumpDir,   recursive = TRUE)

message("--------------------------------------------------")
message(paste0("Output directory:  ", directory))
message(paste0("Dump directory:    ", dumpDir))
message("--------------------------------------------------")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# READ IN DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
mismatches     <- read.table(opt$mismatch, header = TRUE, row.names = NULL, stringsAsFactors = FALSE)
trnatable      <- read.table(opt$trna)
Sampletable.raw <- read.csv(opt$samples, header = TRUE, stringsAsFactors = FALSE)
Sampletable    <- data.frame(
  key     = Sampletable.raw$Sample_ID,
  display = Sampletable.raw$Sample_ID,
  stringsAsFactors = FALSE
)

message(paste0(opt$mismatch, " loaded successfully"))
message(paste0(opt$samples,  " loaded successfully"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# CONSTANTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
mismatchPseudocount <- 10
dotSize             <- 0.4
aminoDotSize        <- 0.8

positionorder <- c(
  '-1','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16',
  '17','18','19','20','21','22','23','24','25','26','27','28','29','30','31',
  '32','33','34','35','36','37','38','39','40','41','42','43','44','45',
  'e1','e2','e3','e4','e5','e6','e7','e8','e9','e10','e11','e12','e13','e14',
  'e15','e16','e17','e18','e19',
  '46','47','48','49','50','51','52','53','54','55','56','57','58','59','60',
  '61','62','63','64','65','66','67','68','69','70','71','72','73','74','75','76'
)

getPalette <- colorRampPalette(brewer.pal(9, "Set1"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# COMPUTE MISMATCH AND DELETION PERCENTAGES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
mismatches$percentmismatch <- mismatches$mismatchedbases / (mismatches$coverage + mismatchPseudocount)
mismatches$percentdelete   <- mismatches$deletedbases    / (mismatches$coverage + mismatchPseudocount)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# IDENTIFY HIGH-MISMATCH POSITIONS (written to DUMP)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
totalmism <- aggregate(
  mismatches$percentmismatch,
  by  = list(position = mismatches$position, Feature = mismatches$Feature),
  FUN = max
)
posmism <- aggregate(
  totalmism$x,
  by  = list(position = totalmism$position),
  FUN = function(mism) { sum(mism > .1) }
)
colnames(posmism) <- c("tRNA_position", "Mismatched_Transcripts")
write.table(posmism, file = paste0(dumpDir, "positionmismatches.txt"), row.names = TRUE)
mismatchpositions <- posmism[posmism$Mismatched_Transcripts > 5, "tRNA_position"]

totaldelete <- aggregate(
  mismatches$deletions / (mismatches$deletions + mismatches$coverage + 30),
  by  = list(position = mismatches$position, Feature = mismatches$Feature),
  FUN = max
)
posdelete <- aggregate(
  totaldelete$x,
  by  = list(position = totaldelete$position),
  FUN = function(mism) { sum(mism > .1) }
)
colnames(posdelete) <- c("tRNA_position", "Deleted_Transcripts")
write.table(posdelete, file = paste0(dumpDir, "positiondeletions.txt"), row.names = TRUE)
deletepositions <- posdelete[posdelete$Deleted_Transcripts > 5, "tRNA_position"]

positions <- union(mismatchpositions, deletepositions)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# GLOBAL POSITION MISMATCH DOTPLOT (DUMP)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
mismatchmeltfilter <- mismatches[
  mismatches$adenines + mismatches$thymines + mismatches$cytosines + mismatches$guanines > 50,
  c("Feature", "Sample", "percentmismatch", "position")
]

mismatchmeltposagg <- aggregate(
  mismatchmeltfilter$percentmismatch,
  by  = list(position = mismatchmeltfilter$position, Feature = mismatchmeltfilter$Feature),
  FUN = max
)
mismatchmeltposagg <- mismatchmeltposagg[mismatchmeltposagg$position %in% positionorder, ]
mismatchmeltposagg$position <- factor(mismatchmeltposagg$position, levels = positionorder)
colnames(mismatchmeltposagg) <- c("Position", "Feature", "percentmismatch")
mismatchmeltposagg$amino    <- trnatable[match(mismatchmeltposagg$Feature, trnatable[, 1]), 3]
mismatchmeltposagg$anticodon <- trnatable[match(mismatchmeltposagg$Feature, trnatable[, 1]), 4]

posname <- paste0(dumpDir, "trnapositionmismatches", outputFmt)
ggplot(data = mismatchmeltposagg, aes(x = Position, y = percentmismatch)) +
  theme_bw() +
  geom_jitter(aes(color = amino), size = 1, width = 0.25) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  ggtitle("Position Mismatches") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Position") +
  ylab("Maximum percent Misincorporation") +
  geom_abline(intercept = .1, slope = 0, linetype = 2)
ggsave(filename = posname, width = 20, height = 7)
message(paste0("\tDUMP: ", basename(posname)))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# GLOBAL FIVE-PRIME END DOTPLOT (DUMP)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
fiveprimemeltfilter <- mismatches[mismatches$tRNAreadstotal > 30,
  c("Feature", "Sample", "position", "readstarts", "tRNAreadstotal", "coverage")]
fiveprimemeltfilter$percentstart <- fiveprimemeltfilter$readstarts / fiveprimemeltfilter$tRNAreadstotal

threeprimemeltfilter <- mismatches[mismatches$tRNAreadstotal > 30,
  c("Feature", "Sample", "position", "readends", "tRNAreadstotal", "coverage")]
threeprimemeltfilter$percentstart <- threeprimemeltfilter$readends / threeprimemeltfilter$tRNAreadstotal

mismatchfilter <- mismatches[mismatches$tRNAreadstotal > 30,
  c("Feature", "Sample", "position", "percentmismatch", "coverage")]
deletefilter <- mismatches[mismatches$tRNAreadstotal > 30,
  c("Feature", "Sample", "position", "percentdelete", "coverage")]

fiveprimeposagg <- aggregate(
  fiveprimemeltfilter$percentstart,
  by  = list(position = fiveprimemeltfilter$position, Sample = fiveprimemeltfilter$Sample),
  FUN = mean
)
fiveprimeposagg <- fiveprimeposagg[fiveprimeposagg$position %in% positionorder, ]
fiveprimeposagg$position <- factor(fiveprimeposagg$position, levels = positionorder)
colnames(fiveprimeposagg) <- c("Position", "Sample", "percentstart")

fiveprimemeltfilter$amino  <- trnatable[match(fiveprimemeltfilter$Feature, trnatable[, 1]), 3]
threeprimemeltfilter$amino <- trnatable[match(threeprimemeltfilter$Feature, trnatable[, 1]), 3]
mismatchfilter$amino       <- trnatable[match(mismatchfilter$Feature, trnatable[, 1]), 3]
deletefilter$amino         <- trnatable[match(deletefilter$Feature, trnatable[, 1]), 3]

posname <- paste0(dumpDir, "trnapositionfiveprime", outputFmt)
ggplot(data = fiveprimeposagg, aes(x = Position, y = percentstart)) +
  theme_bw() +
  geom_jitter(aes(color = Sample), size = dotSize, width = 0.25) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  ggtitle("Position Starts") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Position") +
  ylab("Maximum Starts")
ggsave(filename = posname, width = 20, height = 7)
message(paste0("\tDUMP: ", basename(posname)))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# PER-ISOTYPE LOOP
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
for (curramino in unique(trnatable[, 3])) {

  #----- Subset data for this isotype
  fiveprimeamino  <- fiveprimemeltfilter[fiveprimemeltfilter$amino == curramino, ]
  threeprimeamino <- threeprimemeltfilter[threeprimemeltfilter$amino == curramino, ]
  mismatchamino   <- mismatchfilter[mismatchfilter$amino == curramino, ]
  deleteamino     <- deletefilter[deletefilter$amino == curramino, ]

  if (nrow(fiveprimeamino) < 1) next

  #----- Per-isotype position mismatch dotplot (DUMP)
  posname <- paste0(dumpDir, curramino,"_trnapositionmismatches", outputFmt)
  mismatchmeltposaggamino <- mismatchmeltposagg[mismatchmeltposagg$amino == curramino, ]
  ggplot(data = mismatchmeltposaggamino, aes(x = Position, y = percentmismatch)) +
    theme_bw() +
    geom_jitter(aes(color = anticodon), size = aminoDotSize, width = 0.25) +
    scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
    ggtitle("Position Mismatches") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab("Position") +
    ylab("Maximum percent Misincorporation")
  ggsave(filename = posname, width = 20, height = 7)
  message(paste0("\tDUMP: ", basename(posname)))

  #----- Five-prime end bar chart (DUMP)
  fiveprimeaminoposagg <- aggregate(
    fiveprimeamino$percentstart,
    by  = list(
      Sample   = Sampletable[match(fiveprimeamino$Sample, Sampletable[, 1]), 2],
      Position = fiveprimeamino$position
    ),
    FUN = mean
  )
  fiveprimeaminoposagg <- fiveprimeaminoposagg[fiveprimeaminoposagg$Position %in% positionorder, ]
  fiveprimeaminoposagg$Position <- factor(fiveprimeaminoposagg$Position, levels = positionorder)
  colourCount <- length(unique(fiveprimeaminoposagg$Position)) + 1

  currplot <- ggplot(fiveprimeaminoposagg, aes(x = Sample, y = x, fill = Position, stat = "identity")) +
    theme_bw() +
    theme(panel.border = element_rect(linetype = "blank"), panel.grid = element_line(linetype = "blank")) +
    geom_bar(position = "fill", stat = "identity") +
    geom_bar(position = "fill", stat = "identity", color = "black", show.legend = FALSE) +
    scale_y_continuous(labels = percent_format()) +
    theme(axis.text.x = element_text(size = 5)) +
    ggtitle(paste0(curramino, " Five-prime ends")) +
    xlab("Sample") +
    ylab("Percentage of reads that start at position") +
    labs(fill = "position") +
    scale_fill_manual(values = getPalette(colourCount)) +
    theme(axis.title.x = element_text(face = "bold", size = 15),
          axis.text.x  = element_text(face = "bold", size = 9, angle = 90, vjust = .5))

  posname <- paste0(dumpDir, curramino,"_fiveprimecounts", outputFmt)
  ggsave(posname, currplot)
  message(paste0("\tDUMP: ", basename(posname)))

  #----- Three-prime end bar chart (DUMP)
  threeprimeaminoposagg <- aggregate(
    threeprimeamino$percentstart,
    by  = list(
      Sample   = Sampletable[match(threeprimeamino$Sample, Sampletable[, 1]), 2],
      Position = threeprimeamino$position
    ),
    FUN = mean
  )
  threeprimeaminoposagg <- threeprimeaminoposagg[threeprimeaminoposagg$Position %in% positionorder, ]
  threeprimeaminoposagg$Position <- factor(threeprimeaminoposagg$Position, levels = positionorder)
  colourCount <- length(unique(threeprimeaminoposagg$Position)) + 1

  currplot <- ggplot(threeprimeaminoposagg, aes(x = Sample, y = x, fill = Position, stat = "identity")) +
    theme_bw() +
    theme(panel.border = element_rect(linetype = "blank"), panel.grid = element_line(linetype = "blank")) +
    geom_bar(position = "fill", stat = "identity") +
    geom_bar(position = "fill", stat = "identity", color = "black", show.legend = FALSE) +
    scale_y_continuous(labels = percent_format()) +
    theme(axis.text.x = element_text(size = 5)) +
    ggtitle(paste0(curramino, " Three-prime ends")) +
    xlab("Sample") +
    ylab("Percentage of reads that end at position") +
    labs(fill = "position") +
    scale_fill_manual(values = getPalette(colourCount)) +
    theme(axis.title.x = element_text(face = "bold", size = 15),
          axis.text.x  = element_text(face = "bold", size = 9, angle = 90, vjust = .5))

  posname <- paste0(dumpDir, curramino,"_threeprimecounts", outputFmt)
  ggsave(posname, currplot)
  message(paste0("\tDUMP: ", basename(posname)))

  #----- Five-prime heatmap (DUMP)
  fiveprimeaminoposagg <- aggregate(
    fiveprimeamino$percentstart,
    by  = list(
      Sample   = Sampletable[match(fiveprimeamino$Sample, Sampletable[, 1]), 2],
      position = fiveprimeamino$position,
      Feature  = fiveprimeamino$Feature
    ),
    FUN = mean
  )
  fiveprimeaminoposagg$percentstart <- fiveprimeaminoposagg$x
  fiveprimeaminoposagg <- fiveprimeaminoposagg[fiveprimeaminoposagg$position %in% positionorder, ]
  fiveprimeaminoposagg$position <- factor(fiveprimeaminoposagg$position, levels = positionorder)

  currplot <- ggplot(fiveprimeaminoposagg, aes(x = position, y = Sample, fill = percentstart)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "blue", limits = c(0, 1)) +
    theme_bw() +
    facet_grid(rows = vars(Feature)) +
    theme(panel.border   = element_rect(linetype = "blank"),
          panel.grid     = element_line(linetype = "blank"),
          axis.title.x   = element_blank(),
          axis.text.y    = element_text(colour = "black", size = 20),
          axis.text.x    = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16),
          strip.text.y   = element_text(size = 24, angle = 0)) +
    xlab("Position") +
    ylab("Sample") +
    labs(fill = "five-prime ends")

  posname <- paste0(dumpDir, curramino,"_fiveprimeheatmap", outputFmt)
  ggsave(posname, currplot,
    height = .25 * length(unique(fiveprimeamino$Feature)) * length(unique(fiveprimeamino$Sample)),
    width  = .25 * length(unique(fiveprimeamino$position)),
    limitsize = FALSE
  )
  message(paste0("\tDUMP: ", basename(posname)))

  #----- Three-prime heatmap (DUMP)
  threeprimeaminoposagg <- aggregate(
    threeprimeamino$percentstart,
    by  = list(
      Sample   = Sampletable[match(threeprimeamino$Sample, Sampletable[, 1]), 2],
      position = threeprimeamino$position,
      Feature  = threeprimeamino$Feature
    ),
    FUN = mean
  )
  threeprimeaminoposagg$percentstart <- threeprimeaminoposagg$x
  threeprimeaminoposagg <- threeprimeaminoposagg[threeprimeaminoposagg$position %in% positionorder, ]
  threeprimeaminoposagg$position <- factor(threeprimeaminoposagg$position, levels = positionorder)

  currplot <- ggplot(threeprimeaminoposagg, aes(x = position, y = Sample, fill = percentstart)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "green", limits = c(0, 1)) +
    theme_bw() +
    facet_grid(rows = vars(Feature)) +
    theme(panel.border   = element_rect(linetype = "blank"),
          panel.grid     = element_line(linetype = "blank"),
          axis.title.x   = element_blank(),
          axis.text.y    = element_text(colour = "black", size = 20),
          axis.text.x    = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16),
          strip.text.y   = element_text(size = 24, angle = 0)) +
    xlab("Position") +
    ylab("Sample") +
    labs(fill = "three-prime ends")

  posname <- paste0(dumpDir, curramino,"_threeprimeheatmap", outputFmt)
  ggsave(posname, currplot,
    height = .25 * length(unique(fiveprimeamino$Feature)) * length(unique(fiveprimeamino$Sample)),
    width  = .25 * length(unique(fiveprimeamino$position)),
    limitsize = FALSE
  )
  message(paste0("\tDUMP: ", basename(posname)))

  #----- Mismatch heatmap (KEEP)
  mismatchaminoposagg <- aggregate(
    mismatchamino$percentmismatch,
    by  = list(
      Sample   = Sampletable[match(mismatchamino$Sample, Sampletable[, 1]), 2],
      position = mismatchamino$position,
      Feature  = mismatchamino$Feature
    ),
    FUN = mean
  )
  mismatchaminoposagg$percentmismatch <- mismatchaminoposagg$x
  mismatchaminoposagg <- mismatchaminoposagg[mismatchaminoposagg$position %in% positionorder, ]
  mismatchaminoposagg$position <- factor(mismatchaminoposagg$position, levels = positionorder)

  currplot <- ggplot(mismatchaminoposagg, aes(x = position, y = Sample, fill = percentmismatch)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "red", limits = c(0, 1)) +
    theme_bw() +
    facet_grid(rows = vars(Feature)) +
    theme(panel.border   = element_rect(linetype = "blank"),
          panel.grid     = element_line(linetype = "blank"),
          axis.title.x   = element_blank(),
          axis.text.y    = element_text(colour = "black", size = 20),
          axis.text.x    = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16),
          strip.text.y   = element_text(size = 24, angle = 0)) +
    xlab("Position") +
    ylab("Sample") +
    labs(fill = "mismatches")

  posname <- paste0(directory, "/", curramino, "_mismatchheatmap", outputFmt)
  ggsave(posname, currplot,
    height = .25 * length(unique(fiveprimeamino$Feature)) * length(unique(fiveprimeamino$Sample)),
    width  = .25 * length(unique(fiveprimeamino$position)),
    limitsize = FALSE
  )
  message(paste0("\tPlotted: ", basename(posname)))

  #----- Deletion heatmap (DUMP)
  deleteaminoposagg <- aggregate(
    deleteamino$percentdelete,
    by  = list(
      Sample   = Sampletable[match(deleteamino$Sample, Sampletable[, 1]), 2],
      position = deleteamino$position,
      Feature  = deleteamino$Feature
    ),
    FUN = mean
  )
  deleteaminoposagg$percentdelete <- deleteaminoposagg$x
  deleteaminoposagg <- deleteaminoposagg[deleteaminoposagg$position %in% positionorder, ]
  deleteaminoposagg$position <- factor(deleteaminoposagg$position, levels = positionorder)

  currplot <- ggplot(deleteaminoposagg, aes(x = position, y = Sample, fill = percentdelete)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "purple", limits = c(0, 1)) +
    theme_bw() +
    facet_grid(rows = vars(Feature)) +
    theme(panel.border   = element_rect(linetype = "blank"),
          panel.grid     = element_line(linetype = "blank"),
          axis.title.x   = element_blank(),
          axis.text.y    = element_text(colour = "black", size = 20),
          axis.text.x    = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16),
          strip.text.y   = element_text(size = 24, angle = 0)) +
    xlab("Position") +
    ylab("Sample") +
    labs(fill = "deletions")

  posname <- paste0(dumpDir, curramino,"_deletionheatmap", outputFmt)
  ggsave(posname, currplot,
    height = .25 * length(unique(fiveprimeamino$Feature)) * length(unique(fiveprimeamino$Sample)),
    width  = .25 * length(unique(fiveprimeamino$position)),
    limitsize = FALSE
  )
  message(paste0("\tDUMP: ", basename(posname)))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# PER-POSITION LOOP (ALL TO DUMP)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
for (currpos in positionorder) {

  poslabel <- ifelse(currpos == "-1", "neg1", currpos)

  mismatchmelt  <- mismatches[mismatches$position == currpos, c("Feature", "Sample", "percentmismatch")]
  deletionmelt  <- mismatches[mismatches$position == currpos, c("Feature", "Sample", "percentdelete")]
  fiveprimemelt <- mismatches[mismatches$position == currpos,
    c("Feature", "Sample", "readstarts", "tRNAreadstotal", "coverage")]
  fiveprimemelt$percentstart <- fiveprimemelt$readstarts / fiveprimemelt$tRNAreadstotal

  if (nrow(fiveprimemelt) == 0) next

  #----- Per-position read starts boxplot (DUMP)
  fiveprimeagg <- aggregate(
    fiveprimemelt$percentstart,
    by  = list(Feature = fiveprimemelt$Feature,
               Sample  = Sampletable[match(fiveprimemelt$Sample, Sampletable[, 1]), 2]),
    FUN = mean
  )
  colnames(fiveprimeagg) <- c("Feature", "Sample", "percentstart")

  posname <- paste0(dumpDir, poslabel, "_possamplereadstarts", outputFmt)
  ggplot(data = fiveprimeagg, aes(x = Sample, y = percentstart)) +
    geom_boxplot(aes(fill = Sample), outlier.shape = NA) +
    theme_bw() +
    geom_jitter(aes(fill = Sample), size = dotSize) +
    ggtitle(paste0("Position ", currpos, " Read Starts")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ylim(0, 1) +
    xlab("Sample") +
    ylab("Percent Read Starts")
  ggsave(filename = posname, width = 7, height = 7)

  #----- Per-position sample mismatch boxplot (DUMP)
  mismatchmeltagg <- aggregate(
    mismatchmelt$percentmismatch,
    by  = list(Feature = mismatchmelt$Feature,
               Sample  = Sampletable[match(mismatchmelt$Sample, Sampletable[, 1]), 2]),
    FUN = mean
  )
  colnames(mismatchmeltagg) <- c("Feature", "Sample", "percentmismatch")

  posname <- paste0(dumpDir, poslabel, "_possamplemismatches", outputFmt)
  ggplot(data = mismatchmeltagg, aes(x = Sample, y = percentmismatch)) +
    geom_boxplot(aes(fill = Sample), outlier.shape = NA) +
    theme_bw() +
    geom_jitter(aes(fill = Sample), size = dotSize) +
    ggtitle(paste0("Position ", currpos, " Mismatches")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
    xlab("Sample") +
    ylab("Percent Misincorporation")
  ggsave(filename = posname, width = 7, height = 7)

  #----- Per-position isotype mismatch boxplot (DUMP)
  mismatchmeltaminoagg <- aggregate(
    mismatchmelt$percentmismatch,
    by  = list(Sample = mismatchmelt$Sample,
               amino  = trnatable[match(mismatchmelt$Feature, trnatable[, 1]), 3]),
    FUN = mean
  )
  colnames(mismatchmeltaminoagg) <- c("Sample", "amino", "percentmismatch")

  posname <- paste0(dumpDir, poslabel, "_posaminomismatches", outputFmt)
  ggplot(data = mismatchmeltaminoagg, aes(x = amino, y = percentmismatch)) +
    geom_boxplot(aes(fill = amino), outlier.shape = NA) +
    theme_bw() +
    geom_jitter(aes(fill = amino), size = dotSize) +
    ggtitle(paste0("Position ", currpos, " Mismatches")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
    xlab("Isotype") +
    ylab("Percent Misincorporation")
  ggsave(filename = posname, width = 7, height = 7)

  #----- Per-position isotype deletion boxplot (DUMP)
  deletemeltaminoagg <- aggregate(
    deletionmelt$percentdelete,
    by  = list(Sample = deletionmelt$Sample,
               amino  = trnatable[match(deletionmelt$Feature, trnatable[, 1]), 3]),
    FUN = mean
  )
  colnames(deletemeltaminoagg) <- c("Feature", "amino", "percentdeletions")

  posname <- paste0(dumpDir, poslabel, "_posaminodeletions", outputFmt)
  ggplot(data = deletemeltaminoagg, aes(x = amino, y = percentdeletions)) +
    geom_boxplot(aes(fill = amino), outlier.shape = NA) +
    theme_bw() +
    geom_jitter(aes(fill = amino), size = dotSize) +
    ggtitle(paste0("Position ", currpos, " Deletions")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
    xlab("Isotype") +
    ylab("Percent Skipped")
  ggsave(filename = posname, width = 7, height = 7)

  #----- Per-position sample deletion boxplot (DUMP)
  deletemeltSampleagg <- aggregate(
    deletionmelt$percentdelete,
    by  = list(Feature = deletionmelt$Feature,
               amino   = trnatable[match(deletionmelt$Feature, trnatable[, 1]), 3],
               Sample  = Sampletable[match(deletionmelt$Sample, Sampletable[, 1]), 2]),
    FUN = mean
  )
  colnames(deletemeltSampleagg) <- c("Feature", "amino", "Sample", "percentdeletions")

  posname <- paste0(dumpDir, poslabel, "_possampledeletions", outputFmt)
  ggplot(data = deletemeltSampleagg, aes(x = Sample, y = percentdeletions)) +
    geom_boxplot(aes(fill = Sample), outlier.shape = NA) +
    theme_bw() +
    geom_jitter(aes(fill = Sample), size = dotSize) +
    ggtitle(paste0("Position ", currpos, " Deletions")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
    xlab("Sample") +
    ylab("Percent Skipped")
  ggsave(filename = posname, width = 7, height = 7)
}

message("--------------------------------------------------")
message(paste0("Mismatch heatmaps written to: ", directory))
message(paste0("All other figures written to: ", dumpDir))
message("Done!")
