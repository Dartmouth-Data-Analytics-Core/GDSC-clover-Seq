# Clover-Seq Pipeline Documentation
![status](https://img.shields.io/badge/status-in--development-orange)
![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)
![Snakemake](https://img.shields.io/badge/Snakemake-v7.18.0-red?logo=snakemake&logoColor=white)
![Python Version](https://img.shields.io/badge/python-3.8.5-blue)
![R Version](https://img.shields.io/badge/R-4.4.3-blue)

Modular Snakemake workflows for the comprehensive analyses of mature tRNAs and other small RNAs (smRNAs) from high-throughput sequencing data. 


<img src="img/CloverSeq_CQB_logo.png" alt="Description" width="700" height="450" style="border: none;" />


# Table of Contents
- [Introduction](#introduction)
- [Installation](#installation)
- [Database Build Module](#database-build-module)
- [Database Build Implementation](#database-build-implementation)
- [Preprocessing Module](#preprocessing-module)
- [Preprocessing Implementation](#preprocessing-implementation)
- [Preprocessing Outputs](#preprocessing-outputs)
    - [01 Trimming](#01-trimming)
    - [02 tRNA Alignment and Alignment Stats](#02-trna-alignment-and-alignment-stats)
    - [03 tRNA Counts](#03-trna-counts)
    - [04 smRNA Counts](#04-smrna-counts)
    - [05 Normalized](#05-normalized)
    - [06 PCA](#06-pca)
    - [07 Rds Files](#07-rds-files)
    - [08 Plots](#08-plots)
    - [09 QC](#09-qc)
- [Differential Expression Implementation](#differential-expression-implementation)
    - [10 Differential Expression Outputs](#10-differential-expression-outputs)
    - [11 Clover Seq Report](#11-clover-seq-report)
- [Development Notes](#development-notes)
- [Contact](#contact)
- [Citation and Licensing](#citation-and-licensing)


## Introduction
This pipeline supports the analysis of mature-tRNAs and other small RNAs (smRNAs) for human (hg38), mouse (mm10), and fly(dm6) genomes through 3 main modules. While typical RNA-Seq preprocessing strategies are emplyed in this pipeline, special considerations to handle tRNA biology are included. Custom reference databases encompass putative mature tRNAs as well as unique tRNA isodecoders as well as tRNA loci in the native host genome along with other annotated smRNAs including miRNA, sRNA, siRNA, snRNA, snoRNA, ribozymes and others. Mature tRNA transcripts are modified to include the addition of a 3' CCA tail which are not encoded genomically. Results are reported at two levels: the gene level which included mature tRNAs plus all other native tRNA-loci and smRNAs, and the isodecoder level, which is tRNA-specific.

The pipeline generates a number of helpful figures for analysis including PCA plots, coverage plots, and normalized abundances of mature tRNAs and other smRNAs, viewable as a pdf report generated at the end of the differential expression module.

<img src="img/pipeline_graphic.png" alt="Description" width="750" height="400" style="border: none;" />




## Installation
To install this code, clone the github repository

```shell
#----- Clone repository
git clone https://github.com/mikemartinez99/clover-Seq

```

Several [conda environments](https://anaconda.org/anaconda/conda) are required to run this code successfully. For your convenience, these conda environments have been prebuilt and are hosted publically at the following path:

`/dartfs/rc/nosnapshots/G/GMBSR_refs/envs/GDSC-Clover-Seq`

If you wish to build these environments yourself, the associated yaml file in the [env_config folder](https://github.com/mikemartinez99/clover-Seq/tree/main/workflows/env_config) can be built using the following command:

```shell
#----- Build conda environment
conda env create -f env_config/<name.yaml>
```


## Database Build Module
**NOTE: THIS IS OPTIONAL, BUT IT IS RECOMMENDED TO USE A PRE-BUILT DATABASE**

Reference databases are created through files obtained from [gtRNAdb](https://gtrnadb.ucsc.edu) through [tRNA-scan](https://lowelab.ucsc.edu/tRNAscan-SE/) experiments. These custom databases encompass not only mature-tRNA sequences, but the full native host genome as well, allowing access to all tNRA-loci. Mature tRNA sequences are created with the addition of a 3' CCA tail, removal of introns, and addition of the histidine post-transcriptional 5'G base. Sequences are padded on both ends by 20 bp to allow for extra bases when mapping. Within this module, multiple-sequence alignments (MSAs) are generated for both the mature tRNA sequences and genomic tRNA loci using [cmalign](https://manpages.ubuntu.com/manpages/trusty/man1/cmalign.1.html) and [Infernal](https://github.com/brendanf/inferrnal) with eukaryotic tRNA covariance models from tRNAscan_SE. These MSAs allow the downstream calculation of [Sprinzl-positions](http://polacek.dcbp.unibe.ch/publications/Holmes%20et%20al_tDR%20nomenclature_Nat.Meth_2023.pdf) for each tRNA, facilitating per-base read coverage. Additionally, this module generates a [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) index for the tRNA-genome (mature tRNA sequences + the full reference genome). 

These references have been pre-downloaded, pre-built, and pre-indexed and are hosted by [Genomics and Molecular Biology Shared Resources](https://geiselmed.dartmouth.edu/gsr/) on Discovery for ease of use and efficiency. However, if you wish to build from scratch, or customize the reference, the first module of this pipeline can be used. Pre-built references can be accessed at the following path on Discovery:

`/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/genomic_references/tRAX_databases`

<img src="img/database_rulegraph.png" alt="Description" width="200" height="300"/>

|Snakemake Rule|Purpose|Conda Environent|
|--------------|-------|----------------|
|`generate_gtRNA_db`|Download Ensembl GTF, gtRNA-db data tarball, genome and tRNA fasta file and build resulting tRNA database for downstream analysis in the preprocessing module|`clover-seq`|
|`concat_tRNAs`|Concatenate mature tRNA fasta and genome fasta to get tRNA-genome fasta file for use in profiling all smRNAs|`clover-seq`|
|`tRNA_bt2_index`|Create Bowtie2 index of tRNA-genome for use in alignment step of preprocessing module|`clover-bowtie2`|


## Database Build Implementation
To create a tRNA database from scratch, following the steps below:

1. Update your `module1_job.script.sh` script to point Snakemake to one of the prebuilt configs in `prebuilt_configs` folder. These configs are organism-specific and point to specific URLs on [gtRNAdb](https://gtrnadb.ucsc.edu). 

```shell
#----- Snakemake call
snakemake -s workflows/module-1-build-db.smk \
    --configfile prebuilt_configs/hg38_config.yaml \     # EDIT THIS LINE
    --use-conda \
    --conda-frontend conda \
    --conda-prefix /dartfs/rc/nosnapshots/G/GMBSR_refs/envs/GDSC-Clover-Seq \
    --profile cluster_profile \
    --rerun-incomplete \
    --keep-going 

```

2. Submit your `module1_job.script` script

```shell
#----- Submit snakemake job
sbatch module1_job.script.sh

```

To check the status of your Snakemake job and all child jobs it spawns, run the following (replacing NETID with your Dartmouth NetID)

**TIP** This job checking line can be applied for all 3 modules!

```shell
#----- Check user job status
squeue | grep "NETID"
```

**TIP** You will notice your directory will populate with logs following the convention: log_X_RuleID_JobID.out. For rules that run on a per-sample basis, there should be one log for each sample. For rules that run once for all samples together, there should be only a single log. 

## Preprocessing Module
Clover-Seq is adapted from the [tRAX Pipeline](https://github.com/UCSC-LoweLab/tRAX) to be implemented as a Snakemake workflow, allowing sample parallelization and improved modularity and efficiency. As input, this pipeline takes raw fastq.gz files, a config.yaml, and a sample metadata sheet.

This module starts with basic "best-practices" for RNA-Sequencing by trimming adapters. Currently, the Snakemake workflow only supports single-end, UMI-free reads, but support for paired-end and UMIs is under development. A minimum required read length of 15 bp is required for a read to be retained. Following adapter trimming, reads are mapped using [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) in very-sensitive mode, ignoring quality scores and allowing a maximum of 100 alignments per read. Non-tRNA smRNA annotations are obtained through an Ensembl gtf file (contained in the tRNA database). Read count metrics are calculated for both tRNA-specific alignments (mapping to the mature tRNA sequences) as well as tRNA-isotype mappings (mapping to tRNA-loci).

This module supports the calculation of tRNA fragment types as well. Reads where both the 5' and 3' end lie within 5 nt of their respective ends on the mature tRNAs are categorized as whole-tRNAs. Reads overlapping or aligning closely to the 5' end are categorized as 5' fragments. Likewise, reads aligning closely to the 3' end are categorized as 3' fragments. Fragments not meeting any of these criteria are categorized as "other" (Holmes, et al., 2022)

A gene-level dataset and tRNA-isotype dataset are generated. The former contains mature-tRNA mappings + all other smRNA featrures, and the latter contains tRNA-isotypes from genomic loci. These datasets are used to calculate principal components (separately), normalization size factors through [DESeq2](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html), and Rds files for differental expression in module 3. 

<img src="img/preprocessing_rulegraph.png" alt="Description" width="400" height="400"/>

|Snakemake Rule|Purpose|Conda Environent|
|--------------|-------|----------------|
|`trimming`|Trim fastq files to remove adapters and reads not meeting size threshold|`clover-seq`|
|`tRNA_align`|Align sequencing reads to combined tRNA sequences and host genome, perform length and quality filtering|`clover-bowtie2`|
|`tRNA_mark_duplicates`|Flag sequencing duplicates|`Picard`|
|`tRNA_map_stats`|Collate idxstats and flagstats metrics|`clover-seq`|
|`tRNA_count`|Count tRNA isotype reads as well as gene-level tRNAs/smRNAs in the genome|`clover-seq`|
|`read_length_distribution`|Collate read-length statistics for all reads (tRNA + non-tRNA)|`clover-seq`|
|`count_smRNAs`|Collate smRNA biotype counts across groups and samples|`clover-seq`|
|`normalize_and_PCA`|Calculate size factors and normalized-counts for tRNA-isotypes and gene-level counts, run PCA|`clover-seq`|
|`get_tRNA_coverage`|Calculare normalized coverage at each Sprinzl position for mature tRNAs|`clover-seq`|
|`get_mismatches`|Calculate normalized mismatches for each position in mature tRNAs|`clover-seq`|
|`plot_counts`|Generation of various plots|`clover-seq`|

## Preprocessing Implementation
To run the preprocessing pipeline, modify the `Sample_list_SE.txt` file. This is a comma-delimited three column file (Sample_ID, fastq_1, Group).
Note that modifying header names will cause the pipeline to fail! Simply modify your sample names, path to the raw fastq.gz file, and Group. Below is an example of how this file should look. 

```text
Sample_ID,fastq_1,Group
IB1,data/IB_1_S1_R1_001.fastq.gz,WT
IB2,data/IB_2_S2_R1_001.fastq.gz,WT
IB3,data/IB_3_S3_R1_001.fastq.gz,KO
IB4,data/IB_4_S4_R1_001.fastq.gz,KO

```

The `config.yaml` file provides important file paths that allow the Snakefile to dynamically run. Utilize the config specific to your organism in the [preprocessing prebuilt configs folder](https://github.com/mikemartinez99/clover-Seq/tree/main/preprocessing_prebuilt_configs). A list of case-sensitive parameters as specified in the config are listed below:

|Parameter|Value|
|---------|-----|
|`sample_txt`|Sample_list_SE.txt|
|`layout`|one of `single` or `paired` (currently only supports single)|
|`genome`|one of `hg38`, `mm10`, or `dm6`|
|`refLevel`|A reference level from the group column in Sample_list_SE.txt|
|`trna_db`|Path to publically hosted, pre-built tRNA databases|
|`bt2_index`|Path + bowtie2 index prefix for publically hosted, pre-built tRNA genome Bowtie2 index|
|`adapter_1`|Adapter sequence 1 (Illumina universal)|
|`adapter_2` |Adapter sequence 2 if paired-end (Illumina universal)|
|`maxMaps`|Number of mutli-mappings to accept per read|
|`nPenalty`|Penalty score for ambiguous bases (set to 5 to account for tRNA biology)|
|`log2FC_magnitude_threshold`|Magnitude of effect size (log2 Fold change) at which results will be considered significant|
|`padj_threshold`|Adjusted p-value threshold at which results will be considered statisitcally significant|


Modify the `module2_job.script.sh` script accordingly to point to your config file using the `--configfile` argument. Below is an example using the hg38 tRNA genome and specifying to submit 10 jobs in parallel across different cluster nodes as specified by `--profile cluster_profile`

```shell
#----- Run snakemake workflow
snakemake -s workflows/module-2-preprocess.smk \
    --use-conda \
    --configfile prebuilt_configs/hg38_config.yaml \     # EDIT THIS LINE
    --conda-frontend conda \
    --conda-prefix /dartfs/rc/nosnapshots/G/GMBSR_refs/envs/GDSC-Clover-Seq \
    --profile cluster_profile \
    --rerun-incomplete \
    --keep-going 
```

To submit your job, run the following code. 

```shell
#----- Submit snakemake job script
sbatch module2_job.script.sh
```

## Preprocessing Outputs
The preprocessing module of this workflow contains 9 major output folder. Their contents are exaplined in the following sections. 

### 01 Trimming
Contains trimmed fastq files and trimming logs in the `logs` folder (if SE and using CutAdapt)

### 02 tRNA Alignment and Alignment Stats
Contains multiple alignment files. Additionally, a third folder (02_tRNA_unaligned) holds fastq files containing reads that did not align to the tRNA-genome. 

|Files|Content|Rule|
|-----|-------|----|
|`.srt.bam`|Bam files sorted by coordinate and filtered for reads between 15-90 bp|`rule tRNA_alignment`|
|`.alignment.log.txt`|Bowtie2 alignment log|`rule tRNA_alignment`|
|`.mkdup.bam`|Filtered and sorted bam file with duplicates flagged|`rule tRNA_mark_duplicates`|
|`.mkdup.log.txt`|Picard Markduplicates log file|`rule tRNA_mark_duplicates`|
|`stats/.mkdup.bam.idxstats`|Samtools idxstats metrics|`rule tRNA_map_stats`|
|`stats/.mkdup.bam.flagstat`|Samtools flagstat metrics|`rule tRNA_map_stats`|
|`full_alignment_read_length_distribution.txt`|Per-sample read length distribution information for all reads|`read_length_distribution`|

### 03 tRNA Counts
Contains detailed per-sample tRNA count information at the gene and isoform level as well as coverage and mismatch information.

|Files|Content|Rule|
|-----|-------|----|
|`genetype_counts.txt`||`rule tRNA_count`|
|`tRNA_isotype_counts.txt`|Per sample counts of tRNA isoforms|`rule tRNA_count`|
|`gene_level_counts_detailed.txt`|Per sample counts of tRNAs broken down by read type (5', 3' antisense, other) and genome smRNAs/tRNA-loci|`rule tRNA_count`|
|`gene_level_counts_collapsed.txt`|Per sample counts of tRNAs collapsed to isodecoder level and genome smRNAs|`rule tRNA_count`|
|`tRNA_ends_counts.txt`|Contains counts for CCA, CC, and C ending tRNAs|`rule tRNA_count`|
|`mature_tRNA_coverages.txt`|Contains normalized coverages for tRNAs by Sprinzl position|`rule get_tRNA_coverage`|
|`Z_scored_summed_mismatches_relative_to_*.png`|Contains normalized mismatch information at each position, summed across replicates, relative to your reference group, visualized as a Z-scaled heatmap|`rule get_mismatches`|
|`mismatch_figures`|Folder containing normalized mismatch information at each position, summed across replicates, relative to your reference group, split by isoacceptor family, vlisualized as a Z-scaled and hierarchically clustered heatmap|`rule get_mismatches`|

### 04 smRNA Counts
Contains detailed per-sample smRNA count information at the gene biotype level for tRNAs and other smRNAs.

|Files|Content|Rule|
|-----|-------|----|
|`raw_amino_counts_by_group.txt`|Per-greoup raw total counts for all amino acids|`rule smRNA_count`|
|`read_length_distribution.txt`|Per-sample read length distribution for tRNAs split by pre-tRNAs and mature tRNAs|`rule smRNA_count`|
|`smRNA_raw_counts_by_group.txt`|Per-group raw total counts for all smRNA gene biotypes|`rule smRNA_count`|
|`smRNA_raw_counts_by_sample.txt`|Per-sample raw total counts for all smRNA gene biotypes|`rule smRNA_count`|
|`raw_anticodon_counts_by_sample.txt`|Per-sample raw total counts for all anticodons|`rule smRNA_count`|
|`subgroup_counts.txt`|Per-sample raw total counts grouped by tRNA or non-tRNA|`rule smRNA_count`| 

### 05 Normalized
Contains normalization information, and normalized tRNA counts.

|Files|Content|Rule|
|-----|-------|----|
|`gene_level_counts_size_factors.csv`|Per-sample size factors calculated by median-of-ratios (DESeq2) for tRNA isodecoders + genomic smRNAs|`rule normalize_and_PCA`|
|`normalized_gene_level_counts.csv`|Per-sample normalized counts calculated by rlog transformation (DESeq2) for tRNA isodecoders + genomic smRNAs|`rule normalize_and_PCA`|
|`tRNA_isotype_counts_size_factors`|Per-sample size factors calculated by median-of-ratios (DESeq2) for tRNA isotypes only|`rule normalize_and_PCA`|
|`normalized_tRNA_isotype_counts.csv`|Per-sample normalized counts calculated by rlog transformation (DESeq2) for tRNA isotypes only|`rule normalize_and_PCA`|
|`CCA_ends_normalized.csv`|CCA-end counts normalized using tRNA-isotype dataset size factors|`rule plot_counts`|

### 06 PCA
Contains tabular data for principal component analysis and associated PCA plots.

|Files|Content|Rule|
|-----|-------|----|
|`gene_level_variance_plot.png`|Shows ranked gene variances for rlog-normalized tRNA isodecoders + genomic smRNAs|`rule normalize_and_PCA`|
|`gene_level_loadings.csv`|PCA loadings per sample based on 500 most variable features for rlog-normalized tRNA isodecoders + genomic smRNAs|`rule normalize_and_PCA`|
|`gene_level_PCA.png`|PCA plot based on 500 most variable features for rlog-normalized tRNA isodecoders + genomic smRNAs|`rule normalize_and_PCA`|
|`tRNA_isotype_variance_plot.png`|Shows ranked gene variance for rlog-normalized tRNA isotypes only|`rule normalize_and_PCA`|
|`tRNA_isotype_loadings.csv`|PCA loadings per sample based on 500 most variable features for rlog-normalized tRNA isotypes only|`rule normalize_and_PCA`|
|`tRNA_isotype_PCA.png`|PCA plot based on 500 most variable features for rlog-normalized tRNA isotypes only|`rule normalize_and_PCA`|
|`PCA_Analysis_Summary.png`|Side by side PCA plots for full gene-level and tRNA only analysis|`rule normalize_and_PCA`|
|`PCA_Direct_Comparison.png`|PCA data for both experiments on one plot to highlight shifts|`rule normalize_and_PCA`|


### 07 RDS Files

Contains DESeq2 serialized R-data objects to be used in the differential expression workflow.

|Files|Content|Rule|
|-----|-------|----|
|`gene_level_DESeq2_object.Rds`|DESeq dataset object for tRNA isodecoders + genomic smRNAs dataset with (reference level set in config)|`rule normalize_and_PCA`|
|`tRNA_isotype_DESeq2_object.Rds`|DESeq dataset object for tRNA isotypes only (reference level set in config)|`rule normalize_and_PCA`|

### 08 Plots
Contains various plots to explore counts.

|Files|Content|Rule|
|-----|-------|----|
`CCA_ends_normalized_absolute_abundances.png`|Absolute abundances of CCA end types per sample, normalized using tRNA-isotype size factors with resulting normalized tabular data saved in 05_normalized folder|`rule plot_counts`|
|`CCA_ends_relative_abundances.png`|Relative abundances (calculated from non-normalized CCA-end counts) of CCA end types per sample|`rule plot_counts`|
|`Grouped_boxplot_norm_tRNA_isotypes_by_Sample_and_Anticodon.png`|High detial plot showing counts of each tRNA isotype for specific codons|`rule plot_counts`|
|`Isoacceptor_counts_normalized.png`|Normalized absolute abundance of each isoacceptor type|`rule plot_counts`|
|`Isoacceptor_counts_by_sample_normalized.png`|Absolute composition of amino-acids by Sample|`rule plot_counts`|
|`mature_tRNA_coverages`|Folder containing multiple pngs (1 per mature tRNA feature detected), showing normalized coverage per Sprinzl position, split by sample|`rule plot_counts`|

### 09 QC

Contains general QC summary in html format with associated data tables

|Files|Content|Rule|
|-----|-------|----|
|`tRNA_multi_QC_report_data`|Folder containing multi-qc data|`rule all`|
|`tRNA_multi_QC_report.html`|Web viewable html summary of QC metrics for trimming and alignment|`rule all`|



## Differential Expression Implementation
Differential expression is conducted using [DESeq2](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html). Rds files containing differential expression results are generated in the preprocessing module are analyzed in this module. 
This code expects two Rds files, `gene_level_DESeq2_object.Rds` and `tRNA_isotype_DESeq2_object.Rds` within your `07_rds_files` folder. 

DESeq2 results undergo adaptive shrinkage (ashr) which is a method for empirical Bayes shrinkage estimation. Ashr estimates effect sizes (e.g. log2 fold changes) and adjusted standard errors by shrinking noisy estimates toward a central value with the goal of reducing noise and avoiding over-estimation of effects when counts are low, or high. This approach is adaptive as shrinkage depends on the standard error of each estimate and the estimated global distirbution of true effects. 

Each dataset produces two main plots: an MA plot and a volcano plot. MA plots show the average expression (A) in log10 scale on the X-axis and log2 fold change (M) on the y-axis. Genes on the right side of the plot are highly expressed while those on the left are lowly expressed. Genes are also sized by their significance (-log10 padj). Larger points indicate greater significance. 

The differential expression module consists of a single rule which produces the following outputs:

### 10 Differential Expression Outputs

|Files|Content|Rule|
|-----|-------|----|
|`gene_level_DESeq2_shrink_Results.csv`|Contains DESeq2 results and normalized counts per sample (Median-of-ratios) for gene-level dataset|`rule extract_DESeq2_results`|
|`tRNA_isotype_DESeq2_shrink_Results.csv`|Contains DESeq2 results and normalized counts per sample (Median-of-ratios) for tRNA isotype dataset|`rule extract_DESeq2_results`|
|`Figures`|Folder containing MA and volcano plots for both datasets|`rule extract_DESeq2_results`|

### 11 Clover Seq Report

Upon completion of the differential expression module, a pdf report `clover-seq-report.pdf` is generated. 

<img src="img/report_example.png" alt="Description" width="400" height="400"/>

## Development Notes:

To run custom_mismatches.py (works on a single bam file at a time) and requires the mature tRNA fasta, bed, and stk alignment
Originally, the code that calculates the mismatches happens in `processsamples.py` (original tRAX) which calls `getcoverage.testmain` to calculate mismathces. When I try to just run `getcoverage.testmain` this does not work. `getcoverage.main` DOES work, but the file it outputs is slightly different (just has coverage information and not mismatch information). This is useful to make the coverage plots, but we don't have any mismatch information.
```shell
 python code/mismatches.py \
    --bam 02_tRNA_alignment/IB1.mkdup.bam \
    --fasta /dartfs-hpc/rc/lab/G/GMBSR_bioinfo/genomic_references/tRAX_databases/hg38_db/db-maturetRNAs.fa \
    --bed /dartfs-hpc/rc/lab/G/GMBSR_bioinfo/genomic_references/tRAX_databases/hg38_db/db-maturetRNAs.bed \
    --stk /dartfs-hpc/rc/lab/G/GMBSR_bioinfo/genomic_references/tRAX_databases/hg38_db/db-trnaalign.stk >TEST.out
```

The original `getgenomicmismatches.py` from tRAX github seems to work, but rather than tRNAs, it works on all other smRNAs.
Not enitely sure how to hack it to work for tRNAs. See example file on desktop. This is more robust, but I have no idea what's going on....

```shell
python code/mismatches.py \
    --samplefile=runfile.txt \
    -genomefasta=/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/genomic_references/tRAX_databases/hg38_db/db-maturetRNAs.fa \
    --ensemblgtf=/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/genomic_references/tRAX_databases/hg38_db/genes.gtf \
    --covfile=TEST_coverages.txt \
    --outbed=TEST_out.bed \
    --cores=8 \
    --sizefactors=05_normalized/gene_level_counts_size_factors.csv \
    --bedfile=/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/genomic_references/tRAX_databases/hg38_db/db-maturetRNAs.bed \
    --stk=/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/genomic_references/tRAX_databases/hg38_db/db-trnaalign.stk 


```

```

To make all feature bed file (not sure if this is super relevant...)
```shell
python make_all_feature_bed.py \
    db-maturetRNAs.bed \
    db-trnaloci.bed \
    genes.gtf
```

## Contact
Please address questions to **DataAnalyticsCore@groups.dartmouth.edu** or submit an issue in the GitHub repository.

## Citation and Licensing

**This codebase is adapted from the [original tRAX tool](https://github.com/UCSC-LoweLab/tRAX), licensed under GPL v3.0., modified by Mike Martinez, Dartmouth Genomic Data Science Core** 

**Citation:** [Holmes AD, Howard JM, Chan PP, and Lowe TM.](https://www.biorxiv.org/content/10.1101/2022.07.02.498565v1)

**This pipeline was created with funds from the COBRE grant 1P20GM130454. If you use the pipeline in your own work, please acknowledge the pipeline by citing the grant number in your manuscript.**




