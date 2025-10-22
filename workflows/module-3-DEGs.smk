#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Clover-Seq tRNA differential expression workflow
#
# This code was modified from tRAX (doi: 10.1101/2022.07.02.498565)
# 
# Modified by Mike Martinez (Genomic Data Science Core - Dartmouth)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SET GLOBAL SCOPE PYTHON VARIABLES (EXECUTED BEFORE SNAKEMAKE)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

import pandas as pd 
import csv

#----- Read in the sample data
samples_df = pd.read_table(config["sample_txt"], delimiter = ",").set_index("Sample_ID", drop = False)
sample_list = list(samples_df["Sample_ID"])
genome = config["genome"]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SNAKEMAKE RULES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#---- Final rule
rule all:
    input:
        #----- Rule extract_DESeq2_results outputs
        expand("10_differential_expression/{file}", file = [
            "gene_level_DESeq2_shrink_Results.csv",
            "tRNA_isotype_DESeq2_shrink_Results.csv"]),

    output:
        "11_report/clover-seq-report.pdf"
    conda: "env_config/clover-report.yaml"
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb"
    benchmark: "benchmarks/module_3_all/rule_all_bm.tsv"
    params:
    shell: """
    
        #----- Run dummy command
        Rscript -e 'rmarkdown::render("code/clover-seq-report.Rmd")'
        mv "code/clover-seq-report.pdf" {output}

    
    """

#----- Rule to generate differential expression results
rule extract_DESeq2_results:
    output:
        geneLevelRes = "10_differential_expression/gene_level_DESeq2_shrink_Results.csv",
        isotypeRes = "10_differential_expression/tRNA_isotype_DESeq2_shrink_Results.csv"
    conda: "env_config/clover-seq.yaml"
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb"
    benchmark: "benchmarks/rule_extract_DESeq2_results/extract_DESeq2_results_bm.tsv"
    params:
        extract_DESeq2 = "code/extract_DESeq2_results.R",
        metadata = config["sample_txt"],        
        refLevel = config["refLevel"],
        log2FC_magnitude_threshold = config["log2FC_magnitude_threshold"],
        padj_threshold = config["padj_threshold"]
    shell: """
    
        #----- Run the R script to extract DESeq2 results
        Rscript {params.extract_DESeq2} \
            07_rds_files/ \
            {params.metadata} \
            {params.refLevel} \
            {params.log2FC_magnitude_threshold} \
            {params.padj_threshold}

    """

