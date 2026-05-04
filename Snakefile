#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Clover-Seq version 3.0
# This code was modified from tRAX (doi: 10.1101/2022.07.02.498565)
# and adapted for use on Dartmouth HPC via Snakemake.
# tRNA databases are built and hosted by the Dartmouth 
# Genomic Data Science Core (GDSC)
#
# Author: Mike Martinez
# Date last updated: May 2nd, 2026
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SET GLOBAL SCOPE PYTHON VARIABLES (EXECUTED BEFORE SNAKEMAKE)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# TO DO
# - Optimize count_all_smRNA.py
# - Optimize getCoverage.py
# - Optimize getMismatches.py
# - Make sure DESeq2 script is outputting png and not pdf
# - Custom MultiQC reporting and branding
# - Change clusterProfile and job script back to use standard partition

import pandas as pd
import csv

#----- Read in the sample data
samples_df = pd.read_table(config["sample_txt"], delimiter=",").set_index("Sample_ID", drop=False)
sample_list = list(samples_df["Sample_ID"])
genome = config["genome"]

#----- Optionally include database-building rules
build_database = config.get("build_database", False)
if build_database:
    include: "additional_rules/build_database.smk"

#----- Outputs required when build_database: true
db_build_outputs = [
    f"{genome}_db/genes.gtf",
    f"{genome}_db/genome.fa",
    f"{genome}_db/genome.fa.fai",
    f"{genome}_db/{genome}-tRNAs.fa",
    f"{genome}_db/{genome}-tRNAs.bed",
    f"{genome}_db/{genome}-tRNAs-detailed.ss",
    f"{genome}_db/{genome}-tRNAs-detailed.out",
    f"{genome}_db/{genome}-tRNAs-confidence-set.ss",
    f"{genome}_db/{genome}-tRNAs_name_map.txt",
    f"{genome}_db/{genome}-mature-tRNAs.fa",
    f"{genome}_db/{genome}-filtered-tRNAs.fa",
    f"{genome}_db/db-trnatable.txt",
    f"{genome}_db/db-trnaloci.stk",
    f"{genome}_db/db-trnaloci.bed",
    f"{genome}_db/db-trnaalign.stk",
    f"{genome}_db/db-maturetRNAs.fa",
    f"{genome}_db/db-maturetRNAs.bed",
    f"{genome}_db/db-locusnum.txt",
    f"{genome}_db/db-dbinfo.txt",
    f"{genome}_db/db-alignnum.txt",
    f"{genome}_db/db-tRNAgenome.fa",
    f"{genome}_db/loci_sequences.fa",
    f"{genome}_db/db-tRNAgenome.1.bt2l",
    f"{genome}_db/db-tRNAgenome.2.bt2l",
    f"{genome}_db/db-tRNAgenome.3.bt2l",
    f"{genome}_db/db-tRNAgenome.4.bt2l",
    f"{genome}_db/db-tRNAgenome.rev.1.bt2l",
    f"{genome}_db/db-tRNAgenome.rev.2.bt2l"
] if build_database else []

#----- Generate run script for read counting
def generate_runfile(sample_file):
    with open(sample_file, 'r') as infile, open("runfile.txt", 'w') as outfile:
        reader = csv.DictReader(infile)
        for row in reader:
            sample_id = row["Sample_ID"]
            group     = row['Group']
            outfile.write(f"{sample_id} {group} 02_tRNA_alignment\n")

generate_runfile(config["sample_txt"])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SNAKEMAKE RULES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Final rule
rule all:
    input:
        #----- Rule trimming outputs
        expand("01_trimming/{sample}.R1.trim.fastq.gz", sample=sample_list),

        #----- Rule tRNA_bowtie2 / tRNA_choosemappings outputs
        expand("02_tRNA_alignment/{sample}.srt.bam", sample=sample_list),

        #----- Rule tRNA_map_stats outputs
        expand("02_tRNA_alignment/stats/{sample}.srt.bam.idxstats", sample=sample_list),
        expand("02_tRNA_alignment/stats/{sample}.srt.bam.flagstat", sample=sample_list),

        #----- Rule tRNA_count outputs
        expand("03_Raw_Quant/tRNA_counts/{file}", file=[
            "genetype_counts.txt",
            "tRNA_isotype_counts.txt",
            "gene_level_counts_detailed.txt",
            "gene_level_counts_collapsed.txt",
            "tRNA_ends_counts.txt",
            "unique_tRNA_counts.txt"]),

        #----- Rule get_mismatches outputs
        "05_Mismatches/mature_tRNA_mismatches.txt",
        "05_Mismatches/mature_tRNA_mismatches.bed",
        "05_Mismatches/heatmaps",

        #----- Rule count_smRNAs outputs (includes read length distribution)
        "03_Raw_Quant/raw_amino_counts_by_group.txt",
        "03_Raw_Quant/raw_anticodon_counts_by_sample.txt",
        expand("03_Raw_Quant/other_smRNAs/{file}", file=[
            "read_length_distribution.txt",
            "smRNA_raw_counts_by_group.txt",
            "smRNA_raw_counts_by_sample.txt",
            "subroup_counts.txt"]),

        #----- Rule normalize_and_PCA outputs
        expand("04_Expression/{file}", file=[
            "gene_level_counts_size_factors.csv",
            "normalized_gene_level_counts.csv",
            "tRNA_isotype_counts_size_factors.csv",
            "normalized_tRNA_isotype_counts.csv",
            "gene_level_DESeq2_object.Rds",
            "tRNA_isotype_DESeq2_object.Rds"]),
        expand("07_Plots/PCA/{file}", file=[
            "gene_level_variance_plot.png",
            "gene_level_loadings.csv",
            "gene_level_PCA.png",
            "tRNA_isotype_variance_plot.png",
            "tRNA_isotype_loadings.csv",
            "tRNA_isotype_PCA.png",
            "PCA_Analysis_Summary.png"]),

        #----- Rule get_tRNA_coverage outputs
        "06_Coverages/mature_tRNA_coverages.txt",
        "06_Coverages/pretRNA_locus_coverages.txt",

        #----- Rule plot_counts outputs
        expand("07_Plots/{file}", file=[
            "Grouped_boxplot_norm_tRNA_isotypes_by_Sample_and_Anticodon.png",
            "Isoacceptor_counts_by_sample_normalized.png",
            "Isoacceptor_counts_normalized.png",
            "CCA_ends_Relative_Abundances.png",
            "CCA_ends_normalized_absolute_abundances.png",
            "smRNA_Relative_Abundances.png"]),

        #----- Rule generate_mqc_content outputs
        "08_QC/mqc_custom_content/trna_isotype_abundance_mqc.tsv",
        "08_QC/mqc_custom_content/cca_tail_status_mqc.tsv",
        "08_QC/mqc_custom_content/smrna_biotype_mqc.tsv",

        #----- Database build outputs (only when build_database: true)
        db_build_outputs

    output:
        "08_QC/tRNA_multi_QC_report.html"
    message: "Generating MultiQC report"
    conda: "env_config/clover-seq.yaml"
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb"
    params:
        genome = config["genome"]
    shell: """

        #----- Run MultiQC report
        multiqc \
            01_trimming/logs \
            02_tRNA_alignment \
            02_tRNA_alignment/stats \
            03_Raw_Quant/tRNA_counts \
            08_QC/mqc_custom_content \
            -n 08_QC/tRNA_multi_QC_report.html \
            -c multiqc_config.yaml

        #----- Clean redundant file
        if [ -f 03_Raw_Quant/tRNA_counts/unique_tRNA_counts.txt ]; then
            rm 03_Raw_Quant/tRNA_counts/unique_tRNA_counts.txt
        fi

    """

#----- Rule to trim adapters
rule trimming:
    output:
        trim_1 = "01_trimming/{sample}.R1.trim.fastq.gz",
    log:     "01_trimming/logs/{sample}.cutadapt.log"
    message: "Trimming adapters: {wildcards.sample}"
    conda: "env_config/clover-seq.yaml"
    resources: cpus="8", maxtime="2:00:00", mem_mb="60gb"
    params:
        sample    = lambda wildcards: wildcards.sample,
        fastq_1   = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_1"],
        adapter_1 = config["adapter_1"],
        minlength = config["minlength"]
    shell: """

        #----- Run cutadapt (stdout → report, stderr → log)
        cutadapt \
            -o {output.trim_1} \
            {params.fastq_1} \
            -m {params.minlength} \
            -a {params.adapter_1} \
            -j {resources.cpus} > {log} 
    """

#----- Rule 1: align with bowtie2, write name-ordered raw BAM
rule tRNA_bowtie2:
    input:
        trim_1 = "01_trimming/{sample}.R1.trim.fastq.gz"
    output:
        rawBam = temp("02_tRNA_alignment/{sample}.raw.bam")
    log:     "02_tRNA_alignment/logs/{sample}.bowtie2.log"
    message: "Bowtie2 alignment: {wildcards.sample}"
    conda: "env_config/clover-seq.yaml"
    resources: cpus="10", maxtime="6:00:00", mem_mb="60gb"
    params:
        bt2_index = config["bt2_index"],
        maxMaps   = config["maxMaps"],
        nPenalty  = config["nPenalty"]
    shell: """

        bowtie2 \
            --local \
            -x {params.bt2_index} \
            -U {input.trim_1} \
            -k {params.maxMaps} \
            --very-sensitive \
            --np {params.nPenalty} \
            --ignore-quals \
            --reorder \
            -p {resources.cpus} \
            2> {log} \
        | samtools view -b -o {output.rawBam} -

    """

#----- Rule 2: select best tRNA mappings, coordinate-sort, index
rule tRNA_choosemappings:
    input:
        rawBam = "02_tRNA_alignment/{sample}.raw.bam"
    output:
        srtBam = "02_tRNA_alignment/{sample}.srt.bam"
    log:     "02_tRNA_alignment/logs/{sample}.choosemappings.log"
    message: "Selecting best tRNA mappings: {wildcards.sample}"
    conda: "env_config/clover-seq.yaml"
    resources: cpus="10", maxtime="6:00:00", mem_mb="60gb"
    params:
        sample  = lambda wildcards: wildcards.sample,
        tRNA_db = config["trna_db"]
    shell: """

        python code/choosemappings.py {params.tRNA_db}/db-trnatable.txt \
            --input {input.rawBam} \
            --progname TRAX \
            --fqname {params.sample} \
            --expname {params.sample} \
            --minnontrnasize 20 \
            2> {log} \
        | samtools sort -@ {resources.cpus} -o {output.srtBam} -

        samtools index {output.srtBam}

    """

#----- Rule to collate tRNA mapping statistics
rule tRNA_map_stats:
    input:
        bam = "02_tRNA_alignment/{sample}.srt.bam"
    output:
        idxStats  = "02_tRNA_alignment/stats/{sample}.srt.bam.idxstats",
        flagStats = "02_tRNA_alignment/stats/{sample}.srt.bam.flagstat"
    message: "Collecting alignment statistics: {wildcards.sample}"
    conda: "env_config/clover-seq.yaml"
    resources: cpus="4", maxtime="1:00:00", mem_mb="16gb"
    shell: """

        samtools idxstats {input.bam} > {output.idxStats}
        samtools flagstat {input.bam} > {output.flagStats}

    """

#----- Rule to count tRNAs
rule tRNA_count:
    input:
        expand("02_tRNA_alignment/{sample}.srt.bam", sample=sample_list)
    output:
        genetypeFile        = "03_Raw_Quant/tRNA_counts/genetype_counts.txt",
        tRNA_isotype_counts = "03_Raw_Quant/tRNA_counts/tRNA_isotype_counts.txt",
        trnaCountsDetailed  = "03_Raw_Quant/tRNA_counts/gene_level_counts_detailed.txt",
        trnaCountsCollapsed = "03_Raw_Quant/tRNA_counts/gene_level_counts_collapsed.txt",
        trnaEnds            = "03_Raw_Quant/tRNA_counts/tRNA_ends_counts.txt",
        trnaUniqueCounts    = "03_Raw_Quant/tRNA_counts/unique_tRNA_counts.txt"
    message: "Counting tRNA reads across all samples"
    conda: "env_config/clover-seq.yaml"
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb"
    params:
        countScript = "code/countreads.py",
        runFile     = config["runFile"],
        trna_db     = config["trna_db"]
    shell: """

        #----- Run the countreads.py script
        python {params.countScript} \
            --samplefile={params.runFile} \
            --ensemblgtf={params.trna_db}/genes.gtf \
            --trnaloci={params.trna_db}/db-trnaloci.bed \
            --maturetrnas={params.trna_db}/db-maturetRNAs.bed \
            --trnatable={params.trna_db}/db-trnatable.txt \
            --removepseudo \
            --genetypefile={output.genetypeFile} \
            --countfile={output.trnaCountsDetailed} \
            --trnaends={output.trnaEnds} \
            --trnacounts={output.tRNA_isotype_counts} \
            --trnauniquecounts={output.trnaUniqueCounts} \
            --cores={resources.cpus}

        #----- Collapse detailed tRNA counts to isotype-level counts
        awk '
        NR==1 {{
            print;
            next
        }}
        {{
            split($1, arr, "_");
            base = arr[1];
            if (!(base in seen)) {{
                order[++count] = base;
                seen[base] = 1;
            }}
            for (i=2; i<=NF; i++) {{
                counts[base, i] += $i;
            }}
        }}
        END {{
            for (i=1; i<=count; i++) {{
                tRNA = order[i];
                printf "%s", tRNA;
                for (j=2; j<=NF; j++) {{
                    printf "\t%d", counts[tRNA,j] + 0;
                }}
                print "";
            }}
        }}
        ' {output.trnaCountsDetailed} > {output.trnaCountsCollapsed}
    """

#----- Rule to count read types (smRNAs, amino acids, anticodons) — runs after normalization
#      to use size factors, matching processsamples.py counttypes order
rule count_smRNAs:
    input:
        bams        = expand("02_tRNA_alignment/{sample}.srt.bam", sample=sample_list),
        sizefactors = "04_Expression/gene_level_counts_size_factors.csv"
    output:
        aminoCounts     = "03_Raw_Quant/raw_amino_counts_by_group.txt",
        anticodonCounts = "03_Raw_Quant/raw_anticodon_counts_by_sample.txt",
        readLengths     = "03_Raw_Quant/other_smRNAs/read_length_distribution.txt",
        groupCounts     = "03_Raw_Quant/other_smRNAs/smRNA_raw_counts_by_group.txt",
        counts          = "03_Raw_Quant/other_smRNAs/smRNA_raw_counts_by_sample.txt",
        subGroupFile    = "03_Raw_Quant/other_smRNAs/subroup_counts.txt"
    message: "Counting smRNA reads across all samples"
    conda: "env_config/clover-seq.yaml"
    resources: cpus="12", maxtime="6:00:00", mem_mb="60gb"
    params:
        smRNA_count = "code/count_all_smRNA.py",
        runFile     = config["runFile"],
        trna_db     = config["trna_db"]
    shell: """

        #----- Count all tRNA + smRNA (size-factor normalised, matching tRAX counttypes)
        python {params.smRNA_count} \
            --samplefile={params.runFile} \
            --sizefactors={input.sizefactors} \
            --trnatable={params.trna_db}/db-trnatable.txt \
            --ensemblgtf={params.trna_db}/genes.gtf \
            --trnaloci={params.trna_db}/db-trnaloci.bed \
            --maturetrnas={params.trna_db}/db-maturetRNAs.bed \
            --trnaaminofile={output.aminoCounts} \
            --readlengthfile={output.readLengths} \
            --realcountfile={output.counts} \
            --countfile={output.groupCounts} \
            --mismatchfile={output.subGroupFile} \
            --trnaanticodonfile={output.anticodonCounts} \
            --cores={resources.cpus}

    """

#----- Rule to normalize counts and generate PCA
rule normalize_and_PCA:
    input:
        geneLevelCounts = "03_Raw_Quant/tRNA_counts/gene_level_counts_collapsed.txt",
        isoformCounts   = "03_Raw_Quant/tRNA_counts/tRNA_isotype_counts.txt"
    output:
        "04_Expression/gene_level_counts_size_factors.csv",
        "04_Expression/normalized_gene_level_counts.csv",
        "07_Plots/PCA/gene_level_variance_plot.png",
        "07_Plots/PCA/gene_level_loadings.csv",
        "07_Plots/PCA/gene_level_PCA.png",
        "04_Expression/tRNA_isotype_counts_size_factors.csv",
        "04_Expression/normalized_tRNA_isotype_counts.csv",
        "07_Plots/PCA/tRNA_isotype_variance_plot.png",
        "07_Plots/PCA/tRNA_isotype_loadings.csv",
        "07_Plots/PCA/tRNA_isotype_PCA.png",
        "07_Plots/PCA/PCA_Analysis_Summary.png",
        "04_Expression/gene_level_DESeq2_object.Rds",
        "04_Expression/tRNA_isotype_DESeq2_object.Rds"
    message: "Normalizing counts and generating PCA"
    conda: "env_config/clover-seq.yaml"
    resources: cpus="12", maxtime="6:00:00", mem_mb="60gb"
    params:
        #PCAScript = "code/visualizations/clover-seq-normalize-and-PCA.R",
        DESeqScript = "code/clover-seq-DESeq2.R",
        metadata  = config["sample_txt"],
        refLevel  = config["refLevel"]
    shell: """

        #----- Run normalization and PCA script
        Rscript {params.DESeqScript} \
            {params.metadata} \
            {params.refLevel}
    """

#----- Rule to get tRNA coverages (mature + pre-tRNA locus), matching processsamples.py gettrnacoverage
rule get_tRNA_coverage:
    input:
        sizeFactors = "04_Expression/gene_level_counts_size_factors.csv"
    output:
        coverages     = "06_Coverages/mature_tRNA_coverages.txt",
        lociCoverages = "06_Coverages/pretRNA_locus_coverages.txt"
    message: "Calculating tRNA coverages"
    conda: "env_config/clover-seq.yaml"
    resources: cpus="12", maxtime="6:00:00", mem_mb="60gb"
    params:
        coverageCode = "code/getcoverage.py",
        runFile      = config["runFile"],
        trna_db      = config["trna_db"]
    shell: """

        #----- Run coverage calculation (testmain mode via --trnafasta)
        python {params.coverageCode} \
            --bedfile={params.trna_db}/db-maturetRNAs.bed \
            --locibed={params.trna_db}/db-trnaloci.bed \
            --samplefile={params.runFile} \
            --stkfile={params.trna_db}/db-trnaalign.stk \
            --locistk={params.trna_db}/db-trnaloci.stk \
            --trnafasta={params.trna_db}/db-maturetRNAs.fa \
            --numfile={params.trna_db}/db-alignnum.txt \
            --locinums={params.trna_db}/db-locusnum.txt \
            --sizefactors={input.sizeFactors} \
            --allcoverage={output.coverages} \
            --locicoverage={output.lociCoverages} \
            --lociedgemargin=30 \
            --cores={resources.cpus}
    """

#----- Rule to detect mismatches for mature tRNAs
rule get_mismatches:
    input:
        sizeFactors = "04_Expression/gene_level_counts_size_factors.csv"
    output:
        mismatches = "05_Mismatches/mature_tRNA_mismatches.txt",
        outBed     = "05_Mismatches/mature_tRNA_mismatches.bed",
        heatmapDir = directory("05_Mismatches/heatmaps")
    message: "Detecting tRNA mismatches and generating heatmaps"
    conda: "env_config/clover-seq.yaml"
    resources: cpus="12", maxtime="6:00:00", mem_mb="60gb"
    params:
        mismatchCode = "code/getgenomicmismatches.py",
        heatmapCode  = "code/visualizations/clover-seq-heatmaps.R",
        runFile      = config["runFile"],
        metadata     = config["sample_txt"],
        trna_db      = config["trna_db"],
        refLevel     = config["refLevel"]
    shell: """

        #----- Detect mismatches
        python {params.mismatchCode} \
            --samplefile={params.runFile} \
            --genomefasta={params.trna_db}/db-maturetRNAs.fa \
            --covfile={output.mismatches} \
            --outbed={output.outBed} \
            --cores={resources.cpus} \
            --sizefactors={input.sizeFactors} \
            --bedfile={params.trna_db}/db-maturetRNAs.bed \
            --stkfile={params.trna_db}/db-trnaalign.stk

        #----- Generate mismatch heatmaps
        Rscript {params.heatmapCode} \
            --mismatch={output.mismatches} \
            --trna={params.trna_db}/db-trnatable.txt \
            --samples={params.metadata} \
            --directory={output.heatmapDir}

    """

#----- Rule to generate count plots
rule plot_counts:
    input:
        "04_Expression/normalized_tRNA_isotype_counts.csv",
        "06_Coverages/mature_tRNA_coverages.txt"
    output:
        "07_Plots/Grouped_boxplot_norm_tRNA_isotypes_by_Sample_and_Anticodon.png",
        "07_Plots/Isoacceptor_counts_by_sample_normalized.png",
        "07_Plots/Isoacceptor_counts_normalized.png",
        "07_Plots/CCA_ends_Relative_Abundances.png",
        "07_Plots/CCA_ends_normalized_absolute_abundances.png",
        "07_Plots/smRNA_Relative_Abundances.png"
    message: "Generating count and coverage plots"
    conda: "env_config/clover-seq.yaml"
    resources: cpus="12", maxtime="6:00:00", mem_mb="60gb"
    params:
        plotScript = "code/visualizations/clover-seq-plot-counts.R",
        covPlots   = "code/visualizations/clover-seq-plot-coverages.R",
        metadata   = config["sample_txt"],
        refLevel   = config["refLevel"]
    shell: """

        #----- Run count visualization script
        Rscript {params.plotScript} \
            {params.metadata} \
            {params.refLevel}

        #----- Run coverage visualization script
        Rscript {params.covPlots} \
            06_Coverages/mature_tRNA_coverages.txt
    """

#----- Rule to generate MultiQC custom content files
rule generate_mqc_content:
    input:
        gene_level_counts = "03_Raw_Quant/tRNA_counts/tRNA_isotype_counts.txt",
        ends_counts       = "03_Raw_Quant/tRNA_counts/tRNA_ends_counts.txt",
        biotype_counts    = "03_Raw_Quant/other_smRNAs/smRNA_raw_counts_by_sample.txt"
    output:
        isotype = "08_QC/mqc_custom_content/trna_isotype_abundance_mqc.tsv",
        cca     = "08_QC/mqc_custom_content/cca_tail_status_mqc.tsv",
        biotype = "08_QC/mqc_custom_content/smrna_biotype_mqc.tsv"
    message: "Generating MultiQC custom content"
    conda: "env_config/clover-seq.yaml"
    resources: cpus="1", maxtime="0:30:00", mem_mb="8gb"
    params:
        script = "code/generate_mqc_custom_content.py"
    shell: """

        python {params.script} \
            --gene-level-counts {input.gene_level_counts} \
            --ends-counts       {input.ends_counts} \
            --biotype-counts    {input.biotype_counts} \
            --output-dir        08_QC/mqc_custom_content

    """
