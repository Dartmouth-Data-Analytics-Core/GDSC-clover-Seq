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

        #----- Rule tRNA_align outputs
        expand("02_tRNA_alignment/{sample}.srt.bam", sample=sample_list),

        #----- Rule tRNA_mark_duplicates outputs
        expand("02_tRNA_alignment/duplicates/{sample}.mkdup.bam", sample=sample_list),
        expand("02_tRNA_alignment/duplicates/{sample}.mkdup.log.txt", sample=sample_list),

        #----- Rule tRNA_map_stats outputs
        expand("02_tRNA_alignment/stats/{sample}.mkdup.bam.idxstats", sample=sample_list),
        expand("02_tRNA_alignment/stats/{sample}.mkdup.bam.flagstat", sample=sample_list),

        #----- Rule tRNA_count outputs
        expand("03_Raw_Quant/tRNA_counts/{file}", file=[
            "genetype_counts.txt",
            "tRNA_isotype_counts.txt",
            "gene_level_counts_detailed.txt",
            "gene_level_counts_collapsed.txt",
            "tRNA_ends_counts.txt"]),

        #----- Rule get_mismatches outputs
        "05_Mismatches/mature_tRNA_mismatches.txt",
        "05_Mismatches/mature_tRNA_mismatches.bed",
        "05_Mismatches/heatmaps",

        #----- Rule read_length_distribution outputs
        "02_tRNA_alignment/full_alignment_read_length_distribution.txt",

        #----- Rule count_smRNAs outputs
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

        #----- Rule plot_counts outputs
        expand("07_Plots/{file}", file=[
            "Grouped_boxplot_norm_tRNA_isotypes_by_Sample_and_Anticodon.png",
            "Isoacceptor_counts_by_sample_normalized.png",
            "Isoacceptor_counts_normalized.png",
            "CCA_ends_Relative_Abundances.png",
            "CCA_ends_normalized_absolute_abundances.png",
            "smRNA_Relative_Abundances.png"]),

        #----- Database build outputs (only when build_database: true)
        db_build_outputs

    output:
        "09_QC/tRNA_multi_QC_report.html"
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
            -n 09_QC/tRNA_multi_QC_report.html \
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

#----- Rule to align samples to tRNA database
rule tRNA_align:
    input:
        trim_1 = "01_trimming/{sample}.R1.trim.fastq.gz"
    output:
        srtBam   = "02_tRNA_alignment/{sample}.srt.bam"
    log:     "02_tRNA_alignment/logs/{sample}.tRNA_align.log"
    message: "Aligning to tRNA database: {wildcards.sample}"
    conda: "env_config/clover-bowtie2.yaml"
    resources: cpus="10", maxtime="6:00:00", mem_mb="60gb"
    params:
        sample    = lambda wildcards: wildcards.sample,
        bt2_index = config["bt2_index"],
        maxMaps   = config["maxMaps"],
        nPenalty  = config["nPenalty"]
    shell: """

        #----- Run Bowtie2 alignment (bowtie2 stats → alignLog, pipeline stderr → log)
        bowtie2 \
            --local \
            -x {params.bt2_index} \
            -U {input.trim_1} \
            -k {params.maxMaps} \
            --very-sensitive \
            --np {params.nPenalty} \
            --ignore-quals \
            -p {resources.cpus} \
            -S 02_tRNA_alignment/{params.sample}.aln.sam 2> {log}

        #----- Convert SAM to BAM
        samtools view -bS 02_tRNA_alignment/{params.sample}.aln.sam > 02_tRNA_alignment/{params.sample}.bam

        #----- Sort and index
        samtools sort -@ 4 02_tRNA_alignment/{params.sample}.bam > {output.srtBam}
        samtools index {output.srtBam}

        #----- Remove temp files
        rm -rf 02_tRNA_alignment/{params.sample}.aln.sam
        rm -rf 02_tRNA_alignment/{params.sample}.bam

    """

#----- Rule to mark duplicates
rule tRNA_mark_duplicates:
    input:
        bam = "02_tRNA_alignment/{sample}.srt.bam"
    output:
        mkdup    = "02_tRNA_alignment/duplicates/{sample}.mkdup.bam",
        mkdupLog = "02_tRNA_alignment/duplicates/{sample}.mkdup.log.txt"
    log:     "02_tRNA_alignment/logs/{sample}.picard.log"
    message: "Marking duplicates: {wildcards.sample}"
    conda: "env_config/Picard.yaml"
    resources: cpus="12", maxtime="6:00:00", mem_mb="100gb"
    params:
        sample = lambda wildcards: wildcards.sample
    shell: """

        #----- Run Picard MarkDuplicates
        picard -Xmx16G -Xms16G \
            MarkDuplicates \
            I={input.bam} \
            O={output.mkdup} \
            M={output.mkdupLog} \
            OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 \
            CREATE_INDEX=false \
            MAX_RECORDS_IN_RAM=4000000 \
            ASSUME_SORTED=true \
            MAX_FILE_HANDLES=768 2> {log}

        #----- Index the deduped BAM
        samtools index {output.mkdup}

    """

#----- Rule to collate tRNA mapping statistics
rule tRNA_map_stats:
    input:
        mkdup = "02_tRNA_alignment/duplicates/{sample}.mkdup.bam"
    output:
        idxStats  = "02_tRNA_alignment/stats/{sample}.mkdup.bam.idxstats",
        flagStats = "02_tRNA_alignment/stats/{sample}.mkdup.bam.flagstat"
    message: "Collecting alignment statistics: {wildcards.sample}"
    conda: "env_config/clover-seq.yaml"
    resources: cpus="12", maxtime="6:00:00", mem_mb="60gb"
    params:
        sample = lambda wildcards: wildcards.sample
    shell: """

        #----- Collect alignment metrics
        samtools idxstats {input.mkdup} > {output.idxStats}
        samtools flagstat {input.mkdup} > {output.flagStats}

    """

#----- Rule to count tRNAs
rule tRNA_count:
    input:
        expand("02_tRNA_alignment/duplicates/{sample}.mkdup.bam", sample=sample_list)
    output:
        genetypeFile        = "03_Raw_Quant/tRNA_counts/genetype_counts.txt",
        tRNA_isotype_counts = "03_Raw_Quant/tRNA_counts/tRNA_isotype_counts.txt",
        trnaCountsDetailed  = "03_Raw_Quant/tRNA_counts/gene_level_counts_detailed.txt",
        trnaCountsCollapsed = "03_Raw_Quant/tRNA_counts/gene_level_counts_collapsed.txt",
        trnaEnds            = "03_Raw_Quant/tRNA_counts/tRNA_ends_counts.txt"
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
            --genetypefile={output.genetypeFile} \
            --trnaends={output.trnaEnds} \
            --trnacounts={output.tRNA_isotype_counts} \
            --cores={resources.cpus} > {output.trnaCountsDetailed}

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

#----- Rule to plot read-length distributions for all aligned reads
rule read_length_distribution:
    input:
        expand("02_tRNA_alignment/duplicates/{sample}.mkdup.bam", sample=sample_list)
    output:
        distribution = "02_tRNA_alignment/full_alignment_read_length_distribution.txt"
    message: "Calculating read length distributions"
    conda: "env_config/clover-seq.yaml"
    resources: cpus="12", maxtime="6:00:00", mem_mb="60gb"
    shell: """

        #----- Calculate read length distribution
        echo -e "Length\tSample\tCount" > {output.distribution}

        for bam in {input}; do
            sample=$(basename "$bam" .mkdup.bam)
            samtools view "$bam" | \
                awk -v sample="$sample" '{{
                    len = length($10)
                    counts[len]++
                }}
                END {{
                    for (i = 0; i <= 100; i++) {{
                        printf "%d\t%s\t%d\\n", i, sample, counts[i]+0
                    }}
                }}' >> {output.distribution}
        done
    """

#----- Rule to count other smRNAs
rule count_smRNAs:
    input:
        expand("02_tRNA_alignment/duplicates/{sample}.mkdup.bam", sample=sample_list)
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

        #----- Count all tRNA + smRNA
        python {params.smRNA_count} \
            --samplefile={params.runFile} \
            --trnatable={params.trna_db}/db-trnatable.txt \
            --ensemblgtf={params.trna_db}/genes.gtf \
            --trnaloci={params.trna_db}/db-trnaloci.bed \
            --maturetrnas={params.trna_db}/db-maturetRNAs.bed \
            --trnaaminofile={output.aminoCounts} \
            --readlengthfile={output.readLengths} \
            --realcountfile={output.counts} \
            --countfile={output.groupCounts} \
            --mismatchfile={output.subGroupFile} \
            --trnaanticodonfile={output.anticodonCounts}

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

#----- Rule to get tRNA coverages
rule get_tRNA_coverage:
    input:
        sizeFactors = "04_Expression/gene_level_counts_size_factors.csv"
    output:
        coverages = "06_Coverages/mature_tRNA_coverages.txt"
    message: "Calculating tRNA coverages"
    conda: "env_config/clover-seq.yaml"
    resources: cpus="12", maxtime="6:00:00", mem_mb="60gb"
    params:
        coverageCode = "code/getcoverage.py",
        runFile      = config["runFile"],
        trna_db      = config["trna_db"]
    shell: """

        #----- Run coverage calculation
        python {params.coverageCode} \
            --bedfile={params.trna_db}/db-maturetRNAs.bed \
            --samplefile={params.runFile} \
            --stkfile={params.trna_db}/db-trnaalign.stk \
            --sizefactors={input.sizeFactors} > coverage.tmp.txt

        awk 'NR==1 {{
            for (i = 1; i <= NF; i++) {{
                if ($i ~ /\.gap/ || $i == "-1") skip[i] = 1;
                else keep[++n] = i;
            }}
        }}
        {{
            for (i = 1; i <= n; i++) {{
                printf "%s%s", $(keep[i]), (i < n ? OFS : ORS);
            }}
        }}' OFS="\t" coverage.tmp.txt > {output.coverages}
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
