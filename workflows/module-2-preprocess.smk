#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Clover-Seq tRNA data pre-processing workflow
#
# This code was modified from tRAX (doi: 10.1101/2022.07.02.498565)
# 
# Modified by Mike Martinez (Genomic Data Science Core - Dartmouth)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- TO DO

#- Add code and conda environment to make some cool plots at the read counts step
#- Add code to make some 09_QC plots outside of multi09_QC
#- Add prebuilt configs for different host species
#- Built tRNA-only bowtie2 indices
#- Figure out the mitoDB issue
#- Build tRNA + mitoRNA db if we can figure out the mito issue

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SET GLOBAL SCOPE PYTHON VARIABLES (EXECUTED BEFORE SNAKEMAKE)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

import pandas as pd 
import csv

#----- Read in the sample data
samples_df = pd.read_table(config["sample_txt"], delimiter = ",").set_index("Sample_ID", drop = False)
sample_list = list(samples_df["Sample_ID"])
genome = config["genome"]

#----- Generate run script for read counting
def generate_runfile(sample_file):
    with open(sample_file, 'r') as infile, open("runfile.txt", 'w') as outfile:
        reader = csv.DictReader(infile)
        for row in reader:
            sample_id = row["Sample_ID"]
            group = row['Group']
            outfile.write(f"{sample_id} {group} 02_tRNA_alignment\n")

#----- Run function
generate_runfile(config["sample_txt"])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SNAKEMAKE RULES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Final Rule
rule all:
    input:
        #----- Rule trimming outputs
        expand("01_trimming/{sample}.R1.trim.fastq.gz", sample = sample_list),
        expand("01_trimming/logs/{sample}.cutadapt.report", sample = sample_list),

        #----- Rule tRNA_alignment outputs
        expand("02_tRNA_alignment/{sample}.alignment.log.txt", sample = sample_list),
        expand("02_tRNA_alignment/{sample}.srt.bam", sample = sample_list),

        #----- Rule tRNA_mark_duplicates outputs
        expand("02_tRNA_alignment/{sample}.mkdup.bam", sample = sample_list),
        expand("02_tRNA_alignment/{sample}.mkdup.log.txt", sample = sample_list),

        #----- Rule tRNA_map_stats outputs
        expand("02_tRNA_alignment/stats/{sample}.mkdup.bam.idxstats", sample = sample_list),
        expand("02_tRNA_alignment/stats/{sample}.mkdup.bam.flagstat", sample = sample_list),

        #----- Rule tRNA_count outputs
        expand("03_tRNA_counts/{file}", file = [
            "genetype_counts.txt",
            "tRNA_isotype_counts.txt",
            "gene_level_counts_detailed.txt",
            "gene_level_counts_collapsed.txt",
            "tRNA_ends_counts.txt",
            "mature_tRNA_coverages.txt",
            "mature_tRNA_mismatches.txt",
            "mature_tRNA_mismatches.bed"]),

        #----- Rule read_length_distribution outputs
        "02_tRNA_alignment/full_alignment_read_length_distribution.txt",

        #----- Rule count_smRNAs outputs
        expand("04_smRNA_counts/{file}", file = [
            "raw_amino_counts_by_group.txt",
            "read_length_distribution.txt",
            "smRNA_raw_counts_by_group.txt",
            "smRNA_raw_counts_by_sample.txt",
            "subroup_counts.txt",
            "raw_anticodon_counts_by_sample.txt"]),
        
        #----- Rule normalize_and_PCA outputs
        expand("05_normalized/{file}", file = [
            "gene_level_counts_size_factors.csv",
            "normalized_gene_level_counts.csv",
            "tRNA_isotype_counts_size_factors.csv",
            "normalized_tRNA_isotype_counts.csv"]),
        expand("06_PCA/{file}", file = [
            "gene_level_variance_plot.png",
            "gene_level_loadings.csv",
            "gene_level_PCA.png",
            "tRNA_isotype_variance_plot.png",
            "tRNA_isotype_loadings.csv",
            "tRNA_isotype_PCA.png",
            "PCA_Analysis_Summary.png"]),
        expand("07_rds_files/{file}", file = [
            "gene_level_DESeq2_object.Rds",
            "tRNA_isotype_DESeq2_object.Rds"]),
        
        #----- Rule plot_counts outputs
        expand("08_plots/{file}", file = [
            "Grouped_boxplot_norm_tRNA_isotypes_by_Sample_and_Anticodon.png",
            "Isoacceptor_counts_by_sample_normalized.png",
            "Isoacceptor_counts_normalized.png",
            "CCA_ends_Relative_Abundances.png",
            "CCA_ends_normalized_absolute_abundances.png",
            "smRNA_Relative_Abundances.png"]),

    output:
        "09_QC/tRNA_multi_QC_report.html"
    conda: "env_config/clover-seq.yaml"
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb"
    benchmark: "benchmarks/rule_all_bm.tsv"
    params:
        genome = config["genome"],
    shell:"""
    
        #----- Run Multi_QC Report
        multiqc \
            01_trimming/logs \
            02_tRNA_alignment \
            02_tRNA_alignment/stats \
            03_tRNA_counts \
            -n 09_QC/tRNA_multi_QC_report.html \
            -c multiqc_config.yaml

        #----- Clean redundant file
        if [ -f 03_tRNA_counts/unique_tRNA_counts.txt ]; then
            rm 03_tRNA_counts/unique_tRNA_counts.txt
        fi

        #----- Unlock working directory
        rm -f .snakemake/lock
    
    """

#----- Rule to trim
rule trimming:
    output:
        trim_1 = "01_trimming/{sample}.R1.trim.fastq.gz",
        report = "01_trimming/logs/{sample}.cutadapt.report"
    conda: "env_config/clover-seq.yaml"
    resources: cpus="8", maxtime="2:00:00", mem_mb="60gb"
    benchmark: "benchmarks/rule_01_trimming/{sample}_01_trimming_bm.tsv"
    params:
        sample = lambda wildcards: wildcards.sample,
        fastq_1 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_1"],
        adapter_1 = config["adapter_1"],
        minlength = config["minlength"]
    shell: """
    
        #----- Run cutadapt
        cutadapt \
            -o {output.trim_1} \
            {params.fastq_1} \
            -m {params.minlength} \
            -a {params.adapter_1} \
            -j {resources.cpus} > {output.report} 2>&1
    """

#----- Rule to align samples to tRNA database (need db-trnatable.txt)
rule tRNA_align:
    input:
        trim_1 = "01_trimming/{sample}.R1.trim.fastq.gz",
    output:
        #sam = "alignment/{sample}.aln.sam",
        alignLog = "02_tRNA_alignment/{sample}.alignment.log.txt",
        #unalign = "tRNA_unaligned/{sample}.unalign.fastq",
        srtBam = "02_tRNA_alignment/{sample}.srt.bam"
    conda: "env_config/clover-bowtie2.yaml"
    resources: cpus="10", maxtime="6:00:00", mem_mb="60gb"
    benchmark: "benchmarks/rule_tRNA_align/{sample}_tRNA_align_bm.tsv"
    params:
        sample = lambda wildcards: wildcards.sample,
        bt2_index = config["bt2_index"],
        maxMaps = config["maxMaps"],
        nPenalty = config["nPenalty"],
    shell: """
    
        #----- Run Bowtie2 alignment (removed the --reorder arg from original code to speed up...don't see the need with SE data)
        bowtie2 \
            -x {params.bt2_index} \
            -U {input.trim_1} \
            -D 20 -R 3 -N 1 -L 12 -i S,1,0.50 \
            -k {params.maxMaps} \
            --very-sensitive \
            --np {params.nPenalty} \
            --ignore-qual \
            -p {resources.cpus} \
            -S 02_tRNA_alignment/{params.sample}.aln.sam 2> {output.alignLog}

        #----- subset reads for aligned length > 15 & < 90bp 
        samtools view -h 02_tRNA_alignment/{params.sample}.aln.sam | \
            awk 'BEGIN {{OFS="\t"}} $1 ~ /^@/ || ((length($10) > 15 && length($10) <= 90))' | \
            samtools view -Sb - > 02_tRNA_alignment/{params.sample}.bam

        #----- filter for any reads with MAPQ <=1
        samtools view -h -q 2 02_tRNA_alignment/{params.sample}.bam > 02_tRNA_alignment/{params.sample}.sub.bam

        #----- Sort and filter the bam file
        samtools sort -@ 4 02_tRNA_alignment/{params.sample}.sub.bam > {output.srtBam}
        samtools index {output.srtBam}

        #----- Remove temp files
        rm -rf 02_tRNA_alignment/{params.sample}.aln.sam
        rm -rf 02_tRNA_alignment/{params.sample}.bam
        rm -rf 02_tRNA_alignment/{params.sample}.sub.bam
    """

#----- Rule to mark duplicates
rule tRNA_mark_duplicates:
    input:
        bam = "02_tRNA_alignment/{sample}.srt.bam"
    output:
        mkdup = "02_tRNA_alignment/{sample}.mkdup.bam",
        mkdupLog = "02_tRNA_alignment/{sample}.mkdup.log.txt"
    conda: "env_config/Picard.yaml"
    resources: cpus="12", maxtime="6:00:00", mem_mb="60gb"
    benchmark: "benchmarks/rule_tRNA_mark_duplicates/{sample}_tRNA_mark_duplicates_bm.tsv"
    params:
        sample = lambda wildcards: wildcards.sample
    shell: """
    
        #----- Run Picard mark Duplicates
        picard -Xmx16G -Xms16G  \
            MarkDuplicates \
            I={input.bam} \
            O={output.mkdup} \
            M={output.mkdupLog} \
            OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 \
            CREATE_INDEX=false  \
            MAX_RECORDS_IN_RAM=4000000 \
            ASSUME_SORTED=true \
            MAX_FILE_HANDLES=768

        #----- Index the mkdup bam
        samtools index {output.mkdup}
        
    """

#----- Rule to collate tRNA mapping statistics
rule tRNA_map_stats:
    input:
        mkdup = "02_tRNA_alignment/{sample}.mkdup.bam"
    output:
        idxStats = "02_tRNA_alignment/stats/{sample}.mkdup.bam.idxstats",
        flagStats = "02_tRNA_alignment/stats/{sample}.mkdup.bam.flagstat"
    conda: "env_config/clover-seq.yaml"
    resources: cpus="12", maxtime="6:00:00", mem_mb="60gb"
    benchmark: "benchmarks/rule_tRNA_map_stats/{sample}_tRNA_map_stats_bm.tsv"
    params: 
        sample = lambda wildcards: wildcards.sample
    shell: """
    
        #----- Collect metrics
        samtools idxstats {input.mkdup} > {output.idxStats}
        samtools flagstat {input.mkdup} > {output.flagStats}
  
    """

#----- Rule to count tRNAs
rule tRNA_count:
    input:
        expand("02_tRNA_alignment/{sample}.mkdup.bam", sample = sample_list)
        #expand("02_tRNA_alignment/{sample}.mkdup.bam.bai", sample = sample_list)
    output:
        genetypeFile = "03_tRNA_counts/genetype_counts.txt",
        tRNA_isotype_counts = "03_tRNA_counts/tRNA_isotype_counts.txt",
        trnaCountsDetailed = "03_tRNA_counts/gene_level_counts_detailed.txt",
        trnaCountsCollapsed = "03_tRNA_counts/gene_level_counts_collapsed.txt",
        trnaEnds = "03_tRNA_counts/tRNA_ends_counts.txt"
    conda: "env_config/clover-seq.yaml"
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb",
    benchmark: "benchmarks/rule_tRNA_count/tRNA_count_bm.tsv"
    params:
        countScript = "code/countreads.py",
        runFile = config["runFile"],
        trna_db = config["trna_db"]
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
            --trnacounts={output.tRNA_isotype_counts} > {output.trnaCountsDetailed}

        #----- Collapse the detailed tRNA counts to overall tRNA counts
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


#----- Rule to plot read-length distributions for ALL reads (not just tRNAs)
rule read_length_distribution:
    input:
        expand("02_tRNA_alignment/{sample}.mkdup.bam", sample = sample_list)
    output:
        distribution = "02_tRNA_alignment/full_alignment_read_length_distribution.txt"
    conda: "env_config/clover-seq.yaml"
    resources: cpus="12", maxtime="6:00:00", mem_mb="60gb"
    benchmark: "benchmarks/rule_read_length_distribution/read_length_distribution_bm.tsv"
    params:
    shell: """
    
        #----- Calculate read length distribution information
        echo -e "Length\tSample\tCount" > {output.distribution}

        # Loop through all BAM files
        for bam in {input}; do
            # Extract the sample name from the filename
            sample=$(basename "$bam" .mkdup.bam)

            # Generate read length distribution using samtools and awk
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
        expand("02_tRNA_alignment/{sample}.mkdup.bam", sample = sample_list),
    output:
        aminoCounts = "04_smRNA_counts/raw_amino_counts_by_group.txt",
        readLengths = "04_smRNA_counts/read_length_distribution.txt",
        groupCounts = "04_smRNA_counts/smRNA_raw_counts_by_group.txt",
        counts = "04_smRNA_counts/smRNA_raw_counts_by_sample.txt",
        subGroupFile = "04_smRNA_counts/subroup_counts.txt",
        anticodonCounts = "04_smRNA_counts/raw_anticodon_counts_by_sample.txt"
    conda: "env_config/clover-seq.yaml"
    resources: cpus="12", maxtime="6:00:00", mem_mb="60gb"
    benchmark: "benchmarks/rule_count_smRNAs/count_smRNAs_bm.tsv"
    params:
        smRNA_count = "code/count_all_smRNA.py",
        runFile = config["runFile"],
        trna_db = config["trna_db"]

    shell: """
    
        #----- Run the code to count all tRNA + smRNA
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

#----- Rule to output 09_QC plots
rule normalize_and_PCA:
    input:
        geneLevelCounts = "03_tRNA_counts/gene_level_counts_collapsed.txt",
        isoformCounts = "03_tRNA_counts/tRNA_isotype_counts.txt"
    output:
        "05_normalized/gene_level_counts_size_factors.csv",
        "05_normalized/normalized_gene_level_counts.csv",
        "06_PCA/gene_level_variance_plot.png",
        "06_PCA/gene_level_loadings.csv",
        "06_PCA/gene_level_PCA.png",
        "05_normalized/tRNA_isotype_counts_size_factors.csv",
        "05_normalized/normalized_tRNA_isotype_counts.csv",
        "06_PCA/tRNA_isotype_variance_plot.png",
        "06_PCA/tRNA_isotype_loadings.csv",
        "06_PCA/tRNA_isotype_PCA.png",
        "06_PCA/PCA_Analysis_Summary.png",
        "07_rds_files/gene_level_DESeq2_object.Rds",
        "07_rds_files/tRNA_isotype_DESeq2_object.Rds"
    conda: "env_config/clover-seq.yaml"
    resources: cpus="12", maxtime="6:00:00", mem_mb="60gb"
    benchmark: "benchmarks/rule_normalize_and_PCA/normalize_and_PCA_bm.tsv"
    params:
        PCAScript = "code/visualizations/clover-seq-normalize-and-PCA.R",
        metadata = config["sample_txt"],
        refLevel = config["refLevel"]
    shell: """
    
        #----- Run the script
        Rscript {params.PCAScript} \
            {params.metadata} \
            {params.refLevel}
    """

#----- Rule to get coverages
rule get_tRNA_coverage:
    input:
        sizeFactors = "05_normalized/gene_level_counts_size_factors.csv"
    output:
        coverages = "03_tRNA_counts/mature_tRNA_coverages.txt"
    conda: "env_config/clover-seq.yaml"
    resources: cpus="12", maxtime="6:00:00", mem_mb="60gb"
    benchmark: "benchmarks/rule_get_tRNA_coverage/get_tRNA_coverage_bm.tsv"
    params:
        coverageCode = "code/getcoverage.py",
        runFile = config["runFile"],
        trna_db = config["trna_db"]
    shell: """
    
        #----- Run the getcoverage code
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

#----- Rule to get mismatches for mature tRNAs
rule get_mismatches:
    input:
        sizeFactors = "05_normalized/gene_level_counts_size_factors.csv"
    output:
        mismatches = "03_tRNA_counts/mature_tRNA_mismatches.txt",
        outBed = "03_tRNA_counts/mature_tRNA_mismatches.bed"
    conda: "env_config/clover-seq.yaml"
    resources: cpus="12", maxtime="6:00:00", mem_mb="60gb"
    benchmark: "benchmarks/rule_get_mismatches/get_mismatches_bm.tsv"
    params:
        mismatchCode = "code/getgenomicmismatches.py",
        heatmapCode = "code/visualizations/clover-seq-mismatches.R",
        runFile = config["runFile"],
        metadata = config["sample_txt"],
        trna_db = config["trna_db"],
        refLevel = config["refLevel"]
    shell: """
    
        #----- Run the mismatch code code
        python {params.mismatchCode} \
            --samplefile={params.runFile} \
            --genomefasta={params.trna_db}/db-maturetRNAs.fa \
            --ensemblgtf={params.trna_db}/genes.gtf \
            --covfile={output.mismatches} \
            --outbed={output.outBed} \
            --cores={resources.cpus} \
            --sizefactors={input.sizeFactors} \
            --bedfile={params.trna_db}/db-maturetRNAs.bed \
            --stkfile={params.trna_db}/trnaalign.stk

        #----- Run plotting code
        Rscript {params.heatmapCode} \
            {output.mismatches} \
            {params.metadata} \
            {params.refLevel}

    """

#----- Rule to generate plots
rule plot_counts:
    input:
        "05_normalized/normalized_tRNA_isotype_counts.csv",
        "03_tRNA_counts/mature_tRNA_coverages.txt"
    output:
        "08_plots/Grouped_boxplot_norm_tRNA_isotypes_by_Sample_and_Anticodon.png",
        "08_plots/Isoacceptor_counts_by_sample_normalized.png",
        "08_plots/Isoacceptor_counts_normalized.png",
        "08_plots/CCA_ends_Relative_Abundances.png",
        "08_plots/CCA_ends_normalized_absolute_abundances.png",
        "08_plots/smRNA_Relative_Abundances.png"
    conda: "env_config/clover-seq.yaml"
    resources: cpus="12", maxtime="6:00:00", mem_mb="60gb"
    benchmark: "benchmarks/rule_plot_counts/plot_counts_bm.tsv"
    params:
        plotScript = "code/visualizations/clover-seq-plot-counts.R",
        covPlots = "code/visualizations/clover-seq-plot-coverages.R",
        metadata = config["sample_txt"],
        refLevel = config["refLevel"]
    shell: """
    
        #----- Run plotting script
        Rscript {params.plotScript} \
            {params.metadata} \
            {params.refLevel}

        #----- Run coverage plot script
        Rscript {params.covPlots} \
            03_tRNA_counts/mature_tRNA_coverages.txt
    """

#rule coverage:
 #   input:
 #       sizeFactors = "05_normalized/gene_level_counts_size_factors.csv"
   # conda: "env_config/clover-seq.yaml"
   # resources: cpus="12", maxtime="6:00:00", mem_mb="60gb"
   # benchmark: "benchmarks/rule_coverage_coverage_bm.txt"
   # params:
    #    runFile = config["runFile"],
     #   trna_db = config["trna_db"],
    #    covScript = "code/getcoverage.py"
   # shell: """
    
        #----- Run the coverage script
    #    python {params.covScript} \
     #       --bedfile={params.trna_db}/db-maturetRNAs.bed \
     #       --samplefile={params.runFile} \
     #       --stkfile={params.trna_db}/db-trnaalign.stk
      #      --mincoverage=10 \
      #      --sizefactors={input.sizeFactors} > "10_coverages/all_coverages.txt"

    
   #"""
#This works locally....
    #python code/getcoverage.py --bedfile=/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/genomic_references/tRAX_databases/hg38_db/db-maturetRNAs.bed --samplefile=runfile.txt --stkfile=/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/genomic_references/tRAX_databases/hg38_db/db-trnaalign.stk --mincoverage=10 --sizefactors=05_normalized/trna_isotype_counts_size_factors.csv --uniquegenome=unique > output_test.txt
