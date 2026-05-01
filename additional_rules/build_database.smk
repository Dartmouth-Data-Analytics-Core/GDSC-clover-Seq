#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Clover-Seq database building rules
#
# Included in Snakefile when build_database: true in config.
# Genome variable and config are inherited from the main Snakefile.
#
# This code was modified from tRAX (doi: 10.1101/2022.07.02.498565)
#
# Modified by Mike Martinez (Genomic Data Science Core - Dartmouth)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Rule to build tRNA database from downloads
rule generate_gtRNA_db:
    output:
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
        f"{genome}_db/db-alignnum.txt"
    conda: "../workflows/env_config/clover-seq.yaml"
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb"
    params:
        genome      = config["genome"],
        GTF_URL     = config["GTF_URL"],
        gtRNAdb_URL = config["gtRNAdb_URL"],
        gtRNAdb_OUT = config["gtRNAdb_OUT"],
        gtRNAdb_NAME = config["gtRNAdb_NAME"],
        GENOME_URL  = config["GENOME_URL"],
        FA          = config["FASTA"],
        makeDB      = config["makeDB"]
    shell: """
        set -euo pipefail

        #----- GTF file from Ensembl
        echo "Generating GTF"
        wget -q -O - {params.GTF_URL} | \
            gzip -cd | \
            grep -v '^#' | \
            awk '{{print "chr" $0;}}' | \
            sed 's/chrMT/chrM/g' | \
            grep -e Mt_rRNA -e Mt_tRNA -e miRNA -e misc_RNA -e rRNA -e snRNA -e snoRNA -e ribozyme -e sRNA -e scaRNA \
            > {params.genome}_db/genes.gtf
        echo "GTF done"

        #----- Extract gtRNAdb archive
        echo "Downloading gtRNAdb"
        wget --no-check-certificate {params.gtRNAdb_URL} -O {params.genome}_db/tse.tar.gz
        tar -xzvf {params.genome}_db/tse.tar.gz -C {params.genome}_db
        rm {params.genome}_db/tse.tar.gz
        echo "gtRNAdb done"

        #----- Download genome FASTA
        echo "Downloading genome FASTA"
        if [ {params.FA} == "true" ]
        then
            wget -q -O - {params.GENOME_URL} | gzip -cd > {params.genome}_db/genome.fa
        else
            wget -q -O {params.genome}_db/genome.2bit {params.GENOME_URL}
            twoBitToFa {params.genome}_db/genome.2bit {params.genome}_db/genome.fa
        fi
        echo "FASTA done"

        #----- Build tRAX database
        python {params.makeDB} \
            --databasename={params.genome}_db/db \
            --genomefile={params.genome}_db/genome.fa \
            --trnascanfile={params.genome}_db/{params.gtRNAdb_OUT} \
            --namemapfile={params.genome}_db/{params.gtRNAdb_NAME}
    """

#----- Rule to concatenate mature tRNAs with genome and generate loci FASTA
rule concat_tRNAs:
    input:
        maturetRNAs = f"{genome}_db/db-maturetRNAs.fa",
        genome_fa   = f"{genome}_db/genome.fa",
        loci_bed    = f"{genome}_db/db-trnaloci.bed"
    output:
        tRNAgenome  = f"{genome}_db/db-tRNAgenome.fa",
        loci_fasta  = f"{genome}_db/loci_sequences.fa"
    conda: "../workflows/env_config/clover-seq.yaml"
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb"
    shell: """
        #----- Concatenate mature tRNAs and full genome
        cat {input.maturetRNAs} {input.genome_fa} > {output.tRNAgenome}

        #----- Extract loci sequences
        bedtools getfasta \
            -fi {output.tRNAgenome} \
            -bed {input.loci_bed} \
            -fo {output.loci_fasta}
    """

#----- Rule to build bowtie2 index for tRNA-genome
rule tRNA_bt2_index:
    input:
        tRNAgenome = f"{genome}_db/db-tRNAgenome.fa"
    output:
        f"{genome}_db/db-tRNAgenome.1.bt2l",
        f"{genome}_db/db-tRNAgenome.2.bt2l",
        f"{genome}_db/db-tRNAgenome.3.bt2l",
        f"{genome}_db/db-tRNAgenome.4.bt2l",
        f"{genome}_db/db-tRNAgenome.rev.1.bt2l",
        f"{genome}_db/db-tRNAgenome.rev.2.bt2l"
    conda: "../workflows/env_config/clover-bowtie2.yaml"
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb"
    params:
        genome = config["genome"]
    shell: """
        bowtie2-build \
            {input.tRNAgenome} \
            {params.genome}_db/db-tRNAgenome \
            -p {resources.cpus}
    """
