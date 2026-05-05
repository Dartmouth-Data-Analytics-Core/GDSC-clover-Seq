
# Dartmouth GDSC Clover-Seq Pipeline
![Version](https://img.shields.io/badge/version-3.1.1-blue)

GDSC-Clover-Seq is a Snakemake pipeline for the comprehensive quantitative analysis of tRNA-seq libraries, developed by the Dartmouth Genomic Data Science Core (GDSC). It is adapted from the tRAX framework (Holmes et al., 2022; doi: 10.1101/2022.07.02.498565) and reimplemented as a reproducible, HPC-compatible workflow with per-rule conda environments and SLURM resource management.

The pipeline takes adapter-trimmed short reads and produces a complete tRNA expression atlas. Reads are aligned to a Bowtie2 reference containing mature tRNA sequences downloaded from GtRNADB and the full genome sequence, allowing simultaneous quantification of mature tRNAs and their pre-tRNA precursors. Multi-mapping (common in tRNA-seq because many tRNA gene families share near-identical sequences) is resolved by retaining only the highest-scoring alignment(s) per read. Downstream counting distinguishes uniquely-mapped reads from multi-mapper-resolved reads, and both are carried through parallel DESeq2 differential expression analyses. Per-position coverage and mismatch profiles are computed using the Sprinzl canonical numbering system, enabling comparison of structural positions across tRNAs of different lengths and identification of modification-induced reverse transcriptase stops. CCA tail integrity is quantified as a proxy for tRNA maturation and aminoacylation status.

All outputs are aggregated into an interactive MultiQC report. Testing was performed on human datasets, but prebuilt databases are hosted through the GDSC for human, fly, and mouse.

## Documentation
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Configuring a GDSC-Clover-Seq pipeline run](docs/configuration.md)  
- [Optional Database Building](docs/database.md)  
- [Submitting the Pipeline](docs/submitting.md)  
- [Understanding the Outputs](docs/Understanding_the_Outputs.md)  
- [Contact](#contact)
- [Citation and Licensing](#citation-and-licensing)


## Installation

To install this code, clone the github repository. This will create a new folder in your working directory called `GDSC-clover-Seq`

```shell

#----- Clone repository and move into it.
git clone https://github.com/Dartmouth-Data-Analytics-Core/GDSC-clover-Seq/
cd GDSC-clover-Seq

```

Comprehensive documentation can be obtained at the links above in the Table of Contents.

## Quick Start

**STEP 1:** Make a folder for environments (this only needs to be done once)

```bash
mkdir clover-seq-envs/
```

**STEP 2:** Move your data into the pipeline folder

```bash
mv <path/to/data> GDSC-clover-Seq
```

**STEP 3**: Configure Run

- Ensure Sample_list_SE.txt is correct
- Edit `prebuilt_config/hg38` to use correct reference level
- Edit `job.script.sh` to point to correct config (one of `hg38`, `dm6`, or `mm10`
- Edit `job.script.sh` to point to conda prefix (the path you created in step 1)
- Edit `job.script.sh` to use your Dartmouth email in the SBATCH header

**STEP 4:** Submit the job

```bash
sbatch job.script.sh
```


## Contact
Please address questions to **DataAnalyticsCore@groups.dartmouth.edu** or submit an issue in the GitHub repository.

## Citation and Licensing

This codebase is adapted from the [original tRAX tool](https://github.com/UCSC-LoweLab/tRAX), licensed under GPL v3.0.

Modifications to modernize and adapt this pipeline for use via Snakemake were made by Mike Martinez, Dartmouth Genomic Data Science Core

**Citation:** [Holmes AD, Howard JM, Chan PP, and Lowe TM.](https://www.biorxiv.org/content/10.1101/2022.07.02.498565v1)

**This pipeline was created with funds from the COBRE grant 1P20GM130454. If you use the pipeline in your own work, please acknowledge the pipeline by citing the grant number in your manuscript.**

**Protocols.io DOI:** dx.doi.org/10.17504/protocols.io.n2bvjewz5gk5/v1

