
# Dartmouth GDSC Clover-Seq Pipeline
![Version](https://img.shields.io/badge/version-2.0-blue)

This pipeline provides preprocessing and quality control of tRNA sequencing data. This pipeline has been built and tested using human, mouse, and fly data sets. Major steps of this pipeline include: 

- Trimming of adapters using [*Cutadapt*](https://cutadapt.readthedocs.io/en/stable/)
- Alignment to custom tRNA databases (mature and pre-tRNA loci) using [*bowtie2*](https://github.com/BenLangmead/bowtie2)
- Quantification of mature tRNA and tRNA isotype sequences using custom code
- Mismatch and coverage analyses
- Quality Control and differential expression analyses using [DESeq2](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)

## Documentation
- [Installation](#installation)
- [Comprehensive Documentation](#comprehensive-documentation)
- [Contact](#contact)
- [Citation and Licensing](#citation-and-licensing)


## Installation

To install this code, clone the github repository. This will create a new folder in your working directory called `GDSC-clover-Seq`

```shell

#----- Clone repository and move into it.
git clone https://github.com/Dartmouth-Data-Analytics-Core/GDSC-clover-Seq/
cd GDSC-clover-Seq

```

Several [conda environments](https://anaconda.org/anaconda/conda) are required to run this code successfully. For your convenience, these conda environments have been prebuilt and are hosted publically by the Dartmouth Genomic Data Science Core through the conda-prefix argument in the `job.script.sh`.

## Comprehensive Documentation

For extended information on how to configure a job and submit it on the [Dartmouth high-performance compute cluster, Discovery](https://rc.dartmouth.edu/hpc/discovery-overview/) see the links below:

[Configuring a GDSC-Clover-Seq pipeline run](docs/configuration.md)  
[Optional Database Building](docs/database.md)  
[Submitting the Pipeline](docs/submitting.md)  
[Understanding the Outputs](docs/Understanding_the_Outputs.md)  

## Contact
Please address questions to **DataAnalyticsCore@groups.dartmouth.edu** or submit an issue in the GitHub repository.

## Citation and Licensing

This codebase is adapted from the [original tRAX tool](https://github.com/UCSC-LoweLab/tRAX), licensed under GPL v3.0.

Modifications to modernize and adapt this pipeline for use via Snakemake were made by Mike Martinez, Dartmouth Genomic Data Science Core

**Citation:** [Holmes AD, Howard JM, Chan PP, and Lowe TM.](https://www.biorxiv.org/content/10.1101/2022.07.02.498565v1)

**This pipeline was created with funds from the COBRE grant 1P20GM130454. If you use the pipeline in your own work, please acknowledge the pipeline by citing the grant number in your manuscript.**

**Protocols.io DOI:** dx.doi.org/10.17504/protocols.io.n2bvjewz5gk5/v1

