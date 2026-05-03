# Table of Contents
- [Introduction](#introduction)
- [Installation](#installation)
- [Comprehensive Documentation](#comprehensive-documentation)
- [Contact](#contact)
- [Citation and Licensing](#citation-and-licensing)

## Introduction

Modular Snakemake workflows for the comprehensive analyses of mature tRNAs and other small RNAs (smRNAs) from high-throughput sequencing data. 

## Installation

To install this code, clone the github repository. This will create a new folder in your working directory called `GDSC-clover-Seq`

```shell

#----- Clone repository and move into it.
git clone https://github.com/Dartmouth-Data-Analytics-Core/GDSC-clover-Seq/
cd GDSC-clover-Seq

```

Several [conda environments](https://anaconda.org/anaconda/conda) are required to run this code successfully. For your convenience, these conda environments have been prebuilt and are hosted publically by the Dartmouth Genomic Data Science Core through the conda-prefix argument in the `job.script.sh`.

## Comprehensive Documentation

[Configuring a GDSC-Clover-Seq pipeline run](docs/configuration.md)
[Optional Database Building](docs/database.md)
[Submitting the Pipeline](docs/submitting.md)
[Understanding the Outputs](docs/Understanding_the_Outputs.md)

## Contact
Please address questions to **DataAnalyticsCore@groups.dartmouth.edu** or submit an issue in the GitHub repository.

## Citation and Licensing

**This codebase is adapted from the [original tRAX tool](https://github.com/UCSC-LoweLab/tRAX), licensed under GPL v3.0., modified by Mike Martinez, Dartmouth Genomic Data Science Core** 

**Citation:** [Holmes AD, Howard JM, Chan PP, and Lowe TM.](https://www.biorxiv.org/content/10.1101/2022.07.02.498565v1)

**This pipeline was created with funds from the COBRE grant 1P20GM130454. If you use the pipeline in your own work, please acknowledge the pipeline by citing the grant number in your manuscript.**

**Protocols.io DOI:** dx.doi.org/10.17504/protocols.io.n2bvjewz5gk5/v1

