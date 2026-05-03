# Installation

## Requirements

- Linux or macOS (HPC cluster recommended for large datasets)
- [Conda](https://docs.conda.io/en/latest/) or [Mamba](https://mamba.readthedocs.io/)
- [Snakemake ≥ 7.18](https://snakemake.readthedocs.io/)
- SLURM (for cluster execution via the bundled cluster profile)

---

## 1. Clone the Repository

```shell
git clone https://github.com/Dartmouth-Data-Analytics-Core/GDSC-clover-Seq/
cd GDSC-clover-Seq
```

---

## 2. Set Up Conda Environments

Three Conda environments are required:

| Environment | Used by |
|-------------|---------|
| `clover-seq` | Most rules (trimming, counting, normalization, QC) |
| `clover-bowtie2` | Alignment (`tRNA_align`, `tRNA_bt2_index`) |
| `Picard` | Duplicate marking (`tRNA_mark_duplicates`) |

### Option A — Use pre-built environments (Dartmouth HPC only)

Pre-built environments are hosted on Discovery at:

```
/dartfs/rc/nosnapshots/G/GMBSR_refs/envs/GDSC-Clover-Seq
```

Pass this path to Snakemake via `--conda-prefix`:

```shell
snakemake ... \
    --use-conda \
    --conda-frontend conda \
    --conda-prefix /dartfs/rc/nosnapshots/G/GMBSR_refs/envs/GDSC-Clover-Seq
```

No additional setup is needed — Snakemake will activate the correct environment for each rule.

### Option B — Build environments from YAML

Environment definitions are in `env_config/`. Build each with:

```shell
conda env create -f env_config/clover-seq.yaml
conda env create -f env_config/clover-bowtie2.yaml
conda env create -f env_config/Picard.yaml
```

Key dependencies in `clover-seq`:

- Python ≥ 3.10, pysam, numpy
- samtools, bowtie2, cutadapt, infernal
- R ≥ 4.4: DESeq2, ggplot2, pheatmap, multiQC

---

## 3. tRNA Reference Databases

The pipeline requires a pre-built tRNA database directory containing:

- Mature tRNA FASTA + Bowtie2 index
- Genomic tRNA loci BED file
- Ensembl GTF
- Sprinzl alignment (`.stk`) files for coverage and mismatch analysis

### Pre-built databases (Dartmouth HPC only)

```
/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/genomic_references/tRAX_databases/
├── hg38_db/
├── mm10_db/
└── dm6_db/
```

Set `trna_db` and `bt2_index` in your config file to the appropriate path (see [Configuration](configuration.md)).

### Build your own database

If your organism is not listed above, or you need a custom reference, use Module 1 (see [Database Build](database.md)).

---

## 4. Verify Installation

After cloning and setting up environments, do a dry run to verify the workflow parses correctly:

```shell
snakemake -s workflows/module-2-preprocess.smk \
    --configfile prebuilt_configs/hg38_config.yaml \
    --use-conda \
    --conda-prefix /dartfs/rc/nosnapshots/G/GMBSR_refs/envs/GDSC-Clover-Seq \
    --profile cluster_profile \
    --dry-run
```

A successful dry run prints the list of jobs that would run without executing anything.
