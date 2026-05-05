# Database Build Module

> ⚠️ **These rules are optional and toggleable via prebuilt_configs**   
> Pre-built databases for hg38, mm10, and dm6 are available on Dartmouth HPC at `/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/genomic_references/tRAX_databases/`.
> We recommend using these (as they are specified in the config files) unless you need a custom organism or updated annotation. Default behavior of the > pipeline is to run without building a database.

---

## What the database contains

The tRNA reference database is a combined tRNA-genome reference that enables alignment of reads to both mature tRNA sequences and their genomic loci simultaneously.

**Mature tRNA sequences** are constructed from [gtRNAdb](https://gtrnadb.ucsc.edu) predictions with:

- Introns removed
- 3′ CCA tail appended (not genomically encoded)
- 5′ G added to His-tRNA (post-transcriptional modification)
- 20 bp padding on both ends to accommodate partial reads

**Genomic loci** are full native genome sequences from the appropriate assembly, allowing reads to be assigned to pre-tRNA loci for accurate isotype counting.

Multiple-sequence alignments (MSAs) are generated with [cmalign](https://manpages.ubuntu.com/manpages/trusty/man1/cmalign.1.html) (part of [Infernal](https://github.com/brendanf/inferrnal)) using eukaryotic tRNA covariance models from tRNAscan-SE. These MSAs enable downstream calculation of **Sprinzl positions** — a canonical numbering system that maps structurally equivalent bases across tRNAs of different lengths — used by the coverage and mismatch modules.

A [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) index of the combined tRNA-genome is built as the final step.

---

## Rules

| Rule | Purpose | Conda environment |
|------|---------|-------------------|
| `generate_gtRNA_db` | Download Ensembl GTF, gtRNAdb tarball, genome FASTA, and build the tRNA database | `clover-seq` |
| `concat_tRNAs` | Concatenate mature tRNA FASTA and genome FASTA into the combined tRNA-genome FASTA | `clover-seq` |
| `tRNA_bt2_index` | Build Bowtie2 index of the tRNA-genome | `clover-bowtie2` |

---

## Running the database build

### 1. Choose a prebuilt config

Organism-specific configs in `prebuilt_configs/` point to the correct URLs on gtRNAdb and Ensembl:

```
prebuilt_configs/
├── hg38_config.yaml   # Homo sapiens GRCh38
├── mm10_config.yaml   # Mus musculus GRCm38
└── dm6_config.yaml    # Drosophila melanogaster dm6
```

### 2. Toggle the build rules in your config

Open the `prebuilt_config/` corresponding to your organism and set the `build_database:` parameter to `true`

## Database directory structure

After a successful build, the database directory contains:

```
{genome}_db/
├── {genome}-tRNAs.fa           Mature tRNA sequences (with CCA tails)
├── {genome}-tRNAs-loci.fa      Pre-tRNA genomic locus sequences
├── {genome}-tRNAs.bed          Mature tRNA BED annotation
├── {genome}-tRNAs-loci.bed     Pre-tRNA genomic loci BED annotation
├── {genome}-tRNAs.stk          MSA for Sprinzl position mapping (mature)
├── {genome}-tRNAs-loci.stk     MSA for Sprinzl position mapping (loci)
├── {genome}.gtf                Ensembl GTF for smRNA biotype assignment
├── tRNAgenome.fa               Combined tRNA + genome FASTA
└── bt2_index/
    └── db-tRNAgenome.*         Bowtie2 index files
```

Set `trna_db` and `bt2_index` in your preprocessing config to point to this directory.
