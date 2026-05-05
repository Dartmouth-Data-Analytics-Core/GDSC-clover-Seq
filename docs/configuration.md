# Configuration

All pipeline parameters are controlled by a YAML config file and a sample metadata sheet. Organism-specific starter configs are provided in `prebuilt_configs/`.

---

## Config File

Pass your appropriate config in `job.script.sh`:

```shell
#----- Specify Config (one of "hg38", "mm10", "dm6", case-sensitive and needs to be in quotes.)
CONFIG="hg38"

```

These correspond to 3 prebuilt configs, located in the `prebuilt_configs` folder. Parameters required to run this pipeline can be adjusted there if necessary. **Required** parameters are listed in the general settings section and should be checked and adjusted *before* submitting the pipeline. At minimum, the following three parameters **must** be specified:

- Sample_txt
- genome
- refLevel

## Prebuilt Software Environments

The `job.script.sh` contains a `conda prefix` argument that points to pre-built environments hosted by the GDSC. It is recommended to point this to a folder you create. The first time you run the pipeline, the environments will build to this location and be accessible every time afterwards.


### General settings

| Parameter | Type | Example | Description |
|-----------|------|---------|-------------|
| `sample_txt` | string | `Sample_list_SE.txt` | Path to sample metadata sheet |
| `layout` | string | `single` | Library layout: `single` or `paired` (paired-end support is under development) |
| `genome` | string | `hg38` | Genome build: one of `hg38`, `mm10`, or `dm6` |
| `refLevel` | string | `shGFP` | Reference group for DESeq2 contrasts — must match a value in the `Group` column of your sample sheet |
| `build_database` | bool | `false` | Set `true` to trigger Module 1 (database build); leave `false` when using pre-built databases |

### tRNA database paths

| Parameter | Type | Example | Description |
|-----------|------|---------|-------------|
| `trna_db` | string | `/dartfs-hpc/.../hg38_db` | Path to the pre-built tRNA database directory |
| `bt2_index` | string | `/dartfs-hpc/.../db-tRNAgenome` | Path + index prefix for the Bowtie2 tRNA-genome index |

### Trimming parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `adapter_1` | string | `GAGATCGGAAGAGCACACGTC` | 3′ adapter sequence for read 1 (Illumina TruSeq) |
| `adapter_2` | string | `GATCGTCGGACTGTAGAACTC` | 3′ adapter for read 2, if paired-end |
| `minlength` | string | `15` | Minimum retained read length after trimming (bp) |

### Alignment parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `maxMaps` | string | `100` | Maximum multi-mapping alignments accepted per read. Higher values capture more tRNA paralogs at the cost of ambiguous assignments |
| `nPenalty` | string | `5` | Bowtie2 penalty for ambiguous (N) bases. Set to 5 to tolerate tRNA modification sites that appear as N in the reference |

### Differential expression parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `log2FC_magnitude_threshold` | number | `1` | Minimum |log2FC| to call a result biologically significant |
| `padj_threshold` | number | `0.05` | Adjusted p-value cutoff for statistical significance |

### Example config (hg38)

```yaml
sample_txt: "Sample_list_SE.txt"
layout: "single"
genome: "hg38"
refLevel: "shGFP"
build_database: false

trna_db: "/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/genomic_references/tRAX_databases/hg38_db"
bt2_index: "/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/genomic_references/tRAX_databases/hg38_db/bt2_index/db-tRNAgenome"

adapter_1: "GAGATCGGAAGAGCACACGTC"
adapter_2: "GATCGTCGGACTGTAGAACTC"
minlength: "15"

maxMaps: "100"
nPenalty: "5"

log2FC_magnitude_threshold: 1
padj_threshold: 0.05
```

---

## Sample Sheet

The sample sheet is a comma-delimited file with **exactly three columns**. The header row is required and the column names must match exactly — the pipeline will fail if they are changed.

| Column | Description |
|--------|-------------|
| `Sample_ID` | Unique identifier for each sample |
| `fastq_1` | Path to the raw `.fastq.gz` file |
| `Group` | Experimental group label used for DESeq2 contrasts |

### Example

```text
Sample_ID,fastq_1,Group
shGFP1,data/shGFP_1_R1.fastq.gz,shGFP
shGFP2,data/shGFP_2_R1.fastq.gz,shGFP
shArg1,data/shArg_1_R1.fastq.gz,shArg
shArg2,data/shArg_2_R1.fastq.gz,shArg
shArg3,data/shArg_3_R1.fastq.gz,shArg
```

> ⚠️ **Header names are required**
> Do **not** rename the header columns (`Sample_ID`, `fastq_1`, `Group`). Only edit the data rows.

> ⚠️ **Group labels and refLevel**
> The `refLevel` in your config must exactly match one of the values in the `Group` column.

---


## Cluster Profile

The bundled SLURM cluster profile is at `cluster_profile/config.yaml`. Default settings, which should not need to be editted:

```yaml
jobs: 10           # maximum parallel jobs
executor: slurm
default-resources:
  - partition=standard
  - cpus=1
  - mem_mb=8000
  - time="04:00:00"
```xs
