# Understanding the Outputs

This document describes every output file produced by the Clover-Seq Snakemake pipeline, grouped by output directory. For each file, the rule that produced it, the file format, and the meaning of its contents are described.

---

## Directory Structure Overview

```
01_trimming/          Adapter-trimmed reads
02_tRNA_alignment/    Aligned BAMs, duplicate marking, alignment stats
03_Raw_Quant/         Raw read counts for tRNA and all small RNA biotypes
04_Expression/        DESeq2-normalized counts and serialized R objects
05_Mismatches/        Per-position mismatch profiles and heatmaps
06_Coverages/         Per-position read coverage across mature tRNAs
07_Plots/             Publication-ready figures
09_QC/                Aggregated MultiQC quality report
```

---

## 01_trimming/

### `01_trimming/{sample}.R1.trim.fastq.gz`
**Rule:** `trimming`

Adapter-trimmed reads for each sample, produced by cutadapt. Reads shorter than the configured minimum length (`minlength`) are discarded. This file is the input to alignment and is not typically used directly for downstream analysis.

**Quality check:** Trimming statistics (adapter content, reads passing filter, length distribution) are in `01_trimming/logs/{sample}.cutadapt.log` and are included in the MultiQC report.

---

## 02_tRNA_alignment/

### `02_tRNA_alignment/{sample}.srt.bam`
**Rule:** `tRNA_align`

Sorted, indexed BAM file for each sample aligned to the tRNA genome database using Bowtie2 (`--local --very-sensitive`). The database contains both mature tRNA chromosomes and genomic pre-tRNA loci, allowing reads to be assigned to either. Multi-mapping is permitted up to the `maxMaps` limit configured in the run config.

---

### `02_tRNA_alignment/duplicates/{sample}.mkdup.bam`
**Rule:** `tRNA_mark_duplicates`

Duplicate-marked BAM produced by Picard MarkDuplicates. Duplicate reads are flagged in the BAM flag field but **not** removed — all downstream counting scripts see every read, including duplicates. This is intentional for Clover-Seq, where optical duplicates are low and PCR amplification bias is tracked separately.

### `02_tRNA_alignment/duplicates/{sample}.mkdup.log.txt`
**Rule:** `tRNA_mark_duplicates`

Picard metrics file reporting the number of reads examined, duplicates identified, and estimated library size. Included in the MultiQC report.

---

### `02_tRNA_alignment/stats/{sample}.mkdup.bam.idxstats`
**Rule:** `tRNA_map_stats`

Tab-delimited file from `samtools idxstats`. Each row is one reference sequence (chromosome) in the database. Columns:

| Column | Description |
|--------|-------------|
| Reference | Chromosome/tRNA name (e.g., `tRNA-Val-AAC-5-1`) |
| Length | Reference sequence length (bp) |
| Mapped | Reads mapped to this reference |
| Unmapped | Reads with this reference but unmapped flag |

Useful for checking which tRNAs captured the most reads at the individual-chromosome level.

### `02_tRNA_alignment/stats/{sample}.mkdup.bam.flagstat`
**Rule:** `tRNA_map_stats`

Summary alignment statistics from `samtools flagstat`: total reads, mapped reads, duplicates, paired-end flags. Included in the MultiQC report.

---

### `02_tRNA_alignment/full_alignment_read_length_distribution.txt`
**Rule:** `read_length_distribution`

Tab-delimited table of read-length counts across all samples. All aligned reads (not classified by RNA type) are counted. Lengths 0–100 are reported for every sample.

| Column | Description |
|--------|-------------|
| Length | Read length (nt) |
| Sample | Sample ID |
| Count | Number of reads at this length in this sample |

This reflects the full complexity of the small RNA landscape before biotype classification. tRNA-derived reads typically peak at 18–22 nt and 30–36 nt (mature tRNA length).

---

## 03_Raw_Quant/tRNA_counts/

These five files are all produced by a single run of `rule tRNA_count` (`countreads.py`).

### `03_Raw_Quant/tRNA_counts/genetype_counts.txt`
**Rule:** `tRNA_count` → `countreads.py` (`printtypefile`)

Tab-delimited gene-to-biotype mapping table. Each row describes one feature row that appears in the count files. This file is used by DESeq2 to assign biotype labels for size-factor calculation and differential expression.

| Column | Description |
|--------|-------------|
| Gene name | Feature identifier (e.g., `tRNA-Val-AAC-5-1_wholecounts`) |
| Biotype | Feature type (e.g., `tRNA`, `tRNA_locus`, `protein_coding`) |
| Chromosome | Chromosome the feature comes from; `tRNA` for mature tRNAs |
| Mean read length | Average length of reads mapping to this feature across all samples |

---

### `03_Raw_Quant/tRNA_counts/gene_level_counts_detailed.txt`
**Rule:** `tRNA_count` → `countreads.py` (`printcountfile`)

Tab-delimited raw count table. Rows are individual tRNA features broken out by fragment type; columns are samples. This is the direct output of `countreads.py` before any collapsing.

**Fragment type suffixes per tRNA gene:**

| Suffix | Description |
|--------|-------------|
| `_wholecounts` | Reads spanning ≥ 90% of the full-length tRNA |
| `_fiveprime` | 5′ tRNA fragment reads |
| `_threeprime` | 3′ tRNA fragment reads |
| `_other` | Reads that do not fit whole/5′/3′ categories |
| `_antisense` | Reads mapping to the antisense strand |

Pre-tRNA loci also appear with `_wholeprecounts`, `_partialprecounts`, and `_trailercounts` suffixes.

---

### `03_Raw_Quant/tRNA_counts/gene_level_counts_collapsed.txt`
**Rule:** `tRNA_count` → awk in shell

The detailed count table collapsed to one row per tRNA gene by summing all fragment types for each gene name. Produced by splitting the gene name on `_` and summing all rows sharing the same base name.

This is the primary input to DESeq2 normalization (`rule normalize_and_PCA`). Use this file for differential expression analysis at the gene level.

---

### `03_Raw_Quant/tRNA_counts/tRNA_isotype_counts.txt`
**Rule:** `tRNA_count` → `countreads.py` (`printtrnacountfile`)

Tab-delimited raw count table with one row per mature tRNA gene or pre-tRNA locus, columns are samples. Only features with at least 5 reads in any sample are included. This is the primary input to isotype-level DESeq2 analysis and is the second input to `rule normalize_and_PCA`.

Unlike the detailed file, each row here is a single tRNA identity — no fragment-type breakdown. Pre-tRNA loci (genomic loci) are listed first, followed by mature tRNAs.

---

### `03_Raw_Quant/tRNA_counts/tRNA_ends_counts.txt`
**Rule:** `tRNA_count` → `countreads.py` (`printtrnaendfile`)

Tab-delimited CCA-end type count table. Rows are `{tRNA_gene}\t{end_type}` pairs; columns are samples. End types reflect the 3′ trailer composition of each read and are used to assess CCA end integrity as a proxy for tRNA maturation and aminoacylation.

| End type | Meaning |
|----------|---------|
| `CCA` | Read ends with full CCA trailer (mature, uncharged) |
| `CC` | Read ends with CC (missing the terminal A) |
| `C` | Read ends with only C |
| `Trimmed` | No CCA trailer (fully trimmed or aminoacylated) |

Only tRNAs with ≥ 5 reads in any sample are included.

---

## 03_Raw_Quant/

### `03_Raw_Quant/raw_amino_counts_by_group.txt`
**Rule:** `count_smRNAs` → `count_all_smRNA.py` (`printaminocounts`)

Tab-delimited table of tRNA counts grouped by amino acid. Rows are the 20 standard amino acids (plus suppressor/unknown classes); columns are replicate groups. Counts are size-factor normalized.

A read is counted toward an amino acid only if it maps uniquely to a single amino acid family — multi-mapping reads that could belong to different amino acid families are excluded. This makes the counts conservative but unambiguous.

---

### `03_Raw_Quant/raw_anticodon_counts_by_sample.txt`
**Rule:** `count_smRNAs` → `count_all_smRNA.py` (`printanticodoncounts`)

Tab-delimited table of tRNA counts grouped by anticodon (isoacceptor family). Rows are anticodons; columns are individual samples (not collapsed by replicate). Counts are size-factor normalized.

A read is counted toward an anticodon only if it maps uniquely to a single anticodon — more stringent than the amino acid filter. This is the recommended table for isoacceptor-level analysis.

---

## 03_Raw_Quant/other_smRNAs/

### `03_Raw_Quant/other_smRNAs/smRNA_raw_counts_by_group.txt`
**Rule:** `count_smRNAs` → `count_all_smRNA.py` (`printtypefile`)

Tab-delimited table of **normalized** biotype read counts per replicate group. Rows are RNA biotype categories; columns are replicate group names. Each read is assigned to exactly one category in priority order (tRNA > pre-tRNA > Ensembl biotype > BED feature > other).

**Row order (fixed):**
- `other` — reads not assigned to any known category
- Ensembl biotypes in priority order: `snoRNA`, `snRNA`, `scaRNA`, `sRNA`, `miRNA`, then other biotypes alphabetically, then `Mt_rRNA`, `Mt_tRNA`, `rRNA`
- `pretRNA_antisense` — antisense reads overlapping pre-tRNA loci
- `pretRNA` — reads overlapping pre-tRNA genomic loci
- `tRNA_antisense` — antisense reads on mature tRNA chromosomes
- `tRNA` — sense reads on mature tRNA chromosomes

---

### `03_Raw_Quant/other_smRNAs/smRNA_raw_counts_by_sample.txt`
**Rule:** `count_smRNAs` → `count_all_smRNA.py` (`printrealcounts`)

Same biotype breakdown as `smRNA_raw_counts_by_group.txt` but with **raw (un-normalized) counts** for each individual sample rather than replicate groups. Use this file when you want to inspect per-sample variation before normalization.

---

### `03_Raw_Quant/other_smRNAs/read_length_distribution.txt`
**Rule:** `count_smRNAs` → `count_all_smRNA.py` (`printlengthfile`)

Tab-delimited read-length distribution broken down by RNA class (tRNA, pre-tRNA, other). Unlike `02_tRNA_alignment/full_alignment_read_length_distribution.txt`, this file shows how the size distribution differs across biotypes.

| Column | Description |
|--------|-------------|
| Length | Read length (nt) |
| Sample | Sample ID |
| other | Reads not classified as tRNA or pre-tRNA |
| trnas | Reads classified as mature tRNA |
| pretrnas | Reads classified as pre-tRNA |

---

### `03_Raw_Quant/other_smRNAs/subroup_counts.txt`
**Rule:** `count_smRNAs` → `count_all_smRNA.py` (`printmismatchcounts`)

Tab-delimited mismatch count distributions, reporting how many reads in each sample have 0, 1, 2, … up to 9 mismatches, split into tRNA and non-tRNA reads. Used to assess alignment quality and mismatch tolerance. Note: the filename contains a typo (`subroup` instead of `subgroup`) inherited from the original pipeline.

| Column | Description |
|--------|-------------|
| count | Number of mismatches (0–9) |
| type | `trna` or `nontrna` |
| {sample} … | Normalized read count at this mismatch level per sample |

---

## 04_Expression/

All files in this directory are produced by `rule normalize_and_PCA` (`clover-seq-DESeq2.R`).

### `04_Expression/gene_level_counts_size_factors.csv`
Size factors computed by DESeq2 for gene-level counts. One row per sample with the sample ID and its size factor. Size factors are used to normalize raw counts to a common library depth. This file is also consumed by `rule get_tRNA_coverage` and `rule get_mismatches` to normalize coverage values.

### `04_Expression/normalized_gene_level_counts.csv`
DESeq2 variance-stabilizing-transformed (VST) counts at the gene level. Rows are tRNA genes; columns are samples. These counts are on a log-like scale, suitable for PCA, heatmaps, and other exploratory visualizations but **not** for differential expression (use the raw counts and the `.Rds` object for that).

### `04_Expression/tRNA_isotype_counts_size_factors.csv`
Size factors computed from the tRNA isotype count matrix. Separate from gene-level size factors because the set of features differs.

### `04_Expression/normalized_tRNA_isotype_counts.csv`
DESeq2 VST-normalized counts at the tRNA isotype level. Used directly by `rule plot_counts` to generate isoacceptor count visualizations.

### `04_Expression/gene_level_DESeq2_object.Rds`
Serialized R `DESeqDataSet` object for gene-level analysis. Load in R with `readRDS()` to perform differential expression, extract normalized counts, or re-run statistical tests with different contrasts without re-running the full pipeline.

### `04_Expression/tRNA_isotype_DESeq2_object.Rds`
Serialized R `DESeqDataSet` object for tRNA isotype-level analysis. Same purpose as above but at isotype resolution.

---

## 04_Expression/ & 07_Plots/PCA/ — PCA outputs

All PCA files are produced by `rule normalize_and_PCA`.

### `07_Plots/PCA/gene_level_variance_plot.png`
Scree plot showing the percentage of variance explained by each principal component at the gene level. Used to determine how many PCs to examine.

### `07_Plots/PCA/gene_level_loadings.csv`
PCA loading scores for gene-level features. Each row is a tRNA gene; columns are PC loadings. Identify which genes drive sample separation along each axis.

### `07_Plots/PCA/gene_level_PCA.png`
PCA scatter plot of samples colored by group at the gene level. The primary QC figure for checking batch effects, outliers, and group separation.

### `07_Plots/PCA/tRNA_isotype_variance_plot.png`
Scree plot at the tRNA isotype level.

### `07_Plots/PCA/tRNA_isotype_loadings.csv`
PCA loadings at the isotype level.

### `07_Plots/PCA/tRNA_isotype_PCA.png`
PCA scatter plot at the isotype level.

### `07_Plots/PCA/PCA_Analysis_Summary.png`
Combined summary figure showing both gene-level and isotype-level PCAs side-by-side for quick comparison.

---

## 05_Mismatches/

### `05_Mismatches/mature_tRNA_mismatches.txt`
**Rule:** `get_mismatches` → `getgenomicmismatches.py` (`transcriptcoverage`)

Tab-delimited per-position mismatch table for all mature tRNAs. Only positions where the mismatch fraction exceeds 5% in at least one sample (or a flanking position does) are reported, keeping file size manageable. Rows are position-sample combinations; the file is long-format.

| Column | Description |
|--------|-------------|
| Feature | tRNA name (e.g., `tRNA-Ala-AGC-1-1`) |
| Sample | Sample ID |
| position | Position index within the tRNA sequence (0-based, or Sprinzl number if `--stkfile` supplied) |
| coverage | Normalized read depth at this position |
| readstarts | Normalized count of reads starting at this position |
| readends | Normalized count of reads ending at this position |
| tRNAreadstotal | Total normalized reads for this tRNA in this sample |
| expreadstotal | Total reads across all samples (used for filtering) |
| actualbase | Reference nucleotide at this position (T substituted for U) |
| mismatchedbases | Normalized count of reads with a mismatch at this position |
| deletedbases | Normalized count of reads with a deletion at this position |
| adenines | Count of reads with A at this position |
| thymines | Count of reads with T at this position |
| cytosines | Count of reads with C at this position |
| guanines | Count of reads with G at this position |
| deletions | Count of reads with a deletion gap at this position |

The per-nucleotide columns (adenines–guanines) report the **observed** base in the read, not whether it is a mismatch — use `actualbase` to determine which calls are mismatches.

---

### `05_Mismatches/mature_tRNA_mismatches.bed`
**Rule:** `get_mismatches` → `getgenomicmismatches.py`

BED file marking high-mismatch positions (mismatch fraction > 5% in any sample). Each interval is a single nucleotide. The BED score field contains `mismatch_fraction × 1000` (integer), allowing BED-compatible tools to filter or visualize by mismatch severity.

Format: `chrom  start  end  name  score  strand`
where `name` is `{tRNA_name}_{position}pos`.

---

### `05_Mismatches/heatmaps/`
**Rule:** `get_mismatches` → `clover-seq-heatmaps.R`

Directory of PNG heatmap images, one per amino acid group (e.g., `Ala_heatmap.png`, `Val_heatmap.png`). Each heatmap shows mismatch frequency as a color gradient across tRNA positions (columns) and samples (rows), making it easy to identify modification sites or damage patterns that are consistent across samples or specific to a condition.

---

## 06_Coverages/

### `06_Coverages/mature_tRNA_coverages.txt`
**Rule:** `get_tRNA_coverage` → `getcoverage.py` + awk filtering

Tab-delimited per-position coverage table for all mature tRNAs, with alignment-gap columns (positions labelled `.gap` in the Sprinzl numbering) removed by the awk post-processing step. Long-format: each row is one position-sample combination.

| Column | Description |
|--------|-------------|
| Feature | tRNA name |
| Sample | Sample ID (or replicate group if `--combinereps` used) |
| position | Sprinzl position number (canonical tRNA numbering) |
| coverage | Normalized read depth at this position |
| readstarts | Normalized count of reads starting here |
| readends | Normalized count of reads ending here |
| uniquecoverage | Coverage from reads uniquely mapping to one tRNA |
| multitrnacoverage | Coverage from reads mapping to multiple tRNAs (same anticodon) |
| multianticodoncoverage | Coverage from reads mapping across anticodon families |
| multiaminocoverage | Coverage from reads mapping across amino acid families |
| tRNAreadstotal | Total reads for this tRNA in this sample |
| actualbase | Reference nucleotide at this position |
| mismatchedbases | Normalized mismatch count |
| deletedbases | Normalized deletion count |
| adenines | Reads with A at this position |
| thymines | Reads with T |
| cytosines | Reads with C |
| guanines | Reads with G |
| deletions | Reads with a deletion at this position |

The Sprinzl position system enables comparison across tRNAs of different lengths by mapping equivalent structural positions (e.g., position 34 is always the wobble base of the anticodon).

---

## 07_Plots/

All plots are produced by `rule plot_counts` (`clover-seq-plot-counts.R` and `clover-seq-plot-coverages.R`).

### `07_Plots/Grouped_boxplot_norm_tRNA_isotypes_by_Sample_and_Anticodon.png`
Boxplot of normalized tRNA isotype counts grouped by both sample and anticodon. Useful for visualizing the spread of isoacceptor expression within and between conditions.

### `07_Plots/Isoacceptor_counts_by_sample_normalized.png`
Bar or strip chart showing normalized read counts for each isoacceptor family, faceted or colored by sample. Shows sample-level variation in the tRNA pool composition.

### `07_Plots/Isoacceptor_counts_normalized.png`
Isoacceptor counts averaged or summarized across replicate groups. The primary figure for comparing the tRNA pool composition between conditions.

### `07_Plots/CCA_ends_Relative_Abundances.png`
Relative proportions of CCA end types (CCA / CC / C / Trimmed) per sample or group. A high proportion of `Trimmed` reads can indicate aminoacylation activity or tRNA 3′ end damage. Derived from `03_Raw_Quant/tRNA_counts/tRNA_ends_counts.txt`.

### `07_Plots/CCA_ends_normalized_absolute_abundances.png`
Absolute (but size-factor-normalized) counts of each CCA end type rather than proportions. Useful when the total tRNA pool size also changes between conditions.

### `07_Plots/smRNA_Relative_Abundances.png`
Stacked bar chart showing the relative contribution of each RNA biotype (tRNA, miRNA, snoRNA, rRNA, other, etc.) to the total aligned read count per sample. Derived from `03_Raw_Quant/other_smRNAs/smRNA_raw_counts_by_sample.txt`. A primary QC figure for assessing library composition.

---

## 09_QC/

### `09_QC/tRNA_multi_QC_report.html`
**Rule:** `rule all` → MultiQC

Aggregated HTML quality report consolidating:
- Cutadapt trimming statistics (adapter content, reads passing filter)
- Bowtie2 alignment rates
- Picard duplicate metrics
- SAMtools flagstat and idxstats summaries
- Count table summaries from `03_Raw_Quant/tRNA_counts/`

Open in any web browser. This is the primary deliverable for QC review and should be inspected before proceeding to differential expression analysis.

> **Note:** The intermediate file `03_Raw_Quant/tRNA_counts/unique_tRNA_counts.txt` is removed by this rule after MultiQC finishes.
