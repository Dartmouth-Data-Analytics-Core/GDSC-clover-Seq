# Understanding the Outputs

Clover-Seq v3.0 is a Snakemake pipeline adapted from tRAX (doi: 10.1101/2022.07.02.498565) for tRNA-seq analysis on Dartmouth HPC. This document describes every output file, organized by analysis stage.

```
01_trimming/              Adapter-trimmed reads and trimming logs
02_tRNA_alignment/        Aligned BAMs, mapping logs, alignment statistics
03_Raw_Quant/             Raw read counts — tRNA-specific and all small RNA biotypes
04_Expression/            DESeq2 size factors, normalized counts, DE results, R objects
05_Mismatches/            Per-position mismatch profiles and heatmaps
06_Coverages/             Per-position read coverage across tRNAs and pre-tRNA loci
07_Plots/                 Figures — PCA, isoacceptor counts, CCA ends, coverage
08_QC/                    MultiQC custom content and aggregated HTML report
```

---

## Read Trimming

Adapters are removed with cutadapt. Reads shorter than the configured `minlength` are discarded. Only the trimmed reads proceed to alignment.

| Output | Rule | Contents |
|--------|------|----------|
| `01_trimming/{sample}.R1.trim.fastq.gz` | `trimming` | Adapter-trimmed reads for each sample |
| `01_trimming/logs/{sample}.cutadapt.log` | `trimming` | Trimming stats — adapter content, reads passing filter, length distribution; parsed by MultiQC |

---

## tRNA Alignment

Reads are aligned to a combined reference containing both mature tRNA sequences (with 3′ CCA tails) and the full genomic sequence, using Bowtie2 in local mode (`--very-sensitive -k 100`). Multi-mappers are resolved by `choosemappings.py`, which retains only the highest-scoring alignment(s) and classifies each read as tRNA or non-tRNA. The raw bowtie2 BAM is a temporary file and is deleted after `choosemappings.py` completes.

| Output | Rule | Contents |
|--------|------|----------|
| `02_tRNA_alignment/{sample}.srt.bam` | `tRNA_choosemappings` | Sorted, indexed BAM with best-hit multi-mapper selection applied; input to all downstream counting and coverage rules |
| `02_tRNA_alignment/{sample}.raw.bam` | `tRNA_bowtie2` | Raw bowtie2 output before multi-mapper filtering; temporary file, deleted on completion |
| `02_tRNA_alignment/logs/{sample}.bowtie2.log` | `tRNA_bowtie2` | Bowtie2 alignment summary — overall alignment rate, uniquely/multi-mapped read counts; parsed by MultiQC |
| `02_tRNA_alignment/logs/{sample}.choosemappings.log` | `tRNA_choosemappings` | tRNA read assignment stats — total tRNA reads, multi-transcript/anticodon/amino ambiguity counts, non-tRNA reads, imperfect matches; parsed into the MultiQC report as a table |
| `02_tRNA_alignment/stats/{sample}.srt.bam.idxstats` | `tRNA_map_stats` | Per-reference sequence read counts from `samtools idxstats`; columns: reference name, length, mapped reads, unmapped reads |
| `02_tRNA_alignment/stats/{sample}.srt.bam.flagstat` | `tRNA_map_stats` | Summary alignment flags from `samtools flagstat` — total, mapped, duplicate read counts; parsed by MultiQC |

---

## tRNA Read Counting

All tRNA count files are produced in a single run of `countreads.py`. The pipeline counts reads at multiple levels of resolution: per genomic locus, per mature tRNA gene, and broken out by fragment type. Only features with at least 5 reads in any sample are reported.

The key distinction between count files is what reads they include and how granularly they report them:

- **Unique counts** include only reads that mapped unambiguously to a single tRNA gene (no equal-scoring alternative alignments). These are the most conservative counts and are used as the primary DESeq2 input.
- **Total counts** (isotype file) include all reads assigned to each tRNA — both uniquely mapping and multi-mapping reads that were resolved by `choosemappings.py` to their best-scoring tRNA hit. This is the more sensitive count, equivalent to the primary count used in the original tRAX pipeline.
- **Detailed counts** break every gene into five fragment-type rows. These are not used for DESeq2 directly but are useful for inspecting the composition of reads (e.g., whether signal comes from tRNA halves vs. full-length reads).

| Output | Rule | Contents |
|--------|------|----------|
| `03_Raw_Quant/tRNA_counts/unique_tRNA_counts.txt` | `tRNA_count` | One row per tRNA gene; counts reads that mapped unambiguously to that gene only (YR tag = 1). Primary input to DESeq2 normalization and differential expression. |
| `03_Raw_Quant/tRNA_counts/tRNA_isotype_counts.txt` | `tRNA_count` | One row per mature tRNA gene and per pre-tRNA locus; total reads including resolved multi-mappers. Pre-tRNA loci appear first, then mature tRNAs. Input to the tRNA isotype DESeq2 analysis. |
| `03_Raw_Quant/tRNA_counts/gene_level_counts_detailed.txt` | `tRNA_count` | Five rows per tRNA gene (fragment types: `_wholecounts`, `_fiveprime`, `_threeprime`, `_other`, `_antisense`). Pre-tRNA loci have `_wholeprecounts`, `_partialprecounts`, `_trailercounts`. Use this to inspect whether tRNA signal comes from full-length reads or tRNA-derived fragments (tDRs). |
| `03_Raw_Quant/tRNA_counts/gene_level_counts_collapsed.txt` | `tRNA_count` | Fragment types summed to one row per tRNA gene (all reads, not fragment-specific). Produced by awk collapsing of the detailed file. |
| `03_Raw_Quant/tRNA_counts/tRNA_ends_counts.txt` | `tRNA_count` | CCA tail completeness per tRNA per sample. Rows are `{gene}\t{end_type}` pairs. End types: `CCA` (intact mature 3′ end), `CC`, `C` (partially trimmed), `Trimmed` (3+ nt removed or aminoacylated). |
| `03_Raw_Quant/tRNA_counts/genetype_counts.txt` | `tRNA_count` | Gene-to-biotype mapping table used by DESeq2. Columns: feature name, biotype, chromosome, mean read length. |

---

## Small RNA Counting

`count_all_smRNA.py` classifies every aligned read into a biotype hierarchy (tRNA > pre-tRNA > Ensembl annotation > other) and reports read length distributions, amino acid and anticodon-level tRNA summaries, and per-sample biotype breakdowns. Size factors from DESeq2 are applied for normalization.

| Output | Rule | Contents |
|--------|------|----------|
| `03_Raw_Quant/raw_amino_counts_by_group.txt` | `count_smRNAs` | Normalized tRNA read counts grouped by amino acid (20 standard + suppressor/unknown), one column per replicate group. Reads are only counted if they map unambiguously to a single amino acid family. |
| `03_Raw_Quant/raw_anticodon_counts_by_sample.txt` | `count_smRNAs` | Normalized tRNA counts grouped by anticodon (isoacceptor family), one column per sample. More stringent than amino acid counts — a read must map unambiguously to one anticodon to be counted. |
| `03_Raw_Quant/other_smRNAs/smRNA_raw_counts_by_group.txt` | `count_smRNAs` | Normalized biotype read counts per replicate group. Rows are RNA biotype categories in fixed priority order: tRNA, pre-tRNA, then Ensembl biotypes (miRNA, snoRNA, snRNA, rRNA, etc.), then other. |
| `03_Raw_Quant/other_smRNAs/smRNA_raw_counts_by_sample.txt` | `count_smRNAs` | Same biotype breakdown as above but raw (un-normalized) counts for each individual sample. Use for per-sample QC before normalization. |
| `03_Raw_Quant/other_smRNAs/read_length_distribution.txt` | `count_smRNAs` | Read length counts (0–100 nt) broken down by biotype class (tRNA, pre-tRNA, other) for each sample. Mature tRNA reads typically peak at 30–36 nt; tDRs at 18–22 nt. |
| `03_Raw_Quant/other_smRNAs/subroup_counts.txt` | `count_smRNAs` | Mismatch count distribution — how many reads have 0, 1, 2 … mismatches, split by tRNA vs. non-tRNA. Used to assess alignment quality. (Filename contains a historical typo — `subroup` instead of `subgroup`.) |

---

## Normalization, PCA, and Differential Expression

`clover-seq-DESeq2.R` runs two parallel DESeq2 analyses: one on the gene-level counts (all reads, `gene_level_counts_detailed.txt`) and one on uniquely-mapped tRNA counts (`unique_tRNA_counts.txt`). Both follow the tRAX normalization approach: size factors are estimated separately with `estimateSizeFactors()`, normalization is applied as a manual sweep division (not DESeq2's internal normalized count accessor), and dispersion/DE fitting uses `betaPrior=TRUE`.

Differential expression is run for all pairwise group comparisons (or a supplied comparisons file) and outputs padj tables, log2 fold-change tables, group median counts, combined results tables, and volcano plots — one set per comparison, per analysis.

| Output | Rule | Contents |
|--------|------|----------|
| `04_Expression/gene_level_counts_size_factors.csv` | `normalize_and_PCA` | DESeq2 size factors for gene-level counts. Two-row format: sample names, then size factor values. Consumed by `get_tRNA_coverage` and `get_mismatches` for coverage normalization. |
| `04_Expression/normalized_gene_level_counts.csv` | `normalize_and_PCA` | Sweep-normalized gene-level counts (raw counts divided by size factors). Suitable for visualization; not for direct DE statistics. |
| `04_Expression/gene_level_DESeq2_object.Rds` | `normalize_and_PCA` | Serialized `DESeqDataSet` for gene-level analysis. Load with `readRDS()` to re-run contrasts or extract results without re-running the pipeline. |
| `04_Expression/unique_tRNAs_counts_size_factors.csv` | `normalize_and_PCA` | DESeq2 size factors for uniquely-mapped tRNA counts. Same format as above. |
| `04_Expression/normalized_unique_tRNAs_counts.csv` | `normalize_and_PCA` | Sweep-normalized unique tRNA counts. Used by `plot_counts` for isoacceptor visualizations. |
| `04_Expression/unique_tRNAs_DESeq2_object.Rds` | `normalize_and_PCA` | Serialized `DESeqDataSet` for unique tRNA analysis. |
| `04_Expression/Differential_Expression/{typename}_{comparison}_volcano.pdf` | `normalize_and_PCA` | Volcano plot per comparison per analysis type (`gene_level` or `unique_tRNAs`). Features with \|log2FC\| > 1.5 and padj below the top-10 threshold are labeled. |
| `04_Expression/Differential_Expression/{typename}_padjs.txt` | `normalize_and_PCA` | Matrix of adjusted p-values — rows are features, columns are comparisons. |
| `04_Expression/Differential_Expression/{typename}_logvals.txt` | `normalize_and_PCA` | Matrix of log2 fold-changes — same structure as padjs. |
| `04_Expression/Differential_Expression/{typename}_medians.txt` | `normalize_and_PCA` | Median normalized counts per group for each feature. |
| `04_Expression/Differential_Expression/{typename}_combine.txt` | `normalize_and_PCA` | Combined results table: log2FC columns + padj columns + group median columns side-by-side. The primary table for downstream filtering and reporting. |
| `04_Expression/Differential_Expression/{typename}_dispersions.txt` | `normalize_and_PCA` | DESeq2 dispersion estimates per feature. Useful for diagnosing overdispersion. |
| `07_Plots/PCA/gene_level_variance_plot.png` | `normalize_and_PCA` | Scree plot — variance explained by each PC at the gene level. |
| `07_Plots/PCA/gene_level_loadings.csv` | `normalize_and_PCA` | Sample PC coordinates (loadings) for gene-level rlog-transformed data. |
| `07_Plots/PCA/gene_level_PCA.png` | `normalize_and_PCA` | PCA scatter plot of samples colored by group — gene level. Primary QC figure for batch effects and group separation. |
| `07_Plots/PCA/unique_tRNAs_variance_plot.png` | `normalize_and_PCA` | Scree plot for uniquely-mapped tRNA PCA. |
| `07_Plots/PCA/unique_tRNAs_loadings.csv` | `normalize_and_PCA` | Sample PC coordinates for unique tRNA rlog data. |
| `07_Plots/PCA/unique_tRNAs_PCA.png` | `normalize_and_PCA` | PCA scatter plot — unique tRNA level. |
| `07_Plots/PCA/PCA_Analysis_Summary.png` | `normalize_and_PCA` | Combined figure showing gene-level and unique tRNA PCAs side-by-side. |

---

## Mismatch Analysis

`getgenomicmismatches.py` computes per-position mismatch and base-composition profiles across all mature tRNAs. Positions with mismatch fraction > 5% in any sample are retained. These profiles reflect tRNA modification sites (e.g., m1A at position 58, which causes reverse transcriptase stops) and can distinguish true modifications from sequencing errors by their consistency across samples.

| Output | Rule | Contents |
|--------|------|----------|
| `05_Mismatches/mature_tRNA_mismatches.txt` | `get_mismatches` | Long-format per-position mismatch table. One row per feature × sample × position. Columns include: normalized coverage, read starts/ends, mismatch and deletion counts, per-nucleotide base counts (A/T/C/G/deletion), and the reference base. |
| `05_Mismatches/mature_tRNA_mismatches.bed` | `get_mismatches` | BED file of high-mismatch positions (fraction > 5%). One interval per nucleotide position. Score field = mismatch fraction × 1000. |
| `05_Mismatches/heatmaps/` | `get_mismatches` | One PNG heatmap per amino acid group (e.g., `Ala_heatmap.png`). Color encodes mismatch frequency across tRNA positions (columns) and samples (rows). Consistent high-mismatch sites across samples indicate RNA modifications. |

---

## Coverage Analysis

`getcoverage.py` computes per-position read depth across mature tRNA sequences and pre-tRNA genomic loci. Positions are mapped to the Sprinzl canonical numbering system via Stockholm alignment, enabling comparison across tRNAs of different lengths — position 34 is always the wobble base, position 58 is the T-loop modification site, etc.

| Output | Rule | Contents |
|--------|------|----------|
| `06_Coverages/mature_tRNA_coverages.txt` | `get_tRNA_coverage` | Long-format coverage table for mature tRNAs. One row per feature × sample × Sprinzl position. Columns: Feature, Sample, position, coverage, readstarts, readends, uniquecoverage, multitrnacoverage, and per-nucleotide base counts. Gap positions in the alignment are excluded by awk post-processing. |
| `06_Coverages/pretRNA_locus_coverages.txt` | `get_tRNA_coverage` | Same format as above but for pre-tRNA genomic loci, using locus-specific Sprinzl-equivalent position numbering. Includes trailer coverage (`_trailercounts`) positions downstream of the mature tRNA end. |

---

## Plots

Visualization rules run after normalization and coverage are complete. Coverage plots are produced for every tRNA gene.

| Output | Rule | Contents |
|--------|------|----------|
| `07_Plots/Grouped_boxplot_norm_unique_tRNAss_by_Sample_and_Anticodon.png` | `plot_counts` | Boxplot of normalized unique tRNA counts, grouped by sample and anticodon. Shows spread of isoacceptor expression within and between conditions. |
| `07_Plots/Isoacceptor_counts_by_sample_normalized.png` | `plot_counts` | Normalized read counts for each isoacceptor family, one bar per sample. Shows per-sample variation in tRNA pool composition. |
| `07_Plots/Isoacceptor_counts_normalized.png` | `plot_counts` | Isoacceptor counts summarized across replicate groups. Primary figure for comparing tRNA pool composition between conditions. |
| `07_Plots/CCA_ends_Relative_Abundances.png` | `plot_counts` | Relative proportions of CCA end types (CCA / CC / C / Trimmed) per sample. A high proportion of `Trimmed` indicates aminoacylation or 3′ end damage. |
| `07_Plots/CCA_ends_normalized_absolute_abundances.png` | `plot_counts` | Normalized absolute counts of each CCA end type. Use when total tRNA pool size changes between conditions and proportions alone are misleading. |
| `07_Plots/smRNA_Relative_Abundances.png` | `plot_counts` | Stacked bar chart of RNA biotype contributions (tRNA, miRNA, snoRNA, rRNA, other) per sample. Primary library composition QC figure. |
| `07_Plots/coverage/{tRNA}_normalized_coverage.png` | `plot_counts` | Per-tRNA Sprinzl-position coverage plot, one PNG per mature tRNA gene, faceted by sample. |

---

## Quality Control

MultiQC custom content files are generated by `generate_mqc_custom_content.py` before being aggregated into the final HTML report. The report consolidates outputs from every stage of the pipeline into a single interactive document.

| Output | Rule | Contents |
|--------|------|----------|
| `08_QC/mqc_custom_content/unique_tRNAs_abundance_mqc.tsv` | `generate_mqc_content` | Top-20 tRNA genes by total expression across all samples; formatted as a MultiQC bargraph with % toggle. |
| `08_QC/mqc_custom_content/cca_tail_status_mqc.tsv` | `generate_mqc_content` | CCA tail completeness summed across all tRNAs per sample; formatted as a MultiQC bargraph. |
| `08_QC/mqc_custom_content/smrna_biotype_mqc.tsv` | `generate_mqc_content` | Small RNA biotype distribution per sample; formatted as a MultiQC bargraph. |
| `08_QC/mqc_custom_content/choosemappings_stats_mqc.tsv` | `generate_mqc_content` | Per-sample tRNA alignment assignment statistics parsed from `choosemappings.py` logs. Columns: Total tRNA Reads, Multi-Transcript Reads, Multi-Anticodon Reads, Multi-Amino Reads, Non-tRNA (unique), Non-tRNA (multi-mapped), Imperfect Matches. Displayed as an interactive table in MultiQC. |
| `08_QC/tRNA_multi_QC_report.html` | `rule all` | Aggregated MultiQC HTML report consolidating cutadapt trimming stats, bowtie2 alignment rates, samtools flagstat/idxstats, tRNA count summaries, and all custom content above. Open in any browser. This is the primary deliverable for QC review and should be inspected before interpreting differential expression results. |
