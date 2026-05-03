#!/usr/bin/env python3
"""
Generate MultiQC custom content TSV files for Clover-Seq.

Outputs three bargraph sections to 08_QC/mqc_custom_content/:
  - trna_isotype_abundance_mqc.tsv  : top-20 tRNA isotypes per sample
  - cca_tail_status_mqc.tsv         : CCA tail completeness per sample
  - smrna_biotype_mqc.tsv           : small RNA biotype distribution per sample
"""

import os
import argparse
import pandas as pd


def write_mqc_bargraph(filepath, section_id, section_name, description,
                       plot_id, title, ylab, df):
    """Write a plain TSV; plot config lives in multiqc_config.yaml."""
    df.to_csv(filepath, sep="\t")


def generate_isotype_abundance(counts_file, output_file):
    """Top-20 tRNA genes by total abundance, one row per sample."""
    df = pd.read_table(counts_file, index_col=0)
    # Keep only canonical tRNA entries (exclude rRNA, SNORD, etc.)
    df = df[df.index.str.startswith("tRNA-")]
    top20 = df.sum(axis=1).nlargest(20).index
    # Transpose: rows = samples, cols = gene names
    out = df.loc[top20].T
    write_mqc_bargraph(
        output_file,
        section_id="trna_isotype_abundance",
        section_name="tRNA Gene Abundances (Top 20)",
        description=(
            "Raw read counts for the 20 most abundant tRNA genes per sample. "
            "Toggle to % to compare relative proportions between samples."
        ),
        plot_id="trna_isotype_abundance_plot",
        title="tRNA: Top 20 Gene Abundances",
        ylab="Read Counts",
        df=out,
    )


def generate_cca_tail_status(ends_file, output_file):
    """CCA tail completeness summed across all tRNAs, one row per sample."""
    df = pd.read_table(ends_file, index_col=0)
    # columns: [end, sample1, sample2, ...]
    sample_cols = [c for c in df.columns if c != "end"]
    grouped = df.groupby("end")[sample_cols].sum().T
    # Enforce display order: most-complete first
    col_order = [c for c in ["CCA", "CC", "C", "Trimmed"] if c in grouped.columns]
    out = grouped[col_order]
    write_mqc_bargraph(
        output_file,
        section_id="cca_tail_status",
        section_name="CCA Tail Completeness",
        description=(
            "3&prime; CCA tail completeness summed across all tRNA loci per sample. "
            "<b>CCA</b> = intact 3&prime; end; <b>CC / C</b> = partially trimmed; "
            "<b>Trimmed</b> = 3+ nucleotides removed."
        ),
        plot_id="cca_tail_status_plot",
        title="tRNA: CCA Tail Completeness",
        ylab="Read Counts",
        df=out,
    )


def generate_biotype_distribution(biotype_file, output_file):
    """Small RNA biotype read counts, one row per sample."""
    df = pd.read_table(biotype_file, index_col=0)
    # rows = biotypes, cols = samples → transpose for MultiQC
    out = df.T
    write_mqc_bargraph(
        output_file,
        section_id="smrna_biotype",
        section_name="Small RNA Biotype Distribution",
        description=(
            "Distribution of aligned reads across small RNA biotypes per sample. "
            "Switch to % view to compare relative proportions."
        ),
        plot_id="smrna_biotype_plot",
        title="Small RNA Biotype Distribution",
        ylab="Read Counts",
        df=out,
    )


def main():
    parser = argparse.ArgumentParser(
        description="Generate MultiQC custom content for Clover-Seq"
    )
    parser.add_argument(
        "--gene-level-counts", required=True,
        help="03_Raw_Quant/tRNA_counts/tRNA_isotype_counts.txt"
    )
    parser.add_argument(
        "--ends-counts", required=True,
        help="03_Raw_Quant/tRNA_counts/tRNA_ends_counts.txt"
    )
    parser.add_argument(
        "--biotype-counts", required=True,
        help="03_Raw_Quant/other_smRNAs/smRNA_raw_counts_by_sample.txt"
    )
    parser.add_argument(
        "--output-dir", required=True,
        help="Output directory (08_QC/mqc_custom_content)"
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    generate_isotype_abundance(
        args.gene_level_counts,
        os.path.join(args.output_dir, "trna_isotype_abundance_mqc.tsv"),
    )
    print("Generated: trna_isotype_abundance_mqc.tsv")

    generate_cca_tail_status(
        args.ends_counts,
        os.path.join(args.output_dir, "cca_tail_status_mqc.tsv"),
    )
    print("Generated: cca_tail_status_mqc.tsv")

    generate_biotype_distribution(
        args.biotype_counts,
        os.path.join(args.output_dir, "smrna_biotype_mqc.tsv"),
    )
    print("Generated: smrna_biotype_mqc.tsv")


if __name__ == "__main__":
    main()
