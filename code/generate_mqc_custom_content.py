#!/usr/bin/env python3
"""
Generate MultiQC custom content TSV files for Clover-Seq.

Outputs to 08_QC/mqc_custom_content/:
  - unique_tRNAs_abundance_mqc.tsv   : top-20 tRNA isotypes per sample (bargraph)
  - cca_tail_status_mqc.tsv          : CCA tail completeness per sample (bargraph)
  - smrna_biotype_mqc.tsv            : small RNA biotype distribution per sample (bargraph)
  - choosemappings_stats_mqc.tsv     : per-sample tRNA mapping assignment stats (table)
"""

import os
import glob
import argparse
import pandas as pd


def write_mqc_bargraph(filepath, section_id, section_name, description,
                       plot_id, title, ylab, df):
    """Write a plain TSV; plot config lives in multiqc_config.yaml."""
    df.to_csv(filepath, sep="\t")


def generate_isotype_abundance(counts_file, output_file):
    """Top-20 tRNA genes by total abundance, one row per sample."""
    df = pd.read_table(counts_file, index_col=0)
    df = df[df.index.str.startswith("tRNA-")]
    top20 = df.sum(axis=1).nlargest(20).index
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
    sample_cols = [c for c in df.columns if c != "end"]
    grouped = df.groupby("end")[sample_cols].sum().T
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


def generate_choosemappings_table(log_dir, output_file):
    """Parse per-sample *.choosemappings.log files into a MultiQC table."""
    field_map = {
        "Total tRNA Reads":                    "Total tRNA Reads",
        "tRNA Reads with multiple transcripts": "Multi-Transcript Reads",
        "tRNA Reads with multiple anticodons":  "Multi-Anticodon Reads",
        "tRNA Reads with multiple aminos":      "Multi-Amino Reads",
        "Single mapped non-tRNAs":              "Non-tRNA (unique)",
        "Multiply mapped non-tRNAs":            "Non-tRNA (multi-mapped)",
        "Imperfect matches":                    "Imperfect Matches",
    }
    col_order = list(field_map.values())

    rows = []
    for log_file in sorted(glob.glob(os.path.join(log_dir, "*.choosemappings.log"))):
        sample = os.path.basename(log_file).replace(".choosemappings.log", "")
        row = {"Sample": sample}
        with open(log_file) as fh:
            for line in fh:
                line = line.strip()
                if ":" not in line:
                    continue
                key, val = line.split(":", 1)
                key, val = key.strip(), val.strip()
                if key not in field_map:
                    continue
                col = field_map[key]
                # "Imperfect matches" value is "X/Y" — keep the numerator
                row[col] = int(val.split("/")[0]) if key == "Imperfect matches" else int(val)
        rows.append(row)

    df = pd.DataFrame(rows).set_index("Sample")[col_order]
    df.to_csv(output_file, sep="\t")


def main():
    parser = argparse.ArgumentParser(
        description="Generate MultiQC custom content for Clover-Seq"
    )
    parser.add_argument(
        "--gene-level-counts", required=True,
        help="03_Raw_Quant/tRNA_counts/unique_tRNA_counts.txt"
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
        "--log-dir", required=True,
        help="Directory containing per-sample *.choosemappings.log files (02_tRNA_alignment)"
    )
    parser.add_argument(
        "--output-dir", required=True,
        help="Output directory (08_QC/mqc_custom_content)"
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    generate_isotype_abundance(
        args.gene_level_counts,
        os.path.join(args.output_dir, "unique_tRNAs_abundance_mqc.tsv"),
    )
    print("Generated: unique_tRNAs_abundance_mqc.tsv")

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

    generate_choosemappings_table(
        args.log_dir,
        os.path.join(args.output_dir, "choosemappings_stats_mqc.tsv"),
    )
    print("Generated: choosemappings_stats_mqc.tsv")


if __name__ == "__main__":
    main()
