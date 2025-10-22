#!/usr/bin/env python3

import pysam
import sys
import argparse
from trnasequtils import readrnastk  # Assuming this exists in your environment

# Hardcoded Sprinzl positions for eukaryotic tRNAs
EUK_POSITIONS = [
    -1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,
    '17a',18,19,20,'20a','20b',21,22,23,24,25,26,27,28,
    29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,
    'e1','e2','e3','e4','e5','e6','e7','e8','e9','e10','e11',
    'e12','e13','e14','e15','e16','e17','e18','e19',46,47,48,
    49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,
    67,68,69,70,71,72,73,74,75,76
]

def parse_bed(bed_path):
    with open(bed_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split('\t')
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            name = fields[3] if len(fields) > 3 else f"{chrom}:{start}-{end}"
            yield {"chrom": chrom, "start": start, "end": end, "name": name}

def map_ref_to_alignment(ref_seq, aligned_seq):
    """Map ref base idx to alignment column idx; skips gaps in aligned seq."""
    mapping = []
    ref_pos = 0
    for aln_pos, base in enumerate(aligned_seq):
        if base != '-':
            mapping.append(aln_pos)
            ref_pos += 1
    return mapping

def count_coverage_mismatches(bam_path, fasta_path, bed_path, stk_path,
                              max_mismatches=None, min_coverage=0):
    bam = pysam.AlignmentFile(bam_path, "rb")
    fasta = pysam.FastaFile(fasta_path)
    trnastk = list(readrnastk(open(stk_path)))[0]  # parse first alignment block

    print("\t".join(["Feature", "SprinzlPosition", "ReferenceBase", "Coverage", "Mismatches"]))

    for feature in parse_bed(bed_path):
        chrom, start, end, name = feature["chrom"], feature["start"], feature["end"], feature["name"]
        feature_len = end - start

        ref_seq = fasta.fetch(chrom, start, end).upper()
        aligned_seq = trnastk.aligns.get(name)
        if aligned_seq is None:
            print(f"Warning: No alignment sequence for {name}", file=sys.stderr)
            continue

        ref_to_aln = map_ref_to_alignment(ref_seq, aligned_seq)
        aln_len = len(aligned_seq)
        aln_to_ref = [None] * aln_len
        for ref_idx, aln_idx in enumerate(ref_to_aln):
            aln_to_ref[aln_idx] = ref_idx

        coverage = [0] * feature_len
        mismatches = [0] * feature_len

        for read in bam.fetch(chrom, start, end):
            if read.is_unmapped:
                continue
            try:
                nm = read.get_tag("NM")
                if max_mismatches is not None and nm > max_mismatches:
                    continue
            except KeyError:
                pass

            aligned_pairs = read.get_aligned_pairs(matches_only=False)
            for read_pos, ref_pos in aligned_pairs:
                if read_pos is None or ref_pos is None:
                    continue
                if ref_pos < start or ref_pos >= end:
                    continue
                idx = ref_pos - start
                coverage[idx] += 1
                read_base = read.query_sequence[read_pos].upper()
                ref_base = ref_seq[idx]
                if read_base != ref_base:
                    mismatches[idx] += 1

        if sum(coverage) < min_coverage:
            continue

        for aln_pos in range(aln_len):
            sprinzl_pos = EUK_POSITIONS[aln_pos] if aln_pos < len(EUK_POSITIONS) else str(aln_pos + 1)
            if sprinzl_pos is None or "gap" in str(sprinzl_pos).lower():
                continue

            ref_pos = aln_to_ref[aln_pos]
            if ref_pos is None or ref_pos >= feature_len:
                cov = 0
                mism = 0
                ref_base = '-'
            else:
                cov = coverage[ref_pos]
                mism = mismatches[ref_pos]
                ref_base = ref_seq[ref_pos]

            print(f"{name}\t{sprinzl_pos}\t{ref_base}\t{cov}\t{mism}")

def main():
    parser = argparse.ArgumentParser(description="Count tRNA coverage and mismatches with Sprinzl numbering")
    parser.add_argument("--bam", required=True, help="BAM alignment file")
    parser.add_argument("--fasta", required=True, help="Reference FASTA file")
    parser.add_argument("--bed", required=True, help="BED file with tRNA regions")
    parser.add_argument("--stk", required=True, help="Stockholm alignment file")
    parser.add_argument("--maxmismatches", type=int, default=None, help="Max mismatches per read")
    parser.add_argument("--mincoverage", type=int, default=0, help="Minimum coverage to report")

    args = parser.parse_args()

    count_coverage_mismatches(
        args.bam,
        args.fasta,
        args.bed,
        args.stk,
        max_mismatches=args.maxmismatches,
        min_coverage=args.mincoverage,
    )

if __name__ == '__main__':
    main()

# To run...
    #python code/mismatches.py --bam 02_tRNA_alignment/IB1.mkdup.bam --fasta /dartfs-hpc/rc/lab/G/GMBSR_bioinfo/genomic_references/tRAX_databases/hg38_db/db-maturetRNAs.fa --bed /dartfs-hpc/rc/lab/G/GMBSR_bioinfo/genomic_references/tRAX_databases/hg38_db/db-maturetRNAs.bed --stk /dartfs-hpc/rc/lab/G/GMBSR_bioinfo/genomic_references/tRAX_databases/hg38_db/db-trnaalign.stk >TEST.out
