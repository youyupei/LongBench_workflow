#!/usr/bin/env python3
# Usage: python3 find_longest_reads_in_bam.py <bamfile> <out_fn> [-n <topn>]
# This script collects the top N longest mapped reads from a BAM/SAM file.
# It outputs a TSV file with read names, read lengths, and aligned lengths.

import heapq
import pysam

def collect_top_reads(bamfile, topn=100):
    """Return top N longest mapped reads from one BAM, with aligned length."""
    heap = []
    counter = 0
    bam = pysam.AlignmentFile(bamfile, "rb" if bamfile.endswith(".bam") else "r")
    for aln in bam:
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
            continue
        read_length = aln.query_length           # full read length
        aligned_length = aln.query_alignment_length    # length aligned to reference
        record = (aln.query_name, read_length, aligned_length)
        if len(heap) < topn:
            heapq.heappush(heap, (aligned_length, counter, record))
        else:
            heapq.heappushpop(heap, (aligned_length, counter, record))
        counter += 1
    bam.close()
    return heap

def write_table(records, outfile):
    """Write table of read name, read length, and aligned length."""
    with open(outfile, "w") as out:
        out.write("read_id\tread_length\taligned_length\n")
        for read_length, _, (name, rl, al) in sorted(records, key=lambda x: -x[0]):
            out.write(f"{name}\t{rl}\t{al}\n")

def main(bamfile, outname, topn):
    heap = collect_top_reads(bamfile, topn)
    write_table(heap, outname)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Find top N longest mapped reads from BAM/SAM files")
    parser.add_argument("bamfile", help="Input BAM/SAM file")
    parser.add_argument("out_fn", help="Output TSV file")
    parser.add_argument("-n", type=int, default=100, help="Number of longest reads [default=100]")
    args = parser.parse_args()
    main(args.bamfile, args.out_fn, args.n)