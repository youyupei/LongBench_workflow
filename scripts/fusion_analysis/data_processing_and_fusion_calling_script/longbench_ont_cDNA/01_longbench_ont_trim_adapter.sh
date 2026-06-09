#!/bin/bash

#SBATCH --job-name=longbench_ont_trim_adapter
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=tan.j@wehi.edu.au
#SBATCH --output=sbatch_output/%x.out
#SBATCH --error=sbatch_output/%x.err

set -euo pipefail

module load cutadapt/

sample_name="$1"

if [[ -z "$sample_name" ]]; then
  echo "Usage: $0 SAMPLE_NAME" >&2
  exit 1
fi

input_fastq="/vast/scratch/users/tan.j/longbench_fusion_calling/dataset/longbench_ont/raw_data/${sample_name}.fastq.gz"
trimmed_dir="/vast/scratch/users/tan.j/longbench_fusion_calling/dataset/longbench_ont/1st_trimmed_data"

mkdir -p "$trimmed_dir"

# -O: Requires ≥25 bp overlap; -e: Allows 10% errors; -m: Keeps reads ≥100 nt (-m) after trimming, -j: Uses 16 threads   
cutadapt --report=full -j 16 \
	--revcomp -O 25 -e 0.1 -m 100 \
  -g TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT \
  -g GCAATACGTAACTGAACGAAGTACAGG \
  -g ACGTAACTGAACGAAGTACAGG \
  -g CTTGCCTGTCGCTCTATCTTCAGAGGAG \
  -g TTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG \
  -g CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAG \
  -g CTTGCGGGCGGCGGACTCTCCTC \
  -g ATCGCCTACCGTGACAAGAAAGTTGTCGGTGTCTTTGTGACTTGCCTGTCGCTCTATCTTC \
  -g ATCGCCTACCGTGACAAGAAAGTTGTCGGTGTCTTTGTGTTTCTGTTGGTGCTGATATTGC \
  -a TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT \
  -a GCAATACGTAACTGAACGAAGTACAGG \
  -a ACGTAACTGAACGAAGTACAGG \
  -a CTTGCCTGTCGCTCTATCTTCAGAGGAG \
  -a TTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG \
  -a CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAG \
  -a CTTGCGGGCGGCGGACTCTCCTC \
  -a ATCGCCTACCGTGACAAGAAAGTTGTCGGTGTCTTTGTGACTTGCCTGTCGCTCTATCTTC \
  -a ATCGCCTACCGTGACAAGAAAGTTGTCGGTGTCTTTGTGTTTCTGTTGGTGCTGATATTGC \
  --info-file "${trimmed_dir}/${sample_name}_adapter_hits.txt" \
  -o "${trimmed_dir}/${sample_name}.fastq.gz" \
  "$input_fastq"