#!/bin/bash

#SBATCH --job-name=longbench_ont_dRNA_merged_topup_subset_8.5M
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
shopt -s nullglob  # so the loop doesn’t iterate the literal pattern if no files

sample_name="${1:?Usage: $0 SAMPLE_NAME}"

raw_dir="/vast/scratch/users/tan.j/longbench_fusion_calling/dataset/longbench_ont_dRNA/fastq_q10_filter_merged_with_topup"
subset_dir="/vast/scratch/users/tan.j/longbench_fusion_calling/dataset/longbench_ont_dRNA/merged_topup_subset_8.5M"

mkdir -p "$subset_dir"

in="$raw_dir/${sample_name}.fastq.gz"
out="$subset_dir/${sample_name}_dRNA_ONT_merged_subset_8.5M.fastq.gz"

# Sample exactly 8,500,000 reads with a fixed seed for reproducibility
seqtk sample -s100 "$in" 8500000 | gzip -c > "$out"