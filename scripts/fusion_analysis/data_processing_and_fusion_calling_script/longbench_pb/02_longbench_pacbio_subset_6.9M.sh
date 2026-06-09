#!/bin/bash

#SBATCH --job-name=longbench_pacbio_subset_6.9M
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=tan.j@wehi.edu.au
#SBATCH --output=sbatch_output/%x.out
#SBATCH --error=sbatch_output/%x.err

set -euo pipefail
shopt -s nullglob  # so the loop doesn’t iterate the literal pattern if no files

sample_name="${1:?Usage: $0 SAMPLE_NAME}"

filter_q10_dir="/vast/scratch/users/tan.j/longbench_fusion_calling/dataset/longbench_pacbio/filter_q10"
subset_dir="/vast/scratch/users/tan.j/longbench_fusion_calling/dataset/longbench_pacbio/subset_6.9M"

mkdir -p "$subset_dir"

in="$filter_q10_dir/${sample_name}.fastq.gz"
out="$subset_dir/${sample_name}_subset_6.9M.fastq.gz"

# Sample exactly 6,900,000 reads with a fixed seed for reproducibility
seqtk sample -s100 "$in" 6900000 | gzip -c > "$out"