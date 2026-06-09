#!/bin/bash

#SBATCH --job-name=longbench_ont_trimmed_filter_q10
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH --time=2:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=tan.j@wehi.edu.au
#SBATCH --output=sbatch_output/%x.out
#SBATCH --error=sbatch_output/%x.err

set -euo pipefail

sample_name="${1:?Usage: $0 SAMPLE_NAME}"

trimmed_dir="/vast/scratch/users/tan.j/longbench_fusion_calling/dataset/longbench_ont/1st_trimmed_data"
filter_q10_dir="/vast/scratch/users/tan.j/longbench_fusion_calling/dataset/longbench_ont/trimmed_filter_q10"

mkdir -p "$filter_q10_dir"

input="$trimmed_dir/${sample_name}.fastq.gz"
output="${filter_q10_dir}/${sample_name}.fastq.gz"

gunzip -c "$input" \
  | /vast/scratch/users/tan.j/tools/chopper-linux -q 10 --threads 32 \
  | gzip -c \
  > "$output"