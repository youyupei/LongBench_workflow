#!/bin/bash

#SBATCH --job-name=longbench_pacbio_filter_q10
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=400G
#SBATCH --time=10:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=tan.j@wehi.edu.au
#SBATCH --output=sbatch_output/%x.out
#SBATCH --error=sbatch_output/%x.err

set -euo pipefail

sample_name="${1:?Usage: $0 SAMPLE_NAME}"

raw_dir="/vast/scratch/users/tan.j/longbench_fusion_calling/dataset/longbench_pacbio/raw_data"
filter_q10_dir="/vast/scratch/users/tan.j/longbench_fusion_calling/dataset/longbench_pacbio/filter_q10"

mkdir -p "$filter_q10_dir"

input="$raw_dir/${sample_name}.fastq.gz"
output="${filter_q10_dir}/${sample_name}.fastq.gz"

gunzip -c "$input" \
  | /vast/scratch/users/tan.j/tools/chopper-linux -q 10 --threads 8 \
  | gzip -c \
  > "$output"
