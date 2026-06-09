#!/bin/bash

#SBATCH --job-name=longbench_ont_2nd_trim_adapter_remnants
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
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

input_fastq="/vast/scratch/users/tan.j/longbench_fusion_calling/dataset/longbench_ont/1st_trimmed_data/${sample_name}.fastq.gz"
trimmed_dir="/vast/scratch/users/tan.j/longbench_fusion_calling/dataset/longbench_ont/2nd_trimmed_data_noTrimPolyA"

mkdir -p "$trimmed_dir"

# -O: Requires ≥18 bp overlap; -e: Allows 10% errors; -m: Keeps reads ≥100 nt (-m) after trimming, -j: Uses 48 threads   
cutadapt --report=full -j 48 \
	--revcomp -O 18 -e 0.1 -m 100 \
  -g ^CTTGCGGGCGGCGGACTCTCCTCT \
  -a CTTGCGGGCGGCGGACTCTCCTCT$ \
  -o "${trimmed_dir}/${sample_name}.fastq.gz" \
  "$input_fastq"
