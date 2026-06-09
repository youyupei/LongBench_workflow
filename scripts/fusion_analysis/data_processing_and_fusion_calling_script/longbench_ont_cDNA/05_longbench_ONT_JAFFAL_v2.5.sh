#!/bin/bash

#SBATCH --job-name=longbench_ONT_JAFFAL_v2.5
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=tan.j@wehi.edu.au
#SBATCH --output=script_output/%x_%J.out
#SBATCH --error=script_output/%x_%J.err

set -euo pipefail

module purge
module load micromamba

eval "$(micromamba shell hook --shell=bash)"
micromamba activate /vast/scratch/users/tan.j/condaenvs/jaffal_env
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

sample_name="${1:?Usage: $0 SAMPLE_NAME}"

input_fastq="/vast/scratch/users/tan.j/longbench_fusion_calling/dataset/longbench_ont/trimmed_subset_16.9M/${sample_name}.fastq.gz"
output_dir="/vast/scratch/users/tan.j/longbench_fusion_calling/output/longread_ONT/${sample_name}"

mkdir -p "$output_dir"

if [[ ! -f "$input_fastq" ]]; then
    echo "ERROR: input FASTQ not found: $input_fastq"
    exit 1
fi

cd "$output_dir"

/vast/scratch/users/tan.j/tools/JAFFA/tools/bin/bpipe run \
/vast/scratch/users/tan.j/tools/JAFFA/JAFFAL.groovy \
"$input_fastq"
