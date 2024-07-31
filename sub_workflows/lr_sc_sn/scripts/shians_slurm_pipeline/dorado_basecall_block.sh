#!/bin/bash
#SBATCH --job-name=dorado_basecalling
#SBATCH --partition="gpuq"
#SBATCH --mem=128G
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:A30:4
#SBATCH --time="2:00:00"
#
#SBATCH --output logs/%x-%j_%a.log
#SBATCH --error logs/%x-%j_%a.log
#
#SBATCH --mail-type="ALL"

#usage: sbatch dorado_basecall_block.sh <reference> <kit> <mod> <model>

module load samtools
module load dorado/0.7.0

#REF=$1 # not used because not barcoded
#KIT=$2 # not used because not barcoded
#MOD=$3 # not used
MODEL=$1

BAM_FILE=bam_tmp/bam_$SLURM_ARRAY_TASK_ID.bam
POD5_DIR=pod5_tmp/block_$SLURM_ARRAY_TASK_ID

set -x

# check if the bam file already exists
if [ -f $BAM_FILE ]; then
    exit 0
fi

if [ -f $BAM_FILE.incomplete ]; then
    mv $BAM_FILE.incomplete $BAM_FILE.incomplete.2
    dorado basecaller --recursive \
        # --kit-name $KIT \ # not used because not barcoded
        # --reference $REF \
        # --modified-bases $MOD \ not used
        # --secondary "no" \
        --resume-from $BAM_FILE.incomplete.2 \
        $MODEL $POD5_DIR > $BAM_FILE.incomplete
    rm $BAM_FILE.incomplete.2
else
    dorado basecaller --recursive \
        #--kit-name $KIT \
        #--reference $REF \
        #--modified-bases $MOD \
        #--secondary "no" \
        $MODEL $POD5_DIR > $BAM_FILE.incomplete
fi

mv $BAM_FILE.incomplete $BAM_FILE
