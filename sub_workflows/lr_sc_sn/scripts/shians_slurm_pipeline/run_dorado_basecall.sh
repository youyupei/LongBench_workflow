#!/bin/bash

### PREAMBLE ####
module load dorado/0.7.0
module load samtools/1.19.2

#################
# usage: run_dorado_basecall_6mA.sh <dir> <block_size>
DIR=$1
BLOCK_SIZE=$2

# check if sufficient arguments are provided
if [ -z "$DIR" ] || [ -z "$BLOCK_SIZE" ]; then
    echo "Error: insufficient arguments provided"
    echo "usage: run_dorado_basecall_6mA.sh <dir> <block_size>"
    exit 1
fi

# variables
USER_EMAIL=${USER}@wehi.edu.au
MODEL=$DORADO_MODELS/dna_r10.4.1_e8.2_400bps_sup@v5.0.0

# setup
EMAIL_CONFIG="--mail-user=$USER_EMAIL --mail-type=ALL"
LOG_CONFIG="--output logs/%x-%j_%a.log --error logs/%x-%j_%a.log"
REPORTING_CONFIG="$EMAIL_CONFIG $LOG_CONFIG"

### HELPER FUNCTIONS ###
########################
# logging function, prints to stderr instead of stdout and includes a timestamp
log() {
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] $1" >&2
}

export -f log

setup_tmp_folder() {
    FOLDER=$1
    if [ ! -d $FOLDER ]; then
    	log "creating $FOLDER..."
        mkdir $FOLDER
    elif [ ! -z "$(ls -A $FOLDER)" ]; then
    	log "$FOLDER exists, removing contents..."
        rm -rf $FOLDER/*
    fi
}

export -f setup_tmp_folder

# split files into folders with a maximum number of files
# usage: split_pod5_files <dir> <block_size>
# <dir> is the directory containing the files to be split
# <block_size> is the maximum number of files in each folder
# returns the number of blocks created
split_pod5_files() {
    local dir=$1
    local block_size=$2

    log "splitting pod5 files in $dir..."
    log "block size: $block_size"

    local files=($(find $PWD/$dir -name "*.pod5"))

    local file_number=${#files[@]}
    log "${#files[@]} pod5 files found..."

    local block_number=$((($file_number + $block_size - 1) / $block_size ))
    log "$block_number blocks will be created..."

    # create the block folders, delete if already exists
    setup_tmp_folder pod5_tmp
    
    for ((i=0; i<$block_number; i++)); do
        mkdir pod5_tmp/block_$i
    done

    log "linking files to block folders..."
    # link the files to the block folders
    for ((i=0; i<$file_number; i++)); do
        local block_index=$((i / $block_size))
        local file=${files[$i]}
        local file_name=$(basename $file)
        local new_file=pod5_tmp/block_$block_index/$file_name

        ln -s $file $new_file
    done

    log "done!"

    # print the number of blocks created
    echo $block_number
}

export -f split_pod5_files

# sort_and_index_bam() {
#     local bam=$1
#     local sorted_bam=${bam/.bam/.sorted.bam}
#     samtools sort -@ 8 $bam > $sorted_bam && mv $sorted_bam $bam
#     samtools index -@ 8 $bam
# }

# export -f sort_and_index_bam

### PIPELINE ####
#################
# make a temporary directory for the bam files if none exist
if [ ! -d bam_tmp ]; then
    mkdir bam_tmp
fi
setup_tmp_folder demux

# RUN pod5 split
BLOCK_COUNT=$(split_pod5_files $DIR $BLOCK_SIZE)

# RUN dorado basecaller
# awk '{print $NF}' is used to extract the job id from the sbatch output
log "running dorado on $BLOCK_COUNT blocks..."
DORADO_JOB_ID=$(
    sbatch $REPORTING_CONFIG --array=0-$((BLOCK_COUNT-1)) scripts/dorado_basecall_block.sh $MODEL | \
    awk '{print $NF}'
)

# RUN samtools merge
log "submitting merge job..."
BAM_MERGE_JOB_ID=$(
    sbatch $REPORTING_CONFIG --dependency=afterok:$DORADO_JOB_ID \
        $REPORTING_CONFIG \
        --job-name="samtools_merge" \
        --mem=2G \
        --cpus-per-task=8 \
        --time="48:00:00" \
        --wrap "samtools cat -@ 8 -o bam_tmp/merged_bam.bam bam_tmp/bam_*.bam" | \
        awk '{print $NF}'
)

# CLEANUP temp files
log "submitting cleanup job..."
sbatch $REPORTING_CONFIG --dependency=afterok:$BAM_INDEX_JOB_ID \
    --job-name="cleanup" \
    --mem=2G \
    --cpus-per-task=1 \
    --time="2:00:00" \
    --wrap "rm -rf pod5_tmp bam_tmp"

log "dorado basecall pipeline submitted!"
