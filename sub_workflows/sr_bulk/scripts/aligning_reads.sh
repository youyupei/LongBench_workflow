#!/bin/bash
#SBATCH --job-name=align_readsLB
#SBATCH --cpus-per-task=6
#SBATCH --mem=60GB
#SBATCH --time=10:00:00
#SBATCH --partition=regular
#SBATCH --mail-type=ALL
#SBATCH --mail-user=solano.a@wehi.edu.au

module load R

Rscript --verbose --vanilla align_reads.R