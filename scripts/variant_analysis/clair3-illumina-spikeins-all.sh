#!/bin/bash
#SBATCH --cluster=genius
#SBATCH --partition=batch
#SBATCH --job-name=clair3_illumina_allctgs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --account=lp_ms_genetics
#SBATCH --mail-type=ALL
#SBATCH --mail-user=margaux.david@kuleuven.be

#Add miniconda to path
export PATH=/data/leuven/343/vsc34345/miniconda3/bin:$PATH

source activate /lustre1/project/stg_00093/Margaux/.conda/envs/clair3_illumina
cd /lustre1/project/stg_00093/Margaux/LongBench/Clair3-illumina

# List of samples
declare -a SAMPLES=("H146" "H1975" "H211" "H2228" "H526" "H69" "HCC827" "SHP77")

# Set paths
BAM_DIR=/lustre1/project/stg_00093/Margaux/LongBench/inputs/illumina
REF=/lustre1/project/stg_00093/Margaux/LongBench/inputs/GRCh38.primary_assembly.genome.add.spikein.fa
THREADS=16
PRESET=ilmn

for SAMPLE in "${SAMPLES[@]}"; do
	    ./run_clair3.sh -b ${BAM_DIR}/${SAMPLE}.sorted.bam -f $REF -m /lustre1/project/stg_00093/Margaux/.conda/envs/clair3_illumina/bin/models/ilmn -o spikeins/all_contigs/${SAMPLE} -t $THREADS -p $PRESET --include_all_ctgs 
	    echo "Finished processing sample: $SAMPLE"    
    done
