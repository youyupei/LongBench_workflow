#!/bin/bash
#SBATCH --cluster=genius
#SBATCH --partition=batch
#SBATCH --job-name=clair3_drna_allctgs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --account=lp_ms_genetics
#SBATCH --mail-type=ALL
#SBATCH --mail-user=margaux.david@kuleuven.be

#Add miniconda to path
export PATH=/data/leuven/343/vsc34345/miniconda3/bin:$PATH

source activate /lustre1/project/stg_00093/Margaux/.conda/envs/clair3_rna
cd /lustre1/project/stg_00093/Margaux/LongBench/Clair3-RNA

# List of samples
declare -a SAMPLES=("H146" "H1975" "H211" "H2228" "H526" "H69" "HCC827" "SHP77")

# Set paths
BAM_DIR=/lustre1/project/stg_00093/Margaux/LongBench/inputs
REF=/lustre1/project/stg_00093/Margaux/LongBench/inputs/GRCh38.primary_assembly.genome.add.spikein.fa
THREADS=32
PRESET=ont_dorado_drna004

for SAMPLE in "${SAMPLES[@]}"; do
	    ./run_clair3_rna -B ${BAM_DIR}/dRNA_bulk_${SAMPLE}.sorted.bam -R $REF -o spikeins/all_contigs/drna/${SAMPLE} -t $THREADS -p $PRESET --enable_phasing_model --include_all_ctgs 
	    echo "Finished processing sample: $SAMPLE"    
    done
