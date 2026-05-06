# Detection of CCLE vars

# Create input CCLE file for each cell line

# Separate per sample + remove 0/0 or 0|0 rows
module load bcftools

for SAMPLE in HCC827 NCIH2228 NCIH211 SHP77 NCIH146 NCIH69 NCIH526 NCIH1975
do
  bcftools view -s $SAMPLE -Ov DepmapMutTabConverted.vcf | awk 'BEGIN {FS=OFS="\t"} {if ($0 ~ /^#/ || ($10 !~ /^0\/0/ && $10 !~ /^0\|0/)) print}' > ${SAMPLE}_CCLE.vcf
done

# remove NCI from file name
for FILE in NCI*.vcf
do
  NEW_NAME=$(echo $FILE | sed 's/NCI//')
  mv $FILE $NEW_NAME
done


# Pacbio vs CCLE
cd /vast/scratch/users/david.m/LongBench/compare_vcfs/bcftools/CCLE_pacbio

source activate /vast/scratch/users/david.m/LongBench/conda/rtg-tools
module load bcftools

#bgzip the vcfs
for FILE in /vast/scratch/users/david.m/LongBench/inputs/CCLE/*_CCLE.vcf
do
  bgzip -c $FILE > ${FILE}.gz
done

#index them
for FILE in /vast/scratch/users/david.m/LongBench/inputs/CCLE/*_CCLE.vcf.gz
do
  bcftools index $FILE
done

# Intersection
# Define the base directories for CCLE and PacBio files
CCLE_DIR=/vast/scratch/users/david.m/LongBench/inputs/CCLE
PACBIO_DIR=/vast/scratch/users/david.m/LongBench/Clair3-RNA/spikeins/pacbio

# List of samples
SAMPLES=("H146" "H1975" "H211" "H2228" "H526" "H69" "HCC827" "SHP77")

# Loop through each sample and perform the intersection
for SAMPLE in "${SAMPLES[@]}"
do
  CCLE_FILE=${CCLE_DIR}/${SAMPLE}_CCLE.vcf.gz
  PACBIO_FILE=${PACBIO_DIR}/${SAMPLE}/output.vcf.gz
  OUTPUT_DIR=${SAMPLE}
  bcftools isec -p $OUTPUT_DIR -n=2 $CCLE_FILE $PACBIO_FILE
done

# ONT vs CCLE
cd /vast/scratch/users/david.m/LongBench/compare_vcfs/bcftools/CCLE_ont

source activate /vast/scratch/users/david.m/LongBench/conda/rtg-tools
module load bcftools

# Intersection
# Define the base directories for CCLE and ONT files
CCLE_DIR=/vast/scratch/users/david.m/LongBench/inputs/CCLE
ONT_DIR=/vast/scratch/users/david.m/LongBench/Clair3-RNA/spikeins/ont

# List of samples
SAMPLES=("H146" "H1975" "H211" "H2228" "H526" "H69" "HCC827" "SHP77")

# Loop through each sample and perform the intersection
for SAMPLE in "${SAMPLES[@]}"
do
  CCLE_FILE=${CCLE_DIR}/${SAMPLE}_CCLE.vcf.gz
  ONT_FILE=${ONT_DIR}/${SAMPLE}/output.vcf.gz
  OUTPUT_DIR=${SAMPLE}
  bcftools isec -p $OUTPUT_DIR -n=2 $CCLE_FILE $ONT_FILE
done
 
 #drna vs CCLE
 cd /vast/scratch/users/david.m/LongBench/compare_vcfs/bcftools/CCLE_drna

source activate /vast/scratch/users/david.m/LongBench/conda/rtg-tools
module load bcftools

# Intersection
# Define the base directories for CCLE and ONT files
CCLE_DIR=/vast/scratch/users/david.m/LongBench/inputs/CCLE
DRNA_DIR=/vast/scratch/users/david.m/LongBench/Clair3-RNA/spikeins/drna

# List of samples
SAMPLES=("H146" "H1975" "H211" "H2228" "H526" "H69" "HCC827" "SHP77")

# Loop through each sample and perform the intersection
for SAMPLE in "${SAMPLES[@]}"
do
  CCLE_FILE=${CCLE_DIR}/${SAMPLE}_CCLE.vcf.gz
  DRNA_FILE=${DRNA_DIR}/${SAMPLE}/output.vcf.gz
  OUTPUT_DIR=${SAMPLE}
  bcftools isec -p $OUTPUT_DIR -n=2 $CCLE_FILE $DRNA_FILE
done

# illumina vs CCLE
cd /vast/scratch/users/david.m/LongBench/compare_vcfs/bcftools/CCLE_illumina

source activate /vast/scratch/users/david.m/LongBench/conda/rtg-tools
module load bcftools

# Intersection
# Define the base directories for CCLE and ONT files
CCLE_DIR=/vast/scratch/users/david.m/LongBench/inputs/CCLE
ILL_DIR=/vast/scratch/users/david.m/LongBench/Clair3-illumina/spikeins

# List of samples
SAMPLES=("H146" "H1975" "H211" "H2228" "H526" "H69" "HCC827" "SHP77")

# Loop through each sample and perform the intersection
for SAMPLE in "${SAMPLES[@]}"
do
  CCLE_FILE=${CCLE_DIR}/${SAMPLE}_CCLE.vcf.gz
  ILL_FILE=${ILL_DIR}/${SAMPLE}/merge_output.vcf.gz
  OUTPUT_DIR=${SAMPLE}
  bcftools isec -p $OUTPUT_DIR -n=2 $CCLE_FILE $ILL_FILE
done

#Illumina called by DV vs CCLE
cd /vast/scratch/users/david.m/LongBench/compare_vcfs/bcftools/CCLE_DV_illumina

source activate /vast/scratch/users/david.m/LongBench/conda/rtg-tools
module load bcftools

# Intersection
# Define the base directories for CCLE and ONT files
CCLE_DIR=/vast/scratch/users/david.m/LongBench/inputs/CCLE
DV_DIR=/vast/projects/LongBench/analysis/sr_bulk/result/Deepvariants

# List of samples
SAMPLES=("H146" "H1975" "H211" "H2228" "H526" "H69" "HCC827" "SHP77")

# Loop through each sample and perform the intersection
for SAMPLE in "${SAMPLES[@]}"
do
  CCLE_FILE=${CCLE_DIR}/${SAMPLE}_CCLE.vcf.gz
  DV_FILE=${DV_DIR}/${SAMPLE}/output.vcf.gz
  OUTPUT_DIR=${SAMPLE}
  bcftools isec -p $OUTPUT_DIR -n=2 $CCLE_FILE $DV_FILE
done

# Create summary table with quality filtering (QUAL < 10 -> lowqual)
cd /vast/scratch/users/david.m/LongBench/compare_vcfs/bcftools

CCLE_VCF_DIR="/vast/scratch/users/david.m/LongBench/inputs/CCLE"
PLATFORMS=("illumina" "DV_illumina" "pacbio" "ont" "drna")
SAMPLES=("H146" "H1975" "H2228" "H211" "H69" "H526" "HCC827" "SHP77")

OUTDIR="genotype_comparisons_qual"
mkdir -p "$OUTDIR"

for SAMPLE in "${SAMPLES[@]}"; do
  CCLE_VCF="${CCLE_VCF_DIR}/${SAMPLE}_CCLE.vcf.gz"
  [[ ! -f "$CCLE_VCF" ]] && echo "Missing CCLE VCF for $SAMPLE" && continue

  unset GT_DICT
  declare -A GT_DICT

  while IFS=$'\t' read -r CHROM POS ID REF ALT QUAL FILTER INFO FORMAT GT_CCLE; do
    KEY="${CHROM}_${POS}_${REF}_${ALT}"
    GT_DICT["$KEY"]="$GT_CCLE"
  done < <(bcftools view -H "$CCLE_VCF")

  echo "[$SAMPLE] Loaded ${#GT_DICT[@]} unique variants from $CCLE_VCF"

  OUTPUT="${OUTDIR}/genotype_comparison_${SAMPLE}.tsv"
  HEADER="CHROM\tPOS\tREF\tALT\tGT_CCLE"
  for PLATFORM in "${PLATFORMS[@]}"; do
    HEADER+="\t${PLATFORM}_GT\t${PLATFORM}_QUAL\t${PLATFORM}_STATUS"
  done
  echo -e "$HEADER" > "$OUTPUT"

  for KEY in "${!GT_DICT[@]}"; do
    IFS="_" read -r CHROM POS REF ALT <<< "$KEY"
    GT_CCLE="${GT_DICT[$KEY]}"
    SORTED_GT_CCLE=$(echo "$GT_CCLE" | tr '|' '/' | tr '/' '\n' | sort | paste -sd / -)

    LINE="${CHROM}\t${POS}\t${REF}\t${ALT}\t${GT_CCLE}"

    for PLATFORM in "${PLATFORMS[@]}"; do
      VCF="CCLE_${PLATFORM}/${SAMPLE}/0001.vcf"
      if [[ ! -f "$VCF" ]]; then
        LINE+="\tNA\tNA\tmissing_vcf"
        continue
      fi

      MATCH_LINE=$(bcftools view -H "$VCF" | awk -v c="$CHROM" -v p="$POS" -v r="$REF" -v a="$ALT" '$1==c && $2==p && $4==r && $5==a' | head -n 1)

      if [[ -z "$MATCH_LINE" ]]; then
        LINE+="\tnot_called\tNA\tnot_called"
      else
        QUAL=$(echo "$MATCH_LINE" | awk '{print $6}')
        GT_PLATFORM=$(echo "$MATCH_LINE" | awk '{print $10}' | cut -d ':' -f1)

        if (( $(echo "$QUAL < 10" | bc -l) )); then
          LINE+="\t${GT_PLATFORM}\t${QUAL}\tlowqual"
        else
          SORTED_GT_PLATFORM=$(echo "$GT_PLATFORM" | tr '|' '/' | tr '/' '\n' | sort | paste -sd / -)
          STATUS=$([[ "$SORTED_GT_PLATFORM" == "$SORTED_GT_CCLE" ]] && echo "match" || echo "discordant")
          LINE+="\t${GT_PLATFORM}\t${QUAL}\t${STATUS}"
        fi
      fi
    done

    echo -e "$LINE" >> "$OUTPUT"
  done

  echo "Finished $SAMPLE → $OUTPUT"
done

#Annotate variant type
cd /vast/scratch/users/david.m/LongBench/compare_vcfs/bcftools
IN_DIR="genotype_comparisons_qual"

for FILE in "${IN_DIR}"/genotype_comparison_*.tsv; do
  TMP="${FILE}.tmp"
  
  # Add new header with 'variant_type'
  awk 'NR==1 { print $0 "\tvariant_type"; next } 
       NR>1 {
         len_ref=length($3); len_alt=length($4);
         if (len_ref == 1 && len_alt == 1) type="SNP";
         else if (len_ref != len_alt) type="indel";
         else type="other";
         print $0 "\t" type;
       }' "$FILE" > "$TMP"

  mv "$TMP" "$FILE"
  echo "Annotated $FILE"
done


