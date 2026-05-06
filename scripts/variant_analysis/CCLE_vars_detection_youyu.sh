#!/bin/bash
# TP detection: intersect each platform's clair3 VCFs against CCLE ground truth
# Adapted from Margaux's CCLE_vars_detection.sh
#
# Usage:
#   CCLE_vars_detection_youyu.sh <ccle_source_vcf> <ccle_dir> <outdir> \
#       <pacbio_vcf_dir> <ont_vcf_dir> <drna_vcf_dir> <ill_vcf_dir> <dv_vcf_dir> \
#       [lr_vcf_name]
#
# *_vcf_dir: directory containing per-sample subdirs with VCF files inside.
#   LR platforms (pacbio/ont/drna): <dir>/<sample>/<lr_vcf_name>  (default: output.vcf.gz)
#   illumina clair3:                <dir>/<sample>/merge_output.vcf.gz
#   deepvariant:                    <dir>/<sample>/output.vcf.gz
#
# lr_vcf_name: optional 9th arg, e.g. output_enable_phasing.vcf.gz for phased LR

set -euo pipefail
module load bcftools
module load htslib

if [[ $# -lt 8 ]]; then
    echo "Usage: $0 <ccle_source_vcf> <ccle_dir> <outdir> <pacbio_dir> <ont_dir> <drna_dir> <ill_dir> <dv_dir> [lr_vcf_name]" >&2
    exit 1
fi

CCLE_SOURCE="$1"
CCLE_DIR="$2"
OUTDIR="$3"
PACBIO_DIR="$4"
ONT_DIR="$5"
DRNA_DIR="$6"
ILL_DIR="$7"
DV_DIR="$8"
LR_VCF_NAME="${9:-output.vcf.gz}"

SAMPLES=(H146 H1975 H211 H2228 H526 H69 HCC827 SHP77)

declare -A CCLE_NAME=(
    [H146]=NCIH146 [H1975]=NCIH1975 [H211]=NCIH211 [H2228]=NCIH2228
    [H526]=NCIH526 [H69]=NCIH69 [HCC827]=HCC827 [SHP77]=SHP77
)

declare -A PLATFORM_VCF_DIR=(
    [pacbio]="$PACBIO_DIR"
    [ont]="$ONT_DIR"
    [drna]="$DRNA_DIR"
    [illumina]="$ILL_DIR"
    [deepvariant]="$DV_DIR"
)

declare -A PLATFORM_VCF_NAME=(
    [pacbio]="$LR_VCF_NAME"
    [ont]="$LR_VCF_NAME"
    [drna]="$LR_VCF_NAME"
    [illumina]=merge_output.vcf.gz
    [deepvariant]=output.vcf.gz
)

mkdir -p "$CCLE_DIR" "$OUTDIR"

# --- Step 1: Create per-sample CCLE VCFs (skip if already exist) ---
echo "=== Creating per-sample CCLE VCFs ==="
for SAMPLE in "${SAMPLES[@]}"; do
    CCLE_SAMPLE="${CCLE_NAME[$SAMPLE]}"
    OUT_VCF="${CCLE_DIR}/${SAMPLE}_CCLE.vcf"
    if [[ -f "${OUT_VCF}.gz" ]]; then
        echo "  $SAMPLE already exists, skipping"
        continue
    fi
    echo "  $CCLE_SAMPLE → ${SAMPLE}_CCLE.vcf"
    bcftools view -s "$CCLE_SAMPLE" -Ov "$CCLE_SOURCE" \
        | awk 'BEGIN {FS=OFS="\t"} {if ($0 ~ /^#/ || ($10 !~ /^0\/0/ && $10 !~ /^0\|0/)) print}' \
        > "$OUT_VCF"
    bgzip -f "$OUT_VCF"
    bcftools index -f "${OUT_VCF}.gz"
done

# --- Step 2: bcftools isec for each platform ---
for PLATFORM in pacbio ont drna illumina deepvariant; do
    echo "=== Intersecting CCLE vs $PLATFORM ==="
    ISEC_DIR="${OUTDIR}/CCLE_${PLATFORM}"
    mkdir -p "$ISEC_DIR"

    for SAMPLE in "${SAMPLES[@]}"; do
        CCLE_VCF="${CCLE_DIR}/${SAMPLE}_CCLE.vcf.gz"
        PLATFORM_VCF="${PLATFORM_VCF_DIR[$PLATFORM]}/${SAMPLE}/${PLATFORM_VCF_NAME[$PLATFORM]}"
        SAMPLE_OUTDIR="${ISEC_DIR}/${SAMPLE}"

        if [[ ! -f "$PLATFORM_VCF" ]]; then
            echo "  WARNING: missing $PLATFORM_VCF — skipping $SAMPLE"
            continue
        fi
        echo "  $SAMPLE"
        bcftools isec -p "$SAMPLE_OUTDIR" -n=2 "$CCLE_VCF" "$PLATFORM_VCF"
    done
done

# --- Step 3: Genotype comparison summary table ---
echo "=== Generating genotype comparison table ==="
SUMMARY_DIR="${OUTDIR}/genotype_comparisons_qual"
mkdir -p "$SUMMARY_DIR"

PLATFORMS=(illumina deepvariant pacbio ont drna)

for SAMPLE in "${SAMPLES[@]}"; do
    CCLE_VCF="${CCLE_DIR}/${SAMPLE}_CCLE.vcf.gz"
    [[ ! -f "$CCLE_VCF" ]] && echo "Missing CCLE VCF for $SAMPLE" && continue

    unset GT_DICT
    declare -A GT_DICT

    while IFS=$'\t' read -r CHROM POS ID REF ALT QUAL FILTER INFO FORMAT GT_CCLE; do
        KEY="${CHROM}_${POS}_${REF}_${ALT}"
        GT_DICT["$KEY"]="$GT_CCLE"
    done < <(bcftools view -H "$CCLE_VCF")

    echo "[$SAMPLE] Loaded ${#GT_DICT[@]} unique variants"

    OUTPUT="${SUMMARY_DIR}/genotype_comparison_${SAMPLE}.tsv"
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
            VCF="${OUTDIR}/CCLE_${PLATFORM}/${SAMPLE}/0001.vcf"
            if [[ ! -f "$VCF" ]]; then
                LINE+="\tNA\tNA\tmissing_vcf"
                continue
            fi

            MATCH_LINE=$(bcftools view -H "$VCF" | awk -v c="$CHROM" -v p="$POS" -v r="$REF" -v a="$ALT" \
                '$1==c && $2==p && $4==r && $5==a' | head -n 1)

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

# --- Step 4: Annotate variant type ---
echo "=== Annotating variant types ==="
for FILE in "${SUMMARY_DIR}"/genotype_comparison_*.tsv; do
    TMP="${FILE}.tmp"
    awk 'NR==1 { print $0 "\tvariant_type"; next }
         NR>1 {
           len_ref=length($3); len_alt=length($4);
           if (len_ref == 1 && len_alt == 1) type="SNP";
           else if (len_ref != len_alt) type="indel";
           else type="other";
           print $0 "\t" type;
         }' "$FILE" > "$TMP"
    mv "$TMP" "$FILE"
done

echo "Done. Results in ${SUMMARY_DIR}"
