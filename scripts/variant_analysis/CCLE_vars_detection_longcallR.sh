#!/bin/bash
# TP detection for LongcallR VCFs: intersect pb/ont/dRNA against CCLE ground truth.
# No manual QUAL threshold — records QUAL as-is for downstream filtering in Rmd.
#
# Usage:
#   CCLE_vars_detection_longcallR.sh <ccle_source_vcf> <ccle_dir> <outdir> <longcallR_dir>
#
# longcallR_dir: flat directory containing {pb_bulk|ont_bulk|dRNA_bulk}_{cell_line}.longcallR.vcf
# Output: {outdir}/genotype_comparisons_qual/genotype_comparison_{cell_line}.tsv

set -euo pipefail
module load bcftools
module load htslib  # needed for bgzip/tabix when creating CCLE VCFs (Step 1)

if [[ $# -lt 4 ]]; then
    echo "Usage: $0 <ccle_source_vcf> <ccle_dir> <outdir> <longcallR_dir>" >&2
    exit 1
fi

CCLE_SOURCE="$1"
CCLE_DIR="$2"
OUTDIR="$3"
LCR_DIR="$4"

SAMPLES=(H146 H1975 H211 H2228 H526 H69 HCC827 SHP77)

declare -A CCLE_NAME=(
    [H146]=NCIH146 [H1975]=NCIH1975 [H211]=NCIH211 [H2228]=NCIH2228
    [H526]=NCIH526 [H69]=NCIH69 [HCC827]=HCC827 [SHP77]=SHP77
)

declare -A PLATFORM_PREFIX=(
    [pacbio]=pb_bulk
    [ont]=ont_bulk
    [drna]=dRNA_bulk
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

# --- Step 2: Genotype comparison summary table ---
echo "=== Generating genotype comparison table ==="
SUMMARY_DIR="${OUTDIR}/genotype_comparisons_qual"
mkdir -p "$SUMMARY_DIR"

PLATFORMS=(pacbio ont drna)

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
            VCF="${LCR_DIR}/${PLATFORM_PREFIX[$PLATFORM]}_${SAMPLE}.longcallR.vcf"
            if [[ ! -f "$VCF" ]]; then
                LINE+="\tNA\tNA\tmissing_vcf"
                continue
            fi

            MATCH_LINE=$(awk -v c="$CHROM" -v p="$POS" -v r="$REF" -v a="$ALT" \
                '$1==c && $2==p && $4==r && $5==a' "$VCF" | head -n 1)

            if [[ -z "$MATCH_LINE" ]]; then
                LINE+="\tnot_called\tNA\tnot_called"
            else
                QUAL=$(echo "$MATCH_LINE" | awk '{print $6}')
                GT_PLATFORM=$(echo "$MATCH_LINE" | awk '{print $10}' | cut -d ':' -f1)
                SORTED_GT_PLATFORM=$(echo "$GT_PLATFORM" | tr '|' '/' | tr '/' '\n' | sort | paste -sd / -)
                STATUS=$([[ "$SORTED_GT_PLATFORM" == "$SORTED_GT_CCLE" ]] && echo "match" || echo "discordant")
                LINE+="\t${GT_PLATFORM}\t${QUAL}\t${STATUS}"
            fi
        done

        echo -e "$LINE" >> "$OUTPUT"
    done

    echo "Finished $SAMPLE → $OUTPUT"
done

# --- Step 5: Annotate variant type ---
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
