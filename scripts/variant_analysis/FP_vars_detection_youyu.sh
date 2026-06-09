#!/bin/bash
# FP detection: variants on non-standard contigs with QUAL >= 10
# Adapted from Margaux's FP_vars_detection.sh
#
# Usage:
#   FP_vars_detection_youyu.sh <out_tsv> \
#       <pacbio_vcf_dir> <ont_vcf_dir> <drna_vcf_dir> <ill_vcf_dir> <dv_vcf_dir> \
#       [lr_vcf_name]
#
# *_vcf_dir: directory containing per-sample subdirs with VCF files inside.
#   LR platforms (pacbio/ont/drna): <dir>/<sample>/<lr_vcf_name>  (default: output.vcf.gz)
#   illumina clair3:                <dir>/<sample>/merge_output.vcf.gz
#   deepvariant:                    <dir>/<sample>/output.vcf.gz
# Writes FP_QUAL_only.vcf next to each input VCF as a side effect.
#
# lr_vcf_name: optional 7th arg, e.g. output_enable_phasing.vcf.gz for phased LR
# lr_only:     optional 8th arg, set to "true" to skip illumina and deepvariant

set -euo pipefail

if [[ $# -lt 6 ]]; then
    echo "Usage: $0 <out_tsv> <pacbio_dir> <ont_dir> <drna_dir> <ill_dir> <dv_dir> [lr_vcf_name] [lr_only]" >&2
    exit 1
fi

OUT_TSV="$1"
PACBIO_DIR="$2"
ONT_DIR="$3"
DRNA_DIR="$4"
ILL_DIR="$5"
DV_DIR="$6"
LR_VCF_NAME="${7:-output.vcf.gz}"
LR_ONLY="${8:-false}"

SAMPLES=(H146 H1975 H211 H2228 H526 H69 HCC827 SHP77)
STD_CTGS='^(chr([1-9]|1[0-9]|2[0-2]|X|Y|M)|GL|KI)'

declare -A PLATFORM_DIR=(
    [pacbio]="$PACBIO_DIR"
    [ont]="$ONT_DIR"
    [drna]="$DRNA_DIR"
)

# --- Step 1: Extract FP candidates ---
echo "=== Extracting non-standard-contig variants ==="

for PLATFORM in pacbio ont drna; do
    for SAMPLE in "${SAMPLES[@]}"; do
        IN_VCF="${PLATFORM_DIR[$PLATFORM]}/${SAMPLE}/${LR_VCF_NAME}"
        OUT_VCF="${PLATFORM_DIR[$PLATFORM]}/${SAMPLE}/FP_QUAL_only.vcf"
        if [[ ! -f "$IN_VCF" ]]; then
            echo "  WARNING: missing $IN_VCF"
            continue
        fi
        echo "  $PLATFORM / $SAMPLE"
        zcat "$IN_VCF" | grep "^#" > "$OUT_VCF"
        zcat "$IN_VCF" | awk -v pat="$STD_CTGS" '$1 !~ /^#/ && $1 !~ pat' >> "$OUT_VCF"
    done
done

if [[ "$LR_ONLY" != "true" ]]; then
    for SAMPLE in "${SAMPLES[@]}"; do
        IN_VCF="${ILL_DIR}/${SAMPLE}/merge_output.vcf.gz"
        OUT_VCF="${ILL_DIR}/${SAMPLE}/FP_QUAL_only.vcf"
        if [[ ! -f "$IN_VCF" ]]; then
            echo "  WARNING: missing $IN_VCF"
            continue
        fi
        echo "  illumina / $SAMPLE"
        zcat "$IN_VCF" | grep "^#" > "$OUT_VCF"
        zcat "$IN_VCF" | awk -v pat="$STD_CTGS" '$1 !~ /^#/ && $1 !~ pat' >> "$OUT_VCF"
    done

    for SAMPLE in "${SAMPLES[@]}"; do
        IN_VCF="${DV_DIR}/${SAMPLE}/output.vcf.gz"
        OUT_VCF="${DV_DIR}/${SAMPLE}/FP_QUAL_only.vcf"
        if [[ ! -f "$IN_VCF" ]]; then
            echo "  WARNING: missing $IN_VCF"
            continue
        fi
        echo "  deepvariant / $SAMPLE"
        zcat "$IN_VCF" | grep "^#" > "$OUT_VCF"
        zcat "$IN_VCF" | awk -v pat="$STD_CTGS" '$1 !~ /^#/ && $1 !~ pat' >> "$OUT_VCF"
    done
fi

# --- Step 2: Build summary table ---
echo "=== Building FP summary table ==="
mkdir -p "$(dirname "$OUT_TSV")"

declare -A all_snps

for PLATFORM in pacbio drna ont; do
    for SAMPLE in "${SAMPLES[@]}"; do
        VCF="${PLATFORM_DIR[$PLATFORM]}/${SAMPLE}/FP_QUAL_only.vcf"
        [[ ! -f "$VCF" ]] && continue
        while read -r chrom pos id ref alt qual _; do
            snp_id="${SAMPLE}:${chrom}:${pos}:${ref}:${alt}"
            all_snps["$snp_id,$PLATFORM"]="$qual"
        done < <(grep -v "^#" "$VCF")
    done
done

for SAMPLE in "${SAMPLES[@]}"; do
    VCF="${ILL_DIR}/${SAMPLE}/FP_QUAL_only.vcf"
    [[ ! -f "$VCF" ]] && continue
    while read -r chrom pos id ref alt qual _; do
        snp_id="${SAMPLE}:${chrom}:${pos}:${ref}:${alt}"
        all_snps["$snp_id,illumina"]="$qual"
    done < <(grep -v "^#" "$VCF")
done

for SAMPLE in "${SAMPLES[@]}"; do
    VCF="${DV_DIR}/${SAMPLE}/FP_QUAL_only.vcf"
    [[ ! -f "$VCF" ]] && continue
    while read -r chrom pos id ref alt qual _; do
        snp_id="${SAMPLE}:${chrom}:${pos}:${ref}:${alt}"
        all_snps["$snp_id,deepvariant"]="$qual"
    done < <(grep -v "^#" "$VCF")
done

all_ids=($(printf "%s\n" "${!all_snps[@]}" | cut -d, -f1 | sort -u))

echo -e "SNP_ID\tpacbio\tdrna\tont\tillumina\tdeepvariant" > "$OUT_TSV"
for snp_id in "${all_ids[@]}"; do
    line="$snp_id"
    for method in pacbio drna ont illumina deepvariant; do
        key="${snp_id},${method}"
        line+="\t${all_snps[$key]:-NA}"
    done
    echo -e "$line" >> "$OUT_TSV"
done

echo "Done. FP summary: $OUT_TSV"
