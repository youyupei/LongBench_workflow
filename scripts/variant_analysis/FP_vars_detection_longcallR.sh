#!/bin/bash
# FP detection for LongcallR VCFs: variants on non-standard contigs.
# No QUAL threshold — records QUAL as-is for downstream filtering in Rmd.
#
# Usage:
#   FP_vars_detection_longcallR.sh <out_tsv> <longcallR_dir>
#
# longcallR_dir: flat directory containing {pb_bulk|ont_bulk|dRNA_bulk}_{cell_line}.longcallR.vcf
# Output TSV columns: SNP_ID, pacbio, ont, drna  (QUAL values or NA)

set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <out_tsv> <longcallR_dir>" >&2
    exit 1
fi

OUT_TSV="$1"
LCR_DIR="$2"

SAMPLES=(H146 H1975 H211 H2228 H526 H69 HCC827 SHP77)
STD_CTGS='^(chr([1-9]|1[0-9]|2[0-2]|X|Y|M)|GL|KI)'

declare -A PLATFORM_PREFIX=(
    [pacbio]=pb_bulk
    [ont]=ont_bulk
    [drna]=dRNA_bulk
)

# --- Step 1: Extract FP candidates ---
echo "=== Extracting non-standard-contig variants ==="

declare -A all_snps

for PLATFORM in pacbio ont drna; do
    PREFIX="${PLATFORM_PREFIX[$PLATFORM]}"
    for SAMPLE in "${SAMPLES[@]}"; do
        VCF="${LCR_DIR}/${PREFIX}_${SAMPLE}.longcallR.vcf"
        if [[ ! -f "$VCF" ]]; then
            echo "  WARNING: missing $VCF"
            continue
        fi
        echo "  $PLATFORM / $SAMPLE"
        while read -r chrom pos id ref alt qual _; do
            snp_id="${SAMPLE}:${chrom}:${pos}:${ref}:${alt}"
            all_snps["$snp_id,$PLATFORM"]="$qual"
        done < <(awk -v pat="$STD_CTGS" '$1 !~ /^#/ && $1 !~ pat' "$VCF")
    done
done

# --- Step 2: Build summary table ---
echo "=== Building FP summary table ==="
mkdir -p "$(dirname "$OUT_TSV")"

all_ids=($(printf "%s\n" "${!all_snps[@]}" | cut -d, -f1 | sort -u))

echo -e "SNP_ID\tpacbio\tont\tdrna" > "$OUT_TSV"
for snp_id in "${all_ids[@]}"; do
    line="$snp_id"
    for method in pacbio ont drna; do
        key="${snp_id},${method}"
        line+="\t${all_snps[$key]:-NA}"
    done
    echo -e "$line" >> "$OUT_TSV"
done

echo "Done. FP summary: $OUT_TSV"
