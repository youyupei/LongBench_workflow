
# Detection of FP
# First, for all methods, create a FP_QUAL_only.vcf containing only variants with QUAL ≥ 10 & located on contigs that are not the ‘standard’ ones. 

# illumina Clair3
samples=(H146 H1975 H211 H2228 H526 H69 HCC827 SHP77)
dir="/vast/scratch/users/david.m/LongBench/Clair3-illumina/spikeins/all_contigs"

for sample in "${samples[@]}"; do
    in_vcf="${dir}/${sample}/merge_output.vcf.gz"
    out_vcf="${dir}/${sample}/FP_QUAL_only.vcf"

    zcat "$in_vcf" | grep "^#" > "$out_vcf"
    zcat "$in_vcf" | awk '$1 !~ /^#/ && $6 >= 10 && $1 !~ /^(chr([1-9]|1[0-9]|2[0-2]|X|Y|M)|GL|KI)/' >> "$out_vcf"
done

# drna
for sample in "${samples[@]}"; do
  in_vcf="/vast/scratch/users/david.m/LongBench/Clair3-RNA/spikeins/all_contigs/drna/${sample}/output.vcf.gz"
  out_vcf="/vast/scratch/users/david.m/LongBench/Clair3-RNA/spikeins/all_contigs/drna/${sample}/FP_QUAL_only.vcf"
  zcat "$in_vcf" | grep "^#" > "$out_vcf"
  zcat "$in_vcf" | awk '$1 !~ /^#/ && $6 >= 10 && $1 !~ /^(chr([1-9]|1[0-9]|2[0-2]|X|Y|M)|GL|KI)/' >> "$out_vcf"
done

# pacbio
for sample in "${samples[@]}"; do
  in_vcf="/vast/scratch/users/david.m/LongBench/Clair3-RNA/spikeins/all_contigs/pacbio/${sample}/output.vcf.gz"
  out_vcf="/vast/scratch/users/david.m/LongBench/Clair3-RNA/spikeins/all_contigs/pacbio/${sample}/FP_QUAL_only.vcf"
  zcat "$in_vcf" | grep "^#" > "$out_vcf"
  zcat "$in_vcf" | awk '$1 !~ /^#/ && $6 >= 10 && $1 !~ /^(chr([1-9]|1[0-9]|2[0-2]|X|Y|M)|GL|KI)/' >> "$out_vcf"
done

# ont
for sample in "${samples[@]}"; do
  in_vcf="/vast/scratch/users/david.m/LongBench/Clair3-RNA/spikeins/all_contigs/ont/${sample}/output.vcf.gz"
  out_vcf="/vast/scratch/users/david.m/LongBench/Clair3-RNA/spikeins/all_contigs/ont/${sample}/FP_QUAL_only.vcf"
  zcat "$in_vcf" | grep "^#" > "$out_vcf"
  zcat "$in_vcf" | awk '$1 !~ /^#/ && $6 >= 10 && $1 !~ /^(chr([1-9]|1[0-9]|2[0-2]|X|Y|M)|GL|KI)/' >> "$out_vcf"
done

# DV
for sample in "${samples[@]}"; do
  in_vcf="/vast/scratch/users/david.m/LongBench/Deepvariants/${sample}/output.vcf.gz"
  out_vcf="/vast/scratch/users/david.m/LongBench/Deepvariants/${sample}/FP_QUAL_only.vcf"
  zcat "$in_vcf" | grep "^#" > "$out_vcf"
  zcat "$in_vcf" | awk '$1 !~ /^#/ && $6 >= 10 && $1 !~ /^(chr([1-9]|1[0-9]|2[0-2]|X|Y|M)|GL|KI)/' >> "$out_vcf"
done

# Next, generate a summary table (containing sample info!)
declare -A all_snps

# Loop over Clair3-RNA (LR) methods
for method in pacbio drna ont; do
    for vcf in /vast/scratch/users/david.m/LongBench/Clair3-RNA/spikeins/all_contigs/${method}/*/FP_QUAL_only.vcf; do
        sample=$(basename $(dirname "$vcf"))
        while read -r chrom pos id ref alt _; do
            snp_id="${sample}:${chrom}:${pos}:${ref}:${alt}"
            all_snps["$snp_id,$method"]=TRUE
        done < <(grep -v "^#" "$vcf")
    done
done

# Clair3-illumina
for vcf in /vast/scratch/users/david.m/LongBench/Clair3-illumina/spikeins/all_contigs/*/FP_QUAL_only.vcf; do
    sample=$(basename $(dirname "$vcf"))
    method=illumina
    while read -r chrom pos id ref alt _; do
        snp_id="${sample}:${chrom}:${pos}:${ref}:${alt}"
        all_snps["$snp_id,$method"]=TRUE
    done < <(grep -v "^#" "$vcf")
done

# DeepVariant
for vcf in /vast/scratch/users/david.m/LongBench/Deepvariants/*/FP_QUAL_only.vcf; do
    sample=$(basename $(dirname "$vcf"))
    method=deepvariant
    while read -r chrom pos id ref alt _; do
        snp_id="${sample}:${chrom}:${pos}:${ref}:${alt}"
        all_snps["$snp_id,$method"]=TRUE
    done < <(grep -v "^#" "$vcf")
done

# Get all unique SNP IDs (WITH sample name)
all_ids=($(printf "%s\n" "${!all_snps[@]}" | cut -d, -f1 | sort -u))

# Header
echo -e "SNP_ID\tpacbio\tdrna\tont\tillumina\tdeepvariant" > SNP_FP_QUAL_detail.tsv

# Create table
for snp_id in "${all_ids[@]}"; do
    line="$snp_id"
    for method in pacbio drna ont illumina deepvariant; do
        key="${snp_id},${method}"
        [[ "${all_snps[$key]}" == "TRUE" ]] && line+="\tTRUE" || line+="\tFALSE"
    done
    echo -e "$line" >> SNP_FP_QUAL_detail.tsv
done
