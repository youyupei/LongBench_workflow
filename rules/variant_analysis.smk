# ============================================================
# Variant analysis: TP (CCLE intersection) + FP (spike-in contigs)
#
# Trigger TP analysis:   snakemake variant_tp_all
# Trigger FP analysis:   snakemake variant_fp_all
# Trigger both:          snakemake variant_analysis_all
#
# LongcallR analysis:
# Trigger TP:            snakemake variant_tp_longcallR
# Trigger FP:            snakemake variant_fp_longcallR
# Trigger both:          snakemake variant_analysis_longcallR
# ============================================================

_lr_scratch = sub_wf_config['lr_bulk']['scratch_dir']
_sr_scratch = sub_wf_config['sr_bulk']['scratch_dir']
_sr_results = sub_wf_config['sr_bulk']['output_path']
_cell_lines = sub_wf_config['lr_bulk']['cell_lines']

_CCLE_SOURCE = "/vast/projects/LongBench/analysis/illumina_sc_sn/result/demuxlet/DepmapMutTabConverted.vcf"
_CCLE_DIR    = "/vast/projects/LongBench/analysis/variant_analysis/CCLE"
_VA_DIR      = "/vast/projects/LongBench/analysis/variant_analysis/bcftools"

_pacbio_dir = _lr_scratch + "/Mutation/clair3_rna/pb_bulk"
_ont_dir    = _lr_scratch + "/Mutation/clair3_rna/ont_bulk"
_drna_dir   = _lr_scratch + "/Mutation/clair3_rna/dRNA_bulk"
_ill_dir    = _sr_scratch + "/Mutation/clair3"
_dv_dir     = _sr_results + "/Deepvariants"



rule variant_tp_all:
    """
    TP analysis: intersect each platform's Clair3 VCFs against CCLE ground truth,
    then build a per-cell-line genotype comparison table with match/discordant/not_called status.
    Output: {_VA_DIR}/genotype_comparisons_qual/genotype_comparison_{{cell_line}}.tsv
    """
    input:
        ccle_source = _CCLE_SOURCE,
        pacbio_vcfs = expand(_pacbio_dir + "/{cell_line}/output.vcf.gz",      cell_line=_cell_lines),
        ont_vcfs    = expand(_ont_dir    + "/{cell_line}/output.vcf.gz",      cell_line=_cell_lines),
        drna_vcfs   = expand(_drna_dir   + "/{cell_line}/output.vcf.gz",      cell_line=_cell_lines),
        ill_vcfs    = expand(_ill_dir    + "/{cell_line}/merge_output.vcf.gz", cell_line=_cell_lines),
        dv_vcfs     = expand(_dv_dir     + "/{cell_line}/output.vcf.gz",      cell_line=_cell_lines),
    output:
        tsvs = expand(_VA_DIR + "/genotype_comparisons_qual/genotype_comparison_{cell_line}.tsv", cell_line=_cell_lines)
    params:
        script     = config['main_wf_dir'] + "/scripts/variant_analysis/CCLE_vars_detection_youyu.sh",
        ccle_dir   = _CCLE_DIR,
        outdir     = _VA_DIR,
        pacbio_dir = _pacbio_dir,
        ont_dir    = _ont_dir,
        drna_dir   = _drna_dir,
        ill_dir    = _ill_dir,
        dv_dir     = _dv_dir,
    resources:
        cpus_per_task = 2,
        mem_mb        = 8000
    shell:
        """
        bash {params.script} \
            {input.ccle_source} \
            {params.ccle_dir} \
            {params.outdir} \
            {params.pacbio_dir} \
            {params.ont_dir} \
            {params.drna_dir} \
            {params.ill_dir} \
            {params.dv_dir}
        """


rule variant_tp_phased:
    """
    TP analysis using phased LR VCFs (output_enable_phasing.vcf.gz from clair3_rna).
    Illumina and DeepVariant use the same unphased VCFs as variant_tp_all.
    Output: {_VA_DIR}_phased/genotype_comparisons_qual/genotype_comparison_{{cell_line}}.tsv
    """
    input:
        ccle_source = _CCLE_SOURCE,
        pacbio_vcfs = expand(_pacbio_dir + "/{cell_line}/output_enable_phasing.vcf.gz", cell_line=_cell_lines),
        ont_vcfs    = expand(_ont_dir    + "/{cell_line}/output_enable_phasing.vcf.gz", cell_line=_cell_lines),
        drna_vcfs   = expand(_drna_dir   + "/{cell_line}/output_enable_phasing.vcf.gz", cell_line=_cell_lines),
        ill_vcfs    = expand(_ill_dir    + "/{cell_line}/merge_output.vcf.gz",           cell_line=_cell_lines),
        dv_vcfs     = expand(_dv_dir     + "/{cell_line}/output.vcf.gz",                cell_line=_cell_lines),
    output:
        tsvs = expand(_VA_DIR + "_phased/genotype_comparisons_qual/genotype_comparison_{cell_line}.tsv", cell_line=_cell_lines)
    params:
        script     = config['main_wf_dir'] + "/scripts/variant_analysis/CCLE_vars_detection_youyu.sh",
        ccle_dir   = _CCLE_DIR,
        outdir     = _VA_DIR + "_phased",
        pacbio_dir = _pacbio_dir,
        ont_dir    = _ont_dir,
        drna_dir   = _drna_dir,
        ill_dir    = _ill_dir,
        dv_dir     = _dv_dir,
    resources:
        cpus_per_task = 2,
        mem_mb        = 8000
    shell:
        """
        bash {params.script} \
            {input.ccle_source} \
            {params.ccle_dir} \
            {params.outdir} \
            {params.pacbio_dir} \
            {params.ont_dir} \
            {params.drna_dir} \
            {params.ill_dir} \
            {params.dv_dir} \
            output_enable_phasing.vcf.gz
        """


rule variant_fp_all:
    """
    FP analysis: extract QUAL>=10 variants on non-standard contigs per platform/sample,
    then build a combined presence/absence summary table.
    Output: {_VA_DIR}/SNP_FP_QUAL_detail.tsv
    """
    input:
        pacbio_vcfs = expand(_pacbio_dir + "/{cell_line}/output.vcf.gz",      cell_line=_cell_lines),
        ont_vcfs    = expand(_ont_dir    + "/{cell_line}/output.vcf.gz",      cell_line=_cell_lines),
        drna_vcfs   = expand(_drna_dir   + "/{cell_line}/output.vcf.gz",      cell_line=_cell_lines),
        ill_vcfs    = expand(_ill_dir    + "/{cell_line}/merge_output.vcf.gz", cell_line=_cell_lines),
        dv_vcfs     = expand(_dv_dir     + "/{cell_line}/output.vcf.gz",      cell_line=_cell_lines),
    output:
        tsv = _VA_DIR + "/SNP_FP_QUAL_detail.tsv"
    params:
        script     = config['main_wf_dir'] + "/scripts/variant_analysis/FP_vars_detection_youyu.sh",
        pacbio_dir = _pacbio_dir,
        ont_dir    = _ont_dir,
        drna_dir   = _drna_dir,
        ill_dir    = _ill_dir,
        dv_dir     = _dv_dir,
    resources:
        cpus_per_task = 2,
        mem_mb        = 4000
    shell:
        """
        bash {params.script} \
            {output.tsv} \
            {params.pacbio_dir} \
            {params.ont_dir} \
            {params.drna_dir} \
            {params.ill_dir} \
            {params.dv_dir}
        """


rule variant_fp_phased:
    """
    FP analysis using phased LR VCFs (output_enable_phasing.vcf.gz from clair3_rna).
    Illumina and DeepVariant use the same unphased VCFs as variant_fp_all.
    Output: {_VA_DIR}_phased/SNP_FP_QUAL_detail.tsv
    """
    input:
        pacbio_vcfs = expand(_pacbio_dir + "/{cell_line}/output_enable_phasing.vcf.gz", cell_line=_cell_lines),
        ont_vcfs    = expand(_ont_dir    + "/{cell_line}/output_enable_phasing.vcf.gz", cell_line=_cell_lines),
        drna_vcfs   = expand(_drna_dir   + "/{cell_line}/output_enable_phasing.vcf.gz", cell_line=_cell_lines),
        ill_vcfs    = expand(_ill_dir    + "/{cell_line}/merge_output.vcf.gz",           cell_line=_cell_lines),
        dv_vcfs     = expand(_dv_dir     + "/{cell_line}/output.vcf.gz",                cell_line=_cell_lines),
    output:
        tsv = _VA_DIR + "_phased/SNP_FP_QUAL_detail.tsv"
    params:
        script     = config['main_wf_dir'] + "/scripts/variant_analysis/FP_vars_detection_youyu.sh",
        pacbio_dir = _pacbio_dir,
        ont_dir    = _ont_dir,
        drna_dir   = _drna_dir,
        ill_dir    = _ill_dir,
        dv_dir     = _dv_dir,
    resources:
        cpus_per_task = 2,
        mem_mb        = 4000
    shell:
        """
        bash {params.script} \
            {output.tsv} \
            {params.pacbio_dir} \
            {params.ont_dir} \
            {params.drna_dir} \
            {params.ill_dir} \
            {params.dv_dir} \
            output_enable_phasing.vcf.gz
        """


rule variant_analysis_all:
    input:
        rules.variant_tp_all.output,
        rules.variant_fp_all.output,
        rules.variant_tp_phased.output,
        rules.variant_fp_phased.output
    output:
        touch(os.path.join(config['flag_dir'], "variant_analysis_all.done"))


# ────────────────────────────────────────────────────────────
# LongcallR variant analysis (pb/ont/dRNA only, no QUAL threshold)
# ────────────────────────────────────────────────────────────

_longcallR_dir = sub_wf_config['lr_bulk']['output_path'] + "/LongcallR"
_VA_LCR_DIR    = "/vast/projects/LongBench/analysis/variant_analysis/bcftools_longcallR"

rule variant_tp_longcallR:
    """
    TP analysis for LongcallR VCFs: intersect pb/ont/dRNA against CCLE ground truth.
    Output: {_VA_LCR_DIR}/genotype_comparisons_qual/genotype_comparison_{cell_line}.tsv
    """
    input:
        ccle_source  = _CCLE_SOURCE,
        pacbio_vcfs  = expand(_longcallR_dir + "/pb_bulk_{cell_line}.longcallR.vcf",   cell_line=_cell_lines),
        ont_vcfs     = expand(_longcallR_dir + "/ont_bulk_{cell_line}.longcallR.vcf",  cell_line=_cell_lines),
        drna_vcfs    = expand(_longcallR_dir + "/dRNA_bulk_{cell_line}.longcallR.vcf", cell_line=_cell_lines),
    output:
        tsvs = expand(_VA_LCR_DIR + "/genotype_comparisons_qual/genotype_comparison_{cell_line}.tsv", cell_line=_cell_lines)
    params:
        script      = config['main_wf_dir'] + "/scripts/variant_analysis/CCLE_vars_detection_longcallR.sh",
        ccle_dir    = _CCLE_DIR,
        outdir      = _VA_LCR_DIR,
        longcallR_dir = _longcallR_dir,
    resources:
        cpus_per_task = 2,
        mem_mb        = 8000
    shell:
        """
        bash {params.script} \
            {input.ccle_source} \
            {params.ccle_dir} \
            {params.outdir} \
            {params.longcallR_dir}
        """


rule variant_fp_longcallR:
    """
    FP analysis for LongcallR VCFs: variants on non-standard contigs, QUAL recorded as-is.
    Output: {_VA_LCR_DIR}/SNP_FP_QUAL_detail.tsv
    """
    input:
        pacbio_vcfs  = expand(_longcallR_dir + "/pb_bulk_{cell_line}.longcallR.vcf",   cell_line=_cell_lines),
        ont_vcfs     = expand(_longcallR_dir + "/ont_bulk_{cell_line}.longcallR.vcf",  cell_line=_cell_lines),
        drna_vcfs    = expand(_longcallR_dir + "/dRNA_bulk_{cell_line}.longcallR.vcf", cell_line=_cell_lines),
    output:
        tsv = _VA_LCR_DIR + "/SNP_FP_QUAL_detail.tsv"
    params:
        script        = config['main_wf_dir'] + "/scripts/variant_analysis/FP_vars_detection_longcallR.sh",
        longcallR_dir = _longcallR_dir,
    resources:
        cpus_per_task = 2,
        mem_mb        = 4000
    shell:
        """
        bash {params.script} \
            {output.tsv} \
            {params.longcallR_dir}
        """


rule variant_analysis_longcallR:
    input:
        rules.variant_tp_longcallR.output,
        rules.variant_fp_longcallR.output,
    output:
        touch(os.path.join(config['flag_dir'], "variant_analysis_longcallR.done"))
