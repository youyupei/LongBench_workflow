# ============================================================
# Downsampled variant analysis (clair3 + LongcallR)
# Inherits _CCLE_SOURCE, _CCLE_DIR, _ill_dir, _dv_dir,
# _cell_lines from the main variant_analysis.smk via include chain.
#
# Trigger clair3 ds:      snakemake variant_analysis_ds
# Trigger LongcallR ds:   snakemake variant_analysis_longcallR_ds
# ============================================================

_CCLE_SOURCE = "/vast/projects/LongBench/analysis/illumina_sc_sn/result/demuxlet/DepmapMutTabConverted.vcf"
_CCLE_DIR    = "/vast/projects/LongBench/analysis/variant_analysis/CCLE"
_lr_ds_scratch  = "/vast/scratch/users/you.yu/LongBench/lr_bulk_ds"
_pacbio_ds_dir  = _lr_ds_scratch + "/Mutation/clair3_rna/pb_bulk"
_ont_ds_dir     = _lr_ds_scratch + "/Mutation/clair3_rna/ont_bulk"
_drna_ds_dir    = _lr_ds_scratch + "/Mutation/clair3_rna/dRNA_bulk"
_VA_DS_DIR      = "/vast/projects/LongBench/analysis/variant_analysis/bcftools_ds"
_cell_lines = config['cell_lines']

_longcallR_ds_dir = "/vast/projects/LongBench/analysis/lr_bulk_ds/result/LongcallR_with_spike_in"
_VA_LCR_DS_DIR    = "/vast/projects/LongBench/analysis/variant_analysis/bcftools_longcallR_ds"


_ill_dir    =  "/vast/scratch/users/you.yu/LongBench/sr_bulk/Mutation/clair3"
_dv_dir     =  "/vast/projects/LongBench/analysis/sr_bulk/result/Deepvariants"

# ── clair3 downsampled ───────────────────────────────────────

rule variant_tp_ds:
    input:
        ccle_source = _CCLE_SOURCE,
        pacbio_vcfs = expand(_pacbio_ds_dir + "/{cell_line}/output.vcf.gz",      cell_line=_cell_lines),
        ont_vcfs    = expand(_ont_ds_dir    + "/{cell_line}/output.vcf.gz",      cell_line=_cell_lines),
        drna_vcfs   = expand(_drna_ds_dir   + "/{cell_line}/output.vcf.gz",      cell_line=_cell_lines)
    output:
        tsvs = expand(_VA_DS_DIR + "/genotype_comparisons_qual/genotype_comparison_{cell_line}.tsv", cell_line=_cell_lines)
    params:
        script     = config['main_wf_dir'] + "/scripts/variant_analysis/CCLE_vars_detection_youyu.sh",
        ccle_dir   = _CCLE_DIR,
        outdir     = _VA_DS_DIR,
        pacbio_dir = _pacbio_ds_dir,
        ont_dir    = _ont_ds_dir,
        drna_dir   = _drna_ds_dir,
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


rule variant_tp_phased_ds:
    input:
        ccle_source = _CCLE_SOURCE,
        pacbio_vcfs = expand(_pacbio_ds_dir + "/{cell_line}/output_enable_phasing.vcf.gz", cell_line=_cell_lines),
        ont_vcfs    = expand(_ont_ds_dir    + "/{cell_line}/output_enable_phasing.vcf.gz", cell_line=_cell_lines),
        drna_vcfs   = expand(_drna_ds_dir   + "/{cell_line}/output_enable_phasing.vcf.gz", cell_line=_cell_lines),
        ill_vcfs    = expand(_ill_dir       + "/{cell_line}/merge_output.vcf.gz",           cell_line=_cell_lines),
        dv_vcfs     = expand(_dv_dir        + "/{cell_line}/output.vcf.gz",                cell_line=_cell_lines),
    output:
        tsvs = expand(_VA_DS_DIR + "_phased/genotype_comparisons_qual/genotype_comparison_{cell_line}.tsv", cell_line=_cell_lines)
    params:
        script     = config['main_wf_dir'] + "/scripts/variant_analysis/CCLE_vars_detection_youyu.sh",
        ccle_dir   = _CCLE_DIR,
        outdir     = _VA_DS_DIR + "_phased",
        pacbio_dir = _pacbio_ds_dir,
        ont_dir    = _ont_ds_dir,
        drna_dir   = _drna_ds_dir,
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


rule variant_fp_ds:
    input:
        pacbio_vcfs = expand(_pacbio_ds_dir + "/{cell_line}/output.vcf.gz",      cell_line=_cell_lines),
        ont_vcfs    = expand(_ont_ds_dir    + "/{cell_line}/output.vcf.gz",      cell_line=_cell_lines),
        drna_vcfs   = expand(_drna_ds_dir   + "/{cell_line}/output.vcf.gz",      cell_line=_cell_lines),
        ill_vcfs    = expand(_ill_dir       + "/{cell_line}/merge_output.vcf.gz", cell_line=_cell_lines),
        dv_vcfs     = expand(_dv_dir        + "/{cell_line}/output.vcf.gz",      cell_line=_cell_lines),
    output:
        tsv = _VA_DS_DIR + "/SNP_FP_QUAL_detail.tsv"
    params:
        script     = config['main_wf_dir'] + "/scripts/variant_analysis/FP_vars_detection_youyu.sh",
        pacbio_dir = _pacbio_ds_dir,
        ont_dir    = _ont_ds_dir,
        drna_dir   = _drna_ds_dir,
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


rule variant_fp_phased_ds:
    input:
        pacbio_vcfs = expand(_pacbio_ds_dir + "/{cell_line}/output_enable_phasing.vcf.gz", cell_line=_cell_lines),
        ont_vcfs    = expand(_ont_ds_dir    + "/{cell_line}/output_enable_phasing.vcf.gz", cell_line=_cell_lines),
        drna_vcfs   = expand(_drna_ds_dir   + "/{cell_line}/output_enable_phasing.vcf.gz", cell_line=_cell_lines),
        ill_vcfs    = expand(_ill_dir       + "/{cell_line}/merge_output.vcf.gz",           cell_line=_cell_lines),
        dv_vcfs     = expand(_dv_dir        + "/{cell_line}/output.vcf.gz",                cell_line=_cell_lines),
    output:
        tsv = _VA_DS_DIR + "_phased/SNP_FP_QUAL_detail.tsv"
    params:
        script     = config['main_wf_dir'] + "/scripts/variant_analysis/FP_vars_detection_youyu.sh",
        pacbio_dir = _pacbio_ds_dir,
        ont_dir    = _ont_ds_dir,
        drna_dir   = _drna_ds_dir,
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


rule variant_analysis_ds:
    input:
        rules.variant_tp_ds.output,
        rules.variant_fp_ds.output,
        rules.variant_tp_phased_ds.output,
        rules.variant_fp_phased_ds.output
    output:
        touch(os.path.join(config['flag_dir'], "variant_analysis_ds.done"))


# ── LongcallR downsampled ────────────────────────────────────

rule variant_tp_longcallR_ds:
    """
    TP analysis for downsampled LongcallR VCFs.
    Output: {_VA_LCR_DS_DIR}/genotype_comparisons_qual/genotype_comparison_{cell_line}.tsv
    """
    input:
        ccle_source  = _CCLE_SOURCE,
        pacbio_vcfs  = expand(_longcallR_ds_dir + "/pb_bulk_{cell_line}.longcallR.vcf",   cell_line=_cell_lines),
        ont_vcfs     = expand(_longcallR_ds_dir + "/ont_bulk_{cell_line}.longcallR.vcf",  cell_line=_cell_lines),
        drna_vcfs    = expand(_longcallR_ds_dir + "/dRNA_bulk_{cell_line}.longcallR.vcf", cell_line=_cell_lines),
    output:
        tsvs = expand(_VA_LCR_DS_DIR + "/genotype_comparisons_qual/genotype_comparison_{cell_line}.tsv", cell_line=_cell_lines)
    params:
        script        = config['main_wf_dir'] + "/scripts/variant_analysis/CCLE_vars_detection_longcallR.sh",
        ccle_dir      = _CCLE_DIR,
        outdir        = _VA_LCR_DS_DIR,
        longcallR_dir = _longcallR_ds_dir,
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


rule variant_fp_longcallR_ds:
    """
    FP analysis for downsampled LongcallR VCFs.
    Output: {_VA_LCR_DS_DIR}/SNP_FP_QUAL_detail.tsv
    """
    input:
        pacbio_vcfs  = expand(_longcallR_ds_dir + "/pb_bulk_{cell_line}.longcallR.vcf",   cell_line=_cell_lines),
        ont_vcfs     = expand(_longcallR_ds_dir + "/ont_bulk_{cell_line}.longcallR.vcf",  cell_line=_cell_lines),
        drna_vcfs    = expand(_longcallR_ds_dir + "/dRNA_bulk_{cell_line}.longcallR.vcf", cell_line=_cell_lines),
    output:
        tsv = _VA_LCR_DS_DIR + "/SNP_FP_QUAL_detail.tsv"
    params:
        script        = config['main_wf_dir'] + "/scripts/variant_analysis/FP_vars_detection_longcallR.sh",
        longcallR_dir = _longcallR_ds_dir,
    resources:
        cpus_per_task = 2,
        mem_mb        = 4000
    shell:
        """
        bash {params.script} \
            {output.tsv} \
            {params.longcallR_dir}
        """


rule variant_analysis_longcallR_ds:
    input:
        rules.variant_tp_longcallR_ds.output,
        rules.variant_fp_longcallR_ds.output,
    output:
        touch(os.path.join(config['flag_dir'], "variant_analysis_longcallR_ds.done"))
