# ============================================================
# Base count analysis
#
# Count total bases per FASTQ across all platforms (LR + SR) and
# produce a per-cell-line cross-platform summary TSV.
#
# Trigger: snakemake aggregate_base_counts
# ============================================================

import glob
from os.path import join as _join

# Convenience aliases derived from sub-workflow configs
_lr_cfg = sub_wf_config['lr_bulk']
_sr_cfg = sub_wf_config['sr_bulk']

_lr_results   = _lr_cfg['output_path']
_sr_results   = _sr_cfg['output_path']
_main_results = main_wf_config['output_path']
_gb_scratch   = config['scratch_dir'] + "/gb_subsample"

# ────────────────────────────────────────────────────────────
# Phase A — base counting (delegate to sub-workflow rules)
# ────────────────────────────────────────────────────────────

rule aggregate_base_counts:
    """
    Collect all per-FASTQ base-count files and write a summary TSV plus a
    suggested target bases integer.  Review the TSV, then set
    gb_subsample_target_bases in config/config_main_wf.yaml before Phase B.
    """
    input:
        lr = expand(
            _lr_results + "/qc/base_counts/{sample}_{cell_line}.total_bases",
            sample=_lr_cfg['sample_id'],
            cell_line=_lr_cfg['cell_lines']),
        sr = expand(
            _sr_results + "/qc/base_counts/{sample}_{cell_line}.total_bases",
            sample=_sr_cfg['sample_id'],
            cell_line=_sr_cfg['cell_lines'])
    output:
        summary = _main_results + "/gb_analysis/base_count_summary.tsv",
        target  = _main_results + "/gb_analysis/suggested_target_bases.txt"
    params:
        script = config['main_wf_dir'] + "/scripts/aggregate_base_counts.py"
    localrule: True
    shell:
        """
        mkdir -p $(dirname {output.summary})
        python3 {params.script} {output.summary} {output.target} {input.lr} {input.sr}
        """

# 
# # ────────────────────────────────────────────────────────────
# # Phase B — LR transcript alignment on gb-downsampled FASTQs
# # ────────────────────────────────────────────────────────────
# 
# 
# # ────────────────────────────────────────────────────────────
# # Phase B — LR FASTQ subsampling
# # ────────────────────────────────────────────────────────────
# 
# rule gb_downsample_fastq_lr:
#     """
#     Subsample a LR FASTQ to match the Illumina (ill_bulk) base count for the
#     same cell line.  No global target needed — each cell line gets its own
#     per-cell-line target derived from the SR base-count file.
#     """
#     input:
#         fastq = lambda w: glob.glob(
#             os.path.join(_lr_cfg['samples_fastq_dir'][w.sample],
#                          f"{w.cell_line}.fastq*"))[0],
#         lr_base_count = _lr_results + "/qc/base_counts/{sample}_{cell_line}.total_bases",
#         sr_base_count = _sr_results + "/qc/base_counts/" + _sr_cfg['sample_id'][0] + "_{cell_line}.total_bases"
#     output:
#         fastq = temp(_gb_scratch + "/lr/{sample}/{cell_line}.fastq.gz")
#     params:
#         seed = config['random_seed']
#     resources:
#         cpus_per_task=4,
#         mem_mb=8000
#     shell:
#         """
#         mkdir -p $(dirname {output.fastq})
#         ACTUAL=$(cat {input.lr_base_count})
#         TARGET=$(cat {input.sr_base_count})
#         FRAC=$(python3 -c "print(min(1.0, $TARGET / $ACTUAL))")
#         seqtk sample -s {params.seed} {input.fastq} $FRAC | gzip > {output.fastq}
#         """
# 
# 
# # ────────────────────────────────────────────────────────────
# # Phase B — LR transcript alignment on gb-downsampled FASTQs
# # ────────────────────────────────────────────────────────────
# 
# rule minimap2_transcript_gb_lr:
#     """Align gb-downsampled LR FASTQ to transcriptome."""
#     input:
#         fastq = _gb_scratch + "/lr/{sample}/{cell_line}.fastq.gz",
#         ref   = _lr_cfg['reference']['transcript']
#     output:
#         bam = temp(_gb_scratch + "/lr/TranscriptAlignment/{sample}_{cell_line}.bam")
#     resources:
#         cpus_per_task=16,
#         mem_mb=64000
#     params:
#         minimap2 = config['software']['minimap2'],
#         opts     = lambda w: _lr_cfg['minimap2_trans_options'][w.sample]
#     shell:
#         """
#         mkdir -p $(dirname {output.bam})
#         {params.minimap2} {params.opts} -t {resources.cpus_per_task} \
#             {input.ref} {input.fastq} | samtools view -bS > {output.bam}
#         """
# 
# 
# # ────────────────────────────────────────────────────────────
# # Phase B — LR oarfish quantification on gb-downsampled data
# # ────────────────────────────────────────────────────────────
# 
# rule oarfish_cov_gb_lr:
#     """Run oarfish (coverage model) on the gb-downsampled transcript BAM."""
#     input:
#         bam = _gb_scratch + "/lr/TranscriptAlignment/{sample}_{cell_line}.bam",
#         ref = _lr_cfg['reference']['transcript']
#     output:
#         out_dir = directory(_lr_results + "/gb_subsample/oarfish_cov_output/{sample}/{cell_line}")
#     conda:
#         _lr_cfg['conda']['oarfish']
#     resources:
#         cpus_per_task=16,
#         mem_mb=64000
#     shell:
#         """
#         mkdir -p {output.out_dir}
#         oarfish --alignments {input.bam} \
#                 --threads {resources.cpus_per_task} \
#                 --output {output.out_dir}/ \
#                 --model-coverage -d . \
#                 --filter-group no-filters \
#                 --num-bootstraps 50
#         """
# 
# 
# # ────────────────────────────────────────────────────────────
# # Phase B — SR salmon quantification on gb-downsampled data
# # ────────────────────────────────────────────────────────────
# 
# rule salmon_gb_sr:
#     """Run salmon quasi-mapping on gb-downsampled SR FASTQs."""
#     input:
#         R1    = _gb_scratch + "/sr/{sample}/{cell_line}_R1.fastq.gz",
#         R2    = _gb_scratch + "/sr/{sample}/{cell_line}_R2.fastq.gz",
#         index = _sr_results + "/salmon/index"
#     output:
#         out_dir = directory(_sr_results + "/gb_subsample/salmon_output/{sample}/{cell_line}")
#     conda:
#         _sr_cfg['conda']['main']
#     resources:
#         cpus_per_task=16,
#         mem_mb=32000
#     shell:
#         """
#         mkdir -p {output.out_dir}
#         salmon quant -i {input.index} \
#                      -l A \
#                      -1 {input.R1} \
#                      -2 {input.R2} \
#                      --validateMappings \
#                      -o {output.out_dir} \
#                      -p {resources.cpus_per_task} \
#                      --numBootstraps 50
#         """
# 
# 
# # ────────────────────────────────────────────────────────────
# # Target rules
# # ────────────────────────────────────────────────────────────
# 
# rule run_gb_downsample_lr:
#     """Run full LR gb-downsampling pipeline for all samples and cell lines."""
#     input:
#         expand(
#             rules.oarfish_cov_gb_lr.output.out_dir,
#             sample=_lr_cfg['sample_id'],
#             cell_line=_lr_cfg['cell_lines'])
#     output:
#         touch(config['flag_dir'] + "/gb_downsample_lr.done")
# 
# 
# rule run_gb_downsample_sr:
#     """Run full SR gb-downsampling pipeline for all samples and cell lines."""
#     input:
#         expand(
#             rules.salmon_gb_sr.output.out_dir,
#             sample=_sr_cfg['sample_id'],
#             cell_line=_sr_cfg['cell_lines'])
#     output:
#         touch(config['flag_dir'] + "/gb_downsample_sr.done")
