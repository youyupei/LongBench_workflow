# annotation cell lines of SC / SN short and long read data
# 1. from SR, get the cell line SNP annotation
from os.path import join
lr_config = config_parser.load_configfile(join(config['main_wf_dir'],config["sub_wf_config"]["lr_sc_sn"]))
sr_config = config_parser.load_configfile(join(config['main_wf_dir'],config["sub_wf_config"]["sr_sc_sn"]))
main_conda = join(config['main_wf_dir'], config['conda_config']['main'])

rule get_vcf_from_sr:
    input:
        BAM1='/vast/projects/LongBench/analysis/sr_bulk/result/subjunc/bam/H146.sorted.bam',
        BAM2='/vast/projects/LongBench/analysis/sr_bulk/result/subjunc/bam/H1975.sorted.bam',
        BAM3='/vast/projects/LongBench/analysis/sr_bulk/result/subjunc/bam/H211.sorted.bam',
        BAM4='/vast/projects/LongBench/analysis/sr_bulk/result/subjunc/bam/H2228.sorted.bam',
        BAM5='/vast/projects/LongBench/analysis/sr_bulk/result/subjunc/bam/H526.sorted.bam',
        BAM6='/vast/projects/LongBench/analysis/sr_bulk/result/subjunc/bam/H69.sorted.bam',
        BAM7='/vast/projects/LongBench/analysis/sr_bulk/result/subjunc/bam/HCC827.sorted.bam',
        BAM8='/vast/projects/LongBench/analysis/sr_bulk/result/subjunc/bam/SHP77.sorted.bam'
    output:
        '/home/users/allstaff/you.yu/vast/test/test_cellsnp_bulk/test_common_SNP/cellSNP.cells.vcf'
    resources:
        cpus_per_task=32,
        mem_mb=32000
    conda:
        main_conda
    shell:
        """
        mkdir -p /home/users/allstaff/you.yu/vast/test/test_cellsnp_bulk/test_common_SNP
        cellsnp-lite -s {input.BAM1},{input.BAM2},{input.BAM3},{input.BAM4},{input.BAM5},{input.BAM6},{input.BAM7},{input.BAM8} \
            -I H146,H1975,H211,H2228,H526,H69,HCC827,SHP77 \
            -O /home/users/allstaff/you.yu/vast/test/test_cellsnp_bulk/test_common_SNP \
            -R /home/users/allstaff/you.yu/LongBench/reference_files/GRCh38/commonSNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
            -p 32 \
            --cellTAG None --UMItag None --gzip --genotype 

        gunzip /home/users/allstaff/you.yu/vast/test/test_cellsnp_bulk/test_common_SNP/cellSNP.cells.vcf.gz
        """

rule Cellsnp_lite_rule:
    input:
        bam= join(lr_config['output_path'], 'flames_out/{sample}/align2genome.bam'),
        bam_sorted_idx= join(lr_config['output_path'], 'flames_out/{sample}/align2genome.bam.bai'),
        bc_list = join(lr_config['output_path'], 'flames_out/{sample}/whitelist.csv'),
        vcf = '/home/users/allstaff/you.yu/vast/test/test_cellsnp_bulk/test_common_SNP/cellSNP.cells.vcf'
    output:
        directory(join(lr_config['output_path'],'/vireo/{sample}_cellsnp'))
    resources:
        cpus_per_task=32,
        mem_mb=32000
    conda:
        main_conda
    shell:
        """
        cellsnp-lite \
          -s {input.bam} \
          -b {input.bc_list} \
          -O {output} \
          -R {input.vcf} \
          -p {resources.cpus_per_task} \
          --minMAF 0.1 \
          --minCOUNT 10 \
          --gzip \
          --genotype 
        """

rule vireo_rule:
    """
    Run the demultiplexing tool on the input VCF file.
    """
    input:
        cellsnp_rst = rules.Cellsnp_lite_rule.output,
        vcf = '/home/users/allstaff/you.yu/vast/test/test_cellsnp_bulk/test_common_SNP/cellSNP.cells.vcf'
    output:
        directory(lr_config['output_path'] + '/vireo/{sample}_vireo')
    resources:
        cpus_per_task=16,
        mem_mb=50000,
        slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
    conda:
        main_conda
    shell:
        """
        vireo \
            -c {input.cellsnp_rst} \
            -d {input.vcf} \
            -o {output} \
            -t GT \
            --callAmbientRNAs \
            -p {resources.cpus_per_task}
        """

# SR
use rule Cellsnp_lite_rule as Cellsnp_lite_rule_short_read with:
    input:
        bam= join(sr_config['output_path'], 'cellranger/{sample}/outs/possorted_genome_bam.bam'),
        bam_sorted_idx= join(sr_config['output_path'], 'cellranger/{sample}/outs/possorted_genome_bam.bam.bai'),
        bc_list = join(sr_config['output_path'], 'cellranger/{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz'),
        vcf = '/home/users/allstaff/you.yu/vast/test/test_cellsnp_bulk/test_common_SNP/cellSNP.cells.vcf'
    output:
        directory(join(sr_config['output_path'],'vireo/{sample}_cellsnp'))

use rule vireo_rule as vireo_rule_short_read with:
    input:
        cellsnp_rst = rules.Cellsnp_lite_rule_short_read.output,
        vcf = '/home/users/allstaff/you.yu/vast/test/test_cellsnp_bulk/test_common_SNP/cellSNP.cells.vcf'
    output:
        directory(join(sr_config['output_path'], 'vireo/{sample}_vireo'))


# This is the target pipeline for the workflow
rule run_cell_line_annotation_pipeline:
    """
    Run demuxlet on the single cell nuclei data.
    """
    input:
        expand(join(sr_config['output_path'], 'vireo/{sample}_vireo'),
                sample=sr_config['sample_id']),
        expand(join(lr_config['output_path'], 'vireo/{sample}_vireo'),
                sample=lr_config['sample_id'])
    output:
        touch(os.path.join(config['flag_dir'], "lr_sr_vireo.done"))

del lr_config
del sr_config