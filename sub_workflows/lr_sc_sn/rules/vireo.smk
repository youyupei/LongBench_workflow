rule run_vireo:
    """
    Run demuxlet on the single cell nuclei data.
    """
    input:
        expand(config['output_path'] + '/vireo/{sample}_vireo',
                sample=config['sample_id'])
    output:
        touch(os.path.join(results_dir, ".flag/vireo.done"))

rule prepare_input_vcf:
    """
    Take the mutation table and prepare the input VCF file for the cell lines.
    """
    input:
        config['reference']['cellline_mutations']
    output:
        os.path.join(results_dir, "demuxlet/DepmapMutTabConverted.vcf")
    resources:
        cpus_per_task=1,
        mem_mb=32000
    params:
        script = os.path.join(config['sub_wf_dir'], 'scripts/DepmapMutTab2Vcf.py')
    shell:
        """
        OUT=$(dirname {output})/
        mkdir -p $OUT
        python3 {params.script} {input} {output}
        """



rule get_vcf_from_sr:
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
        BAM1=/vast/projects/LongBench/analysis/sr_bulk/result/subjunc/bam/H146.sorted.bam
        BAM2=/vast/projects/LongBench/analysis/sr_bulk/result/subjunc/bam/H1975.sorted.bam
        BAM3=/vast/projects/LongBench/analysis/sr_bulk/result/subjunc/bam/H211.sorted.bam
        BAM4=/vast/projects/LongBench/analysis/sr_bulk/result/subjunc/bam/H2228.sorted.bam
        BAM5=/vast/projects/LongBench/analysis/sr_bulk/result/subjunc/bam/H526.sorted.bam
        BAM6=/vast/projects/LongBench/analysis/sr_bulk/result/subjunc/bam/H69.sorted.bam
        BAM7=/vast/projects/LongBench/analysis/sr_bulk/result/subjunc/bam/HCC827.sorted.bam
        BAM8=/vast/projects/LongBench/analysis/sr_bulk/result/subjunc/bam/SHP77.sorted.bam

        cellsnp-lite -s $BAM1,$BAM2,$BAM3,$BAM4,$BAM5,$BAM6,$BAM7,$BAM8 \
            -I H146,H1975,H211,H2228,H526,H69,HCC827,SHP77 \
            -O /home/users/allstaff/you.yu/vast/test/test_cellsnp_bulk/test_common_SNP \
            -R /home/users/allstaff/you.yu/LongBench/reference_files/GRCh38/commonSNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
            -p 32 \
            --cellTAG None --UMItag None --gzip --genotype 

        gunzip /home/users/allstaff/you.yu/vast/test/test_cellsnp_bulk/test_common_SNP/cellSNP.cells.vcf.gz
        """

rule Cellsnp_lite_rule:
    input:
        bam= config['output_path'] + '/flames_out/{sample}/align2genome.bam',
        bam_sorted_idx= config['output_path'] + '/flames_out/{sample}/align2genome.bam.bai',
        bc_list = config['output_path'] + '/flames_out/{sample}/whitelist.csv',
        #vcf = config['output_path'] + '/demuxlet/DepmapMutTabConverted.vcf'
        vcf = '/home/users/allstaff/you.yu/vast/test/test_cellsnp_bulk/test_common_SNP/cellSNP.cells.vcf'
    output:
        directory(config['output_path'] + '/vireo/{sample}_cellsnp')
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
        #vcf = config['output_path'] + '/demuxlet/DepmapMutTabConverted.vcf'
        vcf = '/home/users/allstaff/you.yu/vast/test/test_cellsnp_bulk/test_common_SNP/cellSNP.cells.vcf'
    output:
        directory(config['output_path'] + '/vireo/{sample}_vireo')
    resources:
        cpus_per_task=16,
        mem_mb=50000,
        slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
    params:
        demuxlet_excutable = config['software']['demuxlet']
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