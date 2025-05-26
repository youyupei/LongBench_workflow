rule _pseudobulk_read_count_single_run:
    input:
        fastq = results_dir + "/PseudoBulkAlignment/{sample}_{cell_line}_pseudo_bulk.fastq"
    output:
        txt = temp(results_dir + "/PseudoBulkQC/{sample}_{cell_line}_pseudo_bulk_read_count.txt")
    resources:
        cpus_per_task=1
    shell:
        """
        mkdir -p $(dirname {output.txt})
        expr $(( $(wc -l < {input.fastq}) / 4 )) > {output.txt}
        """

rule pseudobulk_read_count:
    input:
        expand(rules._pseudobulk_read_count_single_run.output[0], 
            cell_line = config['cell_line_list'],
            sample = config['sample_id']) 
    output:
        results_dir + "/PseudoBulkQC/pseudo_bulk_read_count.csv"
    resources:
        cpus_per_task=1
    script:
        "../scripts/combine_read_count.R"


rule pseudobulk_Filter_and_Subsample_bam_for_qc:
    input:
        join(results_dir, "PseudoBulkAlignment/Genome/{sample}_{cell_line}.sorted.bam"),
    output:
        bam = temp(join(scratch_dir,"subsample_data/{sample}_{cell_line}/pseudo_bulk_genome_map_subsample_{n_reads}.bam")),
        bai = temp(join(scratch_dir,"subsample_data/{sample}_{cell_line}/pseudo_bulk_genome_map_subsample_{n_reads}.bam.bai"))
    resources:
        cpus_per_task=1,
        mem_mb=8000
    params:
        seed = config['random_seed'],
    shell:
        """
        mkdir -p $(dirname {output.bam})
        samtools view {input} | cut -f1 | sort | uniq > {output.bam}.read_ids
        shuf -n $(numfmt --from=si  {wildcards.n_reads}) --random-source=<(yes {params.seed}) {output.bam}.read_ids  > {output.bam}.read_ids.subsampled
        rm {output.bam}.read_ids
        samtools view -b -N {output.bam}.read_ids.subsampled -F2304 {input} > {output.bam}
        samtools index {output.bam}
        """

rule pseudobulk_alignQC_analysis_subsample:
    input: 
        fa=config['reference']['genome'], 
        anno=config['reference']['gtf_gz'],
        genome_bam = join(scratch_dir,"subsample_data/{sample}_{cell_line}/pseudo_bulk_genome_map_subsample_1M.bam"),
        genome_bai = join(scratch_dir,"subsample_data/{sample}_{cell_line}/pseudo_bulk_genome_map_subsample_1M.bam.bai")
    resources:
        cpus_per_task=4,
        mem_mb=200000,
        slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au" #
    output:
        output = directory(join(results_dir, "PseudoBulkQC/AlignQC/{sample}_{cell_line}/")),
        tmp_dir = temp(directory(join(scratch_dir, "alignQC_pseudobulk_tmp","{sample}_{cell_line}")))
    priority: 10
    container: "docker://vacation/alignqc"
    retries: 3
    shell:
        """
        mkdir -p $(dirname {output.output})
        mkdir -p {output.tmp_dir}
        alignqc analyze {input.genome_bam} \
            -g {input.fa} \
            --gtf {input.anno} \
            --output_folder {output.output} \
            --threads {resources.cpus_per_task} \
            --specific_tempdir {output.tmp_dir}
        """

rule pseudobulk_qc:
    input:
        rules.pseudobulk_read_count.output[0],
        expand([
                rules.pseudobulk_alignQC_analysis_subsample.output[0]
            ], 
            cell_line = config['cell_line_list'],
            sample = config['sample_id']
        )
        
    output:
        touch(join(flag_dir, "pseudobulk_qc.done"))