import textwrap, os

rule Subsample_bam_for_qc:
    input:
        bam = os.path.join(results_dir, "flames_out/{x}/align2genome.bam"),
        bai = os.path.join(results_dir, "flames_out/{x}/align2genome.bam.bai"),
    output:
        bam = os.path.join(scratch_dir,"subsample_data/{x}/align2genome_{n_reads}.bam"),
        bai = os.path.join(scratch_dir,"subsample_data/{x}/align2genome_{n_reads}.bam.bai")
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
        samtools view -b -N {output.bam}.read_ids.subsampled {input} > {output.bam}
        samtools index {output.bam}
        """


rule bam_F2304_filtering:
    input:
        os.path.join(scratch_dir,"subsample_data/{x}/align2genome_3M.bam")
    output:
        os.path.join(scratch_dir,"subsample_data/{x}/align2genome.3M.F2304.bam"),
        os.path.join(scratch_dir,"subsample_data/{x}/align2genome.3M.F2304.bam.bai")
    resources:
        cpus_per_task=1,
        mem_mb=8000
    shell:
        """
        mkdir -p $(dirname {output})
        samtools view -b -F 2304 {input} > {output}
        samtools index {output}

        """

rule _internal_priming_identifier_single_run:
    input: 
        bam = os.path.join(scratch_dir,"subsample_data/{x}/align2genome.3M.F2304.bam"),
        bai = os.path.join(scratch_dir,"subsample_data/{x}/align2genome.3M.F2304.bam.bai"),
        gtf = config['reference']['gtf'],
        genome = config['reference']['genome']
    output:
        summary=join(results_dir, "int_prim_analysis/{x}_summary.txt"),
        bam = join(scratch_dir, "int_prim_analysis/{x}_IP_tag_added.bam"),
        bai = join(scratch_dir, "int_prim_analysis/{x}_IP_tag_added.bam.bai")
    resources:
        cpus_per_task=32,
        mem_mb=32000
    params:
        python_script="/home/users/allstaff/you.yu/github/PrimeSpotter/PrimeSpotter/PrimeSpotter.py"
    conda:
        config["conda"]["PrimeSpotter"]
    shell:
        """
        mkdir -p $(dirname {output.summary})
        mkdir -p $(dirname {output.bam})
        # module load samtools
        
        python3 {params.python_script} --bam_file {input.bam} \
                                        --gtf_file {input.gtf} \
                                        --output-summary {output.summary} \
                                        --genome-ref {input.genome} \
                                        --processes {resources.cpus_per_task}| samtools view -S -b | samtools sort > {output.bam}
        samtools index {output.bam}
        """


rule internal_priming_identifier:
    input:
        expand(join(results_dir, "int_prim_analysis/{x}_summary.txt"),
                x = config['sample_id'])