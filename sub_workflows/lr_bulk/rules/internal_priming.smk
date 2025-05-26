import textwrap, os
from  os.path import  join


rule run_internal_priming_analysis:
    input:
        ".flag/{x}_run_primspotter.done"
    output:
        touch(".flag/{x}_internal_priming.done")

rule _internal_priming_identifier_single_run:
    input: 
        bam =  os.path.join(scratch_dir,"subsample_data/{sample}_{cell_line}/genome_map_3M.F2304.bam"),
        gtf = config['reference']['gtf_hunman'],
        genome = config['reference']['genome']
    output:
        summary=join(results_dir, "int_prim_analysis/{sample}_{cell_line}_summary.txt"),
        bam = join(scratch_dir, "int_prim_analysis/{sample}_{cell_line}_IP_tag_added.bam"),
        bai = join(scratch_dir, "int_prim_analysis/{sample}_{cell_line}_IP_tag_added.bam.bai")
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
        expand(join(results_dir, "int_prim_analysis/{sample}_{cell_line}_summary.txt"),
                sample = config['sample_id'],
                cell_line = config['cell_lines'])

