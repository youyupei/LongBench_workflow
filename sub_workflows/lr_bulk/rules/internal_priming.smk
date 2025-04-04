import textwrap, os
from  os.path import  join


rule run_internal_priming_analysis:
    input:
        ".flag/{x}_run_primspotter.done"
    output:
        touch(".flag/{x}_internal_priming.done")

rule _internal_priming_identifier_single_run:
    input: 
        bam = os.path.join(scratch_dir,"subsample_data/{sample}_{cell_line}/genome_map_subsample_rate_0.05.bam"),
        bai = os.path.join(scratch_dir,"subsample_data/{sample}_{cell_line}/genome_map_subsample_rate_0.05.bam.bai"),
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

# rule internal_priming_site_analysis:
#     """
#     Find the potential internal priming sites using Genome and annotation
#     """
#     output:
#         join(results_dir, "int_prim_analysis/internal_priming_site_analysis/int_prim_site_summary.csv")
#     resources:
#         cpus_per_task=1,
#         mem_mb=32000
#     conda:
#         main_conda
#     params:
#         python_script="/home/users/allstaff/you.yu/github/PrimeSpotter/PrimeSpotter/internal_priming_site_analysis.py"
#     shell:
#         """
#         mkdir -p $(dirname {output[0]})
#         python3 {params.python_script} 
#         """


# rule internal_priming_gene_count:
#     input:
#         join(results_dir, "int_prim_analysis/{sample}/{sample}_{y}.bam"),
#         join(results_dir, "int_prim_analysis/{sample}/{sample}_{y}.bam.bai")
#     output:
#         join(results_dir, "int_prim_analysis/{sample}/{sample}_{y}_gene_count.csv")
#     params:
#         python_dir="/stornext/Home/data/allstaff/y/you.yu/.cache/R/basilisk/1.14.0/FLAMES/1.9.1/flames_env/bin/python3",
#         flames_module_dir = "/home/users/allstaff/you.yu/github/FLAMES/inst/python"
#     resources:
#         cpus_per_task=32,
#         mem_mb=32000,
#         slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
#     shell:
#         textwrap.dedent(
#             """
#             {params.python_dir} -c '''
#             import sys
#             sys.path.append("{params.flames_module_dir}")

#             from count_gene import quantify_gene
#             gene_count_mat, dup_read_lst, umi_lst = quantify_gene("{input[0]}", "{config[reference][gtf]}", {resources.cpus_per_task})
            
#             gene_count_mat.to_csv("{output}")
        
#             '''
#             """
#         )

