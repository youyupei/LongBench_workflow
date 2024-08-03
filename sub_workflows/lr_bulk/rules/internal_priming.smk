import textwrap, os
results_dir = config["output_path"]

rule run_internal_priming_analysis:
    input:
        ".flag/{x}_run_primspotter.done"
    output:
        touch(".flag/{x}_internal_priming.done")

rule _internal_priming_identifier:
    input: 
        bam = os.path.join(results_dir, "flames_out/{x}/align2genome.bam"),
        bai = os.path.join(results_dir, "flames_out/{x}/align2genome.bam.bai"),
        github = ".flag/gitrepo_youyupei_PrimeSpotter.commit",
        gtf = config['reference']['gtf'],
        genome = config['reference']['genome']
    output:
        flag=touch(".flag/{x}_run_primspotter.done"),
        summary=os.path.join(results_dir, "int_prim_analysis/{x}_summary.txt"),
        bam = os.path.join(results_dir, "int_prim_analysis/{x}_IP_tag_added.bam"),
        bai = os.path.join(results_dir, "int_prim_analysis/{x}_IP_tag_added.bam.bai")
    resources:
        cpus_per_task=32,
        mem_mb=32000,
        slurm_extra="--mail-type=END,FAIL --mail-user=you.yu@wehi.edu.au"
    params:
        python_script="git_repo/PrimeSpotter/PrimeSpotter/PrimeSpotter.py",
        python_dir="/stornext/Home/data/allstaff/y/you.yu/.cache/R/basilisk/1.14.0/FLAMES/1.9.1/flames_env/bin/python3"
    shell:
        """
        mkdir -p $(dirname {output.summary})
        module load samtools
        
        {params.python_dir} {params.python_script} --bam_file {input.bam} \
                                        --gtf_file {input.gtf} \
                                        --output-summary {output.summary} \
                                        --genome-ref {input.genome} \
                                        --processes {resources.cpus_per_task}| samtools view -S -b | samtools sort > {output.bam}
        samtools index {output.bam}
        """


# rule internal_priming_site_analysis:
#     """
#     Find the potential internal priming sites using Genome and annotation
#     """
#     output:
#         #".flag/{y}_{x}_internal_priming.done",
#         os.path.join(results_dir, "int_prim_analysis/internal_priming_site_analysis/int_prim_site_summary.csv")
#     resources:
#         cpus_per_task=1,
#         mem_mb=32000,
#         slurm_extra="--mail-type=END,FAIL --mail-user=you.yu@wehi.edu.au"
#     params:
#         python_script="scripts/internal_priming_site_analysis.py",
#         python_dir="/stornext/Home/data/allstaff/y/you.yu/.cache/R/basilisk/1.14.0/FLAMES/1.9.1/flames_env/bin/python3"
#     shell:
#         """
#         mkdir -p $(dirname {output[0]})
#         {params.python_dir} {params.python_script} 
#         """


# rule internal_priming_gene_count:
#     input:
#         os.path.join(results_dir, "int_prim_analysis/{sample}/{sample}_{y}.bam"),
#         os.path.join(results_dir, "int_prim_analysis/{sample}/{sample}_{y}.bam.bai")
#     output:
#         os.path.join(results_dir, "int_prim_analysis/{sample}/{sample}_{y}_gene_count.csv")
#     params:
#         python_dir="/stornext/Home/data/allstaff/y/you.yu/.cache/R/basilisk/1.14.0/FLAMES/1.9.1/flames_env/bin/python3",
#         flames_module_dir = "/home/users/allstaff/you.yu/github/FLAMES/inst/python"
#     resources:
#         cpus_per_task=32,
#         mem_mb=32000,
#         slurm_extra="--mail-type=END,FAIL --mail-user=you.yu@wehi.edu.au"
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

