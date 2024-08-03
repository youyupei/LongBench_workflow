results_dir = config["output_path"]
main_conda_env = config['conda_config_path'] + "/main.yaml"

# PseudoBulk quantificaiton
rule run_salmon:
    input:
        expand(results_dir + "/.flag/{run_id}_{cell_line}.salmon.done",
                run_id = config["run_id"],
                cell_line = config["cell_line_list"])
    output:
        touch(results_dir + "/.flag/run_salmon.done")

# Salmon
rule salmon:
    input:
        bam = results_dir + "/Alignment/ont_bulk_{cell_line}.bam",
        ref = config["reference"]["transcript"]
    output:
        out_dir = directory(os.path.join(results_dir,"salmon_output/{run_id}_{cell_line}")),
        out_flag = touch(results_dir + "/.flag/{run_id}_{cell_line}.salmon.done")
    conda:
        main_conda_env
    resources:
        cpus_per_task=32,
        mem_mb=32000,
        slurm_extra="--mail-type=END,FAIL --mail-user=you.yu@wehi.edu.au"
    params:
        lib_type = "A" # auto
    shell:
        "salmon quant -t {input.ref} -l {params.lib_type} -a {input.bam} -p {resources.cpus_per_task}  -o {output.out_dir} --ont "


# rule picard_SortSam:
#     """
#     sort sam by name
#     """
#     input:
#         "{x}.bam"
#     output:
#         "{x}.name_sorted.bam"
#     params:
#         picard_dir = "/home/users/allstaff/you.yu/project/software/picard.jar",
#         tmp_dir="/home/users/allstaff/you.yu/LongBench/tmp"
#     resources:
#         cpus_per_task=8,
#         mem_mb=32000
#     shell:
#         """
#         module load picard-tools/2.26.11 # this is just for loading the dependency, check the params.picard_dir for the actual path
#         java -jar {params.picard_dir} SortSam \
#                     INPUT={input} \
#                     OUTPUT={output} \
#                     TMP_DIR={params.tmp_dir} \
#                     SORT_ORDER=queryname
#         """
