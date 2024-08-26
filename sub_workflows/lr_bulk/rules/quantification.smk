results_dir = config["output_path"]
main_conda = config["conda"]["main"]

# PseudoBulk quantificaiton
rule run_salmon:
    input:
        expand(results_dir + "/.flag/{run_id}_{cell_line}.salmon.done",
                run_id = config["sample_id"],
                cell_line = config["cell_lines"])
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
        main_conda
    resources:
        cpus_per_task=32,
        mem_mb=32000,
        slurm_extra="--mail-type=END,FAIL --mail-user=you.yu@wehi.edu.au"
    params:
        lib_type = "A", # auto
        gibbs_samples = 40
    shell:
        "salmon quant -t {input.ref} -l {params.lib_type} -a {input.bam} -p {resources.cpus_per_task}  -o {output.out_dir} --ont --numBootstraps {params.gibbs_samples}"
