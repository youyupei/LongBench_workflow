results_dir = config["output_path"]
main_conda = config["conda"]["main"]

# quantificaiton
rule run_quantification:
    input:
        expand([
            os.path.join(results_dir,"salmon_output/{sample}/{cell_line}"),
            os.path.join(results_dir,"oarfish_nocov_output/{sample}/{cell_line}"),
            os.path.join(results_dir,"oarfish_cov_output/{sample}/{cell_line}")
            ],
            sample = config["sample_id"],
            cell_line = config["cell_lines"]),
    output:
        touch(results_dir + "/.flag/run_salmon.done")

# Salmon
rule salmon:
    input:
        bam = results_dir + "/TranscriptAlignment/{sample}_{cell_line}.bam",
        ref = config["reference"]["transcript"]
    output:
        out_dir = directory(os.path.join(results_dir,"salmon_output/{sample}/{cell_line}"))
    conda:
        main_conda
    resources:
        cpus_per_task=32,
        mem_mb=100000,
        slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
    params:
        lib_type = "A", # auto
        gibbs_samples = 50
    shell:
        "salmon quant -t {input.ref} -l {params.lib_type} -a {input.bam} -p {resources.cpus_per_task}  -o {output.out_dir} --ont --numBootstraps {params.gibbs_samples}"



rule oarfish_no_cov:
    input:
        bam = results_dir + "/TranscriptAlignment/{sample}_{cell_line}.bam",
        ref = config["reference"]["transcript"]
    output:
        out_dir_nocov = directory(os.path.join(results_dir,"oarfish_nocov_output/{sample}/{cell_line}")),
    conda:
        config["conda"]["oarfish"]
    resources:
        cpus_per_task=16,
        mem_mb=64000,
        slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
    params:
        lib_type = "A", # auto
        gibbs_samples = 50
    shell:
        """
        mkdir -p {output.out_dir_nocov}
        oarfish --alignments {input.bam} --threads {resources.cpus_per_task} --output {output.out_dir_nocov}/ -d . --filter-group no-filters --num-bootstraps 50
        """

rule oarfish_cov:
    input:
        bam = results_dir + "/TranscriptAlignment/{sample}_{cell_line}.bam",
        ref = config["reference"]["transcript"]
    output:
        out_dir_cov = directory(os.path.join(results_dir,"oarfish_cov_output/{sample}/{cell_line}"))
    conda:
        config["conda"]["oarfish"]
    resources:
        cpus_per_task=16,
        mem_mb=64000,
        slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
    params:
        lib_type = "A", # auto
        gibbs_samples = 50
    shell:
        """
        mkdir -p {output.out_dir_cov}
        oarfish --alignments {input.bam} --threads {resources.cpus_per_task} --output {output.out_dir_cov}/ --model-coverage  -d . --filter-group no-filters --num-bootstraps 50
        """