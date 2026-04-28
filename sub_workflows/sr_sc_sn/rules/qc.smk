# Required variables (set in sr_sc_sn Snakefile):
## config["samples_fastq_dir"]  — dict: sample_name -> fastq directory
## config["sample_name"]        — list of sample names
## results_dir

######## Count bases ########
rule count_bases_in_fastq_sr_sc_sn:
    """Count total bases sequenced from raw FASTQ directory (all .fastq.gz files)."""
    input:
        lambda wildcards: config["samples_fastq_dir"][wildcards.sample_name]
    output:
        os.path.join(results_dir, "qc/base_counts/{sample_name}.total_bases")
    resources:
        cpus_per_task=1,
        mem_mb=4000
    shell:
        """
        mkdir -p $(dirname {output})
        zcat {input}/*.fastq.gz | awk 'NR%4==2{{sum+=length($0)}} END{{print sum}}' > {output}
        """
