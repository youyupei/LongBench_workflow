results_dir = config["output_path"]

rule split_pooled_bam:
    """
    Split the pooled bam file into individual bam files
    """
    input:
        config['samples_raw_bam_dir']['ont_bulk']
    output:
        touch(".flag/{x}.bam.split_by_bc_done")
    resources:
        cpus_per_task=16,
        mem_mb=32000
    params:
        bc_dict = config["barcode_index"].keys()
        cell_line = config["cell_line"].values()
    shell:
        """
        # module load samtools
        cd $(dirname {input})
        samtools split -d BC --threads {resources.cpus_per_task} {input} 
        """

rule split_pooled_bam:
    """
    Split the pooled bam file into individual bam files
    """
    input:
        '{x}.bam'
    output:
        touch(".flag/{x}.bam.split_by_bc_done")
    resources:
        cpus_per_task=16,
        mem_mb=32000
    params:
        bc_dict = config["barcode_index"].keys()
        cell_line = config["cell_line"].values()
    shell:
        """
        # module load samtools
        cd $(dirname {input})
        samtools split -d BC --threads {resources.cpus_per_task} {input} 
        """