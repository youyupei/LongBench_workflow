results_dir = config["output_path"]
barcode_list = config['barcode_list']
cell_line_to_barcode = {cl: bc for d in barcode_list for bc, cl in d.items()}

wildcard_constraints:
    cell_line='|'.join(cell_line_to_barcode.keys())


rule run_all_mapping:
    input:
        expand(results_dir + "/.flag/ont_bulk_{cell_line}_AlignmentDone.flag", cell_line = cell_line_to_barcode.keys()) 
    output:
        touch(results_dir + "/.flag/run_all_mapping.done")


rule lr_bulk_minimap2_transcript:
    priority: 10
    input:
        fastq = lambda w: os.path.join(config['samples_fastq_dir'][w.sample], "{cell_line}.fastq"),
        ref = config['reference']['transcript']
    output:
        bam = results_dir + "/TranscriptAlignment/{sample}_{cell_line}.bam"
    resources:
        cpus_per_task=16,
        mem_mb=64000
    params:
        minimap2 = config["software"]["minimap2"],
        minimap2_trans_options = lambda w: config["minimap2_trans_options"][w.sample]
    shell:
        """
        {params.minimap2} {params.minimap2_trans_options} -t {resources.cpus_per_task} {input.ref}  - | samtools view -bS - > {output.bam}
        """

rule lr_bulk_minimap2_Genome:
    priority: 10
    input:
        fastq = lambda w: os.path.join(config['samples_fastq_dir'][w.sample], "{cell_line}.fastq"),
        ref = config['reference']['genome']
    output:
        bam = results_dir + "/GenomeAlignment/{sample}_{cell_line}.bam"
    resources:
        cpus_per_task=16,
        mem_mb=64000
    params:
        minimap2 = config["software"]["minimap2"],
        minimap2_genome_options = lambda w: config["minimap2_genome_options"][w.sample]
    shell:
        """
        {params.minimap2} {params.minimap2_genome_options} -t {resources.cpus_per_task} {input.ref}  {input.fastq} | samtools view -bS - > {output.bam}
        """

rule ont_bulk_cat_unsorted_bam:
    input:
        [results_dir + f"/TranscriptAlignment/ont_bulk_{pool}_" + "{cell_line}.bam" for pool  in ["pool", "pool2", "pool3"]]
    output:
        results_dir + "/TranscriptAlignment/ont_bulk_{cell_line}.bam"
    wildcard_constraints:
        cell_line='|'.join(cell_line_to_barcode.keys())
    resources:
        cpus_per_task=8,
        mem_gb=64
    shell:
        """
        # module load samtools
        samtools cat -@ {resources.cpus_per_task} -o {output} {input}
        """

rule ont_bulk_clean_up:
    input:
        results_dir + "/Alignment/ont_bulk_{cell_line}.bam"
    output: 
        touch(results_dir + "/.flag/ont_bulk_{cell_line}_AlignmentDone.flag")
    localrule: True
    params:
        bam = lambda w: expand(results_dir + "/Alignment/ont_bulk_{pool}_{x}.bam",
                        pool = ["pool", "pool2", "pool3"], x = [w.cell_line])
    shell:
        """
        rm -f {params.bam}
        """
