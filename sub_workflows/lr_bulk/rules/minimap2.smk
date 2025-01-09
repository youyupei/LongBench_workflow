results_dir = config["output_path"]
barcode_list = config['barcode_list']
cell_line_to_barcode = {cl: bc for d in barcode_list for bc, cl in d.items()}

wildcard_constraints:
    cell_line='|'.join(cell_line_to_barcode.keys())

rule lr_bulk_minimap2_transcript:
    priority: 10
    input:
        fastq = lambda w: glob.glob(os.path.join(config['samples_fastq_dir'][w.sample],f"{w.cell_line}.fastq*"))[0],
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
        {params.minimap2} {params.minimap2_trans_options} -t {resources.cpus_per_task} {input.ref}  {input.fastq}  | samtools view -bS - > {output.bam}
        """

rule lr_bulk_minimap2_Genome:
    priority: 10
    input:
        fastq = lambda w: glob.glob(os.path.join(config['samples_fastq_dir'][w.sample],f"{w.cell_line}.fastq*"))[0],
        ref = config['reference']['genome']
    output:
        bam = temp(results_dir + "/GenomeAlignment/{sample}_{cell_line}.bam")
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

rule lr_bulk_transcript_sam_flagstat:
    input:
        bam = results_dir + "/TranscriptAlignment/{sample}_{cell_line}.bam"
    output:
        flagstat = results_dir + "/TranscriptAlignment/{sample}_{cell_line}.flagstat"
    shell:
        "samtools flagstat {input.bam} > {output.flagstat}"

rule lr_bulk_genome_sam_flagstat:
    input:
        bam = results_dir + "/GenomeAlignment/{sample}_{cell_line}.sorted.bam"
    output:
        flagstat = results_dir + "/GenomeAlignment/{sample}_{cell_line}.flagstat"
    shell:
        "samtools flagstat {input.bam} > {output.flagstat}"
rule run_all_mapping:
    input:
        expand([results_dir + "/TranscriptAlignment/{sample}_{cell_line}.bam", 
                results_dir + "/GenomeAlignment/{sample}_{cell_line}.sorted.bam",
                results_dir + "/TranscriptAlignment/{sample}_{cell_line}.flagstat",
                results_dir + "/GenomeAlignment/{sample}_{cell_line}.flagstat"], 
                sample = config['sample_id'],
                cell_line = config['cell_lines'])
    output:
        touch(results_dir + "/.flag/run_all_mapping.done")