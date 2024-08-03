results_dir = config["output_path"]
barcode_list = config['barcode_list']
cell_line_to_barcode = {cl: bc for d in barcode_list for bc, cl in d.items()}


rule run_all_mapping:
    input:
        expand(results_dir + "/.flag/ont_bulk_{cell_line}_AlignmentDone.flag", cell_line = cell_line_to_barcode.keys())
        #expand(results_dir + "/Alignment/ont_bulk_{cell_line}.bam", cell_line = cell_line_to_barcode.keys())
        
    output:
        touch(results_dir + "/.flag/run_all_mapping.done")


rule ont_bulk_minimap2_transcript:
    priority: 10
    input:
        in_bam = lambda w: os.path.join(config['samples_raw_bam_dir']['ont_bulk'],
                            f"bulk_{w.pool}_SQK-PCB114-24_barcode0{cell_line_to_barcode[w.cell_line]}.bam"),
        ref = config['reference']['transcript']
    output:
        bam = results_dir + "/Alignment/ont_bulk_{pool}_{cell_line}.bam"
    resources:
        cpus_per_task=8,
        mem_gb=64
    params:
        minimap2 = "/home/users/allstaff/you.yu/LongBench/software/minimap2-2.28_x64-linux/minimap2"
    shell:
        """
        module load samtools
        samtools fastq -T* {input.in_bam} | {params.minimap2} -ax lr:hq -t {resources.cpus_per_task} {input.ref}  - | samtools view -bS - > {output.bam}
        """

rule ont_bulk_cat_unsorted_bam:
    input:
        # lambda w: expand(results_dir + "/Alignment/ont_bulk_{p}_{cl}.bam",
        #                 p = ["pool", "pool2", "pool3"], cl = [w.cell_line])
        expand(results_dir + "/Alignment/ont_bulk_{pool}_{cell_line}.bam", 
                        cell_line = cell_line_to_barcode.keys(),
                        pool = ["pool", "pool2", "pool3"])
    output:
        results_dir + "/Alignment/ont_bulk_{cell_line}.bam"
    resources:
        cpus_per_task=8,
        mem_gb=64
    shell:
        """
        module load samtools
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
