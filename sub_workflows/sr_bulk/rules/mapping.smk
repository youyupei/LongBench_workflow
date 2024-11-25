# rule for mapping
from os.path import join

# Step 1: Trim reads using fastp
rule fastp:
    input:
        R1 = join(input_fastq_dirs['ill_bulk'], "{cell_line}_R1.fastq.gz"),
        R2 = join(input_fastq_dirs['ill_bulk'], "{cell_line}_R2.fastq.gz")
    output:
        R1=os.path.join(scratch_dir, "trimmed_fq/{cell_line}_R1.fastq.gz"),
        R2=os.path.join(scratch_dir, "trimmed_fq/{cell_line}_R2.fastq.gz"),
        html=join(results_dir, "qc/fastp/{cell_line}.html"),
        json=join(results_dir, "qc/fastp/{cell_line}.json")
    conda:
        config['conda']['fastp']
    resources:
        cpus_per_task=16,
        mem_mb=32000
    shell:
        """
        mkdir -p $(dirname {output[0]})
        fastp -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2}  --thread {resources.cpus_per_task} -j {output.json} -h {output.html}
        """

rule split_fa:
    input:
        genome = config["reference"]["genome"]
    output:
        temp(join(scratch_dir, "subjunc/genome.fa"))
    resources:
        cpus_per_task=4,
        mem_mb=32000
    shell:
        """
        mkdir -p $(dirname {output})
        awk '
            {{
                if ($0 ~ /^>/) {{
                    # Print header lines as-is
                    print $0
                }} else {{
                    # For sequence lines, split into 1000-character chunks
                    while (length($0) > 500) {{
                        print substr($0, 1, 500)
                        $0 = substr($0, 501)
                    }}
                    print $0  # Print any remaining sequence
                }}
            }}' {input.genome} > {output}
        """

rule subjunc_index:
    input:
        genome = rules.split_fa.output
    output:
        index = directory(join(config["output_path"], "subjunc/genome_index"))
    resources:
        cpus_per_task=4,
        mem_mb=32000
    script: join(config['sub_wf_dir'], 'scripts/build_index.R')
        
rule subjunc_mapping:
    input:
        R1 = rules.fastp.output.R1,
        R2 = rules.fastp.output.R2,
        index = rules.subjunc_index.output.index,
        gtf = config["reference"]["gtf"]
    output:
        bam = temp(join(config["output_path"], "subjunc/bam/{cell_line}.bam"))
    resources:
        cpus_per_task=16,
        mem_mb=64000
    script: join(config['sub_wf_dir'], 'scripts/align_reads.R')

rule sort_and_index_bam:
    input:
        bam = rules.subjunc_mapping.output.bam
    output:
        sorted_bam = join(config["output_path"], "subjunc/bam/{cell_line}.sorted.bam"),
        bai = join(config["output_path"], "subjunc/bam/{cell_line}.sorted.bam.bai")
    resources:
        cpus_per_task=16,
        mem_mb=32000
    shell:
        """
        samtools sort -@ {resources.cpus_per_task} -o {output.sorted_bam} {input.bam}
        samtools index {output.sorted_bam}
        """

rule mapping:
    input:
        expand(
            [
                rules.sort_and_index_bam.output.sorted_bam,
                rules.sort_and_index_bam.output.bai
            ],
            cell_line = config['cell_lines']
        )
    output:
        touch(join(config["output_path"], ".flag/mapping.done"))