# Step 2: Run salmon on the trimmed reads
## Step 2.1: Index the transcriptome
rule salmon_index:
    input: config['reference']['transcript']
    output: directory(join(results_dir, "salmon/index"))
    resources:
        cpus_per_task=16,
        mem_mb=32000
    conda:
        config["conda"]["main"]
    shell:
        """
        salmon index -t {input} -i {output} -k 31
        """

rule salmon_quant:
    input:
        R1 = rules.fastp.output.R1,
        R2 = rules.fastp.output.R2,
        index = rules.salmon_index.output
    output:
         directory(join(results_dir, "salmon/salmon_quant/{cell_line}"))
    conda:
        config["conda"]["main"]
    shell:
        """
        mkdir -p $(dirname {output})
        salmon quant -i {input.index} \
                    -l A \
                    -1 {input.R1} \
                    -2 {input.R2} \
                    --validateMappings \
                    -o {output} \
                    -p 16 \
                    --numBootstraps 50
        """

# Entire quantification worflow head
rule salmon:
    input:
        expand(
            [
                rules.salmon_quant.output[0]
            ],
            cell_line = config['cell_lines']
        )
    output:
        touch(join(results_dir, "qc/.flag/qc.done"))