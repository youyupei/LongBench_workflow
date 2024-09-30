# setup dirs
results_dir = config["output_path"]
main_conda = config["conda"]["main"]

# Main work flow
rule blaze:
    input:
        fastq = lambda wildcards: config["samples_fastq_dir"][wildcards.sample]
    output:
        out_fastq = config["output_path"] + "/flames_out/{sample}/matched_reads.fastq",
        other_output = [
            config["output_path"] + "/flames_out/{sample}/putative_bc.csv",
            config["output_path"] + "/flames_out/{sample}/whitelist.csv",
            config["output_path"] + "/flames_out/{sample}/summary.txt"
        ]
    params:
        manual_count_thres = lambda wildcards: config["count_thres"][wildcards.sample] if config["count_thres"][wildcards.sample] else "",
        expected_cells = lambda wildcards: config["expected_cells"][wildcards.sample] if config["expected_cells"][wildcards.sample] else "",
        whitelist = lambda wildcards: config["whitelist_10x"][wildcards.sample] if config["whitelist_10x"][wildcards.sample] else "",
        force_cell_n = lambda wildcards: config["force_cells"][wildcards.sample] if config["force_cells"][wildcards.sample] else ""
    conda:
        main_conda
    resources:
        cpus_per_task=32,
        mem_mb=100000,
        slurm_extra="--mail-type=END,FAIL --mail-user=you.yu@wehi.edu.au"
    shell:
        """
        out_fn=$(dirname {output.out_fastq})/
        mkdir -p $out_fn

        # check if the whitelist is empty
        if [ "{params.whitelist}" == "" ]; then
            # no whitelist
            full_whitelist_arg=""
        else
            full_whitelist_arg="--full-bc-whitelist {params.whitelist}"
        fi


        # Specify the demultiplexing option (overwrite priority: --force-cells >  --count-threshold > --expect-cells )
        if [ "{params.force_cell_n}" != "" ]; then
            arg_cell_calling="--force-cells {params.force_cell_n}"
        elif [ "{params.manual_count_thres}" != "" ]; then
            arg_cell_calling="--count-threshold {params.manual_count_thres}"
        elif [ "{params.expected_cells}" != "" ]; then
            arg_cell_calling="--expect-cells {params.expected_cells}"
        else
            arg_cell_calling=""
        fi

        blaze \\
        $full_whitelist_arg  \\
        $arg_cell_calling \\
        --output-prefix $out_fn --output-fastq matched_reads.fastq --threads {resources.cpus_per_task} \\
        --minimal_stdout {input.fastq}
        """

rule flames_genome_mapping:
    input:
        fastq = rules.blaze.output.out_fastq,
        config_file = config['flames_config']['genome_mapping'],
        gtf = config['reference']['gtf'],
        genome = config['reference']['genome']
    output: 
        bam = os.path.join(results_dir,"flames_out/{sample}/align2genome.bam"),
        bai = os.path.join(results_dir,"flames_out/{sample}/align2genome.bam.bai")
    params:
        minimap2 = config["software"]["minimap2"],
        k8 = os.path.dirname(config["software"]["minimap2"]) + '/k8'
    resources:
        cpus_per_task=32,
        mem_mb=500000,
        slurm_extra="--mail-type=END,FAIL --mail-user=you.yu@wehi.edu.au"
    shell:
        """
        out_dir=$(dirname {input.fastq})

        Rscript -e "
            library(FLAMES)
            config_file <- '{input.config_file}'
            fastq_flames= '{input.fastq}'
            outdir = '$out_dir'
            GTF = '{input.gtf}'
            genome = '{input.genome}'

            sce <- sc_long_pipeline(fastq=fastq_flames, 
                                    outdir=outdir, 
                                    annot=GTF, 
                                    genome_fa=genome, 
                                    config_file=config_file,
                                    barcodes_file='not used',
                                    minimap2='{params.minimap2}',
                                    k8='{params.k8}')
            "
        """

rule flames_gene_quant:
    input:
        bam = rules.flames_genome_mapping.output,
        fastq = rules.blaze.output.out_fastq,
        config_file = config['flames_config']['gene_quantification'],
        gtf = config['reference']['gtf'],
        genome = config['reference']['genome']
    output:
        fastq = os.path.join(results_dir,"flames_out/{sample}/matched_reads_dedup.fastq"),
        csv = os.path.join(results_dir,"flames_out/{sample}/gene_count.csv")
    resources:
        cpus_per_task=32,
        mem_mb=500000,
        slurm_extra="--mail-type=END,FAIL --mail-user=you.yu@wehi.edu.au"
    params:
        minimap2 = config["software"]["minimap2"],
        k8 = os.path.dirname(config["software"]["minimap2"]) + '/k8'
    shell:
        """
        out_dir=$(dirname {input.fastq})

        Rscript -e "
            library(FLAMES)
            config_file <- '{input.config_file}'
            fastq_flames= '{input.fastq}'
            outdir = '$out_dir'
            GTF = '{input.gtf}'
            genome = '{input.genome}'

            sce <- sc_long_pipeline(fastq=fastq_flames, 
                                    outdir=outdir, 
                                    annot=GTF, 
                                    genome_fa=genome, 
                                    config_file=config_file,
                                    barcodes_file='not used',
                                    minimap2='{params.minimap2}',
                                    k8='{params.k8}')
            "
        """


rule flames_trans_map_and_quant:
    input:
        fastq = rules.flames_gene_quant.output.fastq,
        config_file = config['flames_config']['flames_trans_map_and_quant'],
        gtf = config['reference']['gtf'],
        genome = config['reference']['genome']
    output:
        bam = os.path.join(results_dir,"flames_out/{sample}/realign2transcript.bam"),
        bai = os.path.join(results_dir,"flames_out/{sample}/realign2transcript.bam.bai"),
        trans_quant = os.path.join(results_dir,"flames_out/{sample}/transcript_count.csv.gz"),
        flag = touch(results_dir + "/.flag/flames_{sample}.done")
    resources:
        cpus_per_task=32,
        mem_mb=500000,
        slurm_extra="--mail-type=END,FAIL --mail-user=you.yu@wehi.edu.au"
    params:
        minimap2 = config["software"]["minimap2"],
        k8 = os.path.dirname(config["software"]["minimap2"]) + '/k8'
    shell:
        """
        out_dir=$(dirname {input.fastq})

        Rscript -e "
            library(FLAMES)
            config_file <- '{input.config_file}'
            fastq_flames= '{input.fastq}'
            outdir = '$out_dir'
            GTF = '{input.gtf}'
            genome = '{input.genome}'

            sce <- sc_long_pipeline(fastq=fastq_flames, 
                                    outdir=outdir, 
                                    annot=GTF, 
                                    genome_fa=genome, 
                                    config_file=config_file,
                                    barcodes_file='not used',
                                    minimap2='{params.minimap2}',
                                    k8='{params.k8}')
            "
        """


rule flames:
    input:
        expand(os.path.join(results_dir,"flames_out/{sample}/transcript_count.csv.gz"), sample=config["sample_id"])
    output:
        touch(results_dir + "/.flag/flames.done")