configfile: "config/config.yaml"
results_dir = config["output_path"]

# Main work flow
rule blaze:
    input:
        fastq = lambda wildcards: config["samples_fastq_dir"][wildcards.sample]
    output:
        flag = ".flag/blaze_{sample}.done",
        out_fastq = config["output_path"] + "/flames_out/{sample}/matched_reads.fastq"
    params:
        blaze_path = "/stornext/Home/data/allstaff/y/you.yu/.cache/R/basilisk/1.14.0/FLAMES/1.9.1/flames_env/bin",
        manual_count_thres = lambda wildcards: config["count_thres"][wildcards.sample] if config["count_thres"][wildcards.sample] else "",
        expected_cells = lambda wildcards: config["expected_cells"][wildcards.sample] if config["expected_cells"][wildcards.sample] else "",
        whitelist = lambda wildcards: config["whitelist_10x"][wildcards.sample] if config["whitelist_10x"][wildcards.sample] else "",
        force_cell_n = lambda wildcards: config["force_cells"][wildcards.sample] if config["force_cells"][wildcards.sample] else ""
    resources:
        cpus_per_task=32,
        mem_mb=100000,
        slurm_extra="--mail-type=END,FAIL --mail-user=you.yu@wehi.edu.au"
    shell:
        """
        out_fn=$(dirname {output.out_fastq})/
        mkdir -p $out_fn
        mkdir -p $(dirname {output.flag})

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

        {params.blaze_path}/blaze \\
        $full_whitelist_arg  \\
        $arg_cell_calling \\
        --output-prefix $out_fn --output-fastq matched_reads.fastq --threads {resources.cpus_per_task} \\
        --minimal_stdout {input.fastq}

        touch {output.flag}
        """

rule flame_run_no_identification:
    input:
        flag = ".flag/blaze_{sample}.done",
        config_file = 'config/flames.json',
        fastq = config["output_path"] + "/flames_out/{sample}/matched_reads.fastq"
    output:
        flag = ".flag/flames_{sample}.done",
        output_list = [
            os.path.join(results_dir,"flames_out/{sample}/align2genome.bam"),
            os.path.join(results_dir,"flames_out/{sample}/realign2transcript.bam")
        ]
    resources:
        cpus_per_task=32,
        mem_mb=500000,
        slurm_extra="--mail-type=END,FAIL --mail-user=you.yu@wehi.edu.au"
    shell:
        """
        module load minimap2
        module load samtools

        out_dir=$(dirname {input.fastq})
        mkdir -p $out_dir
        mkdir -p $(dirname {output.flag})

        Rscript -e "
            library(FLAMES)
            config_file <- '{input.config_file}'
            fastq_flames= '$out_dir/matched_reads.fastq'
            outdir = '$out_dir'
            GTF = '{config[reference][gtf]}'
            genome = '{config[reference][genome]}'

            sce <- sc_long_pipeline(fastq=fastq_flames, 
                                    outdir=outdir, 
                                    annot=GTF, 
                                    genome_fa=genome, 
                                    config_file=config_file,
                                    barcodes_file='not used')
            "

        touch {output.flag}
        """


rule flame_run_no_identification_tmp:
    """
    This is a tmp rule for regenerating the flames output with no identification (assuming 
    previous run has generated the gene quantification and deduplicated fastq)
    """
    input:
        config_file = 'config/flames_trans_quant_only.json',
        fastq = config["output_path"] + "/flames_out/{sample}/matched_reads_dedup.fastq"
    output:
        flag = ".flag/flames_{sample}_tmp.done"
    resources:
        cpus_per_task=32,
        mem_mb=500000,
        slurm_extra="--mail-type=END,FAIL --mail-user=you.yu@wehi.edu.au"
    shell:
        """
        module load minimap2
        module load samtools

        out_dir=$(dirname {input.fastq})
        mkdir -p $out_dir
        mkdir -p $(dirname {output.flag})

        Rscript -e "
            library(FLAMES)
            config_file <- '{input.config_file}'
            fastq_flames= '$out_dir/matched_reads.fastq'
            outdir = '$out_dir'
            GTF = '{config[reference][gtf]}'
            genome = '{config[reference][genome]}'

            sce <- sc_long_pipeline(fastq=fastq_flames, 
                                    outdir=outdir, 
                                    annot=GTF, 
                                    genome_fa=genome, 
                                    config_file=config_file,
                                    barcodes_file='not used')
            "

        touch {output.flag}
        """
